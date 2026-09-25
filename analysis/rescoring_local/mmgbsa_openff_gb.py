# -*- coding: utf-8 -*-
"""MM-GBSA corregido — MASIVE-ALS.

Corrige los fallos de la auditoría del 11/09/2026:

  1) DISOLVENTE: añade el término de Born generalizado OBC2 (GBSAOBC2Force,
     radios mbondi2, SA=ACE) a los tres sistemas. El original no tenía
     ningún modelo de disolvente y medía energía potencial en vacío.
     Equivale a lo que hace openmm/app/data/implicit/obc2.xml.

  2) GEOMETRÍA: el ligando NO se re-embebe con ETKDG. Se transfieren las
     coordenadas de la pose acoplada y solo se añaden los hidrógenos, de
     modo que la geometría de partida es idéntica entre corridas.

  3) MINIMIZACIÓN DETERMINISTA: tolerancia estricta (1 kJ/mol/nm).

  4) PLATAFORMA COHERENTE (12/09/2026): la plataforma de minimización la
     fijaba MMGBSA_FINAL_PLATFORM, pero las ENERGÍAS finales las medía
     mmgbsa_openff.potential_energy(), que crea su propio Context con
     get_platform() -> leía MMGBSA_PLATFORM y por defecto elegía CUDA en
     precisión mixta. Se minimizaba en una plataforma y se medía la
     energía en otra. Además get_platform() no aplicaba opciones de
     determinismo: en CPU creaba Contexts con 20 hilos (resultado no
     reproducible y saturación con N workers). Ahora TODA plataforma pasa
     por configurar_plataforma().

  5) RECEPTOR CONGELADO: con --fijo el receptor no se vuelve a preparar
     (PDBFixer colocaba los hidrógenos de forma no determinista). Lo
     prepara una vez preparar_receptores.py con criterio de CAJA de
     docking, y queda como artefacto con md5 registrado.

Ajustes por entorno:
  MMGBSA_MAXITER=8000    iteraciones máximas de minimización
  MMGBSA_TOL=1.0         tolerancia en kJ/mol/nm
  MMGBSA_FINAL_PLATFORM=CPU|CUDA|Reference|OpenCL
  MMGBSA_REPEATS=1       minimizaciones repetidas (se queda con la mejor)
  MMGBSA_CPU_THREADS=1   hilos de la plataforma CPU (1 = determinista)
  MMGBSA_CUDA_DOUBLE=1   fuerza precisión doble en CUDA

Uso:
  python mmgbsa_openff_gb.py --pose <pose.pdbqt> --receptor <rec.pdbqt> \
      --smiles "..." --out json
  python mmgbsa_openff_gb.py --pose <pose.pdbqt> \
      --receptor receptores_fijos/SOD1_fijo.pdb --fijo --smiles "..." --out json
"""
import argparse
import json
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M  # noqa: E402  (helpers de parseo y extracción)

TOLERANCE = float(os.environ.get("MMGBSA_TOL", "1.0"))    # kJ/mol/nm
MAX_ITER = int(os.environ.get("MMGBSA_MAXITER", "8000"))
REPEATS = int(os.environ.get("MMGBSA_REPEATS", "1"))
DEFAULT_PLATFORM = os.environ.get("MMGBSA_FINAL_PLATFORM", "CPU")

# Alinear la plataforma que usa mmgbsa_openff.potential_energy() con la de
# minimización (ese módulo lee MMGBSA_PLATFORM).
if os.environ.get("MMGBSA_PLATFORM", "").strip().upper() != DEFAULT_PLATFORM.upper():
    os.environ["MMGBSA_PLATFORM"] = DEFAULT_PLATFORM

SOLVENT_DIELECTRIC = 78.5
SOLUTE_DIELECTRIC = 1.0
SA_MODEL = "ACE"


def configurar_plataforma(plat):
    """Opciones de determinismo por plataforma (idempotente)."""
    name = plat.getName().upper()
    if name == "CUDA":
        # Mismo resultado entre corridas a coordenadas fijas (medido).
        plat.setPropertyDefaultValue("DeterministicForces", "true")
        if os.environ.get("MMGBSA_CUDA_DOUBLE", "0") == "1":
            plat.setPropertyDefaultValue("Precision", "double")
        if os.environ.get("MMGBSA_CUDA_DEVICE"):
            plat.setPropertyDefaultValue("DeviceIndex",
                                         os.environ["MMGBSA_CUDA_DEVICE"])
    elif name == "CPU":
        # Con 1 hilo las sumas son secuenciales y el resultado reproducible.
        plat.setPropertyDefaultValue("Threads",
                                     os.environ.get("MMGBSA_CPU_THREADS", "1"))
    return plat


# mmgbsa_openff.get_platform() se usa en la minimización Y en cada medición
# de energía, y no aplica opciones de determinismo. Lo envolvemos para que
# toda plataforma quede configurada igual (especialmente Threads=1 en CPU).
_get_platform_base = M.get_platform


def _get_platform_alineada():
    return configurar_plataforma(_get_platform_base())


M.get_platform = _get_platform_alineada


def add_gb_obc2(system, topology, solvent=SOLVENT_DIELECTRIC,
                solute=SOLUTE_DIELECTRIC, sa=SA_MODEL):
    """Añade el término GBSA/OBC2 al System (igual que implicit/obc2.xml)."""
    from openmm import NonbondedForce, CustomNonbondedForce
    from openmm.app.internal.customgbforces import GBSAOBC2Force

    nbs = [f for f in system.getForces() if isinstance(f, NonbondedForce)]
    if len(nbs) != 1:
        raise RuntimeError("se requiere exactamente 1 NonbondedForce (%d)"
                           % len(nbs))
    nb = nbs[0]
    # obc2.xml hace lo mismo: el apantallamiento lo aporta el término GB.
    nb.setReactionFieldDielectric(1.0)

    force = GBSAOBC2Force(solventDielectric=solvent,
                          soluteDielectric=solute, SA=sa)
    # Radios mbondi2 y parámetro de apantallamiento por ELEMENTO, de modo
    # que sirve igual para la proteína y para el ligando.
    params = GBSAOBC2Force.getStandardParameters(topology)
    for i, p in enumerate(params):
        q, _, _ = nb.getParticleParameters(i)
        force.addParticle([q, p[0], p[1]])
    force.finalize()
    force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
    system.addForce(force)
    return force


def prepare_ligand_pose(smiles, pose_path, forcefield):
    """Ligando con la GEOMETRÍA DE LA POSE (determinista)."""
    from openmm import unit
    from rdkit import Chem
    from rdkit.Chem import rdFMCS
    from rdkit.Geometry import Point3D
    from openff.toolkit import Molecule

    docked = M.pose_to_rdmol(pose_path)
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        raise RuntimeError("SMILES invalido: %s" % smiles)
    if rdmol.GetNumAtoms() != docked.GetNumAtoms():
        raise RuntimeError("n atomos no coincide: SMILES %d vs pose %d"
                           % (rdmol.GetNumAtoms(), docked.GetNumAtoms()))
    mcs = rdFMCS.FindMCS([rdmol, docked],
                         bondCompare=rdFMCS.BondCompare.CompareAny,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         timeout=15)
    if mcs.numAtoms != rdmol.GetNumAtoms():
        raise RuntimeError("MCS incompleto: %d/%d"
                           % (mcs.numAtoms, rdmol.GetNumAtoms()))
    pat = Chem.MolFromSmarts(mcs.smartsString)
    m1 = rdmol.GetSubstructMatch(pat)
    m2 = docked.GetSubstructMatch(pat)
    if not m1 or not m2:
        raise RuntimeError("MCS sin correspondencia")

    conf_dock = docked.GetConformer()
    conf = Chem.Conformer(rdmol.GetNumAtoms())
    for i in range(rdmol.GetNumAtoms()):
        conf.SetAtomPosition(i, Point3D(0.0, 0.0, 0.0))
    for i in range(len(m1)):
        p = conf_dock.GetAtomPosition(m2[i])
        conf.SetAtomPosition(m1[i], Point3D(p.x, p.y, p.z))
    rdmol.RemoveAllConformers()
    rdmol.AddConformer(conf, assignId=True)

    # Solo se generan los H (determinista a partir de la geometría de la pose).
    rdmol_h = Chem.AddHs(rdmol, addCoords=True)
    offmol = Molecule.from_rdkit(rdmol_h, allow_undefined_stereo=True)
    try:
        offmol.assign_partial_charges(partial_charge_method="mmff94")
    except Exception:
        pass
    top_openmm = offmol.to_topology().to_openmm()
    vals = np.asarray(offmol.conformers[0].magnitude, dtype=float)
    return top_openmm, unit.Quantity(vals, unit.angstrom), offmol


def _platform_by_name(name):
    """Plataforma concreta por nombre (sin pasar por MMGBSA_PLATFORM)."""
    from openmm import Platform
    for i in range(Platform.getNumPlatforms()):
        p = Platform.getPlatform(i)
        if p.getName().upper() == name.upper():
            return p
    return None


def minimize_on(system, topology, positions, platform_name,
                max_iterations, tolerance):
    """Minimiza en una plataforma concreta con tolerancia explícita."""
    from openmm import LangevinIntegrator, LocalEnergyMinimizer, unit
    from openmm.app import Simulation

    plat = _platform_by_name(platform_name) or M.get_platform()
    configurar_plataforma(plat)
    integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                    2 * unit.femtoseconds)
    sim = Simulation(topology, system, integrator, plat)
    sim.context.setPositions(positions)
    if not M._finite(sim.context.getState(getPositions=True).getPositions()):
        raise RuntimeError("NaN tras setPositions")
    LocalEnergyMinimizer.minimize(
        sim.context,
        tolerance=tolerance * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=max_iterations)
    out = sim.context.getState(getPositions=True).getPositions()
    if not M._finite(out):
        raise RuntimeError("NaN tras minimizar")
    return out


def minimize_best_of(system, topology, positions, repeat=REPEATS,
                     max_iterations=MAX_ITER, tolerance=TOLERANCE,
                     platform_name=None):
    """Minimiza N veces desde la misma geometría y conserva el mínimo más
    bajo."""
    plat = platform_name or DEFAULT_PLATFORM
    best_pos, best_e = None, None
    for _ in range(max(1, repeat)):
        pos = minimize_on(system, topology, positions, plat,
                          max_iterations=max_iterations,
                          tolerance=tolerance)
        e = M.potential_energy(system, topology, pos)
        if best_e is None or e < best_e:
            best_pos, best_e = pos, e
    return best_pos


def rescore_one(pose_path, receptor_pdbqt, smiles):
    from openmm.app import Modeller
    from openff.toolkit import ForceField as OFFForceField, Topology as OFFTopology

    if not os.path.exists(M.FF_PATH):
        raise RuntimeError("no se encuentra openff-2.1.0.offxml en %s" % M.FF_PATH)
    ff = OFFForceField(M.FF_PATH)

    lig_top, lig_pos, offmol = prepare_ligand_pose(smiles, pose_path, ff)
    lig_off = offmol.to_topology()
    lig_center = np.array([[a["x"], a["y"], a["z"]]
                           for a in M.parse_pdbqt_atoms(pose_path)]).mean(axis=0)

    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    try:
        M.extract_receptor_pdb(receptor_pdbqt, lig_center, tmp)
        rec_top, rec_pos = M.prepare_receptor(tmp)
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass

    n_rec = rec_top.getNumAtoms()
    rec_off, rec_mols = M.off_topology_from_pdb(rec_top, rec_pos)

    complex_mod = Modeller(rec_top, rec_pos)
    complex_mod.add(lig_top, lig_pos)
    complex_off = OFFTopology()
    for m in rec_mols:
        complex_off.add_molecule(m)
    complex_off.add_molecule(offmol)

    complex_sys = M.build_system(complex_off, ff)
    add_gb_obc2(complex_sys, complex_mod.topology)
    min_pos = minimize_best_of(complex_sys, complex_mod.topology,
                               complex_mod.positions)
    e_complex = M.potential_energy(complex_sys, complex_mod.topology, min_pos)

    rec_sys = M.build_system(rec_off, ff)
    add_gb_obc2(rec_sys, rec_top)
    e_rec = M.potential_energy(rec_sys, rec_top, min_pos[:n_rec])

    lig_sys = M.build_system(lig_off, ff)
    add_gb_obc2(lig_sys, lig_top)
    e_lig = M.potential_energy(lig_sys, lig_top, min_pos[n_rec:])

    dg = e_complex - e_rec - e_lig
    return {"mmgbsa_dG": round(dg, 2), "e_complex": round(e_complex, 2),
            "e_receptor": round(e_rec, 2), "e_ligand": round(e_lig, 2),
            "n_rec": n_rec}


def preparar_receptor_fijo(receptor_pdbqt, out_pdb, cutoff=12.0,
                           lig_center=None):
    """Prepara el receptor UNA vez y lo guarda en disco (receptor congelado)."""
    from openmm.app import PDBFile
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    try:
        M.extract_receptor_pdb(receptor_pdbqt, lig_center, tmp, cutoff=cutoff)
        rec_top, rec_pos = M.prepare_receptor(tmp)
        with open(out_pdb, "w") as f:
            PDBFile.writeFile(rec_top, rec_pos, f)
        return rec_top.getNumAtoms()
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass


def rescore_one_fijo(pose_path, rec_pdb_preparado, smiles):
    """Rescoring con receptor YA preparado y congelado."""
    from openmm.app import Modeller, PDBFile
    from openff.toolkit import ForceField as OFFForceField, Topology as OFFTopology

    ff = OFFForceField(M.FF_PATH)
    lig_top, lig_pos, offmol = prepare_ligand_pose(smiles, pose_path, ff)
    lig_off = offmol.to_topology()

    rec_pdb = PDBFile(rec_pdb_preparado)
    rec_top, rec_pos = rec_pdb.topology, rec_pdb.positions
    n_rec = rec_top.getNumAtoms()
    rec_off, rec_mols = M.off_topology_from_pdb(rec_top, rec_pos)

    complex_mod = Modeller(rec_top, rec_pos)
    complex_mod.add(lig_top, lig_pos)
    complex_off = OFFTopology()
    for m in rec_mols:
        complex_off.add_molecule(m)
    complex_off.add_molecule(offmol)

    complex_sys = M.build_system(complex_off, ff)
    add_gb_obc2(complex_sys, complex_mod.topology)
    min_pos = minimize_best_of(complex_sys, complex_mod.topology,
                               complex_mod.positions)
    e_complex = M.potential_energy(complex_sys, complex_mod.topology, min_pos)

    rec_sys = M.build_system(rec_off, ff)
    add_gb_obc2(rec_sys, rec_top)
    e_rec = M.potential_energy(rec_sys, rec_top, min_pos[:n_rec])

    lig_sys = M.build_system(lig_off, ff)
    add_gb_obc2(lig_sys, lig_top)
    e_lig = M.potential_energy(lig_sys, lig_top, min_pos[n_rec:])

    dg = e_complex - e_rec - e_lig
    return {"mmgbsa_dG": round(dg, 2), "e_complex": round(e_complex, 2),
            "e_receptor": round(e_rec, 2), "e_ligand": round(e_lig, 2),
            "n_rec": n_rec, "plataforma": DEFAULT_PLATFORM}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pose")
    ap.add_argument("--receptor",
                    help="PDBQT original, o PDB ya preparado si se usa --fijo")
    ap.add_argument("--smiles")
    ap.add_argument("--fijo", action="store_true",
                    help="el receptor ya está preparado (no re-fixear)")
    ap.add_argument("--preparar-receptor", dest="prep_rec", default=None,
                    help="prepara la bolsa una vez y la guarda (requiere "
                         "--pose y --out)")
    ap.add_argument("--cutoff", type=float, default=12.0)
    ap.add_argument("--out", default="json")
    a = ap.parse_args()
    try:
        if a.prep_rec:
            center = np.array([[x["x"], x["y"], x["z"]]
                               for x in M.parse_pdbqt_atoms(a.pose)]).mean(axis=0)
            n = preparar_receptor_fijo(a.prep_rec, a.out, a.cutoff, center)
            print(json.dumps({"receptor_fijo": a.out, "atomos": n}))
            return
        if a.fijo:
            r = rescore_one_fijo(a.pose, a.receptor, a.smiles)
        else:
            r = rescore_one(a.pose, a.receptor, a.smiles)
    except Exception as e:
        r = {"mmgbsa_dG": None, "error": str(e)[:250]}
    print(json.dumps(r))


if __name__ == "__main__":
    main()
