# -*- coding: utf-8 -*-
"""
MM-GBSA con OpenFF (SMIRNOFF 2.1.0) — MASIVE-ALS (PC Windows, sin ambertools).

Sustituye a rescoring_mmgbsa_robusto.py en el PC local:
- Parametriza ligando y receptor con OpenFF 2.1.0 (no necesita GAFF/antechamber).
- Lee poses PDBQT directamente (sin obabel): parseo manual ATOM/HETATM.
- Receptor: extrae el frame del PDBQT de docking y lo repara con pdbfixer.
- MM-GBSA: dG = E(complejo minimizado) - E(receptor) - E(ligando),
  con minimizacion local (L-BFGS via LocalEnergyMinimizer de OpenMM).

Uso:
  python mmgbsa_openff.py --pose <pose.pdbqt> --receptor <receptor.pdbqt>
        --smiles <SMILES> [--out json]

Plataforma: usa la GPU (CUDA) si esta disponible; si no, CPU.
Override: variable de entorno MMGBSA_PLATFORM (Reference|CPU|CUDA|OpenCL).
"""
import argparse
import json
import os
import sys
import tempfile

import numpy as np

from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS
from openff.toolkit import Molecule, ForceField
from openff.units import unit as off_unit


AD_TYPE2ELEM = {
    "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N", "NX": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "P": "P",
    "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "H": "H", "HD": "H", "HS": "H",
}
METAL_ELEMS = {"CU", "ZN", "CA", "FE", "MG", "MN", "CO", "NI", "CD", "NA", "K"}
# Ruta absoluta del offxml (mas robusta: la version yanked de openff-toolkit
# no registra la carpeta data/forcefield por defecto en el registro de nombres).
def _find_ff():
    """Localiza openff-2.1.0.offxml de forma portable.

    Orden: MMGBSA_FF_PATH > paquete openforcefields instalado >
    openff-toolkit instalado > copia dentro del entorno de Windows.
    """
    import site
    cands = [os.environ.get("MMGBSA_FF_PATH", "")]
    roots = []
    try:
        roots += list(site.getsitepackages())
    except Exception:
        pass
    try:
        import openforcefields
        roots.append(os.path.dirname(os.path.abspath(openforcefields.__file__)))
    except Exception:
        pass
    for sp in roots:
        cands.append(os.path.join(sp, "openforcefields", "offxml",
                                  "openff-2.1.0.offxml"))
        cands.append(os.path.join(sp, "openff", "toolkit", "data",
                                  "forcefield", "openff-2.1.0.offxml"))
    cands.append(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                              "..", "..", "rescoring_env", "Lib",
                              "site-packages", "openff", "toolkit", "data",
                              "forcefield", "openff-2.1.0.offxml"))
    for c in cands:
        if c and os.path.exists(c):
            return c
    return "openff-2.1.0"


FF_PATH = _find_ff()
FF_NAME = "openff-2.1.0"


def parse_pdbqt_atoms(pdbqt_path, primer_model=True):
    """Devuelve lista de dicts {name,resname,resnum,x,y,z,elem} de ATOM/HETATM.

    Si `primer_model` es True (por defecto), solo lee el primer bloque MODEL
    (las poses de Vina-GPU tienen num_modes=3).
    """
    atoms = []
    en_model = False
    for line in open(pdbqt_path, encoding="utf-8", errors="replace"):
        if line.startswith("MODEL"):
            if atoms and primer_model:
                break  # ya leido el primer modelo
            en_model = True
            continue
        if line.startswith("ENDMDL"):
            en_model = False
            if primer_model and atoms:
                break
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            name = line[12:16].strip()
            resname = line[17:20].strip()
            resnum = line[22:26].strip()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            adtype = line[77:80].strip().upper()
            elem = AD_TYPE2ELEM.get(adtype, adtype)
        except (ValueError, IndexError):
            continue
        if elem in ("H",) or elem in METAL_ELEMS:
            continue
        atoms.append(dict(name=name, resname=resname, resnum=resnum,
                          x=x, y=y, z=z, elem=elem))
    return atoms


STANDARD_AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
               "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
               "THR", "TRP", "TYR", "VAL"}


def extract_receptor_pdb(receptor_pdbqt, lig_center, out_pdb, cutoff=12.0):
    """Escribe el bolsillo del receptor como PDB simple a partir del PDBQT.

    Los PDBQT de receptores vienen sucios tras el preparado de Vina-GPU:
    - FUS tenia restos de RNA (G, U) que OpenFF no puede leer.
    - SOD1 trae aguas (HOH) y metales (CU, ZN, CA).
    - SOD1 ademas CONCATENA ~18 copias de la cadena con numeracion 1..153
      repetida; sin tratarlas, OpenMM fusiona los residuos con el mismo
      numero en un residuo gigante y OpenFF lo rechaza.

    Limpieza:
    1) Solo aminoacidos estandar (fuera HOH, metales, RNA y raros).
    2) Por cada numero de residuo se conserva la COPIA mas cercana al
       centro de la pose (las copias lejanas se cancelarian en dG).
    3) Se conservan los residuos con algun atomo a menos de `cutoff` A
       del centro de la pose (bolsillo). Si no quedara ninguno, los 60
       mas cercanos.
    4) Reenumeracion secuencial y TER (cadena nueva) entre tramos
       separados, para que OpenMM no enlace residuos que no tocan.

    MM-GBSA local: el resto de la proteina cancela en dG = E_c - E_r y no
    aporta (interacciones de corto alcance), pero multiplicaria el costo
    computacional por el numero de copias.
    """
    atoms = parse_pdbqt_atoms(receptor_pdbqt)
    if not atoms:
        raise RuntimeError("receptor sin atomos: %s" % receptor_pdbqt)

    # 1) solo aminoacidos estandar
    atoms = [a for a in atoms if a["resname"] in STANDARD_AA]
    if not atoms:
        raise RuntimeError("receptor sin aminoacidos estandar: %s" % receptor_pdbqt)

    center = None if lig_center is None else np.asarray(lig_center, dtype=float)

    # 2) agrupar residuos consecutivos y marcar copias (reinicios numericos)
    residues = []  # dicts: resnum, resname, copy, atoms
    copy_id, prev_rn = 0, None
    for a in atoms:
        try:
            rn = int(a["resnum"])
        except (TypeError, ValueError):
            continue
        if prev_rn is not None and rn < prev_rn:
            copy_id += 1
        if not residues or residues[-1]["resnum"] != rn or residues[-1]["copy"] != copy_id:
            residues.append({"resnum": rn, "resname": a["resname"],
                             "copy": copy_id, "atoms": []})
        residues[-1]["atoms"].append(a)
        prev_rn = rn
    if not residues:
        raise RuntimeError("receptor sin residuos validos: %s" % receptor_pdbqt)

    def res_dist(r):
        if center is None:
            return 0.0
        pts = np.array([[x["x"], x["y"], x["z"]] for x in r["atoms"]])
        return float(np.linalg.norm(pts - center, axis=1).min())

    # 3) por resnum, la copia mas cercana al ligando
    best = {}
    for r in residues:
        d = res_dist(r)
        k = r["resnum"]
        if k not in best or d < best[k][0]:
            best[k] = (d, r)
    kept = [(d, r) for k, (d, r) in best.items()
            if center is None or d <= cutoff]
    if not kept:
        kept = sorted(best.values(), key=lambda t: t[0])[:60]
    kept.sort(key=lambda t: (t[1]["copy"], t[1]["resnum"]))  # copias contiguas

    # 4) escribir PDB: reenumerado, cadena nueva (TER) entre tramos
    lines, serial, chain_idx, res_serial = [], 0, 0, 0
    prev = None
    for d, r in kept:
        nueva_cadena = (prev is not None and
                        (r["copy"] != prev["copy"] or
                         r["resnum"] != prev["resnum"] + 1))
        if nueva_cadena:
            lines.append("TER")
            chain_idx += 1
            res_serial = 0
        res_serial += 1
        for a in r["atoms"]:
            serial += 1
            lines.append(
                "ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
                % (serial, a["name"][:4], a["resname"][:3],
                   chr(ord("A") + chain_idx % 26), res_serial,
                   a["x"], a["y"], a["z"], a["elem"][:2]))
        prev = r
    lines.append("TER")
    lines.append("END")
    with open(out_pdb, "w") as f:
        f.write("\n".join(lines) + "\n")
    return out_pdb


def pose_to_rdmol(pose_path):
    """PDBQT -> RDKit mol con conformador (sin obabel, sin Hs).

    Añade enlaces por distancia (ConnectTheDots) para que el MCS con
    ringMatchesRingOnly funcione (la pose PDBQT no trae enlaces).
    """
    atoms = parse_pdbqt_atoms(pose_path)
    if not atoms:
        raise RuntimeError("pose vacia: %s" % pose_path)
    rw = Chem.RWMol()
    conf = Chem.Conformer(len(atoms))
    for i, a in enumerate(atoms):
        idx = rw.AddAtom(Chem.Atom(a["elem"]))
        conf.SetAtomPosition(idx, (a["x"], a["y"], a["z"]))
    rw.AddConformer(conf)
    # RDKit >= 2025 eliminó ConnectTheDots/DetermineBonds: enlaces por
    # distancia (misma heurística: suma de radios covalentes + 0.45 A).
    pt = Chem.GetPeriodicTable().GetRcovalent
    for i in range(rw.GetNumAtoms()):
        ai = rw.GetAtomWithIdx(i)
        for j in range(i + 1, rw.GetNumAtoms()):
            aj = rw.GetAtomWithIdx(j)
            if ai.GetAtomicNum() == 0 or aj.GetAtomicNum() == 0:
                continue
            d = (conf.GetAtomPosition(i) - conf.GetAtomPosition(j)).Length()
            r = pt(ai.GetAtomicNum()) + pt(aj.GetAtomicNum()) + 0.45
            if d < r:
                rw.AddBond(i, j, Chem.BondType.SINGLE)
    return rw.GetMol()


def prepare_ligand(smiles, pose_path, forcefield):
    """Devuelve (topology_openmm, positions_openmm) del ligando en la pose."""
    from openmm import unit
    docked = pose_to_rdmol(pose_path)
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        raise RuntimeError("SMILES invalido: %s" % smiles)
    if rdmol.GetNumAtoms() != docked.GetNumAtoms():
        raise RuntimeError(
            "n atomos no coincide: SMILES %d vs pose %d (%s)"
            % (rdmol.GetNumAtoms(), docked.GetNumAtoms(), pose_path))
    mcs = rdFMCS.FindMCS([rdmol, docked],
                         bondCompare=rdFMCS.BondCompare.CompareAny,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         timeout=15)
    if mcs.numAtoms != rdmol.GetNumAtoms():
        raise RuntimeError("MCS incompleto: %d/%d" % (mcs.numAtoms, rdmol.GetNumAtoms()))
    pat = Chem.MolFromSmarts(mcs.smartsString)
    m1 = rdmol.GetSubstructMatch(pat)
    m2 = docked.GetSubstructMatch(pat)
    if not m1 or not m2:
        raise RuntimeError("MCS sin correspondencia")
    conf = docked.GetConformer()
    coordmap = {}
    for i in range(len(m1)):
        p = conf.GetAtomPosition(m2[i])
        coordmap[m1[i]] = Chem.rdGeometry.Point3D(p.x, p.y, p.z)
    rdmol_h = Chem.AddHs(rdmol)
    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    params.maxIterations = 500
    params.SetCoordMap(coordmap)
    ok = -1
    for seed in (42, 2026, 777, 1234, 9999, 31337):
        params.randomSeed = seed
        ok = AllChem.EmbedMolecule(rdmol_h, params)
        if ok == 0:
            break
    if ok != 0:
        raise RuntimeError("embed fallo (6 semillas) para %s" % pose_path)

    offmol = Molecule.from_rdkit(rdmol_h, allow_undefined_stereo=True)
    # Cargas: sin ambertools en Windows, usamos mmff94 (RDKit).
    try:
        offmol.assign_partial_charges(partial_charge_method="mmff94")
    except Exception:
        pass
    top_openmm = offmol.to_topology().to_openmm()
    # El conformer de OpenFF es pint.Quantity; OpenMM exige unit.Quantity
    # propio (si no, Modeller.add falla por items sin unidades).
    conf = offmol.conformers[0]
    vals = np.asarray(conf.magnitude, dtype=float)
    pos = unit.Quantity(vals, unit.angstrom)
    return top_openmm, pos, offmol


def prepare_receptor(pdb_path):
    from pdbfixer import PDBFixer
    fixer = PDBFixer(filename=pdb_path)
    fixer.missingResidues = []
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer.topology, fixer.positions


def off_topology_from_pdb(top_openmm, positions):
    """openmm.app.Topology + positions -> openff.toolkit.Topology (API 0.18).

    Escribe el sistema a PDB temporal y usa Topology.from_pdb (soporta
    multiples fragmentos de bolsillo, que from_polymer_pdb rechaza).
    Devuelve (topologia_openff, lista_de_moleculas).
    """
    from openff.toolkit import Topology as OFFTopology
    from openmm.app import PDBFile
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    try:
        with open(tmp, "w") as f:
            PDBFile.writeFile(top_openmm, positions, f)
        off = OFFTopology.from_pdb(tmp)
        return off, list(off.molecules)
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass


def build_system(topology_off, forcefield):
    """topology_off: openff.toolkit.Topology -> openmm.System."""
    from openmm import NonbondedForce
    fn = getattr(forcefield, "create_openmm_system",
                 getattr(forcefield, "createOpenMMSystem", None))
    if fn is None:
        raise RuntimeError("ForceField sin create_openmm_system")
    # Sin ambertools, OpenFF no puede asignar cargas AM1BCC: asignamos
    # mmff94 (RDKit) a las moleculas que no tengan cargas parciales y las
    # pasamos via charge_from_molecules (si no, create_openmm_system
    # reintenta AM1BCC y falla).
    chg, seen_sig = [], set()
    for m_ in topology_off.molecules:
        try:
            if m_.partial_charges is None:
                m_.assign_partial_charges(partial_charge_method="mmff94")
            # firma isomorfica (SMILES canonico): moléculas isomorfas —
            # p.ej. las dos cadenas identicas de un homodimero — comparten
            # firma y charge_from_molecules exige unicidad isomorfica;
            # mmff94 solo depende del grafo, no de la geometria
            try:
                sig = ("smi", m_.to_smiles())
            except Exception:
                sig = ("graph",
                       tuple(a.atomic_number for a in m_.atoms),
                       tuple(sorted((b.atom1.molecule_atom_index,
                                     b.atom2.molecule_atom_index)
                                    for b in m_.bonds)))
            if sig not in seen_sig:
                seen_sig.add(sig)
                chg.append(m_)
        except Exception:
            pass
    system = fn(topology_off, charge_from_molecules=chg)
    # OpenFF 0.18 no expone nonbondedMethod: forzamos NoCutoff sobre el
    # NonbondedForce (MM-GBSA sin correccion de corte).
    n_conv = 0
    for i in range(system.getNumForces()):
        f = system.getForce(i)
        if isinstance(f, NonbondedForce):
            f.setNonbondedMethod(NonbondedForce.NoCutoff)
            n_conv += 1
    if n_conv == 0:
        raise RuntimeError("sin NonbondedForce en el sistema")
    return system


def _finite(positions):
    from openmm import unit as u
    try:
        a = positions.value_in_unit(u.angstrom)
    except Exception:
        a = np.asarray(positions)
    return np.isfinite(np.asarray(a, dtype=float)).all()


def get_platform():
    """Plataforma OpenMM: CUDA (GPU) si existe, si no CPU. Override con
    la variable de entorno MMGBSA_PLATFORM (p.ej. MMGBSA_PLATFORM=CPU)."""
    import openmm
    from openmm import Platform
    wanted = os.environ.get("MMGBSA_PLATFORM", "").strip().upper()
    names = [Platform.getPlatform(i).getName().upper()
             for i in range(Platform.getNumPlatforms())]
    if wanted:
        if wanted not in names:
            raise RuntimeError("MMGBSA_PLATFORM=%s no disponible (%s)"
                               % (wanted, ",".join(names)))
        idx = names.index(wanted)
    elif "CUDA" in names:
        idx = names.index("CUDA")
    elif "OPENCL" in names:
        idx = names.index("OPENCL")
    elif "CPU" in names:
        idx = names.index("CPU")
    else:
        idx = 0
    plat = Platform.getPlatform(idx)
    print("[mmgbsa] plataforma OpenMM: %s" % plat.getName(), file=sys.stderr,
          flush=True)
    return plat


def minimize(system, topology, positions, max_iterations=1000):
    from openmm import LangevinIntegrator, LocalEnergyMinimizer, unit
    from openmm.app import Simulation
    integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                    2 * unit.femtoseconds)
    sim = Simulation(topology, system, integrator, get_platform())
    sim.context.setPositions(positions)
    if not _finite(sim.context.getState(getPositions=True).getPositions()):
        raise RuntimeError("NaN tras setPositions")
    LocalEnergyMinimizer.minimize(sim.context, maxIterations=max_iterations)
    out = sim.context.getState(getPositions=True).getPositions()
    if not _finite(out):
        raise RuntimeError("NaN tras minimizar")
    return out


def potential_energy(system, topology, positions):
    from openmm import LangevinIntegrator, unit
    from openmm.app import Simulation
    integrator = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                    2 * unit.femtoseconds)
    sim = Simulation(topology, system, integrator, get_platform())
    sim.context.setPositions(positions)
    state = sim.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def rescore_one(pose_path, receptor_pdbqt, smiles):
    from openmm.app import Modeller
    from openff.toolkit import ForceField as OFFForceField
    if not os.path.exists(FF_PATH):
        raise RuntimeError("no se encuentra openff-2.1.0.offxml en %s" % FF_PATH)
    ff = OFFForceField(FF_PATH)
    lig_top, lig_pos, offmol = prepare_ligand(smiles, pose_path, ff)
    lig_off = offmol.to_topology()
    lig_center = np.array([[a["x"], a["y"], a["z"]]
                           for a in parse_pdbqt_atoms(pose_path)]).mean(axis=0)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    try:
        extract_receptor_pdb(receptor_pdbqt, lig_center, tmp)
        rec_top, rec_pos = prepare_receptor(tmp)
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass
    n_rec = rec_top.getNumAtoms()
    rec_off, rec_mols = off_topology_from_pdb(rec_top, rec_pos)
    complex_mod = Modeller(rec_top, rec_pos)
    complex_mod.add(lig_top, lig_pos)
    # Complejo: Topology OpenFF = [receptor(fragmentos), ligando] (mismo
    # orden que complex_mod: receptor primero, ligando despues).
    from openff.toolkit import Topology as OFFTopology
    complex_off = OFFTopology()
    for m in rec_mols:
        complex_off.add_molecule(m)
    complex_off.add_molecule(offmol)
    complex_sys = build_system(complex_off, ff)
    min_pos = minimize(complex_sys, complex_mod.topology, complex_mod.positions)
    e_complex = potential_energy(complex_sys, complex_mod.topology, min_pos)
    rec_sys = build_system(rec_off, ff)
    e_rec = potential_energy(rec_sys, rec_top, min_pos[:n_rec])
    lig_sys = build_system(lig_off, ff)
    e_lig = potential_energy(lig_sys, lig_top, min_pos[n_rec:])
    dg = e_complex - e_rec - e_lig
    return {"mmgbsa_dG": round(dg, 2), "e_complex": round(e_complex, 2),
            "e_receptor": round(e_rec, 2), "e_ligand": round(e_lig, 2)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pose", required=True)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--out", default="json")
    args = ap.parse_args()
    try:
        r = rescore_one(args.pose, args.receptor, args.smiles)
        r["error"] = None
    except Exception as e:
        r = {"mmgbsa_dG": None, "error": str(e)[:250]}
    if args.out == "json":
        print(json.dumps(r))
    else:
        with open(args.out, "w") as f:
            json.dump(r, f)


if __name__ == "__main__":
    main()