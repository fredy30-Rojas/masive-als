#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rescoring MM-GBSA v7 (AM1-BCC real) — MASIVE-ALS.

v7 (22 sep 2026, noche): igual que v6 pero con cargas AM1-BCC de verdad.

v6 (la de la tabla del informe) pone las cargas del ligando a cero para saltarse
el AM1-BCC (es lento: ~1 min por ligando en CPU). Eso convierte el dG en algo
muy parecido a una energia de van der Waals mas superficie, y su pendiente con
el tamaño (-0.498 kcal/mol por atomo pesado) podria ser consecuencia del atajo.

v7 responde: si el sesgo era del atajo, con cargas reales baja. Si no baja, es
de la medida. Todo lo demas es IDENTICO a v6 (receptor semillado con
SEMILLA_RECEPTOR=20260922, un hilo de OpenMM, GAFF 2.11, OBC2, una pose por
ligando, minimizacion de 200 pasos).

Uso (en Oracle, env conda mmgbsa):
  export OPENMM_CPU_THREADS=1 OMP_NUM_THREADS=1
  python rescoring_mmgbsa_bcc.py --check
  python rescoring_mmgbsa_bcc.py --candidates estratos/candidatos.csv \
      --receptores receptores_estratos.json --row 1
  python rescoring_mmgbsa_bcc.py --candidates estratos/candidatos.csv \
      --receptores receptores_estratos.json --out rescoring_bcc.csv
"""
import argparse
import json
import os
import random
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

REQ = ["openmm", "pdbfixer", "rdkit", "pandas", "numpy",
       "openmmforcefields", "openff.toolkit"]

AD_TYPE2ELEM = {
    "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N", "NX": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "P": "P",
    "F": "F", "CL": "CL", "BR": "BR", "I": "I",
    "H": "H", "HD": "H", "HS": "H",
}
METAL_ELEMS = {"CU", "ZN", "CA", "FE", "MG", "MN", "CO", "NI", "CD", "NA", "K"}

# semilla de la preparacion del receptor: igual que v6
SEMILLA_RECEPTOR = 20260922


def sh(cmd):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"cmd falló ({r.returncode}): {cmd}\n{r.stderr[:400]}")
    return r.stdout


def check_dependencies():
    missing = []
    for m in REQ:
        try:
            __import__(m.split(".")[0])
        except Exception:
            missing.append(m)
    for exe in ["obabel", "antechamber"]:
        if subprocess.run(f"which {exe}", shell=True,
                          capture_output=True).returncode != 0:
            missing.append(exe)
    return missing


def parse_pdbqt(pdbqt_path):
    chains = {}
    for line in open(pdbqt_path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21:22].strip() or " "
        resnum = line[22:26].strip()
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        adtype = line[77:80].strip().upper()
        elem = AD_TYPE2ELEM.get(adtype, adtype)
        if elem in ("H",) or elem in METAL_ELEMS:
            continue
        chains.setdefault(chain, []).append(
            dict(name=name, resname=resname, chain=chain, resnum=resnum,
                 x=x, y=y, z=z, elem=elem))
    return chains


def extract_receptor_pdb(pdbqt_path, lig_center, out_pdb):
    chains = parse_pdbqt(pdbqt_path)
    if not chains:
        raise RuntimeError("receptor sin átomos: " + pdbqt_path)
    keep_chains = list(chains.keys())
    if len(chains) > 4:
        best, best_d = None, 1e18
        for c, atoms in chains.items():
            d = min(((a["x"] - lig_center[0]) ** 2 +
                     (a["y"] - lig_center[1]) ** 2 +
                     (a["z"] - lig_center[2]) ** 2) ** 0.5 for a in atoms)
            if d < best_d:
                best_d, best = d, c
        keep_chains = [best]

    lines, serial = [], 0
    for c in keep_chains:
        for a in chains[c]:
            serial += 1
            lines.append(
                f"ATOM  {serial:5d} {a['name']:>4s} {a['resname']:>3s} "
                f"{a['chain']:1s}{int(a['resnum']):4d}    "
                f"{a['x']:8.3f}{a['y']:8.3f}{a['z']:8.3f}"
                f"  1.00  0.00          {a['elem']:>2s}")
        lines.append("TER")
    lines.append("END")
    with open(out_pdb, "w") as f:
        f.write("\n".join(lines) + "\n")
    return out_pdb


def _recortar_pdb_canonical(pdb_path, lig_center, out_pdb, radio=10.0):
    """Identica a v6: recorte canonico; cae el residuo que queda solo en su cadena."""
    import math
    keep = set()
    rec_atoms = []
    for line in open(pdb_path, errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21:22].strip() or " "
        resnum = line[22:26].strip()
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        d = math.sqrt((x - lig_center[0])**2 + (y - lig_center[1])**2 + (z - lig_center[2])**2)
        if d < radio:
            keep.add((chain, resnum))
        rec_atoms.append(line)

    por_cadena = {}
    for ch, num in keep:
        por_cadena.setdefault(ch, set()).add(num)
    sueltas = {ch for ch, nums in por_cadena.items() if len(nums) == 1}
    if sueltas:
        keep = {k for k in keep if k[0] not in sueltas}
    lines = []
    serial = 0
    for line in rec_atoms:
        chain = line[21:22].strip() or " "
        if (chain, line[22:26].strip()) in keep:
            serial += 1
            lines.append(line[:6] + "%5d" % serial + line[11:].rstrip())
    if not lines:
        raise RuntimeError("recorte vacío (radio %.1f A)" % radio)
    with open(out_pdb, "w") as f:
        f.write("\n".join(lines) + "\nEND\n")
    return out_pdb


def prepare_receptor(pdb_path):
    """Identica a v6: PDBFixer con semilla fija + random/numpy sembrados."""
    from pdbfixer import PDBFixer
    import tempfile

    cleaned = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w")
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("HETATM"):
                resname = line[17:20].strip()
                if resname in ("CA", "CAL", "CAC"):
                    continue
            if line.startswith("ATOM"):
                elem = line[76:78].strip() if len(line) > 77 else ""
                if elem == "CA":
                    continue
            cleaned.write(line)
    cleaned.close()

    fixer = PDBFixer(filename=cleaned.name)
    fixer.missingResidues = []
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=SEMILLA_RECEPTOR)
    random.seed(SEMILLA_RECEPTOR)
    np.random.seed(SEMILLA_RECEPTOR)
    fixer.addMissingHydrogens(7.0)
    os.unlink(cleaned.name)
    return fixer.topology, fixer.positions


def parse_pose_pdbqt(pose_path):
    atoms = []
    for line in open(pose_path):
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            atoms.append((name, x, y, z))
    return atoms


def ligand_center(pose_path):
    """Identica a v6: centro de la PRIMERA pose (MODEL 1)."""
    import pose_pdbqt

    pose = pose_pdbqt.leer_pose(pose_path, modelo=1)
    at = list(pose["atoms"].values())
    if not at:
        raise RuntimeError("pose vacía: " + pose_path)
    return np.array([[a[0], a[1], a[2]] for a in at]).mean(axis=0)


def prepare_ligand(smiles, pose_path, forcefield, cargas_bcc=True):
    """v7: lo unico que cambia — AM1-BCC real en vez de cargas a cero."""
    from openmmforcefields.generators import GAFFTemplateGenerator
    from openff.toolkit import Molecule
    from rdkit import Chem
    from rdkit.Chem import AllChem

    import pose_pdbqt

    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        raise RuntimeError(f"SMILES ilegible: {smiles[:60]}")

    pose = pose_pdbqt.leer_pose(pose_path)
    coordmap = pose_pdbqt.coordmap_de_pose(pose, rdmol)

    rdmol_h = Chem.AddHs(rdmol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    params.maxIterations = 500
    params.SetCoordMap(coordmap)
    ret = AllChem.EmbedMolecule(rdmol_h, params)
    if ret != 0:
        p2 = AllChem.ETKDGv3()
        p2.randomSeed = 43
        p2.useRandomCoords = True
        p2.maxIterations = 2000
        ret = AllChem.EmbedMolecule(rdmol_h, p2)
        if ret == 0:
            _llevar_a_la_pose(rdmol_h, coordmap)
    if ret != 0:
        raise RuntimeError(f"embedding fallo (ret={ret}) para {smiles[:40]}")

    offmol = Molecule.from_rdkit(rdmol_h, allow_undefined_stereo=True)
    if cargas_bcc:
        # v7: cargas AM1-BCC reales. Los wrappers por defecto (RDKit/Built-in) no
        # las implementan: hay que pasar por AmberToolsToolkitWrapper, que usa el
        # antechamber/sqm del env (por eso el runner activa el conda mmgbsa).
        from openff.toolkit.utils.toolkits import AmberToolsToolkitWrapper
        offmol.assign_partial_charges(
            partial_charge_method="am1bcc",
            toolkit_registry=AmberToolsToolkitWrapper())
        if offmol.partial_charges is None or len(offmol.partial_charges) == 0:
            raise RuntimeError("AM1-BCC no devolvio cargas")
    else:
        import numpy as np
        from openff.units import unit
        zero_charges = np.zeros(offmol.n_atoms) * unit.elementary_charge
        offmol.partial_charges = zero_charges
    generator = GAFFTemplateGenerator(molecules=[offmol],
                                      forcefield="gaff-2.11")
    forcefield.registerTemplateGenerator(generator.generator)

    from openmm import unit as omm_unit
    top = offmol.to_topology().to_openmm()
    pos = offmol.conformers[0].magnitude * omm_unit.angstrom
    return top, pos


def _llevar_a_la_pose(rdmol_h, coordmap):
    """Identica a v6."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign

    ref = Chem.Mol(rdmol_h)
    conf = ref.GetConformer()
    for j, p in coordmap.items():
        conf.SetAtomPosition(j, p)
    rdMolAlign.AlignMol(rdmol_h, ref, atomMap=[(j, j) for j in coordmap])


def build_system(topology, forcefield):
    from openmm.app import NoCutoff, HBonds
    return forcefield.createSystem(topology, nonbondedMethod=NoCutoff,
                                   constraints=HBonds)


def _finite(positions):
    from openmm import unit as u
    a = positions.value_in_unit(u.angstrom) if hasattr(positions, "value_in_unit") \
        else np.asarray(positions)
    return np.isfinite(np.asarray(a)).all()


def minimize(system, topology, positions, max_iterations=200):
    from openmm import LangevinIntegrator, LocalEnergyMinimizer, unit
    from openmm.app import Simulation

    integrator = LangevinIntegrator(300 * unit.kelvin,
                                    1 / unit.picosecond,
                                    2 * unit.femtoseconds)
    sim = Simulation(topology, system, integrator)
    sim.context.setPositions(positions)
    if not _finite(sim.context.getState(getPositions=True).getPositions()):
        raise RuntimeError("coordenadas NaN tras setPositions")
    LocalEnergyMinimizer.minimize(sim.context, maxIterations=max_iterations)
    out = sim.context.getState(getPositions=True).getPositions()
    if not _finite(out):
        raise RuntimeError("coordenadas NaN tras minimizar")
    return out


def potential_energy(system, topology, positions):
    from openmm import LangevinIntegrator, unit
    from openmm.app import Simulation

    integrator = LangevinIntegrator(300 * unit.kelvin,
                                    1 / unit.picosecond,
                                    2 * unit.femtoseconds)
    sim = Simulation(topology, system, integrator)
    sim.context.setPositions(positions)
    state = sim.context.getState(getEnergy=True)
    return state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)


def rescore_one(receptor_top, receptor_pos, smiles, pose_path, forcefield):
    """v7: pasa el flag de cargas al preparador del ligando."""
    from openmm.app import Modeller

    lig_top, lig_pos = prepare_ligand(smiles, pose_path, forcefield,
                                      cargas_bcc=True)
    n_rec = receptor_top.getNumAtoms()

    complex_mod = Modeller(receptor_top, receptor_pos)
    complex_mod.add(lig_top, lig_pos)

    complex_sys = build_system(complex_mod.topology, forcefield)
    min_pos = minimize(complex_sys, complex_mod.topology,
                       complex_mod.positions)
    e_complex = potential_energy(complex_sys, complex_mod.topology, min_pos)

    rec_sys = build_system(receptor_top, forcefield)
    e_rec = potential_energy(rec_sys, receptor_top, min_pos[:n_rec])

    lig_sys = build_system(lig_top, forcefield)
    e_lig = potential_energy(lig_sys, lig_top, min_pos[n_rec:])

    return e_complex - e_rec - e_lig, e_complex, e_rec, e_lig


def compute_row(r, receptores, forcefield):
    """Identica a v6 salvo el preparador del ligando."""
    tgt, lig, smiles = r["target"], r["ligand"], r["smiles"]
    pose = r["pose_pdbqt"]
    base = {"ligand": lig, "target": tgt, "vina_affinity": r.get("affinity")}
    if tgt not in receptores:
        return {**base, "mmgbsa_dG": None, "error": "target sin receptor"}
    try:
        c = ligand_center(pose)
        tmp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        rpath = receptores[tgt]
        if rpath.lower().endswith(".pdb") and not rpath.lower().endswith(".pdbqt"):
            _recortar_pdb_canonical(rpath, c, tmp_pdb, radio=14.0)
        else:
            extract_receptor_pdb(rpath, c, tmp_pdb)
        rec_top, rec_pos = prepare_receptor(tmp_pdb)
        dg, ec, er, el = rescore_one(rec_top, rec_pos, smiles, pose, forcefield)
        return {**base, "mmgbsa_dG": round(dg, 2),
                "e_complex": round(ec, 2), "e_receptor": round(er, 2),
                "e_ligand": round(el, 2), "error": None}
    except Exception as e:
        return {**base, "mmgbsa_dG": None, "error": str(e)[:200]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--candidates", default="candidatos_42.csv")
    ap.add_argument("--receptores", default="receptores.json")
    ap.add_argument("--out", default="rescoring_bcc.csv")
    ap.add_argument("--row", type=int, default=None)
    args = ap.parse_args()

    missing = check_dependencies()
    if args.check:
        for m in REQ + ["obabel", "antechamber"]:
            print("  [%s] %s" % ("OK" if m not in missing else "FALTA", m))
        sys.exit(0 if not missing else 1)
    if missing:
        print("Faltan: %s" % ", ".join(missing)); sys.exit(1)

    from openmm.app import ForceField
    forcefield = ForceField("amber14-all.xml", "implicit/obc2.xml")

    with open(args.receptores) as f:
        receptores = json.load(f)

    df = pd.read_csv(args.candidates)

    if args.row is not None:
        r = df.iloc[args.row]
        out = compute_row(r, receptores, forcefield)
        print(json.dumps(out))
        sys.exit(0)

    rows = []
    for idx, r in df.iterrows():
        out = compute_row(r, receptores, forcefield)
        rows.append(out)
        tag = f"dG={out['mmgbsa_dG']:.2f}" if out["mmgbsa_dG"] is not None \
            else f"ERROR: {out['error'][:120]}"
        print(f"[{idx+1}/{len(df)}] {out['ligand']} {out['target']} {tag}",
              flush=True)

    out_df = pd.DataFrame(rows).sort_values("mmgbsa_dG", na_position="last")
    out_df.to_csv(args.out, index=False)
    print("Guardado:", args.out)
    print(out_df.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
