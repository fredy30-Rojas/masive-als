#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rescoring MM-GBSA v3 — MASIVE-ALS.

Usa las POSES acopladas por Vina (PDBQT) como geometría inicial del ligando,
recorta el receptor a la cadena de unión (y quita metales/H) y calcula
MM-GBSA de una sola trayectoria con OpenMM + GBSA-OBC2:

    dG_bind = E(complejo_min) - E(receptor) - E(ligando)

Requisitos (Linux/Conda, Python 3.11): openmm, openmmforcefields,
openff-toolkit, rdkit, pdbfixer, pandas, numpy, ambertools (antechamber) y obabel.

Uso:
    python rescoring_mmgbsa.py --check
    python rescoring_mmgbsa.py --candidates candidatos_42.csv \
        --receptores receptores.json --out rescoring_mmgbsa.csv
"""
import argparse
import json
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

REQ = ["openmm", "pdbfixer", "rdkit", "pandas", "numpy",
       "openmmforcefields", "openff.toolkit"]

# AutoDock atom type -> elemento PDB
AD_TYPE2ELEM = {
    "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N", "NX": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "P": "P",
    "F": "F", "CL": "CL", "BR": "BR", "I": "I",
    "H": "H", "HD": "H", "HS": "H",
}
# metales: OpenMM no los parametriza -> se descartan del receptor
METAL_ELEMS = {"CU", "ZN", "CA", "FE", "MG", "MN", "CO", "NI", "CD", "NA", "K"}


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


# ----------------------------------------------------------------------------
# Receptor
# ----------------------------------------------------------------------------
def parse_pdbqt(pdbqt_path):
    """Átomos del receptor PDBQT -> dicts por cadena (sin H ni metales)."""
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
            continue  # fuera H (se reañade con PDBFixer) y metales
        chains.setdefault(chain, []).append(
            dict(name=name, resname=resname, chain=chain, resnum=resnum,
                 x=x, y=y, z=z, elem=elem))
    return chains


def extract_receptor_pdb(pdbqt_path, lig_center, out_pdb):
    """Recorta el receptor PDBQT a la(s) cadena(s) de unión y escribe PDB limpio.

    Si hay > 4 cadenas (p.ej. el ensamblado de 18 cadenas de SOD1), conserva
    solo la cadena cuyo átomo queda más cerca del ligando.
    """
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


def prepare_receptor(pdb_path):
    """PDBFixer sin rellenar residuos ausentes (el recorte no es contiguo)."""
    from pdbfixer import PDBFixer

    fixer = PDBFixer(filename=pdb_path)
    fixer.missingResidues = []   # no rellenamos residuos
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    return fixer.topology, fixer.positions


# ----------------------------------------------------------------------------
# Ligando (pose acoplada + SMILES autoritativo)
# ----------------------------------------------------------------------------
def parse_pose_pdbqt(pose_path):
    """Coordenadas de los átomos del mejor modelo (MODEL 1)."""
    atoms = []
    for line in open(pose_path):
        if line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip()
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            atoms.append((name, x, y, z))
    return atoms


def ligand_center(pose_path):
    at = parse_pose_pdbqt(pose_path)
    if not at:
        raise RuntimeError("pose vacía: " + pose_path)
    c = np.array([[x, y, z] for (_, x, y, z, _) in at]) if at and len(at[0]) == 5 \
        else np.array([[a[1], a[2], a[3]] for a in at])
    return c.mean(axis=0)


def prepare_ligand(smiles, pose_path, forcefield):
    """SMILES + pose PDBQT -> (topology, positions) con H y cargas GAFF/AM1-BCC."""
    from openmmforcefields.generators import GAFFTemplateGenerator
    from openff.toolkit import Molecule
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS

    # 1. pose PDBQT -> SDF (obabel percibe enlaces, mantiene coordenadas)
    pose_sdf = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False).name
    sh(f"obabel {pose_path} -O {pose_sdf} 2>/dev/null")

    docked = Chem.MolFromMolFile(pose_sdf, removeHs=True, sanitize=True)
    if docked is None:
        raise RuntimeError(f"obabel no leyó la pose: {pose_path}")

    # 2. correspondencia de átomos pesados SMILES <-> pose vía MCS
    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol.GetNumAtoms() != docked.GetNumAtoms():
        raise RuntimeError(
            f"nº átomos no coincide: SMILES {rdmol.GetNumAtoms()} vs "
            f"pose {docked.GetNumAtoms()}")
    mcs = rdFMCS.FindMCS([rdmol, docked],
                         bondCompare=rdFMCS.BondCompare.CompareAny,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         ringMatchesRingOnly=True, completeRingsOnly=True,
                         timeout=15)
    if mcs.numAtoms != rdmol.GetNumAtoms():
        raise RuntimeError(f"MCS incompleto: {mcs.numAtoms}/{rdmol.GetNumAtoms()}")
    pat = Chem.MolFromSmarts(mcs.smartsString)
    m1 = rdmol.GetSubstructMatch(pat)
    m2 = docked.GetSubstructMatch(pat)

    conf = docked.GetConformer()
    coordmap = {}
    for i in range(len(m1)):
        p = conf.GetAtomPosition(m2[i])
        coordmap[m1[i]] = Chem.rdGeometry.Point3D(p.x, p.y, p.z)

    # 3. añadir H con los pesados fijos en las coordenadas acopladas
    rdmol_h = Chem.AddHs(rdmol)
    AllChem.EmbedMolecule(rdmol_h, coordMap=coordmap, randomSeed=42,
                          useRandomCoords=True, maxIterations=500)

    # 4. openff + plantilla GAFF (cargas AM1-BCC vía antechamber)
    offmol = Molecule.from_rdkit(rdmol_h, allow_undefined_stereo=True)
    generator = GAFFTemplateGenerator(molecules=[offmol],
                                      forcefield="gaff-2.11")
    forcefield.registerTemplateGenerator(generator.generator)

    top = offmol.to_topology().to_openmm()
    pos = offmol.conformers[0]
    return top, pos


# ----------------------------------------------------------------------------
# MM-GBSA
# ----------------------------------------------------------------------------
def build_system(topology, forcefield):
    from openmm.app import NoCutoff, HBonds
    return forcefield.createSystem(topology, nonbondedMethod=NoCutoff,
                                   constraints=HBonds)


def minimize(system, topology, positions, max_iterations=1000):
    from openmm import LangevinIntegrator, LocalEnergyMinimizer, unit
    from openmm.app import Simulation

    integrator = LangevinIntegrator(300 * unit.kelvin,
                                    1 / unit.picosecond,
                                    2 * unit.femtoseconds)
    sim = Simulation(topology, system, integrator)
    sim.context.setPositions(positions)
    LocalEnergyMinimizer.minimize(sim.context, maxIterations=max_iterations)
    state = sim.context.getState(getPositions=True)
    return state.getPositions(asNumpy=True)


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
    from openmm.app import Modeller

    lig_top, lig_pos = prepare_ligand(smiles, pose_path, forcefield)
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--candidates", default="candidatos_42.csv")
    ap.add_argument("--receptores", default="receptores.json")
    ap.add_argument("--out", default="rescoring_mmgbsa.csv")
    ap.add_argument("--pilot", action="store_true")
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
    if args.pilot:
        df = df.head(1)

    rows = []
    for idx, r in df.iterrows():
        tgt, lig, smiles = r["target"], r["ligand"], r["smiles"]
        pose = r["pose_pdbqt"]
        if tgt not in receptores:
            continue
        try:
            c = ligand_center(pose)
            tmp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb",
                                                  delete=False).name
            extract_receptor_pdb(receptores[tgt], c, tmp_pdb)
            rec_top, rec_pos = prepare_receptor(tmp_pdb)
            dg, ec, er, el = rescore_one(rec_top, rec_pos, smiles, pose,
                                         forcefield)
            rows.append({"ligand": lig, "target": tgt,
                         "vina_affinity": r.get("affinity"),
                         "mmgbsa_dG": round(dg, 2),
                         "e_complex": round(ec, 2),
                         "e_receptor": round(er, 2),
                         "e_ligand": round(el, 2)})
            print(f"[{idx+1}/{len(df)}] {lig} {tgt} dG={dg:.2f}", flush=True)
        except Exception as e:
            rows.append({"ligand": lig, "target": tgt,
                         "vina_affinity": r.get("affinity"),
                         "mmgbsa_dG": None, "error": str(e)[:200]})
            print(f"[{idx+1}/{len(df)}] {lig} {tgt} ERROR: {str(e)[:200]}",
                  flush=True)

    out = pd.DataFrame(rows).sort_values("mmgbsa_dG", na_position="last")
    out.to_csv(args.out, index=False)
    print("Guardado:", args.out)
    print(out.head(15).to_string(index=False))


if __name__ == "__main__":
    main()
