#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rescoring MM-GBSA v6 (robusto) — MASIVE-ALS.

v6 (22 sep 2026, tarde): la preparacion del receptor va con semilla fija y el
recorte se centra en la primera pose.
Cada proceso preparaba su propio receptor y dos corridas de la misma fila daban
coordenadas distintas (hasta 0,24 A) y dG que se movian varios kcal/mol
(ACT_isoproterenol: -15,65 / -13,52 / -7,69 en tres corridas). El azar estaba
donde menos se esperaba, medido paso a paso: `addMissingAtoms` con semilla si es
determinista (0,000 A), pero **`addMissingHydrogens` coloca los hidrogenos con
el generador aleatorio de Python** (302 de 875 atomos se movian >0,01 A). Se
siembran `random` y `numpy.random` antes de anadir los hidrogenos y la
preparacion pasa a ser identica entre corridas (0,001 A).

Y de paso se arreglo un fallo de bulto en `ligand_center`: los .pdbqt de Vina
traen 3 MODEL y se promediaban los tres, con lo que el centro caia a 7,3 A del
sitio real y el recorte del receptor se hacia alrededor de otro punto. Ahora se
usa solo el MODEL 1, que es la pose buena.

v5: igual que v4 pero con los dos arreglos que sacaron los 8 fallos del
rescoring por estratos del 22 sep 2026:

  1. La pose ya no se convierte con obabel. Vina anota en el propio .pdbqt el
     mapa entre su SMILES y los atomos del fichero (`REMARK SMILES` y
     `REMARK SMILES IDX`), asi que las coordenadas se colocan con ese mapa
     (modulo `pose_pdbqt`). obabel fallaba con los bytes NUL de relleno y no
     acertaba con amonios cuaternarios ni amidinas.
  2. El recorte del receptor deja fuera los residuos que se quedan solos en su
     cadena. Un residuo suelto no tiene terminales validas: PDBFixer le pone a
     la vez los H del N-terminal y el OXT del C-terminal, y ninguna plantilla de
     amber14 encaja ('No template found for residue N (THR)' en la H54 de SOD1).
     Un residuo suelto no se puede arreglar ni renombrando ni quitando atomos:
     le faltan las cadenas vecinas, asi que no se cuela un enlace falso.

Lo demas sigue igual que en v3/v4: guardas anti-NaN en la minimizacion, modo
`--row N` para una sola fila con JSON (lo usa el runner) y la tanda no se cae
entera si una pose falla.

CARGAS — CORRECCION (22 sep 2026, noche): el parrafo viejo de aqui abajo decia que el
ligando iba SIN cargas ("Set zero charges to skip slow AM1-BCC calculation"). Es un
no-op fisico: GAFFTemplateGenerator SIEMPRE calcula AM1-BCC al construir la plantilla
(openmmforcefields 0.16.0, template_generators.py:594) y sobrescribe lo que se le haya
puesto. Verificado: el NonbondedForce montado con "cargas a cero" tiene la misma carga
maxima (0.827 e) que con AM1-BCC explícito, y una variante con AM1-BCC explícito
(rescoring_mmgbsa_bcc.py) reproduce bit a bit tres filas de la tanda. La tabla de los
90 lleva cargas AM1-BCC reales; el "atajo" nunca ahorró el tiempo que decia ahorrar.

SIN CARGAS EN EL LIGANDO (texto historico, FALLO en la fisica): `prepare_ligand` pone
las cargas parciales a cero para saltarse el AM1-BCC. Lo que se calcula es GAFF 2.11
con el ligando neutro. No vale llamarlo MM-GBSA de referencia.
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

# semilla de la preparacion del receptor: PDBFixer.addMissingAtoms y el
# generador aleatorio de Python que usa addMissingHydrogens
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
    """Recorta un PDB canonico a los residuos con algun atomo a < radio del ligando.

    Si de una cadena solo queda un residuo, ese residuo se cae del recorte: una
    cadena de un residuo no tiene terminales validas, PDBFixer le pone a la vez
    los H del N-terminal y el OXT del C-terminal, y ninguna plantilla de amber14
    encaja (el fallo "No template found for residue N (THR)" de la H54 de SOD1).
    Se cae el residuo suelto, no la cadena entera, y solo cuando queda solo.
    """
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
    from pdbfixer import PDBFixer
    import random
    import tempfile, os

    # Strip Ca atoms and Ca-containing HETATM residues before PDBFixer
    # (OpenMM amber14 force field has no template for Ca)
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
    # addMissingHydrogens tira del generador aleatorio de Python: sin sembrarlo,
    # dos corridas de la misma fila dan receptores distintos (hasta 0,24 A) y el
    # dG se mueve varios kcal/mol. Sembrado, la preparacion es identica.
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
    """Centro del ligando en su PRIMERA pose.

    Ojo: un .pdbqt de Vina trae varios MODEL (aqui 3), y el centro se calcula
    con el primero, que es la pose buena. Promediando los tres modelos el
    centro cae a 7,3 A del sitio real (el modelo 2 estaba en otro bolsillo) y el
    recorte del receptor de 14 A se hacia alrededor de un punto equivocado.
    """
    import pose_pdbqt

    pose = pose_pdbqt.leer_pose(pose_path, modelo=1)
    at = list(pose["atoms"].values())
    if not at:
        raise RuntimeError("pose vacía: " + pose_path)
    return np.array([[a[0], a[1], a[2]] for a in at]).mean(axis=0)


def prepare_ligand(smiles, pose_path, forcefield):
    from openmmforcefields.generators import GAFFTemplateGenerator
    from openff.toolkit import Molecule
    from rdkit import Chem
    from rdkit.Chem import AllChem

    import pose_pdbqt

    rdmol = Chem.MolFromSmiles(smiles)
    if rdmol is None:
        raise RuntimeError(f"SMILES ilegible: {smiles[:60]}")

    # El grafo lo pone el SMILES del CSV y las coordenadas la pose, con el mapa
    # de atomos que Vina escribe en el propio .pdbqt (nada de obabel).
    pose = pose_pdbqt.leer_pose(pose_path)
    coordmap = pose_pdbqt.coordmap_de_pose(pose, rdmol)

    rdmol_h = Chem.AddHs(rdmol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    params.maxIterations = 500
    params.SetCoordMap(coordmap)
    ret = AllChem.EmbedMolecule(rdmol_h, params)
    # Si el embedding con las coordenadas de la pose falla (ret != 0, pasa con
    # ligandos flexibles: ACT_adrenalina), se genera la geometria sin
    # restricciones y despues se lleva el ligando al bolsillo superponiendo los
    # atomos mapeados sobre las coordenadas de la pose. Sin esa superposicion el
    # ligando se queda flotando fuera del bolsillo y el dG sale ~0.
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
    # Set zero charges to skip slow AM1-BCC calculation
    import numpy as np
    from openff.units import unit
    n_atoms = offmol.n_atoms
    zero_charges = np.zeros(n_atoms) * unit.elementary_charge
    offmol.partial_charges = zero_charges
    generator = GAFFTemplateGenerator(molecules=[offmol],
                                      forcefield="gaff-2.11")
    forcefield.registerTemplateGenerator(generator.generator)

    from openmm import unit as omm_unit
    top = offmol.to_topology().to_openmm()
    pos = offmol.conformers[0].magnitude * omm_unit.angstrom
    return top, pos


def _llevar_a_la_pose(rdmol_h, coordmap):
    """Superpone los atomos mapeados del ligando a las coordenadas de la pose."""
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


def compute_row(r, receptores, forcefield):
    """Calcula una fila y devuelve dict con dG o error."""
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
    ap.add_argument("--out", default="rescoring_mmgbsa.csv")
    ap.add_argument("--pilot", action="store_true")
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
    if args.pilot:
        df = df.head(1)

    # ---- modo una sola fila (usado por el runner con timeout) ----
    if args.row is not None:
        r = df.iloc[args.row]
        out = compute_row(r, receptores, forcefield)
        # normalizar None -> null para JSON
        print(json.dumps(out))
        sys.exit(0)

    # ---- modo tanda completa ----
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
