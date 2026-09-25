#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Lectura robusta de poses PDBQT de Vina para el MM-GBSA.

POR QUE NO OBABEL
-----------------
Los .pdbqt que escribe Vina traen bytes NUL de relleno y varios MODEL, y ademas
obabel no acierta con la percepcion de enlaces en algunas quimicas (amonio
cuaternario, amidinas, guanidinas). De ahi venian los 6 fallos de ligando del
MM-GBSA por estratos:

  * "obabel no leyo la pose"  -> el N atom #7 con 4 enlaces no pasa sanitizacion
    en RDKit (acetilcolina, [N+](C)(C)C), y la pose con NUL ni se lee.
  * "n atomos no coincide"    -> obabel deja 1 o 2 hidrogenos sin quitar y RDKit
    los cuenta como pesados.

LA SOLUCION
-----------
Vina ya escribe en el propio fichero el mapa entre los atomos de su SMILES de
entrada y los atomos del .pdbqt:

    REMARK SMILES <smiles de entrada>
    REMARK SMILES IDX <smi_idx pdbqt_idx> ...   (pares, solo atomos pesados)
    REMARK H PARENT <pdbqt_idx_h parent_idx> ... (que atomos son hidrogeno)

Con eso se colocan las coordenadas de la pose sobre los atomos del SMILES del
CSV sin depender de obabel ni de percibir enlaces: el grafo lo pone el SMILES y
las coordenadas la pose.
"""
import sys

from rdkit import Chem
from rdkit.Geometry import Point3D

AD_TYPE2ELEM = {
    "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N", "NX": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "P": "P",
    "F": "F", "CL": "CL", "BR": "BR", "I": "I",
    "H": "H", "HD": "H", "HS": "H",
}


def leer_pose(path, modelo=1):
    """Lee un MODEL del .pdbqt ignorando los NUL y devuelve atomos y su mapa.

    Devuelve dict con:
      atoms      {indice pdbqt (1-based): (x, y, z, elemento)} solo pesados
      h_parents  [(indice pdbqt del H, indice pdbqt del padre)]
      vina_smiles  SMILES de entrada que anota Vina
      idx_pairs  [(indice en el SMILES de Vina, indice en el pdbqt)]
                 (Vina lo parte en varias lineas; aqui van acumulados)
      n_modelos  cuantos MODEL trae el fichero
    """
    crudo = open(path, "rb").read().replace(bytes([0]), b"")
    lineas = crudo.decode("latin-1").splitlines()

    n_modelos = sum(1 for l in lineas if l.startswith("MODEL"))
    actual = 0
    leyendo = (modelo is None)
    atoms, h_parents = {}, []
    vina_smiles, idx_pairs = None, []
    for line in lineas:
        if line.startswith("MODEL"):
            actual += 1
            leyendo = (modelo is None or actual == modelo)
            continue
        if not leyendo:
            continue
        if line.startswith(("ATOM", "HETATM")):
            i = int(line[6:11])
            adt = line[77:80].strip().upper()
            elem = AD_TYPE2ELEM.get(adt, adt)
            if elem == "H":
                continue
            atoms[i] = (float(line[30:38]), float(line[38:46]),
                        float(line[46:54]), elem)
        elif line.startswith("REMARK SMILES IDX"):
            # Vina parte esta linea en trozos de ~76 caracteres: hay que sumarlos
            idx_pairs = (idx_pairs or []) + [int(v) for v in line.split()[3:]]
        elif line.startswith("REMARK SMILES"):
            vina_smiles = line.split("REMARK SMILES", 1)[1].strip()
        elif line.startswith("REMARK H PARENT"):
            nums = h_parents + [int(v) for v in line.split()[3:]]
            h_parents = list(zip(nums[0::2], nums[1::2]))
        elif line.startswith("ENDMDL") and modelo is not None:
            break

    if not atoms:
        raise RuntimeError("pose sin atomos pesados: " + path)
    return dict(atoms=atoms, h_parents=h_parents, vina_smiles=vina_smiles,
                idx_pairs=idx_pairs, n_modelos=n_modelos)


def _match_sin_h(molv, rdmol):
    """Empareja los atomos del SMILES de Vina con los del SMILES del CSV.

    Devuelve (match, pos): `match[k]` es el indice en rdmol del k-esimo atomo
    pesado de molv, y `pos[i]` la posicion del atomo i de molv en esa lista (los
    hidrogenos explicitos, si los hay, se quitan para poder emparejar).
    """
    pesados = [a.GetIdx() for a in molv.GetAtoms() if a.GetSymbol() != "H"]
    if len(pesados) == molv.GetNumAtoms():
        mol_sin_h = molv
    else:
        mol_sin_h = Chem.RemoveHs(molv)
        if mol_sin_h.GetNumAtoms() != len(pesados):
            raise RuntimeError("no se pudo quitar los H del SMILES de la pose")

    match = rdmol.GetSubstructMatch(mol_sin_h)
    if not match or len(match) != mol_sin_h.GetNumAtoms():
        raise RuntimeError(
            f"el SMILES de la pose no encaja en el del CSV "
            f"({len(match)}/{mol_sin_h.GetNumAtoms()} atomos)")
    pos = {idx: k for k, idx in enumerate(pesados)}
    return match, pos


def coordmap_de_pose(pose, rdmol):
    """Coloca las coordenadas de la pose sobre los atomos de `rdmol`.

    `rdmol` es la molecula del SMILES del CSV (sin hidrogenos). Devuelve
    {indice en rdmol: Point3D}.
    """
    if not pose["vina_smiles"]:
        raise RuntimeError("la pose no trae la linea REMARK SMILES")
    if not pose["idx_pairs"]:
        raise RuntimeError("la pose no trae la linea REMARK SMILES IDX")

    molv = Chem.MolFromSmiles(pose["vina_smiles"])
    if molv is None:
        raise RuntimeError("SMILES de la pose ilegible: " + pose["vina_smiles"])
    match, pos = _match_sin_h(molv, rdmol)
    n_pesados_v = len(pos)

    pares = list(zip(pose["idx_pairs"][0::2], pose["idx_pairs"][1::2]))
    if len(pares) != n_pesados_v:
        raise RuntimeError(
            f"el mapa de Vina tiene {len(pares)} pares y el SMILES de la pose "
            f"{n_pesados_v} atomos pesados")

    mapa = {}
    for i_smi, i_pdb in pares:
        if not 1 <= i_smi <= molv.GetNumAtoms():
            raise RuntimeError(f"indice de SMILES fuera de rango: {i_smi}")
        if i_pdb not in pose["atoms"]:
            raise RuntimeError(f"indice de pdbqt fuera de rango: {i_pdb}")
        k = pos.get(i_smi - 1)
        if k is None:
            raise RuntimeError(f"el mapa apunta al hidrogeno {i_smi}")
        x, y, z, elem = pose["atoms"][i_pdb]
        j = match[k]
        esp = rdmol.GetAtomWithIdx(j).GetSymbol().upper()
        if esp != elem.upper():
            raise RuntimeError(
                f"elemento distinto en el atomo {i_smi}: pose {elem}, "
                f"CSV {esp}")
        mapa.setdefault(j, Point3D(x, y, z))

    if len(mapa) != rdmol.GetNumAtoms():
        raise RuntimeError(
            f"cobertura incompleta: {len(mapa)}/{rdmol.GetNumAtoms()} atomos")
    return mapa


def _verificar(csv_path, poses_dir=None):
    """Comprueba el mapa de Vina en todas las poses del CSV."""
    import csv
    import os
    import traceback

    filas = list(csv.DictReader(open(csv_path, encoding="utf-8-sig")))
    malos = []
    sin_idx = 0
    con_h = 0
    for r in filas:
        pose_path = r["pose_pdbqt"]
        if poses_dir:
            pose_path = os.path.join(poses_dir, os.path.basename(pose_path))
        try:
            pose = leer_pose(pose_path)
            if not pose["idx_pairs"]:
                sin_idx += 1
                raise RuntimeError("sin REMARK SMILES IDX")
            if any(a.GetSymbol() == "H"
                   for a in Chem.MolFromSmiles(pose["vina_smiles"]).GetAtoms()):
                con_h += 1
            rdmol = Chem.MolFromSmiles(r["smiles"])
            if rdmol is None:
                raise RuntimeError("SMILES del CSV ilegible")
            m = coordmap_de_pose(pose, rdmol)
            n = len(m)
            if n != rdmol.GetNumAtoms():
                raise RuntimeError(f"cobertura {n}/{rdmol.GetNumAtoms()}")
        except Exception as e:
            malos.append((r["ligand"], str(e)[:120]))
            continue
    print(f"poses: {len(filas)}  fallos: {len(malos)}  "
          f"sin IDX: {sin_idx}  con H explicito: {con_h}")
    for lig, err in malos:
        print(f"  {lig}: {err}")
    return malos


if __name__ == "__main__":
    _verificar(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
