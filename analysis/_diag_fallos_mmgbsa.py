#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Diagnóstico de los 8 fallos de preparación del MM-GBSA por estratos.

Mira, para cada ligando que falló:
  - si el .pdbqt de la pose trae bytes NUL (basura de Vina)
  - qué saca obabel al convertirlo (con y sin limpiar los NUL)
  - cuántos átomos pesados tiene la pose frente a los del SMILES del CSV
  - la composición por elemento de los dos, para ver qué sobra
"""
import csv
import os
import subprocess
import sys
import tempfile

from rdkit import Chem

AD_TYPE2ELEM = {
    "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N", "NX": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "P": "P",
    "F": "F", "CL": "CL", "BR": "BR", "I": "I",
    "H": "H", "HD": "H", "HS": "H",
}


def contar_elementos(atoms):
    d = {}
    for e in atoms:
        d[e] = d.get(e, 0) + 1
    return dict(sorted(d.items()))


def elementos_de_pdbqt(path):
    syms = []
    for line in open(path, errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        adt = line[77:80].strip().upper()
        el = AD_TYPE2ELEM.get(adt, adt).capitalize()
        if el == "H":
            continue
        syms.append(el)
    return syms


def obabel_sdf(pose_path, out_sdf, limpiar):
    if limpiar:
        raw = open(pose_path, "rb").read().replace(b"\x00", b"")
        src = tempfile.NamedTemporaryFile(suffix=".pdbqt", delete=False).name
        open(src, "wb").write(raw)
    else:
        src = pose_path
    r = subprocess.run(f"obabel {src} -O {out_sdf} 2>&1",
                       shell=True, capture_output=True, text=True)
    msg = (r.stdout + r.stderr).strip().replace("\n", " | ")[:200]
    if limpiar:
        os.unlink(src)
    return msg


def main(csv_path, poses_dir, csv_smiles=None):
    filas = list(csv.DictReader(open(csv_path, encoding="utf-8-sig")))
    if csv_smiles:
        smis = {r["ligand"]: r["smiles"]
                for r in csv.DictReader(open(csv_smiles, encoding="utf-8-sig"))}
        for r in filas:
            r.setdefault("smiles", smis.get(r["ligand"], ""))
    fallos = [r for r in filas if (r.get("error") or "").strip()]
    print(f"fallos: {len(fallos)}\n")
    for r in fallos:
        lig, smi, err = r["ligand"], r["smiles"], r["error"]
        pose = os.path.join(poses_dir, f"{lig}_out.pdbqt")
        print("=" * 70)
        print(f"{lig}\n  tipo fallo: {err[:70]}")

        raw = open(pose, "rb").read()
        n_nul = raw.count(bytes([0]))
        print(f"  pose: {len(raw)} bytes, NUL = {n_nul}")
        print(f"  pdbqt pesados: {contar_elementos(elementos_de_pdbqt(pose))}")

        rd = Chem.MolFromSmiles(smi)
        print(f"  SMILES pesados={rd.GetNumAtoms() if rd else None} "
              f"{contar_elementos([a.GetSymbol() for a in rd.GetAtoms()]) if rd else smi}")

        for limpiar in (False, True):
            sdf = tempfile.NamedTemporaryFile(suffix=".sdf", delete=False).name
            msg = obabel_sdf(pose, sdf, limpiar)
            m = Chem.MolFromMolFile(sdf, removeHs=True, sanitize=True)
            if m is None:
                print(f"  obabel limpiar={limpiar}: None  ({msg})")
            else:
                print(f"  obabel limpiar={limpiar}: pesados={m.GetNumAtoms()} "
                      f"{contar_elementos([a.GetSymbol() for a in m.GetAtoms()])}")
            os.unlink(sdf)
        print()


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
