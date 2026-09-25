#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Nombres de atomos de la cadena H (residuo 54 THR) tras preparar el receptor."""
import math
import sys
import tempfile

sys.path.insert(0, "/home/ubuntu/mmgbsa")

import pandas as pd  # noqa: E402

from rescoring_mmgbsa_robusto import ligand_center, prepare_receptor  # noqa


def recorte(pdb, center, out, radio):
    keep = set()
    lineas = []
    for line in open(pdb, errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        d = math.dist((float(line[30:38]), float(line[38:46]),
                       float(line[46:54])), tuple(center))
        if d < radio:
            keep.add((line[21], int(line[22:26])))
        lineas.append(line)
    idx, serial = [], 0
    for n, line in enumerate(lineas):
        if (line[21], int(line[22:26])) in keep:
            serial += 1
            idx.append(line[:6] + "%5d" % serial + line[11:].rstrip())
    open(out, "w").write("\n".join(idx) + "\nEND\n")
    return out


def main():
    lig = sys.argv[1] if len(sys.argv) > 1 else "DEC_CHEMBL1414576"
    rows = pd.read_csv("estratos/candidatos.csv")
    fila = rows[rows.ligand == lig].iloc[0]
    c = ligand_center(fila["pose_pdbqt"])
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    recorte("estratos/receptor.pdb", c, tmp, 14.0)
    print(f"{lig}: recorte en {tmp}")

    print("--- residuo H54 THR en el recorte, tal cual:")
    for line in open(tmp):
        if line.startswith("ATOM") and line[21] == "H":
            print("   ", line[:30].rstrip())

    top, pos = prepare_receptor(tmp)
    print("--- tras PDBFixer:")
    for i, r in enumerate(top.residues()):
        if i < len(list(top.residues())) - 3 or r.name == "THR":
            nombres = [a.name for a in r.atoms()]
            print(f"   #{i+1} {r.chain.id}{r.id} {r.name}: {nombres}")


if __name__ == "__main__":
    main()
