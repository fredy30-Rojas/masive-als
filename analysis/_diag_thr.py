#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Reproduce el fallo 'No template found for residue N (THR)' en dos ligandos.

Hace lo mismo que compute_row: centra en la pose, recorta el receptor a 14 A,
prepara el receptor con PDBFixer e imprime la traza completa si falla.
"""
import sys
import tempfile
import traceback

sys.path.insert(0, "/home/ubuntu/mmgbsa")

import pandas as pd  # noqa: E402

from rescoring_mmgbsa_robusto import (  # noqa: E402
    _recortar_pdb_canonical, ligand_center, prepare_receptor)

LIGS = ["DEC_CHEMBL1414576", "DECH_CHEMBL3309988", "ACT_isoproterenol"]


def resumen(pdb):
    res = {}
    for line in open(pdb, errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            res.setdefault((line[21], line[22:26].strip(),
                            line[17:20].strip()), 0)
            res[(line[21], line[22:26].strip(), line[17:20].strip())] += 1
    return res


def main():
    rows = pd.read_csv("estratos/candidatos.csv")
    for lig in LIGS:
        fila = rows[rows.ligand == lig]
        if fila.empty:
            print(f"{lig}: no está en el CSV"); continue
        pose = fila.iloc[0]["pose_pdbqt"]
        c = ligand_center(pose)
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        _recortar_pdb_canonical("estratos/receptor.pdb", c, tmp, radio=14.0)
        r = resumen(tmp)
        print("=" * 70)
        print(f"{lig}  centro=({c[0]:.1f},{c[1]:.1f},{c[2]:.1f})  "
              f"recorte={len(r)} residuos")
        for ch in sorted({k[0] for k in r}):
            nums = sorted({int(k[1]) for k in r if k[0] == ch})
            huecos = [f"{a}->{b}" for a, b in zip(nums, nums[1:]) if b != a + 1]
            print(f"   cadena {ch}: {len(nums)} res, huecos={huecos}")
        for (ch, num, rn), n in sorted(r.items()):
            if rn == "THR":
                print(f"   THR {ch}{num}: {n} atomos")
        try:
            top, pos = prepare_receptor(tmp)
            print(f"   receptor OK -> {top.getNumAtoms()} atomos")
        except Exception:
            print("   FALLO:")
            traceback.print_exc()
        print(f"   recorte guardado en {tmp}")


if __name__ == "__main__":
    main()
