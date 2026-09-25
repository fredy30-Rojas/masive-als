#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿De dónde sale el ruido del MM-GBSA: del recorte o del PDBFixer?

Hace dos veces el recorte y la preparación del receptor para el mismo ligando y
compara: número de átomos, coordenadas, y si PDBFixer admite semilla.
"""
import inspect
import sys
import tempfile

import numpy as np

sys.path.insert(0, "/home/ubuntu/mmgbsa")

from rescoring_mmgbsa_robusto import (  # noqa: E402
    _recortar_pdb_canonical, ligand_center, prepare_receptor)

import pandas as pd  # noqa: E402

LIG = "ACT_isoproterenol"


def dif(a, b):
    a = np.array(a.value_in_unit(a.unit))
    b = np.array(b.value_in_unit(b.unit))
    if a.shape != b.shape:
        return f"formas distintas {a.shape} vs {b.shape}"
    return f"max|d|={np.abs(a - b).max():.3e} A"


def main():
    fila = pd.read_csv("estratos/candidatos.csv")
    fila = fila[fila.ligand == LIG].iloc[0]
    c = ligand_center(fila["pose_pdbqt"])
    print("centro del ligando:", np.round(c, 3), flush=True)

    # --- el recorte, dos veces (mismo centro) ---
    pdb_txt = []
    for rep in (1, 2):
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        _recortar_pdb_canonical("estratos/receptor.pdb", c, tmp, radio=14.0)
        pdb_txt.append(open(tmp).read())
    iguales = pdb_txt[0] == pdb_txt[1]
    print(f"recorte: {len(pdb_txt[0].splitlines())} vs "
          f"{len(pdb_txt[1].splitlines())} lineas -> "
          f"{'IDENTICO' if iguales else 'DISTINTO'}", flush=True)

    # --- la preparacion (PDBFixer), dos veces ---
    tops, poss = [], []
    for rep in (1, 2):
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        _recortar_pdb_canonical("estratos/receptor.pdb", c, tmp, radio=14.0)
        t, p = prepare_receptor(tmp)
        tops.append(t)
        poss.append(p)
        print(f"  rep{rep}: {len(list(t.atoms()))} atomos", flush=True)
    print("receptor preparado:", dif(poss[0], poss[1]), flush=True)

    # --- ¿PDBFixer admite semilla? ---
    from pdbfixer import PDBFixer
    for nombre in ("addMissingAtoms", "addMissingHydrogens", "findMissingAtoms"):
        try:
            print(f"  PDBFixer.{nombre}{inspect.signature(getattr(PDBFixer, nombre))}",
                  flush=True)
        except Exception as e:
            print(f"  PDBFixer.{nombre}: sin firma ({e})", flush=True)


if __name__ == "__main__":
    main()
