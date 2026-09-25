#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿Cambia el receptor preparado de un proceso a otro?

Imprime un resumen numerico de las coordenadas del receptor preparado para
ACT_isoproterenol (atomo, suma de coordenadas). Lanzado varias veces en
procesos distintos: si los numeros cambian, la preparacion depende del proceso
(sospecha: orden de iteracion de conjuntos, PYTHONHASHSEED); si coinciden, la
preparacion es reproducible.
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, "/home/ubuntu/mmgbsa")

from rescoring_mmgbsa_robusto import (  # noqa: E402
    _recortar_pdb_canonical, ligand_center, prepare_receptor)

import pandas as pd  # noqa: E402


def main():
    fila = pd.read_csv("estratos/candidatos.csv")
    fila = fila[fila.ligand == "ACT_isoproterenol"].iloc[0]
    c = ligand_center(fila["pose_pdbqt"])
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    _recortar_pdb_canonical("estratos/receptor.pdb", c, tmp, radio=14.0)
    top, pos = prepare_receptor(tmp)
    a = np.array(pos.value_in_unit(pos.unit))
    print(f"PYTHONHASHSEED={os.environ.get('PYTHONHASHSEED', 'sin fijar')} "
          f"pid={os.getpid()} atomos={len(a)} "
          f"suma={a.sum():.6f} suma2={float((a ** 2).sum()):.6f} "
          f"primeros={np.round(a[0], 6)}", flush=True)


if __name__ == "__main__":
    main()
