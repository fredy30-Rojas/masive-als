#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿Basta con sembrar random/numpy antes de addMissingHydrogens?

Dos pruebas: (a) sin sembrar nada, (b) sembrando `random` y `numpy.random` justo
antes de anadir los hidrogenos. Si (b) sale identico, el arreglo es una semilla
de Python y el numero pasa a ser repetible.
"""
import random
import sys
import tempfile

import numpy as np

sys.path.insert(0, "/home/ubuntu/mmgbsa")

from rescoring_mmgbsa_robusto import (  # noqa: E402
    SEMILLA_RECEPTOR, _recortar_pdb_canonical, ligand_center)

import pandas as pd  # noqa: E402
from pdbfixer import PDBFixer  # noqa: E402


def pos(p):
    return np.array(p.value_in_unit(p.unit))


def prepara(pdb, sembrar):
    fixer = PDBFixer(filename=pdb)
    fixer.missingResidues = []
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms(seed=SEMILLA_RECEPTOR)
    if sembrar is not None:
        random.seed(sembrar)
        np.random.seed(sembrar)
    fixer.addMissingHydrogens(7.0)
    return fixer.topology, fixer.positions


def compara(etiqueta, t1, p1, t2, p2):
    a, b = pos(p1), pos(p2)
    d = np.linalg.norm(a - b, axis=1)
    print(f"{etiqueta}: {len(d)} atomos, {int((d > 0.01).sum())} se mueven "
          f">0,01 A, max={d.max():.3f} A", flush=True)


def main():
    fila = pd.read_csv("estratos/candidatos.csv")
    fila = fila[fila.ligand == "ACT_isoproterenol"].iloc[0]
    c = ligand_center(fila["pose_pdbqt"])
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    _recortar_pdb_canonical("estratos/receptor.pdb", c, tmp, radio=14.0)

    a = prepara(tmp, None)
    b = prepara(tmp, None)
    compara("sin sembrar", a[0], a[1], b[0], b[1])

    c1 = prepara(tmp, 7)
    c2 = prepara(tmp, 7)
    compara("sembrando random y numpy", c1[0], c1[1], c2[0], c2[1])

    # y de propina: sembrado una vez contra sin sembrar
    compara("sembrado vs sin sembrar", a[0], a[1], c1[0], c1[1])


if __name__ == "__main__":
    main()
