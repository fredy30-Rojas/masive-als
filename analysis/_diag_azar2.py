#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿En que paso se vuelve azarosa la preparacion del receptor?

Hace la preparacion dos veces, paso a paso (PDBFixer, buscar faltantes, anadir
atomos con semilla, anadir hidrogenos) y compara las coordenadas despues de cada
paso, diciendo que atomos se mueven mas.
"""
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


def compara(etapa, t1, p1, t2, p2):
    a, b = pos(p1), pos(p2)
    if a.shape != b.shape:
        print(f"{etapa}: FORMAS DISTINTAS {a.shape} vs {b.shape}", flush=True)
        return None
    d = np.linalg.norm(a - b, axis=1)
    n = int((d > 0.01).sum())
    print(f"{etapa}: {len(d)} atomos, {n} se mueven >0,01 A, "
          f"max={d.max():.3f} A", flush=True)
    if n:
        idx = np.argsort(-d)[:6]
        for i in idx:
            at = list(t1.atoms())[i]
            print(f"    {at.residue.name}{at.residue.id} {at.name}: "
                  f"{d[i]:.3f} A", flush=True)
    return d


def prepara(pdb, paso):
    fixer = PDBFixer(filename=pdb)
    fixer.missingResidues = []
    fixer.removeHeterogens(keepWater=False)
    if paso < 1:
        return fixer.topology, fixer.positions
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    if paso < 2:
        return fixer.topology, fixer.positions
    fixer.addMissingAtoms(seed=SEMILLA_RECEPTOR)
    if paso < 3:
        return fixer.topology, fixer.positions
    fixer.addMissingHydrogens(7.0)
    return fixer.topology, fixer.positions


def main():
    fila = pd.read_csv("estratos/candidatos.csv")
    fila = fila[fila.ligand == "ACT_isoproterenol"].iloc[0]
    c = ligand_center(fila["pose_pdbqt"])
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    _recortar_pdb_canonical("estratos/receptor.pdb", c, tmp, radio=14.0)

    nombres = {0: "tras removeHeterogens",
               1: "tras findMissingAtoms",
               2: "tras addMissingAtoms(seed)",
               3: "tras addMissingHydrogens"}
    previos = {}
    for paso in (0, 1, 2, 3):
        t1, p1 = prepara(tmp, paso)
        t2, p2 = prepara(tmp, paso)
        compara(nombres[paso], t1, p1, t2, p2)
        previos[paso] = (t1, p1)
        if paso > 0:
            t_prev, p_prev = previos[paso - 1]
            if len(list(t1.atoms())) == len(list(t_prev.atoms())):
                compara(f"  (paso {paso} contra el paso anterior, una sola vez)",
                        t1, p1, t_prev, p_prev)


if __name__ == "__main__":
    main()
