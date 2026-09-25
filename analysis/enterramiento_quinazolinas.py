#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Mide el ENTERRAMIENTO de las poses acopladas (anomalía de las quinazolinas).

El 20 de septiembre de 2026 se encontró que dos quinazolinas cristalizadas en
Trp32 puntúan -10,0 y -9,8 kcal/mol mientras otra de la misma química, con un
anillo de piperazina en vez de diazepano, puntúa -4,97. Cinco kcal/mol por un
CH2 no es físico. Este script mide, sobre las poses ya calculadas, si la
diferencia se explica por **cuánto receptor rodea al ligando** (enterramiento),
que es la explicación trivial y hay que descartarla con un número.

Cuenta, para la mejor pose de cada ligando:
  - átomos pesados del receptor a <= 4,5 A
  - átomos pesados del receptor a <= 5,5 A
  - residuos distintos del receptor en contacto (<= 4,5 A)
  - afinidad, y afinidad por átomo pesado

Uso:
    python analysis/enterramiento_quinazolinas.py
"""
from __future__ import annotations

import glob
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
RECEPTOR = os.path.abspath(os.path.join(BASE, "..", "gpu_dock", "SOD1.pdbqt"))
OUT = os.path.join(BASE, "validacion_SOD1_v4", "out")
INTERES = [
    ("diazepanoquinazolina_4MQ", "ACT_diazepanoquinazolina_4MQ"),
    ("cf3quinazolina_ZO0", "ACT_cf3quinazolina_ZO0"),
    ("quinazolina_12I", "ACT_quinazolina_12I"),
    ("Lig9_6B3", "ACT_Lig9_6B3"),
    ("5-fluorouridina", "ACT_5-fluorouridina"),
    ("isoproterenol", "ACT_isoproterenol"),
]


def atomos_pesados(ruta, solo_primera_pose=False):
    """[(x, y, z, residuo)] de los atomos pesados (H/HD/HS fuera)."""
    out, en_pose, pose = [], 1, 1
    with open(ruta, encoding="utf-8", errors="replace") as f:
        for l in f:
            if l.startswith("MODEL"):
                pose = int(l.split()[1])
                continue
            if not l.startswith(("ATOM", "HETATM")):
                continue
            if solo_primera_pose and pose != 1:
                continue
            if l.rsplit(None, 1)[-1].upper() in ("H", "HD", "HS", "D", "DD"):
                continue
            try:
                xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            except ValueError:
                continue
            res = l[17:20].strip() + l[22:26].strip()
            out.append((*xyz, res))
    return out


def dist2(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def afinidad(ruta):
    """Mejor afinidad de las que Vina escribe como REMARK en el fichero."""
    with open(ruta, encoding="utf-8", errors="replace") as f:
        for l in f:
            if l.startswith("REMARK VINA RESULT"):
                return float(l.split()[3])
    return None


def main():
    rec = atomos_pesados(RECEPTOR)
    print("átomos pesados del receptor: %d\n" % len(rec))
    print("%-30s %8s %6s %10s %10s %9s" %
          ("ligando", "afin", "pesad", "recep<=4.5A", "recep<=5.5A", "residuos"))
    for nombre, fichero in INTERES:
        p = os.path.join(OUT, fichero + "_out.pdbqt")
        if not os.path.exists(p):
            print("%-30s  (sin pose)" % nombre)
            continue
        lig = atomos_pesados(p, solo_primera_pose=True)
        n45 = n55 = 0
        residuos = set()
        for r in rec:
            d2 = min(dist2(r, l) for l in lig)
            if d2 <= 4.5 ** 2:
                n45 += 1
                residuos.add(r[3])
            if d2 <= 5.5 ** 2:
                n55 += 1
        aff = afinidad(p)
        print("%-30s %8.2f %6d %10d %10d %9d" %
              (nombre, aff if aff else float("nan"), len(lig), n45, n55, len(residuos)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
