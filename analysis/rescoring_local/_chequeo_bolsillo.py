# -*- coding: utf-8 -*-
"""Comprueba que los ligandos de la validación caen DENTRO del bolsillo fijo.

El receptor congelado es una bolsa extraída una sola vez. Si la pose de un
control quedara fuera de esa bolsa, se estaría midiendo la energía de un
ligando flotando en disolvente, no de un complejo: el dG no significaría
nada y la comparación controles/fondo sería inválida.

Se compara el centro de masa del ligando (pose acoplada) con el centro de la
bolsa y con su radio, para controles y fondo por separado.
"""
import csv
import os
import sys

import numpy as np

BASE = r"C:\Users\Fredy\masive-als"
if not os.path.isdir(BASE):
    BASE = "/mnt/c/Users/Fredy/masive-als"
LOCAL = os.path.join(BASE, "analysis", "rescoring_local")
RES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
FIJO = os.path.join(LOCAL, "receptores_fijos")
VAL = os.path.join(LOCAL, "validacion_controles.csv")


def atomos_pdb(path):
    xs = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith(("ATOM", "HETATM")):
                try:
                    xs.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                except ValueError:
                    pass
    return np.array(xs)


def atomos_pdbqt(path):
    xs = []
    with open(path, encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("ENDMDL"):     # solo la primera pose
                break
            if l.startswith(("ATOM", "HETATM")):
                try:
                    xs.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                except ValueError:
                    pass
    return np.array(xs)


def main():
    filas = list(csv.DictReader(open(VAL, encoding="utf-8")))
    recs = {}
    for t in sorted({r["target"] for r in filas}):
        p = os.path.join(FIJO, "%s_fijo.pdb" % t)
        a = atomos_pdb(p)
        centro = a.mean(axis=0)
        radio = float(np.linalg.norm(a - centro, axis=1).max())
        recs[t] = (centro, radio, len(a))
        print("bolsillo %-9s %5d atomos | radio %.1f A" % (t, len(a), radio))

    for t in sorted({r["target"] for r in filas}):
        centro, radio, _ = recs[t]
        print("\n--- %s ---" % t)
        for tipo in ("control", "fondo"):
            d = []
            for r in filas:
                if r["target"] != t or r["tipo"] != tipo:
                    continue
                pose = os.path.join(RES, "results_" + t, r["ligand"] + "_out.pdbqt")
                if not os.path.exists(pose):
                    continue
                a = atomos_pdbqt(pose)
                if len(a) == 0:
                    continue
                d.append(float(np.linalg.norm(a.mean(axis=0) - centro)))
            if not d:
                print("  %-8s sin poses" % tipo)
                continue
            d = np.array(d)
            dentro = int((d <= radio).sum())
            print("  %-8s n=%3d | distancia al centro: media %.1f  min %.1f  max %.1f"
                  " | dentro del radio: %d/%d"
                  % (tipo, len(d), d.mean(), d.min(), d.max(), dentro, len(d)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
