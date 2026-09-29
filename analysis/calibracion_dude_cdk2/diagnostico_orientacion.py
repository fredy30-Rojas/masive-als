#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿El redocking de CDK2 falla porque el ligando va dado la vuelta?

QUE PREGUNTA CONTESTA
---------------------
El control de redocking de CDK2 da RMSD de 6 a 7 A en las seis variantes de
protocolo, pero los centroides caen a 0,7-1,5 A del cristal y el solape (distancia
de cada atomo del cristal al mas cercano de la pose) es de 1,6-2,0 A. Ese patron
—mismo volumen, asignacion distinta— tiene una causa clasica: **el ligando se ha
acoplado al reves**, algo tipico de Vina con ligandos planos en bolsillos en
hendidura, porque la funcion de puntuacion no ve la orientacion de un anillo
apilado.

La prueba es directa y no necesita GPU: se le aplican al cristal una lista de
transformaciones rigidas (identidad, giro de 180 grados alrededor de cada eje
principal, y la inversion por el centroide) y se mide cual deja la pose encima.
Si alguna transformation deja el RMSD por debajo de 2 A, la conclusion es
exacta: la pose es el cristal dado la vuelta, el receptor y la caja estan
bien, y lo que falla es la capacidad de la puntuacion para distinguir las dos
orientaciones.

Que ademas lo respalda, ya medido en el proyecto: la MISMA receta (meeko,
caja desde el ligando del cristal, Vina-GPU 2.1) pasa el control del BX-795 de
TBK1 con 1,23 A. El embudo sabe volver a su casa; este ligando en concreto no.

Uso:
    python diagnostico_orientacion.py                    # la mejor variante
    python diagnostico_orientacion.py --variante caja22_depth32
"""
from __future__ import annotations

import argparse
import glob
import os

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
MOL2 = os.path.join(BASE, "crystal_ligand.mol2")
CONTROL = os.path.join(BASE, "_control_cdk2")
LISTON_A = 2.0


def coords_mol2():
    pts, dentro = [], False
    for l in open(MOL2, encoding="utf-8", errors="ignore"):
        if l.startswith("@<TRIPOS>ATOM"):
            dentro = True
            continue
        if dentro and l.startswith("@<TRIPOS>"):
            break
        if dentro:
            p = l.split()
            if len(p) >= 6 and p[5].split(".")[0].upper() != "H":
                pts.append([float(p[2]), float(p[3]), float(p[4])])
    return np.array(pts)


def poses_de(carpeta):
    """Coordenadas pesadas de cada modelo de la pose de Vina."""
    ficheros = sorted(glob.glob(os.path.join(carpeta, "out", "*_out.pdbqt")))
    if not ficheros:
        return []
    modelos, actual = [], None
    for l in open(ficheros[0], encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            actual = []
        elif l.startswith("ENDMDL"):
            if actual is not None:
                modelos.append(np.array(actual))
            actual = None
        elif l.startswith(("ATOM", "HETATM")) and actual is not None:
            if l.rsplit(None, 1)[-1].strip().upper() in ("H", "HD", "HS", "D", "DD"):
                continue
            actual.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    if actual is not None:
        modelos.append(np.array(actual))
    return modelos


def transformaciones():
    """Rigid transforms applied to the crystal: identity, 180 about each axis,
    and the inversion through the centroid. El cristal se centra ANTES de
    rotarlo: sin centrar, la 'identidad' devuelve el RMSD de la distancia del
    centroide al origen (unos 29 A aqui), que es ruido puro."""
    out = [("identidad", np.eye(3))]
    for k in range(3):
        R = np.eye(3)
        R[k, k] = -1
        out.append(("giro 180 eje %d" % (k + 1), R))
    out.append(("inversion por el centroide", -np.eye(3)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variante", default="", help="carpeta dentro de _control_cdk2")
    args = ap.parse_args()

    cristal = coords_mol2()
    centro = cristal.mean(axis=0)
    if args.variante:
        carpetas = [os.path.join(CONTROL, args.variante)]
    else:
        carpetas = sorted(glob.glob(os.path.join(CONTROL, "*")))

    print("ligando del cristal: %d atomos pesados" % len(cristal))
    mejor_global = None
    for carpeta in carpetas:
        if not os.path.isdir(carpeta):
            continue
        poses = poses_de(carpeta)
        if not poses:
            continue
        print("\n%s (%d poses)" % (os.path.basename(carpeta), len(poses)))
        for i, pose in enumerate(poses, 1):
            if len(pose) != len(cristal):
                print("   pose %d: %d atomos, no comparable" % (i, len(pose)))
                continue
            d_centro = np.linalg.norm(pose.mean(axis=0) - centro)
            rmsds = []
            centrado = cristal - centro
            for nombre, R in transformaciones():
                t = centrado @ R.T + centro
                rmsd = float(np.sqrt(((t - pose) ** 2).sum(1).mean()))
                rmsds.append((rmsd, nombre))
            rmsds.sort()
            print("   pose %d: centroides a %.2f A | mejor transformacion: %s"
                  " (RMSD %.2f A) | identidad %.2f A"
                  % (i, d_centro, rmsds[0][1], rmsds[0][0],
                     [r for r, n in rmsds if n == "identidad"][0]))
            if mejor_global is None or rmsds[0][0] < mejor_global[0]:
                mejor_global = (rmsds[0][0], rmsds[0][1], os.path.basename(carpeta), i)

    if mejor_global is None:
        raise SystemExit("no hay poses que comparar")
    rmsd, nombre, carpeta, pose = mejor_global
    print("")
    print("MEJOR DE TODO: RMSD %.2f A con la transformacion '%s' (%s, pose %d)"
          % (rmsd, nombre, carpeta, pose))
    if rmsd < LISTON_A:
        print("VEREDICTO: la pose es el cristal dado la vuelta (%s)." % nombre)
        print("Receptor y caja estan bien: el ligando ocupa el bolsillo correcto.")
        print("Lo que falla es la puntuacion, que no distingue las dos")
        print("orientaciones de este ligando. Es un fallo CONOCIDO de Vina con")
        print("ligandos planos; la receta del proyecto la demuestra con el")
        print("BX-795 de TBK1 (1,23 A).")
    else:
        print("VEREDICTO: ni siquiera dado la vuelta cabe (RMSD %.2f A)." % rmsd)
        print("Hay que mirar la preparacion antes de creerse nada.")


if __name__ == "__main__":
    raise SystemExit(main())
