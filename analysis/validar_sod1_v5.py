#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validacion de SOD1 v5: la v4 repetida contra el receptor LIMPIO (21 sep 2026).

POR QUE EXISTE
--------------
La v4 midio el AUC de SOD1 contra `gpu_dock/SOD1.pdbqt`, y ese fichero resulto
no ser una proteina limpia: 21.585 atomos, 1.763 aguas cristalograficas, las 18
copias de 1HL5 y los metales. Las aguas no se pueden desplazar y las 18 copias
abren grietas de empaquetamiento cristalino: ahi se fueron las poses de -10
kcal/mol de las quinazolinas, que son las que daban el «positivo» del proyecto
(ver redocking_trp32/INFORME_REDOCKING_QUIMIAS_NUEVAS_2026-09-21.md).

Aqui se repite TODO el calculo de la v4 —los 14 controles y los 483 señuelos de
los tres fondos— contra `gpu_dock/SOD1_limpio.pdbqt` (un solo dimer biologico,
cadenas A y H, sin aguas ni metales; ver construir_receptor_sod1.py), con la
misma caja, la misma exhaustividad y el mismo pipeline. Despues se comparan las
dos tablas lado a lado.

No reimplementa nada: importa `validar_sod1_v4` y le cambia las rutas, para que
la comparacion sea de receptor y no de codigo.

Uso:
    python validar_sod1_v5.py --preparar
    python validar_sod1_v5.py --acoplar [--workers 16]
    python validar_sod1_v5.py --analizar
    python validar_sod1_v5.py --comparar
    python validar_sod1_v5.py            # todo menos comparar
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import shutil
import sys

import validar_sod1_v4 as V4
import validar_senuelos as VS

BASE = V4.BASE
GPU = os.path.abspath(os.path.join(BASE, "..", "gpu_dock"))

# --- se le cambian las rutas a la maquinaria de la v4 -------------------------
V4.RECEPTOR = os.path.join(GPU, "SOD1_limpio.pdbqt")
V4.WORK = os.path.join(BASE, "validacion_SOD1_v5")
V4.LIGDIR = os.path.join(V4.WORK, "ligands")
V4.OUTDIR = os.path.join(V4.WORK, "out")
# Los señuelos de la v3 se RE-ACOPLAN contra el receptor limpio: ahora el fondo
# sale de la carpeta de esta version, no de la de la v3.
V4.V3_OUT = V4.OUTDIR
V4.ANTIGUO = os.path.join(V4.WORK, "antiguos_v5.csv")

FUENTES = [
    os.path.join(BASE, "validacion_SOD1_v4", "ligands"),
    os.path.join(BASE, "validacion_SOD1_v3", "ligands"),
    os.path.join(BASE, "_validacion_SOD1", "ligands"),
]


def log(m):
    print(m, flush=True)


def preparar():
    """Reune en una sola carpeta todos los ligandos de los tres conjuntos.

    El orden importa: el primero que llega manda, y el primero es el de la v4
    (los controles nuevos, preparados con el pipeline actual).
    """
    os.makedirs(V4.LIGDIR, exist_ok=True)
    total, nuevos = 0, 0
    for fuente in FUENTES:
        for p in sorted(glob.glob(os.path.join(fuente, "*.pdbqt"))):
            dst = os.path.join(V4.LIGDIR, os.path.basename(p))
            total += 1
            if os.path.exists(dst):
                continue
            shutil.copy(p, dst)
            nuevos += 1
    log("ligandos del conjunto: %d (copiados %d) en %s"
        % (len(glob.glob(os.path.join(V4.LIGDIR, "*.pdbqt"))), nuevos, V4.LIGDIR))


def escribir_antiguos():
    """CSV con las afinidades NUEVAS de los señuelos DEC_ (los 199 antiguos).

    La v4 leia de `_validacion_SOD1/validacion_SOD1_trp32.csv`, que son las
    afinidades medidas contra el receptor VIEJO. Aqui hay que leer las de esta
    version, o se estaria mezclando receptores en la misma tabla.
    """
    filas = []
    for p in sorted(glob.glob(os.path.join(V4.OUTDIR, "DEC_*_out.pdbqt"))):
        n = os.path.basename(p).replace("_out.pdbqt", "")
        a = VS.parse_affinity(p)
        if a is not None:
            filas.append({"ligand": n, "rol": "decoy", "affinity": a})
    with open(V4.ANTIGUO, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ligand", "rol", "affinity"])
        w.writeheader()
        w.writerows(filas)
    log("señuelos antiguos re-acoplados y anotados: %d -> %s" % (len(filas), V4.ANTIGUO))
    return len(filas)


def comparar():
    """Tabla lado a lado v4 (receptor sucio) frente a v5 (receptor limpio)."""
    def leer(ruta):
        if not os.path.exists(ruta):
            return {}
        with open(ruta, encoding="utf-8") as f:
            return {r["ligando"]: r for r in csv.DictReader(f)}

    a = leer(os.path.join(BASE, "validacion_SOD1_v4", "detalle_controles_sod1_v4.csv"))
    b = leer(os.path.join(BASE, "validacion_SOD1_v5", "detalle_controles_sod1_v4.csv"))
    log("")
    log("%-30s %-9s %-9s %-8s %-12s %-12s"
        % ("control", "v4", "v5", "cambio", "puesto v4", "puesto v5"))
    for nombre in a:
        ra, rb = a.get(nombre), b.get(nombre, {})
        fa, fb = ra.get("afinidad", ""), rb.get("afinidad", "")
        try:
            cambio = "%+.2f" % (float(fb) - float(fa))
        except (TypeError, ValueError):
            cambio = "n/d"
        log("%-30s %-9s %-9s %-8s %-12s %-12s"
            % (nombre, fa, fb or "n/d", cambio,
               ra.get("puesto_frente_a_los_tres_fondos", ""),
               rb.get("puesto_frente_a_los_tres_fondos", "n/d")))

    log("")
    log("AUC por familia:")
    log("%-34s %-19s %-7s %-7s %-7s %s" % ("familia", "fondo", "AUC v4", "AUC v5",
                                           "sin tam v4", "sin tam v5"))
    fa = _auc(os.path.join(BASE, "validacion_SOD1_v4", "analisis_sod1_v4.csv"))
    fb = _auc(os.path.join(BASE, "validacion_SOD1_v5", "analisis_sod1_v4.csv"))
    for clave in fa:
        x, y = fa[clave], fb.get(clave, {})
        log("%-34s %-19s %-7s %-7s %-7s %s"
            % (clave[0], clave[1], x.get("AUC"), y.get("AUC", "n/d"),
               x.get("AUC_sin_tamano"), y.get("AUC_sin_tamano", "n/d")))
    return 0


def _auc(ruta):
    if not os.path.exists(ruta):
        return {}
    with open(ruta, encoding="utf-8") as f:
        return {(r["familia"], r["fondo"]): r for r in csv.DictReader(f)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preparar", action="store_true")
    ap.add_argument("--acoplar", action="store_true")
    ap.add_argument("--analizar", action="store_true")
    ap.add_argument("--comparar", action="store_true")
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()
    os.makedirs(V4.WORK, exist_ok=True)
    os.makedirs(V4.OUTDIR, exist_ok=True)
    todo = not any([args.preparar, args.acoplar, args.analizar, args.comparar])

    if not os.path.exists(V4.RECEPTOR):
        log("FALTA el receptor limpio %s -> corre construir_receptor_sod1.py" % V4.RECEPTOR)
        return 1
    log("receptor: %s" % V4.RECEPTOR)

    if args.preparar or todo:
        log("=== preparando (reuniendo los tres conjuntos) ===")
        preparar()
    if args.acoplar or todo:
        log("=== acoplando contra el receptor limpio ===")
        V4.acoplar(workers=args.workers)
        escribir_antiguos()
    if args.analizar or todo:
        log("=== analizando ===")
        V4.analizar()
    if args.comparar:
        comparar()
    return 0


if __name__ == "__main__":
    sys.exit(main())
