#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Quita de rescoring_estratos.csv las 8 filas con error, para que el runner las repita.

El runner reanuda saltando los ligandos ya hechos, asi que para que vuelva a
calcular los 8 que fallaban hay que sacarlos de la tabla (y antes, copia de
seguridad del CSV).
"""
import csv
import os
import shutil
import time

OUT = "/home/ubuntu/mmgbsa/rescoring_estratos.csv"
FALLOS = ["DECM_CHEMBL978", "DECM_CHEMBL39736", "DECM_CHEMBL170988",
          "DEC_CHEMBL1414576", "DECM_CHEMBL1969867", "DECH_CHEMBL3309988",
          "DECM_CHEMBL20651", "DECM_CHEMBL4435214"]


def main():
    respaldo = OUT + ".bak_" + time.strftime("%Y%m%d_%H%M")
    shutil.copy2(OUT, respaldo)
    with open(OUT, newline="") as f:
        lector = csv.DictReader(f)
        campos = lector.fieldnames
        filas = list(lector)
    fuera = [r for r in filas if r["ligand"] in FALLOS]
    quedan = [r for r in filas if r["ligand"] not in FALLOS]
    tmp = OUT + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(quedan)
    os.replace(tmp, OUT)
    print(f"respaldo: {respaldo}")
    print(f"filas: {len(filas)} -> {len(quedan)} (fuera {len(fuera)})")
    for r in fuera:
        print(f"  fuera: {r['ligand']} (error: {r['error'][:60]})")


if __name__ == "__main__":
    main()
