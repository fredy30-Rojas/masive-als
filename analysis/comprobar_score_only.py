#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""comprobar_score_only.py — ¿es comparable la afinidad de `--score_only` con la del docking?

POR QUE ESTA COMPROBACION
------------------------
`repunuar_vinardo.py` no re-acopla: toma las poses ya guardadas y las vuelve a puntuar
con `--score_only`. Para que las metricas de `metricas_vinardo.py` con la columna `vina`
sean comparables con las que publico `validar_diana_limpia.py`, haria falta que
`--score_only` devolviera EXACTAMENTE la misma afinidad que la que Vina escribio en el
`_out.pdbqt` al acoplar. No tiene por que ser asi: al acoplar, Vina hace una
optimizacion local (BFGS) de la pose y luego puntua la pose ya optimizada; con
`--score_only` se puntua la pose tal cual esta guardada.

Si la diferencia es notable, la columna `vina` del re-puntado es "Vina sobre la pose
guardada", que es una cifra un poco distinta de "Vina en el docking", y las dos no se
pueden mezclar en la misma tabla sin decirlo. Este script mide la diferencia en vez de
suponer que es cero.

Uso:
    python comprobar_score_only.py --target TDP43
    python comprobar_score_only.py --target SOD1
"""
import argparse
import csv
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="TDP43", choices=["TDP43", "SOD1"])
    args = ap.parse_args()

    import validar_sod1_limpia as V
    import validar_diana_limpia as D
    cfg = D.DIANAS[args.target]

    csv_p = os.path.join(BASE, "vinardo_%s.csv" % args.target.lower())
    con_p = {r["ligand"]: r for r in csv.DictReader(open(csv_p, encoding="utf-8"))}

    poses2 = cfg.get("poses2") or os.path.join(os.path.dirname(V.LIGS_FONDO2), "out")
    carpetas = [("fondo blando", cfg["poses"]), ("fondo duro", poses2),
                ("positivos", os.path.join(cfg["salida"], "out"))]

    print("DIANA %s — afinidad del docking frente a la de --score_only" % args.target)
    print("%-12s %6s %10s %10s %10s %10s"
          % ("carpeta", "n", "mediana", "maximo", "p95", "mayor de 0,05"))
    for etiqueta, d in carpetas:
        if not d or not os.path.isdir(d):
            continue
        filas = []
        for p in sorted(os.listdir(d)):
            if not p.endswith("_out.pdbqt"):
                continue
            nombre = p.replace("_out.pdbqt", "")
            if nombre.startswith("ACT_"):
                nombre = nombre[4:]
            fila = con_p.get(nombre)
            if not fila or not fila.get("vina"):
                continue
            aff = V.afinidad(os.path.join(d, p))
            if aff is None:
                continue
            filas.append(abs(float(fila["vina"]) - aff))
        if not filas:
            continue
        filas.sort()
        n = len(filas)
        print("%-12s %6d %10.4f %10.4f %10.4f %10d"
              % (etiqueta, n, filas[n // 2], filas[-1], filas[int(n * 0.95)],
                 sum(1 for x in filas if x > 0.05)))
    print("")
    print("Una diferencia grande NO es un fallo del script: es que al acoplar Vina")
    print("optimiza la pose antes de puntuarla y `--score_only` no lo hace. Por eso la")
    print("columna `vina` del re-puntado es la misma funcion sobre la misma pose, pero")
    print("no necesariamente el mismo numero que el docking, y las dos columnas no se")
    print("mezclan en la misma tabla.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
