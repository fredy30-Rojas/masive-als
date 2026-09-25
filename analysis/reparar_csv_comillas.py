#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""reparar_csv_comillas.py — arregla CSVs con comas SIN comillas en un campo de texto.

EL PROBLEMA
-----------
`controles_sod1_v3.csv` se escribio a mano (o con un writer que no cito el campo) y
una de sus columnas, `redocking_pasa`, lleva texto libre con comas:

    ... ,no (y su pose depositada no es un mínimo del potencial: F···N de 2,07 Å),
    "co-cristalización; ribosa al disolvente (sitio MEDIDO)",Wright et al. 2013 ...

Esa coma de mas parte la fila en dos y **desplaza todas las columnas siguientes una
posicion**: al leerlo con `csv.DictReader`, `tipo_ensayo` recibe la cola del texto y
`cita` recibe el tipo de ensayo. El fichero parece bueno y no lo es — que es la peor
clase de fichero malo.

COMO SE REPARA
--------------
Se lee cada fila en crudo, y si trae mas campos que la cabecera se **vuelven a unir
los sobrantes en la columna de texto libre** (`--campo`, por defecto la ultima
columna que puede tener comas). Se reescribe el fichero con `csv.writer`, que cita
lo que haga falta, y se guarda copia del original en `.comillas_rotas`.

Antes de escribir se comprueba que despues TODAS las filas tienen el mismo numero de
campos que la cabecera. Si alguna no cuadra, no se toca nada.

Uso:
  python reparar_csv_comillas.py controles_sod1_v3.csv
  python reparar_csv_comillas.py fichero.csv --campo redocking_pasa --seco
"""
import argparse
import csv
import os
import shutil
import sys

BASE = os.path.dirname(os.path.abspath(__file__))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ficheros", nargs="+")
    ap.add_argument("--campo", default=None,
                    help="nombre de la columna de texto libre. Si no se dice, se "
                         "reconstruye la fila dejando el sobrante en la ultima "
                         "columna conocida antes del desfase")
    ap.add_argument("--seco", action="store_true", help="solo informar, no escribir")
    args = ap.parse_args()

    for f in args.ficheros:
        ruta = f if os.path.isabs(f) else os.path.join(BASE, f)
        filas = list(csv.reader(open(ruta, encoding="utf-8")))
        if not filas:
            print("%s: vacio" % f)
            continue
        cab = filas[0]
        n = len(cab)
        cuerpos = filas[1:]
        desfasadas = [i for i, r in enumerate(cuerpos) if len(r) != n]
        print("%s: %d columnas, %d filas, %d con desfase"
              % (os.path.basename(ruta), n, len(cuerpos), len(desfasadas)))
        if not desfasadas:
            print("   nada que reparar")
            continue

        # la columna donde va el texto libre: la que se indica, o la que precede
        # siempre al primer desplazamiento observado
        if args.campo and args.campo in cab:
            idx = cab.index(args.campo)
        else:
            # el desfase empieza en la primera fila mala: se toma la columna cuyo
            # contenido en las filas buenas parece texto libre (espacios y signos)
            idx = n - 3
            for i, r in enumerate(cuerpos):
                if len(r) != n:
                    idx = min(len(r), n) - 3
                    break
            print("   (no se dijo --campo: se usa la columna %d = '%s')"
                  % (idx, cab[idx]))

        fuera = []
        for r in cuerpos:
            if len(r) == n:
                fuera.append(r)
                continue
            sobra = len(r) - n
            fila = r[:idx] + [",".join(r[idx:idx + 1 + sobra])] + r[idx + 1 + sobra:]
            fuera.append(fila)

        malas = [i for i, r in enumerate(fuera) if len(r) != n]
        if malas:
            print("   NO se toca: despues de unir siguen %d filas con desfase"
                  % len(malas))
            continue

        for i in desfasadas[:5]:
            print("   fila %d: %s -> %s" % (i + 2, cab[idx], fuera[i][idx][:80]))

        if args.seco:
            print("   (--seco: no se ha escrito nada)")
            continue

        shutil.copy(ruta, ruta + ".comillas_rotas")
        with open(ruta, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(cab)
            for r in fuera:
                w.writerow(r)
        print("   reparado y guardado (original en %s.comillas_rotas)"
              % os.path.basename(ruta))

        # comprobacion final: se lee otra vez como lo leeria cualquier lector
        prueba = list(csv.DictReader(open(ruta, encoding="utf-8")))
        anchos = {len(r) for r in prueba}
        print("   releido: %d filas, %d columnas, campos vacios: %d"
              % (len(prueba), anchos.pop() if anchos else 0,
                 sum(1 for r in prueba for v in r.values() if v is None)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
