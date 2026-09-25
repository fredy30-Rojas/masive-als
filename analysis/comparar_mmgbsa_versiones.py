#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""comparar_mmgbsa_versiones.py — cuanto se movieron los dG al arreglar el metodo.

Lee dos CSV de resultado del MM-GBSA (el viejo, mixto v4/v5, y el nuevo v6 con
la preparacion determinista y el recorte centrado en la primera pose) y saca:

  * cuantos ligandos cambian y cuanto;
  * los que mas cambian, con su papel;
  * si el orden se conserva (correlacion de puestos) y cuantos positivos se
    mueven de mitad de tabla.

Uso:
  python comparar_mmgbsa_versiones.py viejo.csv nuevo.csv
"""
import csv
import sys

import numpy as np


def leer(ruta):
    d = {}
    for r in csv.DictReader(open(ruta, encoding="utf-8")):
        try:
            d[r["ligand"]] = (float(r["mmgbsa_dG"]), r.get("papel", ""))
        except (KeyError, TypeError, ValueError):
            continue
    return d


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 1
    viejo, nuevo = leer(sys.argv[1]), leer(sys.argv[2])
    comunes = sorted(set(viejo) & set(nuevo))
    print(f"viejos {len(viejo)} | nuevos {len(nuevo)} | comunes {len(comunes)}")
    if not comunes:
        return 1

    dv = np.array([viejo[k][0] for k in comunes])
    dn = np.array([nuevo[k][0] for k in comunes])
    dif = dn - dv
    print(f"\ndiferencia (nuevo - viejo): media {dif.mean():+.2f} "
          f"| mediana {np.median(dif):+.2f} | max {dif.max():+.2f} "
          f"| min {dif.min():+.2f}")
    print(f"correlacion entre versiones: {np.corrcoef(dv, dn)[0, 1]:.3f}")
    for umbral in (0.5, 1.0, 2.0, 4.0):
        print(f"  cambian mas de {umbral:.1f} kcal/mol: "
              f"{int((np.abs(dif) > umbral).sum())} de {len(dif)}")

    print("\nlos 10 que mas cambian:")
    for k in sorted(comunes, key=lambda k: -abs(nuevo[k][0] - viejo[k][0]))[:10]:
        print(f"  {k:34s} {viejo[k][0]:8.2f} -> {nuevo[k][0]:8.2f} "
              f"({nuevo[k][0] - viejo[k][0]:+7.2f})  {nuevo[k][1]}")

    # puestos: se conserva el orden?
    pv = {k: i for i, k in enumerate(sorted(comunes, key=lambda k: viejo[k][0]))}
    pn = {k: i for i, k in enumerate(sorted(comunes, key=lambda k: nuevo[k][0]))}
    print(f"\ncorrelacion de puestos: "
          f"{np.corrcoef([pv[k] for k in comunes], [pn[k] for k in comunes])[0, 1]:.3f}")
    pos = [k for k in comunes if nuevo[k][1] == "positivo"]
    if pos:
        print("\npuestos de los positivos (1 = mejor dG):")
        for k in sorted(pos, key=lambda k: pn[k]):
            print(f"  {k:34s} puesto {pn[k] + 1:3d} de {len(comunes)} "
                  f"(antes {pv[k] + 1:3d})  {nuevo[k][0]:8.2f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
