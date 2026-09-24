#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""separar_fondos_sod1.py — parte el fondo de SOD1 en sus dos bloques.

POR QUE
-------
La validacion limpia de SOD1 puntua **un solo fondo de 482** señuelos, pero ese fondo
son tres cosas pegadas: 199 emparejados antiguos (`DEC_`), 140 emparejados nuevos
(`DECM_`) y **143 duros** —quelantes de metales y redox-activos, generados con los
SMARTS de `validar_sod1_v3.py`— (`DECH_`). La regla nueva
(`REGLA_DECISION_2026-09-24.md`) decide sobre el fondo DURO, que es el unico que puede
separar quimia especifica de quimia generica; y en TDP-43 el bloque mezclado fue
justamente el que se cayo por dos milesimas. Para leer SOD1 con la misma vara hay que
separarlos.

El prefijo no es una suposicion: lo pone el propio preparador
(`clave_del_fondo()` en `validar_sod1_v3.py` -> `DECH_%s`, `DECM_` para los emparejados
nuevos, `DEC_` para los antiguos).

No recalcula nada: copia `validar_sod1_limpia.csv` y solo cambia la columna `papel`
(`fondo2` para los `DECH_`, `fondo` para el resto de señuelos). Los positivos se quedan
como estan.

Uso:
    python separar_fondos_sod1.py
Salida: regla_decision/sod1_fondos_separados.csv
"""
import csv
import os

BASE = os.path.dirname(os.path.abspath(__file__))
ENTRADA = os.path.join(BASE, "validar_sod1_limpia.csv")
SALIDA = os.path.join(BASE, "regla_decision", "sod1_fondos_separados.csv")

DURO = "DECH_"          # quelantes de metales y redox-activos
EMPAREJADOS = ("DEC_", "DECM_")


def main():
    with open(ENTRADA, encoding="utf-8") as f:
        filas = list(csv.DictReader(f))
    campos = list(filas[0])
    cuenta = {"fondo": 0, "fondo2": 0, "sin_tocar": 0}
    for r in filas:
        n = r["ligand"]
        if n.startswith(DURO):
            r["papel"] = "fondo2"
            cuenta["fondo2"] += 1
        elif n.startswith(EMPAREJADOS):
            r["papel"] = "fondo"
            cuenta["fondo"] += 1
        else:
            cuenta["sin_tocar"] += 1
    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(filas)
    print("escrito %s" % SALIDA)
    print("   emparejados (DEC_/DECM_) %d | duros (DECH_) %d | %d sin tocar (positivos)"
          % (cuenta["fondo"], cuenta["fondo2"], cuenta["sin_tocar"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
