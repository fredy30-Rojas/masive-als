#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""El banco de prueba de TBK1: los compuestos con afinidad medida de verdad.

POR QUE
-------
La auditoria del 26 de septiembre dejo claro que TDP-43 no se puede calibrar
porque de sus 40.125 actividades en ChEMBL **cero** son una medida de afinidad
(son un ensayo de crecimiento en levadura). TBK1, en cambio, tiene quimica
medida. Ese es todo el cambio.

EL RECUENTO DEL PLAN NO SE PUDO REPRODUCIR, Y ESTE ES EL MOTIVO
---------------------------------------------------------------
El plan (`PLAN_CALIBRACION_TBK1_2026-09-26.md`) dice 561 compuestos y una tabla
de reparto por potencia (82 / 117 / 232 / 110 / 20). Hoy, bajando la diana
entera, salen **2315**. La causa esta localizada y es un limite del script que
hizo el recuento:

  `analysis/buscar_diana_calibracion.py` tiene `def medidas(tid, tope=6)` y
  recorre `while offset // 1000 < tope`: **solo seis paginas**, o sea 6000 de las
  13123 actividades de la diana. Con esas primeras 6000 salen exactamente los
  1175 compuestos de la tabla del plan (comprobado hoy). La diana entera tiene
  2841 medidas de afinidad sobre 2315 compuestos.

Es decir: la tabla del plan es una lectura truncada, y el banco de verdad es
cuatro veces mayor y mejor (464 compuestos a 10 nanomolar o mejor, no 82). No se
cambia el plan por gusto: se cambia porque el numero que habia no era el numero.

QUE CUENTA COMO ACTIVO
----------------------
Solo medidas de afinidad de verdad (IC50, Ki, Kd, EC50) con valor pChEMBL,
contra la diana humana de proteina unica **TBK1 = CHEMBL5408**. De cada compuesto
se queda **el mejor valor** (mas pChEMBL), que es lo que hace
`verdad_de_referencia.py` para las otras dianas, para que las tablas sean
comparables.

Lo que NO entra: porcentajes de inhibicion, ensayos celulares sin valor de
afinidad, y las actividades sin pChEMBL.

QUE ESCRIBE
-----------
    compuestos_tbk1.csv        molecule_chembl_id, smiles, tipo, valor, unidades,
                               pchembl y franja de potencia
    actividades_tbk1.json      las 13123 actividades crudas (cache: el banco se
                               puede rehacer sin volver a pedir nada)
    banco_tbk1.txt             el informe legible: recuento, reparto y comparacion
                               con lo que decia el plan

Uso:
    python bajar_compuestos_tbk1.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
import time
import urllib.request

BASE = os.path.dirname(os.path.abspath(__file__))
CSV_SALIDA = os.path.join(BASE, "compuestos_tbk1.csv")
CACHE = os.path.join(BASE, "actividades_tbk1.json")
TXT_SALIDA = os.path.join(BASE, "banco_tbk1.txt")

API = "https://www.ebi.ac.uk/chembl/api/data/"
DIANA = "CHEMBL5408"          # TBK1 humana, proteina unica
TIPO_DIANA = "SINGLE PROTEIN"
UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) masive-als/banco-tbk1"

TIPOS = {"IC50", "Ki", "Kd", "EC50"}

FRANJAS = [
    (8.0, None, "10 nanomolar o mejor (pChEMBL >= 8)"),
    (7.0, 8.0, "de 10 a 100 nanomolar"),
    (6.0, 7.0, "de 0,1 a 1 micromolar"),
    (5.0, 6.0, "de 1 a 10 micromolar"),
    (None, 5.0, "peor que 10 micromolar"),
]

# Lo que decia el plan, para poder poner las dos columnas juntas.
PLAN = {"compuestos": 561,
        "reparto": {"10 nanomolar o mejor (pChEMBL >= 8)": 82,
                    "de 10 a 100 nanomolar": 117,
                    "de 0,1 a 1 micromolar": 232,
                    "de 1 a 10 micromolar": 110,
                    "peor que 10 micromolar": 20}}
PLAN_TABLA = {"actividades_diana": 13123, "con_afinidad": 1342, "compuestos": 1175}


def pide(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=90) as r:
        return json.load(r)


def franja(p):
    for lo, hi, nombre in FRANJAS:
        if (lo is None or p >= lo) and (hi is None or p < hi):
            return nombre
    return "sin franja"


def bajar_actividades():
    """Todas las actividades de la diana, o la cache si ya esta."""
    if os.path.exists(CACHE) and os.path.getsize(CACHE) > 1000:
        return json.load(open(CACHE, encoding="utf-8"))
    acts, offset = [], 0
    while True:
        d = pide(API + "activity?target_chembl_id=%s&format=json&limit=1000&offset=%d"
                 % (DIANA, offset))
        a = d.get("activities", [])
        if not a:
            break
        acts += a
        offset += 1000
        if offset >= d["page_meta"]["total_count"]:
            break
        time.sleep(0.2)
    with open(CACHE, "w", encoding="utf-8") as f:
        json.dump(acts, f)
    return acts


def mejor_por_compuesto(acts, limite=None):
    """{molecule_chembl_id: {pchembl, ...}} con el mejor valor de afinidad."""
    mejor = {}
    for a in (acts[:limite] if limite else acts):
        if a.get("standard_type") not in TIPOS:
            continue
        p = a.get("pchembl_value")
        if not p:
            continue
        try:
            p = float(p)
        except (TypeError, ValueError):
            continue
        mol = a.get("molecule_chembl_id")
        if not mol:
            continue
        prev = mejor.get(mol)
        if prev is None or p > prev["pchembl"]:
            v = a.get("standard_value")
            try:
                v = float(v) if v is not None else None
            except (TypeError, ValueError):
                v = None
            mejor[mol] = {"molecule_chembl_id": mol,
                          "smiles": (a.get("canonical_smiles") or "").strip(),
                          "tipo": a.get("standard_type"), "valor": v,
                          "unidades": a.get("standard_units"), "pchembl": p}
    return mejor


def main():
    # --- la diana, comprobada ---
    t = pide(API + "target/%s?format=json" % DIANA)
    if t.get("organism") != "Homo sapiens" or t.get("target_type") != TIPO_DIANA:
        raise SystemExit("CHEMBL5408 no es una diana humana de proteina unica: %s"
                         % t.get("target_type"))
    print("diana: %s (%s, %s)" % (t.get("pref_name"), DIANA, t.get("organism")))

    acts = bajar_actividades()
    print("actividades colgando de la diana: %d" % len(acts))

    # --- el reparto que decia el plan, reproducido con las primeras 6000 ---
    parcial = mejor_por_compuesto(acts, limite=6000)

    # --- el banco completo ---
    mejor = mejor_por_compuesto(acts)
    sin_smiles = [m for m, v in mejor.items() if not v["smiles"]]
    for mol in sin_smiles:
        try:
            m = pide(API + "molecule/%s?format=json" % mol)
            s = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if s:
                mejor[mol]["smiles"] = s.strip()
            time.sleep(0.15)
        except Exception:  # noqa: BLE001
            pass
    sin_smiles = [m for m, v in mejor.items() if not v["smiles"]]
    if sin_smiles:
        print("   sin SMILES tras reintentar: %d (se dejan fuera del banco)"
              % len(sin_smiles))

    filas = [v for v in mejor.values() if v["smiles"]]
    for v in filas:
        v["franja"] = franja(v["pchembl"])
    filas.sort(key=lambda v: -v["pchembl"])

    campos = ["molecule_chembl_id", "smiles", "tipo", "valor", "unidades",
              "pchembl", "franja"]
    with open(CSV_SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(filas)

    reparto = {}
    for v in filas:
        reparto[v["franja"]] = reparto.get(v["franja"], 0) + 1

    # --- informe, con las dos columnas juntas ---
    L = []
    L.append("BANCO DE PRUEBA DE TBK1 (%s) — %s" % (DIANA, time.strftime("%Y-%m-%d")))
    L.append("")
    L.append("actividades colgando de la diana: %d" % len(acts))
    L.append("actividades de afinidad de verdad (IC50/Ki/Kd/EC50 con pChEMBL): %d"
             % sum(1 for a in acts if a.get("standard_type") in TIPOS
                   and a.get("pchembl_value")))
    L.append("compuestos del banco (el mejor valor de cada uno): %d" % len(filas))
    L.append("")
    L.append("reparto por potencia, y lo que decia el plan:")
    L.append("   %-38s %8s %8s" % ("franja", "hoy", "plan"))
    for lo, hi, nombre in FRANJAS:
        L.append("   %-38s %8d %8d"
                 % (nombre, reparto.get(nombre, 0), PLAN["reparto"][nombre]))
    L.append("")
    L.append("el numero del plan, explicado:")
    L.append("   buscar_diana_calibracion.py solo pedia 6 paginas (6000 actividades)")
    L.append("   de las %d de la diana. Con esas primeras 6000 salen %d compuestos,"
             % (len(acts), len(parcial)))
    L.append("   que es el 1175 de la tabla del plan. La diana entera trae %d."
             % len(filas))
    L.append("")
    L.append("los mas potentes:")
    for v in filas[:10]:
        L.append("   %-16s %-5s %-10s pChEMBL %.2f  %s"
                 % (v["molecule_chembl_id"], v["tipo"],
                    ("%g %s" % (v["valor"], v["unidades"])) if v["valor"] else "-",
                    v["pchembl"], v["franja"]))
    L.append("")
    L.append("CSV: %s" % os.path.basename(CSV_SALIDA))
    L.append("cache de actividades: %s" % os.path.basename(CACHE))
    texto = "\n".join(L)
    with open(TXT_SALIDA, "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    print(texto)
    return 0


if __name__ == "__main__":
    sys.exit(main())
