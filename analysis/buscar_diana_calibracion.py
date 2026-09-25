#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""buscar_diana_calibracion.py — ¿que diana de ELA tiene quimica conocida?

POR QUE
-------
El proyecto no puede validar su embudo contra TDP-43 porque no hay ligandos
potentes publicados contra el bolsillo de RRM2 (ver
AUDITORIA_ACTIVOS_TDP43_2026-09-26.md). Sin referencia, cualquier AUC es ruido.

La salida no es acoplar mas compuestos: es calibrar el instrumento en una diana
donde SI se pueda saber la respuesta, y despues aplicar el embudo ya medido a
TDP-43. Este script mira, para los genes que causan ELA, cuanta afinidad medida
hay publicada en ChEMBL. Cuanta mas, mejor se puede calibrar.

QUE CUENTA
----------
Solo la diana humana de proteina unica, y solo medidas de afinidad de verdad
(IC50, Ki, Kd, EC50, con valor pChEMBL), no ensayos funcionales ni celulares.

Uso:
    python buscar_diana_calibracion.py
"""
import json
import time
import urllib.parse
import urllib.request

UA = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) masive-als/diana-calibracion"
API = "https://www.ebi.ac.uk/chembl/api/data/"

GENES = [
    "SOD1", "VCP", "TBK1", "FUS", "HDAC6", "OPTN", "C9orf72", "NEK1",
    "KIF5A", "ATXN2", "STMN2", "SQSTM1", "UBQLN2", "PFN1", "ANXA11",
    "TARDBP",
]
TIPOS = {"IC50", "Ki", "Kd", "EC50", "AC50", "Kd(app)"}


def pide(url):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def diana_humana(gen):
    q = urllib.parse.quote(gen)
    d = pide(API + "target/search?q=%s&format=json&limit=50" % q)
    for t in d.get("targets", []):
        if (t.get("organism") == "Homo sapiens"
                and t.get("target_type") == "SINGLE PROTEIN"):
            return t
    return None


def medidas(tid, tope=6):
    """Cuenta actividades de afinidad, paginando un maximo de `tope` paginas."""
    total = 0
    compuestos = set()
    offset = 0
    while offset // 1000 < tope:
        d = pide(API + "activity?target_chembl_id=%s&format=json&limit=1000&offset=%d"
                 % (tid, offset))
        acts = d.get("activities", [])
        if not acts:
            break
        for a in acts:
            if a.get("standard_type") in TIPOS and a.get("pchembl_value"):
                total += 1
                compuestos.add(a.get("molecule_chembl_id"))
        offset += 1000
        if offset >= d["page_meta"]["total_count"]:
            break
    return total, len(compuestos)


print("%-10s %-16s %10s %10s %10s   %s"
      % ("gen", "CHEMBL", "actividades", "afinidad", "compuestos", "nombre"))
print("-" * 100)
filas = []
for gen in GENES:
    try:
        t = diana_humana(gen)
        if not t:
            print("%-10s %s" % (gen, "sin diana humana de proteina unica"))
            continue
        tid = t["target_chembl_id"]
        tot = pide(API + "activity?target_chembl_id=%s&format=json&limit=1" % tid)
        n_tot = tot["page_meta"]["total_count"]
        n_af, n_comp = medidas(tid)
        filas.append((n_af, n_comp, gen, tid, t["pref_name"]))
        print("%-10s %-16s %10d %10d %10d   %s"
              % (gen, tid, n_tot, n_af, n_comp, t["pref_name"]))
    except Exception as e:
        print("%-10s FALLO: %s" % (gen, e))
    time.sleep(0.3)

print()
print("ORDENADO POR COMPUESTOS CON AFINIDAD MEDIDA")
for n_af, n_comp, gen, tid, nombre in sorted(filas, reverse=True):
    print("   %-10s %-16s %5d compuestos, %5d medidas   %s"
          % (gen, tid, n_comp, n_af, nombre))
