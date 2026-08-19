# -*- coding: utf-8 -*-
"""Busca activos conocidos de TDP-43 y FUS en ChEMBL (controles positivos
para la validacion de senuelos del embudo — revision cientifica 18/08).

Usa el endpoint bioactivities del target CHEMBL correspondiente.
TDP-43 = TARDBP, FUS = FUS (proteinas de union a RNA).
"""
import json
import time
import urllib.request
import urllib.parse

BASE = "https://www.ebi.ac.uk/chembl/api/data"


def get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return json.loads(r.read().decode("utf-8"))


def buscar_target(query):
    url = "%s/target/search.json?q=%s" % (BASE, urllib.parse.quote(query))
    d = get(url)
    return d.get("targets", [])[:5]


def bioactividades(target_chembl_id, pchembl_min=6.0):
    """Actividades con pChEMBL >= umbral (potencia reportada)."""
    out = []
    url = ("%s/bioactivity.json?target_chembl_id=%s&pchembl_value__gte=%.1f"
           "&limit=1000" % (BASE, target_chembl_id, pchembl_min))
    d = get(url)
    for b in d.get("bioactivities", []):
        out.append({
            "molecule_chembl_id": b.get("molecule_chembl_id"),
            "pchembl": b.get("pchembl_value"),
            "type": b.get("standard_type"),
            "value": b.get("standard_value"),
            "relation": b.get("standard_relation"),
        })
    return out


for query in ["TARDBP", "TDP-43", "FUS"]:
    print("=== target search: %s ===" % query)
    try:
        for t in buscar_target(query):
            print("  %s | %s | %s" % (t.get("target_chembl_id"),
                                      t.get("pref_name"),
                                      t.get("organism")))
    except Exception as e:
        print("  error:", e)
    time.sleep(1)
