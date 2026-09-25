# -*- coding: utf-8 -*-
"""Descarga bioactividades de TDP-43 (CHEMBL2362981) y FUS (CHEMBL5724679)
desde ChEMBL con pChEMBL >= 5 y obtiene SMILES de cada molecula.
Guarda activos_tdp43_chembl.csv y activos_fus_chembl.csv (name,smiles,target)."""
import json
import time
import urllib.request
import urllib.parse

BASE = "https://www.ebi.ac.uk/chembl/api/data"


def get(url):
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.loads(r.read().decode("utf-8"))


def bioactividades(target_chembl_id, pchembl_min=5.0):
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
            "units": b.get("standard_units"),
            "relation": b.get("standard_relation"),
        })
    return out


def smiles_molecula(mol_id):
    url = "%s/molecule/%s.json" % (BASE, mol_id)
    try:
        d = get(url)
        return d.get("molecule_structures", {}).get("canonical_smiles")
    except Exception as e:
        print("  mol %s error: %s" % (mol_id, str(e)[:80]))
        return None


for tid, tgt, out in [("CHEMBL2362981", "TDP43", "activos_tdp43_chembl.csv"),
                      ("CHEMBL5724679", "FUS", "activos_fus_chembl.csv")]:
    print("=== %s (%s) ===" % (tgt, tid))
    try:
        acts = bioactividades(tid, 5.0)
    except Exception as e:
        print("  error bioactividades:", str(e)[:120])
        continue
    print("  bioactividades pChEMBL>=5:", len(acts))
    # deduplicar por molecula, quedarse con mejor pchembl
    mejor = {}
    for a in acts:
        mid = a["molecule_chembl_id"]
        if not mid:
            continue
        pc = a.get("pchembl")
        try:
            pc = float(pc) if pc else 0.0
        except Exception:
            pc = 0.0
        if mid not in mejor or pc > mejor[mid]["pchembl"]:
            mejor[mid] = a
    print("  moleculas unicas:", len(mejor))
    filas = []
    for mid, a in sorted(mejor.items(), key=lambda kv: -float(kv[1].get("pchembl") or 0)):
        smi = smiles_molecula(mid)
        if smi:
            filas.append((mid, smi, tgt, a.get("pchembl"), a.get("type")))
        time.sleep(0.3)
    print("  con SMILES:", len(filas))
    with open(out, "w", encoding="utf-8") as f:
        f.write("name,smiles,target,pchembl,type\n")
        for mid, smi, t, pc, ty in filas:
            f.write("%s,%s,%s,%s,%s\n" % (mid, smi, t, pc, ty or ""))
    print("  guardado:", out)
    time.sleep(1)
