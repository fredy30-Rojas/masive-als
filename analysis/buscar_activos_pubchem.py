# -*- coding: utf-8 -*-
"""Controles positivos para TDP-43 (TARDBP) y FUS desde PubChem bioassays.
ChEMBL API esta caida (500); PubChem PUG REST responde OK.

1. Para cada gen, obtiene los AIDs de ensayos y los CIDs activos.
2. Descarga SMILES de esos CIDs.
3. Guarda activos_tdp43.csv y activos_fus.csv (name,smiles,target).
"""
import json
import time
import urllib.request
import urllib.parse

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def get(url, accept="application/json"):
    req = urllib.request.Request(url, headers={
        "Accept": accept, "User-Agent": "Mozilla/5.0 (MASIVE-ALS research)"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return r.read()


def cids_activos_por_gen(gen):
    """Devuelve lista de CIDs activos en ensayos del gen."""
    url = "%s/assay/target/genesymbol/%s/cids/JSON" % (BASE, gen)
    d = json.loads(get(url))
    cids = []
    for info in d.get("InformationList", {}).get("Information", []):
        cids.extend(info.get("CID", []))
    return list(dict.fromkeys(cids))  # dedup preservando orden


def smiles_de_cids(cids, lote=50):
    """SMILES canonicos para una lista de CIDs (peticiones por lotes)."""
    mapa = {}
    for i in range(0, len(cids), lote):
        lote_ids = cids[i:i + lote]
        url = ("%s/compound/cid/%s/property/CanonicalSMILES,IsomericSMILES/JSON"
               % (BASE, ",".join(str(c) for c in lote_ids)))
        try:
            raw = get(url)
            d = json.loads(raw)
            for prop in d.get("PropertyTable", {}).get("Properties", []):
                cid = prop.get("CID")
                smi = prop.get("CanonicalSMILES") or prop.get("IsomericSMILES")
                if cid and smi:
                    mapa[str(cid)] = smi
        except Exception as e:
            print("  lote %d error: %s" % (i, e))
        time.sleep(0.5)
    return mapa


for gen, out, tgt in [("TARDBP", "activos_tdp43.csv", "TDP43"),
                      ("FUS", "activos_fus.csv", "FUS")]:
    print("=== %s ===" % gen)
    try:
        cids = cids_activos_por_gen(gen)
    except Exception as e:
        print("  error obteniendo cids:", e)
        continue
    print("  CIDs activos:", len(cids))
    mapa = smiles_de_cids(cids)
    print("  con SMILES:", len(mapa))
    with open(out, "w", encoding="utf-8") as f:
        f.write("name,smiles,target\n")
        for cid in cids:
            smi = mapa.get(str(cid))
            if smi:
                f.write("CID%s,%s,%s\n" % (cid, smi, tgt))
    print("  guardado:", out)
