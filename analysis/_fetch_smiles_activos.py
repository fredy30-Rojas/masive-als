# -*- coding: utf-8 -*-
"""Descarga SMILES de los activos conocidos de TDP-43 y FUS (literatura,
revision Francois-Moutal 2021 y Frontiers 2025) via PubChem y guarda
activos_tdp43.csv / activos_fus.csv."""
import json
import time
import urllib.request

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

ACTIVOS = {
    "TDP43": [
        ("rTRD01", 56739767),        # RRM1/RRM2, desplaza RNA (HSQC-NMR)
        ("nTRD22", 45189578),        # NTD, alosterico sobre union RNA
        ("bis-ANS", 123808),         # CTD, modula LLPS
        ("CongoRed", 11313),         # CTD, modula LLPS
        ("5FUrd", 9427),             # RRM (analogo de nucleosido)
        ("isoproterenol", 3779),     # compuesto conocido del panel SOD1/TDP43
    ],
    "FUS": [
        ("Dehydroxymethylflazine", 5488120),   # inhibitor natural FUS (docking+MD)
        ("CleroindicinC", 10975728),           # inhibitor natural FUS (docking+MD)
    ],
}


def get_json(url):
    req = urllib.request.Request(url, headers={
        "Accept": "application/json", "User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.loads(r.read().decode("utf-8"))


def smiles_de_cid(cid):
    url = "%s/compound/cid/%d/property/ConnectivitySMILES/JSON" % (BASE, cid)
    d = get_json(url)
    return d["PropertyTable"]["Properties"][0].get("ConnectivitySMILES", "")


for target, lista in ACTIVOS.items():
    out = "activos_%s.csv" % ("tdp43" if target == "TDP43" else "fus")
    with open(out, "w", encoding="utf-8") as f:
        f.write("name,smiles,target\n")
        for nom, cid in lista:
            try:
                smi = smiles_de_cid(cid)
                print("%s | %s | %s" % (target, nom, smi))
                f.write("%s,%s,%s\n" % (nom, smi, target))
            except Exception as e:
                print("%s | %s | ERROR %s" % (target, nom, e))
            time.sleep(0.6)
    print("guardado:", out)
