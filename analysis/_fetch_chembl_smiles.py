# -*- coding: utf-8 -*-
"""Descarga SMILES canonicos de ChEMBL para los ligandos del rescoring sin SMILES local."""
import json
import time
import urllib.request
from pathlib import Path
import pandas as pd

DATA = Path(r"C:\Users\Fredy\masive-als\rescoring_datos")
FILES = ["rescoring_mmgbsa_v5.csv", "rescoring_nuevo.csv", "rescoring_tdp43_top100.csv",
         "rescoring_tdp43_top100_v2.csv", "rescoring_tdp43_mmgbsa.csv", "rescoring_mmgbsa.csv"]

valid = set()
for f in FILES:
    d = pd.read_csv(DATA / f)
    d.columns = [str(c).strip() for c in d.columns]
    if "mmgbsa_dG" not in d.columns or "ligand" not in d.columns:
        continue
    v = pd.to_numeric(d["mmgbsa_dG"], errors="coerce")
    valid |= set(d.loc[v.notna(), "ligand"].astype(str).str.strip())

have = set()
for f in ["candidatos_limpios.csv", "candidatos_mmgbsa.csv"]:
    d = pd.read_csv(DATA / f)
    d.columns = [str(c).strip().lower() for c in d.columns]
    if "smiles" in d.columns:
        idc = next((c for c in ["ligand", "molecule_chembl_id", "chembl_id", "compound_id", "id"] if c in d.columns), None)
        if idc:
            have |= set(d[idc].astype(str).str.strip())

missing = sorted(valid - have)
print(f"ligandos sin SMILES local: {len(missing)}")
out = []
errs = []
for i, cid in enumerate(missing):
    url = f"https://www.ebi.ac.uk/chembl/api/data/molecule/{cid}.json"
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "masive-als-rescoring/1.0"})
        with urllib.request.urlopen(req, timeout=20) as r:
            js = json.load(r)
        smi = (js.get("molecule_structures") or {}).get("canonical_smiles")
        if smi:
            out.append({"ligand": cid, "smiles": smi})
        else:
            errs.append((cid, "sin smiles"))
    except Exception as e:
        errs.append((cid, str(e)[:60]))
    time.sleep(0.5)

if out:
    pd.DataFrame(out).to_csv(DATA / "chembl_smiles_extra.csv", index=False)
print(f"obtenidos: {len(out)} | fallos: {len(errs)}")
for cid, e in errs[:8]:
    print("  ", cid, "->", e)
