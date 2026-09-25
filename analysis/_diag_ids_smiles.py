# -*- coding: utf-8 -*-
"""Diagnostico de emparejamiento ligand->SMILES entre rescoring y candidatos."""
import re
import pandas as pd
from pathlib import Path

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

maps = {}
for f in ["candidatos_limpios.csv", "candidatos_mmgbsa.csv"]:
    d = pd.read_csv(DATA / f)
    d.columns = [str(c).strip().lower() for c in d.columns]
    if "smiles" in d.columns:
        idc = next((c for c in ["ligand", "molecule_chembl_id", "chembl_id", "compound_id", "id"] if c in d.columns), None)
        if idc:
            for k, s in zip(d[idc].astype(str).str.strip(), d["smiles"].astype(str).str.strip()):
                maps.setdefault(k, s)

found = valid & set(maps)
missing = sorted(valid - set(maps))
print("rescatados validos:", len(valid), "| con SMILES exacto:", len(found), "| sin match:", len(missing))
print("ejemplos sin match:", missing[:8])

def norm(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())

nvalid = {norm(x) for x in valid}
nmaps = {norm(x) for x in maps}
inter = nvalid & nmaps
print("con normalizacion alfanumerica:", len(inter))
extra = sorted({x for x in missing if norm(x) in inter})[:8]
print("sin match que encajan al normalizar:", extra)
# colisiones: un normalizado con varios SMILES
from collections import defaultdict
g = defaultdict(set)
for k, s in maps.items():
    g[norm(k)].add(s)
col = {k: v for k, v in g.items() if len(v) > 1 and k in nvalid}
print("colisiones en normalizados:", len(col))
for k in list(col)[:3]:
    print("  ", k, "->", list(col[k])[:3])
