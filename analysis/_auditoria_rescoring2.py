# -*- coding: utf-8 -*-
"""Pruebas adicionales: el dG depende del tamano/carga del ligando? (artefacto de vacio)"""
import csv
import math
import statistics as st
from collections import defaultdict

CSV = r"C:\Users\Fredy\masive-als\analysis\rescoring_local\rescoring_lista_corta.csv"


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def spearman(xs, ys):
    n = len(xs)
    if n < 5:
        return None

    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r

    rx, ry = ranks(xs), ranks(ys)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else None


from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors

rows = []
with open(CSV, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if fnum(r.get("mmgbsa_dG")) is not None and r.get("smiles"):
            rows.append(r)

data = []
for r in rows:
    m = Chem.MolFromSmiles(r["smiles"])
    if m is None:
        continue
    d = fnum(r["mmgbsa_dG"])
    q = Chem.GetFormalCharge(m)
    data.append({
        "t": r["target"], "lig": r["ligand"], "dG": d,
        "ha": m.GetNumHeavyAtoms(), "q": q,
        "tpsa": Descriptors.TPSA(m), "logp": Descriptors.MolLogP(m),
        "rot": rdMolDescriptors.CalcNumRotatableBonds(m),
    })

print("n =", len(data))

print("\n--- Correlacion de dG con propiedades del LIGANDO (por target) ---")
for t in sorted(set(d["t"] for d in data)):
    sub = [d for d in data if d["t"] == t]
    print("  %-9s n=%d  rho(dG,atomos)=%s  rho(dG,TPSA)=%s  rho(dG,logP)=%s" % (
        t, len(sub),
        "%.3f" % spearman([d["dG"] for d in sub], [d["ha"] for d in sub]),
        "%.3f" % spearman([d["dG"] for d in sub], [d["tpsa"] for d in sub]),
        "%.3f" % spearman([d["dG"] for d in sub], [d["logp"] for d in sub])))

print("\n--- dG medio por carga formal del ligando ---")
byq = defaultdict(list)
for d in data:
    byq[d["q"]].append(d["dG"])
for q in sorted(byq):
    v = byq[q]
    print("  carga %+d : n=%-5d dG medio=%7.2f" % (q, len(v), st.mean(v)))

print("\n--- dG medio por carga, separando el tipo de soluto ---")
print("  (si el artefacto fuera de vacio, los cationes/aniones serian extremos)")

print("\n--- Tamano del sistema por target (offsets que no cancelan) ---")
by_er = defaultdict(list)
for r in rows:
    er = fnum(r.get("e_receptor"))
    if er is not None:
        by_er[r["target"]].append(er)
for t in sorted(by_er):
    v = by_er[t]
    print("  %-9s e_receptor: min=%9.1f mediana=%9.1f max=%9.1f  (n=%d)"
          % (t, min(v), st.median(v), max(v), len(v)))
