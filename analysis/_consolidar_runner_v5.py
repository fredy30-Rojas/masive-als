# -*- coding: utf-8 -*-
"""Consolida rescoring_mmgbsa_v5_raw.csv -> ranking MM-GBSA definitivo.

Deduplica por (ligand, target) quedándose con la mejor fila válida
(dG mas negativo, ignorando dG=0.0 degenerados y errores).
Genera ranking_mmgbsa_v5.csv ordenado por dG por target y el ranking
consolidado multi-target.
"""
import csv, json, sys, os

SRC = os.path.join(os.path.dirname(__file__), "rescoring_mmgbsa_v5_raw.csv")
OUT = os.path.join(os.path.dirname(__file__), "ranking_mmgbsa_v5.csv")

rows = list(csv.DictReader(open(SRC, encoding="utf-8-sig")))

def parse_dg(r):
    try:
        return float(r["mmgbsa_dG"])
    except (TypeError, ValueError):
        return None

def is_bad(r):
    dg = parse_dg(r)
    err = (r.get("error") or "").strip()
    if err:
        return True
    if dg is None:
        return True
    # dG exactamente 0.0 => complejo == partes separadas (degenerado, no confiable)
    if abs(dg) < 1e-6:
        return True
    return False

# deduplicar: mejor fila valida por (ligand,target)
best = {}
for r in rows:
    key = (r["ligand"], r["target"])
    cur = best.get(key)
    if cur is None:
        best[key] = r
        continue
    # preferir valida sobre invalida
    if is_bad(cur) and not is_bad(r):
        best[key] = r
    elif (not is_bad(cur)) and (not is_bad(r)):
        if (parse_dg(r) or 0) < (parse_dg(cur) or 0):
            best[key] = r

pares = sorted(best.values(), key=lambda r: (r["target"], parse_dg(r) or 0))

# escribir ranking por target
with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["ligand", "target", "vina_affinity", "mmgbsa_dG", "status"])
    w.writeheader()
    for r in pares:
        dg = parse_dg(r)
        status = "error" if is_bad(r) else "ok"
        w.writerow({
            "ligand": r["ligand"], "target": r["target"],
            "vina_affinity": r["vina_affinity"],
            "mmgbsa_dG": "" if dg is None else ("%.2f" % dg),
            "status": status,
        })

print("pares unicos: %d (de %d filas crudas)" % (len(pares), len(rows)))
print("ok: %d | degenerados/error: %d" % (
    sum(1 for r in pares if not is_bad(r)),
    sum(1 for r in pares if is_bad(r)),
))
print()
# top por target
for tgt in ["SOD1", "TDP43", "FUS"]:
    sub = [r for r in pares if r["target"] == tgt and not is_bad(r)]
    sub.sort(key=lambda r: parse_dg(r))
    print("=== %s (%d ok) ===" % (tgt, len(sub)))
    for r in sub[:10]:
        print("  %-16s %8.2f" % (r["ligand"], parse_dg(r)))
    print()
