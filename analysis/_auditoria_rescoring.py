# -*- coding: utf-8 -*-
"""Auditoria cientifica del rescoring MM-GBSA (lista corta).

Comprueba si los dG son fisicamente utilizables:
  1. cobertura y errores por target
  2. consistencia interna dG = e_complex - e_receptor - e_ligand
  3. sanciones fisicas: e_ligand > 0 (imposible), saltos entre corridas
  4. eficiencia de ligando (|dG| / atomos pesados) vs rango realista
  5. correlacion Vina vs MM-GBSA (Spearman) por target
  6. reproducibilidad: mismo ligando en corridas/targets distintos
"""
import csv
import math
import os
import statistics as st
from collections import defaultdict

BASE = r"C:\Users\Fredy\masive-als"
CSV = os.path.join(BASE, "analysis", "rescoring_local", "rescoring_lista_corta.csv")
V5 = os.path.join(BASE, "rescoring_datos", "rescoring_mmgbsa_v5.csv")
T100 = os.path.join(BASE, "rescoring_datos", "rescoring_tdp43_mmgbsa.csv")
SRC = os.path.join(BASE, "analysis", "lista_corta_candidatos.csv")


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def heavy_atoms(smi):
    if not smi:
        return None
    try:
        from rdkit import Chem
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        return m.GetNumHeavyAtoms()
    except Exception:
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


# ---------- cargar ----------
rows = []
with open(CSV, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        rows.append(r)

src = {}
with open(SRC, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        src[(r["target"], r["ligand"])] = r

print("=" * 78)
print("1) COBERTURA Y ERRORES")
print("=" * 78)
by_t = defaultdict(lambda: {"ok": 0, "err": 0, "dgs": []})
for r in rows:
    d = fnum(r.get("mmgbsa_dG"))
    k = r.get("target")
    if d is None:
        by_t[k]["err"] += 1
    else:
        by_t[k]["ok"] += 1
        by_t[k]["dgs"].append(d)
tot_ok = sum(v["ok"] for v in by_t.values())
tot_err = sum(v["err"] for v in by_t.values())
for k in sorted(by_t):
    v = by_t[k]
    d = v["dgs"]
    rng = "%.2f (%.2f..%.2f)" % (st.mean(d), min(d), max(d)) if d else "-"
    print("  %-10s ok=%-5d err=%-5d  media=%s" % (k, v["ok"], v["err"], rng))
print("  TOTAL ok=%d err=%d  (%.1f%% exito)" % (tot_ok, tot_err,
                                                100.0 * tot_ok / (tot_ok + tot_err)))

print()
print("=" * 78)
print("2) CONSISTENCIA INTERNA  dG == e_complex - e_receptor - e_ligand")
print("=" * 78)
bad_internal = 0
checked = 0
for r in rows:
    d, ec, er, el = (fnum(r.get("mmgbsa_dG")), fnum(r.get("e_complex")),
                     fnum(r.get("e_receptor")), fnum(r.get("e_ligand")))
    if None in (d, ec, er, el):
        continue
    checked += 1
    if abs((ec - er - el) - d) > 0.05:
        bad_internal += 1
print("  filas con las 4 energias: %d | inconsistentes: %d" % (checked, bad_internal))

print()
print("=" * 78)
print("3) SANCIONES FISICAS")
print("=" * 78)
elig = [fnum(r.get("e_ligand")) for r in rows]
elig = [v for v in elig if v is not None]
pos = [v for v in elig if v > 0]
print("  e_ligand disponible: %d | e_ligand > 0 (imposible): %d (%.1f%%)"
      % (len(elig), len(pos), 100.0 * len(pos) / len(elig)))
if elig:
    print("  e_ligand: min=%.1f  mediana=%.1f  max=%.1f" % (min(elig), st.median(elig), max(elig)))
erec = [fnum(r.get("e_receptor")) for r in rows]
erec = [v for v in erec if v is not None]
if erec:
    print("  e_receptor: min=%.1f  mediana=%.1f  max=%.1f" % (min(erec), st.median(erec), max(erec)))
    # bimodalidad: cuantos receptores grandes vs pequenos
    big = [v for v in erec if v < -3000]
    print("  e_receptor < -3000 (receptor grande): %d ; >= -3000 (recorte): %d"
          % (len(big), len(erec) - len(big)))

dgs = [fnum(r.get("mmgbsa_dG")) for r in rows]
dgs = [v for v in dgs if v is not None]
dgs_s = sorted(dgs)
print("  dG: min=%.2f  P25=%.2f  mediana=%.2f  P75=%.2f  max=%.2f  desv=%.2f"
      % (min(dgs), dgs_s[len(dgs_s) // 4], st.median(dgs),
         dgs_s[3 * len(dgs_s) // 4], max(dgs), st.pstdev(dgs)))

print()
print("=" * 78)
print("4) EFICIENCIA DE LIGANDO  |dG| / atomos pesados")
print("=" * 78)
eff = []
no_smiles = 0
for r in rows:
    d = fnum(r.get("mmgbsa_dG"))
    key = (r.get("target"), r.get("ligand"))
    smi = r.get("smiles") or (src.get(key, {}) or {}).get("smiles")
    if not smi:
        no_smiles += 1
    if d is None or not smi:
        continue
    ha = heavy_atoms(smi)
    if not ha:
        continue
    eff.append((abs(d) / ha, r.get("ligand"), ha, d))
print("  con SMILES: %d | sin SMILES: %d" % (len(rows) - no_smiles, no_smiles))
if eff:
    vals = sorted(e[0] for e in eff)
    print("  LE: min=%.3f  mediana=%.3f  max=%.3f  kcal/mol/atomo"
          % (vals[0], st.median(vals), vals[-1]))
    for lim in (1.5, 2.0, 2.5, 3.0):
        n = sum(1 for v in vals if v > lim)
        print("     LE > %.1f : %d (%.1f%%)" % (lim, n, 100.0 * n / len(vals)))
    print("  Top-5 mas extremos:")
    for v, lig, ha, d in sorted(eff, reverse=True)[:5]:
        print("     %-14s %2d atomos  dG=%.2f  LE=%.2f" % (lig, ha, d, v))

print()
print("=" * 78)
print("5) CORRELACION VINA vs MM-GBSA (Spearman) POR TARGET")
print("=" * 78)
for t in sorted(by_t):
    xs, ys = [], []
    for r in rows:
        if r.get("target") != t:
            continue
        a, d = fnum(r.get("vina_affinity")), fnum(r.get("mmgbsa_dG"))
        if a is not None and d is not None:
            xs.append(a)
            ys.append(d)
    print("  %-10s n=%-5d rho=%s" % (t, len(xs),
                                     ("%.3f" % spearman(xs, ys)) if spearman(xs, ys) is not None else "-"))

print()
print("=" * 78)
print("6) REPRODUCIBILIDAD ENTRE CORRIDAS")
print("=" * 78)
old = {}
for path, tag in ((V5, "v5"), (T100, "tdp43_real")):
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                d = fnum(r.get("mmgbsa_dG"))
                if d is not None:
                    old.setdefault((r.get("ligand") or "").strip(), []).append((tag, r.get("target"), d))
newmap = {}
for r in rows:
    d = fnum(r.get("mmgbsa_dG"))
    if d is not None:
        newmap.setdefault((r.get("ligand") or "").strip(), []).append(
            (r.get("target"), d))
common = sorted(set(old) & set(newmap))
print("  ligandos en corridas antiguas: %d | en corrida nueva: %d | comunes: %d"
      % (len(old), len(newmap), len(common)))
diffs = []
for lig in common:
    for tag, tgt_old, d_old in old[lig]:
        for tgt_new, d_new in newmap[lig]:
            if tgt_old and tgt_new and tgt_old.replace("_v2", "") == tgt_new.replace("_v2", ""):
                diffs.append((lig, tgt_old, d_old, d_new, d_new - d_old))
if diffs:
    dd = [x[4] for x in diffs]
    print("  pares comparables: %d | delta absoluto medio=%.2f  max=%.2f kcal/mol"
          % (len(dd), st.mean([abs(x) for x in dd]), max(abs(x) for x in dd)))
    for lig, tgt, a, b, d in sorted(diffs, key=lambda x: -abs(x[4]))[:8]:
        print("     %-14s %-9s viejo=%.2f nuevo=%.2f  delta=%+.2f" % (lig, tgt, a, b, d))

print()
print("=" * 78)
print("7) MISMO LIGANDO EN VARIOS TARGETS (dependencia del receptor)")
print("=" * 78)
sizes = [(len(v), k) for k, v in newmap.items() if len(v) > 1]
print("  ligandos con dG en >1 target: %d" % len(sizes))
spread = []
for lig, vals in newmap.items():
    if len(vals) > 1:
        ds = [v[1] for v in vals]
        spread.append(max(ds) - min(ds))
if spread:
    print("  dispersion entre targets: mediana=%.2f  max=%.2f kcal/mol"
          % (st.median(spread), max(spread)))
