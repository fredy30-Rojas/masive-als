# -*- coding: utf-8 -*-
"""Test decisivo: enriquecimiento. Donde caen los ligandos conocidos (controles)
en el ranking de Vina y en el de MM-GBSA local?"""
import csv
import os
from collections import defaultdict

BASE = r"C:\Users\Fredy\masive-als"
CTRL = os.path.join(BASE, "analysis", "controles_calibracion.csv")
RES = os.path.join(BASE, "analysis", "rescoring_local", "rescoring_lista_corta.csv")


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


# controles en libreria
ctrl = []
with open(CTRL, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        if (r.get("en_libreria") or "").strip().lower() == "si":
            ctrl.append(r)

print("Controles positivos presentes en la libreria: %d" % len(ctrl))
for c in ctrl:
    print("   %-8s %-22s lig=%-14s vina=%s"
          % (c["target"], c["control"], c.get("ligand_libreria"),
             c.get("afinidad_vina")))

# rankings
rows = []
with open(RES, encoding="utf-8") as f:
    for r in csv.DictReader(f):
        rows.append(r)

by_t = defaultdict(list)
for r in rows:
    by_t[r["target"]].append(r)

print()
print("=" * 74)
print("POSICION DE LOS CONTROLES EN EL RANKING (1 = mejor)")
print("=" * 74)
for c in ctrl:
    t = (c["target"] or "").replace("_v2", "")
    lig = (c.get("ligand_libreria") or "").strip()
    sub = by_t.get(t, [])
    if not sub:
        # probar variante TDP43_v2
        sub = by_t.get(t + "_v2", [])
    dgs = [(fnum(x["mmgbsa_dG"]), x["ligand"]) for x in sub if fnum(x["mmgbsa_dG"]) is not None]
    vinas = [(fnum(x["vina_affinity"]), x["ligand"]) for x in sub if fnum(x["vina_affinity"]) is not None]
    if not dgs:
        continue
    dgs.sort()
    vinas.sort()
    n = len(dgs)
    pos_d = next((i + 1 for i, (d, l) in enumerate(dgs) if l == lig), None)
    pos_v = next((i + 1 for i, (d, l) in enumerate(vinas) if l == lig), None)
    pct_d = ("%.0f%%" % (100.0 * pos_d / n)) if pos_d else "ausente"
    pct_v = ("%.0f%%" % (100.0 * pos_v / len(vinas))) if pos_v else "ausente"
    print("  %-8s %-20s (%s): MMGBSA %s/%d (%s)   VINA %s/%d (%s)"
          % (c["target"], c["control"], lig, pos_d, n, pct_d, pos_v, len(vinas), pct_v))

print()
print("Si el rescoring fuera valido, los controles (activos conocidos) deberian")
print("caer en el TOP 5-10% de sus listas.")
