# -*- coding: utf-8 -*-
"""Diagnostico de las 237 filas con mmgbsa_dG=0 descartadas por entrenar_modelo_rescoring.py."""
import pandas as pd
from pathlib import Path

DATA = Path(r"C:\Users\Fredy\masive-als\rescoring_datos")
FILES = ["rescoring_mmgbsa_v5.csv", "rescoring_nuevo.csv", "rescoring_tdp43_top100.csv",
         "rescoring_tdp43_top100_v2.csv", "rescoring_tdp43_mmgbsa.csv", "rescoring_mmgbsa.csv"]

def num(s):
    try:
        x = float(str(s).strip().replace(",", "."))
        return x if abs(x) != float("inf") else None
    except (TypeError, ValueError):
        return None

for f in FILES:
    p = DATA / f
    if not p.exists():
        print(f, "NO EXISTE")
        continue
    d = pd.read_csv(p)
    d.columns = [str(c).strip() for c in d.columns]
    dg = next((c for c in ("mmgbsa_dG", "dG", "mmgbsa") if c in d.columns), None)
    err = next((c for c in ("error", "status", "fail") if c in d.columns), None)
    print("=" * 25, f, "| filas:", len(d))
    print("  columnas:", list(d.columns))
    if dg is None:
        print("  sin columna de energia")
        continue
    vals = d[dg].map(num)
    ceros = vals.abs() < 1e-9
    nulos = vals.isna()
    print(f"  dG: ceros={int(ceros.sum())}, nulos={int(nulos.sum())}, validos>0={int((~ceros & ~nulos).sum())}")
    if err:
        print("  valores columna", err, ":", d[err].value_counts(dropna=False).head(8).to_dict())
    if int(ceros.sum()) and ceros.any():
        zr = d[ceros]
        muestra = zr.head(2)
        print("  muestra filas cero:")
        for _, r in muestra.iterrows():
            print("   ", {k: str(r[k])[:60] for k in r.index})
        # los mismos ligandos con dG valido en otros archivos
        lig_cero = set(zr.get("ligand", pd.Series(dtype=str)).astype(str))
        lig_ok = set()
        for f2 in FILES:
            p2 = DATA / f2
            if not p2.exists():
                continue
            d2 = pd.read_csv(p2)
            d2.columns = [str(c).strip() for c in d2.columns]
            if dg in d2.columns and "ligand" in d2.columns:
                v2 = d2[dg].map(num)
                lig_ok |= set(d2.loc[v2.abs() > 1e-9, "ligand"].astype(str))
        inter = lig_cero & lig_ok
        print(f"  ligandos con cero aqui pero con dG valido en otro archivo: {len(inter)} de {len(lig_cero)}")
        if inter:
            print("  ejemplos:", sorted(inter)[:5])
