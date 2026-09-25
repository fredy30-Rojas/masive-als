#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Con la semilla puesta, ¿da el mismo dG dos veces la misma fila?

Lanza dos veces la misma fila en paralelo y compara los tres terminos. Si los
numeros coinciden, el metodo es repetible y la tabla se puede citar.
"""
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

LIG = "ACT_isoproterenol"
REPS = 3


def una(rep):
    import csv
    idx = 0
    with open("estratos/candidatos.csv", newline="") as f:
        for i, r in enumerate(csv.DictReader(f)):
            if r["ligand"] == LIG:
                idx = i
                break
    cmd = [sys.executable, "rescoring_mmgbsa_robusto.py",
           "--candidates", "estratos/candidatos.csv",
           "--receptores", "receptores_estratos.json", "--row", str(idx)]
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    lineas = [l for l in p.stdout.splitlines() if l.startswith("{")]
    if not lineas:
        return rep, None, (p.stderr or "")[-200:], dt
    return rep, json.loads(lineas[-1]), None, dt


def main():
    with ThreadPoolExecutor(max_workers=REPS) as ex:
        res = list(ex.map(una, range(1, REPS + 1)))
    for rep, d, err, dt in res:
        if err:
            print(f"rep{rep}: ERROR {err}", flush=True)
            continue
        print(f"rep{rep}: dG={d['mmgbsa_dG']} E_comp={d['e_complex']} "
              f"E_rec={d['e_receptor']} E_lig={d['e_ligand']} ({dt:.0f} s)",
              flush=True)
    ds = [d["mmgbsa_dG"] for _, d, e, _ in res if d]
    if len(ds) == REPS:
        print(f"\nspread de dG: {min(ds)} a {max(ds)} "
              f"(rango {max(ds) - min(ds):.2f} kcal/mol)", flush=True)
        print("REPETIBLE" if max(ds) - min(ds) == 0 else "NO REPETIBLE", flush=True)


if __name__ == "__main__":
    main()
