# -*- coding: utf-8 -*-
"""Mide coste y dispersion del rescoring con receptor CONGELADO.

Toma un ligando por diana de la lista corta y lo reescala N veces con el
receptor fijo (--fijo) en la plataforma indicada, midiendo tiempo y dG.
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import time

BASE = "/mnt/c/Users/Fredy/masive-als"
LOCAL = os.path.join(BASE, "analysis", "rescoring_local")
HIJO = os.path.join(LOCAL, "mmgbsa_openff_gb.py")
FIJO = os.path.join(LOCAL, "receptores_fijos")
POSES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
LISTA = os.path.join(BASE, "analysis", "lista_corta_candidatos.csv")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=3)
    ap.add_argument("--platform", default="CUDA")
    ap.add_argument("--targets", default="SOD1,TDP43_v2,FUS")
    a = ap.parse_args()

    objetivos = [t.strip() for t in a.targets.split(",") if t.strip()]
    elegidos = {}
    # Dos ligandos por diana: el primero y otro al azar determinista (el 30).
    with open(LISTA, encoding="utf-8") as f:
        filas = list(csv.DictReader(f))
    for t in objetivos:
        del_t = [r for r in filas if r["target"] == t]
        elegidos[t] = [del_t[0], del_t[min(29, len(del_t) - 1)]]

    resumen = []
    env = dict(os.environ, MMGBSA_FINAL_PLATFORM=a.platform)
    for t, ligs in elegidos.items():
        rec = os.path.join(FIJO, "%s_fijo.pdb" % t)
        for r in ligs:
            lig, smi = r["ligand"], r["smiles"]
            pose = os.path.join(POSES, "results_%s" % t, lig + "_out.pdbqt")
            if not os.path.exists(pose):
                print("SIN POSE %s %s" % (t, lig), flush=True)
                continue
            dgs, ts = [], []
            for i in range(a.repeats):
                cmd = [sys.executable, HIJO, "--pose", pose, "--receptor", rec,
                       "--fijo", "--smiles", smi, "--out", "json"]
                t0 = time.time()
                p = subprocess.run(cmd, capture_output=True, text=True, env=env)
                dt = time.time() - t0
                try:
                    row = json.loads(p.stdout.strip().splitlines()[-1])
                except Exception:
                    row = {"error": (p.stderr or "")[-200:]}
                if row.get("mmgbsa_dG") is not None:
                    dgs.append(row["mmgbsa_dG"])
                ts.append(dt)
                print("  %-9s %-14s dG=%s  %.1f s"
                      % (t, lig, row.get("mmgbsa_dG", row.get("error")), dt),
                      flush=True)
            if dgs:
                disp = max(dgs) - min(dgs)
                resumen.append((t, lig, dgs[0], disp, sum(ts) / len(ts)))

    print("\n=== RESUMEN ===")
    print("%-9s %-14s %8s %10s %8s" % ("diana", "ligando", "dG", "disp", "t(s)"))
    for t, lig, dg, disp, tt in resumen:
        print("%-9s %-14s %8.2f %10.2f %8.1f" % (t, lig, dg, disp, tt))
    if resumen:
        print("\nt medio: %.1f s | dispersion media: %.2f kcal/mol"
              % (sum(x[4] for x in resumen) / len(resumen),
                 sum(x[3] for x in resumen) / len(resumen)))


if __name__ == "__main__":
    main()
