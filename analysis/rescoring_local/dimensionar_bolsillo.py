# -*- coding: utf-8 -*-
"""Dimensiona el bolsillo del receptor a partir de las poses reales acopladas.

Para cada diana mide, sobre las poses de la lista corta, la distancia del
atomo de ligando mas lejano al CENTRO DE LA CAJA de docking. El radio del
bolsillo se fija por encima del percentil 99 de esa distancia, con margen:
asi todos los compuestos se evaluan contra un bolsillo que los contiene
completos, y no se incluyen atomos que no aportan nada al dG.
"""
import argparse
import csv
import json
import os
import sys

import numpy as np

BASE = "/mnt/c/Users/Fredy/masive-als"
if not os.path.isdir(BASE):
    BASE = r"C:\Users\Fredy\masive-als"
LOCAL = os.path.join(BASE, "analysis", "rescoring_local")
POSES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
LISTA = os.path.join(BASE, "analysis", "lista_corta_candidatos.csv")

sys.path.insert(0, LOCAL)
import mmgbsa_openff as M  # noqa: E402

CAJAS = {
    "SOD1":     {"centro": (46.5, 80.0, 73.3), "size": 22},
    "TDP43_v2": {"centro": (24.23, 16.89, -15.87), "size": 26},
    "FUS":      {"centro": (-14.5, 15.1, -7.8), "size": 25},
    "TDP43":    {"centro": (16.3, 41.1, 48.5), "size": 24},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--muestra", type=int, default=None,
                    help="limita el numero de poses por diana")
    a = ap.parse_args()

    filas = {}
    with open(LISTA, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            filas.setdefault(r["target"], []).append(r["ligand"])

    salida = {}
    for t, ligs in filas.items():
        if t not in CAJAS:
            continue
        c = np.array(CAJAS[t]["centro"], dtype=float)
        if a.muestra:
            ligs = ligs[:a.muestra]
        dmax, errores = [], 0
        for lig in ligs:
            p = os.path.join(POSES, "results_%s" % t, lig + "_out.pdbqt")
            if not os.path.exists(p):
                errores += 1
                continue
            try:
                at = M.parse_pdbqt_atoms(p, primer_model=True)
            except Exception:
                errores += 1
                continue
            if not at:
                errores += 1
                continue
            xyz = np.array([[x["x"], x["y"], x["z"]] for x in at])
            dmax.append(float(np.linalg.norm(xyz - c, axis=1).max()))
        if not dmax:
            salida[t] = {"error": "sin poses legibles", "n": 0}
            continue
        d = np.array(dmax)
        info = {"n_poses": len(d), "size_caja": CAJAS[t]["size"],
                "mitad": CAJAS[t]["size"] / 2.0, "errores": errores,
                "max": round(float(d.max()), 2),
                "p50": round(float(np.percentile(d, 50)), 2),
                "p90": round(float(np.percentile(d, 90)), 2),
                "p99": round(float(np.percentile(d, 99)), 2),
                "p999": round(float(np.percentile(d, 99.9)), 2)}
        info["radio_recomendado"] = round(info["p999"] + 5.0, 1)
        salida[t] = info
        print("%-9s n=%d  max=%.1f  p50=%.1f  p90=%.1f  p99=%.1f  p99.9=%.1f "
              "-> radio %.1f A (caja mitad=%.1f)"
              % (t, info["n_poses"], info["max"], info["p50"], info["p90"],
                 info["p99"], info["p999"], info["radio_recomendado"],
                 info["mitad"]), flush=True)

    rep = os.path.join(LOCAL, "bolsillo_dimensionado.json")
    with open(rep, "w", encoding="utf-8") as f:
        json.dump(salida, f, ensure_ascii=False, indent=2)
    print("reporte: %s" % rep)


if __name__ == "__main__":
    main()
