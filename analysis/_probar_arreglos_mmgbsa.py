#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prueba los arreglos del MM-GBSA en los 8 ligandos que fallaron.

Por cada ligando lanza el script parcheado en modo --row y saca el dG o el
error. Con el ligando de control (que ya salia bien) lo lanza tambien con el
script viejo, para comprobar que el arreglo no cambia el resultado de los que
ya funcionaban.
"""
import csv
import json
import subprocess
import sys
import time

CASOS = ["DEC_CHEMBL1414576", "DECH_CHEMBL3309988", "DECM_CHEMBL978"]
CONTROL = "ACT_adrenalina"

NUEVO = "rescoring_mmgbsa_robusto.py"
VIEJO = "rescoring_mmgbsa_robusto.py.bak_20260922"


def fila_de(lig):
    with open("estratos/candidatos.csv", encoding="utf-8-sig") as f:
        for i, r in enumerate(csv.DictReader(f)):
            if r["ligand"] == lig:
                return i
    raise SystemExit(f"{lig} no está en el CSV")


def correr(script, idx):
    cmd = [sys.executable, script, "--candidates", "estratos/candidatos.csv",
           "--receptores", "receptores_estratos.json", "--row", str(idx)]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    salida = [l for l in (r.stdout or "").splitlines() if l.startswith("{")]
    if not salida:
        return None, (r.stderr or "")[-300:], dt
    return json.loads(salida[-1]), None, dt


def tarea(nombre, lig, script):
    out, err, dt = correr(script, fila_de(lig))
    if out is None:
        return f"{lig:22s} {nombre:5s} ERROR de script: {err}"
    if out.get("error"):
        return f"{lig:22s} {nombre:5s} error: {str(out['error'])[:80]}"
    return (f"{lig:22s} {nombre:5s} dG={out['mmgbsa_dG']:8.2f}  ({dt:.0f} s)"
            f"\t{out['mmgbsa_dG']}")


def main():
    from concurrent.futures import ThreadPoolExecutor

    trabajos = [("NUEVO", lig, NUEVO) for lig in CASOS + [CONTROL]]
    trabajos.append(("VIEJO", CONTROL, VIEJO))
    with ThreadPoolExecutor(max_workers=3) as ex:
        for res in ex.map(lambda t: tarea(*t), trabajos):
            print(res, flush=True)


if __name__ == "__main__":
    main()
