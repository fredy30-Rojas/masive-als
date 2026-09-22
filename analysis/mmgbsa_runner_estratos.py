#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""mmgbsa_runner_estratos.py — lanza el MM-GBSA por estratos de tamaño en Oracle.

Se ejecuta EN ORACLE (necesita openmm + openmmforcefields + openff + obabel +
antechamber, que aqui no estan y en el PC de Fredy tampoco).

POR QUE PARALELO
----------------
El runner que ya existia en el servidor (`runner_sod1.py`) va de uno en uno y tarda
~83 s por ligando medidos en el piloto. Con 90 ligandos son 2 horas y pico de reloj
para una prueba. La maquina tiene 4 nucleos, asi que aqui se lanzan 3 a la vez
(dejando uno libre: en esa VM corren tambien la memoria semantica y el indexado).

QUE HACE
--------
  * lee `estratos/candidatos.csv` y `receptores_estratos.json`;
  * reanuda: si `rescoring_estratos.csv` ya tiene un ligando, no lo repite;
  * por cada pendiente lanza `rescoring_mmgbsa_robusto.py --row N` y guarda el
    JSON que imprime;
  * los errores no tumban la tanda: se apuntan en la columna `error` y se sigue.

Uso (en Oracle):
  source /home/ubuntu/miniforge3/etc/profile.d/conda.sh && conda activate mmgbsa
  python /home/ubuntu/mmgbsa/mmgbsa_runner_estratos.py            # todos
  python /home/ubuntu/mmgbsa/mmgbsa_runner_estratos.py --workers 3
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

DIR = "/home/ubuntu/mmgbsa"
CAND = os.path.join(DIR, "estratos", "candidatos.csv")
RECEP = os.path.join(DIR, "receptores_estratos.json")
OUT = os.path.join(DIR, "rescoring_estratos.csv")
LOG = os.path.join(DIR, "runner_estratos.log")
SCRIPT = os.path.join(DIR, "rescoring_mmgbsa_robusto.py")
PYTHON = "/home/ubuntu/miniforge3/envs/mmgbsa/bin/python"

TIMEOUT = 1800
COLS = ["ligand", "target", "affinity", "mmgbsa_dG", "e_complex", "e_receptor",
        "e_ligand", "error", "papel", "quimia"]

# Un hilo por ligando (OPENMM_CPU_THREADS=1): con varios hilos la suma de
# fuerzas de la plataforma CPU lleva otro orden cada vez, la minimizacion toma
# caminos distintos y el dG se movia hasta 6 kcal/mol entre corridas iguales.
# Con un hilo el calculo es determinista; el paralelismo se hace con workers.
ENTORNO = {**os.environ, "OPENMM_CPU_THREADS": "1", "OMP_NUM_THREADS": "1"}


def log(m):
    linea = "[%s] %s" % (time.strftime("%H:%M:%S"), m)
    print(linea, flush=True)
    with open(LOG, "a") as f:
        f.write(linea + "\n")


def una(t):
    idx, lig = t
    t0 = time.time()
    try:
        p = subprocess.run([PYTHON, script, "--candidates", CAND,
                            "--receptores", RECEP, "--row", str(idx)],
                           capture_output=True, text=True, timeout=TIMEOUT,
                           env=ENTORNO)
        if p.returncode != 0 or not p.stdout.strip():
            return idx, lig, {"error": (p.stderr or "sin salida")[-300:]},
            time.time() - t0
        lineas = [l for l in p.stdout.splitlines() if l.strip().startswith("{")]
        if not lineas:
            return idx, lig, {"error": "sin JSON en la salida"}, time.time() - t0
        return idx, lig, json.loads(lineas[-1]), time.time() - t0
    except subprocess.TimeoutExpired:
        return idx, lig, {"error": "timeout de %d s" % TIMEOUT}, time.time() - t0
    except Exception as e:
        return idx, lig, {"error": str(e)[:300]}, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--script", default=SCRIPT,
                    help="script de rescoring a lanzar por fila (v6 o v7/bcc)")
    ap.add_argument("--out", default=OUT,
                    help="CSV de salida (uno distinto por variante)")
    args = ap.parse_args()
    script, out_csv = args.script, args.out

    cands = list(csv.DictReader(open(CAND)))
    log("candidatos: %d" % len(cands))

    hechos = set()
    if os.path.exists(out_csv):
        for r in csv.DictReader(open(out_csv)):
            if r.get("error") in ("", None) and r.get("mmgbsa_dG"):
                hechos.add(r["ligand"])
    log("ya hechos (sin error): %d" % len(hechos))

    pend = [(i, r) for i, r in enumerate(cands) if r["ligand"] not in hechos]
    log("pendientes: %d | workers: %d | script: %s"
        % (len(pend), args.workers, os.path.basename(script)))
    if not pend:
        log("TODO HECHO")
        return 0

    nuevo = not os.path.exists(out_csv) or os.path.getsize(out_csv) == 0
    outf = open(out_csv, "a", newline="")
    w = csv.writer(outf)
    if nuevo:
        w.writerow(COLS)

    t0 = time.time()
    n_ok, n_err = 0, 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        for idx, lig, res, dt in ex.map(una, [(i, r["ligand"]) for i, r in pend]):
            fila = cands[idx]
            err = res.get("error") or ""
            w.writerow([fila["ligand"], fila["target"], fila["affinity"],
                        res.get("mmgbsa_dG", ""), res.get("e_complex", ""),
                        res.get("e_receptor", ""), res.get("e_ligand", ""),
                        err, fila["papel"], fila["quimia"]])
            outf.flush()
            if err:
                n_err += 1
                log("%-30s ERROR en %.0f s: %s" % (lig, dt, err[:110]))
            else:
                n_ok += 1
                log("%-30s dG=%8.2f  (%.0f s)  [%d ok / %d err]"
                    % (lig, res["mmgbsa_dG"], dt, n_ok, n_err))

    log("")
    log("terminado: %d ok, %d con error, en %.1f min"
        % (n_ok, n_err, (time.time() - t0) / 60.0))
    log("salida: %s" % out_csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
