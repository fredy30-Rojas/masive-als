#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Acopla el banco CDK2 de DUD-E con Vina-GPU, en la RTX 4080.

ES EL DE TBK1, COPIADO
----------------------
Misma caja, mismos parametros, mismos lotes de 500 con volcado en vivo y el
mismo filtro de tipos de atomo que AutoDock no entiende. No se inventa nada: la
pregunta que responde es si el embudo, tal cual esta, separa activos de senuelos
en un conjunto donde se sabe la respuesta.

LA DIFERENCIA (una): `--decoys N` acota los senuelos con semilla fija, para tener
una lectura en horas en vez de en un dia. Con 474 activos y 4.740 senuelos (1:10)
el margen del AUROC es de unas dos centesimas: suficiente para decidir. El
AUROC no depende de la proporcion (es un rango), asi que la submuestra no
cambia el numero esperado, solo su precision.

Reanudable: salta lo que ya tiene pose, asi que se puede parar y relanzar.

Uso:
    python lanzar_banco_cdk2.py --decoys 4740     # primera lectura (1:10)
    python lanzar_banco_cdk2.py                    # los 27.850 senuelos
    python lanzar_banco_cdk2.py --solo-listar
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import random
import shutil
import subprocess
import time

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(BASE))
GPU = os.path.join(RAIZ, "gpu_dock")
EXE = os.path.join(GPU, "Vina-GPU-2.1-win.exe")
LIGANDS = os.path.join(BASE, "ligands")
SALIDA = os.path.join(BASE, "out")
CSV_SALIDA = os.path.join(BASE, "resultados_cdk2.csv")
CAJA = os.path.join(BASE, "caja_cdk2.json")
LOG = os.path.join(BASE, "banco_cdk2_docking.log")

SEARCH_DEPTH = 20
NUM_MODES = 3
THREAD = 8000
SEMILLA = 20260929

TIPOS_AUTODOCK = {"A", "C", "N", "NA", "OA", "SA", "S", "HD", "H",
                  "F", "Cl", "Br", "I", "P"}


def log(m):
    linea = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M"), m)
    print(linea, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(linea + "\n")


def leer_caja():
    if not os.path.exists(CAJA):
        raise SystemExit("falta caja_cdk2.json: ejecuta antes preparar_receptor_cdk2.py")
    c = json.load(open(CAJA, encoding="utf-8"))
    receptor = os.path.join(BASE, c["receptor_pdbqt"])
    if not os.path.exists(receptor):
        raise SystemExit("falta el receptor %s" % receptor)
    return c, receptor


def escribir_config(receptor, ligdir, outdir, caja, cfg):
    with open(cfg, "w", encoding="utf-8") as f:
        f.write("receptor = %s\n" % receptor.replace("\\", "/"))
        f.write("ligand_directory = %s\n" % ligdir.replace("\\", "/"))
        f.write("output_directory = %s\n" % outdir.replace("\\", "/"))
        f.write("opencl_binary_path = %s\n" % GPU.replace("\\", "/"))
        f.write("center_x = %s\ncenter_y = %s\ncenter_z = %s\n"
                % tuple(caja["centro_caja"]))
        f.write("size_x = %d\nsize_y = %d\nsize_z = %d\n"
                % ((caja["tamano_caja"],) * 3))
        f.write("search_depth = %d\nnum_modes = %d\nthread = %d\n"
                % (SEARCH_DEPTH, NUM_MODES, THREAD))


def tipos_de(ruta):
    out = set()
    with open(ruta, encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                out.add(l[77:79].strip())
    return out


def separables(ruta):
    return tipos_de(ruta) <= TIPOS_AUTODOCK


def volcar(outdir):
    n = 0
    for p in glob.glob(os.path.join(outdir, "*_out.pdbqt")):
        try:
            os.replace(p, os.path.join(SALIDA, os.path.basename(p)))
            n += 1
        except OSError:
            pass
    return n


def afinidad(ruta):
    try:
        with open(ruta, encoding="utf-8", errors="ignore") as f:
            for l in f:
                if "REMARK VINA RESULT" in l:
                    return float(l.split()[3])
    except OSError:
        return None
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--decoys", type=int, default=0,
                    help="maximo de senuelos (0 = todos). Semilla fija.")
    ap.add_argument("--lote", type=int, default=500)
    ap.add_argument("--solo-listar", action="store_true")
    args = ap.parse_args()

    caja, receptor = leer_caja()
    log("receptor: %s" % os.path.basename(receptor))
    log("caja: centro %s, %d A | search_depth %d num_modes %d thread %d"
        % (caja["centro_caja"], caja["tamano_caja"], SEARCH_DEPTH, NUM_MODES,
           THREAD))

    todos = sorted(glob.glob(os.path.join(LIGANDS, "*.pdbqt")))
    if not todos:
        raise SystemExit("no hay ligandos en ligands/ (ejecuta preparar_banco_cdk2.py)")
    activos = [p for p in todos if os.path.basename(p).startswith("ACT_")]
    decoys = [p for p in todos if os.path.basename(p).startswith("DEC_")]
    log("banco: %d activos y %d senuelos" % (len(activos), len(decoys)))

    ligandos, descartados = [], []
    for p in todos:
        (ligandos if separables(p) else descartados).append(p)
    if descartados:
        con_tipo = {}
        for p in descartados:
            for t in tipos_de(p) - TIPOS_AUTODOCK:
                con_tipo[t] = con_tipo.get(t, 0) + 1
        log("descartados %d: tipo de atomo que AutoDock no entiende (%s)"
            % (len(descartados), ", ".join("%s=%d" % kv for kv in sorted(con_tipo.items()))))
        with open(os.path.join(BASE, "descartados_cdk2.txt"), "w",
                  encoding="utf-8") as f:
            for p in descartados:
                f.write("%s\t%s\n" % (os.path.basename(p),
                                      " ".join(sorted(tipos_de(p) - TIPOS_AUTODOCK))))

    if args.decoys and args.decoys < len(decoys):
        rng = random.Random(SEMILLA)
        keep = rng.sample(decoys, args.decoys)
        ligandos = [p for p in ligandos if os.path.basename(p).startswith("ACT_")] + keep
        log("senuelos acotados a %d (semilla %d, proporcion 1:%.1f)"
            % (args.decoys, SEMILLA, args.decoys / max(len(activos), 1)))

    os.makedirs(SALIDA, exist_ok=True)
    ya = set(os.path.basename(p)[:-len("_out.pdbqt")]
             for p in glob.glob(os.path.join(SALIDA, "*_out.pdbqt")))
    pendientes = [p for p in ligandos
                  if os.path.basename(p)[:-len(".pdbqt")] not in ya]
    log("ya con pose: %d | pendientes: %d" % (len(ya), len(pendientes)))
    if args.solo_listar:
        return 0

    t0 = time.time()
    hechos_inicio = len(ya)
    while pendientes:
        lote = pendientes[:args.lote]
        pendientes = pendientes[args.lote:]
        tmp = os.path.join(BASE, "_lote_tmp")
        shutil.rmtree(tmp, ignore_errors=True)
        ligdir = os.path.join(tmp, "ligands")
        outdir = os.path.join(tmp, "out")
        os.makedirs(ligdir)
        os.makedirs(outdir)
        for p in lote:
            shutil.copy2(p, os.path.join(ligdir, os.path.basename(p)))
        cfg = os.path.join(tmp, "config.txt")
        escribir_config(receptor, ligdir, outdir, caja, cfg)
        log("lote de %d ligandos -> Vina-GPU" % len(lote))
        with open(os.path.join(BASE, "motor_stdout.log"), "a",
                  encoding="utf-8", errors="ignore") as f:
            f.write("\n===== lote de %d, %s =====\n" % (len(lote), time.strftime("%Y-%m-%d %H:%M")))
            f.flush()
            proc = subprocess.Popen([EXE, "--config", cfg], cwd=GPU,
                                    stdout=f, stderr=subprocess.STDOUT)
            while proc.poll() is None:
                time.sleep(15)
                volcar(outdir)
            n = volcar(outdir)
        log("   codigo de salida %s | %d poses volcadas" % (proc.returncode, n))
        shutil.rmtree(tmp, ignore_errors=True)
        n = len(glob.glob(os.path.join(SALIDA, "*_out.pdbqt")))
        ritmo = (time.time() - t0) / max(n - hechos_inicio, 1)
        log("   poses guardadas: %d/%d | %.2f s por ligando | faltan ~%.1f h"
            % (n, len(ligandos), ritmo, ritmo * len(pendientes) / 3600.0))
        if not lote:
            break

    filas = []
    for p in sorted(glob.glob(os.path.join(SALIDA, "*_out.pdbqt"))):
        e = afinidad(p)
        if e is not None:
            filas.append({"ligand": os.path.basename(p)[:-len("_out.pdbqt")],
                          "target": "CDK2", "energy": e})
    with open(CSV_SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ligand", "target", "energy"])
        w.writeheader()
        w.writerows(filas)
    log("CSV: %s (%d filas) en %.1f min"
        % (os.path.basename(CSV_SALIDA), len(filas), (time.time() - t0) / 60.0))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
