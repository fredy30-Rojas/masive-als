#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Acopla el banco de TBK1 con Vina-GPU, con el MISMO protocolo del proyecto.

PROTOCOLO: COPIADO, NO INVENTADO
--------------------------------
El unico protocolo de GPU que la regla de decision del proyecto valido es el de
`analysis/control_gpu_tdp43.py`: **search_depth 20, num_modes 3, thread 8000**, y
una caja de 24 A. Aqui se usan esos mismos valores, porque este banco existe
justamente para medir ese protocolo en una diana donde se sabe la respuesta.

  OJO: `search_depth` de Vina-GPU NO es `exhaustiveness` de Vina. Son perillas
  distintas y no hay equivalencia exacta. Por eso este banco mide el ORDEN y los
  veredictos, no valores absolutos de afinidad.

RECEPTOR Y CAJA
---------------
Los dos salen de `caja.json`, que lo escribe `construir_receptor_tbk1.py` con la
caja fijada por el ligando cocristalizado BX-795 de la 4EUU. No hay ninguna caja
que elegir aqui.

Uso:
    python lanzar_banco_tbk1.py --solo-listar     # comprueba y no toca la GPU
    python lanzar_banco_tbk1.py --probar 300      # prueba corta y mide el ritmo
    python lanzar_banco_tbk1.py                   # el banco entero
    python lanzar_banco_tbk1.py --lote 500        # por lotes de 500 (por defecto)

Reanudable: se salta los ligandos cuya pose ya este en `out/`, y las poses se
vuelcan en vivo cada 15 s, de modo que una muerte no cuesta mas que el ligando
que estuviera en curso.

Descarta los ligandos con boro, silicio o sodio: AutoDock 4 no tiene esos tipos
y los perdia en silencio (142 de 35.349 el 28 sep). Los nombres quedan en
`descartados_tbk1.txt` para poder contarlos despues.

Salida: `out/` con las poses y `resultados_tbk1.csv` con ligand,target,energy.
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import shutil
import subprocess
import time

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(BASE))
GPU = os.path.join(RAIZ, "gpu_dock")
EXE = os.path.join(GPU, "Vina-GPU-2.1-win.exe")
LIGANDS = os.path.join(BASE, "ligands")
SALIDA = os.path.join(BASE, "out")
CSV_SALIDA = os.path.join(BASE, "resultados_tbk1.csv")
CAJA = os.path.join(BASE, "caja.json")

SEARCH_DEPTH = 20
NUM_MODES = 3
THREAD = 8000

# Tipos de atomo que AutoDock 4 entiende. Los que quedan fuera producen
# "ATOM syntax incorrect" y el ligando se pierde en silencio dentro del lote.
# El 28 sep 2026 se contaron 142 ligandos del banco con boro, silicio o sodio.
TIPOS_AUTODOCK = {"A", "C", "N", "NA", "OA", "SA", "S", "HD", "H",
                  "F", "Cl", "Br", "I", "P"}


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def leer_caja():
    if not os.path.exists(CAJA):
        raise SystemExit("falta caja.json: ejecuta antes construir_receptor_tbk1.py")
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
                % (caja["tamano_caja"], caja["tamano_caja"], caja["tamano_caja"]))
        f.write("search_depth = %d\nnum_modes = %d\nthread = %d\n"
                % (SEARCH_DEPTH, NUM_MODES, THREAD))


def tipos_de(ruta):
    """Tipos de atomo de un pdbqt, leidos de las columnas 78-79 del formato."""
    out = set()
    with open(ruta, encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("ATOM") or l.startswith("HETATM"):
                out.add(l[77:79].strip())
    return out


def separables(ruta):
    """False si el ligando lleva un tipo de atomo que AutoDock no reconoce."""
    return tipos_de(ruta) <= TIPOS_AUTODOCK


def volcar(outdir):
    """Mueve las poses ya escritas por Vina a la carpeta buena.

    Se llama MIENTRAS corre el motor, no al final del lote. El 28 sep el lote de
    5000 se cayo a las tres horas y se perdio entero, porque las poses solo se
    copiaban cuando el lote terminaba. Con lotes de 500 y volcado cada 15 s, una
    muerte se lleva como mucho el ligando que este en curso.
    """
    n = 0
    for p in glob.glob(os.path.join(outdir, "*_out.pdbqt")):
        try:
            os.replace(p, os.path.join(SALIDA, os.path.basename(p)))
            n += 1
        except OSError:
            pass
    return n


def afinidad(ruta):
    """Afinidad de la mejor pose de un *_out.pdbqt."""
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
    ap.add_argument("--solo-listar", action="store_true")
    ap.add_argument("--probar", type=int, default=0)
    ap.add_argument("--lote", type=int, default=500,
                    help="ligandos por lote. Pequeno a proposito: 500 son unos "
                         "22 min, y un lote de 5000 que se cae se lleva 3 h.")
    args = ap.parse_args()

    caja, receptor = leer_caja()
    log("receptor: %s (%s)" % (os.path.relpath(receptor, RAIZ),
                               "+".join(caja["cadenas_receptor"])))
    log("caja: centro %s, %d A | search_depth %d num_modes %d thread %d"
        % (caja["centro_caja"], caja["tamano_caja"], SEARCH_DEPTH, NUM_MODES, THREAD))

    todos = sorted(glob.glob(os.path.join(LIGANDS, "*.pdbqt")))
    if not todos:
        raise SystemExit("no hay ligandos en %s (ejecuta preparar_banco_tbk1.py)"
                         % os.path.relpath(LIGANDS, RAIZ))

    ligandos, descartados = [], []
    for p in todos:
        (ligandos if separables(p) else descartados).append(p)
    log("ligandos en el banco: %d" % len(todos))
    if descartados:
        con_tipo = {}
        for p in descartados:
            for t in tipos_de(p) - TIPOS_AUTODOCK:
                con_tipo[t] = con_tipo.get(t, 0) + 1
        log("descartados %d: tipo de atomo que AutoDock no entiende (%s)"
            % (len(descartados),
               ", ".join("%s=%d" % kv for kv in sorted(con_tipo.items()))))
        with open(os.path.join(BASE, "descartados_tbk1.txt"), "w",
                  encoding="utf-8") as f:
            for p in descartados:
                f.write("%s\t%s\n" % (os.path.basename(p),
                                      " ".join(sorted(tipos_de(p) - TIPOS_AUTODOCK))))

    os.makedirs(SALIDA, exist_ok=True)
    ya = set(os.path.basename(p)[:-len("_out.pdbqt")]
             for p in glob.glob(os.path.join(SALIDA, "*_out.pdbqt")))
    pendientes = [p for p in ligandos
                  if os.path.basename(p)[:-len(".pdbqt")] not in ya]
    log("ya con pose: %d | pendientes: %d" % (len(ya), len(pendientes)))
    if args.solo_listar:
        log("(--solo-listar: no se lanza nada)")
        return 0
    if args.probar:
        pendientes = pendientes[:args.probar]

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
            f.write("\n===== lote de %d, %s =====\n"
                    % (len(lote), time.strftime("%Y-%m-%d %H:%M")))
            f.flush()
            proc = subprocess.Popen([EXE, "--config", cfg], cwd=GPU,
                                    stdout=f, stderr=subprocess.STDOUT)
            # Volcado en vivo: lo que Vina ya escribio pasa a out/ sin esperar
            # al final del lote. Cada 15 s, para no castigar el disco.
            while proc.poll() is None:
                time.sleep(15)
                n = volcar(outdir)
                if n:
                    log("   %d poses volcadas en vivo" % n)
            n = volcar(outdir)
        log("   codigo de salida %s | %d poses volcadas" % (proc.returncode, n))
        shutil.rmtree(tmp, ignore_errors=True)
        n = len(glob.glob(os.path.join(SALIDA, "*_out.pdbqt")))
        ritmo = (time.time() - t0) / max(n - hechos_inicio, 1)
        log("   poses guardadas: %d/%d | %.2f s por ligando | faltan ~%.1f h"
            % (n, len(ligandos), ritmo, ritmo * len(pendientes) / 3600.0))
        if not lote:
            break

    # --- CSV con el mismo formato que resultados_vinagpu_total.csv ---
    filas = []
    for p in sorted(glob.glob(os.path.join(SALIDA, "*_out.pdbqt"))):
        e = afinidad(p)
        if e is not None:
            filas.append({"ligand": os.path.basename(p)[:-len("_out.pdbqt")],
                          "target": "TBK1", "energy": e})
    with open(CSV_SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ligand", "target", "energy"])
        w.writeheader()
        w.writerows(filas)
    log("CSV: %s (%d filas) en %.1f min"
        % (os.path.relpath(CSV_SALIDA, RAIZ), len(filas), (time.time() - t0) / 60.0))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
