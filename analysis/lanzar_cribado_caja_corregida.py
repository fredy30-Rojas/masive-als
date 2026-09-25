#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""lanzar_cribado_caja_corregida.py — cribar la libreria en la caja que SI esta validada.

POR QUE EXISTE
--------------
El cribado grande que hay en el disco (tandas z001-z011, `resultados_vinagpu_total.csv`,
597.917 poses) se hizo en las cajas VIEJAS: TDP-43 con `gpu_dock/TDP43.pdbqt` en
(28.3, 43.7, 52.5) y SOD1 con `gpu_dock/SOD1.pdbqt` en (27.9, 111.8, 64.4), las dos de
25 A. Esas cajas se midieron despues con el embudo y **no pasan**: la de TDP-43 da un
AUC de 0,465, que es el azar (`REVISION_VALIDACION_BOLSILLOS_2026-08-20.md`), y por eso
el 20 de agosto se corrigio la caja.

La caja CORREGIDA es la que el 24 de septiembre paso la regla de decision: TDP-43 con
`analysis/_tdp43_bolsillo_v2/4BS2_ph74.pdbqt`, centro (24.23, 16.89, -15.87), 26 A, AUC
0,731-0,738 y PASA en las cuatro corridas; SOD1 con `gpu_dock/SOD1_limpio.pdbqt`, centro
(46.5, 80.0, 73.3), 22 A. La libreria, en cambio, **nunca se ha cribado en ellas**: solo
se re-acoplaron los de arriba del ranking viejo (162 en `_redock_tdp43_masivo.py`), y con
la caja corregida puntuan mejor de -0,4 a -1,3 kcal/mol, o sea que el ranking viejo no
ordenaba lo mismo.

Este script es el que cierra ese hueco: la misma libreria preparada, el mismo motor y el
MISMO protocolo que el control de GPU que la regla valido.

PROTOCOLO: NO INVENTADO, COPIADO DEL CONTROL
--------------------------------------------
`analysis/control_gpu_tdp43.py`, que es la corrida de GPU que la regla del 24 uso como
cuarta corrida (AUC 0,733): search_depth 20, num_modes 3, thread 8000, y la caja
corregida. Los ligandos NO se re-preparan: se usan los `.pdbqt` ya preparados de
`gpu_dock/libreria_ligands/`, que es lo que hizo el control.

Uso:
    python lanzar_cribado_caja_corregida.py --solo-listar      # no toca la GPU
    python lanzar_cribado_caja_corregida.py --probar 500        # mide el ritmo real
    python lanzar_cribado_caja_corregida.py                     # la libreria entera
    python lanzar_cribado_caja_corregida.py --diana TDP43       # solo una diana
Resumible: salta los ligandos cuya pose ya este en la carpeta de salida.
Salida: gpu_dock/cribado_caja_corregida/<DIANA>/  y un CSV con el mismo formato que
`resultados_vinagpu_total.csv` (ligand,target,energy) para que lo lea el mismo codigo.
"""
import argparse
import csv
import glob
import os
import subprocess
import time

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
GPU = os.path.join(RAIZ, "gpu_dock")
EXE = os.path.join(GPU, "Vina-GPU-2.1-win.exe")
LIGDIR = os.path.join(GPU, "libreria_ligands")
SALIDA = os.path.join(GPU, "cribado_caja_corregida")

# Las DOS cajas corregidas, con el receptor y el tamano de la validacion que paso la
# regla. El orden importa: TDP-43 es la que pasa el bloque duro (PASA) y la que decide;
# SOD1 esta en SIN EVIDENCIA.
DIANAS = {
    "TDP43": {"receptor": os.path.join(BASE, "_tdp43_bolsillo_v2", "4BS2_ph74.pdbqt"),
              "centro": (24.23, 16.89, -15.87), "tamano": 26},
    "SOD1": {"receptor": os.path.join(GPU, "SOD1_limpio.pdbqt"),
             "centro": (46.5, 80.0, 73.3), "tamano": 22},
}

# Iguales que en control_gpu_tdp43.py: es el unico protocolo de GPU que la regla valido.
SEARCH_DEPTH = 20
NUM_MODES = 3
THREAD = 8000


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def posicion_mejor(ruta):
    """Afinidad de la primera pose (la mejor) de un *_out.pdbqt, o None."""
    try:
        with open(ruta, encoding="utf-8", errors="ignore") as f:
            for linea in f:
                if "REMARK VINA RESULT" in linea:
                    return float(linea.split()[3])
    except OSError:
        return None
    return None


def escribir_config(diana, cfg, ligdir, outdir):
    c = DIANAS[diana]
    with open(cfg, "w", encoding="utf-8") as f:
        f.write("receptor = %s\nligand_directory = %s\noutput_directory = %s\n"
                % (c["receptor"], ligdir, outdir))
        f.write("opencl_binary_path = %s\n" % GPU)
        f.write("center_x = %s\ncenter_y = %s\ncenter_z = %s\n"
                % (c["centro"][0], c["centro"][1], c["centro"][2]))
        f.write("size_x = %d\nsize_y = %d\nsize_z = %d\n"
                % (c["tamano"], c["tamano"], c["tamano"]))
        f.write("search_depth = %d\nnum_modes = %d\nthread = %d\n"
                % (SEARCH_DEPTH, NUM_MODES, THREAD))


def lanzar(diana, ligdir, outdir, cfg):
    subprocess.run([EXE, "--config", cfg], cwd=GPU, capture_output=True, timeout=None)
    hechos = glob.glob(os.path.join(outdir, "*_out.pdbqt"))
    return len(hechos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--diana", choices=sorted(DIANAS), default=None)
    ap.add_argument("--solo-listar", action="store_true",
                    help="comprueba receptor, caja y ligandos, escribe el config y para")
    ap.add_argument("--probar", type=int, default=0,
                    help="acopla solo los primeros N ligandos, para medir el ritmo")
    args = ap.parse_args()

    objetivos = [args.diana] if args.diana else ["TDP43", "SOD1"]
    ligandos = sorted(glob.glob(os.path.join(LIGDIR, "*.pdbqt")))
    log("libreria preparada: %d ligandos en %s" % (len(ligandos), os.path.relpath(LIGDIR, RAIZ)))
    if not ligandos:
        raise SystemExit("no hay ligandos preparados en %s" % LIGDIR)

    for diana in objetivos:
        c = DIANAS[diana]
        if not os.path.exists(c["receptor"]):
            raise SystemExit("falta el receptor de %s: %s" % (diana, c["receptor"]))
        outdir = os.path.join(SALIDA, diana)
        cfg = os.path.join(SALIDA, "config_%s.txt" % diana)
        os.makedirs(outdir, exist_ok=True)
        escribir_config(diana, cfg, LIGDIR, outdir)

        ya = set(os.path.basename(p)[:-len("_out.pdbqt")]
                 for p in glob.glob(os.path.join(outdir, "*_out.pdbqt")))
        pendientes = [p for p in ligandos
                      if os.path.basename(p)[:-len(".pdbqt")] not in ya]
        log("%s: receptor %s" % (diana, os.path.relpath(c["receptor"], RAIZ)))
        log("   caja centro %s tamano %d A | search_depth %d num_modes %d thread %d"
            % (c["centro"], c["tamano"], SEARCH_DEPTH, NUM_MODES, THREAD))
        log("   ya hechos: %d | pendientes: %d" % (len(ya), len(pendientes)))
        log("   config: %s" % os.path.relpath(cfg, RAIZ))
        if args.solo_listar:
            log("   (--solo-listar: no se lanza nada; el comando seria"
                " %s --config %s)" % (os.path.relpath(EXE, RAIZ), os.path.relpath(cfg, RAIZ)))
            continue

        if args.probar:
            # Prueba con copias, para no ensuciar la salida buena con un lote corto.
            import shutil
            probdir = os.path.join(SALIDA, "_prueba_%s" % diana, "ligands")
            os.makedirs(probdir, exist_ok=True)
            salida_prob = os.path.join(SALIDA, "_prueba_%s" % diana, "out")
            os.makedirs(salida_prob, exist_ok=True)
            for p in pendientes[:args.probar]:
                shutil.copy2(p, os.path.join(probdir, os.path.basename(p)))
            cfg_prob = os.path.join(SALIDA, "_prueba_%s" % diana, "config.txt")
            escribir_config(diana, cfg_prob, probdir, salida_prob)
            t0 = time.time()
            n = lanzar(diana, probdir, salida_prob, cfg_prob)
            dt = time.time() - t0
            log("   PRUEBA: %d ligandos en %.0f s -> %.2f s por ligando; la libreria"
                " entera serian %.1f h por diana"
                % (n, dt, dt / max(n, 1), dt / max(n, 1) * len(ligandos) / 3600.0))
            continue

        t0 = time.time()
        n = lanzar(diana, LIGDIR, outdir, cfg)
        log("   %s: %d poses nuevas en %.0f s" % (diana, n, time.time() - t0))

        # El CSV, con el mismo formato que resultados_vinagpu_total.csv.
        filas = []
        for p in sorted(glob.glob(os.path.join(outdir, "*_out.pdbqt"))):
            e = posicion_mejor(p)
            if e is not None:
                filas.append({"ligand": os.path.basename(p)[:-len("_out.pdbqt")],
                              "target": diana, "energy": e})
        csv_out = os.path.join(SALIDA, "resultados_%s.csv" % diana)
        with open(csv_out, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["ligand", "target", "energy"])
            w.writeheader()
            w.writerows(filas)
        log("   guardado: %s (%d filas)" % (os.path.relpath(csv_out, RAIZ), len(filas)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
