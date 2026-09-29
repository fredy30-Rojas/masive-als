#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""acoplar_fondo_rbind_sod1.py — el fondo duro (R-BIND 2.0) contra SOD1.

POR QUE HACE FALTA
-----------------
El validador de SOD1 (`validar_sod1_limpia.py`) tiene mounted el segundo fondo, el
"duro", pero DESACTIVADO: `POSES_FONDO2 = None`. Se dejo asi porque no habia ninguna
pose de unidores de ARN de R-BIND 2.0 en la caja de SOD1: las 152 que hay estan
acopladas contra TDP-43, y una pose en otra caja no sirve para nada.

Sin ese bloque, la validacion de SOD1 responde a UNA sola pregunta: si el embudo
reconoce alguna quimica frente a senuelos emparejados. La pregunta que de verdad
importa en una proteina que se une a ARN es otra: si separa quimica ESPECIFICA de
quimica GENERICA de ARN. Un farmaco de TDP-43 se parece mucho a unPegamento de ARN,
y si el embudo no los separa, no esta(selectividad.

LO QUE NO HACE: MEZCLAR
-----------------------
Los dos fondos NO se suman en una sola lista. Cada bloque se puntua con sus propios
positivos y su propio fondo, y asi el puesto de un positivo no depende de cuantos
pegamentos de ARN haya. Esa decision esta escrita en el docstring de
`validar_sod1_limpia.py` y aqui se respeta.

MISMO RECEPTOR, MISMA CAJA, MISMA EXHAUSTIVIDAD
-----------------------------------------------
Es lo que hace que la comparacion sea legitima. Si se cambiara cualquiera de las
cuatro cosas (receptor, caja, exhaustividad, receta de preparacion), el fondo no
seria comparable con los positivos y la medida no valdria. Los cuatro numeros estan
tomados del propio validador de SOD1, no de memoria:

    receptor       gpu_dock/SOD1_limpio.pdbqt
    caja           (46.5, 80.0, 73.3), tamano 22
    exhaustividad  8

Y el motor es `vina.exe` de CPU, no Vina-GPU: el fondo duro de TDP-43 se acoplo
igual, y asi las dos dianas se puntuan con el mismo motor. No toca la tarjeta, que
ahora mismo esta ocupada con el banco de TBK1.

ES REANUDABLE
-------------
Lo que ya tiene pose con `REMARK VINA RESULT` no se repite. Se puede cortar y volver
a lanzar sin perder nada, que con 159 ligandos y 20 nucleos lleva rato.

Uso:
    python acoplar_fondo_rbind_sod1.py
    python acoplar_fondo_rbind_sod1.py --limite 3    # medir tiempos
"""
import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
VINA = os.path.join(RAIZ, "tools", "vina.exe")

# --- los cuatro numeros, copiados del validador de SOD1, no puestos de memoria ---
RECEPTOR = os.path.join(RAIZ, "gpu_dock", "SOD1_limpio.pdbqt")
CENTRO = (46.5, 80.0, 73.3)
TAMANO = 22
EXHAUSTIVIDAD = 8

SALIDA = os.path.join(BASE, "_validacion_SOD1_rbind")
LIGDIR = os.path.join(SALIDA, "ligands")
OUTDIR = os.path.join(SALIDA, "out")
LOG = os.path.join(SALIDA, "acoplar.log")
CSV_ENTRADA = os.path.join(BASE, "_rbind", "rbind_fondo_sm.csv")
HILOS = 6

sys.path.insert(0, RAIZ)
from preparar_ligando import escribir   # noqa: E402


def log(m, destino=None):
    linea = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M"), m)
    print(linea, flush=True)
    with open(destino or LOG, "a", encoding="utf-8") as f:
        f.write(linea + "\n")


def nombre_de(rbind_id):
    """`R-BIND (SM) 0001` -> `RB_SM_0001`, la convencion que espera el validador."""
    d = rbind_id.replace("(SM)", "").strip().replace(" ", "_")
    return d


def tiene_pose(ruta):
    try:
        with open(ruta, encoding="utf-8", errors="ignore") as f:
            return "REMARK VINA RESULT" in f.read()
    except OSError:
        return False


def acoplar(tarea):
    """Un acoplamiento. Las opciones son DOBLES a proposito: este vina.exe
    (v1.2.3) no acepta `-r`, falla al instante con "unrecognised option", y es un
    fallo silencioso si no se mira el codigo de salida: parece un fallo de tiempo
    o de receptor, y no es ninguna de las dos cosas."""
    exe, rec, lig, out, cx, cy, cz, tam, exh = tarea
    cmd = [exe, "--receptor", rec, "--ligand", lig,
           "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", str(tam), "--size_y", str(tam), "--size_z", str(tam),
           "--exhaustiveness", str(exh), "--num_modes", "3",
           "--out", out, "--cpu", "1"]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        ok = tiene_pose(out)
        detalle = "" if ok else (r.stderr or r.stdout or "")[-160:].replace("\n", " ")
        return os.path.basename(lig), r.returncode, ok, detalle
    except Exception as e:  # noqa: BLE001
        return os.path.basename(lig), -1, False, repr(e)[:120]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--limite", type=int, default=0)
    ap.add_argument("--hilos", type=int, default=HILOS)
    args = ap.parse_args()

    if not os.path.exists(RECEPTOR):
        raise SystemExit("no esta el receptor %s" % RECEPTOR)
    if not os.path.exists(VINA):
        raise SystemExit("no esta el vina de CPU %s" % VINA)
    os.makedirs(LIGDIR, exist_ok=True)
    os.makedirs(OUTDIR, exist_ok=True)
    open(LOG, "w", encoding="utf-8").close()

    log("FONDO DURO R-BIND 2.0 CONTRA SOD1   %s"
        % time.strftime("%Y-%m-%d %H:%M"))
    log("receptor  %s" % os.path.basename(RECEPTOR))
    log("caja      %s tamano %d exhaustividad %d (Vina de CPU, %d hilos)"
        % (CENTRO, TAMANO, EXHAUSTIVIDAD, args.hilos))

    filas = [r for r in csv.DictReader(open(CSV_ENTRADA, encoding="utf-8"))
             if r["estado"].startswith("entra")]
    if args.limite:
        filas = filas[:args.limite]
    log("ligandos del fondo duro: %d" % len(filas))

    tareas, fallos = [], []
    t0 = time.time()
    for i, r in enumerate(filas, 1):
        nombre = nombre_de(r["rbind_id"])
        out = os.path.join(OUTDIR, nombre + "_out.pdbqt")
        if tiene_pose(out):
            continue
        try:
            # `escribir` devuelve la ruta del PDBQT, o None si no pudo.
            p = escribir(nombre, r["smiles"], LIGDIR)
        except Exception as e:  # noqa: BLE001
            p = None
            fallos.append((nombre, "no se pudo preparar: %r" % e))
        if not p or not os.path.exists(p):
            if not any(n == nombre for n, _ in fallos):
                fallos.append((nombre, "no se pudo preparar"))
            continue
        tareas.append((VINA, RECEPTOR, p, out, CENTRO[0], CENTRO[1], CENTRO[2],
                       TAMANO, EXHAUSTIVIDAD))
        if i % 25 == 0:
            log("   preparados %d/%d (%.0f s)" % (i, len(filas), time.time() - t0))
    log("")
    log("a acoplar: %d   (fallos de preparacion: %d)" % (len(tareas), len(fallos)))
    for n, m in fallos:
        log("   %-16s %s" % (n, m))

    t0 = time.time()
    n_ok = 0
    with ProcessPoolExecutor(max_workers=args.hilos) as ex:
        for k, (nombre, code, ok, err) in enumerate(ex.map(acoplar, tareas), 1):
            if ok:
                n_ok += 1
            else:
                log("   FALLO %-16s codigo %s %s" % (nombre, code, err), LOG)
            if k % 20 == 0 or k == len(tareas):
                log("   %d/%d acoplados, %d con pose (%.0f s)"
                    % (k, len(tareas), n_ok, time.time() - t0))

    total = len([f for f in os.listdir(OUTDIR) if f.endswith("_out.pdbqt")])
    log("")
    log("poses en %s: %d" % (OUTDIR, total))
    log("para activarlo en el validador: POSES_FONDO2 y LIGS_FONDO2 a esta carpeta")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
