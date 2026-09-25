#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""acoplar_fondo_rbind.py — acopla el fondo duro (R-BIND 2.0) contra la diana.

Es el mismo acoplamiento que el de los positivos y los señuelos, para que la
comparacion sea legitima: mismo receptor, misma caja, misma exhaustividad, la
misma receta de preparacion (`preparar_ligando`). Si se cambiara cualquiera de
las cuatro, el fondo no seria comparable y la medida no valdria nada.

Los ligandos van a `_validacion_TDP43_rbind/ligands/` con nombre `RB_SM_0001` y
las poses a `_validacion_TDP43_rbind/out/RB_SM_0001_out.pdbqt`, que es
exactamente la convencion que espera `validar_sod1_limpia.py` para leer un fondo.

Es reanudable: lo que ya tiene pose con `REMARK VINA RESULT` no se repite. Se
puede cortar y volver a lanzar sin perder nada.

Uso:
    python acoplar_fondo_rbind.py              # todos
    python acoplar_fondo_rbind.py --limite 3   # solo 3, para medir tiempos

Se puede partir en dos procesos que van uno por cada punta de la lista
(`--reversa --log acoplar_rev.log`), que es lo que se hizo el 23 de septiembre
cuando el acople se freno con los ligandos grandes: la maquina tiene 20 nucleos y
cada vina va con `--cpu 1`, asi que con seis procesos sobraban catorce.
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

# La diana: los mismos cuatro numeros que usa la validacion limpia de TDP-43
# (`validar_diana_limpia.py`, bloque TDP43).
RECEPTOR = os.path.join(BASE, "_tdp43_bolsillo_v2", "4BS2_ph74.pdbqt")
CENTRO = (24.23, 16.89, -15.87)
TAMANO = 26
EXHAUSTIVIDAD = 16

SALIDA = os.path.join(BASE, "_validacion_TDP43_rbind")
LIGDIR = os.path.join(SALIDA, "ligands")
OUTDIR = os.path.join(SALIDA, "out")
LOG = os.path.join(SALIDA, "acoplar.log")
CSV_ENTRADA = os.path.join(BASE, "_rbind", "rbind_fondo_sm.csv")
HILOS = 6

sys.path.insert(0, RAIZ)
from preparar_ligando import escribir   # noqa: E402


def log(m):
    print(m, flush=True)
    os.makedirs(SALIDA, exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")


def nombre_de(rbind_id):
    """'R-BIND (SM) 0001' -> 'RB_SM_0001'."""
    return "RB_" + rbind_id.replace("R-BIND (SM) ", "SM_").replace(" ", "_")


def tiene_pose(ruta):
    if not os.path.exists(ruta):
        return False
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            return True
    return False


def _dock(t):
    vina, rec, lig, out, cx, cy, cz, tam, exh = t
    subprocess.run([vina, "--receptor", rec, "--ligand", lig,
                    "--center_x", str(cx), "--center_y", str(cy),
                    "--center_z", str(cz),
                    "--size_x", str(tam), "--size_y", str(tam), "--size_z", str(tam),
                    "--exhaustiveness", str(exh), "--num_modes", "3",
                    "--out", out, "--cpu", "1"],
                   capture_output=True, timeout=7200)
    return out


def main():
    global LOG, HILOS
    ap = argparse.ArgumentParser()
    ap.add_argument("--limite", type=int, default=0)
    ap.add_argument("--reversa", action="store_true",
                    help="recorrer la lista por el final (segundo proceso)")
    ap.add_argument("--log", default=None, help="fichero de registro alternativo")
    ap.add_argument("--hilos", type=int, default=HILOS)
    args = ap.parse_args()
    if args.log:
        LOG = os.path.join(SALIDA, args.log)
    HILOS = args.hilos

    os.makedirs(LIGDIR, exist_ok=True)
    os.makedirs(OUTDIR, exist_ok=True)
    open(LOG, "w", encoding="utf-8").close()
    log("FONDO DURO R-BIND 2.0 CONTRA TDP-43   %s"
        % time.strftime("%Y-%m-%d %H:%M"))
    log("receptor  %s" % os.path.basename(RECEPTOR))
    log("caja      %s tamano %d exhaustividad %d" % (CENTRO, TAMANO, EXHAUSTIVIDAD))

    filas = [r for r in csv.DictReader(open(CSV_ENTRADA, encoding="utf-8"))
             if r["estado"].startswith("entra")]
    if args.reversa:
        filas = list(reversed(filas))
    if args.limite:
        filas = filas[:args.limite]
    log("ligandos del fondo duro: %d" % len(filas))

    # ------------------------------------------------------------- preparacion
    tareas, fallos = [], []
    t0 = time.time()
    for i, r in enumerate(filas, 1):
        nombre = nombre_de(r["rbind_id"])
        out = os.path.join(OUTDIR, nombre + "_out.pdbqt")
        if tiene_pose(out):
            continue
        p = escribir(nombre, r["smiles"], LIGDIR)
        if p is None:
            fallos.append((nombre, r["nombre"], "no se pudo preparar"))
            continue
        tareas.append((VINA, RECEPTOR, p, out, CENTRO[0], CENTRO[1], CENTRO[2],
                       TAMANO, EXHAUSTIVIDAD))
        if i % 25 == 0:
            log("   preparados %d/%d (%.0f s)" % (i, len(filas), time.time() - t0))
    log("")
    log("a acoplar: %d   (fallos de preparacion: %d)" % (len(tareas), len(fallos)))
    for n, nom, m in fallos:
        log("   %-16s %-30s %s" % (n, nom[:30], m))

    # -------------------------------------------------------------- acoplamiento
    t0 = time.time()
    hechos = 0
    if tareas:
        with ProcessPoolExecutor(max_workers=HILOS) as ex:
            for out in ex.map(_dock, tareas):
                hechos += 1
                if hechos % 10 == 0 or hechos == len(tareas):
                    log("   acoplados %d/%d  (%.1f min, %.0f s por ligando)"
                        % (hechos, len(tareas), (time.time() - t0) / 60.0,
                           (time.time() - t0) / max(1, hechos)))

    # ------------------------------------------------------------------ recuento
    con_pose = 0
    sin_pose = []
    for r in filas:
        nombre = nombre_de(r["rbind_id"])
        if tiene_pose(os.path.join(OUTDIR, nombre + "_out.pdbqt")):
            con_pose += 1
        else:
            sin_pose.append(nombre)
    log("")
    log("RESULTADO: %d/%d con pose" % (con_pose, len(filas)))
    if sin_pose:
        log("   sin pose: %s" % ", ".join(sin_pose[:20]))
    log("poses en %s" % OUTDIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
