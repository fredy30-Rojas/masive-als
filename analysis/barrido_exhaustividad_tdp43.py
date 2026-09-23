#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""barrido_exhaustividad_tdp43.py — el MISMO conjunto de TDP-43 acoplado con Vina
de CPU a varias exhaustividades, para separar el ruido de la perilla del ruido
del motor.

POR QUE
-------
El control con Vina-GPU (INFORME_CONTROL_VINAGPU_TDP43_2026-09-23.md) enseño que
el veredicto de quimiotipos cambia al cambiar de motor: con los señuelos
emparejados, la CPU (exhaustividad 16) dice PASA con 2 de 5 y la GPU
(search_depth 20) dice NO PASA con 1 de 5. Pero ese cambio puede venir de dos
sitios distintos y el control no los separa:

  a) del MOTOR (Vina 1.2.3 y Vina-GPU 2.1 son implementaciones distintas), o
  b) del ESFUERZO DE BUSQUEDA (`exhaustiveness` de Vina no es `search_depth` de
     Vina-GPU: son perillas distintas y no hay equivalencia exacta).

Este barrido mueve SOLO el esfuerzo, dentro del mismo motor: la misma CPU, el
mismo Vina, el mismo receptor, la misma caja, los mismos ficheros de ligando.
Si el veredicto tambien se mueve aqui, entonces lo que no aguanta es el criterio
(que decide cuantos quimiotipos cruzan un corte), no el motor.

QUE SE ACOPLA
-------------
Los 289 ligandos del control: los 7 positivos de union medida, los 122 señuelos
emparejados y el fondo duro de R-BIND 2.0 (152). No se re-prepara nada: se copian
los `.pdbqt` ya preparados con la receta canonica.

El punto medio ya existe y NO se repite: la validacion limpia a exhaustividad 16
(`_validacion_TDP43` + `_validacion_TDP43_rbind`).

Uso:
    python barrido_exhaustividad_tdp43.py --todo
    python barrido_exhaustividad_tdp43.py --acoplar 8      # solo una
    python barrido_exhaustividad_tdp43.py --puntuar 8
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
VINA = os.path.join(RAIZ, "tools", "vina.exe")

RECEPTOR = os.path.join(BASE, "_tdp43_bolsillo_v2", "4BS2_ph74.pdbqt")
CENTRO = (24.23, 16.89, -15.87)
TAMANO = 26
NUM_MODES = 3
HILOS = 14            # la maquina tiene 20 nucleos; se dejan seis libres
EXHAUSTIVIDADES = [8, 32]

SALIDA = os.path.join(BASE, "_barrido_TDP43")
LOG = os.path.join(SALIDA, "barrido.log")
PREFIJO_DURO = "RB_SM_"

ORIGENES = [os.path.join(BASE, "_validacion_TDP43", "ligands"),
            os.path.join(BASE, "_validacion_TDP43_rbind", "ligands"),
            os.path.join(BASE, "validar_tdp43_limpia", "ligands")]


def log(m):
    print(m, flush=True)
    os.makedirs(SALIDA, exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write("[%s] %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), m))


def dirs(exh):
    d = os.path.join(SALIDA, "exh%d" % exh)
    return d, os.path.join(d, "ligands"), os.path.join(d, "out"), os.path.join(d, "out_duro")


def _dock(t):
    vina, rec, lig, out, cx, cy, cz, tam, exh = t
    subprocess.run([vina, "--receptor", rec, "--ligand", lig,
                    "--center_x", str(cx), "--center_y", str(cy),
                    "--center_z", str(cz),
                    "--size_x", str(tam), "--size_y", str(tam), "--size_z", str(tam),
                    "--exhaustiveness", str(exh), "--num_modes", str(NUM_MODES),
                    "--out", out, "--cpu", "1"],
                   capture_output=True, timeout=7200)
    return out


def tiene_resultado(ruta):
    if not ruta or not os.path.exists(ruta):
        return False
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            return True
    return False


def montar(exh):
    """Copia los ligandos ya preparados a la carpeta de esta exhaustividad."""
    _, ligs, _, _ = dirs(exh)
    os.makedirs(ligs, exist_ok=True)
    for origen in ORIGENES:
        if not os.path.isdir(origen):
            continue
        for p in sorted(glob.glob(os.path.join(origen, "*.pdbqt"))):
            destino = os.path.join(ligs, os.path.basename(p))
            if not os.path.exists(destino):
                shutil.copy2(p, destino)
    return len(glob.glob(os.path.join(ligs, "*.pdbqt")))


def acoplar(exh, limite=0):
    d, ligs, out, out_duro = dirs(exh)
    os.makedirs(out, exist_ok=True)
    os.makedirs(out_duro, exist_ok=True)
    n = montar(exh)
    log("exhaustividad %d: %d ligandos en el montaje" % (exh, n))

    tareas = []
    for p in sorted(glob.glob(os.path.join(ligs, "*.pdbqt"))):
        nombre = os.path.basename(p)[:-len(".pdbqt")]
        destino = os.path.join(out_duro if nombre.startswith(PREFIJO_DURO) else out,
                               nombre + "_out.pdbqt")
        if tiene_resultado(destino):
            continue
        tareas.append((VINA, RECEPTOR, p, destino, CENTRO[0], CENTRO[1], CENTRO[2],
                       TAMANO, exh))
    ya = n - len(tareas)
    if limite:
        tareas = tareas[:limite]
    log("exhaustividad %d: %d ya acoplados, %d por acoplar (hilos %d, "
        "num_modes %d, caja %d A)"
        % (exh, ya, len(tareas), HILOS, NUM_MODES, TAMANO))
    if not tareas:
        return 0

    t0 = time.time()
    hechos = 0
    with ProcessPoolExecutor(max_workers=HILOS) as ex:
        for _ in ex.map(_dock, tareas):
            hechos += 1
            if hechos % 10 == 0 or hechos == len(tareas):
                log("   exhaustividad %d: acoplados %d/%d (%.1f min, %.1f min por "
                    "ligando)" % (exh, hechos, len(tareas),
                                  (time.time() - t0) / 60.0,
                                  (time.time() - t0) / 60.0 / hechos))
    faltan = [t[3] for t in tareas if not tiene_resultado(t[3])]
    log("exhaustividad %d: terminado en %.1f min; sin pose: %d"
        % (exh, (time.time() - t0) / 60.0, len(faltan)))
    for f in faltan[:10]:
        log("   sin pose: %s" % os.path.basename(f))
    return 0


def puntuar(exh):
    d, ligs, out, out_duro = dirs(exh)
    log("exhaustividad %d: puntuando (%d poses + %d del fondo duro)"
        % (exh, len(glob.glob(os.path.join(out, "*_out.pdbqt"))),
           len(glob.glob(os.path.join(out_duro, "*_out.pdbqt")))))
    sys.path.insert(0, BASE)
    import validar_sod1_limpia as V

    V.RECEPTOR = RECEPTOR
    V.CENTRO = CENTRO
    V.TAMANO = TAMANO
    V.EXHAUSTIVIDAD = exh
    V.POSES_FONDO = out
    V.LIGS_FONDO = ligs
    V.SALIDA = os.path.join(SALIDA, "validar_exh%d" % exh)
    V.DIANA = "TDP43"
    V.NOTA_FONDO = ("122 señuelos emparejados en propiedades, acoplados con Vina "
                    "CPU a exhaustividad %d" % exh)
    V.ETIQUETA_FONDO = "señuelos emparejados (exh %d)" % exh
    V.POSES_FONDO2 = out_duro
    V.LIGS_FONDO2 = ligs
    V.NOTA_FONDO2 = ("unidores de ARN de R-BIND 2.0, acoplados con Vina CPU a "
                     "exhaustividad %d en la misma caja" % exh)
    V.ETIQUETA_FONDO2 = "R-BIND 2.0 (exh %d)" % exh
    V.log("BARRIDO DE EXHAUSTIVIDAD: %d   %s" % (exh, time.strftime("%Y-%m-%d %H:%M")))
    return V.main()


def resumen():
    """Deja en un fichero los tres bloques de cada exhaustividad, para leerlos de un golpe."""
    filas = []
    for exh in sorted(EXHAUSTIVIDADES):
        p = os.path.join(SALIDA, "validar_exh%d_resumen.csv" % exh)
        if not os.path.exists(p):
            continue
        import csv as _csv
        r = list(_csv.DictReader(open(p, encoding="utf-8")))[0]
        filas.append({
            "exhaustividad": exh,
            "blando_auc": r.get("auc_crudo"), "blando_atomo": r.get("auc_atomo"),
            "blando_residual": r.get("auc_residual"),
            "blando_quimias": r.get("quimias_pareadas_por_tamano"),
            "blando_veredicto": r.get("veredicto"),
            "duro_auc": r.get("duro_auc_crudo"),
            "duro_residual": r.get("duro_auc_residual"),
            "duro_quimias": r.get("duro_quimias_pareadas_por_tamano"),
            "duro_ganan": r.get("duro_quimias_que_ganan"),
            "duro_veredicto": r.get("duro_veredicto"),
            "pendiente": r.get("pendiente_tamano_kcal_por_atomo"),
        })
    if not filas:
        return
    sal = os.path.join(SALIDA, "barrido_resumen.csv")
    with open(sal, "w", newline="", encoding="utf-8") as f:
        w = _csv.DictWriter(f, fieldnames=list(filas[0]))
        w.writeheader()
        w.writerows(filas)
    log("resumen del barrido en %s" % sal)
    for r in filas:
        log("   exh %-3s blando AUC %s quimias %s %s | duro AUC %s quimias %s %s"
            % (r["exhaustividad"], r["blando_auc"], r["blando_quimias"],
               r["blando_veredicto"], r["duro_auc"], r["duro_quimias"],
               r["duro_veredicto"]))


def main():
    global HILOS
    ap = argparse.ArgumentParser()
    ap.add_argument("--acoplar", type=int, default=None, metavar="EXH")
    ap.add_argument("--puntuar", type=int, default=None, metavar="EXH")
    ap.add_argument("--todo", action="store_true")
    ap.add_argument("--limite", type=int, default=0)
    ap.add_argument("--hilos", type=int, default=HILOS)
    args = ap.parse_args()
    HILOS = args.hilos
    if not any([args.acoplar, args.puntuar, args.todo]):
        ap.error("hay que decir que hacer: --acoplar N, --puntuar N o --todo")

    os.makedirs(SALIDA, exist_ok=True)
    log("")
    log("BARRIDO DE EXHAUSTIVIDAD EN CPU — TDP-43, el mismo conjunto  %s"
        % time.strftime("%Y-%m-%d %H:%M"))
    for exh in EXHAUSTIVIDADES:
        if args.todo or args.acoplar == exh:
            if acoplar(exh, args.limite):
                return 1
        if args.todo or args.puntuar == exh:
            puntuar(exh)
    if args.todo:
        resumen()
    return 0


if __name__ == "__main__":
    sys.exit(main())
