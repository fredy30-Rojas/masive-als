#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""reparar_sod1_v3_glue.py — aparta los pseudo-atomos glue del conjunto SOD1 v3.

POR QUE
-------
Meeko con los ajustes por defecto no cierra los anillos de 7 eslabones en
adelante: los abre y pega los extremos con pseudo-atomos de pegamento (CG0/G0),
duplicando dos atomos por anillo. Vina ve el anillo abierto, su energia interna
del ligando libre (REMARK UNBOUND) sale a +8 kcal/mol en vez de ~0, y como el
score reportado es inter + intra - unbound, el compuesto recibe un regalo de
3-6 kcal/mol que no existe.

El conjunto SOD1 v3 tiene 11 afectados, todos del FONDO (5 del fondo duro DECH_
y 6 del emparejado DECM_), y varios estan en la cabeza del ranking. Como el AUC
de esta validacion (0.671 crudo, 0.486 sin tamaño) se calcula contra esos fondos,
el glue contribuye a hundirlo.

QUE HACE
--------
1. Detecta los ligandos con tipos CG0/G0 en `validacion_SOD1_v3/ligands`.
2. Aparta la pose vieja en `out/_antes_glue/` y el ligando viejo en
   `ligands/_roto_glue/` (no borra nada).
3. Vuelve a preparar el ligando con `rigid_macrocycles=True` y comprueba que el
   numero de atomos cuadra con su SMILES.
4. Re-acopla SOLO esos ligandos con el MISMO receptor, caja y exhaustividad que
   uso `validar_sod1_v3.py` (gpu_dock/SOD1.pdbqt, centro 46.5/80.0/73.3,
   tamaño 22, exhaustividad 8), para que el antes y el despues sean comparables.
5. Guarda la copia del analisis previo y avisa de que hay que re-lanzar
   `python validar_sod1_v3.py --analizar` para recalcular las metricas.

Uso:
  python reparar_sod1_v3_glue.py
"""
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor

BASE = os.path.dirname(os.path.abspath(__file__))
GPU = os.path.abspath(os.path.join(BASE, "..", "gpu_dock"))
sys.path.insert(0, BASE)

from reparar_ligandos_glue import contar_fichero, reparar          # noqa: E402
from reparar_validaciones_glue import (_dock, leer_afinidad,       # noqa: E402
                                       leer_unbound, log, VINA_CPU)

TRABAJO = os.path.join(BASE, "validacion_SOD1_v3")
LIGDIR = os.path.join(TRABAJO, "ligands")
OUTDIR = os.path.join(TRABAJO, "out")

# Los mismos parametros que validar_sod1_v3.py (lineas 72-75).
RECEPTOR = os.path.join(GPU, "SOD1.pdbqt")
CENTRO = (46.5, 80.0, 73.3)
TAMANO = 22
EXHAUSTIVIDAD = 8
HILOS = 6


def main():
    log("")
    log("=" * 78)
    log("=== SOD1 v3: reparacion de pseudo-atomos glue")
    log("    ligando   %s" % LIGDIR)
    log("    receptor  %s (el original, para que el antes y el despues sean"
        " comparables)" % os.path.basename(RECEPTOR))
    log("    caja      %s tamano %d exhaustividad %d"
        % (CENTRO, TAMANO, EXHAUSTIVIDAD))

    rotos = []
    for p in sorted(os.listdir(LIGDIR)):
        if not p.endswith(".pdbqt"):
            continue
        n, glue = contar_fichero(os.path.join(LIGDIR, p))
        if glue:
            rotos.append((os.path.splitext(p)[0], n, glue))
    if not rotos:
        log("    ningun ligando con pseudo-atomos. Nada que hacer.")
        return 0
    log("    ligandos con pseudo-atomos: %d" % len(rotos))

    antes = os.path.join(OUTDIR, "_antes_glue")
    rotos_dir = os.path.join(LIGDIR, "_roto_glue")
    os.makedirs(antes, exist_ok=True)
    os.makedirs(rotos_dir, exist_ok=True)

    # copia del analisis previo, para poder comparar
    analisis = os.path.join(TRABAJO, "analisis_sod1_v3.csv")
    copia = os.path.join(TRABAJO, "analisis_sod1_v3_antes_glue.csv")
    if os.path.exists(analisis) and not os.path.exists(copia):
        shutil.copy(analisis, copia)
        log("    analisis previo guardado en %s" % os.path.basename(copia))

    tareas, cambios = [], []
    for fichero, n, glue in rotos:
        p = os.path.join(LIGDIR, fichero + ".pdbqt")
        pose = os.path.join(OUTDIR, fichero + "_out.pdbqt")
        unbound_antes = leer_unbound(pose) if os.path.exists(pose) else None
        if os.path.exists(pose):
            shutil.copy(pose, os.path.join(antes, os.path.basename(pose)))
        txt, n_ok, motivo = reparar(p)
        if txt is None:
            log("    %-24s NO reparado: %s" % (fichero, motivo))
            continue
        shutil.copy(p, os.path.join(rotos_dir, fichero + ".pdbqt"))
        with open(p, "w", encoding="utf-8") as fh:
            fh.write(txt)
        if os.path.exists(pose):
            os.remove(pose)
        tareas.append((VINA_CPU, RECEPTOR, p, pose,
                       CENTRO[0], CENTRO[1], CENTRO[2], TAMANO, EXHAUSTIVIDAD))
        cambios.append((fichero, n, n_ok, glue, unbound_antes))

    log("    re-acoplando %d ligandos (Vina CPU, exhaustividad %d)..."
        % (len(tareas), EXHAUSTIVIDAD))
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=HILOS) as ex:
        list(ex.map(_dock, tareas))
    log("    acoplados en %.0f s" % (time.time() - t0))

    log("")
    log("    %-24s %8s %9s %9s %8s %10s"
        % ("ligando", "fondo", "antes", "despues", "cambio", "UNBOUND"))
    log("    " + "-" * 74)
    for fichero, n, n_ok, glue, unb in cambios:
        pose = os.path.join(OUTDIR, fichero + "_out.pdbqt")
        nuevo = leer_afinidad(pose)
        viejo = leer_afinidad(os.path.join(antes, fichero + "_out.pdbqt"))
        fondo = "duro" if fichero.startswith("DECH_") else "emparejado"
        if nuevo is None or viejo is None:
            log("    %-24s %8s %9s %9s   (sin pose)"
                % (fichero, fondo,
                   "%.2f" % viejo if viejo is not None else "-",
                   "%.2f" % nuevo if nuevo is not None else "-"))
            continue
        log("    %-24s %8s %9.2f %9.2f %+9.2f %10s"
            % (fichero, fondo, viejo, nuevo, nuevo - viejo,
               "%.2f" % unb if unb is not None else "-"))

    log("")
    log("ahora hay que recalcular las metricas:")
    log("    python validar_sod1_v3.py --analizar")
    log("y comparar con %s" % os.path.basename(copia))
    log("originales guardados en out/_antes_glue/ y ligands/_roto_glue/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
