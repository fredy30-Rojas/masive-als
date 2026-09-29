#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Espera al banco de CDK2 y luego hace el control de redocking de ANDR.

POR QUE ESTE SCRIPT Y NO UNO SOLO
---------------------------------
`encadenar_cdk2.py` ya estaba CORRIENDO cuando se decidio anadir el control de `andr`
al final de la cadena. Un proceso que ya ha cargado su codigo no lo recarga solo, asi
que el cambio no le llega. Y no se puede reiniciar a la fuerza: lleva 19:28 esperando
a que TBK1 libere la tarjeta, y matar ese proceso seria tirar la espera.

La solucion es que este script haga SOLO la parte nueva, y esperando por su cuenta a
que la tarjeta este libre. Los dos se turnan sin pisarse: este no toca nada hasta que
`resultados_cdk2.csv` ha aparecido y la GPU esta descargada.

QUE HACE, EN ORDEN
------------------
1. Espera a que exista `resultados_cdk2.csv` y a que la GPU este descargada.
2. Lanza el control de `andr` con la caja del proyecto (24 A, depth 20).
3. Repite con `search_depth 128` en la misma caja, que separa BUSQUEDA de PUNTUACION:
   si la pose del cristal tampoco aparece agotando la busqueda, no es que el motor
   no la encuentre, es que no la quiere.
4. Si con ninguna de las dos pasa el liston de 2 A, hace el barrido completo
   20/22/24 A x depth 20/32, que es lo que se hizo con CDK2 antes de concluir que el
   problema no estaba en el motor.

Uso:
    python encadenar_andr.py
    python encadenar_andr.py --ya        # sin esperar, para lanzarlo a mano
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
CSV_CDK2 = os.path.join(BASE, "resultados_cdk2.csv")
LOG = os.path.join(BASE, "encadenar_andr.log")
STDOUT = os.path.join(BASE, "encadenado_stdout.log")
CONTROL = os.path.join(BASE, "control_redocking_andr.py")
BARRIDO = os.path.join(BASE, "barrido_redocking_andr.txt")

ESPERA = 180
SILENCIO = 1800        # el CSV de CDK2 sin moverse media hora = banco parado
TOPE = 24 * 3600
GPU_MAX_UMBRAL = 2000  # MiB: por encima, la tarjeta sigue ocupada


def log(m):
    linea = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), m)
    print(linea, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(linea + "\n")


def gpu_libre():
    """MiB de VRAM en uso, leidos con nvidia-smi. -1 si no se puede leer."""
    try:
        r = subprocess.run(
            ["nvidia-smi", "--query-gpu=memory.used", "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=30)
        return int(r.stdout.strip().splitlines()[0])
    except Exception:  # noqa: BLE001
        return -1


def esperar():
    log("esperando a que el banco de CDK2 termine y la tarjeta se quede libre")
    t0 = visto = 0.0
    avisado = False
    while True:
        time.sleep(ESPERA)
        ahora = time.time()
        if os.path.exists(CSV_CDK2):
            m = os.path.getmtime(CSV_CDK2)
            if not visto:
                visto = m
                log("   resultados_cdk2.csv ya existe; esperando a que la"
                    " tarjeta quede libre")
            elif ahora - visto > SILENCIO:
                log("   el CSV lleva media hora sin moverse: se da el banco por"
                    " parado y se sigue")
                return
            elif m != visto:
                visto = m
        mem = gpu_libre()
        if mem >= 0 and mem <= GPU_MAX_UMBRAL:
            log("   tarjeta libre (%d MiB): adelante" % mem)
            return
        if not avisado and int(ahora - t0) % 1800 < ESPERA:
            log("   la tarjeta sigue ocupada (%d MiB). esperando desde hace %d min"
                % (mem, int((ahora - t0) / 60)))
            avisado = True
        if ahora - t0 > TOPE:
            log("mas de un dia esperando: se lanza igual")
            return


def correr(orden, etiqueta):
    with open(STDOUT, "a", encoding="utf-8", errors="ignore") as f:
        f.write("\n===== %s %s =====\n" % (etiqueta, time.strftime("%Y-%m-%d %H:%M")))
        f.flush()
        return subprocess.run(orden, cwd=BASE, stdout=f,
                              stderr=subprocess.STDOUT).returncode


def pasado():
    """Si el control ya paso el liston en algun momento."""
    try:
        with open(BARRIDO, encoding="utf-8", errors="ignore") as f:
            return "PASA" in f.read()
    except OSError:
        return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ya", action="store_true",
                    help="no esperar: lanzar el control ahora mismo")
    args = ap.parse_args()

    if not args.ya:
        esperar()
        time.sleep(60)   # un minuto de cortesia para que la tarjeta se asiente

    log("control de andr: el esteroide TES contra su propio cristal (2AM9)")
    for tam, depth in ((24, 20), (24, 128)):
        log("   variante: caja %d A, search_depth %d" % (tam, depth))
        correr([sys.executable, CONTROL, "--caja", str(tam),
                "--depth", str(depth), "--centro", "bolsillo"],
               "control de andr, caja %d depth %d" % (tam, depth))

    if not pasado():
        log("   no pasa con la caja de 24 A: barrido 20/22/24 x depth 20/32")
        correr([sys.executable, CONTROL, "--barrido", "--centro", "bolsillo"],
               "control de andr, barrido completo")
    else:
        log("   PASA con la caja del proyecto: no hace falta el barrido")

    log("control de andr terminado. Resultado en %s"
        % os.path.basename(BARRIDO))
    return 0


if __name__ == "__main__":
    sys.exit(main())
