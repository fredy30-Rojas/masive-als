#!/usr/bin/env python
"""Vigila la preparacion de ligandos en CPU y arranca el docking en GPU.

Fredy: "cuando termine, haslo en la GPU".

La preparacion del banco de TBK1 va con 6 nucleos de CPU
(acelerar_banco_tbk1.py) y la RTX 4080 esta ociosa. Este script espera a que
la preparacion se estabilice y, cuando se queda quieta, lanza
lanzar_banco_tbk1.py, que es el que acopla en la GPU.

Se lanza en segundo plano y sobrevive:

  pythonw vigilar_y_lanzar_gpu.py

No borra nada y no relanza si el docking ya esta vivo. Si se lanza, lo deja
corriendo: lanzar_banco_tbk1.py es reanudable, asi que si se corta, la
siguiente vuelta sigue por donde iba.
"""
import glob
import os
import subprocess
import sys
import time
from datetime import datetime

BASE = os.path.dirname(os.path.abspath(__file__))
LIGANDS = os.path.join(BASE, "ligands")
PREPARADOR = "acelerar_banco_tbk1.py"
LANZADOR = os.path.join(BASE, "lanzar_banco_tbk1.py")
LOG = os.path.join(BASE, "vigilar_gpu.log")
SALIDA = os.path.join(BASE, "out")

# Estabilizacion: cuantos segundos sin escribir ficheros nuevos para dar
# por terminada la preparacion.
SEGURO_S = 240
# Cada cuanto se mira.
CADA_S = 60
# Maximo de horas que se queda esperando antes de rendirse.
MAX_ESPERA_H = 14


def log(mensaje):
    marca = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    linea = "[%s] %s" % (marca, mensaje)
    print(linea, flush=True)
    try:
        with open(LOG, "a", encoding="utf-8", errors="ignore") as fh:
            fh.write(linea + "\n")
    except OSError:
        pass


def cuenta_ligandos():
    try:
        return len(os.listdir(LIGANDS))
    except OSError:
        return 0


def preparando():
    """True si el proceso que prepara ligandos sigue vivo."""
    try:
        salida = subprocess.run(
            ["tasklist", "/FI", "IMAGENAME eq python.exe", "/FO", "CSV"],
            capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return False
    # tasklist no da la linea de comandos, asi que se mira el fichero de
    # salida: si el preparador esta vivo, motor_stdout.log o el propio
    # contador se mueven. Por eso la comprobacion fuerte es estabilizacion.
    return True


def docking_vivo():
    """True si hay un Vina-GPU usando la GPU ahora mismo."""
    try:
        r = subprocess.run(
            ["nvidia-smi", "--query-compute-apps=process_name",
             "--format=csv,noheader"],
            capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired):
        return False
    return "Vina" in (r.stdout or "")


if __name__ == "__main__":
    log("=" * 60)
    log("Vigilante GPU para el banco TBK1")
    log("=" * 60)
    log("Preparador: %s" % PREPARADOR)
    log("Lanzador:   %s" % os.path.basename(LANZADOR))
    log("Criterio: %d s sin escribir nuevos ligandos = preparacion terminada"
        % SEGURO_S)
    log("")

    inicio = time.time()
    anterior = cuenta_ligandos()
    ultimo_cambio = time.time()
    log("Ligandos ahora: %d" % anterior)

    while True:
        time.sleep(CADA_S)
        actual = cuenta_ligandos()
        ahora = time.time()

        if actual != anterior:
            log("Ligandos: %d (+%d)" % (actual, actual - anterior))
            anterior = actual
            ultimo_cambio = ahora
            continue

        if docking_vivo():
            log("Ya hay un Vina-GPU corriendo. No lanzo nada.")
            log("Vigilante terminado.")
            break

        if ahora - ultimo_cambio >= SEGURO_S:
            log("Sin escribir desde hace %d s: la preparacion ha terminado."
                % (ahora - ultimo_cambio))
            log("Ligandos listos: %d" % actual)
            break

        if ahora - inicio > MAX_ESPERA_H * 3600:
            log("Se agotaron las %d h de espera. Me rindo sin lanzar."
                % MAX_ESPERA_H)
            break

        log("... esperando (%d s quieto, faltan %d)"
            % (int(ahora - ultimo_cambio), SEGURO_S))

    if docking_vivo():
        log("Nada que hacer.")
        sys.exit(0)

    os.makedirs(SALIDA, exist_ok=True)
    hechos = len(glob.glob(os.path.join(SALIDA, "*_out.pdbqt")))
    log("Poses ya escritas: %d" % hechos)
    log("Lanzando docking en la RTX 4080...")

    try:
        proceso = subprocess.Popen(
            [sys.executable, LANZADOR],
            cwd=BASE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.STDOUT,
            creationflags=getattr(subprocess, "CREATE_NO_WINDOW", 0),
        )
    except OSError as exc:
        log("NO SE PUDO LANZAR: %s" % exc)
        sys.exit(1)

    log("Lanzo con PID %d" % proceso.pid)
    log("El docking corre en segundo plano. Vigila el progreso con:")
    log("  tail -f %s" % os.path.join(BASE, "motor_stdout.log"))
