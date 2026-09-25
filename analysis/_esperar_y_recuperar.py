# -*- coding: utf-8 -*-
"""Espera a que termine el re-docking masivo principal y luego lanza la recuperacion."""
import os, sys, time

CHECKPOINT = r"C:\Users\Fredy\masive-als\analysis\_redock_tdp43_masivo\checkpoint.csv"
LOG = r"C:\Users\Fredy\masive-als\analysis\_recuperacion_lanzador.log"

def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")

t0 = time.time()
last = 0
quiet = 0
log("Vigilando el pool principal (checkpoint crece o se detiene)...")
while True:
    n = sum(1 for _ in open(CHECKPOINT)) if os.path.exists(CHECKPOINT) else 0
    if n == last:
        quiet += 1
        # el pool escribe cada ~1s; si 5 min sin cambios, termino
        if quiet >= 6 and time.time() - t0 > 600:
            log("Checkpoint estable en %d lineas -> pool principal terminado" % n)
            break
    else:
        quiet = 0
        last = n
    time.sleep(60)

log("Lanzando recuperacion de timeouts...")
sys.exit(os.system("python _recuperar_timeouts_tdp43.py"))
