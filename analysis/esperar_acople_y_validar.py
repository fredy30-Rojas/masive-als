#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""esperar_acople_y_validar.py — cuando TERMINA el acople del fondo duro,
lanza la validacion limpia de TDP-43 con el fondo duro dentro.

Por que existe: el acople de los ultimos ligandos grandes tarda y no se sabe
cuando acaba. En vez de estar mirando, este vigilante espera a que no quede
ningun vina.exe vivo (dos comprobaciones seguidas, para no cortar por un
hueco entre ligandos) y entonces corre `validar_diana_limpia.py --target TDP43`,
que ya lleva el fondo blando (122 señuelos) y el fondo duro (R-BIND 2.0).

No acopla nada ni toca los ligandos: solo espera y lanza.

Uso (oculto):
    wscript.exe ejecutar_bat_oculto.vbs "C:\\Users\\Fredy\\masive-als\\analysis\\acople_y_validar_tdp43.bat"
"""
import os
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "_validacion_TDP43_rbind", "out")
LOG = os.path.join(BASE, "_esperar_acople.log")
LOG_VALIDACION = os.path.join(BASE, "_validar_tdp43_limpia_out.log")
ESPERA = 30          # segundos entre comprobaciones
MAX_HORAS = 6        # tope de seguridad


def log(m):
    linea = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), m)
    print(linea, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(linea + "\n")


def vinas_vivos():
    """Cuenta los vina.exe que quedan (el acople corre con seis a la vez)."""
    try:
        r = subprocess.run(["tasklist", "/FI", "IMAGENAME eq vina.exe", "/NH"],
                           capture_output=True, text=True, timeout=60)
        # OJO: sin nadie que listar, o si tasklist falla, `stdout` puede venir
        # como None. La primera version del vigilante murio ahi el 23 sep 2026
        # (AttributeError: 'NoneType' object has no attribute 'splitlines'),
        # justo en el momento en que ya no quedaba ningun vina.
        txt = r.stdout or ""
        return sum(1 for l in txt.splitlines() if "vina.exe" in l.lower())
    except Exception as e:
        log("no se pudo preguntar por los procesos: %s" % e)
        return 1


def poses():
    if not os.path.isdir(OUT):
        return 0
    return sum(1 for n in os.listdir(OUT) if n.endswith("_out.pdbqt"))


def main():
    log("vigilante del acople arrancado; poses ahora: %d" % poses())
    quietos = 0
    t0 = time.time()
    while True:
        n = vinas_vivos()
        if n == 0:
            quietos += 1
        else:
            quietos = 0
        log("vinas vivos: %d  (poses: %d, seguidos sin nadie: %d)"
            % (n, poses(), quietos))
        if quietos >= 2:
            break
        if time.time() - t0 > MAX_HORAS * 3600:
            log("se agoto el tope de %d horas; se lanza la validacion igual" % MAX_HORAS)
            break
        time.sleep(ESPERA)

    log("el acople ya no corre; poses finales: %d" % poses())
    log("lanzando la validacion limpia de TDP-43 con el fondo duro...")
    with open(LOG_VALIDACION, "w", encoding="utf-8") as f:
        f.write("VALIDACION LIMPIA DE TDP-43 CON FONDO DURO   %s\n"
                % time.strftime("%Y-%m-%d %H:%M"))
    try:
        with open(LOG_VALIDACION, "a", encoding="utf-8") as f:
            r = subprocess.run([sys.executable, "-X", "utf8",
                                "validar_diana_limpia.py", "--target", "TDP43"],
                               cwd=BASE, stdout=f, stderr=subprocess.STDOUT,
                               timeout=6 * 3600)
        log("validacion terminada con codigo %s" % r.returncode)
    except Exception as e:
        log("la validacion fallo al lanzarse: %s" % e)
        return 1
    log("fin. Log de la validacion: %s" % LOG_VALIDACION)
    return 0


if __name__ == "__main__":
    sys.exit(main())
