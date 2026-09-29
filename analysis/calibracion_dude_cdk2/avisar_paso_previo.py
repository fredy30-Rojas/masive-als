#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Avisa a Fredy cuando el paso previo del encadenador de CDK2 haya terminado.

QUE ESTA ESPERANDO
-----------------
Fredy pidio el 29 de septiembre que le avisaran en cuanto acabara el PASO PREVIO del
encadenador, para leer el resultado de la BUSQUEDA AGOTADA. El paso previo es la
prueba que separa BUSQUEDA de PUNTUACION: acopla el ligando del cristal de CDK2 con
`search_depth 128` en la misma caja del banco. Si la pose del cristal tampoco aparece
agotando la busqueda, no es que el motor no la encuentre: es que no la quiere.

COMO SABE QUE TERMINO
---------------------
Mirando `encadenar_cdk2.log`. El paso previo acaba justo despues de la REMEDICION del
barrido, asi que la senal buena es la linea

    remedicion hecha (barrido_redocking_cdk2_corregido.txt)

que es lo ULTIMO que hace el paso previo antes de que arranque el banco de cuatro
horas. Cuando aparece, se lee el barrido remedido y se avisa por voz.

Si el paso previo falla y el proceso se muere, tambien hay que avisar: si no, Fredy se
queda esperando toda la noche un aviso que no va a llegar. Se controla mirando si el
fichero de lock del proceso sigue vivo, o, mas simple, si el log lleva mucho tiempo
quieto Y no ha aparecido la senal.

QUE DICE EL AVISO
-----------------
El resultado de las variantes de depth 128, que son las que responden a la pregunta.
El numero que decide es el `bolsillo` (RMSD sin alinear) con liston de 2 A. El `piso`
es la cota inferior que no depende del orden de atomos.

Lo que se lee es `barrido_redocking_cdk2_corregido.txt` a proposito, no el crudo: el
crudo no lleva la separacion de la medida del bolsillo y el piso, que es la que vale.

Uso:
    python avisar_paso_previo.py
"""
from __future__ import annotations

import os
import re
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
LOG = os.path.join(BASE, "encadenar_cdk2.log")
BARRIDO = os.path.join(BASE, "barrido_redocking_cdk2_corregido.txt")
AVISO = os.path.join(BASE, "aviso_paso_previo.log")
HABLAR = r"C:\Users\Fredy\hablar.py"

ESPERA = 60
SILENCIO = 5400        # el log quieto 90 min sin la senal = algo se rompio
TOPE = 20 * 3600
LISTON_A = 2.0         # el liston del proyecto, el mismo que usa el validador


def log(m):
    linea = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), m)
    print(linea, flush=True)
    with open(AVISO, "a", encoding="utf-8") as f:
        f.write(linea + "\n")


def leer_log():
    try:
        with open(LOG, encoding="utf-8", errors="ignore") as f:
            return f.read()
    except OSError:
        return ""


def variantes_agotadas():
    """Las del barrido remedido, y solo las de depth alto (la busqueda agotada)."""
    filas = []
    try:
        with open(BARRIDO, encoding="utf-8", errors="ignore") as f:
            for l in f:
                if "depth" not in l or "|" not in l:
                    continue
                m = re.search(r"caja\s*(\d+)\s*A\s*\|\s*depth\s*(\d+)", l)
                if not m or int(m.group(2)) < 100:
                    continue
                b = re.search(r"bolsillo\s*([\d.]+)", l)
                p = re.search(r"piso\s*([\d.]+)", l)
                c = re.search(r"centroides\s*([\d.]+)", l)
                v = "PASA" if "PASA" in l and "no pasa" not in l else "no pasa"
                filas.append({"caja": int(m.group(1)),
                              "depth": int(m.group(2)),
                              "bolsillo": float(b.group(1)) if b else None,
                              "piso": float(p.group(1)) if p else None,
                              "solape": (re.search(r"solape\s*([\d.]+)",
                                                   l).group(1)
                                         if re.search(r"solape\s*([\d.]+)", l)
                                         else None),
                              "veredicto": v})
    except OSError:
        pass
    return filas


def hablar(texto):
    """Habla por TTS con la voz del proyecto.

    Si hablar.py falla, el aviso queda igualmente en el log: el aviso nunca depende
    solo de la voz.
    """
    log("VOZ: %s" % texto)
    try:
        subprocess.run([sys.executable, HABLAR, texto], timeout=300,
                       capture_output=True)
    except Exception as e:  # noqa: BLE001
        log("no se pudo hablar: %r" % e)


def resumen_voz(filas):
    """El aviso, en una frase. Sin listas: lo que se oye es lo que se dice."""
    if not filas:
        return ("La busqueda agotada ha terminado, pero no han quedado variantes de "
                "depth alto en el barrido. Hay que mirarlo a mano.")
    con_dato = [f for f in filas if f["bolsillo"] is not None]
    if not con_dato:
        return ("La busqueda agotada ha terminado con %d variantes, pero el "
                "barrido no trae el RMSD del bolsillo. Hay que mirarlo a mano."
                % len(filas))
    mejor = min(con_dato, key=lambda f: f["bolsillo"])
    n_pasa = sum(1 for f in con_dato if f["bolsillo"] < LISTON_A)
    if n_pasa:
        return ("Listo don Fredy, la busqueda agotada ha terminado y si que pasa: "
                "de %d variantes de depth alto, %d %s por debajo de dos "
                "angstroms. La mejor esta a %.2f. Ya puede leerlo."
                % (len(con_dato), n_pasa,
                   "caen" if n_pasa > 1 else "cae", mejor["bolsillo"]))
    return ("Listo don Fredy, la busqueda agotada ha terminado y no pasa ninguna de "
            "las %d variantes de depth alto. La mejor se queda a %.2f angstroms en el "
            "bolsillo, y el piso a %.2f, asi que no es un problema de como se mida. "
            "Con depth ciento veintiocho el motor no encuentra la pose del cristal: "
            "no es que no la encuentre, es que no la quiere. Ya puede leerlo."
            % (len(con_dato), mejor["bolsillo"],
               mejor["piso"] if mejor["piso"] is not None else -1))


def main():
    log("esperando a que el paso previo de CDK2 termine (senal: 'remedicion hecha')")
    t0 = ultimo = time.time()
    visto = len(leer_log())
    avisado = False
    while True:
        time.sleep(ESPERA)
        ahora = time.time()
        texto = leer_log()
        if len(texto) != visto:
            visto = len(texto)
            ultimo = ahora
        if "remedicion hecha" in texto:
            log("SENAL: el paso previo ha terminado (remediacion hecha)")
            # Un momento de margen: el fichero se escribe justo antes del aviso.
            time.sleep(10)
            filas = variantes_agotadas()
            log("variantes de busqueda agotada encontradas: %d" % len(filas))
            for f in filas:
                log("   caja %d depth %d | bolsillo %s | piso %s | %s"
                    % (f["caja"], f["depth"], f["bolsillo"], f["piso"],
                       f["veredicto"]))
            hablar(resumen_voz(filas))
            log("aviso dado. El resultado esta en %s" % os.path.basename(BARRIDO))
            return 0
        if ahora - ultimo > SILENCIO:
            log("AVISO: el log lleva %d min quieto y no ha aparecido la senal."
                " Puede que el paso previo haya fallado."
                % (SILENCIO // 60))
            hablar("Oiga don Fredy, el encadenador lleva hora y media sin moverse y "
                   "el paso previo no ha dado ninguna senal de haber terminado. "
                   "Puede que algo se haya roto, o que este esperando todavia a que "
                   "T B K uno suelte la tarjeta. Conviene mirarlo.")
            return 1
        if ahora - t0 > TOPE:
            log("mas de 20 h esperando sin senal: se avisa igual")
            hablar("Oiga don Fredy, han pasado muchas horas y el paso previo sigue "
                   "sin dar senal de haber terminado.")
            return 1
        if not avisado and int(ahora - t0) % 1800 < ESPERA:
            log("   sigue esperando (lleva %d min)" % int((ahora - t0) / 60))
            avisado = True


if __name__ == "__main__":
    sys.exit(main())
