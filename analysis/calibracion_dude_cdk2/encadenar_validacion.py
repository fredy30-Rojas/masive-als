#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Espera al banco de CDK2 y pasa el validador para leer el AUC de punta a punta.

QUE ESTA ENCadenado
-------------------
El orden de la noche, y por que este script espera en vez de lanzar:

    TBK1 (suelta la tarjeta)
      -> paso previo de CDK2: busqueda agotada depth 128 + remedicion
        -> BANCO DE CDK2 (1:10, unas cuatro horas)
          -> ESTE SCRIPT: validar_banco_cdk2.py

El banco escribe `resultados_cdk2.csv` SOLO cuando ha terminado, asi que la existencia
de ese fichero con un tamaño que ya no cambia es la senal de que se puede medir. No se
mide un CSV a medias: el AUC con el banco incompleto daria una lectura que no es la
que se quiere leer, y ademas el validador espera las 474 activas, que se acoplan al
final.

POR QUE NO BASTA CON ESPERAR AL CSV
----------------------------------
Un banco que se MUERE (alguien lo cerro, se fue la corriente) no escribe el CSV
nunca, y quedarse esperando es perder la noche. Se vigila tambien el registro del
lanzador: si lleva mas de una hora sin moverse, se da por parado y se lanza el
validador igual, que es la forma de perder precision y no perder el dia. Con un banco
muerto, el validador dira que faltan compuesto, y eso ya es un dato.

QUE MIDE EL VALIDADOR
---------------------
El AUC crudo y el AUC por atomo pesado, los EF al 1 % y al 5 %, el BEDROC con la base
al azar MEDIDA (no la teorica), el bootstrap de 2000, los esqueletos de Murcko, y la
comparacion con el 0,791 que publica DUD-E para CDK2 y con el techo de 0,617 de TBK1.

Y mete dentro la seccion del CONTROL, que es lo que hace que el AUC se pueda leer: sin
un control que valga, un AUC de 0,6 no dice si es la diana o el embudo. Con el
control de CDK2 al lado (que no pasa) y el de `andr` encolado, se lee con las dos
caras.

Uso:
    python encadenar_validacion.py
    python encadenar_validacion.py --ya      # medir ahora, sin esperar
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
CSV_CDK2 = os.path.join(BASE, "resultados_cdk2.csv")
LOG_LANZADOR = os.path.join(BASE, "banco_cdk2.log")
LOG = os.path.join(BASE, "encadenar_validacion.log")
INFORME = os.path.join(BASE, "informe_cdk2.txt")
HABLAR = r"C:\Users\Fredy\hablar.py"

ESPERA = 300           # cinco minutos: aqui no hay prisa, el banco dura horas
SILENCIO = 3600        # el registro quieto una hora = banco parado
TOPE = 24 * 3600
N_ACTIVAS = 474        # las activas del banco: si faltan, el AUC no es comparable


def log(m):
    linea = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), m)
    print(linea, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(linea + "\n")


def mtime(ruta):
    try:
        return os.path.getmtime(ruta)
    except OSError:
        return 0.0


def contar_lineas(ruta):
    try:
        with open(ruta, encoding="utf-8", errors="ignore") as f:
            return sum(1 for _ in f)
    except OSError:
        return 0


def hablar(texto):
    log("VOZ: %s" % texto)
    try:
        subprocess.run([sys.executable, HABLAR, texto], timeout=600,
                       capture_output=True)
    except Exception as e:  # noqa: BLE001
        log("no se pudo hablar: %r" % e)


def esperar():
    """Hasta que el CSV exista y se haya estabilizado (o el banco muera)."""
    log("esperando a que el banco de CDK2 termine (senal: resultados_cdk2.csv)")
    log("   %d activas esperadas" % N_ACTIVAS)
    t0 = ultimo_csv = 0.0
    ultimo_log = mtime(LOG_LANZADOR)
    avisado = False
    while True:
        time.sleep(ESPERA)
        ahora = time.time()
        if os.path.exists(CSV_CDK2):
            mc = mtime(CSV_CDK2)
            if mc != ultimo_csv:
                if ultimo_csv:
                    log("   el CSV se ha reescrito: %d lineas" % contar_lineas(CSV_CDK2))
                ultimo_csv = mc
            elif ultimo_csv and ahora - ultimo_csv > 600:
                # Diez minutos quieto con el CSV ya escrito: el lanzador ya no lo
                # esta tocando, o sea que ha terminado.
                log("el CSV lleva 10 min sin moverse: el banco ha terminado")
                return "terminado"
        ml = mtime(LOG_LANZADOR)
        if ml != ultimo_log:
            ultimo_log = ml
        elif ml and ahora - ultimo_log > SILENCIO:
            log("el registro del lanzador lleva mas de 1 h quieto: se da el banco"
                " por parado y se mide igual")
            return "parado"
        if ahora - t0 > TOPE:
            log("mas de un dia esperando: se mide igual")
            return "agotado"
        if not avisado and int(ahora - t0) % 1800 < ESPERA:
            log("   el banco sigue (%d lineas de momento). esperando desde hace"
                " %d min" % (contar_lineas(CSV_CDK2), int((ahora - t0) / 60)))
            avisado = True


def leer_informe():
    """Las cifras que se dicen por voz, leidas del informe y NUNCA inventadas.

    Los patrones se han ajustado al formato que escribe `validar_banco_cdk2.py`:

        banco: 474 activos con pose y 4740 decoys con pose (10.0 por activo)
        * crudo:
          AUC 0.612 (bootstrap 95 %: 0.598 - 0.626) -> NO PASA
          EF1 % = 3.10 | EF5 % = 1.80 | BEDROC(20) = 0.181
        * por atomo pesado:
          AUC 0.625 (...)

    Se cogen las dos AUC por su bloque ("* crudo:" y "* por atomo pesado:") y no por
    la primera coincidencia, porque "AUC" aparece tambien en el texto del control y en
    la comparacion con el 0,791 publicado, y leer el numero equivocado por voz seria
    peor que no decir nada.
    """
    if not os.path.exists(INFORME):
        return None
    with open(INFORME, encoding="utf-8", errors="ignore") as f:
        txt = f.read()

    d = {}
    m = re.search(r"banco:\s*(\d+)\s*activos con pose y\s*(\d+)\s*decoys", txt)
    if m:
        d["act"], d["dec"] = int(m.group(1)), int(m.group(2))

    def bloque(etiqueta, sig):
        m = re.search(r"\*\s*%s\s*:(.*?)%s" % (re.escape(etiqueta), sig), txt,
                      re.S | re.I)
        return m.group(1) if m else ""

    def metricas(bloq):
        out = {}
        for clave, patron in (("auc", r"AUC\s+([\d.]+)\s*\("),
                              ("lo", r"bootstrap 95\s*%:\s*([\d.]+)"),
                              ("hi", r"-\s*([\d.]+)\s*\)"),
                              ("ef1", r"EF1\s*%\s*=\s*([\d.]+)"),
                              ("ef5", r"EF5\s*%\s*=\s*([\d.]+)"),
                              ("bedroc", r"BEDROC\(20\)\s*=\s*([\d.]+)")):
            mm = re.search(patron, bloq, re.I)
            if mm:
                try:
                    out[clave] = float(mm.group(1))
                except ValueError:
                    pass
        out["pasa"] = "PASA" in bloq and "NO PASA" not in bloq
        return out

    crudo = metricas(bloque("crudo", r"\*\s*por\s*atomo"))
    atomo = metricas(bloque("por atomo pesado", r"\n[A-ZÁÉÍÓÚÑ#]"))
    if not atomo:
        atomo = metricas(bloque("por atomo pesado", r"\Z"))
    d["crudo"], d["atomo"] = crudo, atomo
    d["_texto"] = txt
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ya", action="store_true",
                    help="no esperar: pasar el validador ahora mismo")
    args = ap.parse_args()

    if not args.ya:
        motivo = esperar()
        log("motivo: %s" % motivo)

    log("pasando el validador para leer el AUC de punta a punta")
    orden = [sys.executable, os.path.join(BASE, "validar_banco_cdk2.py")]
    with open(os.path.join(BASE, "validacion_stdout.log"), "a",
              encoding="utf-8", errors="ignore") as f:
        f.write("\n===== validacion %s =====\n" % time.strftime("%Y-%m-%d %H:%M"))
        f.flush()
        codigo = subprocess.run(orden, cwd=BASE, stdout=f,
                                stderr=subprocess.STDOUT).returncode
    log("el validador termino con codigo %s" % codigo)

    d = leer_informe()
    if d is None:
        hablar("Oiga don Fredy, el banco de C D K dos ha terminado, pero el validador"
               " no ha dejado informe. Conviene mirarlo a mano.")
        return codigo or 1

    crudo, atomo = d.get("crudo") or {}, d.get("atomo") or {}
    log("banco: %s activos / %s decoys" % (d.get("act", "?"), d.get("dec", "?")))
    log("crudo %s | por atomo %s" % (crudo, atomo))

    # Solo se dicen numeros que el informe trae de verdad. Si el parser no los
    # encuentra, se avisa de que el informe esta y se deja que Fredy lo lea: decir un
    # cero o un nueve por la boca seria peor que callarse.
    if "auc" not in crudo or "auc" not in atomo:
        hablar("Listo don Fredy, el banco de C D K dos ha terminado y el validador ha"
               " escrito el informe, pero no he conseguido leer las cifras del area"
               " bajo la curva y prefiero no inventarle un numero. El informe esta"
               " hecho, puede usted mirarlo.")
        log("informe: %s (cifras no leidas: hay que mirarlas a mano)"
            % os.path.basename(INFORME))
        return codigo

    veredicto = ("el intervalo de confianza queda por encima del azar, asi que el"
                 " AUC se lee como bueno"
                 if crudo.get("pasa") else
                 "el intervalo toca el azar, asi que el AUC se lee como debil")
    hablar("Listo don Fredy, el banco de C D K dos ha terminado y el validador ya ha"
           " leido el AUC de punta a punta. Con %d activos y %d senuelos, el area bajo"
           " la curva sale a %.3f, y por atomo pesado a %.3f. El enriquecimiento al"
           " uno por ciento es de %.2f. Y %s. Ya lo tiene todo en el informe, con la"
           " seccion del control al lado, que es lo que hace que el numero se pueda"
           " leer."
           % (d.get("act", 0), d.get("dec", 0), crudo["auc"], atomo["auc"],
              crudo.get("ef1", 0.0), veredicto))
    log("informe: %s" % os.path.basename(INFORME))
    return codigo


if __name__ == "__main__":
    sys.exit(main())
