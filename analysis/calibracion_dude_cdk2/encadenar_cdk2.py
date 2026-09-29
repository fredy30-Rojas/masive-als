#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Espera a que el banco de TBK1 suelte la GPU y lanza la submuestra de CDK2.

POR QUE EXISTE
--------------
Las dos dianas se acoplan en la MISMA tarjeta, y Vina-GPU no comparte: el banco de
TBK1 esta usando ~11,9 GB de los 12 GB de la RTX 4080. Levantar el de CDK2 a la vez
no lo reparte, lo tumba (o tumba al otro). Asi que uno tiene que ir detras del otro,
y el que va segundo no se puede quedar esperando a que alguien mire el reloj: esto
se encarga.

COMO SABE QUE TBK1 ACABO
------------------------
El lanzador de TBK1 solo reescribe `resultados_tbk1.csv` al terminar, asi que basta
con esperar a que ese fichero cambie. Se comprueban dos senales, por si una falla:
la fecha del CSV, y una linea "CSV:" nueva en `banco_tbk1.log`.

Y por si el banco se paro en vez de terminar (alguien lo mato, se cerro la sesion):
si el registro lleva mas de una hora sin moverse, tambien se lanza. Un banco muerto
no va a soltar la GPU nunca, y quedarse esperandolo es perder la noche.

QUE LANZA
---------
1. Primero, unos minutos de tarjeta para la prueba que separa BUSQUEDA de
   PUNTUACION: el control de redocking con `search_depth 128` en la misma caja
   del banco. Si la pose del cristal no aparece ni agotando la busqueda, no es
   que no la encuentre: es que no la quiere. Va antes que el banco a proposito
   (despues serian cuatro horas de espera para tres minutos de prueba).
2. Despues `lanzar_banco_cdk2.py --decoys 4740`, que es la lectura 1:10 decidida
   con Fredy: 474 activos contra 4.740 señuelos, unas cuatro horas. El banco
   entero son 27.846 señuelos y se lanza despues, sobre lo ya hecho (el lanzador
   reanuda solo).
3. Y al final, el control de `andr`, que es el que DE VERDAD decide si el AUC se
   puede leer: CDK2 no vale como control (su ligando venia dos veces al 50 % y con
   un brazo suelto), asi que se eligio otra diana. Ver `control_andr()`.

Uso:
    python encadenar_cdk2.py
    python encadenar_cdk2.py --decoys 27846     # el banco entero
    python encadenar_cdk2.py --ya                # sin esperar (uso manual)
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
TBK1 = os.path.join(os.path.dirname(BASE), "calibracion_tbk1")
CSV_TBK1 = os.path.join(TBK1, "resultados_tbk1.csv")
LOG_TBK1 = os.path.join(TBK1, "banco_tbk1.log")
LOG = os.path.join(BASE, "encadenar_cdk2.log")

ESPERA = 180            # segundos entre comprobaciones
SILENCIO = 3600         # un registro parado una hora = el banco no va a volver
TOPE = 24 * 3600        # no esperar mas de un dia


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


def lineas_csv_en_log():
    """Cuantas veces ha escrito el lanzador la linea del CSV final."""
    try:
        with open(LOG_TBK1, encoding="utf-8", errors="ignore") as f:
            return sum(1 for l in f if l.startswith("[") and "CSV:" in l)
    except OSError:
        return 0


def esperar_tbk1():
    csv0, cl0 = mtime(CSV_TBK1), lineas_csv_en_log()
    log("esperando a que el banco de TBK1 termine")
    log("   senales iniciales: CSV %s | lineas de CSV en el registro: %d"
        % (time.strftime("%H:%M:%S", time.localtime(csv0)) if csv0 else "no hay",
           cl0))
    t0 = ultimo = time.time()
    while True:
        time.sleep(ESPERA)
        ahora = time.time()
        csv1, cl1 = mtime(CSV_TBK1), lineas_csv_en_log()
        if csv1 != csv0:
            log("el CSV de TBK1 se ha reescrito: el banco ha terminado")
            return "terminado"
        if cl1 > cl0:
            log("hay una linea de CSV nueva en el registro: el banco ha terminado")
            return "terminado"
        if mtime(LOG_TBK1) != ultimo:
            ultimo = mtime(LOG_TBK1)
        elif ahora - ultimo > SILENCIO:
            log("el registro de TBK1 lleva mas de %d min sin moverse: se da por"
                " parado y se sigue" % (SILENCIO // 60))
            return "parado"
        if ahora - t0 > TOPE:
            log("mas de un dia esperando: se lanza igual")
            return "agotado"
        if int(ahora - t0) % 1800 < ESPERA:
            log("   sigue ocupada la GPU (esperando desde hace %d min)"
                % int((ahora - t0) / 60))


def paso_previo():
    """Antes del banco: la prueba que separa BUSQUEDA de PUNTUACION.

    El control falla a 3,67-6,87 A con `search_depth 20`. Eso puede ser que la
    pose no la encuentre (busqueda) o que la encuentre y la descarte (puntuacion).
    Se distinguen agotando la busqueda en la MISMA caja del banco: si con
    `search_depth 128` tampoco aparece, no es que no la encuentre: es que no la
    quiere.

    Cuesta unos minutos de tarjeta y responde la pregunta que Fredy pidio el 29
    de septiembre. Va ANTES del banco a proposito: despues serian cuatro horas de
    espera para tres minutos de prueba.
    """
    log("paso previo: busqueda agotada en la caja del banco (depth 128)")
    for caja in (24, 20):
        orden = [sys.executable, os.path.join(BASE, "control_redocking_cdk2.py"),
                 "--caja", str(caja), "--depth", "128", "--centro", "ligando"]
        with open(os.path.join(BASE, "encadenado_stdout.log"), "a",
                  encoding="utf-8", errors="ignore") as f:
            f.write("\n===== busqueda agotada caja %d =====\n" % caja)
            f.flush()
            subprocess.run(orden, cwd=BASE, stdout=f, stderr=subprocess.STDOUT)
    log("paso previo terminado (el resultado queda en barrido_redocking_cdk2.txt)")

    # Y ahora hay que REMEDIR: el control escribe en el fichero crudo, pero el
    # informe del AUC lee el corregido (el que deduce la correspondencia de atomos y
    # separa el nucleo del brazo). Sin este paso, las variantes de depth 128 no
    # aparecerian en el informe y la seccion del control contaria solo las seis.
    log("remediendo el barrido con la medida buena, para que el informe lo vea")
    with open(os.path.join(BASE, "encadenado_stdout.log"), "a",
              encoding="utf-8", errors="ignore") as f:
        f.write("\n===== remedicion del barrido =====\n")
        f.flush()
        subprocess.run([sys.executable, os.path.join(BASE,
                                                    "diagnostico_orden_poses.py"),
                        "--remedir-barrido"],
                       cwd=BASE, stdout=f, stderr=subprocess.STDOUT)
    log("remedicion hecha (barrido_redocking_cdk2_corregido.txt)")


def lanzar(decoys):
    orden = [sys.executable, os.path.join(BASE, "lanzar_banco_cdk2.py"),
             "--decoys", str(decoys)]
    log("lanzando: %s" % " ".join(orden[1:]))
    # El lanzador escribe su propio registro y su CSV; aqui solo se le deja correr.
    with open(os.path.join(BASE, "encadenado_stdout.log"), "a",
              encoding="utf-8", errors="ignore") as f:
        f.write("\n===== encadenado %s =====\n" % time.strftime("%Y-%m-%d %H:%M"))
        f.flush()
        p = subprocess.Popen(orden, cwd=BASE, stdout=f, stderr=subprocess.STDOUT)
        codigo = p.wait()
    log("el banco de CDK2 termino con codigo %s" % codigo)
    return codigo


def control_andr():
    """EL CONTROL DE VERDAD, en la diana elegida. Al final, a proposito.

    CDK2 resulto ser un control MAL ELEGIDO, no una prueba de que el embudo fallara:
    su ligando venia modelado dos veces al 50 % y el brazo del amonio no tocaba ni la
    proteina ni el agua. Para tener una lectura que signifique algo hace falta OTRO
    control, y el 29 de septiembre se eligio `andr`: su ligando TES es un esteroide
    con 0 enlaces rotatorios, una sola copia en el cristal (B medio 19,2) y los dos
    polares con pareja de proteina a menos de 3,5 A. Ver `elegir_control_dude.py` y
    `verificar_cristal_control.py`.

    Va al FINAL de la cadena, despues de las cuatro horas del banco, por dos
    razones: es lo mas barato de todo (pocas variantes de un minuto cada una), y sus
    numeros no hacen falta para el AUC, que se mide igual. Pero es el que decide si
    el AUC se puede leer o no, asi que tiene que estar SI O SI.

    Con la caja de 24 A primero (la del proyecto) y, si hace falta, el barrido.
    """
    log("control de andr: el esteroide TES contra su propio cristal (2AM9)")
    for tam, depth in ((24, 20), (24, 128)):
        orden = [sys.executable, os.path.join(BASE, "control_redocking_andr.py"),
                 "--caja", str(tam), "--depth", str(depth),
                 "--centro", "bolsillo"]
        with open(os.path.join(BASE, "encadenado_stdout.log"), "a",
                  encoding="utf-8", errors="ignore") as f:
            f.write("\n===== control de andr, caja %d depth %d =====\n"
                    % (tam, depth))
            f.flush()
            subprocess.run(orden, cwd=BASE, stdout=f, stderr=subprocess.STDOUT)
    log("control de andr terminado (barrido_redocking_andr.txt)")

    # Si con la caja del proyecto no pasa, no se abandona: se prueba el barrido
    # completo, que es lo que se hizo con CDK2 antes de concluir que el problema no
    # era el motor.
    try:
        with open(os.path.join(BASE, "barrido_redocking_andr.txt"),
                  encoding="utf-8", errors="ignore") as f:
            texto = f.read()
    except OSError:
        texto = ""
    if "PASA" not in texto:
        log("andr no pasa con la caja de 24 A: se prueba el barrido 20/22/24 x 20/32")
        with open(os.path.join(BASE, "encadenado_stdout.log"), "a",
                  encoding="utf-8", errors="ignore") as f:
            f.write("\n===== control de andr, barrido completo =====\n")
            f.flush()
            subprocess.run([sys.executable,
                            os.path.join(BASE, "control_redocking_andr.py"),
                            "--barrido", "--centro", "bolsillo"],
                           cwd=BASE, stdout=f, stderr=subprocess.STDOUT)
    else:
        log("andr PASA: no hace falta el barrido completo")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--decoys", type=int, default=4740)
    ap.add_argument("--ya", action="store_true",
                    help="no esperar a TBK1 (para lanzarlo a mano)")
    ap.add_argument("--sin-prueba", action="store_true",
                    help="saltarse el paso previo de busqueda agotada")
    ap.add_argument("--sin-andr", action="store_true",
                    help="saltarse el control de andr del final")
    args = ap.parse_args()

    if not args.ya:
        motivo = esperar_tbk1()
        log("TBK1: %s. Se deja un minuto a la tarjeta que se quede libre." % motivo)
        time.sleep(60)
    if not args.sin_prueba:
        paso_previo()
    codigo = lanzar(args.decoys)
    if not args.sin_andr:
        try:
            control_andr()
        except Exception as e:  # noqa: BLE001
            # El banco ya esta hecho: un fallo del control no debe taparlo ni
            # impedir que se escriba el log.
            log("el control de andr fallo: %r" % e)
    return codigo


if __name__ == "__main__":
    sys.exit(main())
