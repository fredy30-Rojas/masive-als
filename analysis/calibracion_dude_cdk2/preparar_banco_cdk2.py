#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara el banco de ligandos CDK2 de DUD-E con la receta canonica.

QUE ES
------
El equivalente CDK2 de `acelerar_banco_tbk1.py`: lee los dos `.ism` de DUD-E
(`actives_final.ism`, 474 activos agrupados, y `decoys_final.ism`, 27.850
senuelos) y convierte cada SMILES a PDBQT con `preparar_ligando.py`, LA receta
canonica del proyecto: macrociclos rigidos, cinco semillas de incrustacion,
sales al fragmento mayor, metales apartados y cero pseudo-atomos.

Los nombres quedan `ACT_<id>` y `DEC_<id>` (con los caracteres raros de Windows
sustituidos). El etiquetado para el AUROC sale del prefijo, como en TBK1, y se
guarda ademas `banco_cdk2.csv` con (nombre, etiqueta, id) para que la lectura
no dependa de recordar nada.

Reanudable: los ficheros que ya estan bien no se tocan (`--verificar` audita).

Uso:
    python preparar_banco_cdk2.py --comprobar
    python preparar_banco_cdk2.py                 # 8 nucleos
    python preparar_banco_cdk2.py --procesos 6
    python preparar_banco_cdk2.py --verificar
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
import time
from multiprocessing import Pool

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(BASE))
LIGANDS = os.path.join(BASE, "ligands")
MANIFIESTO = os.path.join(BASE, "banco_cdk2.csv")
INFORME = os.path.join(BASE, "banco_cdk2_informe.txt")

sys.path.insert(0, RAIZ)
import preparar_ligando as PL  # noqa: E402

LEYENDO = 0


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def esta_bien(ruta):
    """Un PDBQT del banco esta bien si existe, no esta truncado y no lleva
    pseudo-atomos (la misma prueba del banco de TBK1)."""
    if not os.path.exists(ruta) or os.path.getsize(ruta) < 100:
        return False
    try:
        with open(ruta, encoding="utf-8", errors="ignore") as f:
            txt = f.read()
    except OSError:
        return False
    if "TORSDOF" not in txt:
        return False
    return PL.contar_atomos(txt)[1] == 0


def nombre_seguro(cid):
    return re.sub(r'[\\/:*?"<>|]', "_", cid)


def cargar():
    """[(nombre_fichero, etiqueta, id)] desde los dos .ism, activos primero."""
    fuera = []
    for fich, prefijo in (("actives_final.ism", "ACT_"),
                          ("decoys_final.ism", "DEC_")):
        ruta = os.path.join(BASE, fich)
        with open(ruta, encoding="utf-8", errors="ignore") as f:
            for l in f:
                partes = l.split()
                if len(partes) < 2:
                    continue
                cid = partes[1]
                fuera.append((prefijo + nombre_seguro(cid), prefijo.strip("_"),
                              cid, partes[0]))
    return fuera


def _generar(tarea):
    nombre, etiqueta, cid, smi = tarea
    destino = os.path.join(LIGANDS, nombre + ".pdbqt")
    try:
        txt, n, motivo = PL.construir(smi)
        if txt is None:
            return (nombre, "fallo", motivo)
        with open(destino, "w", encoding="utf-8") as f:
            f.write(txt)
    except Exception as e:
        return (nombre, "fallo", "%s: %s" % (type(e).__name__, e))
    return (nombre, "generado", "")


def _auditar(ruta):
    return (os.path.basename(ruta), esta_bien(ruta))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--comprobar", action="store_true",
                    help="solo contar, no tocar nada")
    ap.add_argument("--procesos", type=int, default=8)
    ap.add_argument("--verificar", action="store_true")
    args = ap.parse_args()

    tareas = cargar()
    con_act = sum(1 for t in tareas if t[1] == "ACT")
    log("banco CDK2 de DUD-E: %d activos y %d senuelos (%d en total)"
        % (con_act, len(tareas) - con_act, len(tareas)))

    os.makedirs(LIGANDS, exist_ok=True)
    pendientes = [t for t in tareas if not esta_bien(
        os.path.join(LIGANDS, t[0] + ".pdbqt"))]
    log("ya bien: %d | por generar: %d" % (len(tareas) - len(pendientes),
                                           len(pendientes)))
    if args.comprobar:
        return 0

    if args.verificar:
        with Pool(args.procesos) as p:
            res = p.map(_auditar, [os.path.join(LIGANDS, t[0] + ".pdbqt")
                                   for t in tareas])
        malos = [n for n, ok in res if not ok]
        log("verificados %d, malos %d" % (len(res), len(malos)))
        for n in malos[:20]:
            log("   MAL: %s" % n)
        return 0

    t0 = time.time()
    fallos, apartados = [], []
    with Pool(args.procesos) as p:
        hechos = 0
        for nombre, estado, motivo in p.imap_unordered(_generar, pendientes):
            hechos += 1
            if estado == "fallo":
                fallos.append((nombre, motivo))
            if hechos % 500 == 0:
                ritmo = hechos / (time.time() - t0)
                log("progreso %d/%d | %.1f mol/s | faltan ~%.0f min"
                    % (hechos, len(pendientes), ritmo,
                       (len(pendientes) - hechos) / max(ritmo, 0.01) / 60))

    # manifiesto con etiqueta, para leer el AUROC sin depender del prefijo
    nombres_ok = {t[0] for t in tareas
                  if esta_bien(os.path.join(LIGANDS, t[0] + ".pdbqt"))}
    with open(MANIFIESTO, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["nombre", "etiqueta", "id"])
        for nombre, etiqueta, cid, _smi in tareas:
            if nombre in nombres_ok:
                w.writerow([nombre, etiqueta, cid])

    with open(INFORME, "a", encoding="utf-8") as f:
        f.write("[%s] generados %d de %d | fallos %d\n"
                % (time.strftime("%Y-%m-%d %H:%M"), len(nombres_ok),
                   len(tareas), len(fallos)))
        for nombre, motivo in fallos[:30]:
            f.write("   FALLO %s: %s\n" % (nombre, motivo))
    log("terminado: %d pdbqt buenos, %d fallos | %.0f min"
        % (len(nombres_ok), len(fallos), (time.time() - t0) / 60))
    log("manifiesto: %s" % os.path.basename(MANIFIESTO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
