#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Cosecha el kernel de Kaggle del banco de CDK2 y deja montada la sesion siguiente.

POR QUE EXISTE
--------------
El 4 de octubre de 2026 el banco de CDK2 (28.301 ligandos) se llevo a Kaggle
porque la RTX 4080 iba llena con las 400 tandas de ZINC. El kernel acopla las
23.087 pendientes y escribe su CSV en /kaggle/working.

El relevo estaba contado en `KAGGLE_CDK2_2026-10-04.md` como una receta de tres
pasos, y ese es justo el problema: una receta en prosa se ejecuta mal cuando ya
ha pasado una noche y el que la lee no es el que la escribio. Aqui queda en
codigo, con las guardas en su sitio:

  1. baja la salida del kernel (`kaggle kernels output`);
  2. fusiona su CSV con el de casa, sin pisar lo ya acoplado;
  3. rehace `banco_cdk2.tar.gz` con el CSV nuevo dentro y sube una version del
     dataset (`kaggle datasets version`);
  4. vuelve a empujar el kernel (`kaggle kernels push`), que reanuda solo lo
     que falte.

QUE NO HACE
-----------
No toca la RTX 4080 ni el kernel mientras corre. Si el kernel esta RUNNING, se
niega y lo dice: bajar la salida a medias solo confunde. La guarda de estado
existe para que nadie confunda «ha terminado» con «va por la mitad», que es el
unico error que aqui cuesta cuota.

EL CSV
------
Tanto el de casa como el de Kaggle tienen las MISMAS columnas, porque el kernel
reproduce el protocolo local: `ligand,target,energy`. La fusion es una union por
`ligand`. En los pocos nombres que aparezcan en los dos (si una sesion se repitio),
gana el que ya estaba en casa: da igual el valor porque el protocolo es identico,
y quedarse con el de casa deja el fichero estable y comparable.

Uso:
    python cosechar_cdk2_kaggle.py                 # solo mirar el estado y contar
    python cosechar_cdk2_kaggle.py --bajar         # baja, fusiona y rehace el tar
    python cosechar_cdk2_kaggle.py --bajar --relanzar   # y ademas sube y relanza
    python cosechar_cdk2_kaggle.py --bajar --forzar     # aunque el kernel no haya
                                                        # terminado (a mano, a
                                                        # sabiendas)
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import shutil
import subprocess
import sys
import tarfile
import time

BASE = os.path.dirname(os.path.abspath(__file__))
RUN = os.path.join(BASE, "_kaggle_cdk2_run")
BANCO_DIR = os.path.join(BASE, "_kaggle_banco")
TAR = os.path.join(BANCO_DIR, "banco_cdk2.tar.gz")
LIGANDS = os.path.join(BASE, "ligands")
RECEPTOR = os.path.join(BASE, "receptor_cdk2.pdbqt")
CAJA = os.path.join(BASE, "caja_cdk2.json")
CSV_CASA = os.path.join(BASE, "resultados_cdk2.csv")
CSV_RESPALDO = os.path.join(BASE, "resultados_cdk2_pre_kaggle.csv")
BAJADO = os.path.join(BASE, "_bajado_kaggle")

KERNEL = "yograbotodo/masive-als-cdk2-banco-gpu"
DATASET = "yograbotodo/masive-als-cdk2-banco"
NOMBRE_CSV_KERNEL = "resultados_cdk2_kaggle.csv"

LOG = os.path.join(BASE, "cosechar_cdk2_kaggle.log")


def log(m):
    linea = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), m)
    print(linea, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(linea + "\n")


def kaggle(*args):
    """Llama al CLI como modulo; es como esta instalado en este PC."""
    orden = [sys.executable, "-m", "kaggle", *args]
    return subprocess.run(orden, capture_output=True, text=True)


def estado():
    r = kaggle("kernels", "status", KERNEL)
    if r.returncode != 0:
        raise SystemExit("no se pudo preguntar el estado: %s"
                         % (r.stderr.strip() or r.stdout.strip()))
    salida = (r.stdout + r.stderr).strip()
    for etiqueta in ("COMPLETE", "ERROR", "CANCEL", "RUNNING", "QUEUE"):
        if etiqueta in salida.upper():
            return etiqueta
    return "DESCONOCIDO"


def contar_casa():
    """(hechos, ligandos, pendientes) del banco local."""
    if not os.path.exists(CSV_CASA):
        raise SystemExit("falta %s" % os.path.relpath(CSV_CASA, BASE))
    hechos = set()
    for r in csv.DictReader(open(CSV_CASA, encoding="utf-8")):
        if r.get("ligand") and r.get("energy"):
            hechos.add(r["ligand"].strip())
    ligandos = [os.path.basename(p)[:-len(".pdbqt")]
                for p in glob.glob(os.path.join(LIGANDS, "*.pdbqt"))]
    pendientes = [l for l in ligandos if l not in hechos]
    return hechos, ligandos, pendientes


def bajar():
    """Baja /kaggle/working del kernel a _bajado_kaggle y devuelve su CSV."""
    shutil.rmtree(BAJADO, ignore_errors=True)
    os.makedirs(BAJADO)
    log("bajando la salida del kernel a %s" % os.path.relpath(BAJADO, BASE))
    r = kaggle("kernels", "output", KERNEL, "-p", BAJADO)
    if r.returncode != 0:
        raise SystemExit("kaggle kernels output fallo: %s"
                         % (r.stderr.strip() or r.stdout.strip()))
    hits = glob.glob(os.path.join(BAJADO, "**", NOMBRE_CSV_KERNEL),
                     recursive=True)
    if not hits:
        # Kaggle a veces mete todo bajo un subdirectorio con el nombre del
        # kernel; el glob recursivo ya lo cubre, pero si el CSV se llama de otra
        # forma conviene decir que se ha bajado en vez de morir en silencio.
        listado = glob.glob(os.path.join(BAJADO, "**", "*"), recursive=True)
        raise SystemExit("no encuentro %s en la salida; esto es lo que hay: %s"
                         % (NOMBRE_CSV_KERNEL,
                            ", ".join(os.path.relpath(p, BAJADO)
                                      for p in listado[:20]) or "(vacio)"))
    log("CSV de Kaggle: %s" % os.path.relpath(hits[0], BASE))
    return hits[0]


def fusionar(csv_kaggle):
    """Une el CSV de Kaggle sobre el de casa. Devuelve (añadidos, total)."""
    hechos, ligandos, _ = contar_casa()
    respaldo = not os.path.exists(CSV_RESPALDO)
    if respaldo:
        # Se guarda UNA vez: el respaldo es el estado antes de tocar nada con
        # Kaggle. Volver a copiarlo en cada sesion borraria el original.
        shutil.copy2(CSV_CASA, CSV_RESPALDO)
        log("respaldo del CSV de casa: %s"
            % os.path.relpath(CSV_RESPALDO, BASE))

    filas, vistos = [], set()
    for r in csv.DictReader(open(CSV_CASA, encoding="utf-8")):
        lig = (r.get("ligand") or "").strip()
        if not lig or not r.get("energy"):
            continue
        filas.append({"ligand": lig,
                      "target": r.get("target") or "CDK2",
                      "energy": r["energy"]})
        vistos.add(lig)

    añadidos = 0
    for r in csv.DictReader(open(csv_kaggle, encoding="utf-8")):
        lig = (r.get("ligand") or "").strip()
        if not lig or not r.get("energy") or lig in vistos:
            continue
        filas.append({"ligand": lig,
                      "target": r.get("target") or "CDK2",
                      "energy": r["energy"]})
        vistos.add(lig)
        añadidos += 1

    tmp = CSV_CASA + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ligand", "target", "energy"])
        w.writeheader()
        w.writerows(sorted(filas, key=lambda r: r["ligand"]))
    os.replace(tmp, CSV_CASA)
    log("fusion: +%d de Kaggle | total con energia: %d de %d"
        % (añadidos, len(vistos), len(ligandos)))
    if añadidos == 0:
        log("AVISO: Kaggle no aporto ningun ligando nuevo. O el kernel no llego"
            " a acoplar, o su CSV venia sin energia. Revisar antes de relanzar.")
    return añadidos, len(vistos)


def rehacer_tar():
    """Reconstruye banco_cdk2.tar.gz con el CSV fusionado dentro.

    La estructura es la que el kernel espera al descomprimir: ligands/,
    receptor_cdk2.pdbqt, caja_cdk2.json y resultados_cdk2.csv en la raiz.
    """
    for necesario in (LIGANDS, RECEPTOR, CAJA, CSV_CASA):
        if not os.path.exists(necesario):
            raise SystemExit("falta %s para rehacer el tar"
                             % os.path.relpath(necesario, BASE))
    tmp = TAR + ".tmp"
    log("rehaciendo %s" % os.path.relpath(TAR, BASE))
    with tarfile.open(tmp, "w:gz") as t:
        t.add(LIGANDS, arcname="ligands")
        t.add(RECEPTOR, arcname="receptor_cdk2.pdbqt")
        t.add(CAJA, arcname="caja_cdk2.json")
        t.add(CSV_CASA, arcname="resultados_cdk2.csv")
    os.replace(tmp, TAR)
    log("tar rehecho: %.1f MB" % (os.path.getsize(TAR) / 1e6))


def subir_dataset(csv_nuevos):
    log("subiendo version nueva del dataset %s" % DATASET)
    r = kaggle("datasets", "version", "-p", BANCO_DIR,
               "-m", "CSV de CDK2 con %d acoplados (sesion Kaggle %s)"
               % (csv_nuevos, time.strftime("%Y-%m-%d")))
    if r.returncode != 0:
        raise SystemExit("kaggle datasets version fallo: %s"
                         % (r.stderr.strip() or r.stdout.strip()))
    log("   %s" % (r.stdout.strip() or r.stderr.strip()))
    # Kaggle tarda en procesar; si se empuja el kernel antes de que la version
    # este lista, monta la ANTERIOR y repite trabajo. Se espera a que diga ready.
    for _ in range(30):
        time.sleep(20)
        s = kaggle("datasets", "status", DATASET)
        txt = (s.stdout + s.stderr).strip()
        log("   estado del dataset: %s" % txt)
        if "ready" in txt.lower():
            return
    log("AVISO: el dataset no dijo 'ready' en 10 min; revisar antes de relanzar")


def relanzar():
    log("empujando el kernel %s" % KERNEL)
    r = kaggle("kernels", "push", "-p", RUN)
    if r.returncode != 0:
        raise SystemExit("kaggle kernels push fallo: %s"
                         % (r.stderr.strip() or r.stdout.strip()))
    log("   %s" % (r.stdout.strip() or r.stderr.strip()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bajar", action="store_true",
                    help="baja y fusiona la salida del kernel")
    ap.add_argument("--relanzar", action="store_true",
                    help="sube el dataset y vuelve a empujar el kernel")
    ap.add_argument("--sin-tar", action="store_true",
                    help="fusiona sin rehacer el tar (para mirar antes)")
    ap.add_argument("--forzar", action="store_true",
                    help="no exigir que el kernel haya terminado")
    args = ap.parse_args()

    e = estado()
    hechos, ligandos, pendientes = contar_casa()
    log("kernel %s: %s" % (KERNEL, e))
    log("banco: %d ligandos | hechos %d | pendientes %d"
        % (len(ligandos), len(hechos), len(pendientes)))

    if not args.bajar:
        log("nada que hacer (sin --bajar). Solo se informo del estado.")
        return 0

    if e != "COMPLETE" and not args.forzar:
        raise SystemExit("el kernel esta %s, no COMPLETE: no se baja a medias."
                         " Esperar a que termine, o --forzar si es a proposito."
                         % e)

    csv_kaggle = bajar()
    añadidos, total = fusionar(csv_kaggle)
    _, _, pendientes = contar_casa()
    log("pendientes tras la fusion: %d" % len(pendientes))

    if not args.sin_tar:
        rehacer_tar()
    if args.relanzar:
        if args.sin_tar:
            raise SystemExit("--relanzar exige rehacer el tar: el kernel lee el"
                             " CSV previo del dataset")
        subir_dataset(total)
        relanzar()
        log("relevo hecho: dataset actualizado y kernel relanzado")
    return 0


if __name__ == "__main__":
    sys.exit(main())
