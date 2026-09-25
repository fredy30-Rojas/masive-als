#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Repara los ligandos con pseudo-atomos «glue» de una carpeta, en su sitio.

POR QUE
-------
Meeko, con los ajustes por defecto, no sabe cerrar los anillos de 7 eslabones en
adelante: los parte en dos ramas y pega los extremos con pseudo-atomos de
pegamento (tipos CG0 y G0), DUPLICANDO dos atomos. El proyecto ya lo sabia para
la libreria (`gpu_dock/reparar_libreria_glue.py`, 3.803 afectados), pero los
conjuntos de validacion se prepararon con los ajustes defectuosos.

Lo que hace ese fichero corrupto:
  * Vina ve el anillo abierto, con dos atomos de mas y 14 ramas de torsion en vez
    de 2.
  * La afinidad que Vina reporta es (1) + (2) + (3) - (4), donde (4) es la energia
    interna del ligando LIBRE. Con el anillo abierto esa energia sale ~+8,9
    kcal/mol (en vez de ~0) y al restarla le regala al compuesto ~5 kcal/mol de
    afinidad que no existen.
  * Demostrado en `confirmar_glue_controles.py`: los dos unicos controles de la
    validacion v4 con pseudo-atomos eran justo los dos que daban el «positivo»
    (4MQ -10,0 -> -4,8 y ZO0 -9,8 -> -5,5 al repararlos).

QUE HACE
--------
1. Recorre la carpeta y detecta los ficheros con tipos CG0/G0.
2. Los guarda aparte (carpeta `_roto_glue`) y los vuelve a preparar con
   `rigid_macrocycles=True` y una incrustacion 3D robusta, comprobando que el
   numero de atomos coincide con el SMILES.
3. Borra las poses ya acopladas de esos ligandos, si se le pide, para que se
   re-acoplen.

Uso:
  python reparar_ligandos_glue.py --ligandos validacion_SOD1_v5/ligands \
                                  --poses validacion_SOD1_v5/out
"""
import argparse
import glob
import os
import shutil
import sys

from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))

# La receta de preparacion vive en UN solo sitio: ../preparar_ligando.py.
sys.path.insert(0, os.path.dirname(BASE))
from preparar_ligando import (construir, contar_atomos, contar_fichero,   # noqa: E402
                              smi_del_fichero)


def log(m):
    print(m, flush=True)


def contar_texto(txt):
    """Compatibilidad: ahora vive en preparar_ligando.contar_atomos."""
    return contar_atomos(txt)


def reparar(origen):
    """Devuelve (texto_pdbqt, atomos, motivo) con el anillo cerrado.

    Es la receta canonica aplicada al SMILES que lleva la propia pose o el
    propio ligando en `REMARK SMILES`. Toda la logica (rigid_macrocycles,
    incrustacion robusta y comprobacion de atomos) vive en preparar_ligando.
    """
    return construir(smi_del_fichero(origen))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligandos", required=True)
    ap.add_argument("--poses", default=None,
                    help="carpeta de poses ya acopladas, para borrar las de los afectados")
    args = ap.parse_args()

    ligandos = os.path.abspath(args.ligandos)
    rotos = os.path.join(ligandos, "_roto_glue")
    os.makedirs(rotos, exist_ok=True)

    afectados, reparados, fallos = [], [], []
    for p in sorted(glob.glob(os.path.join(ligandos, "*.pdbqt"))):
        n, glue = contar_fichero(p)
        if not glue:
            continue
        afectados.append((p, n, glue))

    log("ligandos con pseudo-atomos glue: %d de %d"
        % (len(afectados), len(glob.glob(os.path.join(ligandos, "*.pdbqt")))))
    for p, n, glue in afectados:
        nombre = os.path.basename(p)
        txt, n_ok, motivo = reparar(p)
        if txt is None:
            fallos.append((nombre, motivo))
            log("   %-40s NO reparado: %s" % (nombre, motivo))
            continue
        shutil.copy(p, os.path.join(rotos, nombre))   # se conserva el original
        with open(p, "w", encoding="utf-8") as f:
            f.write(txt)
        reparados.append(nombre)
        log("   %-40s %d atomos, %d pseudo -> %d atomos (anillo cerrado)"
            % (nombre, n, glue, n_ok))
        if args.poses:
            pose = os.path.join(os.path.abspath(args.poses),
                                nombre.replace(".pdbqt", "_out.pdbqt"))
            if os.path.exists(pose):
                os.remove(pose)
                log("      pose borrada: %s" % os.path.basename(pose))

    log("")
    log("reparados: %d | no reparables: %d | originales en %s"
        % (len(reparados), len(fallos), rotos))
    if fallos:
        log("los no reparables hay que sacarlos del conjunto, no inventar quimica:")
        for nombre, motivo in fallos:
            log("   %s: %s" % (nombre, motivo))
    return 0


if __name__ == "__main__":
    sys.exit(main())
