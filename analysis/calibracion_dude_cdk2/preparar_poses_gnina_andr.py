#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara las poses de ANDR para puntuarlas con GNINA (score_only).

PARA QUE, SI CDK2 YA LO HIZO
-----------------------------
Con CDK2 el CNN de GNINA tampoco puso las dos conformaciones del cristal por
delante (puestos 29 y 55 de 56, y las ULTIMAS de 56 en afinidad clasica). Se
concluyo que el problema no era solo de Vina: la pose era discutible. Con `andr`
la pregunta es al reves, y es la buena: el ligando TES es un estrofe con 0
enlaces rotatorios, UNA sola copia en el cristal (B medio 19,2) y los dos polares
emparejados a menos de 3,5 A. Si el CNN pone ese cristal por delante de las poses
del motor, tenemos una segunda herramienta, independiente de Vina, que confirma
que la pose cristalografica es la buena. Eso convierte el control en una
validacion y no en una discordancia.

QUE SE PUNTUA
-------------
  * **la pose del cristal** (TES, 2AM9), que es la que tiene que salir buena;
  * **las 9 poses de cada variante** del control de `andr`, cuando el motor las
    haya producido.

OJO: las poses del motor AUN NO EXISTEN. El control corre detras de los bancos, y
la GPU esta ocupada. Este script se puede ejecutar dos veces: la primera prepara
solo el cristal (y sirve para comprobar que GNINA responde), y la segunda, cuando
existan las carpetas del control, anade las poses del motor. Se dice con un
numero en pantalla, no se supone.

Todo en PDBQT de UN SOLO MODELO, porque GNINA en `--score_only` da Parse error
con las lineas MODEL/ENDMDL (lo que se vio el 28 de septiembre con el kernel de
TBK1).

Salida: `_gnina_poses_andr/` con las poses, `receptor_andr.pdbqt` y
`manifiesto.csv`.

Uso:
    python preparar_poses_gnina_andr.py
"""
from __future__ import annotations

import csv
import glob
import os
import shutil

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

import control_redocking_cdk2 as C   # noqa: E402  (la receta, una sola copia)

BASE = C.BASE
PDB_CRISTAL = os.path.join(BASE, "_cache_dude", "2AM9.pdb")
MOL2 = os.path.join(BASE, "_cache_dude", "andr_mol2.mol2")
DESTINO = os.path.join(BASE, "_gnina_poses_andr")
POSES = os.path.join(DESTINO, "out")
NOMBRE_LIG = "TES"
RECEPTOR = os.path.join(BASE, "receptor_andr.pdbqt")


def coords_hetero(nombre, ruta=PDB_CRISTAL):
    """(coords, elementos) del HETATM pedido, en orden de fichero.

    El elemento se saca del NOMBRE del atomo, no de la columna 77-78: el PDB de
    RCSB lo tiene, pero el nombre es lo que no falla nunca, y es lo que ya usan
    el resto de scripts del proyecto.
    """
    from verificar_cristal_control import elemento_de
    coords, elems = [], []
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if not l.startswith("HETATM"):
            continue
        if l[17:20].strip().upper() != nombre:
            continue
        n = l[12:16].strip()
        if not n or n[0] == "H":
            continue
        coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        elems.append(elemento_de(n))
    return np.array(coords), elems


def pdbqt_del_cristal(mol, ref, coord_hetero, elems):
    """PDBQT del cristal, con las coordenadas en el orden del mol2.

    La correspondencia se MIDE (asignacion optima por elemento), no se supone: es
    la leccion del 29 de septiembre, y aqui importa mas que en CDK2 porque el
    esteroide es simetrico en el anillo steroidico y hay mas de una permutacion
    valida.
    """
    perm, rmsd = C.emparejar_orden(elems, coord_hetero, mol, ref)
    en_orden = np.zeros_like(coord_hetero)
    en_orden[perm] = coord_hetero
    m = Chem.Mol(mol)
    cf = m.GetConformer()
    for k, p in enumerate(en_orden):
        cf.SetAtomPosition(k, (float(p[0]), float(p[1]), float(p[2])))
    con_h = Chem.AddHs(m, addCoords=True)
    txt, ok, err = PDBQTWriterLegacy.write_string(
        MoleculePreparation(rigid_macrocycles=True).prepare(con_h)[0])
    if not ok:
        raise SystemExit("meeko no pudo escribir el PDBQT: %s" % err)
    return txt, rmsd


def modelos_de(ruta):
    """Cada pose de un PDBQT de salida, como lineas SIN MODEL/ENDMDL."""
    modelos, actual = [], None
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            actual = []
        elif l.startswith("ENDMDL"):
            if actual:
                modelos.append(actual)
            actual = None
        elif actual is not None:
            actual.append(l.rstrip("\n"))
    return modelos


def main():
    mol = Chem.MolFromMol2File(MOL2, removeHs=True, sanitize=True)
    if mol is None:
        raise SystemExit("no se pudo leer %s" % MOL2)
    conf = mol.GetConformer()
    ref = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])

    shutil.rmtree(DESTINO, ignore_errors=True)
    os.makedirs(POSES)

    manifiesto = [("fichero", "origen", "que es", "kcal/mol", "rmsd_entrada")]

    # --- la pose del cristal: la que tiene que salir buena ---
    coords, elems = coords_hetero(NOMBRE_LIG)
    if len(coords) != mol.GetNumAtoms():
        raise SystemExit("el HETATM %s tiene %d atomos y el mol2 %d: no son el"
                         " mismo grupo" % (NOMBRE_LIG, len(coords),
                                           mol.GetNumAtoms()))
    txt, rmsd = pdbqt_del_cristal(mol, ref, coords, elems)
    fichero = "cristal_TES.pdbqt"
    with open(os.path.join(POSES, fichero), "w", encoding="utf-8") as f:
        f.write(txt)
    manifiesto.append((fichero, "2AM9 %s" % NOMBRE_LIG,
                       "pose del cristal (esteroide, ocupacion 1,00)", "",
                       "%.4f" % rmsd))
    print("cristal: %d atomos, correspondencia con el mol2 a %.4f A"
          % (len(coords), rmsd))

    # --- las poses del motor, si el control ya ha corrido ---
    carpetas = sorted(d for d in glob.glob(os.path.join(BASE, "_control_andr", "*"))
                      if os.path.isdir(d))
    n_motor = 0
    for carpeta in carpetas:
        variante = os.path.basename(carpeta)
        salidas = sorted(glob.glob(os.path.join(carpeta, "out", "*_out.pdbqt")))
        if not salidas:
            continue
        for i, modelo in enumerate(modelos_de(salidas[0]), 1):
            f2 = "vina_%s_m%d.pdbqt" % (variante, i)
            with open(os.path.join(POSES, f2), "w", encoding="utf-8") as f:
                f.write("\n".join(modelo) + "\n")
            energia = ""
            for l in modelo:
                if "REMARK VINA RESULT" in l:
                    energia = l.split()[3]
                    break
            manifiesto.append((f2, "Vina-GPU %s" % variante,
                               "pose %d del motor" % i, energia, ""))
            n_motor += 1

    shutil.copy2(RECEPTOR, os.path.join(DESTINO, "receptor_andr.pdbqt"))
    with open(os.path.join(DESTINO, "manifiesto.csv"), "w", encoding="utf-8",
              newline="") as f:
        csv.writer(f).writerows(manifiesto)

    print("poses preparadas: %d (cristal 1 + motor %d)" % (1 + n_motor, n_motor))
    if n_motor == 0:
        print("   AVISO: el control de andr aun no ha corrido (la GPU esta")
        print("   ocupada), asi que solo esta la pose del cristal. Hay que")
        print("   volver a ejecutar este script cuando existan las carpetas")
        print("   _control_andr/*, o GNINA solo podra ver una pose y no podra")
        print("   decir si el CNN la prefiere: eso no es una segunda opinion.")
    print("receptor: receptor_andr.pdbqt")
    print("carpeta: %s" % os.path.basename(DESTINO))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
