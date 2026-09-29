#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara las poses de CDK2 para puntuarlas con GNINA (score_only).

QUE SE PUNTUA
-------------
La pregunta es si el CNN de GNINA distingue la pose del cristal de las que devuelve
Vina-GPU. Para eso hay que meter en el mismo saco, y con el mismo formato:

  * **las dos conformaciones del cristal**: FAP (el mol2 de DUD-E, altloc A) y FCP
    (altloc B). Ver `cristal_1h00_dos_conformaciones.py`: son el mismo ligando al
    50 % cada una, a 1,41 A;
  * **las 9 poses de cada una de las 6 variantes** del barrido, que son las que
    devolvio el motor.

Todo en PDBQT de UN SOLO MODELO, porque GNINA en modo `--score_only` da Parse error
con las lineas MODEL/ENDMDL (fallo resuelto el 28 de septiembre con el kernel de
TBK1: `_kaggle/masive-als-gnina-score.ipynb`, celda 3).

Salida: `_gnina_poses/` con las poses, `receptor_cdk2.pdbqt`, y `manifiesto.csv` que
dice de donde sale cada fichero. Ese directorio es el que se sube a Kaggle.

Uso:
    python preparar_poses_gnina.py
"""
from __future__ import annotations

import glob
import os
import shutil

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

import control_redocking_cdk2 as C   # noqa: E402  (la receta, una sola copia)

BASE = C.BASE
PDB = os.path.join(BASE, "1H00_cristal_completo.pdb")
DESTINO = os.path.join(BASE, "_gnina_poses")
POSES = os.path.join(DESTINO, "out")


def coords_hetero(nombre):
    """(coords, elementos) de un HETATM del PDB original, en orden de fichero."""
    coords, elems = [], []
    for l in open(PDB, encoding="utf-8", errors="ignore"):
        if not l.startswith("HETATM") or l[17:20].strip() != nombre:
            continue
        coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        elems.append(l[76:78].strip().upper())
    return np.array(coords), elems


def pdbqt_del_cristal(mol, ref, coord_hetero, elems):
    """PDBQT de una conformacion del cristal, con la receta del proyecto.

    Las coordenadas llegan en el orden del PDB y hay que ponerlas en el orden del
    mol2 antes de escribir: la correspondencia se mide, no se supone (la leccion del
    29 de septiembre, ver `diagnostico_orden_poses.py`).
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
    """Las poses de un PDBQT de salida, cada una como lista de lineas SIN MODEL/ENDMDL."""
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
    mol, ref = C.leer_ligando()
    shutil.rmtree(DESTINO, ignore_errors=True)
    os.makedirs(POSES)

    manifiesto = [("fichero", "origen", "que es", "kcal/mol", "rmsd_bolsillo")]

    # --- las dos conformaciones del cristal ---
    for nombre, etiqueta in (("FAP", "cristal_FAP_altlocA"),
                             ("FCP", "cristal_FCP_altlocB")):
        coords, elems = coords_hetero(nombre)
        if len(coords) != 30:
            raise SystemExit("se esperaban 30 atomos en %s; hay %d"
                             % (nombre, len(coords)))
        txt, rmsd = pdbqt_del_cristal(mol, ref, coords, elems)
        fichero = "%s.pdbqt" % etiqueta
        with open(os.path.join(POSES, fichero), "w", encoding="utf-8") as f:
            f.write(txt)
        manifiesto.append((fichero, "1H00 %s" % nombre,
                           "conformacion del cristal (%s)" % nombre, "",
                           "%.3f" % rmsd))

    # --- las poses del motor, un fichero por pose ---
    for carpeta in sorted(d for d in glob.glob(os.path.join(BASE, "_control_cdk2",
                                                            "*"))
                          if os.path.isdir(d)):
        variante = os.path.basename(carpeta)
        salidas = sorted(glob.glob(os.path.join(carpeta, "out", "*_out.pdbqt")))
        if not salidas:
            continue
        for i, modelo in enumerate(modelos_de(salidas[0]), 1):
            fichero = "vina_%s_m%d.pdbqt" % (variante, i)
            with open(os.path.join(POSES, fichero), "w", encoding="utf-8") as f:
                f.write("\n".join(modelo) + "\n")
            energia = ""
            for l in modelo:
                if "REMARK VINA RESULT" in l:
                    energia = l.split()[3]
                    break
            manifiesto.append((fichero, "Vina-GPU %s" % variante,
                               "pose %d del motor" % i, energia, ""))

    shutil.copy2(os.path.join(BASE, "receptor_cdk2.pdbqt"),
                 os.path.join(DESTINO, "receptor_cdk2.pdbqt"))

    with open(os.path.join(DESTINO, "manifiesto.csv"), "w", encoding="utf-8",
              newline="") as f:
        for fila in manifiesto:
            f.write(",".join(str(x) for x in fila) + "\n")

    n_poses = len(glob.glob(os.path.join(POSES, "*.pdbqt")))
    print("poses preparadas: %d" % n_poses)
    print("   del cristal: %d" % len(glob.glob(os.path.join(POSES, "cristal_*"))))
    print("   del motor:   %d" % len(glob.glob(os.path.join(POSES, "vina_*"))))
    print("receptor: receptor_cdk2.pdbqt")
    print("carpeta para Kaggle: %s" % os.path.relpath(DESTINO, C.RAIZ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
