#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CONTROL: se vuelve a acoplar el ligando del cristal de ANDR (TES, 2AM9).

POR QUE ESTE CONTROL Y NO EL DE CDK2
------------------------------------
CDK2 resulto ser un control MAL ELEGIDO. No era una prueba de que el embudo fallara:
su ligando cristalino venia modelado DOS veces al 50 % (FAP y FCP, altloc A y B, a
1,41 A uno del otro, B medio ~47), y el brazo del amonio no tocaba ni la proteina
(4,17 A) ni el agua mas cercana (5,90 A). Con un ligando asi, el crystal no dice
cual es la geometria buena, y que el motor no la encuentre no dice nada del motor.
Ni GNINA coincidia: el CNNscore del cristal (0,230 y 0,112) quedaba en los ultimos
puestos de 56.

ANDR se eligio cribando las 102 dianas de DUD-E con criterios que se pueden MEDIR
sin acoplar (`elegir_control_dude.py`) y luego con la pericia del cristal
(`verificar_cristal_control.py`). La criba fue dura: solo 4 de 101 ligandos tienen
0 o 1 enlace rotatorio. Gano `andr`, y gana en todo lo que importaba:

  - TES es un ESTEROIDE, C19H30O2: **0 enlaces rotatorios**, 4 anillos, 21 pesados.
    No tiene ni un brazo que pueda girar, que era justo el problema de CDK2.
  - una sola copia en el cristal, sin altloc, ocupacion 1,00. El cristal no duda.
  - factor B medio 19,2, el mas bajo de las cuatro candidatas (CDK2: ~47): la
    densidad esta bien definida.
  - los 2 polares tienen N/O de proteina a menos de 3,5 A: **0 de 2 sueltos**.
  - solo 1 de 21 atomos a mas de 4,5 A de cualquier cosa: esta dentro del bolsillo.
  - el mol2 de DUD-E ES el HETATM del PDB, a 0,0000 A al emparejar por elemento:
    receptor y ligando son del mismo cristal, que es lo que hay que comprobar.

QUE MIDE
--------
Lo mismo que el control de CDK2, y por las mismas razones (ver el docstring de
`control_redocking_cdk2.py`, que esta todo.detailed ahi):

  - `bolsillo`: RMSD de la pose al cristal SIN alinear. **Es el que decide** con el
    liston de 2 A, porque responde "esta en el mismo sitio del bolsillo?".
  - `piso`: el mejor RMSD que Podria dar cualquier emparejamiento que respete el
    elemento. Si el piso supera el liston, el veredicto no depende del orden de
    atomos ni de la simetria.
  - `alineado`: GetBestRMS de RDKit, que superpone antes de medir. Es un
    diagnostico de FORMA ("se conserva la forma?"), no un criterio de sitio.

Y la correspondencia de atomos NO se supone: se deduce con la asignacion optima
exigiendo mismo elemento sobre el PDBQT de entrada. Si esa asignacion no da ~0,00 A
el control se para, porque significaria que lo acoplado no era la pose del cristal.

COMO REUTILIZA EL CODIGO DE CDK2
-------------------------------
Importa las funciones ya escritas y probadas de `control_redocking_cdk2` en vez de
copiarlas: si la forma de medir cambia, el cambio se hace una vez y los dos controles
miden igual. Lo unico que se redefine aqui son las rutas, el nombre del ligando y el
`main`. Es un poco incomodo tener que tocar los globales del modulo importado, pero
es explicito y no esconde logica duplicada.

Salida: `_control_andr/<etiqueta>/` y `barrido_redocking_andr.txt`.

Uso:
    python control_redocking_andr.py --barrido      # 20/22/24 A x depth 20/32
    python control_redocking_andr.py --caja 24 --depth 20
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np

import control_redocking_cdk2 as cdk2
from control_redocking_cdk2 import (  # noqa: F401
    LISTON_A,
    TOL_ENTRADA,
    acoplar,
    atomos_pesados_pdbqt,
    emparejar_orden,
    escribir_ligando,
    leer_ligando,
    log,
    medir,
)

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "_cache_dude")

MOL2 = os.path.join(CACHE, "andr_mol2.mol2")
PDB = os.path.join(BASE, "receptor_andr.pdb")
CAJA = os.path.join(BASE, "caja_andr.json")
RESUMEN = os.path.join(BASE, "barrido_redocking_andr.txt")
NOMBRE_PDBQT = "TES_cristal.pdbqt"
CARPETA = os.path.join(BASE, "_control_andr")

# Las funciones importadas usan los globales de su propio modulo, asi que se
# redirigen aqui. Se hace de forma explicita para que quede claro que el control de
# CDK2 no se toca al ejecutar este.
cdk2.BASE = BASE
cdk2.MOL2 = MOL2
cdk2.PDB = PDB
cdk2.CAJA = CAJA
cdk2.RESUMEN = RESUMEN


def escribir_ligando_andr(mol, carpeta):
    """Como la de CDK2 pero nombrando el PDBQT del esteroide.

    El nombre importa poco, pero si se dejara `FAP_cristal.pdbqt` en una carpeta de
    `andr` el informe diria una cosa y el fichero otra.
    """
    import shutil

    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem

    os.makedirs(carpeta, exist_ok=True)
    con_h = Chem.AddHs(Chem.Mol(mol), addCoords=True)
    txt, ok, err = PDBQTWriterLegacy.write_string(
        MoleculePreparation(rigid_macrocycles=True).prepare(con_h)[0])
    if not ok:
        raise SystemExit("meeko no pudo escribir el PDBQT: %s" % err)
    ruta = os.path.join(carpeta, NOMBRE_PDBQT)
    with open(ruta, "w", encoding="utf-8") as f:
        f.write(txt)
    return ruta


def variante(receptor, mol, ref, centro, tamano, depth, etiqueta):
    """Una variante de caja/profundidad, escribiendo su linea en el barrido.

    Logica identica a la de CDK2; se repite aqui solo para que las rutas y el
    nombre del ligando sean los de `andr`.
    """
    import shutil

    carpeta = os.path.join(CARPETA, etiqueta)
    shutil.rmtree(carpeta, ignore_errors=True)
    escribir_ligando_andr(mol, os.path.join(carpeta, "ligands"))

    # Antes de medir nada: la correspondencia, deducida del PDBQT de entrada.
    tipos, coords = atomos_pesados_pdbqt(os.path.join(carpeta, "ligands",
                                                      NOMBRE_PDBQT))
    perm, rmsd_entrada = emparejar_orden(tipos, coords, mol, ref)
    if rmsd_entrada > TOL_ENTRADA:
        raise SystemExit("lo que se acoplo no es la pose del cristal (%.3f A):"
                         " la correspondencia no vale y el control no"
                         " significaria nada" % rmsd_entrada)
    log("  entrada: %d pesados, correspondencia con el cristal a %.4f A"
        % (len(coords), rmsd_entrada))

    modelos, _ = acoplar(carpeta, receptor, centro, tamano, depth)
    filas = medir(mol, ref, modelos, perm, tipos)
    validas = [f for f in filas if f[2] is not None]
    mejor = min(validas, key=lambda f: f[2]) if validas else None
    if mejor:
        log("  pose %d: %+.2f kcal/mol | bolsillo %.2f A | piso %.2f A"
            " | alineado %.2f A | solape %.2f A | centroides %.2f A"
            % (mejor[0], mejor[1], mejor[2], mejor[3], mejor[4], mejor[5],
               mejor[6]))
        if mejor[3] is not None and mejor[3] > LISTON_A:
            log("  ni con el mejor emparejamiento de atomos posible baja de"
                " %.2f A: el veredicto no depende del orden" % mejor[3])
    else:
        log("  sin poses comparables")
    linea = "caja %2d A | depth %2d | entrada %.4f A | %s" % (
        tamano, depth, rmsd_entrada,
        ("mejor pose %d | %+.2f kcal/mol | bolsillo %.2f A | piso %.2f A"
         " | alineado %.2f A | solape %.2f | centroides %.2f | %s"
         % (mejor[0], mejor[1], mejor[2], mejor[3], mejor[4], mejor[5],
            mejor[6], "PASA" if mejor[2] < LISTON_A else "no pasa"))
        if mejor else "sin poses")
    log(linea)
    with open(RESUMEN, "a", encoding="utf-8") as f:
        f.write(linea + "\n")
    return mejor


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--caja", type=int, default=0,
                    help="tamano de caja; 0 = el de caja_andr.json (24)")
    ap.add_argument("--depth", type=int, default=20)
    ap.add_argument("--barrido", action="store_true",
                    help="prueba cajas 20/22/24 con depth 20 y 32")
    ap.add_argument("--centro", default="bolsillo",
                    help="'bolsillo' (centro de residuos <8 A) o 'ligando'")
    args = ap.parse_args()

    for p in (MOL2, PDB, CAJA):
        if not os.path.exists(p):
            raise SystemExit("falta %s: ejecuta antes preparar_receptor_andr.py"
                             % os.path.basename(p))
    caja = json.load(open(CAJA, encoding="utf-8"))
    receptor = os.path.join(BASE, caja["receptor_pdbqt"])

    mol, ref = leer_ligando()
    from rdkit.Chem import rdMolDescriptors
    log("ligando del cristal: %d atomos pesados, %d enlaces rotatorios, %d anillos"
        % (mol.GetNumAtoms(), rdMolDescriptors.CalcNumRotatableBonds(mol),
           rdMolDescriptors.CalcNumRings(mol)))
    log("diana %s | PDB %s | liston %.1f A (medida del bolsillo, sin alinear)"
        % (caja["diana"], caja["pdb"], LISTON_A))

    ligando = np.array(caja["centro_caja"], dtype=float)
    if args.centro == "bolsillo":
        prot = []
        for l in open(PDB, encoding="utf-8", errors="ignore"):
            if l.startswith("ATOM") and l[76:78].strip().upper() != "H":
                prot.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        prot = np.array(prot)
        centro = prot[np.linalg.norm(prot - ligando[None, :], axis=1) < 8.0].mean(axis=0)
        log("centro del ligando %s | centro del bolsillo %s | se separan %.2f A"
            % (np.round(ligando, 2), np.round(centro, 2),
               float(np.linalg.norm(centro - ligando))))
    else:
        centro = ligando

    os.makedirs(CARPETA, exist_ok=True)
    with open(RESUMEN, "a", encoding="utf-8") as f:
        f.write("\n=== barrido %s (centro: %s) ===\n"
                % (time.strftime("%Y-%m-%d %H:%M"), args.centro))

    if args.barrido:
        for tam in (20, 22, 24):
            for depth in (20, 32):
                variante(receptor, mol, ref, centro, tam, depth,
                         "caja%d_depth%d" % (tam, depth))
    else:
        tam = args.caja or caja["tamano_caja"]
        variante(receptor, mol, ref, centro, tam, args.depth,
                 "caja%d_depth%d" % (tam, args.depth))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
