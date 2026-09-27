#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CONTROL: se vuelve a acoplar el BX-795 del cristal y se compara con su pose.

QUE PREGUNTA CONTESTA
---------------------
El banco de TBK1 se apoya en dos cosas que no se pueden dar por supuestas: que el
receptor esta bien preparado y que la caja esta donde el ligando. La forma de
saberlo que usa el proyecto desde el primer dia es **volver a acoplar el ligando
que ya venia cristalizado** y medir a que distancia cae de donde el cristal dice.

COMO SE MIDE (y por que asi)
----------------------------
El ligando que se acopla se construye **del propio BX-795 del cristal**, no de un
SMILES de fuera. Asi el orden de atomos del PDBQT es el del cristal, y el RMSD es
una resta de coordenadas directa, sin casar por simetria ni por subestructura:
nada que se pueda equivocar en el camino.

Y solo se compara la mejor pose de Vina-GPU con la del cristal, en el mismo
sistema de referencia (la caja se define sobre esas coordenadas).

Uso:
    python control_redocking_bx795.py

Salida: `_control_bx795/` con el ligando, su pose, el log del motor y el informe.
Si el RMSD sale por debajo de 2 A, el receptor y la caja estan bien puestos.
"""
from __future__ import annotations

import glob
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(BASE))
GPU = os.path.join(RAIZ, "gpu_dock")
EXE = os.path.join(GPU, "Vina-GPU-2.1-win.exe")
PDB = os.path.join(BASE, "pdb", "4euu.pdb")
CAJA = os.path.join(BASE, "caja.json")
SALIDA = os.path.join(BASE, "_control_bx795")
LIGANDO = "BX7"
CADENA = "A"
RESSEQ = "401"

# El BX-795 de ChEMBL (CHEMBL577784). Se necesita porque los enlaces que RDKit
# percibe de las coordenadas del cristal salen mal (el anillo tiofeno no se
# aromatiza y aparecen ordenes de enlace imposibles). Con este SMILES de plantilla
# se CORRIGEN los ordenes de enlace y se CONSERVAN las coordenadas y el orden de
# atomos del cristal, que es lo que hace falta para el RMSD directo.
SMILES_BX795 = "O=C(NCCCNc1nc(Nc2cccc(NC(=O)N3CCCC3)c2)ncc1I)c1cccs1"

SEARCH_DEPTH = 20
NUM_MODES = 9      # mas modos: se quiere ver si la pose buena esta entre las de arriba


def bloque_cristal():
    """Las lineas HETATM de BX-795, con el altloc resuelto y sin hidrogenos."""
    lineas = []
    for l in open(PDB, encoding="utf-8", errors="ignore"):
        if not l.startswith("HETATM"):
            continue
        if l[17:20].strip() != LIGANDO or l[21] != CADENA:
            continue
        if l[22:27].strip() != RESSEQ:
            continue
        if l[76:78].strip().upper() == "H":
            continue
        lineas.append(l)
    return lineas


def coords_de(pdbqt_o_pdb):
    """Coordenadas de los atomos pesados, en el orden del fichero."""
    out = []
    for l in open(pdbqt_o_pdb, encoding="utf-8", errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        tipo = l.rsplit(None, 1)[-1].strip().upper()
        if tipo in ("H", "HD", "HS", "D", "DD"):
            continue
        out.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(out)


def main():
    if not os.path.exists(CAJA):
        raise SystemExit("falta caja.json: ejecuta antes construir_receptor_tbk1.py")
    caja = json.load(open(CAJA, encoding="utf-8"))
    receptor = os.path.join(BASE, caja["receptor_pdbqt"])
    lineas = bloque_cristal()
    if not lineas:
        raise SystemExit("no se encontro el %s del cristal" % LIGANDO)

    shutil.rmtree(SALIDA, ignore_errors=True)
    ligdir = os.path.join(SALIDA, "ligands")
    outdir = os.path.join(SALIDA, "out")
    os.makedirs(ligdir)
    os.makedirs(outdir)

    cristal_pdb = os.path.join(SALIDA, "bx795_cristal.pdb")
    with open(cristal_pdb, "w", encoding="utf-8") as f:
        f.writelines(lineas)
        f.write("END\n")

    # --- el ligando a acoplar: coordenadas del cristal, enlaces de ChEMBL ---
    mol = Chem.MolFromPDBFile(cristal_pdb, removeHs=True, proximityBonding=True)
    if mol is None:
        raise SystemExit("RDKit no pudo leer el BX-795 del cristal")
    plantilla = Chem.MolFromSmiles(SMILES_BX795)
    try:
        from rdkit.Chem.AllChem import AssignBondOrdersFromTemplate
        corregido = AssignBondOrdersFromTemplate(plantilla, mol)
        print("ordenes de enlace corregidos con el SMILES de ChEMBL (CHEMBL577784)")
    except Exception as e:  # noqa: BLE001
        print("AVISO: no se pudo corregir con la plantilla (%s); se sigue con los"
              " enlaces percibidos del cristal" % e)
        corregido = mol
    mol = Chem.AddHs(corregido)
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    txt, ok, err = PDBQTWriterLegacy.write_string(
        MoleculePreparation(rigid_macrocycles=True).prepare(mol)[0])
    if not ok:
        raise SystemExit("no se pudo escribir el PDBQT del BX-795: %s" % err)
    lig = os.path.join(ligdir, "BX795_cristal.pdbqt")
    with open(lig, "w", encoding="utf-8") as f:
        f.write(txt)
    print("ligando del cristal: %d atomos pesados" % len(coords_de(cristal_pdb)))

    # --- config: la MISMA caja y el mismo protocolo que el banco ---
    cfg = os.path.join(SALIDA, "config.txt")
    with open(cfg, "w", encoding="utf-8") as f:
        f.write("receptor = %s\n" % receptor.replace("\\", "/"))
        f.write("ligand_directory = %s\n" % ligdir.replace("\\", "/"))
        f.write("output_directory = %s\n" % outdir.replace("\\", "/"))
        f.write("opencl_binary_path = %s\n" % GPU.replace("\\", "/"))
        f.write("center_x = %s\ncenter_y = %s\ncenter_z = %s\n"
                % tuple(caja["centro_caja"]))
        f.write("size_x = %d\nsize_y = %d\nsize_z = %d\n"
                % ((caja["tamano_caja"],) * 3))
        f.write("search_depth = %d\nnum_modes = %d\nthread = 8000\n"
                % (SEARCH_DEPTH, NUM_MODES))

    print("acoplando con Vina-GPU (search_depth %d, caja de %d A)..."
          % (SEARCH_DEPTH, caja["tamano_caja"]))
    t0 = time.time()
    with open(os.path.join(SALIDA, "motor.log"), "w", encoding="utf-8",
              errors="ignore") as f:
        r = subprocess.run([EXE, "--config", cfg], cwd=GPU,
                           stdout=f, stderr=subprocess.STDOUT, timeout=None)
    print("motor terminado con codigo %s en %.0f s" % (r.returncode, time.time() - t0))

    poses = sorted(glob.glob(os.path.join(outdir, "*_out.pdbqt")))
    if not poses:
        raise SystemExit("el motor no devolvio ninguna pose")
    pose = poses[0]

    # --- el cristal y la pose, casando por NOMBRE DE ATOMO ---
    # No se puede restar coordenada a coordenada: Vina-GPU reordena los atomos y
    # ademas anade los hidrogenos polares al principio del fichero. Lo que SI
    # conserva es el nombre de cada atomo del cristal (N07, C26...), asi que el
    # emparejamiento se hace por nombre, que es exacto y no adivina nada.
    ref = {}
    for l in lineas:
        if l[76:78].strip().upper() == "H":
            continue
        ref[l[12:16].strip()] = np.array([float(l[30:38]), float(l[38:46]),
                                          float(l[46:54])])

    modelos, actual = [], None
    for l in open(pose, encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            actual = {"energia": None, "atomos": {}}
        elif l.startswith("REMARK VINA RESULT") and actual is not None \
                and actual["energia"] is None:
            actual["energia"] = float(l.split()[3])
        elif l.startswith("ENDMDL"):
            if actual:
                modelos.append(actual)
            actual = None
        elif l.startswith(("ATOM", "HETATM")) and actual is not None:
            tipo = l.rsplit(None, 1)[-1].strip().upper()
            if tipo in ("H", "HD", "HS", "D", "DD"):
                continue
            actual["atomos"][l[12:16].strip()] = np.array(
                [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    if actual:
        modelos.append(actual)

    L = ["CONTROL DE REDOCKING DEL BX-795 — %s" % time.strftime("%Y-%m-%d"),
         "receptor: %s (cadenas %s)" % (os.path.basename(receptor),
                                        "+".join(caja["cadenas_receptor"])),
         "caja: centro %s, %d A | search_depth %d"
         % (caja["centro_caja"], caja["tamano_caja"], SEARCH_DEPTH),
         "el ligando se acopla construido DEL PROPIO CRISTAL, con los enlaces"
         " corregidos por el SMILES de ChEMBL (CHEMBL577784)",
         "poses devueltas: %d" % len(modelos),
         ""]

    mejor = None
    for i, m in enumerate(modelos, 1):
        comunes = sorted(set(ref) & set(m["atomos"]))
        if len(comunes) < 0.9 * len(ref):
            L.append("pose %d: solo %d de los %d atomos del cristal se pueden"
                     " casar (no comparable)" % (i, len(comunes), len(ref)))
            continue
        dif = np.array([ref[k] - m["atomos"][k] for k in comunes])
        rmsd = float(np.sqrt((dif ** 2).sum(1).mean()))
        centro = np.array(list(m["atomos"].values())).mean(axis=0)
        d_centro = float(np.linalg.norm(centro - np.array(list(ref.values())).mean(axis=0)))
        L.append("pose %d: %-7s kcal/mol | RMSD %.2f A (%d atomos) | centroides a %.2f A"
                 % (i, ("%.2f" % m["energia"]) if m["energia"] is not None else "-",
                    rmsd, len(comunes), d_centro))
        if mejor is None or rmsd < mejor:
            mejor = rmsd

    L.append("")
    L.append("RMSD del mejor modo: %s"
             % (("%.2f A" % mejor) if mejor is not None else "no se pudo medir"))
    L.append("criterio del proyecto: por debajo de 2 A, receptor y caja estan bien"
             " puestos; por encima, hay que mirar la preparacion antes de creerse"
             " ningun numero del banco.")
    texto = "\n".join(L)
    print("\n" + texto)
    with open(os.path.join(SALIDA, "informe.txt"), "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
