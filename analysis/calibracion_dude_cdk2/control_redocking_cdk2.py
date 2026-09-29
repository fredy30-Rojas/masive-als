#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CONTROL: se vuelve a acoplar el ligando del cristal de CDK2 (1h00).

POR QUE
-------
El banco CDK2 de DUD-E se apoya en las dos mismas cosas que el de TBK1: que el
receptor esta bien preparado y que la caja esta sobre el bolsillo de ATP. La
forma de comprobarlo es la del proyecto: **acoplar el ligando que ya venia
cocristalizado y ver a que distancia cae**. Liston: RMSD por debajo de 2 A.

EL ORDEN DE ATOMOS NO SE SUPONE: SE MIDE (29 sep 2026)
----------------------------------------------------
La primera version de este control daba por fallado el redocking con RMSD 6,6 A.
Estaba mal medido, y por esto: daba por hecho que "meeko respeta el orden de atomos
de la molecula y Vina-GPU mantiene ese orden en las poses", asi que restaba las
coordenadas de la pose contra las del cristal atomo con atomo, por indice.

No es verdad. Meeko escribe los atomos pesados en OTRO orden: el PDBQT del ligando
del cristal, leido en orden de fichero, empieza C,C,C,C,C,C,O,C,C,C,N,... y el mol2
empieza C,C,C,C,C,C,C,C,C,... Restar por indice mete cada coordenada en el atomo
equivocado, y lo que se mide no es la pose: es una molecula con la geometria
desordenada. Contra eso no salva ni el GetBestRMS, porque GetBestRMS casa por
simetria del grafo, no arregla una permutacion arbitraria de las COORDENADAS.

El arreglo es no suponer: la correspondencia se DEDUCE. El PDBQT que se acopla es
el ligando del cristal tal cual (o deberia serlo), asi que se busca la asignacion
optima entre sus atomos y los del mol2 exigiendo que cada atomo caiga en uno del
mismo elemento. Si de verdad es la pose del cristal, esa asignacion da ~0,00 A y
devuelve la permutacion exacta. Con esa permutacion, la resta por indice vuelve a
significar lo que dice y el GetBestRMS se aplica sobre una molecula bien montada.

Si la asignacion NO da ~0,00 A, el control se para: significa que lo que se acoplo
no era la pose del cristal y el problema es anterior al motor.

QUE MIDE, Y POR QUE DOS RMSD
----------------------------
1. `RMSD indice`: resta heavy-atom a heavy-atom, ya con la permutacion medida.
2. `RMSD mol` (GetBestRMS de RDKit): el mismo calculo casando por estructura y
   simetria, sobre un mol RENUMERADO al orden del fichero. Si los dos coinciden,
   el numero es de fiar.
3. `piso`: el RMSD MAS PEQUENO que puede dar cualquier emparejamiento de atomos
   que respete el elemento (asignacion optima libre). No es la medida del
   proyecto, es un SUELO: si hasta el mejor emparejamiento posible esta por
   encima del liston, el veredicto no se puede discutir.
4. `solape`: distancia de cada atomo del cristal al MAS CERCANO de la pose. Si
   el RMSD sale alto y el solape bajo, la pose esta en la zona pero con los
   atomos cambiados de sitio (ligando dado la vuelta, caso clasico de un
   bolsillo plano como el de ATP).

OJO CON EL GetBestRMS: SUPERPONE ANTES DE MEDIR (29 sep 2026)
------------------------------------------------------------
`GetBestRMS` de RDKit alinea la pose sobre el cristal y DESPUES mide. Eso responde
"se conserva la forma?", no responde "esta en el mismo sitio del bolsillo?", que es
lo que se pregunta un redocking. En CDK2 la diferencia es enorme y hubo que verla
para no equivocarse: la mejor pose esta a 6,87 A en el bolsillo y a 1,55 A tras
alinearla. La forma se conserva; lo que no se conserva es el sitio.

Asi que la medida del proyecto (la que decide con el liston de 2 A) es la del
BOLSILLO, sin alinear. El `RMSD alineado` se deja al lado, pero etiquetado como lo
que es: un diagnostico de forma, no el criterio.

EL LIGANDO
----------
Viene en `crystal_ligand.mol2`, con hidrogenos explicitos (24) y un amonio
cuaternario. Se quitan los hidrogenos ANTES de anadirlos de nuevo con meeko: si
se dejan los viejos, el PDBQT sale con hidroguenas polares congelados en la
pose del cristal y el RMSD por GetBestRMS (que compara moleculas completas)
mide una resta de atomos que no se movieron.

BARRIDO DE PROTOCOLO
--------------------
Con la caja de 24 A de la regla del proyecto el control falla. Antes de culpar
al motor se prueba lo que se prueba siempre: caja ajustada al bolsillo y mas
exhaustividad. Cada variante escribe su carpeta y una linea en
`barrido_redocking_cdk2.txt`.

Uso:
    python control_redocking_cdk2.py --barrido            # 20/22/24 A x depth 20/32
    python control_redocking_cdk2.py --caja 20 --depth 32
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import time
from collections import OrderedDict

import numpy as np
from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem, RDLogger
from scipy.optimize import linear_sum_assignment

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(BASE))
GPU = os.path.join(RAIZ, "gpu_dock")
EXE = os.path.join(GPU, "Vina-GPU-2.1-win.exe")
MOL2 = os.path.join(BASE, "crystal_ligand.mol2")
PDB = os.path.join(BASE, "receptor_cdk2.pdb")
CAJA = os.path.join(BASE, "caja_cdk2.json")
RESUMEN = os.path.join(BASE, "barrido_redocking_cdk2.txt")

NUM_MODES = 9
LISTON_A = 2.0
H_TIPOS = ("H", "HD", "HS", "D", "DD")
# El tipo de AutoDock dice el elemento y algo mas (A = carbono aromatico, NA =
# nitrogeno aceptor, OA = oxigeno aceptor). Para emparejar por elemento hay que
# quedarse con el elemento, que es lo unico que el PDBQT y el mol2 comparten.
GRUPO = {"A": "C", "C": "C", "N": "N", "NA": "N", "OA": "O", "O": "O",
         "SA": "S", "S": "S"}
# Lo que tiene que dar la asignacion para que se crea la correspondencia.
TOL_ENTRADA = 0.05


def log(m):
    print(m, flush=True)


def leer_ligando():
    """(mol sin H con coordenadas del cristal, coords de atomos pesados)."""
    mol = Chem.MolFromMol2File(MOL2, removeHs=True, sanitize=False)
    if mol is None:
        raise SystemExit("no se pudo leer %s" % MOL2)
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                     ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    sin_h = Chem.RemoveHs(mol, sanitize=False)
    Chem.SanitizeMol(sin_h)
    conf = sin_h.GetConformer()
    coords = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                        conf.GetAtomPosition(i).z]
                       for i in range(sin_h.GetNumAtoms())])
    return sin_h, coords


def atomos_pesados_pdbqt(ruta):
    """(tipos, coords) de los atomos pesados de un PDBQT, EN ORDEN DE FICHERO."""
    tipos, coords = [], []
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        tipo = l.rsplit(None, 1)[-1].strip().upper()
        if tipo in H_TIPOS:
            continue
        tipos.append(tipo)
        coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return tipos, np.array(coords)


def emparejar_orden(tipos, coords, mol_ref, ref):
    """Que atomo del mol2 ocupa cada posicion del fichero, y con cuanto error.

    Devuelve (permutacion, rmsd). Es una asignacion optima (Hungarian) con la
    distancia como coste y con los elementos distintos prohibidos, que es la
    unica forma de resolver el orden sin suponerlo. Con un ligando que se
    simetrico puede haber varias permutaciones validas: da igual, cualquiera de
    ellas sirve para medir, porque todas llevan los atomos al atomo correcto (o a
    su simetrico).
    """
    g_fich = [GRUPO.get(t.upper(), t.upper()) for t in tipos]
    g_ref = [GRUPO.get(a.GetSymbol().upper(), a.GetSymbol().upper())
             for a in mol_ref.GetAtoms()]
    if sorted(g_fich) != sorted(g_ref):
        raise SystemExit("el PDBQT y el mol2 no tienen los mismos elementos:"
                         " %d pesados contra %d" % (len(g_fich), len(g_ref)))
    coste = np.linalg.norm(np.asarray(coords)[:, None, :] - ref[None, :, :],
                           axis=2)
    for i, gi in enumerate(g_fich):
        for j, gj in enumerate(g_ref):
            if gi != gj:
                coste[i, j] = 1e6
    _, col = linear_sum_assignment(coste)
    rmsd = float(np.sqrt(((np.asarray(coords) - ref[col]) ** 2).sum(1).mean()))
    return col, rmsd


def piso_rmsd(pos, ref, tipos, mol_ref):
    """El RMSD menor posible entre cualquier emparejamiento que respete el elemento.

    Asignacion optima (Hungarian) con la distancia como coste y los elementos
    distintos prohibidos. Es un SUELO del RMSD, no el RMSD: sirve para decir "ni
    con el mejor emparejamiento imaginable baja del liston", que es la forma mas
    fuerte de decir que una pose no esta donde el cristal.
    """
    g_fich = [GRUPO.get(t.upper(), t.upper()) for t in tipos]
    g_ref = [GRUPO.get(a.GetSymbol().upper(), a.GetSymbol().upper())
             for a in mol_ref.GetAtoms()]
    coste = np.linalg.norm(np.asarray(pos)[:, None, :] - ref[None, :, :], axis=2)
    for i, gi in enumerate(g_fich):
        for j, gj in enumerate(g_ref):
            if gi != gj:
                coste[i, j] = 1e6
    _, col = linear_sum_assignment(coste)
    return float(np.sqrt(((np.asarray(pos) - ref[col]) ** 2).sum(1).mean()))


def escribir_ligando(mol, carpeta):
    """PDBQT del ligando del cristal, atomos pesados y H polares de meeko."""
    os.makedirs(carpeta, exist_ok=True)
    con_h = Chem.AddHs(Chem.Mol(mol), addCoords=True)
    txt, ok, err = PDBQTWriterLegacy.write_string(
        MoleculePreparation(rigid_macrocycles=True).prepare(con_h)[0])
    if not ok:
        raise SystemExit("meeko no pudo escribir el PDBQT: %s" % err)
    ruta = os.path.join(carpeta, "FAP_cristal.pdbqt")
    with open(ruta, "w", encoding="utf-8") as f:
        f.write(txt)
    return ruta


def acoplar(carpeta, receptor, centro, tamano, depth):
    ligdir = os.path.join(carpeta, "ligands")
    outdir = os.path.join(carpeta, "out")
    os.makedirs(ligdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)
    cfg = os.path.join(carpeta, "config.txt")
    with open(cfg, "w", encoding="utf-8") as f:
        f.write("receptor = %s\n" % receptor.replace("\\", "/"))
        f.write("ligand_directory = %s\n" % ligdir.replace("\\", "/"))
        f.write("output_directory = %s\n" % outdir.replace("\\", "/"))
        f.write("opencl_binary_path = %s\n" % GPU.replace("\\", "/"))
        f.write("center_x = %.3f\ncenter_y = %.3f\ncenter_z = %.3f\n" % tuple(centro))
        f.write("size_x = %d\nsize_y = %d\nsize_z = %d\n" % (tamano, tamano, tamano))
        f.write("search_depth = %d\nnum_modes = %d\nthread = 8000\n"
                % (depth, NUM_MODES))
    log("acoplando: caja %d A, centro %s, search_depth %d"
        % (tamano, np.round(centro, 2), depth))
    t0 = time.time()
    with open(os.path.join(carpeta, "motor.log"), "w", encoding="utf-8",
              errors="ignore") as f:
        r = subprocess.run([EXE, "--config", cfg], cwd=GPU,
                           stdout=f, stderr=subprocess.STDOUT, timeout=None)
    log("motor: codigo %s en %.0f s" % (r.returncode, time.time() - t0))
    import glob
    poses = sorted(glob.glob(os.path.join(outdir, "*_out.pdbqt")))
    if not poses:
        return [], None
    modelos, actual = [], None
    for l in open(poses[0], encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            actual = {"energia": None, "pos": []}
        elif l.startswith("REMARK VINA RESULT") and actual is not None \
                and actual["energia"] is None:
            actual["energia"] = float(l.split()[3])
        elif l.startswith("ENDMDL"):
            if actual:
                modelos.append(actual)
            actual = None
        elif l.startswith(("ATOM", "HETATM")) and actual is not None:
            tipo = l.rsplit(None, 1)[-1].strip().upper()
            if tipo in H_TIPOS:
                continue
            actual["pos"].append([float(l[30:38]), float(l[38:46]),
                                  float(l[46:54])])
    if actual:
        modelos.append(actual)
    return modelos, poses[0]


def medir(mol_ref, ref, modelos, perm, tipos_fichero=None):
    """RMSD del bolsillo, piso, RMSD alineado y solape, pose por pose.

    `perm` es la correspondencia fichero -> mol2, medida en el PDBQT de entrada,
    no supuesta. `ref[perm]` pone las coordenadas del cristal en el orden del
    fichero, que es cuando la resta por indice significa algo.

    La medida que decide es la del BOLSILLO (sin alinear). El `rmsd_alineado` se
    calcula tambien, pero es otra cosa: responde si la forma coincide, no si la
    pose esta en el mismo sitio.
    """
    ref_en_orden = ref[perm]
    # El mol de referencia, renumerado al orden del fichero: asi las coordenadas
    # de la pose y el grafo vuelven a ser la MISMA molecula y el RMSD alineado
    # mide lo que tiene que medir.
    mol_en_orden = Chem.RenumberAtoms(mol_ref, [int(x) for x in perm])
    filas = []
    for i, m in enumerate(modelos, 1):
        pos = np.array(m["pos"])
        if len(pos) != len(ref):
            filas.append((i, m["energia"], None, None, None, None, len(pos)))
            continue
        rmsd_bolsillo = float(np.sqrt(((ref_en_orden - pos) ** 2).sum(1).mean()))
        piso = (piso_rmsd(pos, ref, tipos_fichero, mol_ref)
                if tipos_fichero else None)
        pose_mol = Chem.Mol(mol_en_orden)
        cf = pose_mol.GetConformer()
        for k, p in enumerate(pos):
            cf.SetAtomPosition(k, (float(p[0]), float(p[1]), float(p[2])))
        try:
            rmsd_alineado = float(Chem.rdMolAlign.GetBestRMS(pose_mol, mol_ref))
        except Exception:  # noqa: BLE001
            rmsd_alineado = float("nan")
        d = np.linalg.norm(ref[:, None, :] - pos[None, :, :], axis=2)
        solape = float(np.sqrt(d.min(axis=1).mean() ** 2))
        d_centro = float(np.linalg.norm(pos.mean(axis=0) - ref.mean(axis=0)))
        filas.append((i, m["energia"], rmsd_bolsillo, piso, rmsd_alineado,
                      solape, d_centro))
    return filas


def variante(receptor, mol, ref, centro, tamano, depth, etiqueta):
    carpeta = os.path.join(BASE, "_control_cdk2", etiqueta)
    shutil.rmtree(carpeta, ignore_errors=True)
    escribir_ligando(mol, os.path.join(carpeta, "ligands"))
    # Antes de medir nada: la correspondencia. Se saca del PDBQT de entrada, que
    # tiene que ser la pose del cristal; si no lo es, aqui se para.
    tipos, coords = atomos_pesados_pdbqt(os.path.join(carpeta, "ligands",
                                                      "FAP_cristal.pdbqt"))
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
                    help="tamano de caja; 0 = el de caja_cdk2.json")
    ap.add_argument("--depth", type=int, default=20)
    ap.add_argument("--barrido", action="store_true",
                    help="prueba cajas 20/22/24 con depth 20 y 32")
    ap.add_argument("--centro", default="bolsillo",
                    help="'bolsillo' (centro de residuos <8 A) o 'ligando'")
    args = ap.parse_args()

    if not os.path.exists(CAJA):
        raise SystemExit("falta caja_cdk2.json: ejecuta antes preparar_receptor_cdk2.py")
    caja = json.load(open(CAJA, encoding="utf-8"))
    receptor = os.path.join(BASE, caja["receptor_pdbqt"])
    mol, ref = leer_ligando()
    log("ligando del cristal: %d atomos pesados, enlaces aromaticos %d"
        % (mol.GetNumAtoms(), sum(1 for b in mol.GetBonds() if b.GetIsAromatic())))

    ligando = np.array(caja["centro_caja"], dtype=float)
    if args.centro == "bolsillo":
        # centro de los atomos de proteina a menos de 8 A del ligando: la regla
        # del proyecto (construir_receptor_tbk1.py)
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

    os.makedirs(os.path.dirname(RESUMEN), exist_ok=True)
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
