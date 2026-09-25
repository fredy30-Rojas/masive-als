#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Corrige la estereoquimica de los 15 PDBQT rescatados (version robusta).

Aprendizajes de la primera pasada:
  - MMFFOptimizeMolecule puede INVERTIR estereocentros (cruzando la barrera)
    en estos sistemas fusionados; por eso aqui NO se minimiza con MMFF despues
    del embedding. Las coordenadas ETKDG son un punto de partida valido para
    Vina, que luego genera sus propias conformaciones.
  - Reflejar un atomo de anillo a traves del centro quiral distorsiona el
    esqueleto y puede crear choques (0.35 A). Solo se reflejan atomos
    TERMINALES (H o atomo pesado con grado 1), y se verifica que no haya
    choques; si el primer candidato choca, se prueba el siguiente.

Protocolo:
  1) Embed sin quiralidad forzada (seed fija, ETKDGv3, anillos pequenos/macro).
  2) Asignar R/S desde 3D (firma de esta version RDKit: sin kwarg force).
  3) Para cada centro invertido: reflejar un vecino terminal a traves del centro.
  4) Verificar coincidencia total + distancia minima no-enlazada > 1.2 A.
  5) Escribir PDBQT corregido (o dejar el original si ya era correcto).

Salida:
  - _redock_corregidos/pdbqt/<ligand>.pdbqt  (sobrescritos si se corrigio)
  - _redock_corregidos/estereoquimica_auditoria.csv
  - _redock_corregidos/estereo_img/<ligand>.png
"""
import csv
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from rdkit.Geometry import Point3D

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = r"C:\Users\Fredy\masive-als"
CLASIF = os.path.join(BASE, r"analysis\ligandos_fallidos_clasificados.csv")
PDBQT_DIR = os.path.join(BASE, r"analysis\_redock_corregidos\pdbqt")
IMG_DIR = os.path.join(BASE, r"analysis\_redock_corregidos\estereo_img")
OUT_CSV = os.path.join(BASE, r"analysis\_redock_corregidos\estereoquimica_auditoria.csv")

FALLOS = ["CHEMBL113882", "CHEMBL1186247", "CHEMBL1186304", "CHEMBL1186306",
          "CHEMBL1186312", "CHEMBL1186318", "CHEMBL1186333", "CHEMBL1186344",
          "CHEMBL1186345", "CHEMBL1186395", "CHEMBL1187666", "CHEMBL1187708",
          "CHEMBL1187732", "CHEMBL1189155", "CHEMBL1196224"]

MIN_NOBOND = 1.2  # A; por debajo se considera choque


def centros_deseados(mol2d):
    return {a.GetIdx(): ("R" if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW else "S")
            for a in mol2d.GetAtoms()
            if a.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)}


def centros_3d(mol3d):
    try:
        Chem.AssignStereochemistryFrom3D(mol3d)
    except Exception as exc:
        print("AVISO AssignStereochemistryFrom3D:", exc)
    return {a.GetIdx(): ("R" if a.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW else "S")
            for a in mol3d.GetAtoms()
            if a.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)}


def min_dist_no_enlazada(mol):
    conf = mol.GetConformer()
    dmin = 99.0
    for i in range(mol.GetNumAtoms()):
        for j in range(i + 1, mol.GetNumAtoms()):
            if mol.GetBondBetweenAtoms(i, j):
                continue
            p1 = conf.GetAtomPosition(i)
            p2 = conf.GetAtomPosition(j)
            dd = (p1.x - p2.x) ** 2 + (p1.y - p2.y) ** 2 + (p1.z - p2.z) ** 2
            if dd < dmin * dmin:
                dmin = dd ** 0.5
    return dmin


def es_terminal(mol, idx):
    a = mol.GetAtomWithIdx(idx)
    return a.GetDegree() <= 1  # H o atomo pesado colgante


def invertir_centro(mol, conf, centro, vecino):
    pc = conf.GetAtomPosition(centro)
    pv = conf.GetAtomPosition(vecino)
    nuevo = Point3D(2.0 * pc.x - pv.x, 2.0 * pc.y - pv.y, 2.0 * pc.z - pv.z)
    conf.SetAtomPosition(vecino, nuevo)


def embed_sin_quiralidad(smiles, seed=2026):
    base = Chem.MolFromSmiles(smiles)
    if base is None:
        return None, None
    m = Chem.AddHs(base)
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = seed
    ps.useSmallRingTorsions = True
    ps.useMacrocycleTorsions = True
    ps.enforceChirality = False
    if AllChem.EmbedMolecule(m, ps) != 0:
        return None, None
    return m, base


def escribir_pdbqt(mol3d, salida):
    tmp = salida + ".sdf"
    w = Chem.SDWriter(tmp)
    w.write(mol3d)
    w.close()
    with open(tmp, "r", encoding="utf-8", errors="replace") as fh:
        sdf_text = fh.read()
    os.remove(tmp)
    from openbabel import openbabel as ob
    conv = ob.OBConversion()
    conv.SetInFormat("sdf")
    omol = ob.OBMol()
    if not conv.ReadString(omol, sdf_text) or not omol.Has3D():
        return False
    conv.SetOutFormat("pdbqt")
    out_s = conv.WriteString(omol)
    if not out_s or len(out_s) < 300:
        return False
    with open(salida, "w", encoding="utf-8") as fh:
        fh.write(out_s)
    return os.path.getsize(salida) > 300


def dibujar(mol2d, centros, nombre, salida):
    from rdkit.Chem.Draw import rdMolDraw2D
    mol = Chem.Mol(mol2d)
    for idx in centros:
        a = mol.GetAtomWithIdx(idx)
        a.SetProp("atomLabel", "%s%s" % (a.GetSymbol(), centros[idx]))
    d = rdMolDraw2D.MolDraw2DCairo(600, 500)
    opts = d.drawOptions()
    opts.addAtomIndices = False
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    with open(salida, "wb") as fh:
        fh.write(d.GetDrawingText())


def main():
    os.makedirs(IMG_DIR, exist_ok=True)
    with open(CLASIF, newline="", encoding="utf-8", errors="replace") as h:
        rows = list(csv.DictReader(h))
    por_lig = {r["ligand"]: r for r in rows}

    filas = []
    resumen = {"OK_sin_cambio": 0, "CORREGIDO": 0, "NO_CORREGIBLE": 0}
    for nombre in FALLOS:
        r = por_lig.get(nombre, {})
        smi = (r.get("smiles_corregido") or "").strip() or (r.get("smiles_libreria") or "").strip()
        m3d, m2d = embed_sin_quiralidad(smi)
        if m3d is None:
            filas.append({"ligand": nombre, "estado": "NO_EMBED", "centros": 0,
                          "coinciden": 0, "fallan": "-", "metodo": "-"})
            resumen["NO_CORREGIBLE"] += 1
            continue

        des = centros_deseados(m2d)
        conf = m3d.GetConformer()
        actual = centros_3d(m3d)
        fallan = [i for i in des if actual.get(i) != des[i]]
        metodo = "embed_sin_quiralidad"

        for intento in range(4):
            if not fallan:
                break
            cambiado = False
            for idx in list(fallan):
                vecinos = [a.GetIdx() for a in m3d.GetAtomWithIdx(idx).GetNeighbors()]
                candidatos = [v for v in vecinos if es_terminal(m3d, v)]
                for cand in candidatos:
                    invertir_centro(m3d, conf, idx, cand)
                    actual = centros_3d(m3d)
                    nuevos_fallan = [i for i in des if actual.get(i) != des[i]]
                    if idx not in nuevos_fallan:
                        cambiado = True
                        fallan = nuevos_fallan
                        metodo = "inversion_geometrica"
                        break
                    else:
                        # deshacer (reflejar de vuelta)
                        invertir_centro(m3d, conf, idx, cand)
            if not cambiado:
                break

        n_centros = len(des)
        n_ok = n_centros - len(fallan)
        dmin = min_dist_no_enlazada(m3d)
        if not fallan and dmin >= MIN_NOBOND:
            if metodo == "inversion_geometrica":
                estado = "CORREGIDO"
                resumen["CORREGIDO"] += 1
            else:
                estado = "OK_sin_cambio"
                resumen["OK_sin_cambio"] += 1
            ok_pdbqt = escribir_pdbqt(m3d, os.path.join(PDBQT_DIR, nombre + ".pdbqt"))
            if not ok_pdbqt:
                estado += "_PDBQT_FALLO"
        elif not fallan:
            estado = "CORREGIDO_SIN_MMFF_CHOQUE(%.2fA)" % dmin
            resumen["NO_CORREGIBLE"] += 1
            # aun asi escribimos el PDBQT (mejor que el invertido)
            escribir_pdbqt(m3d, os.path.join(PDBQT_DIR, nombre + ".pdbqt"))
        else:
            estado = "NO_CORREGIBLE"
            resumen["NO_CORREGIBLE"] += 1

        dibujar(m2d, des, nombre, os.path.join(IMG_DIR, nombre + ".png"))
        filas.append({
            "ligand": nombre, "estado": estado, "centros": n_centros,
            "coinciden": n_ok,
            "fallan": ";".join(str(i) for i in fallan) if fallan else "-",
            "metodo": metodo,
            "min_dist_A": round(dmin, 2),
        })
        print("%-18s %s centros=%d ok=%d/%d dmin=%.2f" % (
            nombre, estado, n_centros, n_ok, n_centros, dmin))

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as h:
        writer = csv.DictWriter(h, fieldnames=["ligand", "estado", "centros", "coinciden", "fallan", "metodo", "min_dist_A"])
        writer.writeheader()
        writer.writerows(filas)
    print("resumen:", resumen)
    print("escrito:", OUT_CSV)
    print("imagenes:", IMG_DIR)


if __name__ == "__main__":
    main()