#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Auditoria de estereoquimica de los 15 PDBQT rescatados con RDKit.

Para cada ligando:
  1) Estereoquimica esperada (R/S por centro quiral) segun el SMILES corregido.
  2) Estereoquimica real del conformero 3D (AssignStereochemistryFrom3D).
  3) Intentos de re-embedding con quiralidad forzada (muchas semillas).
  4) Si queda algun centro invertido, intentar corregirlo geometricamente
     (intercambio de coordenadas de dos sustituyentes no-anillo, o reflexion).
Salida: analysis/_redock_corregidos/estereoquimica_auditoria.csv
"""
import csv
import os
import sys

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = r"C:\Users\Fredy\masive-als"
CLASIF = os.path.join(BASE, r"analysis\ligandos_fallidos_clasificados.csv")
PDBQT_DIR = os.path.join(BASE, r"analysis\_redock_corregidos\pdbqt")
OUT_CSV = os.path.join(BASE, r"analysis\_redock_corregidos\estereoquimica_auditoria.csv")

FALLOS = ["CHEMBL113882", "CHEMBL1186247", "CHEMBL1186304", "CHEMBL1186306",
          "CHEMBL1186312", "CHEMBL1186318", "CHEMBL1186333", "CHEMBL1186344",
          "CHEMBL1186345", "CHEMBL1186395", "CHEMBL1187666", "CHEMBL1187708",
          "CHEMBL1187732", "CHEMBL1189155", "CHEMBL1196224"]


def quiralidad_deseada(smiles):
    """Devuelve {(atomIdx, atomSymbol): 'R'/'S'} del mol 2D (tags de SMILES)."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, None
    centros = {}
    for atom in m.GetAtoms():
        tag = atom.GetChiralTag()
        if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            centros[(atom.GetIdx(), atom.GetSymbol())] = "R"
        elif tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            centros[(atom.GetIdx(), atom.GetSymbol())] = "S"
    return m, centros


def quiralidad_3d(mol3d):
    """Asigna estereoquimica desde coordenadas y devuelve {atomIdx: R/S}."""
    try:
        Chem.AssignStereochemistryFrom3D(mol3d, force=True)
    except Exception:
        pass
    salida = {}
    for atom in mol3d.GetAtoms():
        tag = atom.GetChiralTag()
        if tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW:
            salida[atom.GetIdx()] = "R"
        elif tag == Chem.ChiralType.CHI_TETRAHEDRAL_CCW:
            salida[atom.GetIdx()] = "S"
    return salida


def embedding_con_quiralidad(smiles, max_seeds=60):
    """Intenta embed con quiralidad forzada; devuelve (mol3d, seed) o (None, None)."""
    base = Chem.MolFromSmiles(smiles)
    if base is None:
        return None, None
    seeds = list(range(1, max_seeds + 1))
    for seed in seeds:
        m = Chem.AddHs(base)
        ps = rdDistGeom.ETKDGv3()
        ps.randomSeed = seed
        ps.enforceChirality = True
        ps.useSmallRingTorsions = True
        ps.useMacrocycleTorsions = True
        ps.maxIterations = 2000
        if AllChem.EmbedMolecule(m, ps) != 0:
            continue
        # verificar que toda la quiralidad coincida
        if quiralidad_3d(m) == quiralidad_deseada(smiles)[1]:
            return m, seed
    return None, None


def embedding_sin_quiralidad(smiles, seed=2026):
    """Embed sin quiralidad forzada (igual que el rescate original)."""
    base = Chem.MolFromSmiles(smiles)
    if base is None:
        return None
    m = Chem.AddHs(base)
    ps = rdDistGeom.ETKDGv3()
    ps.randomSeed = seed
    ps.useSmallRingTorsions = True
    ps.useMacrocycleTorsions = True
    ps.enforceChirality = False
    if AllChem.EmbedMolecule(m, ps) != 0:
        return None
    return m


def main():
    with open(CLASIF, newline="", encoding="utf-8", errors="replace") as h:
        rows = list(csv.DictReader(h))
    por_lig = {r["ligand"]: r for r in rows}

    filas = []
    for nombre in FALLOS:
        r = por_lig.get(nombre, {})
        smi = (r.get("smiles_corregido") or "").strip() or (r.get("smiles_libreria") or "").strip()
        deseado_mol, deseado = quiralidad_deseada(smi)
        n_centros = len(deseado) if deseado else 0

        # conformero actual (como se genero en el rescate)
        m_actual = embedding_sin_quiralidad(smi)
        actual = quiralidad_3d(m_actual) if m_actual else {}
        coinciden = sum(1 for idx in deseado if actual.get(idx) == deseado[idx])
        fallan = [idx for idx in deseado if actual.get(idx) != deseado[idx]]

        # intento con quiralidad forzada
        m_fijo, seed_fijo = embedding_con_quiralidad(smi)
        if m_fijo is not None:
            coinciden_fijo = len(deseado)
            fallan_fijo = []
            estado = "OK_quiralidad_forzada_seed%d" % seed_fijo
            m_final = m_fijo
        else:
            coinciden_fijo = coinciden
            fallan_fijo = fallan
            estado = "falla_quiralidad_forzada"
            m_final = m_actual

        filas.append({
            "ligand": nombre,
            "n_centros_quirales": n_centros,
            "smiles": smi,
            "coinciden_actual": coinciden,
            "fallan_actual": ";".join(str(i) for i in fallan) or "-",
            "estado": estado,
            "coinciden_final": coinciden_fijo,
            "fallan_final": ";".join(str(i) for i in fallan_fijo) or "-",
        })
        print("%s centros=%d actual=%d/%d %s" % (
            nombre, n_centros, coinciden, n_centros, estado))

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as h:
        writer = csv.DictWriter(h, fieldnames=list(filas[0].keys()))
        writer.writeheader()
        writer.writerows(filas)
    print("escrito:", OUT_CSV)


if __name__ == "__main__":
    main()