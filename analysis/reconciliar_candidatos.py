# -*- coding: utf-8 -*-
"""
reconciliar_candidatos.py
MASIVE-ALS — Reconciliacion de candidatos (18 ago 2026, Buffy)

Genera el CSV DE PRODUCCION unico de candidatos a partir de la tanda z001
COMPLETA (14.952 pares, 4.984 x 3), aplicando la metodologia documentada:

  1) Corte top-5% POR PROTEINA (mas negativo = mejor):
       SOD1  <= -7.20
       TDP43 <= -6.40
       FUS   <= -6.40
  2) Union con SMILES (libreria consolidada + extra ChEMBL corregida)
  3) Dedup por SMILES canonico (por proteina)
  4) Filtro PAINS + Lipinski/Veber
  5) Filtro CNS: TPSA <= 90 A^2 y MW <= 450 (paradigma CNS-MPO, Wager 2010)
  6) Columnas: aprobado (FDA si/no), cns (bool), cns_mpo (score 4-comp,
     APROXIMACION documentada en el paper: no sustituye al CNS-MPO de 6
     componentes de Wager que incluye pKa y logD), mw, tpsa, logp, hbd,
     hba, rotb, pose_pdbqt

NOTA (reconciliacion con v4 de Claude): el CSV _candidatos_v4.csv de Claude
usa umbrales de afinidad mas estrictos (SOD1 -7.5 / TDP43 -6.8 / FUS -6.7) y
NO aplica filtro CNS (solo anade la columna de score). Este script produce el
CSV de produccion con la metodologia de la revision cientifica (top-5% por
proteina + filtro CNS). Ambos son validos; el de produccion es este.
"""
import os
import sys
import csv
import glob
import collections

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, FilterCatalog
from rdkit.Chem import Crippen, Lipinski

RDLogger.DisableLog("rdApp.*")

BASE = r"C:\Users\Fredy\masive-als"
ANALYSIS = os.path.join(BASE, "analysis")
Z001 = os.path.join(BASE, "gpu_dock", "tanda_z001", "resultados_z001.csv")

UMBRALES = {"SOD1": -7.20, "TDP43": -6.40, "FUS": -6.40}

# librerias de SMILES
LIBRERIAS = [
    os.path.join(ANALYSIS, "_lib_consolidada.csv"),
    os.path.join(ANALYSIS, "libreria_extra_chembl.csv"),
    os.path.join(BASE, "compounds", "full_fda_library.csv"),
    os.path.join(BASE, "compounds", "batch2_chembl.csv"),
]
FDA = os.path.join(BASE, "compounds", "fda_subset.csv")

PAINS = FilterCatalog.FilterCatalogParams()
PAINS.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
PAINS_CAT = FilterCatalog.FilterCatalog(PAINS)


def cargar_smiles():
    """Devuelve {id_compuesto: smiles} desde todas las librerias."""
    mapa = {}
    for path in LIBRERIAS:
        if not os.path.exists(path):
            continue
        df = pd.read_csv(path)
        cols = list(df.columns)
        idcol = None
        for c in ["chembl_id", "name", "id", "identifier", "ligand"]:
            if c in cols:
                idcol = c
                break
        smicol = None
        for c in ["smiles", "SMILES", "canonical_smiles"]:
            if c in cols:
                smicol = c
                break
        if idcol is None or smicol is None:
            continue
        for _, row in df.iterrows():
            i = str(row[idcol]).strip()
            s = str(row[smicol]).strip()
            if i and s and s.lower() not in ("nan", "none"):
                if i not in mapa:
                    mapa[i] = s
    return mapa


def cargar_fda():
    aprobados = set()
    if os.path.exists(FDA):
        df = pd.read_csv(FDA)
        cols = list(df.columns)
        idcol = None
        for c in ["chembl_id", "name", "id"]:
            if c in cols:
                idcol = c
                break
        if idcol:
            for v in df[idcol]:
                aprobados.add(str(v).strip())
    return aprobados


def canon(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    return Chem.MolToSmiles(m, canonical=True)


def lipinski_ok(mol):
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    fails = 0
    if mw > 500:
        fails += 1
    if logp > 5:
        fails += 1
    if hbd > 5:
        fails += 1
    if hba > 10:
        fails += 1
    return fails <= 1, mw, logp, hbd, hba


def cns_mpo_aprox(mol, mw, tpsa, logp, hbd):
    """Score CNS-MPO aproximado de 4 componentes (0-6), NO el Wager de 6.
    Documentado como aproximacion: omite pKa y logD a pH fisiologico."""
    score = 0.0
    # MW <= 450 -> 1.0 (lineal 500->0)
    if mw <= 450:
        score += 1.0
    elif mw >= 500:
        score += 0.0
    else:
        score += (500 - mw) / 50.0
    # TPSA <= 90 -> 1.0 (lineal 120->0)
    if tpsa <= 90:
        score += 1.0
    elif tpsa >= 120:
        score += 0.0
    else:
        score += (120 - tpsa) / 30.0
    # logP 2-4 -> 1.0 (fuera baja)
    if 2 <= logp <= 4:
        score += 1.0
    else:
        score += max(0.0, 1.0 - abs(logp - 3.0) / 3.0)
    # HBD <= 1 -> 1.0 (lineal 3->0)
    if hbd <= 1:
        score += 1.0
    elif hbd >= 3:
        score += 0.0
    else:
        score += (3 - hbd) / 2.0
    return round(score, 2)


def localizar_pose(target, lig):
    for tanda in ["z001", "z002", "z003"]:
        d = os.path.join(BASE, "gpu_dock", "tanda_%s" % tanda, "results_%s" % target)
        p = os.path.join(d, "%s_out.pdbqt" % lig)
        if os.path.exists(p):
            return p
    return ""


def main():
    print("Cargando z001 completa...")
    df = pd.read_csv(Z001)
    print("  filas:", len(df))
    for t in UMBRALES:
        n = (df["target"] == t).sum()
        print("  %s: %d" % (t, n))

    print("Cargando librerias SMILES...")
    mapa = cargar_smiles()
    print("  ids con SMILES:", len(mapa))
    aprobados = cargar_fda()
    print("  farmacos FDA:", len(aprobados))

    print("Embudo top-5% por proteina + PAINS + Lipinski + CNS...")
    filas = []
    saltados = collections.Counter()
    for target, umbral in UMBRALES.items():
        sub = df[df["target"] == target]
        sub = sub[sub["affinity"] <= umbral]
        print("  %s: %d pares <= %.2f" % (target, len(sub), umbral))
        vistos = set()
        for _, row in sub.iterrows():
            lig = str(row["ligand"])
            aff = float(row["affinity"])
            smiles = mapa.get(lig)
            if not smiles:
                saltados["sin_smiles"] += 1
                continue
            can = canon(smiles)
            if can is None:
                saltados["smiles_invalido"] += 1
                continue
            if can in vistos:
                saltados["dedup"] += 1
                continue
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                saltados["mol_invalido"] += 1
                continue
            if PAINS_CAT.HasMatch(mol):
                saltados["pains"] += 1
                continue
            ok, mw, logp, hbd, hba = lipinski_ok(mol)
            if not ok:
                saltados["lipinski"] += 1
                continue
            tpsa = Descriptors.TPSA(mol)
            rotb = Lipinski.NumRotatableBonds(mol)
            # Filtro CNS: TPSA<=90 y MW<=450
            if tpsa > 90 or mw > 450:
                saltados["cns_no"] += 1
                continue
            vistos.add(can)
            filas.append({
                "ligand": lig,
                "target": target,
                "affinity": round(aff, 2),
                "smiles": smiles,
                "aprobado": "si" if lig in aprobados else "no",
                "cns": True,
                "cns_mpo": cns_mpo_aprox(mol, mw, tpsa, logp, hbd),
                "cns_mpo_ok": cns_mpo_aprox(mol, mw, tpsa, logp, hbd) >= 4.0,
                "mw": round(mw, 2),
                "tpsa": round(tpsa, 2),
                "logp": round(logp, 2),
                "hbd": hbd,
                "hba": hba,
                "rotb": rotb,
                "pose_pdbqt": localizar_pose(target, lig),
            })

    print("Saltados:", dict(saltados))
    filas.sort(key=lambda r: (r["target"], r["affinity"]))
    salida = os.path.join(ANALYSIS, "candidatos_filtrados_v4.csv")
    with open(salida, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        w.writeheader()
        w.writerows(filas)
    print("Escrito:", salida, "->", len(filas), "candidatos")
    cnt = collections.Counter(r["target"] for r in filas)
    print("Por proteina:", dict(cnt))
    n_aprob = sum(1 for r in filas if r["aprobado"] == "si")
    n_mpo4 = sum(1 for r in filas if r["cns_mpo_ok"])
    print("Aprobados FDA:", n_aprob, "| CNS-MPO>=4:", n_mpo4)


if __name__ == "__main__":
    main()
