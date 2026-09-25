# -*- coding: utf-8 -*-
"""Audita los fallos de conversion/docking del cribado Kaggle v4 (TDP-43).

Extrae del CSV de resultados los ligandos que fallaron, cruza con la libreria
(_lib_consolidada.csv y full_library_solo.smi), clasifica el motivo y genera:
  1) ligandos_fallidos_kaggle_v4.csv  -> verdad del terreno (lo que fallo en Kaggle)
  2) libreria_problematicos_rdkit.csv  -> barrido RDKit de toda la libreria
     (metales, sales, multi-modelo/poses) con SMILES corregido o motivo de exclusion.
"""
import csv
import os
import re
import sys
from collections import Counter, defaultdict

from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover

BASE = r"C:\Users\Fredy\masive-als"
CSV_V4 = r"C:\Users\Fredy\kaggle_outputs\masive-als-cribado-master-tdp43-v4\resultados_master_tdp43v4.csv"
LIB_CONSOLIDADA = os.path.join(BASE, r"analysis\_lib_consolidada.csv")
LIB_SOLO = os.path.join(BASE, r"analysis\full_library_solo.smi")
OUT_FALLIDOS = os.path.join(BASE, r"analysis\ligandos_fallidos_kaggle_v4.csv")
OUT_PROBLEMATICOS = os.path.join(BASE, r"analysis\libreria_problematicos_rdkit.csv")

# Metales y metaloides que AutoDock no acepta (no tienen tipo AD)
METALES = set("Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po Fr Ra Ac Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Th Pa U Np Pu Am".split())
# Iones/cationes comunes de sales (para el stripper de RDKit y documentacion)
SALES_PATRON = re.compile(r"(SULFATE|SULPHATE|HYDROCHLORIDE|MALEATE|FUMARATE|CITRATE|TARTRATE|PHOSPHATE|ACETATE|SODIUM|POTASSIUM|CALCIUM|MAGNESIUM|ZINC|BROMIDE|CHLORIDE|NITRATE|SUCCINATE|MESYLATE|TOSYLATE|BESYLATE|EDISYLATE|BROMHYDRATE|HYDROBROMIDE|NAPRADISYLATE|OXALATE|LACTATE|MALATE|SALICYLATE|TRIHYDRATE|DIHYDRATE|MONOHYDRATE|HYDRATE)", re.I)


def read_csv_rows(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as handle:
        return list(csv.DictReader(handle))


def load_lib_smiles():
    """Libreria consolidada: ligand -> smiles"""
    mapping = {}
    if os.path.exists(LIB_CONSOLIDADA):
        for row in read_csv_rows(LIB_CONSOLIDADA):
            lig = (row.get("ligand") or "").strip()
            smi = (row.get("smiles") or "").strip()
            if lig:
                mapping[lig] = smi
    # full_library_solo.smi: "SMILES\tID" (puede haber mas IDs)
    if os.path.exists(LIB_SOLO):
        with open(LIB_SOLO, encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    smi, lig = parts[0].strip(), parts[1].strip()
                    if lig and lig not in mapping:
                        mapping[lig] = smi
    return mapping


def clasificar_razon(texto):
    t = texto
    if "multi-MODEL" in t:
        return "multi_modelo_pose"
    if "not a valid AutoDock type" in t:
        m = re.search(r'"([A-Z][a-z]?)" is not a valid', t)
        atom = m.group(1) if m else "?"
        return "metal_tipo_AD:%s" % atom
    if "Unknown or inappropriate tag" in t:
        return "etiqueta_desconocida"
    if "No valid ligands" in t:
        return "sin_ligandos_validos"
    if "affinity not found" in t:
        return "afinidad_no_encontrada"
    if "timeout" in t.lower() or "timed out" in t.lower():
        return "timeout"
    return "otro"


def razon_corta(texto):
    """Primera linea significativa del error para documentacion."""
    lines = [ln.strip() for ln in texto.splitlines() if ln.strip()]
    for ln in lines:
        if ln.startswith("Parse error"):
            return ln[:160]
        if ln.startswith("No valid"):
            return ln[:120]
        if ln.startswith("ATOM syntax"):
            return ln[:140]
    if lines:
        return lines[0][:160]
    return texto[:160]


def tiene_metal(mol):
    if mol is None:
        return None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in METALES:
            return atom.GetSymbol()
    return None


def analizar_smiles(smiles, nombre):
    """Devuelve (ok, mol, smiles_limpio, metal, motivo)."""
    if not smiles:
        return False, None, None, None, "sin_smiles"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, None, None, None, "smiles_invalido_rdkit"
    metal = tiene_metal(mol)
    # SaltRemover con definicion por defecto (quita iones comunes Na+ K+ Cl- etc.)
    remover = SaltRemover()
    try:
        limpio = remover.StripMol(mol, dontRemoveEverything=True)
    except Exception:
        limpio = None
    # Seleccion del fragmento organico mas grande sin metal (para sales multiples)
    mejor = None
    if limpio is not None and limpio.GetNumAtoms() > 0:
        frags = Chem.GetMolFrags(limpio, asMols=True, sanitizeFrags=False)
        organicos = []
        for fr in frags:
            if fr.GetNumAtoms() > 0 and not any(a.GetSymbol() in METALES for a in fr.GetAtoms()):
                organicos.append((fr.GetNumHeavyAtoms(), fr))
        if organicos:
            mejor = max(organicos, key=lambda x: x[0])[1]
        elif frags:
            mejor = max(frags, key=lambda x: x.GetNumHeavyAtoms())
    smiles_limpio = None
    if mejor is not None and mejor.GetNumAtoms() > 0:
        try:
            smiles_limpio = Chem.MolToSmiles(mejor)
        except Exception:
            smiles_limpio = None
    if metal:
        return False, mol, smiles_limpio, metal, "metal_%s" % metal
    if smiles_limpio and smiles_limpio != Chem.MolToSmiles(mol):
        return True, mol, smiles_limpio, None, "sal_strippable"
    return True, mol, smiles_limpio, None, "ok"


def main():
    if not os.path.exists(CSV_V4):
        print("NO_EXISTE_CSV", CSV_V4)
        sys.exit(1)

    rows = read_csv_rows(CSV_V4)
    print("filas_csv:", len(rows))

    # ---- 1) Verdad del terreno: ligandos que fallaron en Kaggle ----
    fallidos = {}          # ligand -> {razones:[], textos:[], n:int}
    ok_set = set()
    for row in rows:
        lig = (row.get("ligand") or "").strip()
        status = row.get("status") or ""
        if not lig:
            continue
        if status.strip() == "ok":
            ok_set.add(lig)
            continue
        if status.startswith("error:"):
            detalle = status[len("error:"):]
            cat = clasificar_razon(detalle)
            info = fallidos.setdefault(lig, {"razones": Counter(), "textos": [], "n": 0})
            info["razones"][cat] += 1
            info["n"] += 1
            corta = razon_corta(detalle)
            if corta not in info["textos"]:
                info["textos"].append(corta)

    # tambien puede haber ligandos que fallaron en las 3 semillas -> siempre error
    siempre_error = {lig for lig, info in fallidos.items() if lig not in ok_set}
    print("ok:", len(ok_set), "fallidos_unicos:", len(fallidos), "siempre_error:", len(siempre_error))

    lib = load_lib_smiles()
    print("smiles_en_libreria:", len(lib))

    with open(OUT_FALLIDOS, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ligand", "veces_error", "categoria", "smiles_libreria", "smiles_ok", "metal", "motivo", "detalle"])
        for lig in sorted(siempre_error):
            info = fallidos[lig]
            cat = info["razones"].most_common(1)[0][0]
            smi = lib.get(lig, "")
            ok, mol, limpio, metal, motivo = analizar_smiles(smi, lig)
            writer.writerow([
                lig, info["n"], cat, smi,
                limpio if limpio else "",
                metal if metal else "",
                motivo,
                " | ".join(info["textos"])[:400],
            ])
    print("escrito:", OUT_FALLIDOS)

    # ---- 2) Barrido RDKit de TODA la libreria (no solo lo que corrio en Kaggle) ----
    problematicos = []
    total = 0
    metal_counter = Counter()
    for lig, smi in sorted(lib.items()):
        total += 1
        ok, mol, limpio, metal, motivo = analizar_smiles(smi, lig)
        if not ok or metal:
            if metal:
                metal_counter[metal] += 1
            problematicos.append((lig, smi, metal or "", motivo, limpio or ""))
    # ademas: nombres con patron de pose/multi-modelo (archivos _out de acoplamientos previos)
    patron_pose = re.compile(r"_(out|SOD1|TDP43|TDP43_v2|FUS)(_out)?$|\.pdbqt_out$", re.I)
    poses = [(lig, lib[lig], "", "pose_acoplamiento_previo", "")
             for lig, smi in lib.items()
             if patron_pose.search(lig) and lig not in {p[0] for p in problematicos}]
    problematicos.extend(poses)

    with open(OUT_PROBLEMATICOS, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["ligand", "smiles", "metal", "motivo", "smiles_corregido"])
        for lig, smi, metal, motivo, limpio in problematicos:
            writer.writerow([lig, smi, metal, motivo, limpio])

    print("libreria_total:", total)
    print("problematicos_libreria:", len(problematicos), "| por metal:", dict(metal_counter.most_common(15)))
    print("escrito:", OUT_PROBLEMATICOS)

    # resumen consola
    cats = Counter()
    for info in fallidos.values():
        cats.update(info["razones"])
    print("categorias_fallos:", dict(cats.most_common(12)))


if __name__ == "__main__":
    main()