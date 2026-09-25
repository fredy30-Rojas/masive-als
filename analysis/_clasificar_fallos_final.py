# -*- coding: utf-8 -*-
"""Clasificacion final de los ligandos que fallan en conversion/docking.

Salida: analysis/ligandos_fallidos_clasificados.csv con una accion por ligando:
  - CORREGIR_SMILES  -> sal strippable, se regenera el PDBQT desde el SMILES limpio
  - EXCLUIR_METAL    -> metal/metaloide organico, AutoDock no lo puede tipar
  - EXCLUIR_POSE     -> archivo de pose de un acoplamiento previo, no es ligando
  - REVISAR          -> etiqueta desconocida sin SMILES en libreria (revisar origen)
"""
import csv
import os
import re
from collections import Counter

from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover

BASE = r"C:\Users\Fredy\masive-als"
IN_FALLIDOS = os.path.join(BASE, r"analysis\ligandos_fallidos_kaggle_v4.csv")
LIB_CONSOLIDADA = os.path.join(BASE, r"analysis\_lib_consolidada.csv")
OUT = os.path.join(BASE, r"analysis\ligandos_fallidos_clasificados.csv")

METALES = set("Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po Fr Ra Ac Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Th Pa U Np Pu Am".split())

remover = SaltRemover()


def analizar_smiles(smiles, nombre):
    """Replica la logica del auditor: metal? sal strippable?"""
    if not smiles:
        return False, None, None, None, "sin_smiles"
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, None, None, None, "smiles_invalido_rdkit"
    metal = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() in METALES:
            metal = atom.GetSymbol()
            break
    try:
        limpio = remover.StripMol(mol, dontRemoveEverything=True)
    except Exception:
        limpio = None
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


def leer(archivo):
    with open(archivo, newline="", encoding="utf-8", errors="replace") as h:
        return list(csv.DictReader(h))

# Sufijos de sal/forma farmaceutica para resolver el nombre base del farmaco
SUFIJOS_SAL = re.compile(
    r"_(SULFATE|SULPHATE|HYDROCHLORIDE|MALEATE|FUMARATE|CITRATE|TARTRATE|PHOSPHATE|ACETATE|"
    r"SODIUM|POTASSIUM|CALCIUM|MAGNESIUM|ZINC|BROMIDE|CHLORIDE|NITRATE|SUCCINATE|MESYLATE|"
    r"TOSYLATE|BESYLATE|EDISYLATE|BROMHYDRATE|HYDROBROMIDE|NAPRADISYLATE|OXALATE|LACTATE|"
    r"MALATE|SALICYLATE|TRIHYDRATE|DIHYDRATE|MONOHYDRATE|HYDRATE|HEMIHYDRATE|HCL|NA|K|CA|MG|ZN)$",
    re.I,
)


def main():
    fallidos = leer(IN_FALLIDOS)
    lib_rows = leer(LIB_CONSOLIDADA)
    lib = {r["ligand"].strip(): r["smiles"].strip() for r in lib_rows if r.get("ligand")}
    # anadir full_library_solo.smi (nombres base como ABACAVIR sin sufijo)
    solo = os.path.join(BASE, r"analysis\full_library_solo.smi")
    if os.path.exists(solo):
        with open(solo, encoding="utf-8", errors="replace") as h:
            for line in h:
                partes = line.rstrip("\n").split("\t")
                if len(partes) >= 2:
                    smi, nom = partes[0].strip(), partes[1].strip()
                    if nom and nom not in lib:
                        lib[nom] = smi

    filas = []
    resumen = Counter()
    con_smiles = 0
    sin_smiles = 0
    for f in fallidos:
        lig = f["ligand"]
        categoria = f["categoria"]
        smi_lib = (f.get("smiles_libreria") or "").strip()
        smi_limpio = (f.get("smiles_ok") or "").strip()
        metal = (f.get("metal") or "").strip()
        motivo = f.get("motivo") or ""
        detalle = f.get("detalle") or ""

        # Fuente del SMILES: si no estaba, intentar resolver el nombre base (ABACAVIR_SULFATE -> ABACAVIR)
        smi_resuelto = ""
        if not smi_lib:
            m = SUFIJOS_SAL.search(lig)
            if m:
                base = lig[:m.start()]
                smi_resuelto = lib.get(base, "")
        if smi_lib:
            con_smiles += 1
        elif smi_resuelto:
            con_smiles += 1
        else:
            sin_smiles += 1

        # Si no habia SMILES y resolvimos el nombre base, usar el del base
        smi_efectivo = smi_lib or smi_resuelto
        if smi_resuelto and not smi_limpio:
            # recomputar SMILES limpio desde el base
            ok2, mol2, limpio2, metal2, motivo2 = analizar_smiles(smi_resuelto, lig)
            smi_limpio = limpio2 or ""
            if not metal:
                metal = metal2 or ""

        if categoria.startswith("metal_tipo_AD") or metal:
            accion = "EXCLUIR_METAL"
            razon = "atomos %s no tipables por AutoDock (organometalico/sal metalica)" % (metal or categoria.split(":")[-1])
            resumen["EXCLUIR_METAL"] += 1
        elif categoria == "multi_modelo_pose" or "multi-MODEL" in detalle:
            accion = "EXCLUIR_POSE"
            razon = "archivo de pose de acoplamiento previo (multi-MODEL), no es un ligando"
            resumen["EXCLUIR_POSE"] += 1
        elif smi_limpio and smi_limpio != smi_efectivo and smi_efectivo:
            accion = "CORREGIR_SMILES"
            razon = "sal/cuaternario: quitar contraion -> %s" % smi_limpio
            resumen["CORREGIR_SMILES"] += 1
        elif categoria == "etiqueta_desconocida" and not smi_efectivo:
            accion = "REVISAR"
            razon = "etiqueta desconocida en PDBQT y sin SMILES resoluble: revisar origen del archivo"
            resumen["REVISAR"] += 1
        elif categoria == "etiqueta_desconocida":
            accion = "CORREGIR_SMILES"
            razon = "PDBQT con etiqueta desconocida (sal); regenerar desde SMILES limpio: %s" % (smi_limpio or smi_efectivo)
            resumen["CORREGIR_SMILES"] += 1
        elif categoria == "afinidad_no_encontrada":
            accion = "REVISAR"
            razon = "Vina corrio pero no devolvio afinidad: revisar config/pose"
            resumen["REVISAR"] += 1
        else:
            accion = "REVISAR"
            razon = "otro motivo: %s" % (motivo or categoria)
            resumen["REVISAR"] += 1

        filas.append({
            "ligand": lig,
            "accion": accion,
            "categoria": categoria,
            "metal": metal,
            "smiles_libreria": smi_efectivo,
            "smiles_corregido": smi_limpio if smi_limpio else "",
            "razon": razon,
            "detalle": detalle[:200],
        })

    with open(OUT, "w", newline="", encoding="utf-8") as h:
        writer = csv.DictWriter(h, fieldnames=["ligand", "accion", "categoria", "metal", "smiles_libreria", "smiles_corregido", "razon", "detalle"])
        writer.writeheader()
        writer.writerows(filas)

    print("total_fallidos:", len(filas))
    print("con_smiles_libreria:", con_smiles, "| sin_smiles:", sin_smiles)
    print("acciones:", dict(resumen))
    print("escrito:", OUT)

    # Top ligandos CORREGIR con su smiles nuevo
    corregir = [f for f in filas if f["accion"] == "CORREGIR_SMILES"]
    print("\n=== ejemplo CORREGIR_SMILES (primeros 10) ===")
    for f in corregir[:10]:
        print("%s | %s | %s" % (f["ligand"], f["smiles_libreria"][:60], f["smiles_corregido"][:60]))


if __name__ == "__main__":
    main()