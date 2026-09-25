# -*- coding: utf-8 -*-
"""Exporta los SMILES corregidos (accion CORREGIR_SMILES) como .smi listo para regenerar PDBQT."""
import csv
import os

BASE = r"C:\Users\Fredy\masive-als"
IN = os.path.join(BASE, r"analysis\ligandos_fallidos_clasificados.csv")
OUT = os.path.join(BASE, r"analysis\ligandos_a_reacoplar.smi")

with open(IN, newline="", encoding="utf-8", errors="replace") as h:
    rows = list(csv.DictReader(h))

corregir = [r for r in rows if r["accion"] == "CORREGIR_SMILES"]
sin_smiles = [r for r in corregir if not (r.get("smiles_corregido") or "").strip()]
con_smiles = [r for r in corregir if (r.get("smiles_corregido") or "").strip()]

with open(OUT, "w", encoding="utf-8") as h:
    for r in con_smiles:
        h.write("%s\t%s\n" % (r["smiles_corregido"].strip(), r["ligand"].strip()))

print("CORREGIR_SMILES:", len(corregir))
print("  con SMILES corregido exportados:", len(con_smiles))
print("  sin SMILES (requieren fetch):", len(sin_smiles))
print("escrito:", OUT)
# sanity: SMILES validos con RDKit
from rdkit import Chem
invalidos = [r["ligand"] for r in con_smiles if Chem.MolFromSmiles(r["smiles_corregido"].strip()) is None]
print("SMILES corregidos invalidos para RDKit:", len(invalidos), invalidos[:10])