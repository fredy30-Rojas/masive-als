#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""buscar_activos_chembl.py — Recolecta compuestos activos conocidos contra SOD1.

Consulta ChEMBL (target CHEMBL2354 = Superoxide dismutase [Cu-Zn] humano),
filtra los compuestos con actividad documentada (pChEMBL >= 5, IC50 <= 10 uM),
obtiene su SMILES canónico y escribe un CSV de activos para validar_senuelos.py.
Además añade los controles positivos de la literatura: LCS-1 y PRG-A01 (PubChem).

Salida: activos_sod1.csv (name, smiles, target)
"""
import csv

import requests

TARGET = "CHEMBL2354"  # SOD1 [Cu-Zn] humano
BASE = "https://www.ebi.ac.uk/chembl/api/data"

# Controles positivos de la literatura (PubChem)
LITERATURA = [
    ("LCS-1", "CC1=CC(=CC=C1)N2C(=O)C(=C(C=N2)Cl)Cl", "SOD1"),
    ("PRG-A01", "CC1(C(CC2=C(O1)C=C3C(=C2)C=CC(=O)O3)OC(=O)C=CC4=CC(=C(C=C4)O)OC)C", "SOD1"),
]


def main():
    compuestos = {}  # mol_id -> pchembl

    for std_type in ("IC50", "Ki", "EC50", "Kd"):
        url = (f"{BASE}/activity.json?target_chembl_id={TARGET}"
               f"&standard_type={std_type}&limit=500")
        try:
            resp = requests.get(url, timeout=60).json()
        except Exception as e:
            print(f"  {std_type}: error {e}")
            continue
        acts = resp.get("activities", [])
        print(f"{std_type}: {len(acts)} actividades")
        for a in acts:
            mol_id = a.get("molecule_chembl_id")
            pch = a.get("pchembl_value")
            if not mol_id:
                continue
            try:
                p = float(pch)
            except (TypeError, ValueError):
                continue
            if p >= 5.0:  # IC50/Ki/EC50/Kd <= 10 uM
                compuestos[mol_id] = max(compuestos.get(mol_id, 0.0), p)

    print(f"Compuestos activos (pChEMBL >= 5): {len(compuestos)}")

    filas = list(LITERATURA)
    for mol_id, p in sorted(compuestos.items(), key=lambda x: -x[1]):
        try:
            m = requests.get(f"{BASE}/molecule/{mol_id}.json", timeout=60).json()
            smi = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if smi:
                filas.append((mol_id, smi, "SOD1"))
            else:
                print(f"  {mol_id}: sin SMILES")
        except Exception as e:
            print(f"  {mol_id}: error {e}")

    with open("activos_sod1.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["name", "smiles", "target"])
        for r in filas:
            w.writerow(r)
    print(f"CSV escrito: activos_sod1.csv con {len(filas)} activos "
          f"({len(filas)-len(LITERATURA)} de ChEMBL + {len(LITERATURA)} de literatura)")


if __name__ == "__main__":
    main()
