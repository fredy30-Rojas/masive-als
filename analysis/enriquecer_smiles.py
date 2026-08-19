# -*- coding: utf-8 -*-
"""
enriquecer_smiles.py
MASIVE-ALS — Completa los SMILES que faltan en las librerias locales.

Lee un CSV de resultados de docking, encuentra los ligandos que NO tienen
SMILES en las librerias locales (compounds/), los busca en la API publica
de ChEMBL y guarda el resultado en compounds/libreria_extra_chembl.csv,
que filtro_rescoring.py detecta automaticamente en su proxima corrida.

Uso:
  python enriquecer_smiles.py ../gpu_dock/resultados_vinagpu_total.csv
  python enriquecer_smiles.py ../gpu_dock/tanda_z001/resultados_z001.csv
"""

import argparse
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from filtro_rescoring import cargar_smiles_desde_librerias, CARPETA_LIBRERIAS

BASE_CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
SALIDA = os.path.join(CARPETA_LIBRERIAS, "libreria_extra_chembl.csv")


def obtener_smiles_por_id(chembl_id):
    """Descarga el SMILES de una molecula individual (endpoint directo)."""
    url = f"{BASE_CHEMBL}/molecule/{chembl_id}.json"
    try:
        with urllib.request.urlopen(url, timeout=30) as r:
            datos = json.load(r)
        estructuras = datos.get("molecule_structures") or {}
        smiles = (estructuras.get("canonical_smiles")
                  or estructuras.get("standard_smiles"))
        return (chembl_id, smiles)
    except Exception:
        return (chembl_id, None)


def obtener_smiles_multi(ids_chembl, workers=4):
    """Descarga SMILES con 4 hilos, respetando el limite de la API publica."""
    resultados = {}
    total = len(ids_chembl)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for i, (chembl_id, smiles) in enumerate(pool.map(obtener_smiles_por_id,
                                                         ids_chembl), 1):
            if smiles:
                resultados[chembl_id] = smiles
            if i % 250 == 0 or i == total:
                print(f"  {i}/{total} procesados, {len(resultados)} SMILES")
    return resultados


def buscar_smiles_por_nombre(nombre):
    """Busca un nombre en ChEMBL y devuelve (chembl_id, smiles) del primer hit."""
    url = f"{BASE_CHEMBL}/molecule/search?q={urllib.parse.quote(nombre)}&format=json"
    try:
        with urllib.request.urlopen(url, timeout=60) as r:
            datos = json.load(r)
        for mol in datos.get("molecules", []):
            estructuras = mol.get("molecule_structures") or {}
            smiles = (estructuras.get("canonical_smiles")
                      or estructuras.get("standard_smiles"))
            if smiles:
                return (mol.get("molecule_chembl_id"), smiles)
    except Exception as e:
        print(f"  Error buscando '{nombre}': {e}")
    return (None, None)


def main():
    parser = argparse.ArgumentParser(description="Completa SMILES via API de ChEMBL")
    parser.add_argument("csv_entrada", help="CSV de resultados de docking")
    parser.add_argument("--top-pct", type=float, default=None,
                        help="Solo ligandos del percentil superior (ej. 5.0). "
                             "Por defecto: todos los del CSV")
    parser.add_argument("--max-nombres", type=int, default=200,
                        help="Limite de busquedas por nombre (default: 200)")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_entrada)
    col_lig = "ligand" if "ligand" in df.columns else "compound_id"
    if args.top_pct:
        col_af = next((c for c in ("affinity", "energy") if c in df.columns), None)
        if col_af is None:
            print("ERROR: no encuentro columna de afinidad (affinity/energy).")
            sys.exit(1)
        umbral = df[col_af].quantile(args.top_pct / 100)
        df = df[df[col_af] <= umbral]
        print(f"Filtrado al top {args.top_pct}% (umbral {umbral:.2f} kcal/mol): "
              f"{len(df)} filas")
    ligandos = df[col_lig].astype(str).str.strip().unique()
    print(f"Cargado '{args.csv_entrada}': {len(ligandos)} ligandos unicos")

    mapa_local = cargar_smiles_desde_librerias(CARPETA_LIBRERIAS)
    faltan = [lig for lig in ligandos if lig.upper() not in mapa_local]
    print(f"{len(faltan)} ligandos sin SMILES en las librerias locales")

    ids_chembl = sorted({lig for lig in faltan if lig.upper().startswith("CHEMBL")})
    nombres = [lig for lig in faltan if not lig.upper().startswith("CHEMBL")]
    print(f"  CHEMBL a buscar por API: {len(ids_chembl)}")
    print(f"  Nombres a buscar por API: {len(nombres)}")

    encontrados = {}

    if ids_chembl:
        encontrados.update(obtener_smiles_multi(ids_chembl))

    for nombre in nombres[:args.max_nombres]:
        chembl_id, smiles = buscar_smiles_por_nombre(nombre)
        if smiles:
            encontrados[nombre.upper()] = smiles
        time.sleep(0.25)

    if not encontrados:
        print("No se encontraron SMILES nuevos. Nada que guardar.")
        return

    filas = [{"name": k, "chembl_id": (k if k.startswith("CHEMBL") else ""),
              "canonical_smiles": v} for k, v in encontrados.items()]
    extra = pd.DataFrame(filas)
    extra.to_csv(SALIDA, index=False)
    print(f"{len(extra)} SMILES guardados en '{SALIDA}'")
    print(f"Faltantes tras enriquecer: {len(faltan) - len(encontrados)}")


if __name__ == "__main__":
    main()
