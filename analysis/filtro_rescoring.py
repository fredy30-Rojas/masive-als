# -*- coding: utf-8 -*-
"""
filtro_rescoring.py
MASIVE-ALS — Embudo de filtrado y re-scoring en 3 pasos (adaptado 17 ago 2026)

Reduce miles de resultados de docking (Vina/Vina-GPU) a un subconjunto
pequeño y de alta confianza, listo para MM-GBSA / dinamica molecular.

Pasos:
  1) Corte por percentil de afinidad (mas negativo = mejor)
  2) Union con SMILES desde las librerias de compuestos, dedup por
     SMILES canonico (por proteina), filtro PAINS + Lipinski/Veber
  3) Exporta la lista final para preparar corridas de MM-GBSA

CSV de resultados soportados (se detecta la columna de afinidad sola):
  - gpu_dock/tanda_z001/resultados_z001.csv        -> target,ligand,affinity
  - gpu_dock/resultados_vinagpu_total.csv          -> ligand,target,energy

Las librerias de compuestos se buscan solas en ../compounds (CSV con
columna smiles e identificador name/chembl_id). Requiere pandas y rdkit:
  pip install rdkit-pypi pandas --break-system-packages

Uso:
  python filtro_rescoring.py ../gpu_dock/resultados_vinagpu_total.csv \
      --top-pct 5 --salida candidatos_filtrados.csv
  python filtro_rescoring.py ../gpu_dock/tanda_z001/resultados_z001.csv \
      --proteina TDP-43 --top-pct 2 --salida candidatos_tdp43.csv
"""

import argparse
import os
import sys

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, FilterCatalog

RDLogger.DisableLog("rdApp.*")

CARPETA_LIBRERIAS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "compounds")
COL_COMPUESTO = "ligand"
COL_TARGET = "target"
COL_AFINIDAD_POSIBLES = ("affinity", "energy")
COL_SMILES_SAIDA = "smiles"
# Nombres reales en las librerias locales de compuestos:
COL_ID_POSIBLES = ("molecule_chembl_id", "chembl_id", "name", "molecule_id", "zinc_id")
COL_SMILES_POSIBLES = ("canonical_smiles", "smiles")


def detectar_col_afinidad(df: pd.DataFrame) -> str:
    """Devuelve la columna de afinidad que exista en el CSV real."""
    for col in COL_AFINIDAD_POSIBLES:
        if col in df.columns:
            return col
    print("ERROR: no encuentro columna de afinidad (affinity/energy) en el CSV.")
    sys.exit(1)


def cargar_smiles_desde_librerias(carpeta: str) -> dict:
    """Escanea los CSV de la carpeta compounds y arma {id_compuesto: smiles}.

    Acepta como identificador las columnas name, chembl_id, molecule_id o
    zinc_id (la primera que exista), siempre que el CSV tenga columna smiles.
    """
    mapa = {}
    if not os.path.isdir(carpeta):
        print(f"Aviso: no existe la carpeta de librerias '{carpeta}'.")
        return mapa
    for raiz, _, archivos in os.walk(carpeta):
        for archivo in archivos:
            if not archivo.lower().endswith(".csv"):
                continue
            ruta = os.path.join(raiz, archivo)
            try:
                df = pd.read_csv(ruta, nrows=200000)
            except Exception:
                continue
            col_smiles = next((c for c in COL_SMILES_POSIBLES if c in df.columns), None)
            if col_smiles is None:
                continue
            cols_id = [c for c in COL_ID_POSIBLES if c in df.columns]
            if not cols_id:
                continue
            for _, fila in df.iterrows():
                sonrisas = fila[col_smiles]
                if not isinstance(sonrisas, str):
                    continue
                # Registrar TODAS las claves disponibles (name y chembl_id),
                # porque los CSV de resultados usan nombres y/o IDs
                for col_id in cols_id:
                    ident = str(fila[col_id]).strip().upper()
                    if ident and ident != "NAN":
                        mapa.setdefault(ident, sonrisas.strip())
    print(f"Librerias: {len(mapa)} compuestos con SMILES cargados desde '{carpeta}'")
    return mapa


def paso1_corte_por_afinidad(df: pd.DataFrame, col_af: str, top_pct: float) -> pd.DataFrame:
    """Se queda solo con el percentil superior (mas negativo) de afinidad."""
    n_total = len(df)
    umbral = df[col_af].quantile(top_pct / 100)
    filtrado = df[df[col_af] <= umbral].copy()
    print(f"[Paso 1] {n_total} -> {len(filtrado)} filas "
          f"(top {top_pct}%, umbral = {umbral:.2f} kcal/mol)")
    return filtrado


def paso2_join_smiles_dedup_pains(df: pd.DataFrame, col_af: str,
                                  mapa_smiles: dict) -> pd.DataFrame:
    """Une SMILES desde la libreria, deduplica y aplica PAINS + Lipinski/Veber."""
    # 2a) Unir SMILES y reportar los que no se encuentran
    def obtener_smiles(ident):
        clave = str(ident).strip().upper()
        return mapa_smiles.get(clave, None)

    df = df.copy()
    df[COL_SMILES_SAIDA] = df[COL_COMPUESTO].apply(obtener_smiles)
    sin_smiles = df[df[COL_SMILES_SAIDA].isna()]
    con_smiles = df[df[COL_SMILES_SAIDA].notna()].copy()
    print(f"[Paso 2] {len(sin_smiles)} filas sin SMILES en las librerias "
          f"(descartadas, no se pueden filtrar ni preparar para MM-GBSA)")

    # 2b) Dedup por SMILES canonico y proteina, conservando la mejor afinidad
    def canonizar(sonrisas):
        mol = Chem.MolFromSmiles(sonrisas)
        return Chem.MolToSmiles(mol) if mol else None

    con_smiles["smiles_canonico"] = con_smiles[COL_SMILES_SAIDA].apply(canonizar)
    con_smiles = con_smiles.dropna(subset=["smiles_canonico"])
    antes = len(con_smiles)
    con_smiles = (con_smiles.sort_values(col_af)
                            .drop_duplicates(subset=[COL_TARGET, "smiles_canonico"])
                            .sort_values(col_af))
    print(f"[Paso 2] dedup por SMILES+proteina: {antes} -> {len(con_smiles)} filas")

    # 2c) PAINS + Lipinski/Veber
    catalogo_params = FilterCatalog.FilterCatalogParams()
    catalogo_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalogo = FilterCatalog.FilterCatalog(catalogo_params)

    filas_validas = []
    for _, fila in con_smiles.iterrows():
        mol = Chem.MolFromSmiles(fila[COL_SMILES_SAIDA])
        if mol is None:
            continue
        if catalogo.HasMatch(mol):
            continue
        peso = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        donores = Descriptors.NumHDonors(mol)
        aceptores = Descriptors.NumHAcceptors(mol)
        rotables = Descriptors.NumRotatableBonds(mol)
        cumple_lipinski = (peso <= 500 and logp <= 5 and donores <= 5 and aceptores <= 10)
        cumple_veber = rotables <= 10
        if cumple_lipinski and cumple_veber:
            filas_validas.append(fila)

    resultado = pd.DataFrame(filas_validas)
    print(f"[Paso 2] {len(con_smiles)} -> {len(resultado)} compuestos "
          f"(pasan PAINS + Lipinski/Veber)")
    return resultado


def paso3_preparar_para_mmgbsa(df: pd.DataFrame, col_af: str, ruta_salida: str) -> None:
    """Exporta la lista final para preparar corridas de MM-GBSA en GROMACS."""
    columnas_salida = [c for c in [COL_COMPUESTO, COL_TARGET, col_af, COL_SMILES_SAIDA]
                       if c in df.columns]
    df[columnas_salida].to_csv(ruta_salida, index=False)
    print(f"[Paso 3] {len(df)} candidatos exportados a '{ruta_salida}' "
          f"para preparar MM-GBSA (ver prepare_md.py)")


def main():
    parser = argparse.ArgumentParser(description="Embudo de filtrado MASIVE-ALS")
    parser.add_argument("csv_entrada", help="CSV de resultados de docking")
    parser.add_argument("--proteina", default=None,
                        help="Filtrar solo esta proteina (ej. TDP-43, SOD1, FUS)")
    parser.add_argument("--top-pct", type=float, default=5.0,
                        help="Percentil superior de afinidad a conservar (default: 5.0)")
    parser.add_argument("--librerias", default=CARPETA_LIBRERIAS,
                        help="Carpeta con los CSV de librerias de compuestos")
    parser.add_argument("--salida", default="candidatos_filtrados.csv",
                        help="Nombre del CSV de salida")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_entrada)
    print(f"Cargado '{args.csv_entrada}': {len(df)} filas totales")
    col_af = detectar_col_afinidad(df)

    if args.proteina and COL_TARGET in df.columns:
        df = df[df[COL_TARGET] == args.proteina]
        print(f"Filtrado a proteina '{args.proteina}': {len(df)} filas")

    df = paso1_corte_por_afinidad(df, col_af, args.top_pct)
    mapa_smiles = cargar_smiles_desde_librerias(args.librerias)
    df = paso2_join_smiles_dedup_pains(df, col_af, mapa_smiles)
    paso3_preparar_para_mmgbsa(df, col_af, args.salida)


if __name__ == "__main__":
    main()
