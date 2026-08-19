# -*- coding: utf-8 -*-
"""
filtro_rescoring.py  (v3 — corrige bloqueantes detectados por Claude CLI 18/08)
MASIVE-ALS — Embudo de filtrado y re-scoring en 3 pasos

Reduce miles de resultados de docking (Vina/Vina-GPU) a un subconjunto
pequeno y de alta confianza, listo para MM-GBSA / dinamica molecular.

CAMBIOS v3 respecto a v2 (bloqueantes reportados 18/08):
1. Columnas mapeadas automaticamente contra los nombres reales de los CSV
   del proyecto (target/ligand/affinity, ligand/target/energy), en lugar
   de asumir compound_id/protein/affinity_kcal_mol fijos.
2. Restaurada la union de SMILES desde una libreria externa (--lib-smiles):
   los CSV de resultados NO traen smiles, hay que unirlos por id de ligando.
3. Restaurada la deduplicacion por SMILES canonico DENTRO de cada proteina,
   quedandose con la mejor (mas negativa) afinidad cuando el mismo
   compuesto aparece con distintos identificadores (name vs chembl_id).
4. El enrichment de controles positivos ahora se calcula sobre la seleccion
   CRUDA del top-pct por proteina (justo despues del paso 1), ANTES de los
   filtros de quimica medicinal — si se calculara despues, un control que
   cae por PAINS/Lipinski/CNS se contaria como "no encontrado" y mentiria.
5. Las filas sin SMILES (o SMILES invalido) ahora se CUENTAN y se avisan
   explicitamente en cada paso, en vez de descartarse en silencio.
6. Restaurado --proteina para poder correr una sola proteina si hace falta.

ADAPTAR antes de usar:
- Ajusta ALIAS_COLUMNAS si tus CSV usan nombres distintos a los ya
  contemplados (target/ligand/affinity/energy/protein/compound_id/smiles).
- Requiere: pandas, rdkit (pip install rdkit-pypi pandas --break-system-packages)

Uso:
    # Embudo completo, une SMILES desde la libreria, por cada proteina
    python filtro_rescoring.py resultados_z001.csv --lib-smiles compounds/libreria_smiles.csv --top-pct 5

    # Solo una proteina, con controles positivos para enrichment
    python filtro_rescoring.py resultados_z001.csv --lib-smiles compounds/libreria_smiles.csv \
        --proteina TDP43 --controles controles_conocidos.csv
"""

import argparse
import sys

import pandas as pd

# ---- Alias de columnas: cubre los nombres reales usados en el proyecto ----
# Cada lista son los nombres posibles que puede tener esa columna en el CSV
# de entrada; se detecta automaticamente cual esta presente.
ALIAS_COLUMNAS = {
    "compuesto": ["ligand", "compound_id", "compound", "ligand_id"],
    "afinidad": ["affinity", "energy", "affinity_kcal_mol"],
    "proteina": ["target", "protein"],
    "smiles": ["smiles", "SMILES"],
    "pose": ["pose_pdbqt_path", "pose_path", "pdbqt_path"],
}

CNS_TPSA_MAX = 90.0
CNS_MW_MAX = 450.0


def detectar_columnas(df: pd.DataFrame) -> dict:
    """Detecta que nombre real tiene cada columna logica en este CSV."""
    columnas = {}
    for logica, alias in ALIAS_COLUMNAS.items():
        encontrada = next((a for a in alias if a in df.columns), None)
        columnas[logica] = encontrada
    faltantes_criticas = [k for k in ("compuesto", "afinidad", "proteina") if columnas[k] is None]
    if faltantes_criticas:
        print(f"ERROR: no se encontraron columnas para: {faltantes_criticas}")
        print(f"Columnas disponibles en el CSV: {list(df.columns)}")
        print("Ajusta ALIAS_COLUMNAS en el script si tu CSV usa otros nombres.")
        sys.exit(1)
    print(f"Columnas detectadas -> compuesto: '{columnas['compuesto']}', "
          f"afinidad: '{columnas['afinidad']}', proteina: '{columnas['proteina']}', "
          f"smiles: '{columnas['smiles']}', pose: '{columnas['pose']}'")
    return columnas


def unir_smiles_desde_libreria(df: pd.DataFrame, cols: dict, ruta_lib: str) -> pd.DataFrame:
    """
    Une SMILES a los resultados desde una libreria externa (--lib-smiles),
    ya que los CSV de resultados de docking normalmente no traen SMILES.

    La libreria debe tener al menos: una columna de id de compuesto
    (mismo valor que aparece en la columna 'compuesto' de resultados) y
    una columna de smiles.
    """
    lib = pd.read_csv(ruta_lib)
    id_lib = next((a for a in ALIAS_COLUMNAS["compuesto"] if a in lib.columns), None)
    smiles_lib = next((a for a in ALIAS_COLUMNAS["smiles"] if a in lib.columns), None)
    if id_lib is None or smiles_lib is None:
        print(f"ERROR: la libreria '{ruta_lib}' no tiene columnas de id/smiles "
              f"reconocibles. Columnas encontradas: {list(lib.columns)}")
        sys.exit(1)

    lib_reducida = lib[[id_lib, smiles_lib]].drop_duplicates(subset=id_lib)
    lib_reducida = lib_reducida.rename(columns={id_lib: cols["compuesto"], smiles_lib: "smiles"})

    antes = len(df)
    df = df.merge(lib_reducida, on=cols["compuesto"], how="left")
    sin_smiles = df["smiles"].isna().sum()
    print(f"[Union SMILES] {antes} filas de resultados unidas contra '{ruta_lib}' "
          f"({len(lib_reducida)} compuestos con smiles). "
          f"{sin_smiles} filas quedaron SIN smiles tras la union (no se descartan aqui, "
          f"se cuentan explicitamente en los pasos siguientes).")
    cols["smiles"] = "smiles"
    return df


def deduplicar_por_smiles_canonico(df: pd.DataFrame, cols: dict) -> pd.DataFrame:
    """
    Deduplica por SMILES canonico DENTRO de cada proteina, quedandose con
    la fila de mejor (mas negativa) afinidad. Evita que el mismo compuesto
    registrado con distintos identificadores (nombre vs ChEMBL ID) cuente
    dos veces e infle el top-pct / enrichment.
    """
    try:
        from rdkit import Chem
    except ImportError:
        print("RDKit no esta instalado. Instala con:")
        print("  pip install rdkit-pypi --break-system-packages")
        sys.exit(1)

    def canonizar(smiles):
        if not isinstance(smiles, str):
            return None
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)

    df = df.copy()
    df["_smiles_canonico"] = df[cols["smiles"]].apply(canonizar)

    invalidos = df["_smiles_canonico"].isna().sum()
    con_smiles = df.dropna(subset=["_smiles_canonico"])

    antes = len(con_smiles)
    con_smiles = con_smiles.sort_values(cols["afinidad"])  # mas negativo primero
    deduplicado = con_smiles.drop_duplicates(subset=[cols["proteina"], "_smiles_canonico"], keep="first")

    print(f"[Dedup] {invalidos} filas con SMILES ausente/invalido (excluidas de dedup, "
          f"quedan fuera del embudo). De las {antes} filas con SMILES valido: "
          f"{len(deduplicado)} tras deduplicar por SMILES canonico + proteina "
          f"(se descartaron {antes - len(deduplicado)} duplicados, quedandose siempre "
          f"con la mejor afinidad).")
    return deduplicado.drop(columns=["_smiles_canonico"])


def paso1_corte_por_afinidad_por_proteina(df: pd.DataFrame, cols: dict, top_pct: float,
                                            solo_proteina: str = None) -> pd.DataFrame:
    """Corte por percentil superior de afinidad, SIEMPRE calculado por proteina."""
    if solo_proteina:
        df = df[df[cols["proteina"]] == solo_proteina]
        print(f"[Paso 1] Restringido a proteina '{solo_proteina}': {len(df)} filas")

    partes = []
    for proteina, grupo in df.groupby(cols["proteina"]):
        umbral = grupo[cols["afinidad"]].quantile(top_pct / 100)
        filtrado = grupo[grupo[cols["afinidad"]] <= umbral].copy()
        print(f"[Paso 1] {proteina}: {len(grupo)} -> {len(filtrado)} "
              f"(top {top_pct}%, umbral propio = {umbral:.2f})")
        partes.append(filtrado)

    resultado = pd.concat(partes, ignore_index=True) if partes else df.iloc[0:0]
    print(f"[Paso 1] Total combinado: {len(resultado)} compuestos "
          f"(umbral independiente por proteina)")
    return resultado


def calcular_enrichment(df_top_crudo: pd.DataFrame, cols: dict, ruta_controles: str) -> None:
    """
    Calcula enrichment SOBRE LA SELECCION CRUDA del top-pct (justo despues
    del paso 1, antes de PAINS/Lipinski/CNS). Calcularlo despues de esos
    filtros haria que un control valido que cae por quimica medicinal se
    contara como "no encontrado", falseando la metrica de validacion.
    """
    controles = pd.read_csv(ruta_controles)
    id_ctrl = next((a for a in ALIAS_COLUMNAS["compuesto"] if a in controles.columns), None)
    prot_ctrl = next((a for a in ALIAS_COLUMNAS["proteina"] if a in controles.columns), None)
    if id_ctrl is None or prot_ctrl is None:
        print(f"ERROR: el CSV de controles no tiene columnas reconocibles de "
              f"compuesto/proteina. Columnas: {list(controles.columns)}")
        return

    ids_top = set(zip(df_top_crudo[cols["proteina"]], df_top_crudo[cols["compuesto"]]))
    print("\n[Validacion] Enrichment de controles positivos "
          "(calculado ANTES de filtros de quimica medicinal):")
    for proteina, grupo in controles.groupby(prot_ctrl):
        total = len(grupo)
        encontrados = sum(1 for _, fila in grupo.iterrows()
                           if (proteina, fila[id_ctrl]) in ids_top)
        print(f"  {proteina}: {encontrados}/{total} controles conocidos "
              f"aparecen en el top crudo por proteina")
    print("  (si esto da 0 para alguna proteina, revisar caja/bolsillo de esa "
          "proteina antes de confiar en su ranking)\n")


def paso2_drug_likeness_pains(df: pd.DataFrame, cols: dict) -> pd.DataFrame:
    """Filtra por reglas de Lipinski/Veber y catalogo PAINS via RDKit."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, FilterCatalog

    catalogo_params = FilterCatalog.FilterCatalogParams()
    catalogo_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
    catalogo = FilterCatalog.FilterCatalog(catalogo_params)

    filas_validas = []
    sin_smiles = 0
    rechazadas_pains = 0
    rechazadas_lipinski = 0

    for _, fila in df.iterrows():
        smiles = fila.get(cols["smiles"])
        if not isinstance(smiles, str):
            sin_smiles += 1
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            sin_smiles += 1
            continue

        if catalogo.HasMatch(mol):
            rechazadas_pains += 1
            continue

        peso_molecular = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        donores_h = Descriptors.NumHDonors(mol)
        aceptores_h = Descriptors.NumHAcceptors(mol)
        enlaces_rotables = Descriptors.NumRotatableBonds(mol)

        cumple = (peso_molecular <= 500 and logp <= 5 and donores_h <= 5
                  and aceptores_h <= 10 and enlaces_rotables <= 10)
        if cumple:
            filas_validas.append(fila)
        else:
            rechazadas_lipinski += 1

    resultado = pd.DataFrame(filas_validas)
    print(f"[Paso 2] {len(df)} -> {len(resultado)} compuestos. "
          f"Descartadas: {sin_smiles} sin SMILES valido, "
          f"{rechazadas_pains} por PAINS, {rechazadas_lipinski} por Lipinski/Veber.")
    return resultado


def paso2b_filtro_cns_bbb(df: pd.DataFrame, cols: dict) -> pd.DataFrame:
    """Filtro orientativo de penetracion al sistema nervioso central (TPSA, MW)."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors

    filas_validas = []
    sin_smiles = 0
    rechazadas_cns = 0

    for _, fila in df.iterrows():
        smiles = fila.get(cols["smiles"])
        if not isinstance(smiles, str):
            sin_smiles += 1
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            sin_smiles += 1
            continue

        tpsa = Descriptors.TPSA(mol)
        peso_molecular = Descriptors.MolWt(mol)
        if tpsa <= CNS_TPSA_MAX and peso_molecular <= CNS_MW_MAX:
            filas_validas.append(fila)
        else:
            rechazadas_cns += 1

    resultado = pd.DataFrame(filas_validas)
    print(f"[Paso 2b - CNS/BBB] {len(df)} -> {len(resultado)} compuestos. "
          f"Descartadas: {sin_smiles} sin SMILES valido, "
          f"{rechazadas_cns} por no cumplir TPSA<={CNS_TPSA_MAX} / MW<={CNS_MW_MAX}.")
    return resultado


def paso3_preparar_para_mmgbsa(df: pd.DataFrame, cols: dict, ruta_salida: str) -> None:
    """Exporta la lista final para preparar corridas de MM-GBSA en GROMACS."""
    columnas_salida = [c for c in
                        [cols["compuesto"], cols["proteina"], cols["afinidad"],
                         cols["smiles"], cols["pose"]] if c]
    if not cols["pose"]:
        print("[Aviso] No hay columna de ruta a la pose acoplada en el CSV de entrada: "
              "el resultado exportado NO estara listo para MM-GBSA directamente, "
              "solo para revision de candidatos.")

    df[columnas_salida].to_csv(ruta_salida, index=False)
    print(f"[Paso 3] {len(df)} compuestos exportados a '{ruta_salida}'")


def main():
    parser = argparse.ArgumentParser(description="Embudo de filtrado MASIVE-ALS (v3)")
    parser.add_argument("csv_entrada", help="CSV de resultados de docking")
    parser.add_argument("--lib-smiles", default=None,
                         help="CSV de libreria con id de compuesto + smiles, para unir "
                              "si el CSV de entrada no trae smiles")
    parser.add_argument("--proteina", default=None,
                         help="Restringir a una sola proteina (ej. TDP43, SOD1, FUS)")
    parser.add_argument("--top-pct", type=float, default=5.0,
                         help="Percentil superior de afinidad a conservar POR PROTEINA "
                              "(default: 5.0)")
    parser.add_argument("--salida", default="candidatos_filtrados.csv",
                         help="Nombre del CSV de salida")
    parser.add_argument("--controles", default=None,
                         help="CSV opcional con ligandos conocidos para calcular enrichment")
    args = parser.parse_args()

    df = pd.read_csv(args.csv_entrada)
    print(f"Cargado '{args.csv_entrada}': {len(df)} filas totales")

    cols = detectar_columnas(df)

    if not cols["smiles"]:
        if not args.lib_smiles:
            print("ERROR: el CSV de entrada no tiene columna de smiles, y no se paso "
                  "--lib-smiles para unirla. El script no puede continuar sin SMILES.")
            sys.exit(1)
        df = unir_smiles_desde_libreria(df, cols, args.lib_smiles)

    df = deduplicar_por_smiles_canonico(df, cols)

    df_top_crudo = paso1_corte_por_afinidad_por_proteina(df, cols, args.top_pct, args.proteina)

    if args.controles:
        calcular_enrichment(df_top_crudo, cols, args.controles)

    df_filtrado = paso2_drug_likeness_pains(df_top_crudo, cols)
    df_filtrado = paso2b_filtro_cns_bbb(df_filtrado, cols)

    paso3_preparar_para_mmgbsa(df_filtrado, cols, args.salida)


if __name__ == "__main__":
    main()
