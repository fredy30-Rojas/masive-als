# -*- coding: utf-8 -*-
"""cns_filtro.py — Post-proceso del embudo MASIVE-ALS (revisión científica 18/08).

Une los candidatos POR PROTEÍNA (umbrales top-5% propios) y aplica:
  1. Descriptores RDKit: TPSA, MW, logP, HBD, HBA, RotB, CNS-MPO (Wager 2010, simplificado).
  2. Filtro CNS: TPSA <= 90 A^2 Y MW <= 450 (criterio orientativo CNS+;
     alternativamente CNS-MPO >= 4).
  3. Columna 'aprobado' (si/no) comparando contra la librería FDA completa.
  4. Columna 'cns' (estimado: si/no) y 'cns_mpo'.
Salida: candidatos_filtrados.csv (reemplaza la versión sesgada global).
"""
import argparse
import os
import sys

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

CARPETA_LIBRERIAS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "compounds")
FDA_CSV = os.path.join(CARPETA_LIBRERIAS, "full_fda_library.csv")


def descriptores(smiles: str) -> dict:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    return {
        "mw": Descriptors.MolWt(mol),
        "tpsa": Descriptors.TPSA(mol),
        "logp": Descriptors.MolLogP(mol),
        "hbd": Descriptors.NumHDonors(mol),
        "hba": Descriptors.NumHAcceptors(mol),
        "rotb": Descriptors.NumRotatableBonds(mol),
    }


def cns_mpo(mw: float, logp: float, tpsa: float, hbd: int) -> float:
    """CNS-MPO simplificado (Wager et al., ACS Chem Neurosci 2010), 0-6.
    Componentes: clogP, MW, TPSA, HBD, pKa (omitido -> 1.0 neutral).
    Valores >= 4 se consideran favorables para penetración CNS.
    """
    def score_lin(valor, min_ok, max_ok, min_no, max_no):
        # puntúa 1 en [min_ok,max_ok], 0 en [min_no,max_no], lineal entre ambos
        if min_no <= valor <= max_no:
            return 0.0
        if min_ok <= valor <= max_ok:
            return 1.0
        if valor < min_no:
            return 1.0 - (min_no - valor) / max(min_no - min_ok, 1e-9) * 0.0 + 0.0
        return 1.0 - (valor - max_no) / max(max_ok - max_no, 1e-9) * 0.0 + 0.0

    # Wager: ClogP 2-4 ok; MW 250-350 ok; TPSA 40-90 ok; HBD 0-2 ok
    s_logp = 1.0 - max(0.0, min(1.0, (logp - 4.0) / 2.0)) - max(0.0, min(1.0, (2.0 - logp) / 2.0))
    s_mw = 1.0 - max(0.0, min(1.0, (mw - 350.0) / 100.0)) - max(0.0, min(1.0, (250.0 - mw) / 100.0))
    s_tpsa = 1.0 - max(0.0, min(1.0, (tpsa - 90.0) / 40.0)) - max(0.0, min(1.0, (40.0 - tpsa) / 40.0))
    s_hbd = 1.0 - max(0.0, min(1.0, (hbd - 2.0) / 2.0))
    return round(s_logp + s_mw + s_tpsa + s_hbd + 1.0, 2)  # +1 por pKa asumido neutro


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cns-mpo-umbral", type=float, default=4.0,
                        help="Umbral CNS-MPO (default 4.0)")
    parser.add_argument("--salida", default="candidatos_filtrados.csv")
    parser.add_argument("--prefijo", default="candidatos",
                        help="Prefijo de los CSVs por proteína (candidatos / candidatos_total)")
    args = parser.parse_args()

    # 1) cargar candidatos por proteína
    partes = []
    for prot in ["TDP43", "SOD1", "FUS"]:
        ruta = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            "%s_%s.csv" % (args.prefijo, prot))
        if not os.path.exists(ruta):
            print("Falta %s: ejecuta primero filtro_rescoring.py --proteina %s" % (ruta, prot))
            sys.exit(1)
        df = pd.read_csv(ruta)
        partes.append(df)
        print("Cargado %s: %d candidatos" % (ruta, len(df)))
    df = pd.concat(partes, ignore_index=True)
    print(f"Total por proteína (top-5% propio + PAINS + Lipinski/Veber): {len(df)}")

    # 2) descriptores
    desc = df["smiles"].apply(lambda s: descriptores(s) if isinstance(s, str) else {})
    for col in ["mw", "tpsa", "logp", "hbd", "hba", "rotb"]:
        df[col] = [d.get(col) if isinstance(d, dict) else None for d in desc]
    df["cns_mpo"] = [
        cns_mpo(r.mw, r.logp, r.tpsa, r.hbd)
        if pd.notna(r.mw) else None
        for r in df.itertuples()
    ]

    # 3) filtro CNS (TPSA <= 90 y MW <= 450) — criterio principal
    n_antes = len(df)
    df["cns"] = (df["tpsa"] <= 90) & (df["mw"] <= 450)
    df = df[df["cns"]].copy()
    print(f"Filtro CNS (TPSA<=90, MW<=450): {n_antes} -> {len(df)}")
    df["cns_mpo_ok"] = df["cns_mpo"] >= args.cns_mpo_umbral

    # 4) columna aprobado (si/no)
    aprobados = set()
    if os.path.exists(FDA_CSV):
        fda = pd.read_csv(FDA_CSV)
        for col in ["chembl_id", "name"]:
            if col in fda.columns:
                for v in fda[col].dropna():
                    aprobados.add(str(v).strip().upper())
    df["aprobado"] = df["ligand"].apply(
        lambda x: "si" if str(x).strip().upper() in aprobados else "no")
    n_aprobados = (df["aprobado"] == "si").sum()
    print(f"Candidatos aprobados (FDA): {n_aprobados} de {len(df)}")

    # 5) ordenar: aprobados primero, luego mejor afinidad; guardar
    col_af = "affinity" if "affinity" in df.columns else "energy"
    df = df.sort_values(["aprobado", col_af], ascending=[True, True])
    cols = ["ligand", "target", col_af, "smiles", "aprobado",
            "cns", "cns_mpo", "cns_mpo_ok", "mw", "tpsa", "logp", "hbd", "hba", "rotb"]
    df[cols].to_csv(args.salida, index=False)
    print(f"Exportado: {args.salida} ({len(df)} candidatos)")


if __name__ == "__main__":
    main()