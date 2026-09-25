# -*- coding: utf-8 -*-
"""
Lista corta de candidatos MASIVE-ALS.

Embudo sobre analysis/ranking_top5_consolidado.csv (top-5% por proteina):
  1) SMILES obligatorio
  2) PAINS (pan-assay interference)  -> elimina falsos positivos promiscuos
  3) Lipinski (MW<=500, logP<=5, HBD<=5, HBA<=10)
  4) Veber (rotb<=10, TPSA<=140)
  5) CNS: cns_mpo_ok (MPO>=4)  [regla CNS de quimica medicinal]
  6) Ranking final por afinidad dentro de cada proteina

Salida: analysis/lista_corta_candidatos.csv
"""
import os
import csv

from rdkit import Chem, RDLogger
from rdkit.Chem import FilterCatalog

RDLogger.DisableLog("rdApp.*")

BASE = r"C:\Users\Fredy\masive-als"
ANALYSIS = os.path.join(BASE, "analysis")
SRC = os.path.join(ANALYSIS, "ranking_top5_consolidado.csv")
OUT = os.path.join(ANALYSIS, "lista_corta_candidatos.csv")

# catalogo PAINS de RDKit
_params = FilterCatalog.FilterCatalogParams()
_params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_CAT = FilterCatalog.FilterCatalog(_params)


def pasa_pains(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return False, "smiles_invalido"
    if _CAT.HasMatch(m):
        return False, "PAINS"
    return True, ""


def main():
    rows = list(csv.DictReader(open(SRC, encoding="utf-8")))
    print("top5 total:", len(rows))

    ok = []
    stats = {"sin_smiles": 0, "pains": 0, "lipinski": 0, "veber": 0, "cns": 0}
    for r in rows:
        smi = (r.get("smiles") or "").strip()
        if not smi:
            stats["sin_smiles"] += 1
            continue
        pasa, motivo = pasa_pains(smi)
        if not pasa:
            stats["pains"] += 1
            continue
        try:
            mw = float(r["mw"]); logp = float(r["logp"])
            hbd = int(r["hbd"]); hba = int(r["hba"])
            rotb = int(r["rotb"]); tpsa = float(r["tpsa"])
        except (ValueError, KeyError):
            continue
        if not (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10):
            stats["lipinski"] += 1
            continue
        if not (rotb <= 10 and tpsa <= 140):
            stats["veber"] += 1
            continue
        if r.get("cns_mpo_ok") != "True":
            stats["cns"] += 1
            continue
        ok.append(r)

    print("stats:", stats)
    print("lista corta:", len(ok))

    ok.sort(key=lambda r: (r["target"], float(r["afinidad"])))
    with open(OUT, "w", newline="", encoding="utf-8") as f:
        cols = ["target", "ligand", "afinidad", "smiles", "mw", "tpsa",
                "logp", "hbd", "hba", "rotb", "cns_mpo"]
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(ok)

    # resumen por proteina
    from collections import Counter
    c = Counter(r["target"] for r in ok)
    for t in ["TDP43_v2", "SOD1", "FUS"]:
        print("  %-10s %d" % (t, c.get(t, 0)))
    print("CSV:", OUT)


if __name__ == "__main__":
    main()