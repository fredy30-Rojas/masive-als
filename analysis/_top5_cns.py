# -*- coding: utf-8 -*-
"""
Top-5% por proteina + propiedades moleculares + filtro CNS/BBB.

Fuente: analysis/ranking_afinidades.csv (target, ligand, afinidad)
generado desde los *_out.pdbqt del GPU.

Salidas:
  analysis/ranking_top5_consolidado.csv  -> top-5% con props y filtros
  analysis/ranking_consolidado_resumen.json
"""
import os
import csv
import json
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

BASE = r"C:\Users\Fredy\masive-als"
ANALYSIS = os.path.join(BASE, "analysis")
AFI = os.path.join(ANALYSIS, "ranking_afinidades.csv")
EXCL_PATH = os.path.join(BASE, "gpu_dock", "ligandos_excluidos.txt")

TARGETS = ["TDP43_v2", "SOD1", "FUS"]


def props(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return {}
    mw = Descriptors.MolWt(m)
    tpsa = Descriptors.TPSA(m)
    logp = Crippen.MolLogP(m)
    hbd = Lipinski.NumHDonors(m)
    hba = Lipinski.NumHAcceptors(m)
    rotb = Lipinski.NumRotatableBonds(m)
    mpo = 0.0
    if logp <= 5: mpo += 1
    if 40 <= tpsa <= 90: mpo += 1
    if mw <= 450: mpo += 1
    if hbd <= 3: mpo += 1
    if rotb <= 10: mpo += 1
    cns_mpo_ok = mpo >= 4.0
    bbb = (0 <= logp <= 5) and (150 <= mw <= 500) and (tpsa < 90) and (hbd <= 3)
    return {
        "smiles": smi, "mw": round(mw, 1), "tpsa": round(tpsa, 1),
        "logp": round(logp, 2), "hbd": hbd, "hba": hba, "rotb": rotb,
        "cns_mpo": round(mpo, 1), "cns_mpo_ok": cns_mpo_ok, "bbb_ok": bbb,
    }


def main():
    excl = set()
    if os.path.exists(EXCL_PATH):
        with open(EXCL_PATH, encoding="utf-8") as f:
            excl = {l.strip() for l in f if l.strip()}

    smiles_idx = {}
    cl = os.path.join(BASE, "rescoring_datos", "candidatos_limpios.csv")
    if os.path.exists(cl):
        with open(cl, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                smiles_idx.setdefault(r["ligand"], r["smiles"])
    smi_path = os.path.join(ANALYSIS, "full_library_solo.smi")
    if os.path.exists(smi_path):
        with open(smi_path, encoding="utf-8", errors="replace") as f:
            for ln in f:
                ln = ln.strip()
                if "\t" in ln:
                    smi, nombre = ln.rsplit("\t", 1)
                    smiles_idx.setdefault(nombre, smi)

    rows = []
    with open(AFI, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append((r["target"], r["ligand"], float(r["afinidad"])))

    por_target = defaultdict(int)
    for t, _, _ in rows:
        por_target[t] += 1

    top5 = []
    cortes = {}
    for t in TARGETS:
        sub = sorted([r for r in rows if r[0] == t], key=lambda r: r[2])
        n = max(1, int(round(len(sub) * 0.05)))
        cortes[t] = {"n_top5": n, "corte_afinidad": sub[n - 1][2],
                     "mejor": sub[0][2], "total": len(sub)}
        top5.extend(sub[:n])

    out = os.path.join(ANALYSIS, "ranking_top5_consolidado.csv")
    cols = ["target", "ligand", "afinidad", "smiles", "mw", "tpsa",
            "logp", "hbd", "hba", "rotb", "cns_mpo", "cns_mpo_ok", "bbb_ok"]
    n_con_smiles = 0
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for t, lig, aff in sorted(top5, key=lambda r: (r[0], r[2])):
            smi = smiles_idx.get(lig, "")
            pr = props(smi) if smi else {}
            if smi:
                n_con_smiles += 1
            w.writerow({"target": t, "ligand": lig, "afinidad": round(aff, 3),
                        "smiles": smi, **pr})

    resumen = {
        "generado": "2026-09-09",
        "outputs_por_target": dict(por_target),
        "excluidos_documentados": len(excl),
        "top5_por_proteina": cortes,
        "top5_con_smiles": n_con_smiles,
        "top5_total": len(top5),
        "csv": out,
    }
    with open(os.path.join(ANALYSIS, "ranking_consolidado_resumen.json"), "w", encoding="utf-8") as f:
        json.dump(resumen, f, indent=2, ensure_ascii=False)
    print(json.dumps(resumen, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()