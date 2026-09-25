# -*- coding: utf-8 -*-
"""
Ranking consolidado MASIVE-ALS desde los resultados del GPU (disco).

- Lee results_<target>/*_out.pdbqt de las 3 proteinas vigentes
  (TDP43_v2, SOD1, FUS). TDP43 viejo (4IUF) queda descartado cientificamente.
- Documenta los ligandos excluidos (imposibles).
- Propiedades moleculares + filtro CNS SOLO para el top-5% por proteina
  (evita RDKit sobre los ~107k ligandos).
- Escribe analysis/ranking_consolidado.csv (afinidades) y
  ranking_top5_consolidado.csv (top-5% con props y filtros CNS).
"""
import glob
import os
import csv
import json
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

BASE = r"C:\Users\Fredy\masive-als"
RESDIR = os.path.join(BASE, "gpu_dock", "resultados_libreria")
ANALYSIS = os.path.join(BASE, "analysis")
EXCL_PATH = os.path.join(BASE, "gpu_dock", "ligandos_excluidos.txt")

TARGETS = ["TDP43_v2", "SOD1", "FUS"]
ETIQUETA = {
    "TDP43_v2": "TDP-43 (4BS2 RRM1-RRM2)",
    "SOD1": "SOD1 (Trp32 anti-agregacion)",
    "FUS": "FUS (6G99)",
}


def leer_afinidad(pdbqt_path):
    try:
        with open(pdbqt_path, encoding="utf-8", errors="replace") as f:
            for l in f:
                if l.startswith("REMARK VINA RESULT"):
                    return float(l.split()[3])
    except Exception:
        pass
    return None


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
    print("excluidos:", len(excl))

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
    print("smiles indexados:", len(smiles_idx))

    rows = []
    por_target = defaultdict(int)
    for t in TARGETS:
        outdir = os.path.join(RESDIR, "results_" + t)
        n = 0
        for p in glob.glob(os.path.join(outdir, "*_out.pdbqt")):
            lig = os.path.basename(p).replace("_out.pdbqt", "")
            aff = leer_afinidad(p)
            if aff is None:
                continue
            n += 1
            rows.append((t, lig, aff))
        por_target[t] = n
        print("  %s: %d" % (t, n))

    # CSV completo (afinidades + smiles cuando exista)
    out_csv = os.path.join(ANALYSIS, "ranking_consolidado.csv")
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["target", "ligand", "afinidad", "smiles", "excluido"])
        for t, lig, aff in sorted(rows, key=lambda r: (r[0], r[2])):
            w.writerow([t, lig, round(aff, 3), smiles_idx.get(lig, ""), "si" if lig in excl else ""])
    print("csv completo:", out_csv, "filas:", len(rows))

    # top-5% por proteina
    top5 = []
    for t in TARGETS:
        sub = sorted([r for r in rows if r[0] == t], key=lambda r: r[2])
        n = max(1, int(round(len(sub) * 0.05)))
        top5.extend(sub[:n])

    out_top = os.path.join(ANALYSIS, "ranking_top5_consolidado.csv")
    with open(out_top, "w", newline="", encoding="utf-8") as f:
        cols = ["target", "ligand", "afinidad", "smiles", "mw", "tpsa",
                "logp", "hbd", "hba", "rotb", "cns_mpo", "cns_mpo_ok", "bbb_ok"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for t, lig, aff in sorted(top5, key=lambda r: (r[0], r[2])):
            smi = smiles_idx.get(lig, "")
            pr = props(smi) if smi else {}
            w.writerow({"target": t, "ligand": lig, "afinidad": round(aff, 3),
                        "smiles": smi, **pr})
    print("csv top5:", out_top, "filas:", len(top5))

    resumen = {
        "generado": "2026-09-09",
        "outputs_por_target": dict(por_target),
        "excluidos_documentados": len(excl),
        "top5_por_proteina": {
            t: {
                "n": sum(1 for r in top5 if r[0] == t),
                "mejor": min((r[2] for r in rows if r[0] == t), default=None),
                "corte_5pct": sorted([r[2] for r in rows if r[0] == t])[max(0, int(round(len([r for r in rows if r[0] == t]) * 0.05)) - 1)],
            }
            for t in TARGETS
        },
    }
    with open(os.path.join(ANALYSIS, "ranking_consolidado_resumen.json"), "w", encoding="utf-8") as f:
        json.dump(resumen, f, indent=2, ensure_ascii=False)
    print(json.dumps(resumen, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()