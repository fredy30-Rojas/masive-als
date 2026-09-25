# -*- coding: utf-8 -*-
"""
Lanzador del rescoring MM-GBSA (OpenFF) sobre la LISTA CORTA (2.973).

- Entrada: analysis/lista_corta_candidatos.csv (target, ligand, afinidad, smiles)
- Poses:   gpu_dock/resultados_libreria/results_<target>/<ligand>_out.pdbqt
- Receptores: gpu_dock/<TDP43_v2|SOD1|FUS>.pdbqt
- Salida:  analysis/rescoring_local/rescoring_lista_corta.csv (reanudable)
- Usa:     rescoring_env (conda-forge, python 3.11, openmm+openff)
- Reanuda: salta pares (target,ligand) ya presentes con mmgbsa_dG valido.
- Paralelismo: --workers N (subprocesos), default 4.

Uso (Windows, desde bash):
  export PATH="/c/Users/Fredy/masive-als/rescoring_venv/micromamba_bin/Library/bin:$PATH"
  export MAMBA_ROOT_PREFIX="/c/Users/Fredy/masive-als/rescoring_env_root"
  micromamba run -p C:/Users/Fredy/masive-als/rescoring_env python \
      analysis/rescoring_local/runner_mmgbsa_local.py --workers 4
"""
import argparse
import csv
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

BASE = r"C:\Users\Fredy\masive-als"
ANALYSIS = os.path.join(BASE, "analysis")
OUTDIR = os.path.join(ANALYSIS, "rescoring_local")
CAND = os.path.join(ANALYSIS, "lista_corta_candidatos.csv")
OUT = os.path.join(OUTDIR, "rescoring_lista_corta.csv")
LOG = os.path.join(OUTDIR, "runner_lista_corta.log")
RES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
RECEPTORES = {
    "TDP43_v2": os.path.join(BASE, "gpu_dock", "TDP43_v2.pdbqt"),
    "SOD1": os.path.join(BASE, "gpu_dock", "SOD1.pdbqt"),
    "FUS": os.path.join(BASE, "gpu_dock", "FUS.pdbqt"),
}
COLS = ["target", "ligand", "vina_affinity", "smiles", "mmgbsa_dG",
        "e_complex", "e_receptor", "e_ligand", "error"]
TIMEOUT = 900  # 15 min por candidato


def log(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def load_candidates():
    rows = []
    with open(CAND, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append({"target": r["target"], "ligand": r["ligand"],
                         "affinity": r["afinidad"], "smiles": r["smiles"]})
    return rows


def load_done():
    done = {}
    if os.path.exists(OUT):
        with open(OUT, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                done[(r["target"], r["ligand"])] = r
    return done


def calcular_fila(r):
    tgt, lig, smi = r["target"], r["ligand"], r["smiles"]
    pose = os.path.join(RES, "results_" + tgt, lig + "_out.pdbqt")
    receptor = RECEPTORES.get(tgt)
    if not receptor or not os.path.exists(receptor):
        return {"target": tgt, "ligand": lig, "mmgbsa_dG": None,
                "error": "sin receptor"}
    if not os.path.exists(pose):
        return {"target": tgt, "ligand": lig, "mmgbsa_dG": None,
                "error": "sin pose"}
    cmd = [sys.executable, os.path.join(OUTDIR, "mmgbsa_openff.py"),
           "--pose", pose, "--receptor", receptor, "--smiles", smi, "--out", "json"]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
        if p.returncode == 0 and p.stdout.strip():
            try:
                row = json.loads(p.stdout.strip().splitlines()[-1])
            except json.JSONDecodeError:
                row = None
            if row is not None and row.get("mmgbsa_dG") is not None:
                return {"target": tgt, "ligand": lig,
                        "vina_affinity": r.get("affinity"), "smiles": smi,
                        "mmgbsa_dG": row["mmgbsa_dG"],
                        "e_complex": row.get("e_complex"),
                        "e_receptor": row.get("e_receptor"),
                        "e_ligand": row.get("e_ligand"), "error": ""}
        err = (p.stderr or p.stdout or "salida vacia")[-250:]
        return {"target": tgt, "ligand": lig, "vina_affinity": r.get("affinity"),
                "smiles": smi, "mmgbsa_dG": None, "error": err.strip()}
    except subprocess.TimeoutExpired:
        return {"target": tgt, "ligand": lig, "vina_affinity": r.get("affinity"),
                "smiles": smi, "mmgbsa_dG": None, "error": "TIMEOUT>%ds" % TIMEOUT}
    except Exception as e:
        return {"target": tgt, "ligand": lig, "vina_affinity": r.get("affinity"),
                "smiles": smi, "mmgbsa_dG": None, "error": str(e)[:200]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max", type=int, default=None, help="max filas (test)")
    args = ap.parse_args()
    os.makedirs(OUTDIR, exist_ok=True)
    cands = load_candidates()
    done = load_done()
    log("candidatos: %d | ya hechos: %d" % (len(cands), len(done)))

    pend = [r for r in cands if (r["target"], r["ligand"]) not in done]
    if args.max:
        pend = pend[:args.max]
    # Orden: primero los targets que SÍ llegan a la GPU (SOD1, TDP43_v2);
    # FUS al final porque falla en la preparación (receptor no estándar)
    # y nunca llega a OpenMM, dejando la GPU ociosa si va primero.
    _PRIO = {"SOD1": 0, "TDP43_v2": 1, "FUS": 2}
    pend.sort(key=lambda r: _PRIO.get(r["target"], 5))
    from collections import Counter
    log("pendientes: %d (orden: %s)" % (
        len(pend),
        Counter(r["target"] for r in pend).most_common()))

    if not pend:
        log("Nada que hacer. DONE.")
        return

    new_file = not os.path.exists(OUT) or os.path.getsize(OUT) == 0
    outf = open(OUT, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(outf, fieldnames=COLS, extrasaction="ignore")
    if new_file:
        w.writeheader()

    n_ok = n_err = 0
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        fut = {ex.submit(calcular_fila, r): r for r in pend}
        for i, f in enumerate(fut):
            r = fut[f]
            row = f.result()
            w.writerow(row)
            outf.flush()
            done[(r["target"], r["ligand"])] = row
            if row.get("mmgbsa_dG") is not None:
                n_ok += 1
                log("[%d/%d] %s %s dG=%.2f"
                    % (i + 1, len(pend), r["target"], r["ligand"], row["mmgbsa_dG"]))
            else:
                n_err += 1
                log("[%d/%d] %s %s ERROR %s"
                    % (i + 1, len(pend), r["target"], r["ligand"],
                       str(row.get("error", ""))[:100]))
    outf.close()
    log("FIN: ok=%d err=%d | CSV: %s" % (n_ok, n_err, OUT))


if __name__ == "__main__":
    main()