# -*- coding: utf-8 -*-
"""Runner v5 del rescoring MM-GBSA (MASIVE-ALS).

Igual que runner_rescoring.py pero:
- Usa receptores_canon.json (PDB canónicos reordenados en el frame del docking).
- Escribe a rescoring_mmgbsa_v5.csv.
- Reanuda: salta pares (ligand,target) ya calculados con error vacío.
"""
import csv
import json
import os
import subprocess
import sys
import time

CAND = "rescoring_pkg/candidatos_42.csv"
REC = "receptores_canon.json"
OUT = "rescoring_mmgbsa_v5.csv"
LOG = "runner_v5.log"
TIMEOUT = 600  # 10 min por candidato

COLS = ["ligand", "target", "vina_affinity", "mmgbsa_dG",
        "e_complex", "e_receptor", "e_ligand", "error"]


def log(msg):
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")


def load_candidates():
    return list(csv.DictReader(open(CAND)))


def load_done():
    done = {}
    if os.path.exists(OUT):
        with open(OUT) as f:
            for r in csv.DictReader(f):
                done[(r["ligand"], r["target"])] = r
    return done


def main():
    cands = load_candidates()
    done = load_done()
    log(f"candidatos: {len(cands)} | ya hechos: {len(done)}")

    new_file = not os.path.exists(OUT) or os.path.getsize(OUT) == 0
    outf = open(OUT, "a", newline="")
    w = csv.writer(outf)
    if new_file:
        w.writerow(COLS)

    pend = [r for r in cands if (r["ligand"], r["target"]) not in done]
    log(f"pendientes: {len(pend)}")

    for i, r in enumerate(pend):
        key = (r["ligand"], r["target"])
        cmd = [sys.executable, "rescoring_mmgbsa_robusto.py",
               "--candidates", CAND, "--receptores", REC, "--row",
               str(cands.index(r))]
        log(f"[{i+1}/{len(pend)}] {r['ligand']} {r['target']}")
        try:
            p = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=TIMEOUT)
            row = None
            if p.returncode == 0 and p.stdout.strip():
                lines = [l for l in p.stdout.splitlines() if l.strip()]
                if lines:
                    try:
                        row = json.loads(lines[-1])
                    except json.JSONDecodeError:
                        row = None
            if row is None:
                err = (p.stderr or p.stdout or "salida vacía")[-300:]
                row = {"ligand": r["ligand"], "target": r["target"],
                       "vina_affinity": r.get("affinity"),
                       "mmgbsa_dG": None, "error": err}
        except subprocess.TimeoutExpired:
            row = {"ligand": r["ligand"], "target": r["target"],
                   "vina_affinity": r.get("affinity"), "mmgbsa_dG": None,
                   "error": f"TIMEOUT>{TIMEOUT}s"}
            log("   TIMEOUT")
        except Exception as e:
            row = {"ligand": r["ligand"], "target": r["target"],
                   "vina_affinity": r.get("affinity"), "mmgbsa_dG": None,
                   "error": str(e)[:300]}
            log(f"   excepción runner: {e}")

        w.writerow([row.get(c, "") for c in COLS])
        outf.flush()
        if row.get("mmgbsa_dG") is not None:
            log(f"   dG={row['mmgbsa_dG']:.2f}")
        else:
            log(f"   error: {str(row.get('error', ''))[:120]}")
        done[key] = row

    outf.close()
    log("FIN runner v5")


if __name__ == "__main__":
    main()
