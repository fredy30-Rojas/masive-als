#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Re-acopla en CPU (Vina 1.1.2) los ligandos faltantes para completar las
poses de los 231 candidatos. No toca la GPU (el cribado sigue en paralelo).
Guarda la pose como {lig}_{target}_out.pdbqt en _rescoring_stage2/poses/.
"""
import os
import shutil
import subprocess
import sys
import tempfile
import glob

BASE = r"C:/Users/Fredy/masive-als"
GD = os.path.join(BASE, "gpu_dock")
STAGE_POSES = os.path.join(BASE, "analysis", "_rescoring_stage2", "poses")
VINA = os.path.join(BASE, "tools", "vina.exe")
CANDIDATOS = os.path.join(BASE, "analysis", "_nuevos_cns_total3.csv")

RECEPTORES = {
    "TDP43": {"receptor": os.path.join(GD, "TDP43.pdbqt"), "centro": (28.3, 43.7, 52.5), "size": 25},
    "SOD1":  {"receptor": os.path.join(GD, "SOD1.pdbqt"),  "centro": (46.5, 80.0, 73.3), "size": 22},
    "FUS":   {"receptor": os.path.join(GD, "FUS.pdbqt"),   "centro": (-14.5, 15.1, -7.8), "size": 25},
}

LIB_DIRS = [os.path.join(BASE, "compounds", d) for d in
            ["batch2_pdbqt", "fda_full_pdbqt", "lotes_v3_pdbqt",
             "lotes_v4_pdbqt", "lotes_100_pdbqt", "next_batches_pdbqt"]]

import pandas as pd


def localizar_ligando(lig):
    for d in LIB_DIRS:
        p = os.path.join(d, lig + ".pdbqt")
        if os.path.exists(p):
            return p
    return None


def main():
    df = pd.read_csv(CANDIDATOS)
    faltan = []
    for _, r in df.iterrows():
        p = os.path.join(STAGE_POSES, f"{r['ligand']}_{r['target']}_out.pdbqt")
        if not os.path.exists(p):
            faltan.append((r["ligand"], r["target"]))
    print(f"Faltan {len(faltan)} poses (CPU Vina, sin tocar GPU)")
    if not faltan:
        print("Nada que hacer.")
        return

    ok, fail = 0, []
    for lig, tgt in faltan:
        dst = os.path.join(STAGE_POSES, f"{lig}_{tgt}_out.pdbqt")
        src_lig = localizar_ligando(lig)
        if src_lig is None:
            fail.append((lig, tgt, "sin ligando PDBQT"))
            continue
        tmp = tempfile.mkdtemp(prefix=f"vina_cpu_")
        lig_pdbqt = os.path.join(tmp, lig + ".pdbqt")
        out_pdbqt = os.path.join(tmp, lig + "_out.pdbqt")
        shutil.copy(src_lig, lig_pdbqt)
        info = RECEPTORES[tgt]
        cfg = os.path.join(tmp, "config.txt")
        with open(cfg, "w") as f:
            f.write("receptor = %s\n" % info["receptor"].replace("\\", "/"))
            f.write("ligand = %s\n" % lig_pdbqt.replace("\\", "/"))
            f.write("out = %s\n" % out_pdbqt.replace("\\", "/"))
            f.write("center_x = %s\n" % info["centro"][0])
            f.write("center_y = %s\n" % info["centro"][1])
            f.write("center_z = %s\n" % info["centro"][2])
            f.write("size_x = %s\n" % info["size"])
            f.write("size_y = %s\n" % info["size"])
            f.write("size_z = %s\n" % info["size"])
            f.write("num_modes = 1\n")
            f.write("cpu = 4\n")
        try:
            r = subprocess.run([VINA, "--config", cfg],
                               capture_output=True, text=True, timeout=600)
            if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 100:
                shutil.copy(out_pdbqt, dst)
                ok += 1
                print(f"  OK {lig} {tgt}")
            else:
                tail = (r.stderr or r.stdout or "")[-150:]
                fail.append((lig, tgt, tail.replace("\n", " ")))
                print(f"  FAIL {lig} {tgt}: {tail[:80]}")
        except subprocess.TimeoutExpired:
            fail.append((lig, tgt, "TIMEOUT"))
            print(f"  TIMEOUT {lig} {tgt}")
        except Exception as e:
            fail.append((lig, tgt, str(e)[:100]))
            print(f"  ERROR {lig} {tgt}: {e}")
        finally:
            shutil.rmtree(tmp, ignore_errors=True)

    print(f"\nOK: {ok} | fallos: {len(fail)}")
    for lig, tgt, err in fail:
        print(f"  {lig} {tgt}: {err[:120]}")

    # recontar
    n_total = len(df)
    n_con = sum(1 for _, r in df.iterrows()
                if os.path.exists(os.path.join(STAGE_POSES, f"{r['ligand']}_{r['target']}_out.pdbqt")))
    print(f"TOTAL con pose: {n_con}/{n_total}")


if __name__ == "__main__":
    main()
