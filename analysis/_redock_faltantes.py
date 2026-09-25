#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Re-acopla los ligandos faltantes para completar las poses de los 187 candidatos.
Busca el ligando PDBQT en compounds/*_pdbqt y corre Vina-GPU contra la diana
correspondiente, guardando la pose como {lig}_{target}_out.pdbqt en
analysis/_rescoring_stage2/poses/ (mismo formato que el resto del paquete).
"""
import os
import shutil
import subprocess
import sys
import tempfile
import time

BASE = r"C:/Users/Fredy/masive-als"
GD = os.path.join(BASE, "gpu_dock")
STAGE_POSES = os.path.join(BASE, "analysis", "_rescoring_stage2", "poses")
EXE = os.path.join(GD, "Vina-GPU-2.1-win.exe")

RECEPTORES = {
    "TDP43": {"receptor": os.path.join(GD, "TDP43.pdbqt"), "centro": (28.3, 43.7, 52.5), "size": 25},
    "SOD1":  {"receptor": os.path.join(GD, "SOD1.pdbqt"),  "centro": (46.5, 80.0, 73.3), "size": 22},
    "FUS":   {"receptor": os.path.join(GD, "FUS.pdbqt"),   "centro": (-14.5, 15.1, -7.8), "size": 25},
}

CANDIDATOS = os.path.join(BASE, "analysis", "_nuevos_cns_total3.csv")
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
    print(f"Faltan {len(faltan)} poses por generar")
    if not faltan:
        print("Nada que hacer.")
        return

    # agrupar por diana
    por_diana = {}
    for lig, tgt in faltan:
        por_diana.setdefault(tgt, []).append(lig)

    for tgt, ligs in por_diana.items():
        info = RECEPTORES[tgt]
        tmp = tempfile.mkdtemp(prefix=f"redock_{tgt}_")
        tmp_lig = os.path.join(tmp, "ligands")
        tmp_out = os.path.join(tmp, "out")
        os.makedirs(tmp_lig)
        os.makedirs(tmp_out)
        for lig in ligs:
            src = localizar_ligando(lig)
            if src is None:
                print(f"  !! {lig}: ligando PDBQT no encontrado en compounds/")
                continue
            shutil.copy(src, os.path.join(tmp_lig, lig + ".pdbqt"))

        cfg = os.path.join(tmp, "config.txt")
        with open(cfg, "w") as f:
            f.write("receptor = %s\n" % info["receptor"].replace("\\", "/"))
            f.write("ligand_directory = %s\n" % tmp_lig.replace("\\", "/"))
            f.write("output_directory = %s\n" % tmp_out.replace("\\", "/"))
            f.write("center_x = %s\n" % info["centro"][0])
            f.write("center_y = %s\n" % info["centro"][1])
            f.write("center_z = %s\n" % info["centro"][2])
            f.write("size_x = %s\n" % info["size"])
            f.write("size_y = %s\n" % info["size"])
            f.write("size_z = %s\n" % info["size"])
            f.write("num_modes = 1\n")

        print(f"Acoplando {len(ligs)} ligandos contra {tgt}...")
        with open(os.path.join(tmp, "vina.log"), "w") as lf:
            subprocess.run([EXE, "--config", cfg], cwd=GD,
                           stdout=lf, stderr=subprocess.STDOUT, timeout=3600)

        # copiar resultados con el nombre del paquete
        n = 0
        for p in os.listdir(tmp_out):
            base = p.replace("_out.pdbqt", "").replace(".pdbqt", "")
            dst = os.path.join(STAGE_POSES, f"{base}_{tgt}_out.pdbqt")
            shutil.copy(os.path.join(tmp_out, p), dst)
            n += 1
        print(f"  {tgt}: {n} poses generadas")
        shutil.rmtree(tmp, ignore_errors=True)

    # recontar
    df2 = pd.read_csv(CANDIDATOS)
    con_pose = sum(1 for _, r in df2.iterrows()
                   if os.path.exists(os.path.join(STAGE_POSES, f"{r['ligand']}_{r['target']}_out.pdbqt")))
    print(f"TOTAL con pose: {con_pose}/{len(df2)}")

if __name__ == "__main__":
    main()
