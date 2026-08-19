#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara el paquete para el rescoring MM-GBSA en Oracle:
   - candidatos_42.csv con ruta de pose por fila
   - poses/  (42 PDBQT acoplados)
   - receptores/ (SOD1.pdbqt, TDP43.pdbqt, FUS.pdbqt)
   Todo comprimido en rescoring_pkg.tar.gz
"""
import os
import shutil
import tarfile

import pandas as pd

ROOT = r"C:/Users/Fredy/masive-als"
GD = os.path.join(ROOT, "gpu_dock")
STAGE = os.path.join(ROOT, "analysis", "_rescoring_stage")

# 1. candidatos
df = pd.read_csv(os.path.join(ROOT, "analysis", "candidatos_filtrados.csv"))
print("candidatos:", len(df), "targets:", df["target"].value_counts().to_dict())

# 2. limpiar stage
if os.path.exists(STAGE):
    shutil.rmtree(STAGE)
os.makedirs(os.path.join(STAGE, "poses"))
os.makedirs(os.path.join(STAGE, "receptores"))

# 3. copiar poses y rellenar ruta relativa
pose_rel = []
for _, r in df.iterrows():
    src = os.path.join(GD, "tanda_z001", f"results_{r['target']}",
                       f"{r['ligand']}_out.pdbqt")
    if not os.path.exists(src):
        raise FileNotFoundError(f"falta pose: {src}")
    dst = os.path.join(STAGE, "poses", f"{r['ligand']}_{r['target']}_out.pdbqt")
    shutil.copy(src, dst)
    pose_rel.append(f"poses/{r['ligand']}_{r['target']}_out.pdbqt")

df["pose_pdbqt"] = pose_rel

# 4. receptores (los mismos PDBQT usados en docking)
for tgt in ["SOD1", "TDP43", "FUS"]:
    shutil.copy(os.path.join(GD, f"{tgt}.pdbqt"),
                os.path.join(STAGE, "receptores", f"{tgt}.pdbqt"))

# 5. CSV dentro del stage
csv_path = os.path.join(STAGE, "candidatos_42.csv")
df.to_csv(csv_path, index=False)
print("CSV:", csv_path)

# 6. tarball
tar_path = os.path.join(ROOT, "analysis", "rescoring_pkg.tar.gz")
with tarfile.open(tar_path, "w:gz") as tf:
    tf.add(STAGE, arcname="rescoring_pkg")
print("tar:", tar_path, os.path.getsize(tar_path), "bytes")
print("n poses:", len(df), "| receptores:", os.listdir(os.path.join(STAGE, "receptores")))
