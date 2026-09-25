#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara el paquete para la tanda 2 del rescoring MM-GBSA en Oracle:
   - candidatos_231.csv con ruta de pose por fila
   - poses/  (PDBQT acoplados, buscando en TODAS las tandas)
   - receptores/ (los PDBQT usados en docking)
   Comprimido en rescoring_pkg2.tar.gz
"""
import os
import shutil
import tarfile
import glob

import pandas as pd

ROOT = r"C:/Users/Fredy/masive-als"
GD = os.path.join(ROOT, "gpu_dock")
STAGE = os.path.join(ROOT, "analysis", "_rescoring_stage2")

def localizar_pose(lig, target):
    """Busca {lig}_{target}_out.pdbqt o {lig}_out.pdbqt en cualquier tanda."""
    patrones = [
        os.path.join(GD, "tanda_*", f"results_{target}", f"{lig}_{target}_out.pdbqt"),
        os.path.join(GD, "tanda_*", f"results_{target}", f"{lig}_out.pdbqt"),
        os.path.join(GD, "tanda_*", f"{lig}_{target}_out.pdbqt"),
        os.path.join(GD, "tanda_*", f"{lig}_out.pdbqt"),
    ]
    for p in patrones:
        hits = glob.glob(p)
        if hits:
            return hits[0]
    return None

# 1. candidatos nuevos (187 que pasan CNS con los datos actualizados)
df = pd.read_csv(os.path.join(ROOT, "analysis", "_nuevos_cns_total3.csv"))
print("candidatos:", len(df), "| targets:", df["target"].value_counts().to_dict())

# 2. crear stage sin borrar poses ya generadas (reacopladas en CPU)
os.makedirs(os.path.join(STAGE, "poses"), exist_ok=True)
os.makedirs(os.path.join(STAGE, "receptores"), exist_ok=True)

# 3. copiar poses y rellenar ruta relativa (respetando las ya presentes en el stage)
pose_rel = []
faltan = []
for _, r in df.iterrows():
    dst = os.path.join(STAGE, "poses", f"{r['ligand']}_{r['target']}_out.pdbqt")
    if os.path.exists(dst):
        pose_rel.append(f"poses/{r['ligand']}_{r['target']}_out.pdbqt")
        continue
    src = localizar_pose(r["ligand"], r["target"])
    if src is None:
        faltan.append(f"{r['ligand']} {r['target']}")
        pose_rel.append("")
        continue
    shutil.copy(src, dst)
    pose_rel.append(f"poses/{r['ligand']}_{r['target']}_out.pdbqt")

df["pose_pdbqt"] = pose_rel
print("poses encontradas:", (df["pose_pdbqt"] != "").sum(), "| faltan:", len(faltan))
if faltan:
    print("FALTAN (primeros 20):")
    for f in faltan[:20]:
        print("  ", f)

# 4. receptores (los mismos PDBQT usados en docking)
for tgt in ["SOD1", "TDP43", "FUS"]:
    src = os.path.join(GD, f"{tgt}.pdbqt")
    if not os.path.exists(src):
        raise FileNotFoundError(f"falta receptor: {src}")
    shutil.copy(src, os.path.join(STAGE, "receptores", f"{tgt}.pdbqt"))

# 5. CSV dentro del stage (solo filas con pose)
df = df[df["pose_pdbqt"] != ""].copy()
csv_path = os.path.join(STAGE, "candidatos_231.csv")
df.to_csv(csv_path, index=False)
print("CSV:", csv_path, "| filas con pose:", len(df))

# 6. tarball
tar_path = os.path.join(ROOT, "analysis", "rescoring_pkg2.tar.gz")
with tarfile.open(tar_path, "w:gz") as tf:
    tf.add(STAGE, arcname="rescoring_pkg")
print("tar:", tar_path, os.path.getsize(tar_path), "bytes")
print("receptores:", os.listdir(os.path.join(STAGE, "receptores")))
