#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿El embedding con las coordenadas de la pose sale a la primera?

Para cada ligando dice el codigo de retorno y cuanto se aparta de la pose, que
es lo que decide si el ligando acaba dentro del bolsillo o flotando fuera.
"""
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/ubuntu/mmgbsa")

import pose_pdbqt  # noqa: E402
from rdkit import Chem  # noqa: E402
from rdkit.Chem import AllChem  # noqa: E402

LIGS = ["ACT_adrenalina", "ACT_isoproterenol", "DEC_CHEMBL1414576",
        "DECM_CHEMBL978", "DECM_CHEMBL39736", "DECM_CHEMBL4435214",
        "DECH_CHEMBL3309988"]


def rmsd(rdmol_h, coordmap):
    conf = rdmol_h.GetConformer()
    d = []
    for j, p in coordmap.items():
        q = conf.GetAtomPosition(j)
        d.append((q.x - p.x) ** 2 + (q.y - p.y) ** 2 + (q.z - p.z) ** 2)
    return float(np.sqrt(np.mean(d)))


def main():
    filas = pd.read_csv("estratos/candidatos.csv")
    for lig in LIGS:
        f = filas[filas.ligand == lig]
        if f.empty:
            print(f"{lig}: no está"); continue
        f = f.iloc[0]
        rdmol = Chem.MolFromSmiles(f["smiles"])
        pose = pose_pdbqt.leer_pose(f["pose_pdbqt"])
        coordmap = pose_pdbqt.coordmap_de_pose(pose, rdmol)
        rdmol_h = Chem.AddHs(rdmol)
        pr = AllChem.ETKDGv3()
        pr.randomSeed = 42
        pr.useRandomCoords = True
        pr.maxIterations = 500
        pr.SetCoordMap(coordmap)
        ret = AllChem.EmbedMolecule(rdmol_h, pr)
        extra = ""
        if ret == 0:
            extra = f"rmsd_a_la_pose={rmsd(rdmol_h, coordmap):.3f} A"
        print(f"{lig:22s} pesados={rdmol.GetNumAtoms():3d} ret={ret:3d}  {extra}")


if __name__ == "__main__":
    main()
