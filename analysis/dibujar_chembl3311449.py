# -*- coding: utf-8 -*-
"""Dibuja CHEMBL3311449 (SOD1, dG MM-GBSA -15.79) en 2D y 3D con RDKit."""
import os
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem import rdMolTransforms
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

SMI = "O=C(Nc1cccc(Nc2nc(-c3ccccc3)c3cc[nH]c3n2)c1)N1CCCC1"
OUTDIR = r"C:/Users/Fredy/masive-als/paper/figures"
os.makedirs(OUTDIR, exist_ok=True)

mol = Chem.MolFromSmiles(SMI)
mol = Chem.AddHs(mol)

# ---- 2D (sin hidrógenos) ----
mol2d = Chem.RemoveHs(mol)
rdDepictor.Compute2DCoords(mol2d)
drawer = rdMolDraw2D.MolDraw2DCairo(900, 700)
opts = drawer.drawOptions()
opts.bondLineWidth = 2
opts.padding = 0.1
rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol2d)
drawer.FinishDrawing()
png2d = os.path.join(OUTDIR, "chembl3311449_2d.png")
with open(png2d, "wb") as f:
    f.write(drawer.GetDrawingText())
print("2D ->", png2d)

# ---- 3D (conformero ETKDG + MMFF, vista inclinada) ----
params = AllChem.ETKDGv3()
params.randomSeed = 0xABCDEF
status = AllChem.EmbedMolecule(mol, params)
print("embed status =", status)
if status == 0:
    AllChem.MMFFOptimizeMolecule(mol, maxIters=1000)

    def rot_mat(axis, ang):
        c, s = np.cos(ang), np.sin(ang)
        if axis == 'x':
            return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
        if axis == 'y':
            return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])

    R = rot_mat('y', 0.6) @ rot_mat('x', -0.35)
    rdMolTransforms.TransformConformer(mol.GetConformer(), R)

    mol_noh = Chem.RemoveHs(mol)
    d3 = rdMolDraw2D.MolDraw2DCairo(900, 700)
    o3 = d3.drawOptions()
    o3.bondLineWidth = 2
    o3.padding = 0.12
    d3.DrawMolecule(mol_noh, confId=0)  # usa las coords 3D rotadas
    d3.FinishDrawing()
    png3d = os.path.join(OUTDIR, "chembl3311449_3d.png")
    with open(png3d, "wb") as f:
        f.write(d3.GetDrawingText())
    print("3D ->", png3d)
else:
    print("ERROR: fallo al embeber conformero 3D")

# ---- Info ----
m_light = Chem.MolFromSmiles(SMI)
print("MW =", round(Descriptors.MolWt(m_light), 1))
print("TPSA =", round(rdMolDescriptors.CalcTPSA(m_light), 1))
print("logP =", round(Crippen.MolLogP(m_light), 2))
print("HBD =", rdMolDescriptors.CalcNumHBD(m_light))
print("HBA =", rdMolDescriptors.CalcNumHBA(m_light))
print("rotB =", rdMolDescriptors.CalcNumRotatableBonds(m_light))
