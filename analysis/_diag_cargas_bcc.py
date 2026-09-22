#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""_diag_cargas_bcc.py — ¿llegan las cargas AM1-BCC al System de OpenMM?

Se ejecuta en /home/ubuntu/mmgbsa (env conda mmgbsa). Reproduce el preparador
de v7 sobre la fila 1 y comprueba la carga de un átomo del ligando en el
NonbondedForce del sistema final. Correrlo con el env activado.
"""
import pandas as pd

df = pd.read_csv("estratos/candidatos.csv")
r = df.iloc[1]
print("fila 1:", r["ligand"])

from openmm.app import ForceField
import rescoring_mmgbsa_bcc as bcc

ff = ForceField("amber14-all.xml", "implicit/obc2.xml")

from openff.toolkit import Molecule
from rdkit import Chem
from rdkit.Chem import AllChem
import pose_pdbqt

rdmol = Chem.MolFromSmiles(r["smiles"])
pose = pose_pdbqt.leer_pose(r["pose_pdbqt"])
cm = pose_pdbqt.coordmap_de_pose(pose, rdmol)
rh = Chem.AddHs(rdmol)
p = AllChem.ETKDGv3()
p.randomSeed = 42
p.useRandomCoords = True
p.maxIterations = 500
p.SetCoordMap(cm)
AllChem.EmbedMolecule(rh, p)
off = Molecule.from_rdkit(rh, allow_undefined_stereo=True)
from openff.toolkit.utils.toolkits import AmberToolsToolkitWrapper
off.assign_partial_charges(partial_charge_method="am1bcc",
                           toolkit_registry=AmberToolsToolkitWrapper())
q = [float(x.m) for x in off.partial_charges]
print("cargas AM1-BCC: suma=%.6f max_abs=%.4f" % (sum(q), max(abs(x) for x in q)))

# el flujo REAL de v7 (cargas_bcc=True dentro de prepare_ligand)
lig_top, lig_pos = bcc.prepare_ligand(r["smiles"], r["pose_pdbqt"], ff,
                                      cargas_bcc=True)
from openmm.app import Modeller
mod = Modeller(lig_top, lig_pos)
s = bcc.build_system(mod.topology, ff)
from openmm import unit
nbf = [f for f in s.getForces() if f.__class__.__name__ == "NonbondedForce"][0]
cargas = [nbf.getParticleParameters(i)[0].value_in_unit(unit.elementary_charge)
          for i in range(nbf.getNumParticles())]
print("carga max_abs en el System:", round(max(abs(c) for c in cargas), 4))
print("num particulas:", nbf.getNumParticles())
