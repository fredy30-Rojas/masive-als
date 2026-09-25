# -*- coding: utf-8 -*-
"""Comprueba que las poses reacopladas pasan los mismos controles que el hijo.

Hace exactamente lo que mmgbsa_openff_gb.prepare_ligand_pose() antes de calcular:
  1. lee la pose con M.pose_to_rdmol(),
  2. compara el numero de atomos con el del SMILES,
  3. exige que el MCS cubra todos los atomos.

Se corre con el entorno del rescoring (necesita openff/openmm):
    /root/rescoring_env/bin/python verificar_poses_macro.py
"""
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M                                     # noqa: E402

from rdkit import Chem, RDLogger                              # noqa: E402
from rdkit.Chem import rdFMCS                                 # noqa: E402
RDLogger.DisableLog("rdApp.*")

AQUI = os.path.dirname(os.path.abspath(__file__))
LISTA = os.path.join(AQUI, "lista_focalizada_rescoring.csv")
# La carpeta se escribe de una forma desde Windows y de otra desde WSL.
POSES = r"C:\Users\Fredy\masive-als\gpu_dock\rescoring_caja_arn\results_TDP43_v2"
if not os.path.isdir(POSES):
    POSES = "/mnt/c/Users/Fredy/masive-als/gpu_dock/rescoring_caja_arn/results_TDP43_v2"
FIJOS = ["CHEMBL152275", "CHEMBL436025", "CHEMBL154281", "CHEMBL282984",
         "CHEMBL403578", "CHEMBL357525", "CHEMBL28028", "CHEMBL449782",
         "CHEMBL5810661", "CHEMBL5870313", "CHEMBL89479", "CHEMBL275006",
         "CHEMBL276462", "CHEMBL312862"]


def main():
    smi = {r["ligand"]: r["smiles"]
           for r in csv.DictReader(open(LISTA, encoding="utf-8"))}
    todo = True
    for lig in FIJOS:
        pose = os.path.join(POSES, lig + "_out.pdbqt")
        docked = M.pose_to_rdmol(pose)
        rd = Chem.MolFromSmiles(smi[lig])
        mcs = rdFMCS.FindMCS([rd, docked],
                             bondCompare=rdFMCS.BondCompare.CompareAny,
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             timeout=30)
        ok = (rd.GetNumAtoms() == docked.GetNumAtoms()
              and mcs.numAtoms == rd.GetNumAtoms())
        todo = todo and ok
        print("%-14s SMILES %2d pose %2d | MCS %2d/%2d -> %s"
              % (lig, rd.GetNumAtoms(), docked.GetNumAtoms(), mcs.numAtoms,
                 rd.GetNumAtoms(), "OK" if ok else "REVISAR"))
    print("TODO OK" if todo else "HAY ALGO MAL")
    return 0 if todo else 1


if __name__ == "__main__":
    sys.exit(main())
