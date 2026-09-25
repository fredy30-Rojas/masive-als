#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿La conformacion del cristal es alcanzable por el acoplamiento? (20 sep 2026).

Antes de culpar al acoplamiento de un redocking fallido hay que descartar una
causa tonta: que el ligando de entrada este en una conformacion tan distinta a
la del cristal que ninguna colocacion pueda acercarse. Como Vina mueve enlaces
rotables pero NO cambia la conformacion de los anillos, el limite real lo ponen
los anillos.

Se generan 300 conformeros con ETKDG, se superpone cada uno al ligando del
cristal y se toma el minimo. Si el minimo es pequeno (< 0,5 A), la conformacion
cristalina es representable y un fallo del redocking es del acoplamiento, no de
la preparacion. Si es grande, la preparacion impone un suelo.

Uso: python conformacion_cristal.py
"""
import os
import sys

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from redock_trp32 import BASE, SISTEMAS, aplanar

RDLogger.DisableLog("rdApp.*")


def kabsch_rmsd(P, Q):
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return np.sqrt(((((P - Pc) @ R.T) + Qc - Q) ** 2).sum(1).mean())


def minimo_alineado(mol_confs, cristal, n_max_matches=2000):
    m = aplanar(Chem.RemoveHs(Chem.Mol(mol_confs)))
    c = aplanar(Chem.RemoveHs(cristal))
    cc = c.GetConformer().GetPositions()
    matches = c.GetSubstructMatches(m, uniquify=False, maxMatches=n_max_matches)
    if not matches:
        return None
    mejor = float("inf")
    for cid in range(mol_confs.GetNumConformers()):
        cp = m.GetConformer(cid).GetPositions()
        for match in matches:
            v = kabsch_rmsd(cp[np.array(match)], cc)
            if v < mejor:
                mejor = v
    return mejor


def main():
    print("Diferencia conformacional entre el ligando y su pose cristalografica")
    print("(minimo sobre 300 conformeros generados; mide si el espacio de busqueda")
    print(" puede representar la conformacion del cristal)\n")
    for entrada, (codigo, _, nombre) in SISTEMAS.items():
        sdf = os.path.join(BASE, "ligands", "%s_ideal.sdf" % codigo)
        pdb = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
        if not (os.path.exists(sdf) and os.path.exists(pdb)):
            print("  %s (%s): faltan datos" % (entrada, nombre))
            continue
        mol = Chem.AddHs(Chem.MolFromMolFile(sdf, removeHs=False))
        AllChem.EmbedMultipleConfs(mol, numConfs=300, randomSeed=42)
        cristal = Chem.MolFromPDBBlock(open(pdb).read(), removeHs=False,
                                       proximityBonding=True)
        val = minimo_alineado(mol, cristal)
        # y la misma medida con un solo conformero (lo que se le dio a Vina)
        uno = Chem.AddHs(Chem.MolFromMolFile(sdf, removeHs=False))
        AllChem.EmbedMolecule(uno, AllChem.ETKDGv3())
        val_uno = minimo_alineado(uno, cristal)
        print("  %s (%s): 300 conformeros %.2f A | un solo conformero %.2f A"
              % (entrada, nombre, val if val else float("nan"),
                 val_uno if val_uno else float("nan")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
