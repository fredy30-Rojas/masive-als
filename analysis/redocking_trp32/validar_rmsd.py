#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Comprueba el medidor de RMSD contra una implementacion independiente.

El 20 de septiembre se encontro un fallo en `redock_trp32.rmsd_en_sitio`: los
atomos se emparejaban al reves, y eso inflaba el RMSD (isoproterenol de 4A7T de
1,88 a 0,46 A). Este script no arregla nada; comprueba que la version corregida
coincide con `rdkit.Chem.rdMolAlign.CalcRMS`, que es la implementacion de RDKit
del RMSD en el sitio CON correccion de simetria.

Como CalcRMS necesita que las dos moleculas tengan el mismo orden de atomos, la
pose se renumera con el emparejamiento encontrado antes de llamarla. Para
comparar tambien se calcula la version vieja (emparejamiento invertido).

Uso: python validar_rmsd.py
"""
import glob
import os

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolAlign

import redock_trp32 as R

RDLogger.DisableLog("rdApp.*")
BASE = R.BASE


def version_vieja(pdb_lig, pose, conf_id=0):
    """La version con el error: cp[i] contra cc[match[i]]."""
    cry = Chem.MolFromPDBBlock(open(pdb_lig).read(), removeHs=False,
                               proximityBonding=True)
    p = R.aplanar(Chem.RemoveHs(Chem.Mol(pose)))
    c = R.aplanar(Chem.RemoveHs(cry))
    cp = p.GetConformer(conf_id).GetPositions()
    cc = c.GetConformer().GetPositions()
    mejor = None
    for match in p.GetSubstructMatches(c, uniquify=False, maxMatches=500,
                                       useChirality=False):
        d2 = sum((cp[i][k] - cc[match[i]][k]) ** 2
                 for i in range(len(match)) for k in range(3))
        v = (d2 / len(match)) ** 0.5
        mejor = v if mejor is None else min(mejor, v)
    return mejor


def rmsd_rdkit(pdb_lig, pose, conf_id=0):
    """CalcRMS de RDKit, renumerando la pose al orden del cristal."""
    cry = Chem.MolFromPDBBlock(open(pdb_lig).read(), removeHs=False,
                               proximityBonding=True)
    # Se quitan los hidrogenos ANTES de elegir conformero, para que el numero de
    # atomos del conformero que se copia coincida con el de la molecula.
    ph = R.aplanar(Chem.RemoveHs(Chem.Mol(pose)))
    c = R.aplanar(Chem.RemoveHs(cry))
    if ph.GetNumAtoms() != c.GetNumAtoms():
        return None
    if ph.GetNumConformers() > 1:
        p = Chem.Mol(ph)
        p.RemoveAllConformers()
        p.AddConformer(ph.GetConformer(conf_id))
    else:
        p = ph
    matches = p.GetSubstructMatches(c, uniquify=False, maxMatches=500,
                                    useChirality=False)
    if not matches:
        return None
    m = matches[0]
    pose2 = Chem.RenumberAtoms(p, list(m))   # ahora el orden es el del cristal
    return rdMolAlign.CalcRMS(pose2, c)


def main():
    casos = []
    for entrada in R.SISTEMAS:
        cry = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
        for patron in ("%s_caja24_e8_s*.pdbqt", "%s_vina_caja24_e8_s*.pdbqt",
                       "%s_flex21-22-23-30-100_vina_caja24_e8_s*.pdbqt"):
            casos += [(cry, f) for f in sorted(glob.glob(
                os.path.join(BASE, "out", patron % entrada)))]

    print("%-58s %-10s %-10s %-10s" % ("fichero", "vieja", "corregida", "RDKit"))
    peor = 0.0
    for cry, f in casos:
        pose = R.molecula_dockeada(f)
        if pose is None:
            continue
        v = version_vieja(cry, pose)
        n = R.rmsd_en_sitio(cry, pose)
        r = rmsd_rdkit(cry, pose)
        if n is not None and r is not None:
            peor = max(peor, abs(n - r))
        print("%-58s %-10s %-10s %-10s"
              % (os.path.basename(f),
                 "%.2f" % v if v is not None else "n/d",
                 "%.2f" % n if n is not None else "n/d",
                 "%.2f" % r if r is not None else "n/d"))
    print("\ndiferencia maxima entre la version corregida y CalcRMS de RDKit: %.4f A"
          % peor)
    print("(si es 0, el medidor corregido es exactamente el de RDKit)")


if __name__ == "__main__":
    main()
