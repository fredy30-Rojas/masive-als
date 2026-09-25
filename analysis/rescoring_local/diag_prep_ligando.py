# -*- coding: utf-8 -*-
"""Diagnostico: por que algunas poses de la lista focalizada tienen 2 atomos de mas.

Compara, para los ligandos que se le pasen:
  1. lo que dice el SMILES de la lista (RDKit: atomos pesados y elementos),
  2. lo que hay en la pose acoplada (results_TDP43_v2/<lig>_out.pdbqt),
  3. lo que produce la receta de preparado (convertir(): AddHs + ETKDG + MMFF +
     meeko + sanear_tipos), que es la que se uso para acoplar la focalizada.

Uso: python diag_prep_ligando.py CHEMBL28028 CHEMBL3921127 ...
"""
import collections
import csv
import os
import sys
import tempfile

sys.path.insert(0, r"C:\Users\Fredy\masive-als\gpu_dock")
from acoplar_controles_tdp43 import convertir          # noqa: E402

from rdkit import Chem, RDLogger                       # noqa: E402
RDLogger.DisableLog("rdApp.*")

AQUI = os.path.dirname(os.path.abspath(__file__))
LISTA = os.path.join(AQUI, "lista_focalizada_rescoring.csv")
POSES = r"C:\Users\Fredy\masive-als\gpu_dock\rescoring_caja_arn\results_TDP43_v2"
SIN_H = ("H", "HD", "HS")


def leer_pose(path):
    tipos = collections.Counter()
    nombres = []
    n = 0
    for line in open(path, encoding="utf-8", errors="replace"):
        line = line.rstrip()
        if line.startswith("MODEL") and n:
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        t = line[77:80].strip().upper()
        tipos[t] += 1
        if t not in SIN_H:
            n += 1
            nombres.append((line[12:16].strip(), t))
    return n, tipos, nombres


def contar_pdbqt_texto(txt):
    tipos = collections.Counter()
    nombres = []
    n = 0
    for line in txt.splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        t = line.rstrip().rsplit(" ", 1)[-1]      # el tipo es el ultimo token
        tipos[t] += 1
        if t not in SIN_H:
            n += 1
            nombres.append((line[12:16].strip(), t))
    return n, tipos, nombres


def main(nombres_pedidos):
    filas = {r["ligand"]: r for r in csv.DictReader(open(LISTA, encoding="utf-8"))}
    tmp = tempfile.mkdtemp(prefix="diagprep_")
    for lig in nombres_pedidos:
        r = filas.get(lig)
        print("=" * 70)
        if not r:
            print(lig, ": no esta en la lista")
            continue
        smi = r["smiles"]
        mol = Chem.MolFromSmiles(smi)
        elems = collections.Counter(a.GetSymbol() for a in mol.GetAtoms())
        print("%s: SMILES de la lista -> %d atomos pesados %s"
              % (lig, mol.GetNumAtoms(), dict(elems)))

        pose = os.path.join(POSES, lig + "_out.pdbqt")
        if os.path.exists(pose):
            n, tipos, nom = leer_pose(pose)
            extra = [x for x in nom if x[1] not in elems or x[0] == "G"]
            print("   pose acoplada -> %d atomos pesados | tipos %s" % (n, dict(tipos)))
            print("   nombres raros en la pose: %s" % extra[:6])
        else:
            print("   sin pose en %s" % POSES)

        ruta, err = convertir("TEST_" + lig, smi, tmp)
        if not ruta:
            print("   convertir() FALLO: %s" % err)
            continue
        txt = open(ruta, encoding="utf-8").read()
        n2, tipos2, nom2 = contar_pdbqt_texto(txt)
        print("   receta convertir() -> %d atomos pesados | tipos %s" % (n2, dict(tipos2)))
        print("   nombres raros del preparado: %s"
              % [x for x in nom2 if x[1] not in elems or x[0] == "G"][:6])
        malos = [l for l in txt.splitlines() if l.startswith(("ATOM", "HETATM"))
                 and l.rstrip().rsplit(" ", 1)[-1] not in
                 {"C", "A", "N", "NA", "NS", "NX", "O", "OA", "OS", "F", "Cl",
                  "Br", "I", "P", "S", "SA", "HD", "HS"}]
        print("   atomos con tipo no valido: %d %s" % (len(malos), malos[:3]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:] or ["CHEMBL28028", "CHEMBL3921127"]))
