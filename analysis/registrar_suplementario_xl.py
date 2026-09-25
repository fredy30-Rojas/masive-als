#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""registrar_suplementario_xl.py — guarda con procedencia los SMILES de XL21 a XL27.

POR QUE EXISTE
--------------
La Figura Suplementaria 1(a) de Gao 2026 (Nature Aging 6:1667) dibuja los ocho hits del
cribado virtual sobre el CR de TDP-43 (XL20 a XL27) y la tabla del panel (b) da el codigo
de catalogo de Asinex de cada uno. De los ocho, solo XL20 tiene union medida y ya vive en
`controles_tdp43_xl20.csv`; los otros siete no son positivos, pero su estructura y su
codigo conviene tenerlos guardados con la procedencia, para no volver a leerlos y para que
nadie los cuente de memoria.

Los SMILES no se escriben aqui a mano: salen del OCR quimico ya hecho
(`_xl20_fig/smiles_crudos.csv`, DECIMER 2.7.2 sobre la figura) y este script solo los
comprueba (RDKit), les saca formula, esqueleto de Murcko y quimiotipo, y los une con los
codigos de Asinex y con la evidencia funcional que el articulo da de cada uno.

QUE NO SON
----------
No son positivos de union medida: entran en la tabla de hits de un acoplamiento. Lo que el
articulo tiene de ellos es funcional (inhibicion de la agregacion del LCD in vitro a
100 uM para XL21 y XL23; neuroproteccion solo a 100 uM para XL27; XL24 empeoro las cosas)
o nada. Por eso NO entran en `verdad_de_referencia.csv`: esa verdad es de union medida.

Uso:
    python analysis/registrar_suplementario_xl.py

Salida: analysis/suplementario_tdp43_xl21_27.csv
"""
import csv
import os
import sys

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

BASE = os.path.dirname(os.path.abspath(__file__))
CRUDOS = os.path.join(BASE, "_xl20_fig", "smiles_crudos.csv")
SALIDA = os.path.join(BASE, "suplementario_tdp43_xl21_27.csv")

CITA = "Gao et al. 2026, Nature Aging 6:1667 (Figura Suplementaria 1a y tabla del panel b)"
FIGURA = "Nature Aging 6:1667, Figura Suplementaria 1a"
SIN_UNION = ("cribado virtual sobre el CR (Asinex BioDesign): acoplamiento, "
             "sin union medida")
NOTA = ("SMILES leido por OCR quimico (DECIMER 2.7.2) sobre la imagen del panel (a); dos "
        "recortes distintos dan el mismo canonico. Codigo de Asinex leido por OCR de la "
        "tabla del panel (b), en dos pasadas coincidentes, y el de XL20 lo confirma el "
        "texto del articulo. No son positivos: no tienen union medida, asi que no entran "
        "en verdad_de_referencia.csv")

# Codigo de Asinex de cada uno (tabla del panel b) y la evidencia funcional que da el
# articulo. Nada de esto se deduce: esta en el texto o en la tabla.
DATOS = {
    "XL21": ("BDE32163394", "inhibe la agregacion del LCD in vitro a 100 uM; sin union medida"),
    "XL22": ("BDE32165749", "sin actividad reportada en el articulo"),
    "XL23": ("LAS51502065", "inhibe la agregacion del LCD in vitro a 100 uM; sin union medida"),
    "XL24": ("BDE31265220", "empeora la muerte neuronal en el ensayo de neuroproteccion; sin union medida"),
    "XL25": ("BDF33572449", "sin actividad reportada en el articulo"),
    "XL26": ("BDE32053214", "sin actividad reportada en el articulo"),
    "XL27": ("BDE32081987", "neuroprotege solo a 100 uM (XL20 lo hace a 6,25 uM); sin union medida"),
}


def quimiotipo(mol, esqueleto: str) -> str:
    """Etiqueta legible del quimiotipo, a partir de lo que dice el propio esqueleto.

    Los patrones estan comprobados contra las siete estructuras: el de la adenina solo
    engancha XL23 y el de la diona engancha los otros seis (XL25 con la diona fundida a un
    benceno, es decir una quinazolinadiona). Se escribe la familia y no "esqueleto propio"
    porque siete etiquetas distintas para dos familias no dicen nada.
    """
    if mol.HasSubstructMatch(Chem.MolFromSmarts("n1cnc2c(N)ncnc12")):
        return "adenina-aminociclohexanol (como XL20)"
    if mol.HasSubstructMatch(Chem.MolFromSmarts("O=c1[nH]c(=O)c2ccccc2[nH]1")):
        return "quinazolinadiona-sulfonamida"
    if mol.HasSubstructMatch(Chem.MolFromSmarts("O=c1[nH]c(=O)c(c[nH]1)")):
        return "uracilo-carboxamida"
    if esqueleto:
        return "esqueleto propio: %s" % esqueleto[:40]
    return "sin clasificar"


def main() -> int:
    if not os.path.exists(CRUDOS):
        sys.exit("falta %s (corre antes leer_figura_xl20.py)" % CRUDOS)
    leidos = {r["etiqueta"]: r["smiles_crudo"] for r in
              csv.DictReader(open(CRUDOS, encoding="utf-8"))}

    filas = []
    for etiqueta in sorted(DATOS):
        asinex, funcional = DATOS[etiqueta]
        smi = leidos.get(etiqueta)
        mol = Chem.MolFromSmiles(smi) if smi else None
        if mol is None:
            print("  %s: SMILES ilegible o invalido, no entra" % etiqueta)
            continue
        esc = MurckoScaffold.MurckoScaffoldSmiles(mol=mol)
        filas.append({
            "ligand": etiqueta,
            "smiles": Chem.MolToSmiles(mol),
            "formula": rdMolDescriptors.CalcMolFormula(mol),
            "pesados": mol.GetNumHeavyAtoms(),
            "anillos": rdMolDescriptors.CalcNumRings(mol),
            "quimia": quimiotipo(mol, esc),
            "esqueleto_murcko": esc,
            "asinex_id": asinex,
            "origen": SIN_UNION,
            "tipo_ensayo": funcional,
            "apto": "no",
            "cita": CITA,
            "figura": FIGURA,
            "nota": NOTA,
        })

    with open(SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        w.writeheader()
        w.writerows(filas)

    print("%d compuestos guardados en %s" % (len(filas), SALIDA))
    for r in filas:
        print("  %-5s %-14s %-9s %2d anillos  %-38s %s"
              % (r["ligand"], r["asinex_id"], r["formula"], r["anillos"],
                 r["quimia"], r["tipo_ensayo"][:44]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
