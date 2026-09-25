#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""comparar_ocr_figura.py — comprobacion visual de lo que devolvio el OCR.

Para cada hit monta una imagen con el recorte ORIGINAL de la figura a la izquierda
y el dibujo de RDKit del SMILES leido a la derecha, para poder compararlos. No
decide nada por si solo: es material para mirar (o para que lo mire el modelo de
vision qwen3-vl:8b) antes de dar el SMILES por bueno.

Uso:
    ocsr_env/Scripts/python.exe comparar_ocr_figura.py
"""
import os

from PIL import Image
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

BASE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(BASE, "_nshogoza_fig1c")

SMILES_LEIDOS = {
    "hit_1": "C1CC(CN(C1)C2=NC=CC=N2)N",
    "hit_2": "CC(C1=CN2C(=N1)SC=N2)N",
    "hit_3": "C1=CN=C(C(=C1)C#N)N2CCC(CC2)N",
}


def main():
    for nombre, smi in SMILES_LEIDOS.items():
        original = os.path.join(DIR, nombre + ".png")
        if not os.path.exists(original):
            print("falta %s" % original)
            continue
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            print("%s: SMILES invalido" % nombre)
            continue
        AllChem.Compute2DCoords(mol)
        dibujo = os.path.join(DIR, nombre + "_rdkit.png")
        Draw.MolToFile(mol, dibujo, size=(700, 700))

        a = Image.open(original).convert("RGB").resize((700, 700))
        b = Image.open(dibujo).convert("RGB").resize((700, 700))
        comp = Image.new("RGB", (1420, 700), "white")
        comp.paste(a, (0, 0))
        comp.paste(b, (720, 0))
        salida = os.path.join(DIR, nombre + "_comparacion.png")
        comp.save(salida)
        print("guardado %s" % salida)


if __name__ == "__main__":
    main()
