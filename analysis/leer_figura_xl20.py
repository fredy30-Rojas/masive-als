#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""leer_figura_xl20.py — leer por OCR quimico las ocho estructuras de la Figura Suplementaria 1(a).

POR QUE EXISTE
--------------
XL20 es, hasta hoy, el unico ligando con union medida del dominio C-terminal de TDP-43
(SPR y CETSA, dependiente del Trp334; Gao et al. 2026, Nature Aging 6:1667). Su
estructura NO esta en el texto del articulo ni en PubChem: aparece dibujada en la Figura
Suplementaria 1(a), y la tabla del panel (b) da su codigo de catalogo de Asinex
(BDF34019555), que PubChem no reconoce.

Asi que la unica via de tener el SMILES sin inventarlo es leer el dibujo, igual que se
hizo con los tres fragmentos de Nshogoza 2019.

LO QUE DEVUELVE NO ES UN DATO VERIFICADO
----------------------------------------
Un OCR de estructuras se equivoca. Lo que sale de aqui es un CANDIDATO a SMILES que hay
que comprobar (validez en RDKit, dibujo de comparacion, segunda lectura con recorte
distinto y, si se puede, contraste con una fuente externa) antes de que entre en
`verdad_de_referencia.csv`.

COMO ENCUENTRA LAS CELDAS
-------------------------
El panel (a) es una fila de ocho dibujos, etiquetados XL20..XL27 de izquierda a derecha
(comprobado leyendo la tira de etiquetas). Aqui no se adivina el reparto: se corta la
tira de dibujo por los huecos de tinta y se asignan los bloques en orden a las etiquetas,
que es mas robusto que repartir por anchos iguales cuando un dibujo es mas ancho.

Uso:
    ocsr_env/Scripts/python.exe leer_figura_xl20.py <figura.png> [y0 y1]

Salida: analysis/_xl20_fig/ con las celdas recortadas y smiles_crudos.csv.
"""
import csv
import os
import sys

import numpy as np
from PIL import Image

BASE = os.path.dirname(os.path.abspath(__file__))
SALIDA = os.path.join(BASE, "_xl20_fig")
ETIQUETAS = ["XL%d" % n for n in range(20, 28)]


def recortar_celdas(ruta_fig, y0=55, y1=300, umbral=200, hueco_min=25):
    """Devuelve [(etiqueta, imagen)] de los dibujos del panel, en orden."""
    im = Image.open(ruta_fig).convert("L")
    panel = im.crop((0, y0, im.width, y1))
    tinta = np.array(panel) < umbral

    activo = tinta.sum(axis=0) > 1
    bloques, ini, hueco = [], None, 0
    for x, act in enumerate(activo):
        if act:
            if ini is None:
                ini = x
            hueco = 0
        elif ini is not None:
            hueco += 1
            if hueco > hueco_min:
                bloques.append((ini, x - hueco))
                ini = None
    if ini is not None:
        bloques.append((ini, len(activo) - 1))
    bloques = sorted(bloques, key=lambda b: b[1] - b[0], reverse=True)[:len(ETIQUETAS)]
    bloques.sort()

    salida = []
    for etiqueta, (x0, x1) in zip(ETIQUETAS, bloques):
        sub = tinta[:, x0 : x1 + 1]
        filas = np.where(sub.sum(axis=1) > 1)[0]
        if len(filas) == 0:
            continue
        a, b = int(filas.min()), int(filas.max())
        margen = 18
        caja = panel.crop(
            (max(0, x0 - margen), max(0, a - margen),
             min(panel.width, x1 + margen), min(panel.height, b + margen))
        )
        print("   %s: bloque x %d-%d, filas %d-%d (alto %d)"
              % (etiqueta, x0, x1, a, b, b - a + 1))
        lado = max(caja.size) + 80
        lienzo = Image.new("L", (lado, lado), 255)
        lienzo.paste(caja, ((lado - caja.width) // 2, (lado - caja.height) // 2))
        if lado < 1000:
            f = 1000 / lado
            lienzo = lienzo.resize((int(lado * f), int(lado * f)), Image.LANCZOS)
        salida.append((etiqueta, lienzo.convert("RGB")))
    return salida


def comparar(etiqueta, img, smi):
    """Recorte original a la izquierda y dibujo de RDKit a la derecha, para mirarlo.

    No decide nada: es el material con el que se da (o no) el visto bueno. Se hizo asi
    con los fragmentos de Nshogoza, y ahi el cotejo visual fue el que confirmo.
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw

    mol = Chem.MolFromSmiles(smi) if smi else None
    if mol is None:
        print("  %s: SMILES invalido, no hay dibujo de comparacion" % etiqueta)
        return None
    AllChem.Compute2DCoords(mol)
    dibujo = os.path.join(SALIDA, "%s_rdkit.png" % etiqueta)
    Draw.MolToFile(mol, dibujo, size=(700, 700))

    a = img.convert("RGB").resize((700, 700))
    b = Image.open(dibujo).convert("RGB").resize((700, 700))
    comp = Image.new("RGB", (1420, 700), "white")
    comp.paste(a, (0, 0))
    comp.paste(b, (720, 0))
    salida = os.path.join(SALIDA, "%s_comparacion.png" % etiqueta)
    comp.save(salida)
    return salida


def main():
    if len(sys.argv) < 2:
        sys.exit("uso: leer_figura_xl20.py <figura.png> [y0 y1]")
    y0 = int(sys.argv[2]) if len(sys.argv) > 2 else 55
    y1 = int(sys.argv[3]) if len(sys.argv) > 3 else 300
    os.makedirs(SALIDA, exist_ok=True)

    celdas = recortar_celdas(sys.argv[1], y0, y1)
    print("celdas encontradas: %d" % len(celdas))
    if len(celdas) != len(ETIQUETAS):
        print("AVISO: se esperaban 8 (XL20..XL27)")

    from DECIMER import predict_SMILES

    filas = []
    for etiqueta, img in celdas:
        png = os.path.join(SALIDA, "%s.png" % etiqueta)
        img.save(png)
        smi, conf = None, None
        try:
            salida = predict_SMILES(png, confidence=True)
            if isinstance(salida, tuple):
                smi, conf = salida[0], salida[1]
            else:
                smi = salida
        except Exception as e:
            print("  %s: ERROR %s" % (etiqueta, str(e)[:200]))
        print("  %s -> %s   (confianza %s)" % (etiqueta, smi, conf))
        comparar(etiqueta, img, smi)
        filas.append({"etiqueta": etiqueta, "png": os.path.basename(png),
                      "smiles_crudo": smi, "confianza": conf})

    with open(os.path.join(SALIDA, "smiles_crudos.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["etiqueta", "png", "smiles_crudo", "confianza"])
        w.writeheader()
        w.writerows(filas)
    print("guardado: %s" % os.path.join(SALIDA, "smiles_crudos.csv"))


if __name__ == "__main__":
    main()
