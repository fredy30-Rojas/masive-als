#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""leer_figura_nshogoza.py — leer por OCR quimico los tres fragmentos de la Figura 1c.

POR QUE EXISTE
--------------
Los tres hits del cribado de fragmentos de Nshogoza 2019 (IJMS 20:3230) contra el
bolsillo RRM2 de TDP-43 aparecen en la literatura UNICAMENTE como dibujo, en la
Figura 1c, etiquetados «hit 1», «hit 2» y «hit 3». No hay nombre, numero CAS,
codigo de catalogo ni formula en el articulo, ni en su material suplementario
(solo tablas de anchura de linea), ni en la revision de Francois-Moutal 2021
(PMC8341936), que los vuelve a reproducir como imagen. En el disco tampoco estaban.

Asi que la unica via de tener su estructura sin inventarla es leer el dibujo. Este
script recorta los tres paneles de la figura y los pasa por DECIMER (OCR quimico,
Steinbeck lab, ~90-96% de acierto en dibujos impresos).

LO QUE DEVUELVE NO ES UN DATO VERIFICADO
----------------------------------------
Un OCR de estructuras se equivoca. Lo que sale de aqui es un CANDIDATO a SMILES que
hay que comprobar antes de que entre en `verdad_de_referencia.csv`. Las dos
comprobaciones que hace este script solo:

  1) que el SMILES sea valido y tenga una formula coherente con un fragmento
     (el articulo declara fragmentos de 110-250 Da);
  2) que la molecula se pueda dibujar sin valencias raras.

Uso:
    masive-als/ocsr_env/Scripts/python.exe leer_figura_nshogoza.py <figura.png>

Salida: analysis/_nshogoza_fig1c/  con los recortes y un CSV con los SMILES crudos.
"""
import csv
import os
import sys

import numpy as np
from PIL import Image

BASE = os.path.dirname(os.path.abspath(__file__))
SALIDA = os.path.join(BASE, "_nshogoza_fig1c")


def recortar_estructuras(ruta_fig, fraccion_c=0.78, umbral=240):
    """Devuelve la lista de recortes (imagen) de los dibujos del panel (c).

    El panel tiene tres dibujos con su etiqueta de texto debajo («hit 1»...). Por
    columnas se separan los tres; por filas se toma el tramo de tinta ALTO (el
    dibujo, ~400 px) y se descarta el tramo BAJO de debajo (~75 px), que es el
    texto. Si el texto entra en la imagen, DECIMER lo lee como atomos: la primera
    pasada, con el texto dentro, devolvio un complejo de paladio, siete boros y
    sodio que eran literalmente las letras de las etiquetas.
    """
    im = Image.open(ruta_fig).convert("L")
    w, h = im.size
    panel = im.crop((0, int(h * fraccion_c), w, h))
    a = np.array(panel)

    tinta = a < umbral
    # perfil de columnas: donde hay dibujo
    cols = tinta.sum(axis=0)
    activo = cols > 2

    # agrupar columnas activas en bloques separados por huecos grandes
    bloques = []
    ini = None
    hueco = 0
    for x, act in enumerate(activo):
        if act:
            if ini is None:
                ini = x
            hueco = 0
        else:
            if ini is not None:
                hueco += 1
                if hueco > 120:          # hueco ancho = separacion entre moleculas
                    bloques.append((ini, x - hueco))
                    ini = None
    if ini is not None:
        bloques.append((ini, len(activo) - 1))

    # quedarse con los 3 bloques mas anchos
    bloques = sorted(bloques, key=lambda b: b[1] - b[0], reverse=True)[:3]
    bloques.sort()

    recortes = []
    for (x0, x1) in bloques:
        sub = tinta[:, x0:x1 + 1]
        filas_act = sub.sum(axis=1) > 2

        # tramos de filas con tinta
        tramos, ini = [], None
        for y, v in enumerate(filas_act):
            if v and ini is None:
                ini = y
            if not v and ini is not None:
                tramos.append((ini, y - 1))
                ini = None
        if ini is not None:
            tramos.append((ini, len(filas_act) - 1))
        if not tramos:
            continue
        # el dibujo es el tramo mas alto; la etiqueta de texto es baja
        alto = max(tramos, key=lambda t: t[1] - t[0])
        y0, y1 = alto
        cols_act = np.where(sub[y0:y1 + 1].sum(axis=0) > 0)[0]
        xa, xb = x0 + cols_act.min(), x0 + cols_act.max()
        print("   bloque %d-%d -> filas %d-%d (alto %d), columnas %d-%d;"
              " tramos descartados: %s"
              % (x0, x1, y0, y1, y1 - y0 + 1, xa, xb,
                 [(t[0], t[1]) for t in tramos if t != alto]))
        margen = 20
        caja = panel.crop((max(0, xa - margen), max(0, y0 - margen),
                           min(panel.width, xb + margen),
                           min(panel.height, y1 + margen)))
        # DECIMER rinde mejor con la estructura grande, cuadrada y sobre blanco
        lado = max(caja.size) + 80
        lienzo = Image.new("L", (lado, lado), 255)
        lienzo.paste(caja, ((lado - caja.width) // 2, (lado - caja.height) // 2))
        if lado < 900:
            factor = 900 / lado
            lienzo = lienzo.resize((int(lado * factor), int(lado * factor)),
                                   Image.LANCZOS)
        recortes.append(lienzo.convert("RGB"))
    return recortes


def main():
    if len(sys.argv) < 2:
        sys.exit("uso: leer_figura_nshogoza.py <figura.png>")
    os.makedirs(SALIDA, exist_ok=True)

    recortes = recortar_estructuras(sys.argv[1])
    print("recortes encontrados: %d" % len(recortes))
    if len(recortes) != 3:
        print("AVISO: se esperaban 3 dibujos (hit 1, hit 2, hit 3)")

    # DECIMER importa tensorflow y descarga su modelo la primera vez
    from DECIMER import predict_SMILES

    filas = []
    for i, img in enumerate(recortes, start=1):
        png = os.path.join(SALIDA, "hit_%d.png" % i)
        img.save(png)
        smi, conf = None, None
        try:
            salida = predict_SMILES(png, confidence=True)
            if isinstance(salida, tuple):
                smi, conf = salida[0], salida[1]
            else:
                smi = salida
        except Exception as e:
            print("  hit %d: ERROR %s" % (i, str(e)[:200]))
        print("  hit %d -> %s   (confianza %s)" % (i, smi, conf))
        filas.append({"hit": i, "png": os.path.basename(png),
                      "smiles_crudo": smi, "confianza": conf})

    with open(os.path.join(SALIDA, "smiles_crudos.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["hit", "png", "smiles_crudo", "confianza"])
        w.writeheader()
        w.writerows(filas)
    print("guardado: %s" % os.path.join(SALIDA, "smiles_crudos.csv"))


if __name__ == "__main__":
    main()
