# -*- coding: utf-8 -*-
"""Texto a PDF, para los documentos del entregable (20 sep 2026).

No hay pandoc en este equipo, asi que se hace con fpdf2 y una fuente TrueType
del sistema, para que los acentos y las comillas salgan bien.

Uso: python txt_a_pdf.py entrada.txt [salida.pdf]
"""
import os
import sys

from fpdf import FPDF

FUENTES = [r"C:\Windows\Fonts\arial.ttf",
           r"C:\Windows\Fonts\consola.ttf",
           r"C:\Windows\Fonts\calibri.ttf"]
TAMANO = 9.5
INTERLINEA = 4.6


def fuente_disponible():
    for f in FUENTES:
        if os.path.exists(f):
            return f
    return None


def convertir(entrada, salida=None):
    salida = salida or os.path.splitext(entrada)[0] + ".pdf"
    ruta_fuente = fuente_disponible()
    pdf = FPDF(format="A4")
    pdf.set_auto_page_break(auto=True, margin=14)
    pdf.add_page()
    if ruta_fuente:
        pdf.add_font("doc", "", ruta_fuente)
        pdf.set_font("doc", size=TAMANO)
    else:
        pdf.set_font("Courier", size=TAMANO)
    ancho = pdf.epw
    texto = open(entrada, encoding="utf-8").read().replace("\r\n", "\n")
    for linea in texto.split("\n"):
        # fpdf2 no traga la cadena vacia: se manda un espacio
        pdf.multi_cell(ancho, INTERLINEA, linea if linea.strip() else " ")
    pdf.output(salida)
    print("Generado: %s (%d paginas)" % (salida, pdf.pages_count))
    return salida


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    convertir(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
    return 0


if __name__ == "__main__":
    sys.exit(main())
