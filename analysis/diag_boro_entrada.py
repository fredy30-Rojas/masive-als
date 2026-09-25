# -*- coding: utf-8 -*-
"""¿Cuántas moléculas se acoplaron sin ser ellas mismas? Mirada directa.

El intento anterior emparejaba la librería de SMILES contra las poses por
CHEMBL ID, y fallaba: los ficheros de ligando se llaman por NOMBRE
("ABACAVIR.pdbqt"), no por ID, así que sólo se emparejaron 24 de 349.

Este script no empareja nada. Va directo a los PDBQT DE ENTRADA, que son
literalmente las moléculas que se acoplaron, y pregunta por cada una:

    ¿tiene algún átomo NOMBRADO "B"?  y si lo tiene, ¿con qué TIPO?

En AutoDock el tipo (la última columna) es lo que manda al puntuar, no el
nombre. Un átomo llamado "B" pero tipado "C" es el desastre que se busca: la
malla se construyó como si fuera carbono, y la molécula se acopló y puntuó
como su análogo de carbono. Eso importa especialmente en los inhibidores de
proteasoma (bortezomib y compañía), donde el boro ES el grupo que se une.

Los PDBQT de entrada conservan el nombre del átomo igual que la pose, así que
esto cubre TODA la librería, no una muestra.

Uso:  python diag_boro_entrada.py
"""
import os
import re
import sys
from collections import Counter

BASE = r"C:\Users\Fredy\masive-als"
if not os.path.isdir(BASE):
    BASE = "/mnt/c/Users/Fredy/masive-als"
LIGDIR = os.path.join(BASE, "gpu_dock", "libreria_ligands")

# Un átomo "con boro en el nombre": B suelto, o BB (diboro), o B con número.
NOMBRE_BORO = re.compile(r"^B[BD0-9]?$")


def revisar(ruta):
    """Devuelve (nombres_B, tipos_de_esos_B)."""
    nombres, tipos = [], []
    try:
        with open(ruta, encoding="utf-8", errors="replace") as f:
            for linea in f:
                if not linea.startswith(("ATOM", "HETATM")):
                    continue
                nombre = linea[12:16].strip().upper()
                if NOMBRE_BORO.match(nombre):
                    nombres.append(nombre)
                    tipos.append(linea[77:79].strip().upper())
    except OSError:
        return [], []
    return nombres, tipos


def main():
    if not os.path.isdir(LIGDIR):
        print("no encuentro %s" % LIGDIR)
        return 1

    ficheros = sorted(n for n in os.listdir(LIGDIR) if n.endswith(".pdbqt"))
    print("Revisando %d PDBQT de entrada (las moléculas que se acoplaron)..."
          % len(ficheros))

    con_boro = 0
    bien = []       # nombrado B y tipado B
    mal = []        # nombrado B pero tipado C (u otra cosa)
    cuantos_tipos = Counter()

    for i, nombre in enumerate(ficheros):
        if i and i % 20000 == 0:
            print("  ...%d de %d" % (i, len(ficheros)), flush=True)
        nombres, tipos = revisar(os.path.join(LIGDIR, nombre))
        if not nombres:
            continue
        con_boro += 1
        for t in tipos:
            cuantos_tipos[t] += 1
        if all(t == "B" for t in tipos):
            bien.append(nombre[:-len(".pdbqt")])
        else:
            mal.append((nombre[:-len(".pdbqt")],
                        "%d átomos B, tipos %s" % (len(tipos), sorted(set(tipos)))))

    print()
    print("=" * 74)
    print("RESULTADO")
    print("=" * 74)
    print("  PDBQT de entrada revisados      : %d" % len(ficheros))
    print("  con algún átomo nombrado B      : %d" % con_boro)
    print("  B bien tipado como B            : %d" % len(bien))
    print("  B MAL tipado (no es boro)       : %d   <-- el problema" % len(mal))
    print("  tipos encontrados en esos átomos: %s" % dict(cuantos_tipos))
    print()
    if mal:
        print("  COMPUESTOS AFECTADOS:")
        for n, detalle in mal:
            print("     %-40s %s" % (n[:40], detalle))
    if bien:
        print()
        print("  (los que están bien, para comparar:)")
        for n in bien[:10]:
            print("     %s" % n[:60])
    return 0


if __name__ == "__main__":
    sys.exit(main())
