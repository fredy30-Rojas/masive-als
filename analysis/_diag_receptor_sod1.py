#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Por qué falla la preparación del receptor en algunos recortes.

Mira el receptor canónico de SOD1: residuos incompletos (THR sobre todo), y
compara el recorte de 14 A de los dos ligandos que fallaron con el que sí va.
"""
import math
import sys
from collections import OrderedDict

REF = {"ALA": 5, "ARG": 11, "ASN": 8, "ASP": 8, "CYS": 6, "GLN": 9, "GLU": 9,
       "GLY": 4, "HIS": 10, "ILE": 8, "LEU": 8, "LYS": 9, "MET": 8, "PHE": 11,
       "PRO": 7, "SER": 6, "THR": 7, "TRP": 14, "TYR": 12, "VAL": 7}


def leer(pdb):
    res = OrderedDict()
    for line in open(pdb, errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        key = (line[21:22], line[22:26].strip())
        res.setdefault(key, []).append(line)
    return res


def main(pdb):
    res = leer(pdb)
    print(f"{pdb}: {len(res)} residuos, {sum(len(v) for v in res.values())} atomos\n")
    malos = []
    for (chain, num), lines in res.items():
        rn = lines[0][17:20].strip()
        esperado = REF.get(rn)
        nombres = [l[12:16].strip() for l in lines]
        if esperado and len(nombres) != esperado:
            malos.append((chain, num, rn, len(nombres), esperado, nombres))
    print(f"residuos con nº de atomos distinto al estandar: {len(malos)}")
    for chain, num, rn, n, esp, nombres in malos:
        print(f"  cadena {chain!r} res {num} {rn}: {n} atomos (esperado {esp}) "
              f"-> {','.join(nombres)}")

    # reparto de cadenas
    cadenas = {}
    for (chain, _), lines in res.items():
        rn = lines[0][17:20].strip()
        cadenas.setdefault(chain, []).append(int(res[ (chain, _)][0][22:26]))
    print("\ncadenas:")
    for ch, nums in cadenas.items():
        faltas = []
        nums_s = sorted(nums)
        for a, b in zip(nums_s, nums_s[1:]):
            if b != a + 1:
                faltas.append(f"{a}->{b}")
        print(f"  {ch!r}: {len(nums)} res, {min(nums)}-{max(nums)}, huecos={faltas}")
    return res


if __name__ == "__main__":
    main(sys.argv[1])
