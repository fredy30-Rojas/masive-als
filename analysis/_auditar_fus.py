# -*- coding: utf-8 -*-
"""Auditoría de poses MM-GBSA extremas (FUS −30 a −32).

Comprueba colisiones estéricas: distancia mínima ligando-receptor por átomo
y si alguna pareja está por debajo de la suma de radios de van der Waals
(colisión real). Una pose con colisiones fuertes da energías MM-GBSA
extremas y no es válida.
"""
import glob, math, os

REC = r"C:\Users\Fredy\masive-als\analysis\_rescoring_stage2\receptores\FUS.pdbqt"
POSES = r"C:\Users\Fredy\masive-als\analysis\_rescoring_stage2\poses"

# radios vdW (A) por elemento
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80,
       "F": 1.47, "CL": 1.75, "BR": 1.85, "I": 1.98, "H": 1.20}

def elem_from_atom(line):
    # pdbqt: columna 77-80 = tipo AD; coger el elemento de ahi o de 76-78
    ad = line[77:80].strip().upper()
    # los tipos AD vienen con mayusculas; extraer elemento principal
    if ad.startswith("CL"): return "CL"
    if ad.startswith("BR"): return "BR"
    if ad.startswith("NA") or ad.startswith("NS") or ad.startswith("NX"): return "N"
    if ad.startswith("OA") or ad.startswith("OS"): return "O"
    if ad.startswith("SA"): return "S"
    if ad.startswith("HD") or ad.startswith("HS"): return "H"
    if ad: return ad[0]
    return "C"

def parse(pdbqt):
    atoms = []
    for line in open(pdbqt):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        e = elem_from_atom(line)
        atoms.append((x, y, z, e))
    return atoms

rec = parse(REC)
print(f"Receptor FUS: {len(rec)} átomos")

targets = ["CHEMBL7563_FUS", "CHEMBL9010_FUS", "CHEMBL8905_FUS",
           "CHEMBL3309490_FUS", "CHEMBL7652_FUS", "CHEMBL3310332_FUS",
           "CHEMBL264055_FUS"]

# rango de contactos: minima distancia por par
for t in targets:
    p = os.path.join(POSES, t + "_out.pdbqt")
    if not os.path.exists(p):
        print(f"{t}: pose NO EXISTE")
        continue
    lig = parse(p)
    # distancia mínima lig-receptor
    dmin = 1e9; worst = None
    clashes = 0
    for (lx, ly, lz, le) in lig:
        for (rx, ry, rz, re_) in rec:
            d = math.sqrt((lx-rx)**2 + (ly-ry)**2 + (lz-rz)**2)
            if d < dmin:
                dmin = d; worst = (le, re_)
            # colision si d < 0.75 * (r1+r2)
            rsum = (VDW.get(le, 1.7) + VDW.get(re_, 1.7)) * 0.75
            if d < rsum:
                clashes += 1
    print(f"{t}: {len(lig)} átomos lig | d_min_receptor={dmin:.2f} A "
          f"({worst[0]}-{worst[1]}) | contactos<0.75vdW={clashes}")
