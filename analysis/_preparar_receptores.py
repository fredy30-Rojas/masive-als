#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara receptores simplificados para el rescoring MM-GBSA.

- TDP43: cadena A de 6b1n (aproximación; lo ideal es el dominio RRM1).
- SOD1: cadenas A y B de 1hl5 (dímero biológico; lo ideal es el sitio anti-agregación).
- FUS: modelo 1 del ensamble NMR 6g99.

Nota: estos son receptores SIMPLIFICADOS para montar el pipeline. El refinado
a los dominios exactos del docking (RRM1, sitio anti-agregación) queda pendiente.
"""
import os

OUT = os.path.dirname(os.path.abspath(__file__))

def extract(src, dst, chains=None, model=None):
    with open(src, 'r', encoding='utf-8', errors='replace') as f:
        lines = f.readlines()
    out = []
    in_model = True
    current_model = 0
    for line in lines:
        if model is not None:
            if line.startswith('MODEL'):
                current_model = int(line.split()[1])
                in_model = (current_model == model)
                continue
            if line.startswith('ENDMDL'):
                if in_model:
                    break
                continue
            if not in_model:
                continue
        if line.startswith('ATOM') or line.startswith('HETATM'):
            if chains is None:
                out.append(line)
            else:
                chain = line[21]
                if chain in chains:
                    out.append(line)
        elif line.startswith('TER'):
            out.append(line)
    with open(dst, 'w', encoding='utf-8') as f:
        f.write(''.join(out))
    print(f"{os.path.basename(src)} -> {os.path.basename(dst)} ({len(out)} lineas)")

base = r"C:\Users\Fredy\masive-als\proteins"

extract(os.path.join(base, "TDP43", "PDB-6b1n.pdb"),
        os.path.join(OUT, "receptor_TDP43.pdb"), chains=('A',))
extract(os.path.join(base, "SOD1", "PDB-1hl5.pdb"),
        os.path.join(OUT, "receptor_SOD1.pdb"), chains=('A', 'B'))
extract(os.path.join(base, "FUS", "PDB-6g99.pdb"),
        os.path.join(OUT, "receptor_FUS.pdb"), model=1)

print("Receptores preparados en", OUT)
