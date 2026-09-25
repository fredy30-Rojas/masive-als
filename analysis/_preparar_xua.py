# -*- coding: utf-8 -*-
"""Prepara el ligando nativo XUA de 4IUF (control gold-standard de TDP-43).

1. Extrae HETATM XUA del PDB 4IUF a un PDB limpio.
2. Lo convierte a PDBQT (obabel + meeko) en el mismo frame.
"""
import os, shutil, subprocess, sys, tempfile

SRC = r"C:\Users\Fredy\masive-als\analysis\_4iuf_full.pdb"
OUT_PDB = r"C:\Users\Fredy\masive-als\analysis\_xua_nativo.pdb"
OUT_PDBQT = r"C:\Users\Fredy\masive-als\analysis\_xua_nativo.pdbqt"

if not os.path.exists(SRC):
    print("descarga 4iuf primero"); sys.exit(1)

lines = []
for line in open(SRC):
    if line.startswith("HETATM") and "XUA" in line[17:20]:
        lines.append(line.rstrip())
if not lines:
    print("XUA no encontrado"); sys.exit(1)

with open(OUT_PDB, "w") as f:
    f.write("\n".join(lines) + "\nEND\n")
print("XUA extraido: %d atomos -> %s" % (len(lines), OUT_PDB))

# a PDBQT con obabel (mismo frame, hidrogenos polares + cargas)
obabel = r"C:\Users\Fredy\AppData\Local\Python\pythoncore-3.14-64\Scripts\obabel.exe"
cmd = [obabel, OUT_PDB, "-O", OUT_PDBQT, "-p", "7.4", "--partialcharge", "gasteiger"]
r = subprocess.run(cmd, capture_output=True, text=True)
print("obabel rc:", r.returncode)
if r.returncode != 0:
    print(r.stderr[:500])
else:
    n = sum(1 for l in open(OUT_PDBQT) if l.startswith(("ATOM", "HETATM")))
    print("PDBQT XUA: %d atomos -> %s" % (n, OUT_PDBQT))
