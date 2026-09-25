# -*- coding: utf-8 -*-
"""Corrige _recortar_pdb_canonical para agrupar por (cadena, residuo)."""
import io

P = "/home/ubuntu/mmgbsa/rescoring_mmgbsa_robusto.py"
with io.open(P, "r", encoding="utf-8") as f:
    txt = f.read()

old = '''    import math
    keep = set()
    rec_atoms = []
    for line in open(pdb_path, errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resnum = line[22:26].strip()
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        d = math.sqrt((x - lig_center[0])**2 + (y - lig_center[1])**2 + (z - lig_center[2])**2)
        if d < radio:
            keep.add(resnum)
        rec_atoms.append(line)
    lines = []
    serial = 0
    for line in rec_atoms:
        if line[22:26].strip() in keep:
            serial += 1
            lines.append(line[:6] + "%5d" % serial + line[11:].rstrip())'''

new = '''    import math
    keep = set()
    rec_atoms = []
    for line in open(pdb_path, errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        chain = line[21:22].strip() or " "
        resnum = line[22:26].strip()
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        d = math.sqrt((x - lig_center[0])**2 + (y - lig_center[1])**2 + (z - lig_center[2])**2)
        if d < radio:
            keep.add((chain, resnum))
        rec_atoms.append(line)
    lines = []
    serial = 0
    for line in rec_atoms:
        chain = line[21:22].strip() or " "
        if (chain, line[22:26].strip()) in keep:
            serial += 1
            lines.append(line[:6] + "%5d" % serial + line[11:].rstrip())'''

assert old in txt, "no encontrado recorte"
txt = txt.replace(old, new, 1)

with io.open(P, "w", encoding="utf-8") as f:
    f.write(txt)
print("recorte por cadena aplicado")
