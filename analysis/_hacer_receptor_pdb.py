# -*- coding: utf-8 -*-
"""Genera receptor_TDP43_rrm1_canon.pdb desde el pdbqt del docking.

El pdbqt del docking tiene las coordenadas correctas (frame del cribado) y
nombres de átomo canónicos, pero AutoDock los reordena (los átomos de cada
residuo aparecen dispersos y en orden no estándar), lo que rompe PDBFixer
(residuos UNK -> NaN/timeouts).

Solución: reagrupar los átomos por residuo y reordenarlos en el orden
canónico de aminoácidos (N, CA, C, O, luego cadenas laterales). Resultado:
PDB válido, frame del docking, química reconocible por PDBFixer.
"""
import collections

PQT = "rescoring_pkg2/receptores/TDP43.pdbqt"
OUT = "receptor_TDP43_rrm1_canon.pdb"

# orden canónico de los átomos del backbone y preferencia para el resto
BACKBONE_ORDER = {"N": 0, "CA": 1, "C": 2, "O": 3, "OXT": 4}

def elem_of(line):
    ad = line[77:80].strip().upper()
    if ad.startswith("CL"): return "CL"
    if ad.startswith(("NA", "NS", "NX")): return "N"
    if ad.startswith(("OA", "OS")): return "O"
    if ad.startswith(("SA",)): return "S"
    if ad.startswith(("HD", "HS")): return "H"
    return ad[0] if ad else "C"

# agrupar por residuo
residues = collections.OrderedDict()   # resnum -> list of atom dicts
for line in open(PQT):
    if not line.startswith(("ATOM", "HETATM")):
        continue
    resnum = int(line[22:26].strip())
    name = line[12:16].strip()
    resname = line[17:20].strip()
    x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
    elem = elem_of(line)
    residues.setdefault(resnum, []).append(
        dict(name=name, resname=resname, x=x, y=y, z=z, elem=elem))

# reordenar cada residuo
out = []
serial = 0
for resnum in sorted(residues.keys()):
    atoms = residues[resnum]
    resname = atoms[0]["resname"] if atoms else "UNK"
    # hidrógenos al final; backbone primero
    def sort_key(a):
        if a["elem"] == "H":
            return (1, 0, a["name"])
        return (0, BACKBONE_ORDER.get(a["name"], 50), a["name"])
    atoms.sort(key=sort_key)
    for a in atoms:
        serial += 1
        name = a["name"]
        if name[0].isdigit():  # nombres tipo 1HG -> HG1
            name = name[1:] + name[0]
        out.append("ATOM  %5d %-4s %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
                   % (serial, name, resname, resnum, a["x"], a["y"], a["z"], a["elem"]))

with open(OUT, "w") as f:
    f.write("\n".join(out) + "\nEND\n")
print("residuos: %d | atomos: %d -> %s" % (len(residues), serial, OUT))
