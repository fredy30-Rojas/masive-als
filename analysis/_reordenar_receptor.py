# -*- coding: utf-8 -*-
"""Reordena los pdbqt del docking a PDB canónico (mismo frame).

- Reagrupa átomos por (cadena, residuo) en orden canónico (AutoDock los
  intercala por torsión, lo que rompe PDBFixer).
- Elimina residuos que NO son aminoácidos estándar (RNA, aguas, ligandos),
  porque PDBFixer no puede parametrizarlos (p.ej. el RNA de 6G99 en FUS).
- Preserva las coordenadas del frame del docking.
"""
import collections

AA = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
      "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
      "HID", "HIE", "HIP", "CYX", "ASH", "GLH", "LYN"}

ORDER = {
    "N": 0, "CA": 1, "C": 2, "O": 3, "OXT": 4, "OT1": 4, "OT2": 5,
    "CB": 10, "CG": 11, "CG1": 12, "CG2": 13, "CD": 14, "CD1": 15, "CD2": 16,
    "CE": 17, "CE1": 18, "CE2": 19, "CE3": 20, "CZ": 21, "CZ2": 22, "CZ3": 23,
    "CH2": 24, "ND": 25, "ND1": 26, "ND2": 27, "NE": 28, "NE1": 29, "NE2": 30,
    "NH1": 31, "NH2": 32, "NZ": 33, "OD1": 34, "OD2": 35, "OE1": 36, "OE2": 37,
    "OG": 38, "OG1": 39, "OH": 40, "SG": 41, "SD": 42,
    "H": 50, "HA": 51, "HB2": 52, "HB3": 53, "HG": 54, "HD1": 55, "HD2": 56,
    "HE": 57, "HE1": 58, "HE2": 59, "HH": 60, "HH2": 61, "HZ": 62, "HZ2": 63,
}

def norm_name(n):
    return n.split("_")[0]

def reordenar(pqt, out):
    atoms = collections.defaultdict(list)  # (chain, resnum) -> list
    for line in open(pqt):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        if resname not in AA:
            continue
        try:
            resnum = int(line[22:26].strip())
        except ValueError:
            continue
        chain = line[21:22].strip() or " "
        name = norm_name(line[12:16].strip())
        atoms[(chain, resnum)].append((ORDER.get(name, 99), name, line))

    out_lines = []
    serial = 0
    for (chain, resnum) in sorted(atoms):
        for _, name, line in sorted(atoms[(chain, resnum)],
                                    key=lambda t: (t[0], t[1])):
            serial += 1
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            resname = line[17:20].strip()
            elem = line[76:78].strip() or name[0]
            out_lines.append(
                "ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
                % (serial, name, resname, chain, resnum, x, y, z, elem))
    with open(out, "w") as f:
        f.write("\n".join(out_lines) + "\nEND\n")
    print("%s: %d residuos, %d atomos -> %s"
          % (pqt, len(atoms), serial, out))

reordenar("rescoring_pkg2/receptores/TDP43.pdbqt", "receptor_TDP43_rrm1_canon.pdb")
reordenar("rescoring_pkg2/receptores/FUS.pdbqt", "receptor_FUS_canon.pdb")
reordenar("rescoring_pkg2/receptores/SOD1.pdbqt", "receptor_SOD1_canon.pdb")
