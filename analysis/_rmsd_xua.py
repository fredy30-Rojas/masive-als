# -*- coding: utf-8 -*-
"""RMSD ligando nativo: pose dockeada (Vina) vs pose cristalográfica (4IUF).

El RMSD se calcula sobre los átomos pesados comunes tras superponer por
correspondencia de nombres de átomo (mismo frame de 4IUF). Como la pose de
Vina sale en el mismo frame del receptor (TDP43.pdbqt == 4IUF), se pueden
comparar directamente las coordenadas.
"""
import math, sys

def parse(path, heavy=True):
    atoms = []
    for line in open(path):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16].strip()
        if heavy and (name.startswith("H") and name[1:2].isdigit() or
                      name in ("H", "HA", "HB", "HG", "HD", "HE", "HZ",
                               "H1", "H2", "H3", "HO", "HXT", "H5T", "H3T")):
            continue
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
        elem = line[76:78].strip() or name[0]
        atoms.append((name, elem, x, y, z))
    return atoms

def parse_pdbqt(path):
    """Pose de Vina: MULTIMODEL. Devuelve lista de conformeros."""
    modes, cur = [], []
    for line in open(path):
        if line.startswith("MODEL"):
            cur = []
        elif line.startswith(("ATOM", "HETATM")):
            name = line[12:16].strip().split("_")[0]
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            cur.append((name, x, y, z))
        elif line.startswith("ENDMDL"):
            if cur:
                modes.append(cur)
    return modes

# pose cristalografica (solo atomos pesados)
crist = [(n, x, y, z) for (n, e, x, y, z) in parse("_xua_nativo.pdb")
         if e != "H"]
crist_names = set(n for (n, x, y, z) in crist)
print("atomos pesados cristalograficos:", len(crist))

for mi, mode in enumerate(parse_pdbqt("_xua_redock_out.pdbqt"), 1):
    dock = [(n, x, y, z) for (n, x, y, z) in mode
            if n in crist_names]
    if not dock:
        continue
    # pares comunes por nombre de atomo (orden estable)
    byname = {}
    for n, x, y, z in crist:
        byname.setdefault(n, (x, y, z))
    pairs = [(byname[n], (x, y, z)) for (n, x, y, z) in dock if n in byname]
    if not pairs:
        continue
    rmsd = math.sqrt(sum(
        (a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2 for a, b in pairs
    ) / len(pairs))
    print("modo %d: %d atomos comunes, RMSD = %.2f A" % (mi, len(pairs), rmsd))
