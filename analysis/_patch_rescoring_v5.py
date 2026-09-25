# -*- coding: utf-8 -*-
"""Parchea rescoring_mmgbsa_robusto.py (v5):
- Si el receptor es un PDB canónico (.pdb, nombres/orden estándar), lo usa
  directo con PDBFixer, recortado al entorno del ligando (radio 10 A) para
  acelerar la minimización y evitar la geometría rota del pdbqt.
- Si es .pdbqt, mantiene el flujo anterior (extract_receptor_pdb).
"""
import io

P = "/home/ubuntu/mmgbsa/rescoring_mmgbsa_robusto.py"
with io.open(P, "r", encoding="utf-8") as f:
    txt = f.read()

# 1) nueva función: recortar PDB canónico al entorno del ligando
helper = '''

def _recortar_pdb_canonical(pdb_path, lig_center, out_pdb, radio=10.0):
    """Recorta un PDB canonico a los residuos con algun atomo a < radio del ligando."""
    import math
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
            lines.append(line[:6] + "%5d" % serial + line[11:].rstrip())
    if not lines:
        raise RuntimeError("recorte vacío (radio %.1f A)" % radio)
    with open(out_pdb, "w") as f:
        f.write("\\n".join(lines) + "\\nEND\\n")
    return out_pdb

'''

# insertar la funcion despues de extract_receptor_pdb
anchor = "def prepare_receptor(pdb_path):"
assert anchor in txt, "no encontrado anchor prepare_receptor"
txt = txt.replace(anchor, helper + anchor, 1)

# 2) en compute_row: detectar .pdb canonico vs .pdbqt
old = '''        c = ligand_center(pose)
        tmp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        extract_receptor_pdb(receptores[tgt], c, tmp_pdb)
        rec_top, rec_pos = prepare_receptor(tmp_pdb)'''
new = '''        c = ligand_center(pose)
        tmp_pdb = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        rpath = receptores[tgt]
        if rpath.lower().endswith(".pdb") and not rpath.lower().endswith(".pdbqt"):
            _recortar_pdb_canonical(rpath, c, tmp_pdb, radio=10.0)
        else:
            extract_receptor_pdb(rpath, c, tmp_pdb)
        rec_top, rec_pos = prepare_receptor(tmp_pdb)'''
assert old in txt, "no encontrado bloque compute_row"
txt = txt.replace(old, new, 1)

with io.open(P, "w", encoding="utf-8") as f:
    f.write(txt)
print("parche v5 aplicado")
