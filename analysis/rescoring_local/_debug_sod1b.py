# -*- coding: utf-8 -*-
import sys
sys.path.insert(0, r"C:\Users\Fredy\masive-als\analysis\rescoring_local")
import mmgbsa_openff as m
import tempfile, os
from collections import Counter, defaultdict

rec = r"C:\Users\Fredy\masive-als\gpu_dock\SOD1.pdbqt"
tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
m.extract_receptor_pdb(rec, None, tmp)

from pdbfixer import PDBFixer
fixer = PDBFixer(filename=tmp)
fixer.missingResidues = []
fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)
top = fixer.topology

res_atoms = defaultdict(list)
for atom in top.atoms():
    res_atoms[atom.residue.id].append(atom.name)
for rid in ["153", "154", "306", "307"]:
    names = res_atoms.get(rid, [])
    print("res", rid, "->", names)

bonds = []
for b in top.bonds():
    a1, a2 = b
    r1, r2 = a1.residue.id, a2.residue.id
    if {r1, r2} in ({"153", "154"}, {"306", "307"}):
        bonds.append((r1, a1.name, r2, a2.name))
print("enlaces entre 153/154 y 306/307:", bonds)
print("chains:", Counter(r.chain.id for r in top.residues()))
print("num residuos:", top.getNumResidues())
os.remove(tmp)
