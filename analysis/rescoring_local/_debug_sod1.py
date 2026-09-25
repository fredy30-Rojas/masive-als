# -*- coding: utf-8 -*-
import sys, traceback
sys.path.insert(0, r"C:\Users\Fredy\masive-als\analysis\rescoring_local")
import mmgbsa_openff as m

pose = r"C:\Users\Fredy\masive-als\gpu_dock\resultados_libreria\results_SOD1\CHEMBL4584906_out.pdbqt"
rec = r"C:\Users\Fredy\masive-als\gpu_dock\SOD1.pdbqt"
smi = "N=C1N/C(=N/NC(=O)c2ccc(CN3C(=O)c4cccc5cccc3c45)cc2)c2ccccc21"

import tempfile, os
print("1) extrayendo receptor PDB...", flush=True)
atoms = m.parse_pdbqt_atoms(rec)
from collections import Counter
print("   residuos en receptor:", Counter(a["resname"] for a in atoms).most_common(10), flush=True)
tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
m.extract_receptor_pdb(rec, None, tmp)
print("   pdb escrito:", tmp, flush=True)

print("2) pdbfixer...", flush=True)
try:
    rec_top, rec_pos = m.prepare_receptor(tmp)
    print("   ok, atomos:", rec_top.getNumAtoms(), flush=True)
except Exception:
    traceback.print_exc()
    sys.exit(1)

print("3) off topology from pdb...", flush=True)
try:
    rec_off, rec_mols = m.off_topology_from_pdb(rec_top, rec_pos)
    print("   ok, moleculas:", len(rec_mols), flush=True)
except Exception:
    traceback.print_exc()
    sys.exit(1)

print("4) ligando...", flush=True)
from openff.toolkit import ForceField
ff = OFFForceField(m.FF_PATH)
try:
    lig_top, lig_pos, offmol = m.prepare_ligand(smi, pose, ff)
    print("   ok", flush=True)
except Exception:
    traceback.print_exc()
    sys.exit(1)

print("5) build system complejo (aqui fallaba)...", flush=True)
from openmm.app import Modeller
complex_mod = Modeller(rec_top, rec_pos)
complex_mod.add(lig_top, lig_pos)
from openff.toolkit import Topology as OFFTopology
complex_off = OFFTopology()
for mol in rec_mols:
    complex_off.add_molecule(mol)
complex_off.add_molecule(offmol)
try:
    system = m.build_system(complex_off, ff)
    print("   OK sistema complejo creado", flush=True)
except Exception:
    traceback.print_exc()
    sys.exit(1)
print("TODO OK", flush=True)
os.remove(tmp)
