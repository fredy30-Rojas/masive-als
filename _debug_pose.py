import sys
sys.path.insert(0, r"C:\Users\Fredy\masive-als\analysis\rescoring_local")
import mmgbsa_openff as M

pose = r"C:\Users\Fredy\masive-als\gpu_dock\resultados_libreria\results_TDP43_v2\CHEMBL4571895_out.pdbqt"
atoms = M.parse_pdbqt_atoms(pose)
print("atoms parseados:", len(atoms))
from collections import Counter
print(Counter(a["elem"] for a in atoms))
mol = M.pose_to_rdmol(pose)
print("mol atoms:", mol.GetNumAtoms())
print("names:", [a["name"] for a in atoms[:6]])

# SMILES
from rdkit import Chem
smi = "O=C(COc1ccc(-c2nc3ccccc3c(=O)[nH]2)cc1)Nc1cccc2ccccc12"
rd = Chem.MolFromSmiles(smi)
print("smiles atoms:", rd.GetNumAtoms())

# MCS
from rdkit.Chem import rdFMCS
mcs = rdFMCS.FindMCS([rd, mol], bondCompare=rdFMCS.BondCompare.CompareAny,
                     atomCompare=rdFMCS.AtomCompare.CompareElements,
                     ringMatchesRingOnly=True, completeRingsOnly=True, timeout=15)
print("MCS:", mcs.numAtoms, "/", rd.GetNumAtoms())
print("MCS smarts:", mcs.smartsString)