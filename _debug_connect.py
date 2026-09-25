import sys
sys.path.insert(0, r"C:\Users\Fredy\masive-als\analysis\rescoring_local")
import mmgbsa_openff as M
from rdkit import Chem
from rdkit.Chem import rdFMCS

pose = r"C:\Users\Fredy\masive-als\gpu_dock\resultados_libreria\results_TDP43_v2\CHEMBL4571895_out.pdbqt"
docked = M.pose_to_rdmol(pose)
smi = "O=C(COc1ccc(-c2nc3ccccc3c(=O)[nH]2)cc1)Nc1cccc2ccccc12"
rd = Chem.MolFromSmiles(smi)

for label, kw in [
    ("estricto", dict(ringMatchesRingOnly=True, completeRingsOnly=True)),
    ("solo anillos", dict(ringMatchesRingOnly=True, completeRingsOnly=False)),
    ("relajado", dict(ringMatchesRingOnly=False, completeRingsOnly=False)),
]:
    mcs = rdFMCS.FindMCS([rd, docked], bondCompare=rdFMCS.BondCompare.CompareAny,
                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                         timeout=20, **kw)
    print(label, "->", mcs.numAtoms, "/", rd.GetNumAtoms())