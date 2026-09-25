import sys, tempfile, os
sys.path.insert(0, r"C:\Users\Fredy\masive-als\analysis\rescoring_local")
import mmgbsa_openff as M
import numpy as np
from openmm.app import PDBFile
from openff.toolkit import Topology as OFFTopology, ForceField

pose = r"C:\Users\Fredy\masive-als\gpu_dock\resultados_libreria\results_TDP43_v2\CHEMBL4571895_out.pdbqt"
lig_center = np.array([[a["x"], a["y"], a["z"]]
                       for a in M.parse_pdbqt_atoms(pose)]).mean(axis=0)
tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
try:
    M.extract_receptor_pdb(r"C:\Users\Fredy\masive-als\gpu_dock\TDP43_v2.pdbqt",
                           lig_center, tmp)
    rec_top, rec_pos = M.prepare_receptor(tmp)
    tmp2 = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    with open(tmp2, "w") as f:
        PDBFile.writeFile(rec_top, rec_pos, f)
    off = OFFTopology.from_pdb(tmp2)
    print("TOPOLOGY OK mols:", off.n_molecules, "atoms:", off.n_atoms)
    ff = ForceField(M.FF_PATH)
    sys2 = M.build_system(off, ff)
    print("SYSTEM OK forces:", sys2.getNumForces())
finally:
    for p in (tmp, tmp2 if "tmp2" in dir() else None):
        if p:
            try:
                os.remove(p)
            except OSError:
                pass