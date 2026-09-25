# -*- coding: utf-8 -*-
"""¿De dónde viene la variación? Compara la energía del complejo ANTES de
minimizar (depende solo de la entrada) con la energía DESPUÉS (depende del
minimizador y de la plataforma)."""
import hashlib
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M
import mmgbsa_openff_gb as G

POSE = "/mnt/c/Users/Fredy/masive-als/gpu_dock/resultados_libreria/results_SOD1/CHEMBL4584906_out.pdbqt"
REC = "/mnt/c/Users/Fredy/masive-als/gpu_dock/SOD1.pdbqt"
SM = "N=C1N/C(=N/NC(=O)c2ccc(CN3C(=O)c4cccc5cccc3c45)cc2)c2ccccc21"

print("FF_PATH:", M.FF_PATH)
print("plataforma por defecto:", G.DEFAULT_PLATFORM)
print("tolerancia:", G.TOLERANCE, "maxiter:", G.MAX_ITER)

from openff.toolkit import ForceField as OFFForceField
from openmm.app import Modeller

ff = OFFForceField(M.FF_PATH)

for intento in (1, 2):
    lig_top, lig_pos, offmol = G.prepare_ligand_pose(SM, POSE, ff)
    import numpy as _np
    _arr = _np.asarray(lig_pos.value_in_unit(lig_pos.unit), dtype=float)
    lighash = hashlib.md5(_arr.round(3).tobytes()).hexdigest()[:12]
    print("  forma de las coordenadas:", _arr.shape)
    print("\n--- intento %d ---" % intento)
    print("  hash de coordenadas del ligando:", lighash)
    print("  atomos ligando:", lig_top.getNumAtoms())

    import numpy as np
    lig_center = np.array([[a["x"], a["y"], a["z"]]
                           for a in M.parse_pdbqt_atoms(POSE)]).mean(axis=0)
    import tempfile
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    M.extract_receptor_pdb(REC, lig_center, tmp)
    rec_top, rec_pos = M.prepare_receptor(tmp)
    os.remove(tmp)
    print("  atomos receptor:", rec_top.getNumAtoms())

    from openff.toolkit import Topology as OFFTopology
    rec_off, rec_mols = M.off_topology_from_pdb(rec_top, rec_pos)
    complex_mod = Modeller(rec_top, rec_pos)
    complex_mod.add(lig_top, lig_pos)
    complex_off = OFFTopology()
    for m in rec_mols:
        complex_off.add_molecule(m)
    complex_off.add_molecule(offmol)
    sysx = M.build_system(complex_off, ff)
    G.add_gb_obc2(sysx, complex_mod.topology)
    e0 = M.potential_energy(sysx, complex_mod.topology, complex_mod.positions)
    print("  E complejo ANTES de minimizar: %.2f" % e0)

    for v in (1, 2):
        min_pos = G.minimize_on(sysx, complex_mod.topology,
                                complex_mod.positions,
                                G.DEFAULT_PLATFORM, G.MAX_ITER, G.TOLERANCE)
        e1 = M.potential_energy(sysx, complex_mod.topology, min_pos)
        print("     minimizacion %d -> E = %.2f" % (v, e1))
