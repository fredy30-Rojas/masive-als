# -*- coding: utf-8 -*-
"""¿Dónde se va el tiempo? Reparte el coste por fases."""
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M
import mmgbsa_openff_gb as G

POSE = "/mnt/c/Users/Fredy/masive-als/gpu_dock/resultados_libreria/results_SOD1/CHEMBL4584906_out.pdbqt"
REC = "/mnt/c/Users/Fredy/masive-als/gpu_dock/SOD1.pdbqt"
SM = "N=C1N/C(=N/NC(=O)c2ccc(CN3C(=O)c4cccc5cccc3c45)cc2)c2ccccc21"
FIJO = "/tmp/SOD1_fijo.pdb"

from openff.toolkit import ForceField as OFFForceField, Topology as OFFTopology
from openmm.app import Modeller, PDBFile
import tempfile

T = {}


def crono(nombre, f):
    t0 = time.time()
    r = f()
    T[nombre] = T.get(nombre, 0) + time.time() - t0
    return r


ff = OFFForceField(M.FF_PATH)
os.environ.setdefault("MMGBSA_FINAL_PLATFORM", "CUDA")

lig_top, lig_pos, offmol = crono(
    "1 ligando (geometría + cargas)",
    lambda: G.prepare_ligand_pose(SM, POSE, ff))
lig_off = offmol.to_topology()

lig_center = np.array([[a["x"], a["y"], a["z"]]
                       for a in M.parse_pdbqt_atoms(POSE)]).mean(axis=0)
if not os.path.exists(FIJO):
    crono("0 receptor fijo (solo 1 vez)",
          lambda: G.preparar_receptor_fijo(REC, FIJO, 12.0, lig_center))

rec_pdb = PDBFile(FIJO)
rec_top, rec_pos = rec_pdb.topology, rec_pdb.positions
n_rec = rec_top.getNumAtoms()
rec_off, rec_mols = crono("2 topología OpenFF del receptor",
                          lambda: M.off_topology_from_pdb(rec_top, rec_pos))

complex_mod = Modeller(rec_top, rec_pos)
complex_mod.add(lig_top, lig_pos)
complex_off = OFFTopology()
for m in rec_mols:
    complex_off.add_molecule(m)
complex_off.add_molecule(offmol)

complex_sys = crono("3 parametrizar COMPLEJO (SMIRNOFF)",
                    lambda: M.build_system(complex_off, ff))
crono("4 GB complejo", lambda: G.add_gb_obc2(complex_sys, complex_mod.topology))
min_pos = crono("5 minimizar complejo (CUDA)",
                lambda: G.minimize_best_of(complex_sys, complex_mod.topology,
                                           complex_mod.positions))
crono("6 energía complejo",
      lambda: M.potential_energy(complex_sys, complex_mod.topology, min_pos))
rec_sys = crono("7 parametrizar RECEPTOR (SMIRNOFF)",
                lambda: M.build_system(rec_off, ff))
crono("8 GB receptor", lambda: G.add_gb_obc2(rec_sys, rec_top))
crono("9 energía receptor",
      lambda: M.potential_energy(rec_sys, rec_top, min_pos[:n_rec]))
lig_sys = crono("10 parametrizar LIGANDO (SMIRNOFF)",
                lambda: M.build_system(lig_off, ff))
crono("11 GB ligando", lambda: G.add_gb_obc2(lig_sys, lig_top))
crono("12 energía ligando",
      lambda: M.potential_energy(lig_sys, lig_top, min_pos[n_rec:]))

total = sum(T.values())
print("\n=== reparto del tiempo (total %.1f s) ===" % total)
for k in sorted(T, key=lambda x: -T[x]):
    print("  %-38s %7.2f s  %5.1f%%" % (k, T[k], 100 * T[k] / total))
