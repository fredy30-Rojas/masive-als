# -*- coding: utf-8 -*-
"""Compromiso tolerancia / tiempo / valor, con receptor FIJO y fuerzas
deterministas (así cada diferencia es del criterio de convergencia)."""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M
import mmgbsa_openff_gb as G

POSE = "/mnt/c/Users/Fredy/masive-als/gpu_dock/resultados_libreria/results_SOD1/CHEMBL4584906_out.pdbqt"
SM = "N=C1N/C(=N/NC(=O)c2ccc(CN3C(=O)c4cccc5cccc3c45)cc2)c2ccccc21"
FIJO = "/tmp/SOD1_fijo.pdb"

import numpy as _np
REC = "/mnt/c/Users/Fredy/masive-als/gpu_dock/SOD1.pdbqt"
if not os.path.exists(FIJO):
    c = _np.array([[a["x"], a["y"], a["z"]]
                   for a in M.parse_pdbqt_atoms(POSE)]).mean(axis=0)
    G.preparar_receptor_fijo(REC, FIJO, 12.0, c)

print("plataforma: %s | receptor fijo: %s"
      % (G.DEFAULT_PLATFORM, os.path.exists(FIJO)))
print()
print("%10s %8s %10s   %s" % ("tol", "tiempo", "dG", "E_complejo"))
base = None
for tol in (1.0, 2.0, 5.0, 10.0, 50.0):
    G.TOLERANCE = tol
    t0 = time.time()
    r = G.rescore_one_fijo(POSE, FIJO, SM)
    dt = time.time() - t0
    if base is None:
        base = r["mmgbsa_dG"]
    print("%10.1f %7.1fs %10.2f   %.1f   (Δ vs tol=1: %+.2f)"
          % (tol, dt, r["mmgbsa_dG"], r["e_complex"], r["mmgbsa_dG"] - base))
