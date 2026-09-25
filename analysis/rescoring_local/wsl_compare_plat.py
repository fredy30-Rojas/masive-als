# -*- coding: utf-8 -*-
"""Compara CUDA (determinista) con Reference usando un receptor FIJADO.

Receptor preparado una sola vez → todas las corridas parten de exactamente
los mismos átomos y coordenadas, así que cualquier diferencia es del cálculo.
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M
import mmgbsa_openff_gb as G
import numpy as np

POSE = "/mnt/c/Users/Fredy/masive-als/gpu_dock/resultados_libreria/results_SOD1/CHEMBL4584906_out.pdbqt"
REC = "/mnt/c/Users/Fredy/masive-als/gpu_dock/SOD1.pdbqt"
SM = "N=C1N/C(=N/NC(=O)c2ccc(CN3C(=O)c4cccc5cccc3c45)cc2)c2ccccc21"
FIJO = "/tmp/SOD1_fijo.pdb"

lig_center = np.array([[a["x"], a["y"], a["z"]]
                       for a in M.parse_pdbqt_atoms(POSE)]).mean(axis=0)
n = G.preparar_receptor_fijo(REC, FIJO, cutoff=12.0, lig_center=lig_center)
import hashlib
h = hashlib.md5(open(FIJO, "rb").read()).hexdigest()[:12]
print("receptor fijado: %d átomos | hash %s | %s" % (n, h, FIJO))
print("tolerancia %s kJ/mol/nm | maxiter %s" % (G.TOLERANCE, G.MAX_ITER))

resultados = {}
for plat in ("CUDA", "CUDA", "Reference", "Reference"):
    os.environ["MMGBSA_FINAL_PLATFORM"] = plat
    t0 = time.time()
    r = G.rescore_one_fijo(POSE, FIJO, SM)
    dt = time.time() - t0
    print("  %-9s -> dG=%8.2f  (E_c=%.1f E_r=%.1f E_l=%.1f)  %5.1f s"
          % (plat, r["mmgbsa_dG"], r["e_complex"], r["e_receptor"],
             r["e_ligand"], dt))
    resultados.setdefault(plat, []).append(r["mmgbsa_dG"])

print()
for plat, vals in resultados.items():
    print("%-9s: %s | dispersión %.2f kcal/mol"
          % (plat, " ".join("%.2f" % v for v in vals),
             max(vals) - min(vals)))
d_cuda = sum(resultados["CUDA"]) / len(resultados["CUDA"])
d_ref = sum(resultados["Reference"]) / len(resultados["Reference"])
print("\ndiferencia de medias CUDA vs Reference: %+.2f kcal/mol" % (d_cuda - d_ref))
