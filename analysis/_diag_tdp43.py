# -*- coding: utf-8 -*-
"""Diagnóstico: dónde se cuelga el rescoring MM-GBSA de TDP43."""
import csv, os, sys, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import rescoring_mmgbsa_robusto as R

rows = list(csv.DictReader(open("rescoring_pkg/candidatos_42.csv")))
r = rows[6]  # caso TDP43
print("target:", r["target"], "lig:", r["ligand"], flush=True)

t0 = time.time()
c = R.ligand_center(r["pose_pdbqt"])
print("[%.1fs] centro %s" % (time.time() - t0, c), flush=True)

t0 = time.time()
tmp = R.tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
R.extract_receptor_pdb("rescoring_pkg2/receptores/TDP43.pdbqt", c, tmp)
n = sum(1 for l in open(tmp) if l.startswith("ATOM"))
print("[%.1fs] receptor extraido: %d atomos (%s)"
      % (time.time() - t0, n, tmp), flush=True)

from openmm.app import ForceField
ff = ForceField("amber14-all.xml", "implicit/obc2.xml")

t0 = time.time()
rec_top, rec_pos = R.prepare_receptor(tmp)
print("[%.1fs] receptor preparado OK (%d atomos)"
      % (time.time() - t0, rec_top.getNumAtoms()), flush=True)

t0 = time.time()
lig_top, lig_pos = R.prepare_ligand(r["smiles"], r["pose_pdbqt"], ff)
print("[%.1fs] ligando preparado OK" % (time.time() - t0,), flush=True)

t0 = time.time()
dg, ec, er, el = R.rescore_one(rec_top, rec_pos, r["smiles"], r["pose_pdbqt"], ff)
print("[%.1fs] dG=%.2f e_complex=%.2f e_rec=%.2f e_lig=%.2f"
      % (time.time() - t0, dg, ec, er, el), flush=True)
