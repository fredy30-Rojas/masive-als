# -*- coding: utf-8 -*-
"""Aísla la fuente de variación en la preparación del receptor."""
import hashlib
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M

POSE = "/mnt/c/Users/Fredy/masive-als/gpu_dock/resultados_libreria/results_SOD1/CHEMBL4584906_out.pdbqt"
REC = "/mnt/c/Users/Fredy/masive-als/gpu_dock/SOD1.pdbqt"


def h(txt):
    return hashlib.md5(txt.encode()).hexdigest()[:12]


def hpos(pos):
    arr = np.asarray(pos.value_in_unit(pos.unit), dtype=float)
    return hashlib.md5(arr.round(4).tobytes()).hexdigest()[:12], arr.shape


lig_center = np.array([[a["x"], a["y"], a["z"]]
                       for a in M.parse_pdbqt_atoms(POSE)]).mean(axis=0)
print("centro del ligando:", " ".join("%.4f" % v for v in lig_center))
print("hash del PDBQT del receptor:",
      h(open(REC, "rb").read().hex()[:200000]))

for it in (1, 2, 3):
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    try:
        M.extract_receptor_pdb(REC, lig_center, tmp)
        pdb_txt = open(tmp).read()
        print("\n--- intento %d ---" % it)
        print("  PDB extraído -> hash %s | %d líneas | %d residuos"
              % (h(pdb_txt), pdb_txt.count("\n"),
                 len({l[17:26] for l in pdb_txt.splitlines()
                      if l.startswith(("ATOM", "HETATM"))})))
        rec_top, rec_pos = M.prepare_receptor(tmp)
        hh, shape = hpos(rec_pos)
        print("  tras PDBFixer -> %d átomos, hash de posiciones %s %s"
              % (rec_top.getNumAtoms(), hh, shape))
        resnames = "".join(r.name for r in rec_top.residues())
        print("  cadena de residuos (hash):", h(resnames), "| n residuos:",
              sum(1 for _ in rec_top.residues()))
        rec_off, rec_mols = M.off_topology_from_pdb(rec_top, rec_pos)
        print("  moléculas OpenFF: %d" % len(rec_mols))
        q = 0.0
        for m in rec_mols:
            if m.partial_charges is None:
                m.assign_partial_charges(partial_charge_method="mmff94")
            q += float(np.sum(np.asarray(m.partial_charges.magnitude)))
        print("  carga total OpenFF (mmff94): %.4f" % q)
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass
