# -*- coding: utf-8 -*-
"""Aisla el origen de la dispersion residual con receptor congelado.

Repite TODAS las etapas del pipeline dos veces (misma sesion y sesiones
distintas) y hashea cada artefacto intermedio: coordenadas del ligando
tras anadir H, cargas parciales, energia inicial del complejo, sistema
(por sus parametros) y energia minimizada.
"""
import hashlib
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M        # noqa: E402
import mmgbsa_openff_gb as G     # noqa: E402
from openmm.app import Modeller, PDBFile  # noqa: E402
from openff.toolkit import ForceField as OFFForceField, Topology as OFFTopology  # noqa: E402


def h(arr):
    return hashlib.md5(np.asarray(arr, dtype=np.float64).tobytes()).hexdigest()[:10]


def hu(q):
    return h(q.value_in_unit(q.unit))


def etapas(pose, rec_pdb, smiles, etiqueta):
    ff = OFFForceField(M.FF_PATH)
    out = {"etiqueta": etiqueta}

    lig_top, lig_pos, offmol = G.prepare_ligand_pose(smiles, pose, ff)
    out["h_lig_pos"] = hu(lig_pos)
    out["n_lig"] = lig_top.getNumAtoms()
    ch = offmol.partial_charges
    out["h_lig_cargas"] = h(ch.magnitude) if ch is not None else "None"

    rec = PDBFile(rec_pdb)
    rec_top, rec_pos = rec.topology, rec.positions
    n_rec = rec_top.getNumAtoms()
    out["n_rec"] = n_rec
    out["h_rec_pos"] = hu(rec_pos)

    rec_off, rec_mols = M.off_topology_from_pdb(rec_top, rec_pos)
    out["n_rec_mols"] = len(rec_mols)

    complex_mod = Modeller(rec_top, rec_pos)
    complex_mod.add(lig_top, lig_pos)
    complex_off = OFFTopology()
    for m in rec_mols:
        complex_off.add_molecule(m)
    complex_off.add_molecule(offmol)

    sys_c = M.build_system(complex_off, ff)
    G.add_gb_obc2(sys_c, complex_mod.topology)
    pos0 = complex_mod.positions
    out["h_pos0"] = hu(pos0)
    out["e0_complejo"] = round(M.potential_energy(sys_c, complex_mod.topology, pos0), 3)

    mp = G.minimize_best_of(sys_c, complex_mod.topology, pos0)
    out["h_min"] = hu(mp)
    out["e_min"] = round(M.potential_energy(sys_c, complex_mod.topology, mp), 3)

    sys_r = M.build_system(rec_off, ff)
    G.add_gb_obc2(sys_r, rec_top)
    e_r = M.potential_energy(sys_r, rec_top, mp[:n_rec])
    out["e_receptor"] = round(e_r, 3)

    lig_off = offmol.to_topology()
    sys_l = M.build_system(lig_off, ff)
    G.add_gb_obc2(sys_l, lig_top)
    e_l = M.potential_energy(sys_l, lig_top, mp[n_rec:])
    out["e_ligando"] = round(e_l, 3)
    out["dG"] = round(out["e_min"] - e_r - e_l, 3)
    return out


if __name__ == "__main__":
    BASE = "/mnt/c/Users/Fredy/masive-als"
    LOC = os.path.join(BASE, "analysis", "rescoring_local")
    rec = os.path.join(LOC, "receptores_fijos", "SOD1_fijo.pdb")
    pose = os.path.join(BASE, "gpu_dock", "resultados_libreria", "results_SOD1",
                        "CHEMBL4584906_out.pdbqt")
    sm = ("N=C1N/C(=N/NC(=O)c2ccc(CN3C(=O)c4cccc5cccc3c45)cc2)c2ccccc21")
    rep = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    res = [etapas(pose, rec, sm, "r%d" % i) for i in range(rep)]
    for r in res:
        print(json.dumps(r, ensure_ascii=False), flush=True)
    claves = [k for k in res[0] if k.startswith("h_")]
    print("\nclaves de hash que cambian:")
    for k in claves:
        vals = {r[k] for r in res}
        print("  %-14s %s" % (k, "IDENTICO" if len(vals) == 1 else "DIFIERE %s" % vals))
