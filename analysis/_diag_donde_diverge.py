#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿En que paso de la fila completa se separan dos corridas de lo mismo?

Repite el mismo camino que `rescore_one` pero imprimiendo los pasos intermedios:
tamaños, donde queda el ligando respecto a la pose, energía del complejo antes y
despues de minimizar, cuanto se mueve el receptor, y los tres terminos.
"""
import os
import sys
import tempfile

import numpy as np

sys.path.insert(0, "/home/ubuntu/mmgbsa")

from rescoring_mmgbsa_robusto import (  # noqa: E402
    _recortar_pdb_canonical, build_system, ligand_center, minimize,
    potential_energy, prepare_ligand, prepare_receptor)

import pandas as pd  # noqa: E402

LIG = "ACT_isoproterenol"


def main():
    from openmm.app import ForceField, Modeller
    from openmm import unit

    fila = pd.read_csv("estratos/candidatos.csv")
    fila = fila[fila.ligand == LIG].iloc[0]
    c = ligand_center(fila["pose_pdbqt"])

    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    _recortar_pdb_canonical("estratos/receptor.pdb", c, tmp, radio=14.0)
    r_top, r_pos = prepare_receptor(tmp)
    n_rec = r_top.getNumAtoms()

    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
    l_top, l_pos = prepare_ligand(fila["smiles"], fila["pose_pdbqt"], ff)
    n_lig = l_top.getNumAtoms()

    rpos = np.array(r_pos.value_in_unit(unit.angstrom))
    lpos = np.array(l_pos.value_in_unit(unit.angstrom))
    import pose_pdbqt
    p = pose_pdbqt.leer_pose(fila["pose_pdbqt"], modelo=1)
    pose = np.array([[a[0], a[1], a[2]] for a in p["atoms"].values()])
    # distancia minima ligando-receptor: si el ligando esta metido en el bolsillo
    # son ~2-3 A; si quedo fuera del recorte, sale grande
    dmin = float(np.linalg.norm(
        lpos[:, None, :] - rpos[None, :, :], axis=2).min())
    print(f"pid={os.getpid()} receptor={n_rec} ligando={n_lig} "
          f"pose_modelo1={len(pose)} menor_distancia_al_receptor={dmin:.2f} A",
          flush=True)
    print(f"  centro receptor {np.round(rpos.mean(0), 3)} | "
          f"ligando {np.round(lpos.mean(0), 3)} | "
          f"pose {np.round(pose.mean(0), 3)} | "
          f"desplazamiento ligando-pose "
          f"{np.linalg.norm(lpos.mean(0) - pose.mean(0)):.3f} A", flush=True)

    mod = Modeller(r_top, r_pos)
    mod.add(l_top, l_pos)
    complejo = build_system(mod.topology, ff)
    e0 = potential_energy(complejo, mod.topology, mod.positions)
    min_pos = minimize(complejo, mod.topology, mod.positions)
    e1 = potential_energy(complejo, mod.topology, min_pos)

    mpos = np.array(min_pos.value_in_unit(unit.angstrom))
    mov_rec = np.linalg.norm(mpos[:n_rec] - rpos, axis=1)
    mov_lig = np.linalg.norm(mpos[n_rec:] - lpos, axis=1)
    print(f"  E complejo inicial {e0:.2f} -> minimizado {e1:.2f}", flush=True)
    print(f"  movido al minimizar: receptor max {mov_rec.max():.3f} A "
          f"(media {mov_rec.mean():.4f}) | ligando max {mov_lig.max():.3f} A "
          f"(media {mov_lig.mean():.4f})", flush=True)

    e_rec = potential_energy(build_system(r_top, ff), r_top, min_pos[:n_rec])
    e_lig = potential_energy(build_system(l_top, ff), l_top, min_pos[n_rec:])
    print(f"  E_rec={e_rec:.2f} E_lig={e_lig:.2f} dG={e1 - e_rec - e_lig:.2f}",
          flush=True)


if __name__ == "__main__":
    main()
