#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Paso a paso del ligando que falla con el residuo THR, con traza completa."""
import sys
import tempfile
import traceback

sys.path.insert(0, "/home/ubuntu/mmgbsa")

import pandas as pd  # noqa: E402

from rescoring_mmgbsa_robusto import (  # noqa: E402
    _recortar_pdb_canonical, ligand_center, prepare_receptor,
    prepare_ligand, build_system, minimize, potential_energy)

LIG = sys.argv[1] if len(sys.argv) > 1 else "DEC_CHEMBL1414576"


def lista_residuos(top):
    vistos = []
    for r in top.residues():
        vistos.append((r.index + 1, r.chain.id if r.chain else "?",
                       r.id, r.name, len(list(r.atoms()))))
    return vistos


def main():
    rows = pd.read_csv("estratos/candidatos.csv")
    fila = rows[rows.ligand == LIG].iloc[0]
    pose, smiles = fila["pose_pdbqt"], fila["smiles"]
    print(f"{LIG}\n  smiles={smiles}\n  pose={pose}")

    c = ligand_center(pose)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    _recortar_pdb_canonical("estratos/receptor.pdb", c, tmp, radio=14.0)
    print(f"  centro=({c[0]:.2f},{c[1]:.2f},{c[2]:.2f})  recorte={tmp}")
    rec_top, rec_pos = prepare_receptor(tmp)
    res = lista_residuos(rec_top)
    print(f"  receptor: {len(res)} residuos")
    for i, ch, rid, name, n in res:
        marca = "  <-- nombrado en el error" if i in (35, 47) else ""
        if i > 30 or name == "THR":
            print(f"    #{i:3d} {ch}{rid} {name} ({n} atomos){marca}")

    from openmm.app import ForceField
    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")

    print("  --- prepare_ligand")
    try:
        lig_top, lig_pos = prepare_ligand(smiles, pose, ff)
        print(f"      OK ligando: {lig_top.getNumAtoms()} atomos")
    except Exception:
        traceback.print_exc(); return

    print("  --- complejo")
    from openmm.app import Modeller
    try:
        mod = Modeller(rec_top, rec_pos)
        mod.add(lig_top, lig_pos)
        sysc = build_system(mod.topology, ff)
        mp = minimize(sysc, mod.topology, mod.positions)
        ec = potential_energy(sysc, mod.topology, mp)
        print(f"      OK complejo: E={ec:.2f}")
    except Exception:
        traceback.print_exc(); return


if __name__ == "__main__":
    main()
