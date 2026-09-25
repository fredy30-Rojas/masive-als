#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""De donde sale el ruido del MM-GBSA: ¿de la preparacion o de OpenMM?

Compara dos preparaciones (receptor y ligando) y dos minimizaciones desde el
mismo punto de partida, en el mismo proceso.
"""
import sys
import tempfile

import numpy as np
import pandas as pd

sys.path.insert(0, "/home/ubuntu/mmgbsa")

from rescoring_mmgbsa_robusto import (  # noqa: E402
    _recortar_pdb_canonical, ligand_center, prepare_receptor, prepare_ligand,
    build_system, minimize, potential_energy)

LIG = "ACT_isoproterenol"


def dif(a, b):
    a = np.array(a.value_in_unit(a.unit))
    b = np.array(b.value_in_unit(b.unit))
    if a.shape != b.shape:
        return f"formas distintas {a.shape} vs {b.shape}"
    return f"max|d|={np.abs(a - b).max():.2e}"


def main():
    fila = pd.read_csv("estratos/candidatos.csv")
    fila = fila[fila.ligand == LIG].iloc[0]
    c = ligand_center(fila["pose_pdbqt"])

    from openmm.app import ForceField

    tops, poss = [], []
    for rep in (1, 2):
        tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
        _recortar_pdb_canonical("estratos/receptor.pdb", c, tmp, radio=14.0)
        t, p = prepare_receptor(tmp)
        tops.append(t)
        poss.append(p)
    print(f"receptor: {len(list(tops[0].atoms()))} vs {len(list(tops[1].atoms()))} "
          f"atomos -> {dif(poss[0], poss[1])}")

    ltops, lposs = [], []
    for rep in (1, 2):
        ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
        t, p = prepare_ligand(fila["smiles"], fila["pose_pdbqt"], ff)
        ltops.append(t)
        lposs.append({rep: (ff, t, p)})
    print(f"ligando: {len(list(ltops[0].atoms()))} vs {len(list(ltops[1].atoms()))} "
          f"atomos -> {dif(lposs[0][1][2], lposs[1][2][2])}")

    # --- minimizacion: mismo sistema, mismo punto de partida, dos veces
    ff1, lt1, lp1 = lposs[0][1]
    from openmm.app import Modeller
    mod = Modeller(tops[0], poss[0])
    mod.add(lt1, lp1)
    print(f"complejo: {len(list(mod.topology.atoms()))} atomos")
    sys1 = build_system(mod.topology, ff1)
    for rep in (1, 2):
        start = mod.positions
        e0 = potential_energy(sys1, mod.topology, start)
        fin = minimize(sys1, mod.topology, start)
        e1 = potential_energy(sys1, mod.topology, fin)
        print(f"  rep{rep}: E_inicial={e0:.3f}  E_minimizada={e1:.3f}  "
              f"({dif(start, fin)})")


if __name__ == "__main__":
    main()
