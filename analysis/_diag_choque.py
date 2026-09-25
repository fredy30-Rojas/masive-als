#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿Con que choca la pose: con que residuo, a que distancia, y es hidrogeno?

Recorta y prepara el receptor como lo hace el rescoring, prepara el ligando y
saca las distancias minimas atomo-atomo (pesados y con hidrogenos), diciendo
residuo, cadena y elemento de los peores choques.
"""
import sys
import tempfile

import numpy as np

sys.path.insert(0, "/home/ubuntu/mmgbsa")

from rescoring_mmgbsa_robusto import (  # noqa: E402
    _recortar_pdb_canonical, ligand_center, prepare_ligand, prepare_receptor)

import pandas as pd  # noqa: E402
import pose_pdbqt  # noqa: E402


def main():
    from openmm.app import ForceField
    from openmm import unit

    fila = pd.read_csv("estratos/candidatos.csv")
    fila = fila[fila.ligand == "ACT_isoproterenol"].iloc[0]
    c = ligand_center(fila["pose_pdbqt"])
    print("centro del recorte (pose 1):", np.round(c, 3), flush=True)

    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    _recortar_pdb_canonical("estratos/receptor.pdb", c, tmp, radio=14.0)
    print("lineas del recorte:", sum(1 for _ in open(tmp)), flush=True)
    r_top, r_pos = prepare_receptor(tmp)
    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
    l_top, l_pos = prepare_ligand(fila["smiles"], fila["pose_pdbqt"], ff)

    rat = list(r_top.atoms())
    lat = list(l_top.atoms())
    rp = np.array(r_pos.value_in_unit(unit.angstrom))
    lp = np.array(l_pos.value_in_unit(unit.angstrom))
    print(f"receptor {len(rat)} atomos | ligando {len(lat)} atomos", flush=True)

    d = np.linalg.norm(lp[:, None, :] - rp[None, :, :], axis=2)
    pares = [(d[i, j], i, j) for i in range(len(lat)) for j in range(len(rat))]
    pares.sort()
    for dist, i, j in pares[:10]:
        ri, rj = rat[j].residue, rat[j]
        print(f"  {dist:.2f} A  ligando {lat[i].name}({lat[i].element.symbol}) "
              f"- receptor {ri.name}{ri.id} cadena {ri.chain.id} "
              f"{rj.name}({rj.element.symbol})", flush=True)

    # solo atomos pesados
    ih = [i for i, a in enumerate(lat) if a.element.symbol != "H"]
    jh = [j for j, a in enumerate(rat) if a.element.symbol != "H"]
    d2 = d[np.ix_(ih, jh)]
    k = np.unravel_index(np.argmin(d2), d2.shape)
    print(f"minimo sin hidrogenos: {d2[k]:.2f} A "
          f"({lat[ih[k[0]]].name} - {rat[jh[k[1]]].residue.name}"
          f"{rat[jh[k[1]]].residue.id} {rat[jh[k[1]]].name})", flush=True)

    # cadenas presentes y cuantos atomos aporta cada una
    por_cadena = {}
    for a in rat:
        por_cadena[a.residue.chain.id] = por_cadena.get(a.residue.chain.id, 0) + 1
    print("atomos por cadena:", por_cadena, flush=True)

    # el ligando, ¿esta dentro del recorte?
    print(f"distancia del ligando al centro del recorte: "
          f"{np.linalg.norm(lp.mean(0) - c):.2f} A", flush=True)
    p = pose_pdbqt.leer_pose(fila["pose_pdbqt"], modelo=1)
    pes = np.array([[a[0], a[1], a[2]] for a in p["atoms"].values()])
    print(f"pose 1: {len(pes)} atomos pesados, centro {np.round(pes.mean(0), 3)}",
          flush=True)


if __name__ == "__main__":
    main()
