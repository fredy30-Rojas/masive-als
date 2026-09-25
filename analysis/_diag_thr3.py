#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Mecanismo exacto del fallo 'No template found for residue N (THR)'.

Prueba variantes del recorte y dice cuál deja montar el sistema del complejo.
"""
import math
import sys
import tempfile
import traceback

sys.path.insert(0, "/home/ubuntu/mmgbsa")

import pandas as pd  # noqa: E402

from rescoring_mmgbsa_robusto import (  # noqa: E402
    ligand_center, prepare_receptor, prepare_ligand, build_system)

LIGS = ["DEC_CHEMBL1414576", "DECH_CHEMBL3309988"]


def leer_atomos(pdb):
    out = []
    for line in open(pdb, errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            out.append((line[21], int(line[22:26]), line[17:20].strip(),
                        line[12:16].strip()))
    return out


def fragmentos(atomos):
    """Bloques de residuos consecutivos (mismo numero) dentro de cada cadena."""
    res = []
    for ch, num, rn, an in atomos:
        if not res or (res[-1][0], res[-1][1]) != (ch, num):
            res.append([ch, num, rn, [an]])
        else:
            res[-1][3].append(an)
    return res


def recorte(pdb, center, out, radio, modo):
    """Recorta a `radio` A; modo 'crudo' = tal cual, 'stubs' quita bloques de 1."""
    keep = set()
    lineas = []
    for line in open(pdb, errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            continue
        d = math.dist((float(line[30:38]), float(line[38:46]), float(line[46:54])),
                      tuple(center))
        if d < radio:
            keep.add((line[21], int(line[22:26])))
        lineas.append(line)
    if modo == "stubs":
        bloques = fragmentos(leer_atomos(pdb))
        for ch, num, rn, an in bloques:
            if len(an) == 1:
                keep.discard((ch, num))
    out_lines, serial = [], 0
    for line in lineas:
        if (line[21], int(line[22:26])) in keep:
            serial += 1
            out_lines.append(line[:6] + "%5d" % serial + line[11:].rstrip())
    open(out, "w").write("\n".join(out_lines) + "\nEND\n")
    return out


def cola(top, n=4):
    rs = list(top.residues())
    return [(r.index, r.chain.id if r.chain else "-", r.id, r.name,
             len(list(r.atoms()))) for r in rs[-n:]]


def prueba(lig, fila, modo):
    pose = fila["pose_pdbqt"]
    c = ligand_center(pose)
    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    recorte("estratos/receptor.pdb", c, tmp, 14.0, modo)
    print(f"  -- modo={modo}")
    try:
        top, pos = prepare_receptor(tmp)
    except Exception:
        print("     prepare_receptor FALLO"); traceback.print_exc(); return
    print(f"     receptor: {len(list(top.residues()))} res, cola={cola(top)}")

    from openmm.app import ForceField, Modeller
    ff = ForceField("amber14-all.xml", "implicit/obc2.xml")
    lig_top, lig_pos = prepare_ligand(fila["smiles"], pose, ff)
    mod = Modeller(top, pos)
    mod.add(lig_top, lig_pos)
    print(f"     complejo: {len(list(mod.topology.residues()))} res, "
          f"cola={cola(mod.topology, 5)}")
    try:
        build_system(mod.topology, ff)
        print("     SISTEMA OK")
    except Exception as e:
        print(f"     SISTEMA FALLO: {str(e)[:120]}")


def main():
    rows = pd.read_csv("estratos/candidatos.csv")
    for lig in LIGS:
        fila = rows[rows.ligand == lig].iloc[0]
        print("=" * 70)
        print(lig)
        for modo in ("crudo", "stubs"):
            prueba(lig, fila, modo)


if __name__ == "__main__":
    main()
