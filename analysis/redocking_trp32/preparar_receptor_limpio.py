#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Receptores preparados SIN que se muevan los atomos del bolsillo (20 sep 2026).

EL PROBLEMA
-----------
En 4A7S, al comparar `receptores/4A7S.pdbqt` con los atomos del PDB, 13 de los
31 atomos de proteina que rodean al ligando NO estaban en su sitio (hasta 2,3 A
de desplazamiento), y entre ellos estaba la **Trp32**, con el anillo indol
girado (CD1 intercambiado con CD2). Eso importa porque la 5-fluorouridina se
ancla justo apilando con la Trp32: con el bolsillo reconstruido, la pose del
cristal cae en una zona de choque (mover la pose +2 A en x la lleva de +0,04 a
+14,8 kcal/mol) y puntua ~0 en vez de puntuar como un union real. El acoplamiento
se va a 10 A. Los resultados de 4A7S de los controles anteriores no son fiables.

CAUSA PROBABLE
--------------
Estas entradas traen conformaciones alternativas (Lys30 A/B al 50 %, Trp32 CA/CB
70/30). Meeko decide con `--default_altloc`, y con ese ajuste reconstruyo parte
de las cadenas laterales del bolsillo en vez de conservar los atomos medidos.

QUE HACE ESTE SCRIPT
--------------------
1. Escribe un PDB "limpio" por entrada conservando SOLO los atomos con
   conformacion alterna vacia o A (los que la mayoria de los articulos usan) y
   descartando aguas, iones y hetero-atomos.
2. Prepara el receptor con `mk_prepare_receptor` (mismo meeko, mismos ajustes).
3. Comprueba, atomo por atomo, que los atomos de proteina que tocan al ligando
   del cristal siguen donde el PDB dice, y puntua la pose del cristal con Vina.

Uso: python preparar_receptor_limpio.py [--entradas 4A7S,4A7T]
"""
import argparse
import os
import subprocess
import sys

import numpy as np

import protocolo_alternativo as P
import redock_trp32 as R

BASE = R.BASE
PDB = os.path.join(BASE, "pdb")


def escribir_limpio(entrada, destino):
    """Atomos de protein con altloc vacio o A; fuera aguas, iones y ligandos."""
    n, fuera = 0, 0
    with open(destino, "w", encoding="utf-8") as f:
        for l in open(os.path.join(PDB, entrada + ".pdb"), encoding="utf-8",
                      errors="ignore"):
            if l.startswith("ATOM"):
                altloc = l[16]
                if altloc not in (" ", "A"):
                    fuera += 1
                    continue
                f.write(l)
                n += 1
        f.write("END\n")
    return n, fuera


def atomos(f, prefijos=("ATOM  ", "HETATM")):
    out = []
    for l in open(f, encoding="utf-8", errors="ignore"):
        if l.startswith(prefijos):
            nombre = l[12:16].strip()
            elem = (l[76:78].strip() or nombre[0]).upper()
            out.append((l[17:20].strip(), l[21], l[22:27].strip(), nombre, elem,
                        np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


def desplazamientos(rec_pdb, rec_pdbqt, lig_pdb, radio=4.5):
    """Atomos del bolsillo que en el PDBQT no estan donde el PDB dice."""
    lig = np.array(P.coord_pesados(lig_pdb))
    pdb = [a for a in atomos(rec_pdb, ("ATOM",)) if a[4] != "H"]
    q = [a for a in atomos(rec_pdbqt) if a[4] != "H"]
    if not q:
        return None, []
    Q = np.array([a[5] for a in q])
    cerca = [a for a in pdb if np.sqrt(((lig - a[5]) ** 2).sum(-1)).min() < radio]
    malos = []
    for a in cerca:
        d = np.sqrt(((Q - a[5]) ** 2).sum(-1))
        j = int(d.argmin())
        if float(d[j]) > 0.05:
            malos.append((a[0], a[2], a[3], float(d[j])))
    return len(cerca), malos


def score(rec, lig, c, lado=24):
    txt = P.correr_vina(["--receptor", rec, "--ligand", lig, "--scoring", "vina",
                         "--score_only", "--center_x", "%.3f" % c[0],
                         "--center_y", "%.3f" % c[1], "--center_z", "%.3f" % c[2],
                         "--size_x", str(lado), "--size_y", str(lado),
                         "--size_z", str(lado)], timeout=1800)
    return P.afinidad_stdout(txt)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entradas", default=",".join(R.SISTEMAS))
    args = ap.parse_args()

    for entrada in [e.strip() for e in args.entradas.split(",") if e.strip()]:
        limpio = os.path.join(BASE, "receptores", "%s_limpio.pdb" % entrada)
        n, fuera = escribir_limpio(entrada, limpio)
        base = os.path.join(BASE, "receptores", "%s_limpio" % entrada)
        rec = base + ".pdbqt"
        if not os.path.exists(rec):
            r = subprocess.run([R.PYTHON, R.MK_RECEPTOR, "--read_pdb", limpio,
                                "-o", base, "-p", "-j",
                                "--default_altloc", "A"],
                               capture_output=True, text=True, timeout=1800)
            if not os.path.exists(rec):
                print("%s: FALLO preparando el receptor limpio: %s"
                      % (entrada, ((r.stdout or "") + (r.stderr or ""))[-300:]))
                continue
        lig = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
        c = R.centroide(lig)
        cry_pdbqt = os.path.join(BASE, "ligands", "%s_cristal_vina.pdbqt" % entrada)

        for etiqueta, ruta in (("de antes", os.path.join(BASE, "receptores",
                                                        entrada + ".pdbqt")),
                               ("limpio", rec)):
            total, malos = desplazamientos(limpio, ruta, lig)
            e = score(ruta, cry_pdbqt, c)
            print("%-5s %-8s atomos del bolsillo: %-3s desplazados: %-3d | E(pose del cristal) %s kcal/mol"
                  % (entrada, etiqueta, total, len(malos),
                     "%.2f" % e if e is not None else "n/d"))
            if malos:
                print("        peores: %s"
                      % ", ".join("%s%s %s a %.2f A" % m for m in
                                  sorted(malos, key=lambda x: -x[3])[:4]))
        print("   PDB limpio: %d atomos (%d descartados por conformacion alterna)\n"
              % (n, fuera))
    return 0


if __name__ == "__main__":
    sys.exit(main())
