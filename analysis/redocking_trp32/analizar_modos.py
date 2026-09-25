#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿El fallo del redocking es de puntuacion o de muestreo? (20 sep 2026).

El control de redocking mira el modo MEJOR PUNTUADO: si su RMSD supera 2 A, el
protocolo falla. Pero ese unico numero no dice por que. Hay dos causas posibles
y se distinguen mirando los nueve modos que Vina guarda:

  * FALLO DE PUNTUACION: la pose correcta SI se genero, pero la funcion de
    puntuacion la coloco en un puesto malo. Se arregla con otra funcion.
  * FALLO DE MUESTREO: la pose correcta NO esta entre los nueve modos. No hay
    nada que puntuar; el problema es el receptor, la caja o la rigidez del
    ligando (por ejemplo, un anillo fijado en una conformacion distinta a la
    del cristal, que el acoplamiento rigido no puede cambiar).

Uso: python analizar_modos.py [--caja 24]
"""
import argparse
import glob
import os
import sys

from redock_trp32 import SISTEMAS, molecula_dockeada, rmsd_en_sitio, BASE

from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--caja", type=int, default=24)
    ap.add_argument("--exhaustividad", type=int, default=8)
    args = ap.parse_args()

    print("=== RMSD de CADA modo de Vina (caja de %d A, exhaustividad %d) ==="
          % (args.caja, args.exhaustividad))
    for entrada in SISTEMAS:
        nombre = SISTEMAS[entrada][2]
        lig_pdb = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
        salidas = sorted(glob.glob(os.path.join(
            BASE, "out", "%s_caja%d_e%d_s*.pdbqt"
            % (entrada, args.caja, args.exhaustividad))))
        if not salidas or not os.path.exists(lig_pdb):
            print("\n%s: sin datos" % entrada)
            continue
        print("\n%s (%s)" % (entrada, nombre))
        for f in salidas:
            semilla = os.path.basename(f).rsplit("_s", 1)[1].replace(".pdbqt", "")
            # todos los modos del fichero: meeko devuelve una molecula con
            # tantas conformaciones como poses traiga el PDBQT
            from meeko import PDBQTMolecule, RDKitMolCreate
            todos = RDKitMolCreate.from_pdbqt_mol(
                PDBQTMolecule.from_file(f, skip_typing=True))
            vals = []
            for m in todos:
                for cid in range(m.GetNumConformers()):
                    v = rmsd_en_sitio(lig_pdb, m, conf_id=cid)
                    vals.append(v if v is not None else float("nan"))
            mejor = min(vals) if vals else None
            orden = " ".join("%.1f" % v for v in vals[:9])
            sobre2 = sum(1 for v in vals if v <= 2.0)
            print("   semilla %-5s modos: %s | mejor %.2f A | modos <=2 A: %d/%d"
                  % (semilla, orden, mejor, sobre2, len(vals)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
