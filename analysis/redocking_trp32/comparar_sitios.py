#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿Los ligandos co-cristalizados ocupan de verdad el mismo sitio? (20 sep 2026).

El proyecto llama «ligandos de Trp32» a 5-fluorouridina, isoproterenol,
dopamina y adrenalina, que vienen de cuatro entradas PDB distintas (4A7S, 4A7T,
4A7V, 4A7U). Cada entrada trae, ademas, varias copias del ligando repartidas por
la superficie, y no todas estan en Trp32. Aqui se elige la copia buena con
criterio explicito (la que mas residuos del bolsillo toca) y se responde:

  1. ¿A que distancia esta cada ligando de los demas, ya superpuestas las
     proteinas (153 CA de la cadena A, Kabsch)?
  2. ¿Que residuos del receptor toca cada uno (a 4,5 A)?
  3. ¿Cuanto se superponen esos conjuntos de residuos entre pares de ligandos?
  4. ¿A que distancia esta cada uno del Trp32?

Uso: python comparar_sitios.py
"""
import os
import sys

import numpy as np

from redock_trp32 import PDB, SISTEMAS, atomos_receptor, elegir_copia

BASE = os.path.dirname(os.path.abspath(__file__))
REFERENCIA = "4A7T"


def log(m):
    print(m, flush=True)


def ca_y_ligando(entrada, codigo, cadena, resseq):
    ca, lig = [], []
    for l in open(os.path.join(PDB, entrada + ".pdb"), encoding="utf-8",
                  errors="ignore"):
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == cadena:
            ca.append((l[22:27].strip(),
                       [float(l[30:38]), float(l[38:46]), float(l[46:54])]))
        if l.startswith("HETATM") and l[17:20].strip() == codigo \
                and l[21] == cadena and l[22:27].strip() == resseq:
            lig.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return dict(ca), np.array(lig)


def kabsch(P, Q):
    """Matriz de rotacion que lleva P sobre Q (ambas N x 3 ya emparejadas)."""
    Pc, Qc = P.mean(0), Q.mean(0)
    H = (P - Pc).T @ (Q - Qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, Pc, Qc


def main():
    log("Eleccion de la copia del ligando (la del bolsillo de Trp32):")
    datos = {}
    for e, (cod, cad, nombre) in SISTEMAS.items():
        rs = elegir_copia(os.path.join(PDB, e + ".pdb"), cod, cad)
        ca, lig = ca_y_ligando(e, cod, cad, rs)
        datos[e] = (ca, lig)
        log("   %s (%s): %d CA, %d atomos de ligando en resSeq %s"
            % (e, nombre, len(ca), len(lig), rs))

    ref = datos[REFERENCIA][0]
    comun = sorted([r for r in ref if all(r in datos[e][0] for e in SISTEMAS)],
                   key=lambda x: int(x))
    log("\nresiduos CA comunes para superponer: %d" % len(comun))
    Q = np.array([ref[r] for r in comun])
    rec = atomos_receptor(REFERENCIA)

    centros, lig_ref = {}, {}
    for e, (ca, lig) in datos.items():
        if e == REFERENCIA:
            lig_r = lig
        else:
            P = np.array([ca[r] for r in comun])
            R, Pc, Qc = kabsch(P, Q)
            rmsd_ca = np.sqrt(((((P - Pc) @ R.T) + Qc - Q) ** 2).sum(1).mean())
            lig_r = ((lig - Pc) @ R.T) + Qc
            log("%s  RMSD de CA tras superponer: %.2f A" % (e, rmsd_ca))
        centros[e] = lig_r.mean(0)
        lig_ref[e] = lig_r

    log("\n=== 1. DISTANCIA ENTRE CENTROIDES DE LIGANDO (mismo sistema) ===")
    claves = list(centros)
    for i in range(len(claves)):
        for j in range(i + 1, len(claves)):
            d = np.linalg.norm(centros[claves[i]] - centros[claves[j]])
            log("   %s vs %s : %6.1f A%s" % (claves[i], claves[j], d,
                                             "   <-- otro sitio" if d > 6 else ""))

    log("\n=== 2. RESIDUOS DEL RECEPTOR A MENOS DE 4,5 A ===")
    contactos = {}
    for e in SISTEMAS:
        lig = lig_ref[e]
        contactos[e] = {(ch, rs) for ch, rs, _, xyz in rec
                        if np.linalg.norm(lig - xyz, axis=1).min() <= 4.5}
        log("   %s (%s): %d residuos -> %s"
            % (e, SISTEMAS[e][2], len(contactos[e]),
               " ".join("%s%s" % (c, r) for c, r in sorted(
                   contactos[e], key=lambda x: (x[0], int(x[1]))))))

    log("\n=== 3. SUPERPOSICION DE RESIDUOS CONTACTADOS (Jaccard) ===")
    for i in range(len(claves)):
        for j in range(i + 1, len(claves)):
            a, b = contactos[claves[i]], contactos[claves[j]]
            jac = len(a & b) / len(a | b)
            log("   %s vs %s : %.3f" % (claves[i], claves[j], jac))

    log("\n=== 4. DISTANCIA MINIMA AL Trp32 ===")
    for e in SISTEMAS:
        lig = lig_ref[e]
        d = min(np.linalg.norm(lig - xyz, axis=1).min()
                for ch, rs, rn, xyz in rec if rs == "32" and rn == "TRP")
        log("   %s (%s): %.1f A" % (e, SISTEMAS[e][2], d))
    return 0


if __name__ == "__main__":
    sys.exit(main())
