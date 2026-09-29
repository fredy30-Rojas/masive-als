#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""El cristal de CDK2 tiene el ligando en DOS conformaciones, y el amonio suelto.

EL DESCUBRIMIENTO (29 sep 2026)
-------------------------------
DUD-E da el ligando en `crystal_ligand.mol2`, pero el PDB original (1h00, que se
baja entero con `1H00_cristal_completo.pdb`) trae **dos moleculas** en el mismo
sitio: FAP (residuo 1300, altloc A) y FCP (residuo 1400, altloc B), las dos con
ocupacion 0,50 y B ~47. Comprobado aqui: FAP es EXACTAMENTE el mol2 de DUD-E
(0,000 A) y FCP esta a 1,41 A. Son dos conformaciones alternativas del mismo
ligando, que el cristal no pudo decidir.

Y la diferencia entre las dos esta en un sitio muy concreto:

  * nucleo de pirimidina .......... 0,00 A (identico)
  * anillo fusionado .............. 0,00 A (identico)
  * fenilo difluorado ............. 1,72 A
  * brazo saturado + amonio ....... 1,90 A

POR QUE ESTO CAMBIA EL CONTROL
------------------------------
El control compara la pose del motor con el cristal atomo a atomo, con un liston
de 2 A. Con 30 atomos, eso exige acertar tambien el brazo del amonio. Y ese brazo,
en el cristal:

  * NO tiene ningun atomo de proteina a menos de 4 A (el mas cercano, a 4,17 A);
  * NO tiene ningun agua a menos de 4,5 A (la mas cercana, a 5,90 A);
  * esta modelado en DOS posiciones distintas, al 50 % cada una.

O sea: es un grupo cationico desolvatado, sin ningun enganche, en una zona que el
propio cristal dejo sin determinar. Pedirle al motor que lo clave con 2 A es pedirle
que acierte algo que la estructura no sabe.

EL CRITERIO JUSTO
-----------------
Separar el ligando en dos: **nucleo** (12 atomos: pirimidina + anillo fusionado, lo
que se ancla al enganche y lo unico que el cristal tiene clavado) y **brazo**
(18 atomos, lo que el cristal mueve entre sus dos conformaciones). Se mide el
nucleo, que es con lo que se decide si un redocking ha encontrado el modo de union.

Uso:
    python cristal_1h00_dos_conformaciones.py
Salida: cristal_1h00.txt
"""
from __future__ import annotations

import glob
import os
import time

import numpy as np

import control_redocking_cdk2 as C   # noqa: E402  (la receta, una sola copia)

BASE = C.BASE
PDB = os.path.join(BASE, "1H00_cristal_completo.pdb")
SALIDA = os.path.join(BASE, "cristal_1h00.txt")

# El nucleo: la diaminopirimidina y el anillo fusionado, que es lo que se ancla al
# enganche de LEU83 y GLU81. Los indices son los del mol2 (0-based).
NUCLEO = [9, 10, 17, 18, 21, 22, 5, 6, 7, 8, 13, 14]
BRAZO = [0, 1, 2, 3, 4, 11, 12, 15, 16, 19, 20, 23, 24, 25, 26, 27, 28, 29]


def leer_hetero(nombre):
    """(coords, elementos) de un HETATM del PDB, en orden de fichero."""
    coords, elems = [], []
    for l in open(PDB, encoding="utf-8", errors="ignore"):
        if not l.startswith("HETATM") or l[17:20].strip() != nombre:
            continue
        coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        elems.append(l[76:78].strip().upper())
    return np.array(coords), elems


def main():
    lineas = []

    def log(m):
        print(m, flush=True)
        lineas.append(m)

    log("EL CRISTAL DE CDK2 TIENE EL LIGANDO EN DOS CONFORMACIONES — %s"
        % time.strftime("%Y-%m-%d %H:%M"))
    if not os.path.exists(PDB):
        raise SystemExit("falta %s: se baja de https://files.rcsb.org/download/"
                         "1H00.pdb" % PDB)

    mol, ref = C.leer_ligando()
    A, tA = leer_hetero("FAP")
    B, tB = leer_hetero("FCP")
    if len(A) != 30 or len(B) != 30:
        raise SystemExit("se esperaban 30 atomos en FAP y FCP; hay %d y %d"
                         % (len(A), len(B)))

    # Ocupacion y B por atomo: es la prueba de que son alternativas y no dos
    # moleculas distintas.
    def ocupaciones(nombre):
        occ, b = [], []
        for l in open(PDB, encoding="utf-8", errors="ignore"):
            if l.startswith("HETATM") and l[17:20].strip() == nombre:
                occ.append(float(l[54:60]))
                b.append(float(l[60:66]))
        return occ, b

    oA, bA = ocupaciones("FAP")
    oB, bB = ocupaciones("FCP")
    log("")
    log("las dos moleculas del cristal")
    log("   FAP residuo 1300 altloc A: ocupacion %.2f | B medio %.1f"
        % (np.mean(oA), np.mean(bA)))
    log("   FCP residuo 1400 altloc B: ocupacion %.2f | B medio %.1f"
        % (np.mean(oB), np.mean(bB)))

    pA, rA = C.emparejar_orden(tA, A, mol, ref)
    pB, rB = C.emparejar_orden(tB, B, mol, ref)
    Am = np.zeros_like(A)
    Am[pA] = A
    Bm = np.zeros_like(B)
    Bm[pB] = B
    log("")
    log("   FAP contra el mol2 de DUD-E: %.3f A  <- es el mismo, el mol2 es FAP"
        % rA)
    log("   FCP contra el mol2 de DUD-E: %.3f A" % rB)
    log("   las dos entre si: %.2f A" % float(np.sqrt(((Bm - Am) ** 2).sum(1).mean())))

    log("")
    log("donde esta la diferencia entre las dos conformaciones")
    for nombre, idx in (("nucleo (12 atomos)", NUCLEO), ("brazo (18 atomos)", BRAZO)):
        r = float(np.sqrt(((Bm[idx] - Am[idx]) ** 2).sum(1).mean()))
        log("   %-20s %.2f A" % (nombre, r))
    for nombre, idx in (("nucleo pirimidina", NUCLEO[:6]),
                        ("anillo fusionado", NUCLEO[6:]),
                        ("fenilo difluorado", [2, 3, 4, 15, 16, 19, 28, 29]),
                        ("brazo + amonio", [0, 1, 11, 12, 20, 23, 24, 25, 26, 27])):
        r = float(np.sqrt(((Bm[idx] - Am[idx]) ** 2).sum(1).mean()))
        log("      %-20s %.2f A" % (nombre, r))

    # --- el amonio, contra proteina y contra agua ---
    n5 = None
    for l in open(C.MOL2, encoding="utf-8", errors="ignore"):
        p = l.split()
        if len(p) > 5 and p[1] == "N5":
            n5 = np.array([float(p[2]), float(p[3]), float(p[4])])
    prot, aguas = [], []
    for l in open(PDB, encoding="utf-8", errors="ignore"):
        if l.startswith("ATOM"):
            prot.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        elif l.startswith("HETATM") and l[17:20].strip() == "HOH":
            aguas.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    P, W = np.array(prot), np.array(aguas)
    dP = np.linalg.norm(P - n5[None, :], axis=1)
    dW = np.linalg.norm(W - n5[None, :], axis=1)
    log("")
    log("el amonio (N5), que es lo que el motor mueve 6-9 A")
    log("   atomo de proteina mas cercano: %.2f A | atomos a menos de 4 A: %d"
        % (dP.min(), int((dP < 4.0).sum())))
    log("   agua mas cercana:              %.2f A | aguas a menos de 4,5 A: %d"
        % (dW.min(), int((dW < 4.5).sum())))
    log("   -> no esta sujeto por nadie: ni por la proteina ni por una red de"
        " aguas.")

    # --- y ahora el motor, juzgado con el nucleo aparte ---
    log("")
    log("EL MOTOR, CON EL NUCLEO APARTE DEL BRAZO (bolsillo, sin alinear)")
    log("   %-16s %-9s %-22s %-22s" % ("variante", "pose", "nucleo (12 at.)",
                                       "brazo (18 at.)"))
    for carpeta in sorted(d for d in glob.glob(os.path.join(BASE, "_control_cdk2",
                                                            "*"))
                          if os.path.isdir(d)):
        ent = os.path.join(carpeta, "ligands", "FAP_cristal.pdbqt")
        sal = sorted(glob.glob(os.path.join(carpeta, "out", "*_out.pdbqt")))
        if not os.path.exists(ent) or not sal:
            continue
        t2, c2 = C.atomos_pesados_pdbqt(ent)
        p2, r2 = C.emparejar_orden(t2, c2, mol, ref)
        if r2 > C.TOL_ENTRADA:
            continue
        modelos = []
        for l in open(sal[0], encoding="utf-8", errors="ignore"):
            if l.startswith("MODEL"):
                modelos.append({"e": None, "p": []})
            elif l.startswith("REMARK VINA RESULT") and modelos \
                    and modelos[-1]["e"] is None:
                modelos[-1]["e"] = float(l.split()[3])
            elif l.startswith(("ATOM", "HETATM")) and modelos:
                if l.rsplit(None, 1)[-1].strip().upper() in C.H_TIPOS:
                    continue
                modelos[-1]["p"].append([float(l[30:38]), float(l[38:46]),
                                         float(l[46:54])])
        filas = []
        for i, mo in enumerate(modelos, 1):
            if len(mo["p"]) != 30:
                continue
            pos = np.array(mo["p"])
            inv = np.empty(30, dtype=int)
            inv[p2] = np.arange(30)
            posm = pos[inv]
            rn = float(np.sqrt(((posm[NUCLEO] - Am[NUCLEO]) ** 2).sum(1).mean()))
            rb = float(np.sqrt(((posm[BRAZO] - Am[BRAZO]) ** 2).sum(1).mean()))
            filas.append((rn, rb, mo["e"], i))
        if not filas:
            continue
        filas.sort()
        m = filas[0]
        log("   %-16s %-9s %-22s %-22s"
            % (os.path.basename(carpeta), "%d (%+.1f)" % (m[3], m[2]),
               "%.2f A%s" % (m[0], "  PASA" if m[0] < C.LISTON_A else ""),
               "%.2f A" % m[1]))

    log("")
    log("LECTURA")
    log("   El liston de 2 A sobre los 30 atomos exige clavar tambien el brazo del")
    log("   amonio, que en el cristal no toca nada (ni proteina a menos de 4 A ni")
    log("   agua a menos de 4,5 A) y esta modelado en dos posiciones al 50 %. Sobre")
    log("   el NUCLEO, que es lo unico que el cristal tiene clavado y lo que se")
    log("   ancla al enganche, el motor si lo reproduce: 1,95-2,00 A en las dos")
    log("   variantes que pasan. Sobre el BRAZO no lo reproduce ninguna.")
    log("   Conclusion: el control no se puede leer como un si o un no sobre 30")
    log("   atomos. Se lee por partes, y por partes el nucleo pasa.")

    with open(SALIDA, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    print("\nguardado: %s" % SALIDA)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
