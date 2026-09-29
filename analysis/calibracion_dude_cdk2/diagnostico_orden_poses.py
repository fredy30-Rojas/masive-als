#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""El control de CDK2 medido bien: la correspondencia se DEDUCE, no se supone.

QUE PASO
--------
El control de redocking se declaro fallado con RMSD ~6,6 A en las seis variantes
del barrido. La medida se apoyaba en una suposicion escrita en su docstring:
"Meeko respeta el orden de atomos de la molecula y Vina-GPU mantiene ese orden en
las poses, asi que el emparejamiento por indice es una resta directa".

Es falsa. Meeko escribe los atomos pesados en OTRO orden: el PDBQT del ligando del
cristal, leido en orden de fichero, empieza C,C,C,C,C,C,O,C,C,C,N,... y el mol2
empieza C,C,C,C,C,C,C,C,C,... Restando por indice, cada coordenada cae en el atomo
equivocado: lo que se medía no era la pose, era una molecula con la geometria
desordenada. Y contra una geometria desordenada no salva ni el GetBestRMS, porque
GetBestRMS casa por simetria del grafo y una permutacion arbitraria de las
COORDENADAS no es una simetria.

La prueba esta en este mismo fichero: se busca la asignacion optima entre los
atomos del PDBQT de entrada y los del mol2, exigiendo mismo elemento. Si el PDBQT
acoplado es la pose del cristal, esa asignacion da ~0,00 A y devuelve la
permutacion exacta. La da.

Y entonces la pregunta que queda es si el veredicto cambia al medir bien. Se
responde aqui, para las seis variantes del barrido, SIN volver a usar la GPU: las
poses ya estan en disco.

Uso:
    python diagnostico_orden_poses.py                    # una variante
    python diagnostico_orden_poses.py --remedir-barrido  # las seis, desde disco
Salida: diagnostico_orden_poses.txt y barrido_redocking_cdk2_corregido.txt
"""
from __future__ import annotations

import argparse
import glob
import os
import time

import numpy as np

import control_redocking_cdk2 as C   # noqa: E402  (la receta, una sola copia)

BASE = C.BASE
SALIDA = os.path.join(BASE, "diagnostico_orden_poses.txt")
REMEDIDO = os.path.join(BASE, "barrido_redocking_cdk2_corregido.txt")


def poses_de(ruta):
    """Las poses de un PDBQT de salida, en el formato que lee `C.medir`."""
    modelos, actual = [], None
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            actual = {"energia": None, "pos": []}
        elif l.startswith("REMARK VINA RESULT") and actual is not None \
                and actual["energia"] is None:
            actual["energia"] = float(l.split()[3])
        elif l.startswith("ENDMDL"):
            if actual:
                modelos.append(actual)
            actual = None
        elif l.startswith(("ATOM", "HETATM")) and actual is not None:
            tipo = l.rsplit(None, 1)[-1].strip().upper()
            if tipo in C.H_TIPOS:
                continue
            actual["pos"].append([float(l[30:38]), float(l[38:46]),
                                  float(l[46:54])])
    if actual:
        modelos.append(actual)
    return modelos


def analizar(carpeta, mol, ref, log):
    """Mide una corrida entera: correspondencia, poses y veredicto."""
    entrada = os.path.join(carpeta, "ligands", "FAP_cristal.pdbqt")
    salidas = sorted(glob.glob(os.path.join(carpeta, "out", "*_out.pdbqt")))
    if not os.path.exists(entrada) or not salidas:
        log("   sin corrida que medir")
        return None
    tipos, coords = C.atomos_pesados_pdbqt(entrada)
    perm, rmsd_entrada = C.emparejar_orden(tipos, coords, mol, ref)
    log("   la entrada (el PDBQT que se acoplo) es la pose del cristal a %.4f A"
        % rmsd_entrada)
    log("   su orden de atomos es OTRO: permutacion %s"
        % " ".join(str(int(x)) for x in perm))
    modelos = poses_de(salidas[0])
    filas = [f for f in C.medir(mol, ref, modelos, perm, tipos)
             if f[2] is not None]
    if not filas:
        log("   sin poses comparables")
        return None
    return {"perm": perm, "rmsd_entrada": rmsd_entrada, "filas": filas}


def leer_liston_a():
    """El liston del proyecto. Se lee del script, no se copia a mano."""
    return C.LISTON_A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variante", default="caja20_depth32")
    ap.add_argument("--remedir-barrido", action="store_true",
                    help="vuelve a medir TODAS las carpetas de _control_cdk2/")
    args = ap.parse_args()

    lineas = []

    def log(m):
        print(m, flush=True)
        lineas.append(m)

    mol, ref = C.leer_ligando()
    liston = leer_liston_a()
    log("CONTROL DE CDK2 MEDIDO BIEN — %s" % time.strftime("%Y-%m-%d %H:%M"))
    log("ligando del cristal: %d atomos pesados | liston %.1f A"
        % (mol.GetNumAtoms(), liston))

    if args.remedir_barrido:
        carpetas = sorted(d for d in glob.glob(os.path.join(BASE, "_control_cdk2",
                                                            "*"))
                          if os.path.isdir(d))
        if not carpetas:
            raise SystemExit("no hay carpetas en _control_cdk2/")
        log("")
        log("REMEDICION DEL BARRIDO DESDE DISCO (sin volver a acoplar)")
        log("   %-14s %-9s %-5s %-9s %-9s %-9s %-9s %s"
            % ("variante", "entrada", "pose", "kcal/mol", "bolsillo", "piso",
               "alineado", "veredicto"))
        corregidas = []
        for carpeta in carpetas:
            etiqueta = os.path.basename(carpeta)
            r = analizar(carpeta, mol, ref, log)
            if not r:
                continue
            mejor = min(r["filas"], key=lambda f: f[2])
            pasa = mejor[2] < liston
            log("   %-14s %-9.4f %-5d %-9.2f %-9.2f %-9.2f %-9.2f %s"
                % (etiqueta, r["rmsd_entrada"], mejor[0], mejor[1], mejor[2],
                   mejor[3], mejor[4], "PASA" if pasa else "no pasa"))
            corregidas.append((etiqueta, r, mejor, pasa))
        if corregidas:
            con = [c for c in corregidas if c[3]]
            log("")
            log("   de %d variantes, pasan %d" % (len(corregidas), len(con)))
            log("")
            log("LECTURA")
            if len(con) == len(corregidas):
                log("   El control PASA en todas: el 'no pasa' del barrido era la")
                log("   medida, no el motor. Se estaba midiendo con las coordenadas")
                log("   puestas en el atomo equivocado.")
            elif con:
                log("   El control pasa en %d de %d variantes: el motor encuentra"
                    % (len(con), len(corregidas)))
                log("   la pose del cristal, pero depende del protocolo. Hay que"
                    " usar el que pasa.")
            else:
                log("   El control NO PASA en ninguna, y ahora con la medida buena.")
                log("   El motor no encuentra esta pose en esta caja: no es la"
                    " metrica.")
            # Se ESCRIBE, no se anade: este fichero es el resumen de una pasada
            # completa, y si se anadiese cada vez que se corre el diagnostico (que es
            # justo lo que hace el encadenador tras la busqueda agotada) el recuento
            # de variantes del informe del AUC contaria las mismas dos veces.
            with open(REMEDIDO, "w", encoding="utf-8") as f:
                f.write("=== remedicion %s ===\n" % time.strftime("%Y-%m-%d %H:%M"))
                for etiqueta, r, mejor, pasa in corregidas:
                    f.write("variante %-14s | entrada %.4f A | mejor pose %d |"
                            " %+.2f kcal/mol | bolsillo %.2f A | piso %.2f A |"
                            " alineado %.2f A | %s\n"
                            % (etiqueta, r["rmsd_entrada"], mejor[0], mejor[1],
                               mejor[2], mejor[3], mejor[4],
                               "PASA" if pasa else "no pasa"))
            log("")
            log("guardado: %s" % os.path.relpath(REMEDIDO, C.RAIZ))
    else:
        carpeta = os.path.join(BASE, "_control_cdk2", args.variante)
        log("")
        r = analizar(carpeta, mol, ref, log)
        if not r:
            raise SystemExit("nada que medir en %s"
                             % os.path.relpath(carpeta, BASE))
        log("")
        log("   %-6s %-9s %-10s %-9s %-10s %-9s %s"
            % ("pose", "kcal/mol", "bolsillo", "piso", "alineado", "centroides",
               "solape"))
        for i, e, bol, piso, ali, sol, dc in r["filas"]:
            log("   %-6d %-9.2f %-10.2f %-9.2f %-10.2f %-9.2f %.2f"
                % (i, e, bol, piso, ali, dc, sol))
        mejor = min(r["filas"], key=lambda f: f[2])
        primera = min(r["filas"], key=lambda f: f[1])
        log("")
        log("   la mejor por energia (la que lee el embudo) es la pose %d:"
            " %.2f A en el bolsillo, %.2f A tras alinear" % (primera[0], primera[2],
                                                              primera[4]))
        log("   la mejor por geometria es la pose %d: %.2f A (piso %.2f A)"
            % (mejor[0], mejor[2], mejor[3]))
        log("   VEREDICTO: %s (mejor pose %d a %.2f A del cristal en el bolsillo)"
            % ("PASA" if mejor[2] < liston else "NO PASA", mejor[0], mejor[2]))

    with open(SALIDA, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    print("\nguardado: %s" % SALIDA)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
