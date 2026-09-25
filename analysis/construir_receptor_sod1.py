#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Receptor de SOD1 limpio: UN dimer biologico y sin disolvente (21 sep 2026).

POR QUE
-------
El receptor con el que se cribo SOD1 (`gpu_dock/SOD1.pdbqt`) no era una proteina
limpia. Medido hoy:

    21.585 atomos pesados | 1.763 aguas cristalograficas | 18 cadenas | Ca, Cu, Zn

Tres consecuencias, todas malas:

1. Las aguas. Vina las trata como parte RIGIDA del receptor: no se pueden
   desplazar. Cinco caen dentro del bolsillo de Trp32. Un ligando que se aprieta
   entre la proteina y aguas congeladas cobra contactos que en la realidad no
   existirian, porque en la realidad esas aguas se van.
2. Las 18 copias. La caja de 22 A centrada en Trp32 no solo ve la proteina: ve lo
   que tenga al lado, y lo que tenia al lado era el contacto de empaquetamiento
   cristalino entre las cadenas A y J. Ahi se iban a meter las poses de -10
   kcal/mol de las quinazolinas, en una grieta que no es un sitio de union de
   nada (ver INFORME_REDOCKING_QUIMIAS_NUEVAS_2026-09-21.md).
3. Los iones. Cu y Zn son los metales cataliticos, estan al otro extremo de la
   molecula (a mas de 20 A del bolsillo de Trp32) y Vina no los parametriza bien.
   Los controles de redocking que SI pasan (las catecolaminas) se hicieron con
   receptores sin metales, asi que se dejan fuera por coherencia con el protocolo
   ya validado.

QUE HACE ESTE SCRIPT
--------------------
1. Identifica la unidad biologica que contiene la cadena A. El REMARK 350 de
   1HL5 dice: BIOMOLECULE 1, DIMERICA, cadenas A y H (sin operador de simetria).
   Se genera exactamente eso.
2. Escribe el receptor solo con los atomos de proteina de esas dos cadenas,
   resolviendo las conformaciones alternas por ocupacion (misma rutina que los
   controles: `redock_trp32.extraer_receptor`).
3. Lo prepara con meeko (`mk_prepare_receptor -p -j -a`).
4. COMPRUEBA que el bolsillo no se ha movido: distancia de cada atomo del
   bolsillo de Trp32 al receptor viejo y al nuevo, y distancia del centro de la
   caja al atomo mas cercano. La caja del cribado (46,5 / 80,0 / 73,3; 22 A) solo
   sigue siendo valida si la cadena A no se ha tocado.

Uso: python construir_receptor_sod1.py
Salida: ../gpu_dock/SOD1_limpio.pdb y ../gpu_dock/SOD1_limpio.pdbqt
"""
import os
import subprocess
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
GPU = os.path.abspath(os.path.join(BASE, "..", "gpu_dock"))
ORIGEN = os.path.abspath(os.path.join(BASE, "..", "proteins", "SOD1", "PDB-1hl5.pdb"))
CADENAS = {"A", "H"}
CENTRO_CAJA = (46.5, 80.0, 73.3)
TAMANO_CAJA = 22

sys.path.insert(0, os.path.join(BASE, "redocking_trp32"))
import redock_trp32 as R  # noqa: E402


def log(m):
    print(m, flush=True)


def atomos_pesados(pdbqt):
    """(residuo, cadena, num, xyz) de los atomos pesados de un PDBQT."""
    out = []
    for l in open(pdbqt, encoding="utf-8", errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        tipo = l.rsplit(None, 1)[-1].upper()
        if tipo in ("H", "HD", "HS", "D", "DD"):
            continue
        out.append((l[17:20].strip(), l[21], l[22:27].strip(),
                    np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


def por_nombre(pdbqt, cadena):
    """(num_residuo, nombre_atomo) -> coordenada, para una cadena."""
    d = {}
    for l in open(pdbqt, encoding="utf-8", errors="ignore"):
        if not l.startswith("ATOM") or l[21] != cadena:
            continue
        d[(l[22:27].strip(), l[12:16].strip())] = np.array(
            [float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return d


def descripcion(pdbqt):
    tot = hoh = 0
    cadenas = set()
    for l in open(pdbqt, encoding="utf-8", errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        tot += 1
        if l[17:20].strip() in ("HOH", "DOD"):
            hoh += 1
        if l.startswith("ATOM"):
            cadenas.add(l[21])
    return tot, hoh, sorted(cadenas)


def main():
    limpio = os.path.join(GPU, "SOD1_limpio.pdb")
    base = os.path.join(GPU, "SOD1_limpio")
    pdbqt = base + ".pdbqt"

    if not os.path.exists(ORIGEN):
        log("falta el PDB de origen: %s" % ORIGEN)
        return 1

    n = R.extraer_receptor(ORIGEN, limpio, cadenas=CADENAS)
    log("receptor limpio escrito: %d atomos de proteina (cadenas %s)"
        % (n, "+".join(sorted(CADENAS))))

    r = subprocess.run([R.PYTHON, R.MK_RECEPTOR, "--read_pdb", limpio,
                        "-o", base, "-p", "-j", "-a", "--default_altloc", "A"],
                       capture_output=True, text=True, timeout=1800)
    if not os.path.exists(pdbqt):
        log("FALLO preparando con meeko:\n%s"
            % ((r.stdout or "") + (r.stderr or ""))[-800:])
        return 1

    viejo = os.path.join(GPU, "SOD1.pdbqt")
    for etiqueta, ruta in (("antes", viejo), ("limpio", pdbqt)):
        if not os.path.exists(ruta):
            continue
        tot, hoh, cad = descripcion(ruta)
        log("  %-7s %6d atomos | %5d aguas | %2d cadenas (%s)"
            % (etiqueta, tot, hoh, len(cad), "".join(cad)))

    # --- comprobacion de que el bolsillo no se ha movido ---
    # Se compara la cadena A del receptor viejo con la del limpio, atomo por atomo
    # por (residuo, nombre). Compararlo con las 18 cadenas no dice nada; lo que
    # importa es que la cadena A del cribado siga donde estaba, porque la caja de
    # 22 A se definio sobre sus coordenadas.
    a = por_nombre(viejo, "A")
    b = por_nombre(pdbqt, "A")
    comunes = sorted(set(a) & set(b))
    b_faltan = sorted(set(b) - set(a))
    desv = np.array([float(np.linalg.norm(a[k] - b[k])) for k in comunes])
    log("")
    log("cadena A del receptor viejo (%d atomos) frente al limpio (%d atomos):"
        % (len(a), len(b)))
    log("   %d atomos comunes | desplazamiento mediano %.4f A | maximo %.4f A"
        % (len(comunes), float(np.median(desv)), float(desv.max())))
    log("   atomos que solo estan en el limpio (altloc resuelto de otra forma): %d %s"
        % (len(b_faltan), b_faltan[:6]))
    b = atomos_pesados(pdbqt)
    A = np.array([x[3] for x in atomos_pesados(viejo) if x[0] not in ("HOH", "DOD")])
    B = np.array([x[3] for x in b])

    centro = np.array(CENTRO_CAJA)
    d_centro = np.sqrt(((B - centro) ** 2).sum(-1))
    log("   atomo mas cercano al centro de la caja del cribado: %.2f A (caja de %d A)"
        % (d_centro.min(), TAMANO_CAJA))
    cerca = np.sqrt(((B - centro) ** 2).sum(-1)) <= TAMANO_CAJA / 2.0
    log("   atomos dentro de la caja: %d (antes %d)"
        % (int(cerca.sum()),
           int((np.sqrt(((A - centro) ** 2).sum(-1))
                <= TAMANO_CAJA / 2.0).sum())))
    return 0


if __name__ == "__main__":
    sys.exit(main())
