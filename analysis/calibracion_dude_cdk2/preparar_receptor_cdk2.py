#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara el receptor y la caja del banco CDK2 de DUD-E.

QUE HACE
--------
DUD-E trae `receptor.pdb` (1h00, proteina sola, sin hidrogenos: meeko los
anhade) y `crystal_ligand.mol2`, el aminopurina del cristal. La receta es la
misma que la del banco de TBK1 (`construir_receptor_tbk1.py`):

  1. Se limpia el PDB: solo ATOM de proteina, sin aguas ni iones.
  2. Se conservan las cadenas que tocan el bolsillo (atomos a menos de
     RADIO_BOLSILLO del ligando del cristal) mas la cadena del ligando entera.
  3. La caja es de 24 A centrada en el ligando del cristal, la convencion de
     DUD-E y la misma del proyecto.
  4. El receptor se prepara con meeko (`mk_prepare_receptor -p -j -a
     --default_altloc A`), igual que TBK1.

Salida: `receptor_cdk2.pdbqt`, `caja_cdk2.json` y este log con los residuos del
bolsillo, para poder comparar con la bisagra (Leu83) y la puerta (Glu81, Phe80)
que publica DUD-E para CDK2.

Uso:
    python preparar_receptor_cdk2.py
"""
from __future__ import annotations

import json
import os
import subprocess

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
PDB = os.path.join(BASE, "receptor.pdb")
MOL2 = os.path.join(BASE, "crystal_ligand.mol2")
LIMPIO = os.path.join(BASE, "receptor_cdk2.pdb")
PDBQT = os.path.join(BASE, "receptor_cdk2.pdbqt")
CAJA = os.path.join(BASE, "caja_cdk2.json")
PYDIR = r"C:\Users\Fredy\AppData\Local\Python\pythoncore-3.14-64"
MK_RECEPTOR = os.path.join(PYDIR, "Scripts", "mk_prepare_receptor.exe")

TAMANO_CAJA = 24         # igual que TBK1: 4-5 A de margen alrededor del ligando
RADIO_BOLSILLO = 8.0     # para listar residuos y elegir cadenas

RESIDUOS_CLAVE = {"Leu83": 83, "Glu81": 81, "Phe80": 80, "Ile10": 10}


def log(m):
    print(m, flush=True)


def coord(linea):
    return np.array([float(linea[30:38]), float(linea[38:46]),
                     float(linea[46:54])])


def ligando_del_mol2():
    """Centro del ligando del cristal, desde el bloque TRIPOS>ATOM."""
    pts = []
    dentro = False
    with open(MOL2, encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("@<TRIPOS>ATOM"):
                dentro = True
                continue
            if l.startswith("@<TRIPOS>") and dentro:
                break
            if dentro:
                parts = l.split()
                if len(parts) >= 6:
                    pts.append([float(parts[2]), float(parts[3]),
                                float(parts[4])])
    if not pts:
        raise SystemExit("no pude leer el ligando de %s" % MOL2)
    return np.array(pts)


def normalizar():
    """Escribe el PDB en formato moderno: meeko exige la columna de elemento
    (77-78) y el `receptor.pdb` de DUD-E es formato antiguo de 54 columnas, sin
    ella. Se aprovecha para quitar los hidrogenos de nomenclatura vieja (HN, HA,
    ...): meeko regenera los polares desde la geometria.

    Elemento de un PDB legacy de aminoacidos: la primera letra del nombre del
    atomo (N, CA, C, O, CB, HN -> N/C/C/O/C/H). Los simbolos de dos letras solo
    aparecen en selenometioninas y aqui no las hay.

    Ademas, el fichero de 2012 trae un residuo `LEV` que es exactamente una
    leucina (sus 8 atomos pesados, N CA C O CB CG CD1 CD2): es la Leu83 de la
    bisagra con el nombre equivocado. Meeko se para ahi ("unknown residues:
    {'LEV'}"); HID/HIE/HIP si estan en sus plantillas y se quedan como estan.
    """
    n = 0
    hidrogenos = 0
    with open(LIMPIO, "w", encoding="utf-8") as g:
        for l in open(PDB, encoding="utf-8", errors="ignore"):
            if not l.startswith("ATOM"):
                continue
            if l[17:20] == "LEV":  # leucina mal escrita en el fichero de DUD-E
                l = l[:17] + "LEU" + l[20:]
            nombre = l[12:16].strip()
            elemento = nombre[0].upper()
            if elemento == "H":
                hidrogenos += 1
                continue
            # formato completo: ocupacion 1.00 y B 0.00, elemento en 77-78
            linea = (l[:54] + "  1.00  0.00          " + " %1s" % elemento)
            g.write(linea[:78].ljust(78) + "\n")
            n += 1
        g.write("END\n")
    log("receptor normalizado: %d atomos pesados, %d hidrogenos fuera"
        % (n, hidrogenos))
    return n


def main():
    lig = ligando_del_mol2()
    centro_lig = lig.mean(axis=0)
    log("ligando del cristal: %d atomos pesados, centro %.2f, %.2f, %.2f"
        % (len(lig), *centro_lig))

    atomos = []
    with open(PDB, encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("ATOM"):
                atomos.append((l[21], int(l[22:26]), l[17:20].strip(),
                               coord(l)))
    log("receptor.pdb: %d atomos de proteina en cadenas %s"
        % (len(atomos), ", ".join(sorted({a[0] for a in atomos}))))

    # cadenas que tocan el bolsillo + la cadena del ligando
    dist = np.linalg.norm(np.array([a[3] for a in atomos]) - centro_lig,
                          axis=1)
    cadenas_bolsillo = sorted({atomos[i][0]
                               for i in np.where(dist < RADIO_BOLSILLO)[0]})
    log("cadenas con atomos a menos de %.1f A del ligando: %s"
        % (RADIO_BOLSILLO, ", ".join(cadenas_bolsillo) or "ninguna"))
    cadenas = set(cadenas_bolsillo) | {atomos[int(np.argmin(dist))][0]}
    log("cadenas conservadas: %s" % ", ".join(sorted(cadenas)))

    with open(PDB, encoding="utf-8", errors="ignore") as f:
        cadenas_pdb = {l[21] for l in f if l.startswith("ATOM")}
    if cadenas_pdb == cadenas:
        n = normalizar()
    else:
        # si sobrara alguna cadena, se filtra primero y se normaliza despues
        filtrado = PDB + ".filtrado"
        with open(PDB, encoding="utf-8", errors="ignore") as f, \
                open(filtrado, "w", encoding="utf-8") as g:
            for l in f:
                if l.startswith("ATOM") and l[21] in cadenas:
                    g.write(l)
        os.replace(filtrado, PDB + ".solo")
        globals()["PDB"] = PDB + ".solo"
        n = normalizar()
    log("receptor limpio: %s (%d atomos)" % (os.path.basename(LIMPIO), n))

    # residuos del bolsillo, para comparar con lo publicado
    nombres = {}
    for i in np.where(dist < RADIO_BOLSILLO)[0]:
        cad, num, res, _ = atomos[i]
        nombres["%s%d" % (res[0] + res[1:].lower(), num)] = dist[i]
    log("residuos a menos de %.1f A: %d" % (RADIO_BOLSILLO, len(nombres)))
    for etiqueta, num in RESIDUOS_CLAVE.items():
        log("   %-8s residuo %d: %s"
            % (etiqueta, num, "EN EL BOLSILLO" if etiqueta in nombres
               else "NO detectado"))

    if not os.path.exists(PDBQT):
        cmd = [MK_RECEPTOR, "--read_pdb", LIMPIO, "-o", PDBQT[:-6],
               "-p", "-j", "-a", "--default_altloc", "A",
               "--box_center", "%.3f" % centro_lig[0], "%.3f" % centro_lig[1],
               "%.3f" % centro_lig[2],
               "--box_size", str(TAMANO_CAJA), str(TAMANO_CAJA),
               str(TAMANO_CAJA)]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if not os.path.exists(PDBQT):
            log("FALLO preparando el receptor con meeko:\n%s"
                % ((r.stdout or "") + (r.stderr or ""))[-1200:])
            return 1
    n_rec = sum(1 for l in open(PDBQT, encoding="utf-8", errors="ignore")
                if l.startswith(("ATOM", "HETATM")))
    log("receptor preparado con meeko: %s (%d atomos)"
        % (os.path.basename(PDBQT), n_rec))

    json.dump({"centro_caja": [round(float(x), 2) for x in centro_lig],
               "tamano_caja": TAMANO_CAJA,
               "receptor_pdbqt": "receptor_cdk2.pdbqt",
               "pdb": "1h00 (DUD-E)",
               "ligando_cristal": "crystal_ligand.mol2 (54 atomos pesados)"},
              open(CAJA, "w", encoding="utf-8"), indent=1)
    log("caja escrita: %s" % os.path.basename(CAJA))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
