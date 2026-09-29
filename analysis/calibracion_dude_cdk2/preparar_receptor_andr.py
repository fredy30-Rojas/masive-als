#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara el receptor y la caja de ANDR (receptor de androgenos) como control.

POR QUE ESTA DIANA Y NO CDK2
-----------------------------
CDK2 resulto ser un control MAL ELEGIDO, no una prueba de que el embudo falle: su
ligando cristalino venia modelado DOS veces al 50 % (FAP y FCP, altloc A y B), con
B medio ~47, y el brazo del amonio no tocaba ni proteina ni agua. Un control asi no
puede decir nada del embudo.

ANDR se eligio cribando las 102 dianas de DUD-E con criterios medibles SIN acoplar
(`elegir_control_dude.py`) y luego con la pericia del cristal
(`verificar_cristal_control.py`). Gano en todo:

  - ligando TES, un ESTEROIDE (C19H30O2, 21 atomos pesados, 4 anillos);
  - **0 enlaces rotatorios**: no tiene un solo brazo que pueda girar. Es lo mas
    rgido que se puede pedir, y la razon por la que se descarto CDK2 era justo la
    flexibilidad del brazo del amonio;
  - una sola copia en el cristal, sin altloc y con ocupacion 1,00 (CDK2: dos al 50 %);
  - factor B medio 19,2, el mas bajo de las cuatro candidatas (CDK2: ~47);
  - los 2 polares tienen pareja de proteina a menos de 3,5 A: **0 de 2 sueltos**
    (CDK2: el amonio a 4,17 A de la proteina y a 5,90 A del agua mas cercana);
  - solo 1 de 21 atomos a mas de 4,5 A de cualquier cosa: el ligando esta
    realmente dentro del bolsillo, no asomando;
  - el mol2 de DUD-E ES el HETATM del PDB, con RMSD 0,0000 A al emparejar por
    elemento, asi que receptor y ligando son del mismo cristal.

QUE HACE
--------
La misma receta que `preparar_receptor_cdk2.py`, que a su vez es la de
`construir_receptor_tbk1.py`:

  1. El `receptor.pdb` de DUD-E viene en formato antiguo de 54 columnas, SIN columna
     de elemento, y meeko la exige. Se reescribe a formato moderno. Los residuos no
     estandar se detectan y se avisan, porque meeko se para en ellos.
  2. La caja es de 24 A centrada en el centroide del ligando del cristal: la
     convencion de DUD-E y la del proyecto.
  3. El receptor se prepara con meeko, como el de CDK2 y el de TBK1.

OJO CON EL ELEMENTO
-------------------
Igual que en `verificar_cristal_control.py`, el elemento se saca del NOMBRE del atomo
y hay que mirar DOS letras antes que una: `CL1` es cloro. En el receptor de andr hay
cisteinas y ningun cloro, pero el parser queda correcto por si acaso.

Salida: `receptor_andr.pdbqt`, `receptor_andr.pdb`, `caja_andr.json` y este log.

Uso:
    python preparar_receptor_andr.py
"""
from __future__ import annotations

import json
import os
import subprocess

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "_cache_dude")
PDB_ORIGEN = os.path.join(CACHE, "andr_receptor.pdb")
MOL2_ORIGEN = os.path.join(CACHE, "andr_mol2.mol2")
LIMPIO = os.path.join(BASE, "receptor_andr.pdb")
PDBQT = os.path.join(BASE, "receptor_andr.pdbqt")
CAJA = os.path.join(BASE, "caja_andr.json")
PYDIR = r"C:\Users\Fredy\AppData\Local\Python\pythoncore-3.14-64"
MK_RECEPTOR = os.path.join(PYDIR, "Scripts", "mk_prepare_receptor.exe")

TAMANO_CAJA = 24         # el mismo que TBK1 y CDK2
RADIO_BOLSILLO = 8.0
PDB_ORIG = "2AM9"
NOMBRE_LIGANDO = "TES"

# Las 20 aminoacidas estandar. Meeko se para con cualquier otra.
ESTANDAR = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    # variantes de histidina que meeko trae en sus plantillas
    "HID", "HIE", "HIP", "CYX", "CYM",
}


def log(m):
    print(m, flush=True)


def coord(linea):
    return np.array([float(linea[30:38]), float(linea[38:46]),
                     float(linea[46:54])])


def elemento_de(nombre):
    """Elemento desde el nombre del atomo (los PDB antiguos no traen columna).

    Dos letras antes que una: `CL1` es cloro, no carbono.
    """
    dos = nombre[:2].upper()
    if dos in ("CL", "BR", "SE", "NA", "MG", "ZN", "FE", "MN", "CA"):
        return dos.capitalize()
    return nombre[0].upper()


def ligandos_del_mol2():
    """Coords de los atomos PESADOS del cristal, leidas del bloque TRIPOS>ATOM.

    Se descartan los H explicitos porque el centro de la caja debe salir de los
    atomos que se van a acoplar, no de los que meeko va a regenerar.
    """
    pts = []
    dentro = False
    with open(MOL2_ORIGEN, encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("@<TRIPOS>ATOM"):
                dentro = True
                continue
            if l.startswith("@<TRIPOS>") and dentro:
                break
            if dentro:
                parts = l.split()
                if len(parts) >= 6:
                    if elemento_de(parts[1].split(".")[0]) in ("H", "D"):
                        continue
                    pts.append([float(parts[2]), float(parts[3]),
                                float(parts[4])])
    if not pts:
        raise SystemExit("no pude leer el ligando de %s" % MOL2_ORIGEN)
    return np.array(pts)


def normalizar():
    """Reescribe el PDB a formato moderno, quitando los H de nomenclatura vieja.

    Meeko exige la columna de elemento (77-78) y el fichero de DUD-E no la tiene.
    Los H (HN, HA, H1...) se quitan porque meeko los regenera desde la geometria.
    """
    n = hidrogenos = 0
    raros = {}
    with open(LIMPIO, "w", encoding="utf-8") as g:
        for l in open(PDB_ORIGEN, encoding="utf-8", errors="ignore"):
            if not l.startswith("ATOM"):
                continue
            res = l[17:20].strip().upper()
            if res not in ESTANDAR:
                raros[res] = raros.get(res, 0) + 1
            nombre = l[12:16].strip()
            elemento = elemento_de(nombre)
            if elemento in ("H", "D"):
                hidrogenos += 1
                continue
            linea = l[:54] + "  1.00  0.00          " + " %1s" % elemento
            g.write(linea[:78].ljust(78) + "\n")
            n += 1
        g.write("END\n")
    log("receptor normalizado: %d atomos pesados, %d hidrogenos fuera"
        % (n, hidrogenos))
    if raros:
        log("AVISO: residuos fuera de la tabla de meeko: %s"
            % ", ".join("%s (%d at.)" % (k, v) for k, v in sorted(raros.items())))
        log("       meeko puede parar en ellos; si lo hace, se avisa aqui.")
    else:
        log("residuos: todos en la tabla de meeko")
    return n


def main():
    for p in (PDB_ORIGEN, MOL2_ORIGEN):
        if not os.path.exists(p):
            log("FALTA %s: corre antes elegir_control_dude.py" % p)
            return 1

    lig = ligandos_del_mol2()
    centro = lig.mean(axis=0)
    log("ligando del cristal (%s): %d atomos pesados, centro %.2f, %.2f, %.2f"
        % (NOMBRE_LIGANDO, len(lig), *centro))

    atomos = []
    for l in open(PDB_ORIGEN, encoding="utf-8", errors="ignore"):
        if l.startswith("ATOM"):
            atomos.append((l[21], int(l[22:26]), l[17:20].strip(), coord(l)))
    log("receptor.pdb: %d atomos de proteina" % len(atomos))
    if not atomos:
        log("el receptor no tiene ATOM: algo va mal con la descarga")
        return 1

    dist = np.linalg.norm(np.array([a[3] for a in atomos]) - centro, axis=1)
    residuos = {}
    for i in np.where(dist < RADIO_BOLSILLO)[0]:
        cad, num, res, _ = atomos[i]
        etiqueta = "%s%d" % (res[0] + res[1:].lower(), num)
        residuos[etiqueta] = min(residuos.get(etiqueta, 1e9), dist[i])
    log("residuos a menos de %.1f A del ligando: %d"
        % (RADIO_BOLSILLO, len(residuos)))
    for etiqueta in sorted(residuos, key=lambda k: residuos[k])[:12]:
        log("   %-8s %.2f A" % (etiqueta, residuos[etiqueta]))

    n = normalizar()
    if n == 0:
        return 1

    if not os.path.exists(PDBQT):
        cmd = [MK_RECEPTOR, "--read_pdb", LIMPIO, "-o", PDBQT[:-6],
               "-p", "-j", "-a", "--default_altloc", "A",
               "--box_center", "%.3f" % centro[0], "%.3f" % centro[1],
               "%.3f" % centro[2],
               "--box_size", str(TAMANO_CAJA), str(TAMANO_CAJA),
               str(TAMANO_CAJA)]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if not os.path.exists(PDBQT):
            log("FALLO preparando el receptor con meeko:\n%s"
                % ((r.stdout or "") + (r.stderr or ""))[-1500:])
            return 1
        log("meeko: %s" % ((r.stdout or "").strip()[-300:] or "sin salida"))
    n_rec = sum(1 for l in open(PDBQT, encoding="utf-8", errors="ignore")
                if l.startswith(("ATOM", "HETATM")))
    log("receptor preparado con meeko: receptor_andr.pdbqt (%d atomos)" % n_rec)

    json.dump({"centro_caja": [round(float(x), 2) for x in centro],
               "tamano_caja": TAMANO_CAJA,
               "receptor_pdbqt": "receptor_andr.pdbqt",
               "pdb": "%s (DUD-E)" % PDB_ORIG,
               "diana": "andr",
               "ligando_cristal": "TES (%d atomos pesados, C19H30O2, "
                                  "0 enlaces rotatorios)" % len(lig),
               "residuos_bolsillo": sorted(residuos)},
              open(CAJA, "w", encoding="utf-8"), indent=1)
    log("caja escrita: caja_andr.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
