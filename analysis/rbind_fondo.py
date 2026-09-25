#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""rbind_fondo.py — el fondo duro de unidores de ARN (R-BIND 2.0).

POR QUE ESTE FONDO Y NO OTRO
----------------------------
El fondo de las validaciones son 122 señuelos **emparejados en propiedades**: se
parecen a los positivos en tamaño, logP y demás, pero nadie ha demostrado que se
unan a nada. Eso mide si el embudo reconoce *alguna* quimia.

Lo que NO mide es lo que importa en TDP-43: el plan lo dice desde el 20 de
septiembre (`PLAN_TDP43_2026-09-20.md`) y el paper lo repite en su limitación 2
(`paper/paper_masive_als.md`, línea 216). Varios de los positivos mejor
documentados de TDP-43 son cationes planos que se unen a ácido nucleico; el
fondo que decide es un conjunto de **moléculas pequeñas que sí se unen a ARN de
verdad**. Si el embudo no separa un ligando de TDP-43 de un unidor genérico de
ARN, no ha demostrado nada.

De donde sale
-------------
R-BIND 2.0 (Donlic et al., ACS Chem Biol 2022, 17:1556; DOI
10.1021/acschembio.2c00224) es la base curada de ligandos bioactivos de ARN no
ribosómico. Son 188 miembros: 153 moléculas pequeñas (SM, <= 700 Da) y 35
grandes (LM). El artículo publica los datos **solo** en su material
suplementario, y la web que los servía (`rbind.chem.duke.edu`) está en
transición desde que el laboratorio Hargrove se mudó a Toronto
(`hargrovelab.org/rbind`). El ACS devuelve 403 a la descarga automática.

La copia buena y verificable es el ZIP que servía la propia web del proyecto,
que quedó archivado:

    https://web.archive.org/web/20220524151659id_/
        https://rbind.chem.duke.edu/data/RBIND_V2.0_ALL_CSV.zip

Se descarga y se descomprime en `analysis/_rbind/`. Trae los dos ficheros SM
(A = una fila por ligando con SMILES y propiedades; B = una fila por ensayo).

Que se hace con ellos
---------------------
  * se toma `SMILES (NC)` (la forma sin carga neta de la tabla; la columna
    `SMILES (Co)` trae estados de protonación que RDKit lee con cargas de hasta
    +6 y no sirve para acoplar). Si `(NC)` no parsea, se cae a `(Co)`;
  * se quita todo lo que ya está en `verdad_de_referencia.csv` (cualquier diana,
    `apto` = si o no): el fondo no puede contener positivos ni controles del
    propio banco, o la medida se contamina;
  * se deduplica por SMILES canónico (la tabla tiene pares A/B de estereoisómeros
    y entradas repetidas);
  * se guardan los SMILES preparables y los descartes con su motivo.

Salida: `_rbind/rbind_fondo_sm.csv` (id, nombre, smiles, mw, pesados, estado) y
el contador en pantalla. Los que entran van a acoplar con `acoplar_fondo_rbind.py`.

Uso:
    python rbind_fondo.py
"""
import csv
import os
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RBIND = os.path.join(BASE, "_rbind")
A_SM = os.path.join(RBIND, "RBIND_v2.0_A_SM.csv")
VERDAD = os.path.join(BASE, "verdad_de_referencia.csv")
SALIDA = os.path.join(RBIND, "rbind_fondo_sm.csv")

URL_ZIP = ("https://web.archive.org/web/20220524151659id_/"
           "https://rbind.chem.duke.edu/data/RBIND_V2.0_ALL_CSV.zip")


def inchi_clave(smi):
    """Primer bloque del InChIKey: la identidad sin estereoquimica ni protones."""
    mol = Chem.MolFromSmiles(smi) if smi else None
    if mol is None:
        return None
    try:
        return Chem.MolToInchiKey(mol).split("-")[0]
    except Exception:
        return None


def cargar_del_banco():
    """Claves del banco de verdad: el fondo no puede llevar nada de aqui."""
    claves = {}
    for r in csv.DictReader(open(VERDAD, encoding="utf-8")):
        k = inchi_clave(r["smiles"])
        if k:
            claves.setdefault(k, []).append(r["target"] + "/" + r["ligand"])
    return claves


def cargar_rbind():
    if not os.path.exists(A_SM):
        sys.exit("Falta %s. Descargar el ZIP de %s y descomprimir en %s"
                 % (A_SM, URL_ZIP, RBIND))
    filas = list(csv.reader(open(A_SM, encoding="utf-8", errors="replace")))
    hdr = filas[1]
    i_id = hdr.index("Database ID")
    i_nom = hdr.index("Name")
    i_nc = hdr.index("SMILES (NC)")
    i_co = hdr.index("SMILES (Co)")
    i_mw = hdr.index("MW")
    out = []
    for r in filas[2:]:
        if len(r) <= max(i_id, i_nc, i_co) or not r[i_id].strip():
            continue
        out.append({"id": r[i_id].strip(), "nombre": r[i_nom].strip(),
                    "nc": r[i_nc].strip(), "co": r[i_co].strip(),
                    "mw": r[i_mw].strip()})
    return out


def main():
    del_banco = cargar_del_banco()
    print("banco de verdad: %d identidades" % len(del_banco))

    entradas = cargar_rbind()
    print("R-BIND 2.0 (SM), filas con ID: %d" % len(entradas))

    guardadas, vistos, salida = {}, {}, []
    for e in entradas:
        mol = Chem.MolFromSmiles(e["nc"])
        origen = "NC"
        if mol is None:
            mol = Chem.MolFromSmiles(e["co"])
            origen = "Co"
        if mol is None:
            salida.append((e["id"], e["nombre"], "", e["mw"], "", "sin SMILES valido"))
            continue
        smi = Chem.MolToSmiles(mol)
        k = inchi_clave(smi)
        n = mol.GetNumHeavyAtoms()

        if k in del_banco:
            salida.append((e["id"], e["nombre"], smi, e["mw"], n,
                           "en el banco de verdad: " + ", ".join(del_banco[k])))
            continue
        if k in vistos:
            salida.append((e["id"], e["nombre"], smi, e["mw"], n,
                           "duplicado de " + vistos[k]))
            continue
        vistos[k] = e["id"]
        guardadas[e["id"]] = smi
        salida.append((e["id"], e["nombre"], smi, e["mw"], n, "entra (%s)" % origen))

    entra = [s for s in salida if s[5].startswith("entra")]
    fuera = [s for s in salida if not s[5].startswith("entra")]
    print("")
    print("entran en el fondo duro : %d" % len(entra))
    print("descartados             : %d" % len(fuera))
    for s in fuera:
        print("    %-20s %-28s %s" % (s[0], s[1][:28], s[5]))
    entradas_solas = sorted(set(s[4] for s in entra))
    print("")
    print("atomos pesados: min %d  mediana %d  max %d"
          % (entradas_solas[0], entradas_solas[len(entradas_solas) // 2],
             entradas_solas[-1]))

    with open(SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["rbind_id", "nombre", "smiles", "mw", "pesados", "estado"])
        for s in salida:
            w.writerow(s)
    print("guardado: %s" % SALIDA)
    return 0


if __name__ == "__main__":
    sys.exit(main())
