#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""XL20 frente a XL23: la pareja de la misma familia, medida y no contada de memoria.

POR QUE EXISTE
--------------
El articulo del CR (Gao et al. 2026, Nature Aging 6:1667) deja **un solo compuesto
con union medida** —XL20— y varios companeros de tabla con actividad funcional.
De todos ellos, uno solo comparte el quimiotipo de XL20: **XL23**, adenina unida a
un aminociclohexanol y ademas cinco anillos. Los dos son, por tanto, la pareja
natural para mirar **relacion estructura-actividad dentro de la misma familia sin
salir de la fuente**: mismo esqueleto de union (la adenina y el aminociclohexanol),
distinta cola.

Este script escribe esa comparacion como numeros reproducibles, leyendo los
SMILES de donde ya viven (nada de copiarlos a mano aqui):

  * XL20 -> analysis/controles_tdp43_xl20.csv
  * XL23 -> analysis/suplementario_tdp43_xl21_27.csv

QUE MIDE
--------
1. Propiedades calculadas de los dos (formula, masa, logP, TPSA, donantes,
   aceptores, enlaces rotables, anillos, anillos aromaticos, centroesteros).
2. Parecido 2D: **Tanimoto ECFP4** (radio 2, 2048 bits), el mismo descriptor con
   **estereoquimica** y **claves MACCS**.
3. **Esqueleto de Murcko**: si coinciden o no (coincidir de esqueleto es lo que
   hace que sean familia, no solo parecerse).
4. **Nucleo comun (MCS)** con anillos completos: cuantos atomos y enlaces
   comparten, cual es ese nucleo en SMARTS, y si ese nucleo contiene la adenina.
5. **Lo que cada uno pone de su parte**: los atomos que quedan FUERA del nucleo
   comun, con su formula; es decir, la cola que diferencia a XL20 de XL23.
6. Comprobaciones de subestructura explicitas (la adenina en los dos, el
   esqueleto entero de XL20 dentro de XL23, etc.).

Uso:
    python analysis/comparar_xl20_xl23.py

Salida en analysis/_xl20_xl23/:
    pareja.csv                  una fila por propiedad, con los dos valores
    pareja.txt                  lo mismo, legible
    XL20_XL23_comparacion.png   los dos dibujos con el nucleo comun resaltado
"""
from __future__ import annotations

import csv
import os
import sys

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, Draw, rdFMCS, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
SALIDA = os.path.join(BASE, "_xl20_xl23")
CSV_XL20 = os.path.join(BASE, "controles_tdp43_xl20.csv")
CSV_XL23 = os.path.join(BASE, "suplementario_tdp43_xl21_27.csv")

TANIMOTO_ECFP4 = 2          # radio, el mismo descriptor ECFP4 que usa el proyecto
N_BITS = 2048

# La adenina como consulta de subestructura: purina con el amino en la posicion 6.
# Los atomos aromaticos en minuscula de SMARTS admiten sustituyentes encima, que es
# justo lo que hace falta: en XL20 el N9 lleva el ciclohexano y en XL23 tambien.
ADENINA = "Nc1ncnc2c1ncn2"


def log(m):
    print(m, flush=True)


def leer_smiles(ruta, ligando):
    """El SMILES de un ligando, leido de su CSV (sin copiarlo a mano)."""
    with open(ruta, encoding="utf-8-sig", newline="") as f:
        for fila in csv.DictReader(f):
            if (fila.get("ligand") or "").strip().upper() == ligando:
                return (fila.get("smiles") or "").strip()
    raise SystemExit("no encontrado %s en %s" % (ligando, ruta))


def uno(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise SystemExit("SMILES invalido: %s" % smi)
    return m


def properties(mol, etiqueta):
    estereo = Chem.FindMolChiralCenters(mol, includeUnassigned=True,
                                        useLegacyImplementation=False)
    definidos = sum(1 for _, c in estereo if c != "?")
    return {
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "masa": round(Descriptors.MolWt(mol), 1),
        "atomos_pesados": mol.GetNumHeavyAtoms(),
        "logp": round(Descriptors.MolLogP(mol), 2),
        "tpsa": round(Descriptors.TPSA(mol), 1),
        "donantes_h": Descriptors.NumHDonors(mol),
        "aceptores_h": Descriptors.NumHAcceptors(mol),
        "rotables": Descriptors.NumRotatableBonds(mol),
        "anillos": Descriptors.RingCount(mol),
        "anillos_aromaticos": rdMolDescriptors.CalcNumAromaticRings(mol),
        "centroesteros": len(estereo),
        "centroesteros_definidos": definidos,
        "heteroatomos": rdMolDescriptors.CalcNumHeteroatoms(mol),
        "etiqueta": etiqueta,
    }


def tanimoto(a, b):
    """Parecido 2D: ECFP4, ECFP4 con estereoquimica y claves MACCS."""
    out = {}
    out["ecfp4"] = round(DataStructs.TanimotoSimilarity(
        AllChem.GetMorganFingerprintAsBitVect(a, TANIMOTO_ECFP4, N_BITS),
        AllChem.GetMorganFingerprintAsBitVect(b, TANIMOTO_ECFP4, N_BITS)), 3)
    out["ecfp4_con_estereo"] = round(DataStructs.TanimotoSimilarity(
        rdMolDescriptors.GetMorganFingerprint(a, TANIMOTO_ECFP4, useChirality=True),
        rdMolDescriptors.GetMorganFingerprint(b, TANIMOTO_ECFP4, useChirality=True)), 3)
    out["maccs"] = round(DataStructs.TanimotoSimilarity(
        rdMolDescriptors.GetMACCSKeysFingerprint(a),
        rdMolDescriptors.GetMACCSKeysFingerprint(b)), 3)
    return out


def nucleo_comun(a, b):
    """Nucleo maximo comun con anillos completos y correspondencia de elementos y enlaces."""
    mcs = rdFMCS.FindMCS([a, b], atomCompare=rdFMCS.AtomCompare.CompareElements,
                        bondCompare=rdFMCS.BondCompare.CompareOrder,
                        ringMatchesRingOnly=True, completeRingsOnly=True,
                        timeout=30)
    return mcs


def submol(mol, indices):
    """Molecula REAL (no una consulta) con los atomos indicados de otra.

    El nucleo comun sale del MCS en forma de consulta SMARTS (`[#6]`, `[#7]`...), y
    una consulta no se puede comparar con un SMARTS de subestructura: pedirle si
    contiene la adenina da falso. Con los atomos que ese nucleo empareja en XL20 se
    reconstruye la molecula de verdad, aromaticidad y ordenes de enlace incluidos, y
    entonces las comprobaciones si dicen algo.
    """
    rw = Chem.RWMol(mol)
    for idx in sorted(set(range(mol.GetNumAtoms())) - set(indices), reverse=True):
        rw.RemoveAtom(idx)
    trozo = rw.GetMol()
    Chem.SanitizeMol(trozo, catchErrors=True)
    return trozo


def ramas(mol, nucleo):
    """Lo que el ligando pone de su parte: lo que queda fuera del nucleo comun.

    Se corta el nucleo con `ReplaceCore`, que deja el punto de corte marcado con un
    atomo comodin, para que se vea por donde cuelga cada cola y no una formula que
    depende de como RDKit rellene las valencias truncadas.
    """
    try:
        trozo = Chem.ReplaceCore(mol, nucleo, labelByIndex=True)
    except Exception:  # noqa: BLE001
        trozo = None
    if trozo is None:
        return 0, "(no se pudo cortar)"
    piezas = []
    for f in Chem.GetMolFrags(trozo, asMols=True, sanitizeFrags=False):
        pesados = sum(1 for at in f.GetAtoms() if at.GetAtomicNum() > 1)
        piezas.append((pesados, Chem.MolToSmiles(f)))
    piezas.sort(key=lambda p: -p[0])
    return (sum(p[0] for p in piezas),
            "; ".join("%s (%d atomos pesados)" % (s, n) for n, s in piezas))


def dibujar(a, b, nucleo, destino):
    """Los dos con el nucleo comun resaltado, para poder mirarlo a ojo."""
    try:
        coincidencias = [m.GetSubstructMatch(nucleo) for m in (a, b)]
        imagen = Draw.MolsToGridImage(
            [a, b], molsPerRow=2, subImgSize=(460, 340),
            legends=["XL20 (union medida)", "XL23 (sin union medida)"],
            highlightAtomLists=coincidencias)
        imagen.save(destino)
        return True
    except Exception as e:  # noqa: BLE001
        log("   (no se pudo dibujar: %s)" % e)
        return False


def main():
    os.makedirs(SALIDA, exist_ok=True)
    smi20 = leer_smiles(CSV_XL20, "XL20")
    smi23 = leer_smiles(CSV_XL23, "XL23")
    x20, x23 = uno(smi20), uno(smi23)

    p20, p23 = properties(x20, "XL20"), properties(x23, "XL23")
    tan = tanimoto(x20, x23)
    mcs = nucleo_comun(x20, x23)
    nucleo_consulta = Chem.MolFromSmarts(mcs.smartsString) if mcs.numAtoms else None
    # el nucleo como molecula real, construida con los atomos que empareja en XL20
    nucleo = (submol(x20, x20.GetSubstructMatch(nucleo_consulta))
              if nucleo_consulta is not None else None)

    esc20 = MurckoScaffold.MurckoScaffoldSmiles(mol=x20)
    esc23 = MurckoScaffold.MurckoScaffoldSmiles(mol=x23)

    adelina = Chem.MolFromSmarts(ADENINA)
    comprobaciones = {
        "esqueleto_murcko_igual": esc20 == esc23,
        "adenina_en_XL20": x20.HasSubstructMatch(adelina),
        "adenina_en_XL23": x23.HasSubstructMatch(adelina),
        "esqueleto_de_XL20_dentro_de_XL23": x23.HasSubstructMatch(x20),
        "esqueleto_de_XL23_dentro_de_XL20": x20.HasSubstructMatch(x23),
        "nucleo_en_los_dos": bool(nucleo) and x20.HasSubstructMatch(nucleo)
                              and x23.HasSubstructMatch(nucleo),
        "nucleo_contiene_la_adenina": bool(nucleo)
                                      and nucleo.HasSubstructMatch(adelina),
    }
    fuera20 = ramas(x20, nucleo) if nucleo is not None else (0, "")
    fuera23 = ramas(x23, nucleo) if nucleo is not None else (0, "")

    filas = []
    for clave in ("formula", "masa", "atomos_pesados", "logp", "tpsa", "donantes_h",
                  "aceptores_h", "rotables", "anillos", "anillos_aromaticos",
                  "centroesteros", "centroesteros_definidos", "heteroatomos"):
        filas.append(("propiedad:" + clave, str(p20[clave]), str(p23[clave])))
    for clave, valor in tan.items():
        filas.append(("tanimoto:" + clave, "%.3f" % valor, ""))
    filas.append(("esqueleto_murcko_XL20", esc20, ""))
    filas.append(("esqueleto_murcko_XL23", esc23, ""))
    filas.append(("nucleo_comun_smarts", mcs.smartsString, ""))
    filas.append(("nucleo_comun_tamano", "%d atomos, %d enlaces"
                  % (mcs.numAtoms, mcs.numBonds), ""))
    filas.append(("fuera_del_nucleo_XL20", "%d atomos: %s" % fuera20, ""))
    filas.append(("fuera_del_nucleo_XL23", "%d atomos: %s" % fuera23, ""))
    for clave, valor in comprobaciones.items():
        filas.append((clave, str(valor), ""))

    with open(os.path.join(SALIDA, "pareja.csv"), "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["medida", "XL20", "XL23"])
        w.writerows(filas)

    lineas = []
    lineas.append("XL20 frente a XL23 (misma familia: adenina-aminociclohexanol)")
    lineas.append("")
    lineas.append("  XL20 = %s" % smi20)
    lineas.append("  XL23 = %s" % smi23)
    lineas.append("")
    lineas.append("  %-34s %-28s %s" % ("propiedad", "XL20", "XL23"))
    for clave in ("formula", "masa", "atomos_pesados", "logp", "tpsa", "donantes_h",
                  "aceptores_h", "rotables", "anillos", "anillos_aromaticos",
                  "centroesteros_definidos", "heteroatomos"):
        lineas.append("  %-34s %-28s %s" % (clave, p20[clave], p23[clave]))
    lineas.append("")
    lineas.append("  parecido 2D (Tanimoto):")
    for clave, valor in tan.items():
        lineas.append("    %-30s %.3f" % (clave, valor))
    lineas.append("")
    lineas.append("  esqueleto de Murcko igual: %s" % comprobaciones["esqueleto_murcko_igual"])
    lineas.append("  nucleo comun: %d atomos, %d enlaces"
                  % (mcs.numAtoms, mcs.numBonds))
    lineas.append("     %s" % mcs.smartsString)
    lineas.append("  el nucleo comun contiene la adenina: %s"
                  % comprobaciones["nucleo_contiene_la_adenina"])
    lineas.append("")
    lineas.append("  lo que pone de su parte cada uno (fuera del nucleo comun):")
    lineas.append("    XL20: %d atomos: %s" % fuera20)
    lineas.append("    XL23: %d atomos: %s" % fuera23)
    lineas.append("")
    lineas.append("  comprobaciones de subestructura:")
    for clave, valor in comprobaciones.items():
        lineas.append("    %-42s %s" % (clave, valor))
    dibujar(x20, x23, nucleo, os.path.join(SALIDA, "XL20_XL23_comparacion.png"))
    lineas.append("")
    lineas.append("  dibujo con el nucleo comun resaltado: XL20_XL23_comparacion.png")
    texto = "\n".join(lineas)
    with open(os.path.join(SALIDA, "pareja.txt"), "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    log(texto)
    return 0


if __name__ == "__main__":
    sys.exit(main())
