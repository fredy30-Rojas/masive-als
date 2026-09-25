#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Segunda ronda del control de redocking: quimias nuevas (21 sep 2026).

QUE MIDE Y POR QUE
------------------
`redock_trp32.py` mide el RMSD del modo MEJOR PUNTUADO, que es el criterio
estandar. Pero cuando ese modo falla hay que distinguir dos cosas que se
confunden todo el rato en docking:

  * fallo de MUESTREO: la pose cristalina no aparece en ninguna de las poses
    generadas. El buscador no la encuentra.
  * fallo de PUNTUACION: la pose cristalina SI esta entre las generadas, pero la
    funcion prefiere otra. El buscador la encuentra y la descarta.

En el control de las catecolaminas el proyecto ya vio los dos casos (la dopamina
cae en el segundo; la 5-fluorouridina, en el primero). Este script repite la
misma medida para las siete quimias nuevas, leyendo las poses ya calculadas: no
vuelve a acoplar nada.

Ademas comprueba de que tamano es el ligando del cristal frente al ligando ideal
de RCSB, porque un ligando a medio resolver no se puede comparar con el ideal y
el RMSD sale «no calculable» (es el caso de 6A9O).

Uso: python analizar_modos_quimias.py
Salida: modos_quimias_nuevas.csv y modos_quimias_nuevas.log
"""
import csv
import glob
import json
import os
import re
import sys

from rdkit import Chem

import redock_trp32 as R

BASE = R.BASE


def rmsd_sobre_resueltos(lig_pdb, pose, conf_id=0):
    """RMSD cuando el ligando del cristal esta a medio resolver (caso de 6A9O).

    En 6A9O el deposito solo modela 24 de los 33 atomos pesados del Lig9, asi que
    la comparacion atomo a atomo no existe y el RMSD sale «no calculable». Aqui se
    empareja el cristal como SUBCONJUNTO de la pose acoplada y se mide el RMSD
    solo sobre los atomos que el cristal si tiene. No es el mismo numero que un
    redocking limpio (mide menos atomos, y los que faltan son justo la cola
    flexible), pero dice si el esqueleto resuelto esta en su sitio.
    """
    cry = Chem.MolFromPDBBlock(open(lig_pdb).read(), removeHs=False,
                              proximityBonding=True)
    if cry is None or pose is None:
        return None
    c = R.aplanar(Chem.RemoveHs(cry))
    p = R.aplanar(Chem.RemoveHs(Chem.Mol(pose)))
    if c.GetNumAtoms() >= p.GetNumAtoms():
        return None
    match = p.GetSubstructMatch(c, useChirality=False)
    if not match:
        return None
    cp = p.GetConformer(conf_id).GetPositions()
    cc = c.GetConformer().GetPositions()
    d2 = 0.0
    for i, j in enumerate(match):
        d2 += float(((cp[j] - cc[i]) ** 2).sum())
    return (d2 / len(match)) ** 0.5
CRISTAL_POR_ENTRADA = {
    "4A7Q": "4A7Q_cristal.pdb", "4A7G": "4A7G_cristal.pdb",
    "2WZ6": "2WZ6_cristal.pdb", "2WZ0": "2WZ0_cristal.pdb",
    "8GSQ": "8GSQ_cristal.pdb", "6A9O": "6A9O_cristal.pdb",
    "5YTO": "5YTO_cristal.pdb",
    # ronda 1, para comparar en la misma tabla
    "4A7S": "4A7S_cristal.pdb", "4A7T": "4A7T_cristal.pdb",
    "4A7U": "4A7U_cristal.pdb", "4A7V": "4A7V_cristal.pdb",
}


def afinidades_de_pdbqt(path):
    """Afinidades de todos los modos, en el orden del fichero."""
    out = []
    for l in open(path, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            m = re.search(r"([-+]?\d+\.\d+)", l)
            if m:
                out.append(float(m.group(1)))
    return out


def analizar(entrada):
    """Todas las poses de todas las semillas, ordenadas por afinidad."""
    lig = os.path.join(BASE, "ligands", CRISTAL_POR_ENTRADA[entrada])
    if not os.path.exists(lig):
        return None
    at_cristal = sum(1 for l in open(lig) if l.startswith(("ATOM", "HETATM")))
    filas = []
    for f in sorted(glob.glob(os.path.join(BASE, "out", "%s_caja24_e8_s*.pdbqt" % entrada))):
        semilla = re.search(r"_s(\d+)\.pdbqt$", f)
        semilla = semilla.group(1) if semilla else "?"
        aps = afinidades_de_pdbqt(f)
        mol = R.molecula_dockeada(f)
        if mol is None:
            continue
        n_at_pose = mol.GetNumAtoms()
        for i in range(mol.GetNumConformers()):
            val = R.rmsd_en_sitio(lig, mol, conf_id=i)
            parcial = False
            if val is None and at_cristal != n_at_pose:
                val = rmsd_sobre_resueltos(lig, mol, conf_id=i)
                parcial = val is not None
            filas.append({
                "entrada": entrada, "semilla": semilla, "modo": i + 1,
                "afinidad": aps[i] if i < len(aps) else None,
                "rmsd": round(val, 2) if val is not None else None,
                "parcial": parcial,
                "atomos_cristal": at_cristal, "atomos_pose": n_at_pose,
            })
    return filas


def main():
    log_path = os.path.join(BASE, "modos_quimias_nuevas.log")
    lineas = []
    todas = []
    for entrada in CRISTAL_POR_ENTRADA:
        filas = analizar(entrada)
        if not filas:
            lineas.append("%-5s sin poses calculadas" % entrada)
            continue
        todas.extend(filas)
        con_rmsd = [f for f in filas if f["rmsd"] is not None]
        bajo2 = [f for f in con_rmsd if f["rmsd"] <= 2.0]
        mejor = min(con_rmsd, key=lambda f: f["rmsd"]) if con_rmsd else None
        top = min((f for f in filas if f["afinidad"] is not None),
                  key=lambda f: f["afinidad"], default=None)
        # ¿el modo mejor puntuado es el que coincide con el cristal?
        coincide_top = (top is not None and top["rmsd"] is not None
                        and top["rmsd"] <= 2.0)
        lineas.append(
            "%-5s atomos cristal=%2d pose=%2d | poses=%2d | modos<=2A: %2d | "
            "mejor RMSD %.2f A%s (modo %s, %.2f kcal/mol) | mejor puntuado: modo %s "
            "%.2f kcal/mol RMSD %s -> %s"
            % (entrada, filas[0]["atomos_cristal"], filas[0]["atomos_pose"],
               len(filas), len(bajo2),
               mejor["rmsd"] if mejor else float("nan"),
               " (solo atomos resueltos)" if mejor and mejor.get("parcial") else "",
               mejor["modo"] if mejor else "-",
               mejor["afinidad"] if mejor and mejor["afinidad"] is not None else float("nan"),
               top["modo"] if top else "-",
               top["afinidad"] if top and top["afinidad"] is not None else float("nan"),
               ("%.2f A" % top["rmsd"]) if top and top["rmsd"] is not None else "n/d",
               "RECUPERA" if bajo2 else "NO RECUPERA"))

    texto = "\n".join(lineas)
    print(texto)
    with open(log_path, "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    csv_path = os.path.join(BASE, "modos_quimias_nuevas.csv")
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["entrada", "semilla", "modo", "afinidad",
                                          "rmsd", "parcial", "atomos_cristal",
                                          "atomos_pose"])
        w.writeheader()
        for fila in todas:
            w.writerow(fila)
    print("CSV: %s" % csv_path)
    return 0


if __name__ == "__main__":
    sys.exit(main())
