#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Acopla XL20 y XL23 en el CR del C-terminal de TDP-43 (primer acoplamiento).

QUE HACE
--------
1. Prepara los dos ligandos con **la receta canonica del proyecto**
   (`preparar_ligando.escribir`, con su comprobacion de pseudo-atomos de pegamento),
   leyendo los SMILES de donde ya viven en el repositorio, no de memoria:
   XL20 de `controles_tdp43_xl20.csv`, XL23 de `suplementario_tdp43_xl21_27.csv`.
2. Los acopla con Vina en **cada modelo** del conjunto de RMN del CR que eligio
   `construir_receptor_cr.py` (2N2C, 6 modelos, caja de 22 A sobre el anillo del
   Trp334 de cada modelo). Mismos parametros que el resto del proyecto:
   exhaustividad 8, 9 modos, semilla 42.
3. Mide, para cada pose, **cuantos atomos del receptor tiene a menos de 4 A** y
   cuantos de sus propios atomos estan en contacto, ademas de si toca el Trp334 y a
   que distancia esta del anillo. Sin eso no se puede distinguir una pose apoyada en
   la proteina de otra flotando en el hueco de la caja, que Vina tambien puntua bien.

LO QUE NO HACE
--------------
No llama a esto «validacion» ni saca una conclusion de dos ligandos. Dos compuestos
no son una relacion estructura-actividad: son un par. El informe dice lo que se puede
afirmar (donde cae cada uno, con que contactos, si el sitio del Trp334 los admite) y
lo que no (que XL23 no se una: eso no esta medido, y un acoplamiento no lo mide).

Uso:
    python analysis/acoplar_xl20_xl23_cr.py [--exhaustividad 8] [--lado 22]

Salida en analysis/_cr_receptor/:
    ligands/XL20.pdbqt, ligands/XL23.pdbqt
    out/XL20_modelo1.pdbqt ...            poses por ligando y modelo
    acoplamiento_cr.csv                   una fila por pose
    acoplamiento_cr.txt                   el resumen legible
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
SALIDA = os.path.join(BASE, "_cr_receptor")
LIGANDS = os.path.join(SALIDA, "ligands")
OUT = os.path.join(SALIDA, "out")
MODELOS = os.path.join(SALIDA, "modelos")

# La receta de preparacion vive en UN solo sitio: ../preparar_ligando.py
sys.path.insert(0, os.path.dirname(BASE))
sys.path.insert(0, os.path.join(BASE, "redocking_trp32"))
from preparar_ligando import escribir as escribir_ligando  # noqa: E402
import redock_trp32 as R  # noqa: E402

CORTE = 4.0            # A, el mismo corte de contactos que usa el resto del proyecto
CR = (320, 340)
TRP = 334
ANILLO_TRP = ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2")

FUENTES = {
    "XL20": ("controles_tdp43_xl20.csv", "XL20", "union medida (SPR + CETSA)"),
    "XL23": ("suplementario_tdp43_xl21_27.csv", "XL23", "sin union medida"),
}


def log(m):
    print(m, flush=True)


def leer_smiles(archivo, ligando):
    with open(os.path.join(BASE, archivo), encoding="utf-8-sig", newline="") as f:
        for fila in csv.DictReader(f):
            if (fila.get("ligand") or "").strip().upper() == ligando:
                return (fila.get("smiles") or "").strip()
    raise SystemExit("no encontrado %s en %s" % (ligando, archivo))


def atomos(pdb, solo_pesados=True):
    """(residuo, nombre, elemento, xyz) de un PDB o PDBQT."""
    out = []
    for l in open(pdb, encoding="utf-8", errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        nombre = l[12:16].strip()
        elemento = (l[76:78].strip() or nombre[0]).upper()
        if solo_pesados and (elemento == "H" or nombre.startswith("H")):
            continue
        out.append((int(l[22:26]), nombre, elemento,
                    np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


def poses_del_pdbqt(ruta):
    """[(modo, afinidad, [(elemento, xyz)])] de la salida de Vina."""
    poses, actual, af = [], None, None
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            actual, af = [], None
        elif l.startswith("REMARK VINA RESULT:"):
            m = re.search(r"(-?\d+\.\d+)", l)
            af = float(m.group(1)) if m else None
        elif l.startswith(("ATOM", "HETATM")) and actual is not None:
            nombre = l[12:16].strip()
            elemento = (l[76:78].strip() or nombre[0]).upper()
            if elemento == "H" or nombre.startswith("H"):
                continue
            actual.append((elemento, np.array([float(l[30:38]), float(l[38:46]),
                                               float(l[46:54])])))
        elif l.startswith("ENDMDL"):
            if actual:
                poses.append((len(poses) + 1, af, actual))
            actual = None
    if actual:
        poses.append((len(poses) + 1, af, actual))
    return poses


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustividad", type=int, default=8)
    ap.add_argument("--semilla", type=int, default=42)
    ap.add_argument("--lado", type=int, default=22)
    args = ap.parse_args()

    os.makedirs(LIGANDS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(SALIDA, "caja.json"), encoding="utf-8") as f:
        caja = json.load(f)

    log("Receptor: %s (%s) | caja de %d A sobre el anillo del Trp334 de cada modelo"
        % (caja["pdb"], caja["residuos"], args.lado))
    log("")

    # --- ligandos: la receta canonica, con su comprobacion de pseudo-atomos ---
    pdbqt_ligando = {}
    for nombre, (archivo, clave, evidencia) in FUENTES.items():
        smi = leer_smiles(archivo, clave)
        ruta = escribir_ligando(nombre, smi, LIGANDS)
        if ruta is None:
            log("%s NO se pudo preparar" % nombre)
            return 1
        pdbqt_ligando[nombre] = ruta
        log("  %-5s %-46s %s" % (nombre, smi, evidencia))
    log("")

    filas = []
    for nombre, ligando in pdbqt_ligando.items():
        for clave in sorted(caja["modelos"], key=int):
            info = caja["modelos"][clave]
            receptor = os.path.join(BASE, info["pdbqt"])
            receptor_pdb = os.path.join(MODELOS, "cr_modelo%s.pdb" % clave)
            out = os.path.join(OUT, "%s_modelo%s.pdbqt" % (nombre, clave))
            cx, cy, cz = info["centro"]
            cmd = [R.VINA, "--receptor", receptor, "--ligand", ligando,
                   "--center_x", "%.3f" % cx, "--center_y", "%.3f" % cy,
                   "--center_z", "%.3f" % cz,
                   "--size_x", str(args.lado), "--size_y", str(args.lado),
                   "--size_z", str(args.lado),
                   "--exhaustiveness", str(args.exhaustividad),
                   "--num_modes", "9", "--seed", str(args.semilla), "--out", out]
            if not (os.path.exists(out) and os.path.getsize(out) > 100):
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
                if not (os.path.exists(out) and os.path.getsize(out) > 100):
                    log("  %s modelo %s: FALLO Vina: %s"
                        % (nombre, clave, (r.stderr or r.stdout)[-200:]))
                    continue

            rec = atomos(receptor_pdb)
            anillo = np.array([xyz for (res, nom, el, xyz) in rec
                               if res == TRP and nom in ANILLO_TRP])
            centro_anillo = anillo.mean(axis=0) if len(anillo) else None
            rec_xyz = np.array([a[3] for a in rec])
            for modo, afinidad, lig in poses_del_pdbqt(out):
                L = np.array([a[1] for a in lig])
                d = np.sqrt(((L[:, None, :] - rec_xyz[None, :, :]) ** 2).sum(-1))
                cerca = d.min(axis=0) < CORTE
                residuos = sorted({rec[i][0] for i in np.nonzero(cerca)[0]})
                en_contacto = int((d.min(axis=1) < CORTE).sum())
                toca_trp = TRP in residuos
                # ojo con np.linalg.norm sin eje: sobre una matriz devuelve la norma
                # de Frobenius (un solo numero), no la distancia de cada atomo.
                dist_centro = dist_anillo_atomo = None
                if centro_anillo is not None:
                    dist_centro = float(np.linalg.norm(L - centro_anillo, axis=1).min())
                    dist_anillo_atomo = float(np.linalg.norm(
                        L[:, None, :] - anillo[None, :, :], axis=-1).min())
                del_cr = [r for r in residuos if CR[0] <= r <= CR[1]]
                filas.append({
                    "ligando": nombre, "modelo": clave, "modo": modo,
                    "afinidad": afinidad,
                    "contactos": len(residuos),
                    "atomos_del_ligando_en_contacto": en_contacto,
                    "fraccion_ligando_en_contacto": round(en_contacto / float(len(L)), 2),
                    "residuos": " ".join(str(r) for r in residuos),
                    "toca_trp334": toca_trp,
                    "dist_min_al_anillo_trp334": (round(dist_centro, 2)
                                                  if dist_centro is not None else None),
                    "dist_min_a_un_atomo_del_anillo": (
                        round(dist_anillo_atomo, 2)
                        if dist_anillo_atomo is not None else None),
                    "residuos_del_cr": " ".join(str(r) for r in del_cr),
                })
            log("  %-5s modelo %s: %d poses acopladas" % (nombre, clave, len(
                [f for f in filas if f["ligando"] == nombre and f["modelo"] == clave])))

    campos = ["ligando", "modelo", "modo", "afinidad", "contactos",
              "atomos_del_ligando_en_contacto", "fraccion_ligando_en_contacto",
              "residuos", "toca_trp334", "dist_min_al_anillo_trp334",
              "dist_min_a_un_atomo_del_anillo", "residuos_del_cr"]
    with open(os.path.join(SALIDA, "acoplamiento_cr.csv"), "w", encoding="utf-8",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(filas)

    lineas = []
    lineas.append("XL20 y XL23 acoplados en el CR del C-terminal de TDP-43")
    lineas.append("receptor %s (%s), 6 modelos de RMN, caja de %d A sobre el anillo"
                  " del Trp334" % (caja["pdb"], caja["residuos"], args.lado))
    lineas.append("exhaustividad %d, 9 modos, semilla %d | contactos a %.1f A"
                  % (args.exhaustividad, args.semilla, CORTE))
    lineas.append("")
    lineas.append("  %-5s %-6s %-7s %-9s %-9s %-8s %-8s %-8s %s"
                  % ("lig", "modelo", "afin", "contactos", "at.lig", "Trp334",
                     "d.aneo", "d.atomo", "residuos"))
    for fila in filas:
        lineas.append("  %-5s %-6s %-7.2f %-9d %-9d %-8s %-8s %-8s %s"
                      % (fila["ligando"], fila["modelo"], fila["afinidad"],
                         fila["contactos"], fila["atomos_del_ligando_en_contacto"],
                         fila["toca_trp334"], fila["dist_min_al_anillo_trp334"],
                         fila["dist_min_a_un_atomo_del_anillo"], fila["residuos"]))

    # --- resumen: la comparacion que se ha venido a hacer ---
    # Una pose cuenta como «apoyada en el Trp334» si el ligando tiene atomos a menos
    # de 4 A de un atomo del anillo indol. Todo lo demas es posarse en otro sitio de
    # la helice, que con una caja de 22 A sobre un peptido de 43 residuos es posible.
    lineas.append("")
    lineas.append("  RESUMEN (pose apoyada en el Trp334 = algun atomo a menos de %.1f A"
                  " de un atomo del anillo)" % CORTE)
    lineas.append("  %-5s %-7s %-9s %-11s %-9s %s"
                  % ("lig", "mejor", "mejor_Trp", "sobre_Trp334", "n_poses",
                     "residuos que toca la mejor pose"))
    resumen = {}
    for nombre in pdbqt_ligando:
        suyas = [f for f in filas if f["ligando"] == nombre and f["afinidad"] is not None]
        if not suyas:
            continue
        sobre = [f for f in suyas
                 if (f["dist_min_a_un_atomo_del_anillo"] or 99) < CORTE]
        mejor = min(suyas, key=lambda f: f["afinidad"])
        mejor_trp = min(sobre, key=lambda f: f["afinidad"]) if sobre else None
        pesados = len(atomos(pdbqt_ligando[nombre]))
        # Eficiencia de ligando: la afinidad repartida por atomo pesado. Es lo que
        # permite comparar un compuesto de 26 atomos con otro de 34 sin que el
        # segundo gane solo por ser mas grande, que es lo que hace Vina.
        resumen[nombre] = {"mejor": float(mejor["afinidad"]), "pesados": pesados,
                           "le": float(mejor["afinidad"]) / pesados, "sobre": len(sobre),
                           "n": len(suyas)}
        lineas.append("  %-5s %-7.2f %-9s %-11s %-9s %s"
                      % (nombre, mejor["afinidad"],
                         "%.2f" % mejor_trp["afinidad"] if mejor_trp else "-",
                         "%d de %d" % (len(sobre), len(suyas)),
                         len(suyas), mejor["residuos"]))
    lineas.append("")
    lineas.append("  POR MODELO: mejor afinidad de cada uno en cada modelo de la RMN")
    lineas.append("  (el CR desordenado no es una foto: si el orden cambiara de un modelo a otro,"
                  " no habria nada que decir)")
    lineas.append("  %-8s %-9s %-9s %-9s %s"
                  % ("modelo", "XL20", "XL23", "diferencia", "gana"))
    modelos = sorted({f["modelo"] for f in filas}, key=int)
    gana = {n: 0 for n in pdbqt_ligando}
    for mod in modelos:
        valores = {}
        for nombre in pdbqt_ligando:
            suyas = [f for f in filas if f["ligando"] == nombre
                     and f["modelo"] == mod and f["afinidad"] is not None]
            if suyas:
                valores[nombre] = min(float(f["afinidad"]) for f in suyas)
        if len(valores) < 2:
            continue
        mejor_modelo = min(valores, key=lambda n: valores[n])
        gana[mejor_modelo] += 1
        lineas.append("  %-8s %-9.2f %-9.2f %-9.2f %s"
                      % (mod, valores["XL20"], valores["XL23"],
                         valores["XL23"] - valores["XL20"], mejor_modelo))
    lineas.append("")
    lineas.append("  eficiencia de ligando (afinidad por atomo pesado) de la mejor pose:")
    for nombre in sorted(resumen, key=lambda n: -resumen[n]["mejor"]):
        r = resumen[nombre]
        lineas.append("    %-5s %d atomos pesados | mejor %.2f kcal/mol | %.3f kcal/mol por atomo"
                      % (nombre, r["pesados"], r["mejor"], r["le"]))
    lineas.append("    ganados por modelo: %s"
                  % ", ".join("%s %d de %d" % (n, gana[n], len(modelos))
                              for n in sorted(gana)))
    texto = "\n".join(lineas)
    with open(os.path.join(SALIDA, "acoplamiento_cr.txt"), "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    log("")
    log(texto)
    return 0


if __name__ == "__main__":
    sys.exit(main())
