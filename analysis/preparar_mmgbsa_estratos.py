#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""preparar_mmgbsa_estratos.py — la muestra para probar MM-GBSA por estratos de tamaño.

QUE SE QUIERE SABER
-------------------
En la validacion de SOD1 contra el receptor limpio se midio que la ordenacion por
energia de Vina es, en la practica, ordenacion por tamaño:

  * dentro del fondo, la afinidad correlaciona -0,720 con los atomos pesados;
  * la cabeza de la lista (25 mejores) tiene 27,1 atomos de media frente a 19,7
    del fondo;
  * el AUC crudo es 0,439 (por debajo del azar) mientras que **emparejando cada
    positivo contra señuelos de SU MISMO TAMAÑO, 5 de 8 quimiotipos ganan**.

La pregunta es si el MM-GBSA, que no es una funcion de superficie de contacto, quita
ese sesgo o lo tiene igual. Y para responderla no hace falta pasar los 493 ligandos:
hace falta una muestra en la que **cada positivo tenga señuelos de su mismo tamaño**
y en la que haya señuelos repartidos por todo el rango de tamaños.

COMO SE ELIGE LA MUESTRA
------------------------
1. Los 11 positivos de union medida de `verdad_de_referencia.csv`.
2. Para cada positivo, hasta `--por-positivo` señuelos con diferencia de atomos
   pesados <= 2. Es la comparacion justa y es la que decide.
3. Señuelos repartidos por cuantiles de tamaño (hasta `--repartidos`), para poder
   ajustar la recta afinidad-tamaño sobre el fondo y para que haya fondo en los
   tamaños extremos.
4. Nada se solapa: un señuelo entra una sola vez.

Salida (para subir a Oracle):
  mmgbsa_estratos/candidatos.csv   ligand,target,affinity,smiles,pose_pdbqt
  mmgbsa_estratos/poses/*.pdbqt    las poses de la muestra
  mmgbsa_estratos/receptor.pdb     el receptor LIMPIO (copia de SOD1_limpio.pdb)

Uso:
  python preparar_mmgbsa_estratos.py
"""
import csv
import os
import shutil
import sys

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
POSES = os.path.join(BASE, "validacion_SOD1_v5", "out")
LIGS = os.path.join(BASE, "validacion_SOD1_v5", "ligands")
RECEPTOR = os.path.join(RAIZ, "gpu_dock", "SOD1_limpio.pdb")
TRABAJO = os.path.join(BASE, "mmgbsa_estratos")

# El runner de Oracle resuelve `pose_pdbqt` relativo a SU directorio de trabajo,
# asi que en el CSV va la ruta absoluta del servidor (si no, falla con
# "No such file or directory: 'poses/...'" en 0,7 s sin haberse quejado de nada).
RUTA_ORACLE = "/home/ubuntu/mmgbsa/estratos"

SEMILLA = 20260921


def pesados_de(smi):
    """Atomos pesados. `GetNumHeavyAtoms`: los [2H] explicitos no cuentan."""
    m = Chem.MolFromSmiles(smi) if smi else None
    return None if m is None else m.GetNumHeavyAtoms()


def smi_de_fichero(ruta):
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK SMILES "):
            return l[len("REMARK SMILES "):].strip()
        if l.startswith(("ATOM", "HETATM")):
            break
    return None


def afinidad(ruta):
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            try:
                return float(l.split()[3])
            except (IndexError, ValueError):
                return None
    return None


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--por-positivo", type=int, default=6)
    ap.add_argument("--repartidos", type=int, default=45)
    ap.add_argument("--max-ligandos", type=int, default=90)
    args = ap.parse_args()

    posdir = os.path.join(TRABAJO, "poses")
    os.makedirs(posdir, exist_ok=True)

    # ------------------------------------------------------------------ positivos
    positivos = [r for r in csv.DictReader(
        open(os.path.join(BASE, "verdad_de_referencia.csv"), encoding="utf-8"))
        if r["target"] == "SOD1" and r["apto"] == "si"]

    # ---------------------------------------------------------------- el fondo
    fondo = []
    for p in sorted(os.listdir(POSES)):
        if not p.endswith("_out.pdbqt") or p.startswith("ACT_"):
            continue
        nombre = p.replace("_out.pdbqt", "")
        lig = os.path.join(LIGS, nombre + ".pdbqt")
        if not os.path.exists(lig):
            continue
        n = pesados_de(smi_de_fichero(lig))
        aff = afinidad(os.path.join(POSES, p))
        if n is None or aff is None:
            continue
        fondo.append({"nombre": nombre, "n": n, "aff": aff, "smi": smi_de_fichero(lig)})
    print("fondo disponible: %d señuelos" % len(fondo))

    rng = np.random.default_rng(SEMILLA)
    elegidos = {}
    razon = {}

    # 1) emparejados por tamaño, seis por positivo
    for r in positivos:
        n = pesados_de(r["smiles"])
        par = [f for f in fondo if abs(f["n"] - n) <= 2 and f["nombre"] not in elegidos]
        par.sort(key=lambda f: abs(f["n"] - n))
        # de los mas cercanos en tamaño, se coge una muestra al azar para no
        # quedarse siempre con los mismos
        cabeza = par[:max(args.por_positivo * 3, args.por_positivo)]
        rng.shuffle(cabeza)
        for f in cabeza[:args.por_positivo]:
            elegidos[f["nombre"]] = f
            razon[f["nombre"]] = "emparejado con %s (%d atomos)" % (r["ligand"], n)

    # 2) repartidos por cuantiles de tamaño, para poder ajustar la recta del fondo
    libres = [f for f in fondo if f["nombre"] not in elegidos]
    if libres and args.repartidos > 0:
        tam = np.array([f["n"] for f in libres])
        cortes = np.quantile(tam, np.linspace(0, 1, args.repartidos + 1))
        for i in range(args.repartidos):
            lo, hi = cortes[i], cortes[i + 1]
            grupo = [f for f in libres
                     if lo <= f["n"] <= hi and f["nombre"] not in elegidos]
            if not grupo:
                continue
            f = grupo[int(rng.integers(len(grupo)))]
            elegidos[f["nombre"]] = f
            razon[f["nombre"]] = "repartido por tamaño (cuantil %d/%d)" % (
                i + 1, args.repartidos)

    # ------------------------------------------------------------------- seleccion
    filas = []
    for r in positivos:
        nombre = "ACT_" + r["ligand"]
        src = os.path.join(POSES, nombre + "_out.pdbqt")
        if not os.path.exists(src):
            print("   sin pose: %s" % r["ligand"])
            continue
        filas.append({"ligand": nombre, "target": "SOD1",
                      "affinity": "%.3f" % afinidad(src), "smiles": r["smiles"],
                      "pose_pdbqt": "%s/poses/%s_out.pdbqt" % (RUTA_ORACLE, nombre),
                      "papel": "positivo", "quimia": r["quimia"], "_src": src})

    for nombre, f in sorted(elegidos.items(), key=lambda kv: kv[1]["n"]):
        src = os.path.join(POSES, nombre + "_out.pdbqt")
        if not os.path.exists(src):
            continue
        filas.append({"ligand": nombre, "target": "SOD1",
                      "affinity": "%.3f" % f["aff"], "smiles": f["smi"],
                      "pose_pdbqt": "%s/poses/%s_out.pdbqt" % (RUTA_ORACLE, nombre),
                      "papel": "fondo", "quimia": razon.get(nombre, ""), "_src": src})

    if len(filas) > args.max_ligandos:
        # se recorta quitando fondo repartido, nunca positivos ni emparejados
        sobra = len(filas) - args.max_ligandos
        rep = [f for f in filas if f["papel"] == "fondo" and "repartido" in f["quimia"]]
        quitar = {f["ligand"] for f in rep[:sobra]}
        filas = [f for f in filas if f["ligand"] not in quitar]
        print("recortados %d señuelos repartidos para no pasar de %d"
              % (len(quitar), args.max_ligandos))

    # las poses se copian DESPUES de recortar, para que la carpeta y el CSV
    # tengan exactamente los mismos ligandos y nadie se pregunte por los sueltos
    for f in os.listdir(posdir):
        os.remove(os.path.join(posdir, f))
    for r in filas:
        shutil.copy(r["_src"], os.path.join(posdir, r["ligand"] + "_out.pdbqt"))
        del r["_src"]

    campos = ["ligand", "target", "affinity", "smiles", "pose_pdbqt", "papel",
              "quimia"]
    with open(os.path.join(TRABAJO, "candidatos.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for r in filas:
            w.writerow(r)

    shutil.copy(RECEPTOR, os.path.join(TRABAJO, "receptor.pdb"))

    n_pos = sum(1 for r in filas if r["papel"] == "positivo")
    n_fon = len(filas) - n_pos
    tam_pos = np.mean([pesados_de(r["smiles"]) for r in filas if r["papel"] == "positivo"])
    tam_fon = np.mean([pesados_de(r["smiles"]) for r in filas if r["papel"] == "fondo"])
    emp = sum(1 for r in filas if "emparejado" in r["quimia"])

    print("")
    print("muestra: %d ligandos (%d positivos + %d fondo)" % (len(filas), n_pos, n_fon))
    print("   de los del fondo, %d estan emparejados por tamaño con un positivo" % emp)
    print("   atomos pesados: positivos %.1f | fondo %.1f" % (tam_pos, tam_fon))
    print("   cada positivo tiene %s señuelos dentro de +-2 atomos"
          % ", ".join(str(sum(1 for f in filas if f["papel"] == "fondo"
                             and abs(pesados_de(f["smiles"])
                                     - pesados_de(r["smiles"])) <= 2))
                      for r in filas if r["papel"] == "positivo"))
    print("")
    print("   escrito: %s" % os.path.join(TRABAJO, "candidatos.csv"))
    print("   poses:   %s" % posdir)
    print("   receptor limpio: %s" % os.path.join(TRABAJO, "receptor.pdb"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
