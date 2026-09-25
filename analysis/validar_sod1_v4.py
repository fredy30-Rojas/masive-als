#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validacion de SOD1 v4: por fin hay QUIMIAS INDEPENDIENTES en el bolsillo Trp32.

POR QUE EXISTE ESTA VERSION
---------------------------
El criterio que el proyecto se fijo el 20 de septiembre exige positivos de AL
MENOS DOS QUIMIAS INDEPENDIENTES para poder reportar un ranking. Hasta ahora
solo habia dos candidatos independientes (LCS-1, PRG-A01) y ninguno con
estructura. La via para conseguirlos no era buscar mas articulos al azar, sino
recorrer las 156 estructuras de SOD1 humana depositadas y medir que moleculas
pequenas estan a menos de 6 A del anillo de Trp32
(`ligandos_cristal_trp32.py`, resultado en `trp32_cristal/`).

Salieron SIETE quimias nuevas, todas con evidencia estructural directa:
quinazolina (4A7G), diazepano-quinazolina (4A7Q), trifluorometil-quinazolina
(2WZ6), anilina (2WZ0), paliperidona (8GSQ, farmaco aprobado, con MST),
fenantridinona Lig9 (6A9O) y aminoalcohol naftalenico (5YTO). Mas las tres
catecolaminas y la 5-fluorouridina que ya estaban.

QUE HACE
--------
1. Acopla TODOS los controles de `controles_sod1_v4.csv` en la MISMA caja Trp32
   del cribado, con el MISMO receptor y el mismo pipeline que los señuelos.
2. Reutiliza los fondos ya acoplados en la validacion v3 (emparejado nuevo y
   duro) y los 199 señuelos antiguos, para que los numeros sean comparables.
3. Calcula AUC bruto y AUC sin tamaño, EF1% y EF5%, y el puesto de cada control.
4. Comprueba el criterio: cuantas QUIMIAS distintas (Murcko) hay por delante del
   fondo, y si el resultado se puede reportar.

Uso:
    python validar_sod1_v4.py --preparar
    python validar_sod1_v4.py --acoplar [--workers 12]
    python validar_sod1_v4.py --analizar
    python validar_sod1_v4.py            # las tres cosas
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold

import validar_senuelos as VS
import validar_sod1_v3 as V3

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(BASE, "validacion_SOD1_v4")
LIGDIR = os.path.join(WORK, "ligands")
OUTDIR = os.path.join(WORK, "out")
CONTROLES = os.path.join(BASE, "controles_sod1_v4.csv")
LIBRERIA = os.path.abspath(os.path.join(BASE, "..", "compounds", "decoys_library.smi"))
RECEPTOR = os.path.abspath(os.path.join(BASE, "..", "gpu_dock", "SOD1.pdbqt"))
V3_OUT = os.path.join(BASE, "validacion_SOD1_v3", "out")
V3_LIG = os.path.join(BASE, "validacion_SOD1_v3", "ligands")
ANTIGUO = os.path.join(BASE, "_validacion_SOD1", "validacion_SOD1_trp32.csv")
CENTRO = (46.5, 80.0, 73.3)
TAMANO = 22
EXHAUSTIVIDAD = 8

# nombre en el csv -> nombre del fichero pdbqt
FICHERO = {
    "isoproterenol": "ACT_isoproterenol",
    "adrenalina": "ACT_adrenalina",
    "dopamina": "ACT_dopamina",
    "5-fluorouridina": "ACT_5-fluorouridina",
    "LCS-1": "ACT_LCS-1",
    "PRG-A01": "ACT_PRG-A01",
    "CHEMBL2165613": "ACT_CHEMBL2165613",
}

# Familias: lo que se pregunta con cada una es distinto.
FAMILIAS = {
    # una sola quimia (3 compuestos, mismo esqueleto): valida el SITIO
    "catecolaminas (1 quimia)": ["isoproterenol", "adrenalina", "dopamina"],
    # las quimias nuevas del PDB, cada una de una familia distinta
    "cristalizadas, quimias nuevas": [
        "5-fluorouridina", "quinazolina_12I", "diazepanoquinazolina_4MQ",
        "cf3quinazolina_ZO0", "anilina_ZZT", "paliperidona_K4I",
        "Lig9_6B3", "naftalenoaminoalcohol_946",
    ],
    "Lig9 sola (manjula)": ["Lig9_6B3"],
    "paliperidona sola (farmaco)": ["paliperidona_K4I"],
    "sin estructura (actividad celular)": ["LCS-1", "PRG-A01"],
    "serie pirazolona": ["CHEMBL2165613"],
    "cristalizadas + sin estructura": [
        "5-fluorouridina", "quinazolina_12I", "diazepanoquinazolina_4MQ",
        "cf3quinazolina_ZO0", "anilina_ZZT", "paliperidona_K4I", "Lig9_6B3",
        "naftalenoaminoalcohol_946", "LCS-1", "PRG-A01",
    ],
    "todos los controles": [
        "isoproterenol", "adrenalina", "dopamina", "5-fluorouridina",
        "quinazolina_12I", "diazepanoquinazolina_4MQ", "cf3quinazolina_ZO0",
        "anilina_ZZT", "paliperidona_K4I", "Lig9_6B3",
        "naftalenoaminoalcohol_946", "LCS-1", "PRG-A01", "CHEMBL2165613",
    ],
}


def log(m):
    print(m, flush=True)


def leer_controles():
    with open(CONTROLES, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def mapa_ficheros():
    """Rellena FICHERO para TODO control (nombre del csv -> fichero pdbqt).

    Las tres fases (preparar, acoplar, analizar) pueden correr en procesos
    distintos, asi que el mapa no puede depender de haber preparado antes.
    """
    for r in leer_controles():
        FICHERO.setdefault(r["ligand"], "ACT_" + r["ligand"])
    return FICHERO


def preparar():
    os.makedirs(LIGDIR, exist_ok=True)
    con = mapa_ficheros()
    for r in leer_controles():
        nombre = r["ligand"]
        fichero = con[nombre]
        destino = os.path.join(LIGDIR, fichero + ".pdbqt")
        if os.path.exists(destino) and os.path.getsize(destino) > 100:
            continue
        de_v3 = os.path.join(V3_LIG, fichero + ".pdbqt")
        if os.path.exists(de_v3):  # mismos parametros, mismo fichero
            shutil.copy(de_v3, destino)
            log("  %-30s copiado de la v3" % nombre)
            continue
        if VS.convertir_pdbqt(fichero, r["smiles"], LIGDIR):
            log("  %-30s preparado (%s)" % (nombre, r["quimia"]))
        else:
            log("  %-30s FALLO al preparar" % nombre)


def acoplar(workers=12):
    mapa_ficheros()
    os.makedirs(OUTDIR, exist_ok=True)
    tareas = []
    for p in sorted(glob.glob(os.path.join(LIGDIR, "*.pdbqt"))):
        nombre = os.path.basename(p).replace(".pdbqt", "")
        out = os.path.join(OUTDIR, nombre + "_out.pdbqt")
        if os.path.exists(out) and os.path.getsize(out) > 100:
            continue
        tareas.append((VS.VINA_CPU, RECEPTOR, p, out, CENTRO[0], CENTRO[1],
                       CENTRO[2], TAMANO, EXHAUSTIVIDAD))
    log("pendientes de acoplar: %d (caja %.1f,%.1f,%.1f; %d A; ex %d)"
        % (len(tareas), *CENTRO, TAMANO, EXHAUSTIVIDAD))
    if tareas:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for i, _ in enumerate(ex.map(V3._dock_one, tareas), 1):
                log("   ... %d/%d" % (i, len(tareas)))
    log("acoplados en total: %d" % len(glob.glob(os.path.join(OUTDIR, "*_out.pdbqt"))))


def afinidades():
    res = {}
    for p in glob.glob(os.path.join(OUTDIR, "*_out.pdbqt")):
        n = os.path.basename(p).replace("_out.pdbqt", "")
        a = VS.parse_affinity(p)
        if a is not None:
            res[n] = a
    return res


def tamanos_control():
    n = {}
    for r in leer_controles():
        m = Chem.MolFromSmiles(r["smiles"])
        if m is not None:
            n[r["ligand"]] = m.GetNumHeavyAtoms()
    return n


def pesados_pdbqt(ruta):
    """Atomos pesados contando el propio PDBQT (vale para señuelos sin SMILES).

    Contar cabria esperar que la libreria los tuviera todos, pero los señuelos
    nuevos se generaron con nombres que no estan en el .smi: leer el fichero que
    se acoplo es la medida directa y no depende de emparejar nombres.
    """
    n = 0
    try:
        with open(ruta, encoding="utf-8", errors="replace") as f:
            for l in f:
                if not l.startswith(("ATOM", "HETATM")):
                    continue
                # En PDBQT el ultimo campo es el TIPO de atomo de Vina (C, A, N,
                # NA, OA, SA, H, HD...). La columna 76:78 no es el elemento
                # quimico, asi que se lee el tipo, que es lo que distingue H.
                tipo = l.rsplit(None, 1)[-1].upper()
                if tipo not in ("H", "HD", "HS", "D", "DD"):
                    n += 1
    except OSError:
        return None
    return n or None


def tamanos_de_ficheros():
    """nombre de fichero (sin _out) -> atomos pesados, leidos del PDBQT."""
    n = {}
    for d in (LIGDIR, V3_LIG, os.path.join(BASE, "_validacion_SOD1", "ligands")):
        for p in glob.glob(os.path.join(d, "*.pdbqt")):
            k = os.path.basename(p).replace(".pdbqt", "")
            v = pesados_pdbqt(p)
            if v and k not in n:
                n[k] = v
    return n


def quimias_murcko():
    q = {}
    for r in leer_controles():
        m = Chem.MolFromSmiles(r["smiles"])
        q[r["ligand"]] = MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m else ""
    return q


def analizar():
    mapa_ficheros()
    aff = afinidades()
    antiguo = {}
    with open(ANTIGUO, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["affinity"]:
                antiguo[r["ligand"]] = float(r["affinity"])

    # fondos: los mismos que en la v3 (mismos ficheros, misma caja)
    fondo_emp, fondo_duro = {}, {}
    for p in glob.glob(os.path.join(V3_OUT, "*_out.pdbqt")):
        n = os.path.basename(p).replace("_out.pdbqt", "")
        a = VS.parse_affinity(p)
        if a is None:
            continue
        if n.startswith("DECM_"):
            fondo_emp[n] = a
        elif n.startswith("DECH_"):
            fondo_duro[n] = a
    fondo_antiguo = {k: v for k, v in antiguo.items() if k.startswith("DEC_")}
    fondo_todo = {**fondo_emp, **fondo_duro, **fondo_antiguo}
    log("fondos: emparejado nuevo %d | duro %d | emparejado antiguo %d"
        % (len(fondo_emp), len(fondo_duro), len(fondo_antiguo)))

    # tamaño: de los señuelos y de los controles
    tam = {}
    for nombre, smi in VS.leer_libreria(LIBRERIA):
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            tam[nombre] = m.GetNumHeavyAtoms()
    tam.update(tamanos_control())

    def tam_de(k):
        base = k.split("_", 1)[1] if k.startswith(("DEC_", "DECM_", "DECH_", "ACT_")) else k
        return tam.get(base) or tam.get(k)

    # el tamaño se lee del propio fichero acoplado: no depende de emparejar nombres
    tam_fichero = tamanos_de_ficheros()

    def tam_efectivo(k):
        return tam_fichero.get(k) or tam_de(k)

    # OJO: los fondos NO estan en `aff` (que solo trae los controles de esta
    # carpeta); hay que meterlos en `scores` o el ajuste por tamaño se queda sin
    # señuelos y el sesgo sale mal.
    scores = {**antiguo, **aff, **fondo_emp, **fondo_duro}
    tam_lig = {k: tam_efectivo(k) for k in scores if tam_efectivo(k)}
    sin_tam = [k for k in fondo_todo if k not in tam_lig]
    resid, a_aj, b_aj = V3.normalizar_por_tamano(scores, tam_lig, list(fondo_todo))
    if b_aj is not None:
        log("sesgo de tamaño (ajustado SOLO con el fondo, %d de %d con tamaño): "
            "afinidad = %.2f %+.3f * pesados"
            % (len(fondo_todo) - len(sin_tam), len(fondo_todo), a_aj, b_aj))
        if sin_tam:
            log("   sin tamaño conocido (fuera del ajuste): %d, p.ej. %s"
                % (len(sin_tam), ", ".join(sin_tam[:5])))

    def val(nombre, d):
        return d.get(FICHERO[nombre])

    filas, detalle = [], []
    for nombre in FAMILIAS["todos los controles"]:
        f = FICHERO[nombre]
        v = aff.get(f, antiguo.get(f))
        detalle.append({"ligando": nombre, "fichero": f, "afinidad": v})

    quimias = quimias_murcko()
    for fam, miembros in FAMILIAS.items():
        act = [v for v in (val(m, aff) for m in miembros) if v is not None]
        act_r = [resid[FICHERO[m]] for m in miembros if FICHERO[m] in resid]
        for etiq, fondo in (("emparejado nuevo", fondo_emp),
                            ("emparejado antiguo", fondo_antiguo),
                            ("duro", fondo_duro),
                            ("los tres juntos", fondo_todo)):
            f_r = {k: resid[k] for k in fondo if k in resid}
            auc = V3.roc_auc(act, list(fondo.values()))
            auc_r = V3.roc_auc(act_r, list(f_r.values())) if act_r and f_r else None
            filas.append({
                "familia": fam, "n_activos": len(act), "fondo": etiq,
                "n_fondo": len(fondo),
                "media_activos": round(float(np.mean(act)), 2) if act else None,
                "media_fondo": round(float(np.mean(list(fondo.values()))), 2) if fondo else None,
                "AUC": round(auc, 3) if auc is not None else None,
                "AUC_sin_tamano": round(auc_r, 3) if auc_r is not None else None,
                "EF1%": round(V3.ef(act, list(fondo.values()), 1.0), 2) if act and fondo else None,
                "EF5%": round(V3.ef(act, list(fondo.values()), 5.0), 2) if act and fondo else None,
            })

    for d in detalle:
        if d["afinidad"] is None:
            continue
        d["quimia_murcko"] = quimias.get(d["ligando"], "")
        d["pesados"] = tam_fichero.get(d["fichero"]) or tam.get(d["ligando"], "")
        d["puesto_frente_a_los_tres_fondos"] = V3.puesto(
            [d["afinidad"]], list(fondo_todo.values()), d["afinidad"])
        d["n_fondo"] = len(fondo_todo)
        d["cuantil"] = round(100.0 * d["puesto_frente_a_los_tres_fondos"]
                             / (len(fondo_todo) + 1), 1)
        d["residuo_sin_tamano"] = round(resid.get(d["fichero"], float("nan")), 3)

    with open(os.path.join(WORK, "analisis_sod1_v4.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["familia", "n_activos", "fondo", "n_fondo",
                                          "media_activos", "media_fondo", "AUC",
                                          "AUC_sin_tamano", "EF1%", "EF5%"],
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(filas)
    with open(os.path.join(WORK, "detalle_controles_sod1_v4.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ligando", "fichero", "afinidad", "pesados",
                                          "quimia_murcko", "residuo_sin_tamano",
                                          "puesto_frente_a_los_tres_fondos", "n_fondo",
                                          "cuantil"], extrasaction="ignore")
        w.writeheader()
        w.writerows(detalle)

    log("")
    log("%-34s %-19s %-5s %-6s %-9s %-14s %s"
        % ("familia", "fondo", "n_act", "AUC", "sin tamaño", "EF1%", "EF5%"))
    for r in filas:
        log("%-34s %-19s %-5s %-6s %-9s %-14s %s"
            % (r["familia"], r["fondo"], r["n_activos"], r["AUC"],
               r["AUC_sin_tamano"], r["EF1%"], r["EF5%"]))
    log("")
    log("controles frente a los tres fondos juntos (N=%d):" % len(fondo_todo))
    for d in sorted(detalle, key=lambda x: (x["afinidad"] is None, x["afinidad"])):
        log("   %-30s %-9s puesto %-4s de %d   mejor que el %5.1f%%   %s"
            % (d["ligando"],
               "%.2f" % d["afinidad"] if d["afinidad"] is not None else "n/d",
               d.get("puesto_frente_a_los_tres_fondos", "n/d"), len(fondo_todo),
               100.0 - (d.get("cuantil") or 100.0), d.get("quimia_murcko", "")[:38]))
    log("")
    log("quimias (Murcko) distintas entre los controles: %d"
        % len({q for q in quimias.values() if q}))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preparar", action="store_true")
    ap.add_argument("--acoplar", action="store_true")
    ap.add_argument("--analizar", action="store_true")
    ap.add_argument("--workers", type=int, default=12)
    args = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    todo = not any([args.preparar, args.acoplar, args.analizar])
    if args.preparar or todo:
        log("=== preparando controles ===")
        preparar()
    if args.acoplar or todo:
        log("=== acoplando ===")
        acoplar(workers=args.workers)
    if args.analizar or todo:
        log("=== analizando ===")
        analizar()
    return 0


if __name__ == "__main__":
    sys.exit(main())
