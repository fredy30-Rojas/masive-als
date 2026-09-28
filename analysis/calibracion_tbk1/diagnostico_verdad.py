"""Construye varias verdades de referencia y mide Vina y MM-GBSA contra cada una.

QUE PREGUNTA CONTESTA
---------------------
El AUROC de 0,602 de Vina contra la lista actual no dice si el docking es malo:
puede que la lista sea mala. Este script rehace la lista con criterios mas
estrictos y mide contra cada una, para separar las dos cosas.

LOS TRES CRITERIOS
  1. EL DEL BANCO: por compuesto, SU MEJOR pChEMBL. 2172 activos, 141 inactivos.
  2. LA MEDIANA: la mediana de sus medidas, que es lo que representa de verdad
     al compuesto. 1552 activos, 141 inactivos. De aqui sale que un 29% de los
     "activos" dejan de serlo.
  3. EL ESTRICTO: solo medidas de UNION (IC50/Ki/Kd/EC50, que son las unicas
     que tienen pChEMBL), sin limites, y con al menos DOS medidas por clase.
     346 activos, 17 inactivos.

Ademas se separa por cuantas medidas tiene el compuesto, porque un unico dato
publicado no es una verdad de laboratorio.

LO QUE HAY QUE MIRAR
Si el AUROC de Vina sube al usar la lista estricta, el problema era la lista.
Si se queda igual o baja, el problema es el metodo, y la lista buena no lo
salva. Las dos cosas son posibles y tienen Crack y consecuencias distintas.
"""
import argparse
import collections
import csv
import json
import math
import os
import statistics
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
ACTIVIDADES = os.path.join(BASE, "actividades_tbk1.json")
CSV_BANCO = os.path.join(BASE, "compuestos_tbk1.csv")

UNION = ("IC50", "Ki", "Kd", "EC50")
LIMITE = ("<", ">", ">=", "<=", "~")


def cargar_actividades():
    acts = json.load(open(ACTIVIDADES, encoding="utf-8"))
    por = collections.defaultdict(list)
    for a in acts:
        try:
            p = float(a.get("pchembl_value"))
        except (TypeError, ValueError):
            continue
        rel = a.get("standard_relation") or "="
        por[a.get("molecule_chembl_id")].append(
            (a.get("standard_type"), p, rel, a.get("target_organism")))
    return por


def afinidad_vina(pdbqt):
    dentro = False
    try:
        fh = open(pdbqt, errors="replace")
    except OSError:
        return None
    with fh:
        for linea in fh:
            if linea.startswith("MODEL"):
                dentro = True
                continue
            if linea.startswith("ENDMDL") and dentro:
                break
            if dentro and linea.startswith("REMARK VINA RESULT"):
                try:
                    return float(linea.split()[3])
                except (IndexError, ValueError):
                    return None
    return None


def auroc(activo, inactivo):
    if not activo or not inactivo:
        return None, None
    n = len(activo) * len(inactivo)
    suma = sum(1.0 if a < i else (0.5 if a == i else 0.0)
               for a in activo for i in inactivo)
    a = suma / n
    na, ni = len(activo), len(inactivo)
    q1 = a / (2.0 - a)
    q2 = 2.0 * a * a / (1.0 + a)
    var = (a * (1 - a) + (na - 1) * (q1 - a * a) + (ni - 1) * (q2 - a * a)) / (na * ni)
    return a, math.sqrt(var) if var > 0 else None


def verdad_por_criterio(por, tipo="IC50"):
    """(id -> True/False) para cada uno de los tres criterios."""
    banco, mediana, estricto = {}, {}, {}
    for cid, vals in por.items():
        uniones = [p for t, p, rel, _ in vals if t in UNION]
        sin_limite = [p for t, p, rel, _ in vals if t in UNION and rel == "="]
        banco[cid] = max(p for _, p, _, _ in vals) >= 6.0
        if uniones:
            mediana[cid] = statistics.median(uniones) >= 6.0
        if len(sin_limite) >= 2:
            n_act = sum(1 for p in sin_limite if p >= 6.0)
            n_ina = sum(1 for p in sin_limite if p < 6.0)
            if n_act >= 2:
                estricto[cid] = True
            elif n_ina >= 2:
                estricto[cid] = False
    return {"banco (mejor valor)": banco,
            "mediana de las medidas": mediana,
            "estricto (2 medidas, sin limites)": estricto}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--criba", default=os.path.join(BASE, "criba_250.csv"))
    ap.add_argument("--dir-poses", default=os.path.join(BASE, "out"))
    args = ap.parse_args()

    por = cargar_actividades()
    dg = {}
    if os.path.exists(args.criba):
        for r in csv.DictReader(open(args.criba, encoding="utf-8")):
            if r["dG_union"]:
                dg[r["ligand"]] = float(r["dG_union"])
        print("dG del MM-GBSA: %d ligandos" % len(dg))

    # afinidades de Vina para todos los compuestos con actividad
    vina = {}
    for cid in por:
        p = os.path.join(args.dir_poses, "ACT_%s_out.pdbqt" % cid)
        v = afinidad_vina(p)
        if v is not None:
            vina[cid] = v
    print("con afinidad de Vina: %d de %d" % (len(vina), len(por)))

    n_med = collections.Counter(min(5, len(v)) for v in por.values())
    print("medidas por compuesto: %s"
          % {("%d+" % k if k == 5 else str(k)): n for k, n in sorted(n_med.items())})

    print("\n%-34s %7s %7s %9s %9s"
          % ("criterio de verdad", "activos", "inact.", "Vina", "MM-GBSA"))
    for nombre, tabla in verdad_por_criterio(por).items():
        comunes = [c for c in tabla if c in vina]
        a = [c for c in comunes if tabla[c]]
        i = [c for c in comunes if not tabla[c]]
        if len(i) < 5:
            print("%-34s %7d %7d   (pocos inactivos: no se puede medir)"
                  % (nombre, len(a), len(i)))
            continue
        av, ev = auroc([vina[c] for c in a], [vina[c] for c in i])
        con = [c for c in comunes if c in dg]
        if len(con) > 10 and sum(1 for c in con if not tabla[c]) >= 5:
            am, em = auroc([dg[c] for c in con if tabla[c]],
                           [dg[c] for c in con if not tabla[c]])
            txt = "%6.3f ±%.3f" % (am, em if em else 0)
        else:
            txt = "     --  "
        print("%-34s %7d %7d %9.3f %9s"
              % (nombre, len(a), len(i), av, txt))
        print("%-34s %7s %7s   error de Vina ±%.3f" % ("", "", "", ev or 0))

    # y ademas: como rinde Vina segun cuantas medidas tenga el compuesto
    print("\nVina segun cuantas medidas tiene el compuesto (criterio banco):")
    banco = verdad_por_criterio(por)["banco (mejor valor)"]
    for etiq, condicion in (("1 sola medida", lambda n: n == 1),
                           ("2 o mas", lambda n: n >= 2)):
        a, i = [], []
        for cid, vals in por.items():
            if cid not in vina or cid not in banco:
                continue
            if not condicion(len(vals)):
                continue
            (a if banco[cid] else i).append(vina[cid])
        if len(i) >= 5 and a:
            av, ev = auroc(a, i)
            print("  %-14s %5d activos %4d inactivos  AUROC %.3f ±%.3f"
                  % (etiq, len(a), len(i), av, ev or 0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
