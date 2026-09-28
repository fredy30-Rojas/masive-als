"""Compara el MM-GBSA con lo que ya da Vina, sobre los MISMOS ligandos.

POR QUE HAY QUE COMPARAR Y NO SOLO MIRAR EL MM-GBSA (28 sep 2026): el
rescoring solo tiene sentido si MEJORA el docking. Vina ya da una afinidad
por ligando, gratis, en 3 segundos. Si el MM-GBSA ordena igual de bien que
Vina, no aporta nada y su coste (36 s por ligando) es puro gasto. Si lo
ordena PEOR, hay que decirlo y buscar otra via.

Por eso el AUROC se calcula para los tres:
  - Vina sola (la referencia)
  - MM-GBSA solo
  - la suma de los dos, con el peso que mejor clasifique
Y se mide el margen de error, porque con veinte ligandos el AUROC de 0,567
no significa nada y hay que decirlo.

La afinidad de Vina se lee del propio pdbqt de salida, en la linea
`REMARK VINA RESULT`, que es la primera del modelo 1. No hace falta volver a
correr el docking.
"""
import argparse
import csv
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))


def afinidad_vina(pdbqt):
    """La afinidad del primer modelo del pdbqt, o None."""
    dentro = False
    for linea in open(pdbqt, errors="replace"):
        if linea.startswith("MODEL"):
            dentro = True
            continue
        if linea.startswith("ENDMDL") and dentro:
            break
        if dentro and linea.startswith("REMARK VINA RESULT"):
            partes = linea.split()
            try:
                return float(partes[3])
            except (IndexError, ValueError):
                return None
    return None


def rango(xs):
    orden = sorted(range(len(xs)), key=lambda i: xs[i])
    r = [0.0] * len(xs)
    i = 0
    while i < len(orden):
        j = i
        while j + 1 < len(orden) and xs[orden[j + 1]] == xs[orden[i]]:
            j += 1
        medio = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            r[orden[k]] = medio
        i = j + 1
    return r


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxy = sum((a - mx) * (b - my) for a, b in zip(xs, ys))
    sxx = sum((a - mx) ** 2 for a in xs)
    syy = sum((b - my) ** 2 for b in ys)
    if sxx <= 0 or syy <= 0:
        return None
    return sxy / (math.sqrt(sxx) * math.sqrt(syy))


def auroc(activo, inactivo):
    """Fraccion de pares (activo mejor unido que inactivo). Mas bajo = mejor."""
    if not activo or not inactivo:
        return None, None
    n = len(activo) * len(inactivo)
    suma = 0.0
    for a in activo:
        for i in inactivo:
            suma += 1.0 if a < i else (0.5 if a == i else 0.0)
    a = suma / n
    # error estandar de Hanley-McNeil, que con una sola clase compuesta no
    # vale; con dos clases y n skewed, se usa la aproximacion normal
    na, ni = len(activo), len(inactivo)
    q1 = a / (2.0 - a)
    q2 = 2.0 * a * a / (1.0 + a)
    var = (a * (1 - a) + (na - 1) * (q1 - a * a) + (ni - 1) * (q2 - a * a)) / (na * ni)
    return a, math.sqrt(var) if var > 0 else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--muestra", default=os.path.join(BASE, "muestra_250.csv"))
    ap.add_argument("--dir-poses", default=os.path.join(BASE, "out"))
    ap.add_argument("--criba", default=os.path.join(BASE, "criba_250.csv"),
                    help="salida del cribado, con dG_union por ligando")
    args = ap.parse_args()

    ref = {r["ligand"]: r for r in csv.DictReader(open(args.muestra,
                                                      encoding="utf-8"))}
    filas = []
    for lig, r in ref.items():
        p = os.path.join(args.dir_poses, "ACT_%s_out.pdbqt" % lig)
        if not os.path.exists(p):
            continue
        v = afinidad_vina(p)
        if v is None:
            continue
        filas.append([lig, r["clase"], float(r["pchembl"]), v, None])
    print("muestra: %d, con afinidad de Vina leida: %d" % (len(ref), len(filas)))
    if not filas:
        return 1

    if os.path.exists(args.criba):
        for r in csv.DictReader(open(args.criba, encoding="utf-8")):
            if not r["dG_union"]:
                continue
            for f in filas:
                if f[0] == r["ligand"]:
                    f[4] = float(r["dG_union"])
        print("con dG del MM-GBSA: %d" % sum(1 for f in filas if f[4] is not None))
    else:
        print("aun no hay criba con dG: %s" % args.criba)

    def auroc_de(sel):
        a = [sel(f) for f in filas if f[1] == "activo"]
        i = [sel(f) for f in filas if f[1] == "inactivo"]
        return auroc(a, i)

    print("\nclases: %d activos, %d inactivos"
          % (sum(1 for f in filas if f[1] == "activo"),
             sum(1 for f in filas if f[1] == "inactivo")))

    v, ev = auroc_de(lambda f: f[3])
    print("Vina sola      AUROC %.3f  (error +-%.3f)  <- la referencia"
          % (v, ev if ev else 0.0))

    con = [f for f in filas if f[4] is not None]
    if con:
        m, em = auroc([f[4] for f in con if f[1] == "activo"],
                      [f[4] for f in con if f[1] == "inactivo"])
        print("MM-GBSA solo   AUROC %.3f  (error +-%.3f)"
              % (m, em if em else 0.0))
        # combinado: se prueban pesos y se queda con el mejor, diciendo cual
        mejor, mejor_pes = None, None
        for paso in range(0, 21):
            pes = paso / 4.0
            a, _ = auroc([f[4] - pes * f[3] for f in con if f[1] == "activo"],
                         [f[4] - pes * f[3] for f in con if f[1] == "inactivo"])
            if mejor is None or a > mejor:
                mejor, mejor_pes = a, pes
        print("combinado      AUROC %.3f  con peso %.2f sobre Vina"
              % (mejor, mejor_pes))
        p = pearson([f[2] for f in con], [f[4] for f in con])
        s = pearson(rango([f[2] for f in con]), rango([f[4] for f in con]))
        print("MM-GBSA contra pChEMBL: Pearson %+.3f  Spearman %+.3f"
              % (p if p else 0, s if s else 0))
        pv = pearson([f[2] for f in con], [f[3] for f in con])
        sv = pearson(rango([f[2] for f in con]), rango([f[3] for f in con]))
        print("Vina    contra pChEMBL: Pearson %+.3f  Spearman %+.3f"
              % (pv if pv else 0, sv if sv else 0))
    else:
        print("\n todavia no hay dG del MM-GBSA para comparar")
    return 0


if __name__ == "__main__":
    sys.exit(main())
