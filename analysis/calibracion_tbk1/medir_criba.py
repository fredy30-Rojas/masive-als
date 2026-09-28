"""Dice si el dG del MM-GBSA separa activos de inactivos, y si correlaciona.

Un solo cuadro, sin optimizacion: eso sesga el VALOR, pero no deberia sesgar
el RANKING, que es lo que se usa. Por eso lo que se mide aqui no es si el
numero es bueno, sino si ORDENA bien.
"""
import csv
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))


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
    return sxy / (sxx ** 0.5 * syy ** 0.5)


def rango(xs):
    """Rangos con empates por la media."""
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


def spearman(xs, ys):
    return pearson(rango(xs), rango(ys))


def auroc(activo, inactivo):
    """Proporcion de pares (activo mas negativo que inactivo) bien ordenados.
    Con signo: mas bajo = mejor unido."""
    if not activo or not inactivo:
        return None
    n = len(activo) * len(inactivo)
    suma = 0.0
    for a in activo:
        for i in inactivo:
            suma += 1.0 if a < i else (0.5 if a == i else 0.0)
    return suma / n


def main(ruta):
    filas = []
    for r in csv.DictReader(open(ruta, encoding="utf-8")):
        if not r["dG_union"]:
            print("  (sin dG: %s  %s)" % (r["ligand"], r["nota"]))
            continue
        filas.append((r["ligand"], r["clase"], float(r["pchembl"]),
                      float(r["dG_union"])))
    if not filas:
        print("sin datos")
        return 1
    print("con dG: %d de los que hay en el csv" % len(filas))
    act = [f for f in filas if f[1] == "activo"]
    ina = [f for f in filas if f[1] == "inactivo"]
    print("activos %d, inactivos %d\n" % (len(act), len(ina)))

    for nombre, grupo in (("activos", act), ("inactivos", ina)):
        if not grupo:
            continue
        dgs = sorted(g[3] for g in grupo)
        print("%-10s dG  min %8.1f  mediana %8.1f  max %8.1f"
              % (nombre, dgs[0],
                 dgs[len(dgs) // 2] if len(dgs) % 2
                 else (dgs[len(dgs) // 2 - 1] + dgs[len(dgs) // 2]) / 2,
                 dgs[-1]))

    a = auroc([g[3] for g in act], [g[3] for g in ina])
    if a is not None:
        print("\nAUROC (activo mejor que inactivo): %.3f" % a)
        print("  0,5 es azar. Por encima de 0,7 ya ordena algo.")
    if len(act) + len(ina) > 2:
        p = pearson([g[2] for g in filas], [g[3] for g in filas])
        s = spearman([g[2] for g in filas], [g[3] for g in filas])
        if p is not None:
            print("\nPearson  dG contra pChEMBL: %+.3f" % p)
            print("Spearman dG contra pChEMBL: %+.3f" % s)
            print("  Si fuera negativo, el mas activo seria el mas negativo.")
    print("\nLos diez mas unidos (dG mas bajo):")
    for f in sorted(filas, key=lambda x: x[3])[:10]:
        print("  %-16s %-9s pChEMBL %5.2f  dG %8.1f"
              % (f[0], f[1], f[2], f[3]))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1
                  else os.path.join(BASE, "criba_20_resultado.csv")))
