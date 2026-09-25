# -*- coding: utf-8 -*-
"""Analiza la validación del rescoring: ¿distingue lo activo de lo que no?

Diseño: 20 controles positivos conocidos (18 SOD1, 2 TDP-43) contra 180
compuestos de fondo tomados al azar de la propia lista corta (todos
acoplados, todos pasando los mismos filtros). Si el rescoring ordena bien,
los controles deben recibir dG más negativos que el fondo.

PROBLEMA QUE HAY QUE CONTROLAR: el dG de un MM-GBSA de una sola trayectoria
premia el tamaño molecular (más átomos, más contactos, más negativo). Y en
SOD1 los controles son MÁS PEQUEÑOS que el fondo (19,5 vs 31,1 átomos
pesados de media): el sesgo de tamaño juega CONTRA los controles. Por eso
se reportan tres números y hay que leerlos juntos:

  1. AUC bruto — controles vs fondo tal cual.
  2. AUC con fondo de tamaño comparable — solo los compuestos de fondo que
     caen en el rango de tamaño de los controles.
  3. AUC sobre residuales — se ajusta dG ~ átomos pesados con el fondo y se
     comparan los residuos: mide si el control es mejor de lo que su tamaño
     predice. Es la prueba que separa "discrimina" de "solo premia tamaño".

Un resultado convincente es que el control gane por su propio mérito y no
por tamaño: AUC bruto alto, y residual positivo (dG mejor que el previsto
por su tamaño).

Uso:  python analizar_validacion.py [csv]
"""
import csv
import os
import sys

import numpy as np

BASE = r"C:\Users\Fredy\masive-als"
if not os.path.isdir(BASE):
    BASE = "/mnt/c/Users/Fredy/masive-als"
LOCAL = os.path.join(BASE, "analysis", "rescoring_local")
RES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
VAL = os.path.join(LOCAL, "validacion_controles.csv")
DEFAULT = os.path.join(LOCAL, "rescoring_validacion.csv")


def atomos_pesados(pose):
    n = 0
    with open(pose, encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("ENDMDL"):
                break
            if l.startswith(("ATOM", "HETATM")) and l[76:79].strip().upper() != "H":
                n += 1
    return n


def auc(mejores, peores):
    """P(mejor puntúa antes que peor), dG más negativo = mejor."""
    a, b = np.asarray(mejores, float), np.asarray(peores, float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    gana = (a[:, None] < b[None, :]).sum() + 0.5 * (a[:, None] == b[None, :]).sum()
    return float(gana) / (len(a) * len(b))


def mannwhitney_p(a, b):
    try:
        from scipy.stats import mannwhitneyu
        return float(mannwhitneyu(a, b, alternative="two-sided").pvalue)
    except Exception:
        x = np.concatenate([a, b])
        r = np.argsort(np.argsort(x)) + 1.0
        n1, n2 = len(a), len(b)
        u1 = r[:n1].sum() - n1 * (n1 + 1) / 2.0
        mu, sd = n1 * n2 / 2.0, np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12.0)
        if sd == 0:
            return 1.0
        from math import erfc, sqrt
        return float(erfc(abs((u1 - mu) / sd) / sqrt(2)))


def spearman(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else 0.0


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    if not os.path.exists(path):
        print("no existe %s" % path)
        return 1
    filas = [r for r in csv.DictReader(open(path, encoding="utf-8"))
             if r.get("mmgbsa_dG")]
    tipos = {r["ligand"]: r["tipo"] for r in
             csv.DictReader(open(VAL, encoding="utf-8"))}
    for r in filas:
        r["tipo"] = tipos.get(r["ligand"], "?")
        pose = os.path.join(RES, "results_" + r["target"], r["ligand"] + "_out.pdbqt")
        r["ha"] = atomos_pesados(pose) if os.path.exists(pose) else 0
        r["dG"] = float(r["mmgbsa_dG"])

    print("filas con dG: %d | plataformas: %s"
          % (len(filas), sorted({r.get("plataforma", "?") for r in filas})))
    print("hash receptor (debe ser único por diana): %s"
          % sorted({(r["target"], r.get("receptor_md5", "")[:12]) for r in filas}))

    for t in sorted({r["target"] for r in filas}):
        ctrl = [r for r in filas if r["target"] == t and r["tipo"] == "control"]
        fond = [r for r in filas if r["target"] == t and r["tipo"] == "fondo"]
        if not fond:
            continue
        dc = np.array([r["dG"] for r in ctrl])
        df = np.array([r["dG"] for r in fond])
        print("\n=== %s ===  controles n=%d | fondo n=%d" % (t, len(dc), len(df)))
        print("  dG controles: mediana %.1f  (%.1f a %.1f)"
              % (np.median(dc), dc.max(), dc.min()) if len(dc) else "  sin controles")
        print("  dG fondo:     mediana %.1f  (%.1f a %.1f)"
              % (np.median(df), df.max(), df.min()))
        hc = [r["ha"] for r in ctrl]
        hf = [r["ha"] for r in fond]
        print("  tamaño: controles %.1f atomos pesados | fondo %.1f"
              % (np.mean(hc) if hc else 0, np.mean(hf)))
        rho = spearman([r["dG"] for r in filas if r["target"] == t],
                       [r["ha"] for r in filas if r["target"] == t])
        print("  Spearman dG vs tamaño = %+.3f" % rho)
        if not ctrl:
            print("  (sin controles con pose: no se puede calcular AUC)")
            continue
        print("  1) AUC bruto = %.3f   p = %.4f" % (auc(dc, df), mannwhitney_p(dc, df)))
        # 2) fondo de tamaño comparable
        lo, hi = min(hc), max(hc)
        comp = [r["dG"] for r in fond if lo <= r["ha"] <= hi]
        if len(comp) >= 5:
            print("  2) AUC tamaño-comparable (%d-%d atomos, n=%d) = %.3f"
                  % (lo, hi, len(comp), auc(dc, comp)))
        else:
            print("  2) tamaño-comparable: solo %d del fondo, insuficiente" % len(comp))
        # 3) residuales: ajuste dG ~ atomos pesados con el fondo
        if len(df) >= 8:
            coef = np.polyfit([r["ha"] for r in fond], df, 1)
            res_c = dc - np.polyval(coef, hc)
            res_f = df - np.polyval(coef, hf)
            print("  3) AUC sobre residuales (quitando el tamaño) = %.3f"
                  % auc(res_c, res_f))
            print("     pendiente dG/tamaño = %.2f kcal/mol por atomo pesado"
                  % coef[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())
