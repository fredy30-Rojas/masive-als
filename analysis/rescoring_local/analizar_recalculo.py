# -*- coding: utf-8 -*-
"""Analiza el recálculo del rescoring hecho con el protocolo validado.

Contexto: el 11 sep 2026 se validó el protocolo MM-GBSA y se vio que el dG
bruto NO sirve para ordenar candidatos: premia el tamaño molecular. En SOD1
los controles son más pequeños que el fondo y el AUC bruto sale 0.253 (peor
que el azar); sólo al quitar el tamaño (residual) se ve la señal real (0.754).
En TDP43_v2 ni siquiera con residuales (0.073): no discrimina.

Este script repite ese análisis sobre `rescoring_recalculo_validado.csv`
(242 compuestos pasados con el protocolo validado) y comprueba tres cosas:

  1. Que el recálculo es comparable: misma plataforma y mismo receptor que
     el snapshot congelado del 11 sep (si no, no se puede comparar nada).
  2. Que los controles siguen comportándose igual con el protocolo bueno.
  3. El ranking de candidatos por residual y por eficiencia de ligando,
     que es como SÍ se puede ordenar. Nunca por dG bruto.

Uso:  python analizar_recalculo.py [csv]
"""
import csv
import os
import sys

import numpy as np
from rdkit import Chem

BASE = r"C:\Users\Fredy\masive-als"
if not os.path.isdir(BASE):
    BASE = "/mnt/c/Users/Fredy/masive-als"
LOCAL = os.path.join(BASE, "analysis", "rescoring_local")
DEFAULT = os.path.join(LOCAL, "rescoring_recalculo_validado.csv")
VAL = os.path.join(LOCAL, "validacion_controles.csv")
MANIFEST = os.path.join(LOCAL, "receptores_fijos", "snapshot_20260911",
                        "manifest.json")


def atomos_pesados(smiles):
    """Átomos pesados del fragmento mayor (las sales y los contadores sueltos
    no son el ligando)."""
    if not smiles:
        return 0
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    trozos = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if not trozos:
        return 0
    return max(m.GetNumHeavyAtoms() for m in trozos)


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
    if len(a) < 3:
        return float("nan")
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    return float((ra * rb).sum() / d) if d else 0.0


def cargar(path):
    filas, errores = [], []
    for r in csv.DictReader(open(path, encoding="utf-8")):
        if r.get("mmgbsa_dG"):
            r["dG"] = float(r["mmgbsa_dG"])
            r["ha"] = atomos_pesados(r.get("smiles", ""))
            filas.append(r)
        elif r.get("error"):
            errores.append(r)
    return filas, errores


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    if not os.path.exists(path):
        print("no existe %s" % path)
        return 1

    filas, errores = cargar(path)
    tipos = {(r["target"], r["ligand"]): r["tipo"]
             for r in csv.DictReader(open(VAL, encoding="utf-8"))}
    for r in filas:
        r["tipo"] = tipos.get((r["target"], r["ligand"]), "fondo")

    # ---- 1. ¿Es comparable con la validación? ----
    print("=" * 72)
    print("1. COMPROBACION DE COMPARABILIDAD")
    print("=" * 72)
    print("filas con dG: %d   |   filas con error: %d" % (len(filas), len(errores)))
    print("plataformas: %s" % sorted({r.get("plataforma", "?") for r in filas}))
    for t in sorted({r["target"] for r in filas}):
        hs = sorted({r.get("receptor_md5", "")[:16] for r in filas if r["target"] == t})
        print("  %-10s receptor md5: %s" % (t, hs))
    if os.path.exists(MANIFEST):
        import json
        man = json.load(open(MANIFEST, encoding="utf-8"))
        print("  snapshot 11 sep:    %s" % man)

    # ---- 2. Controles vs fondo ----
    print()
    print("=" * 72)
    print("2. CONTROLES vs FONDO (lo que decide si el orden sirve)")
    print("=" * 72)
    for t in sorted({r["target"] for r in filas}):
        ctrl = [r for r in filas if r["target"] == t and r["tipo"] == "control"]
        fond = [r for r in filas if r["target"] == t and r["tipo"] == "fondo"]
        print("\n=== %s ===  controles n=%d | fondo n=%d" % (t, len(ctrl), len(fond)))
        if not ctrl:
            print("  sin controles en este recálculo: no se puede evaluar aquí")
            continue
        dc = np.array([r["dG"] for r in ctrl])
        df = np.array([r["dG"] for r in fond])
        hc = [r["ha"] for r in ctrl]
        hf = [r["ha"] for r in fond]
        print("  dG controles: mediana %.1f  (%.1f a %.1f)"
              % (np.median(dc), dc.max(), dc.min()))
        print("  dG fondo:     mediana %.1f  (%.1f a %.1f)"
              % (np.median(df), df.max(), df.min()))
        print("  tamaño: controles %.1f atomos pesados | fondo %.1f"
              % (np.mean(hc), np.mean(hf)))
        print("  Spearman dG vs tamaño = %+.3f"
              % spearman([r["dG"] for r in filas if r["target"] == t],
                         [r["ha"] for r in filas if r["target"] == t]))
        print("  1) AUC bruto = %.3f   p = %.4f"
              % (auc(dc, df), mannwhitney_p(dc, df)))
        lo, hi = min(hc), max(hc)
        comp = [r["dG"] for r in fond if lo <= r["ha"] <= hi]
        if len(comp) >= 5:
            print("  2) AUC tamaño-comparable (%d-%d atomos, n=%d) = %.3f"
                  % (lo, hi, len(comp), auc(dc, comp)))
        else:
            print("  2) tamaño-comparable: solo %d del fondo, insuficiente" % len(comp))
        if len(df) >= 8:
            coef = np.polyfit(hf, df, 1)
            res_c = dc - np.polyval(coef, hc)
            res_f = df - np.polyval(coef, hf)
            print("  3) AUC sobre residuales (quitando el tamaño) = %.3f"
                  % auc(res_c, res_f))
            print("     pendiente dG/tamaño = %.2f kcal/mol por atomo pesado" % coef[0])
            print("  4) AUC sobre eficiencia de ligando (dG/atomo) = %.3f"
                  % auc([r["dG"] / r["ha"] for r in ctrl if r["ha"]],
                        [r["dG"] / r["ha"] for r in fond if r["ha"]]))
            # En qué puesto del ranking caen los controles, que es lo que de
            # verdad se pregunta uno: ¿los buenos quedan arriba?
            for nombre, clave in (("dG bruto", "dG"), ("residual", None),
                                  ("eficiencia", "efic")):
                if clave is None:
                    vals = [(r, r["dG"] - np.polyval(coef, r["ha"])) for r in
                            filas if r["target"] == t]
                else:
                    vals = [(r, r["dG"] / r["ha"] if clave == "efic" else r["dG"])
                            for r in filas if r["target"] == t and r["ha"]]
                vals.sort(key=lambda x: x[1])
                puestos = [i + 1 for i, (r, _) in enumerate(vals)
                           if r["tipo"] == "control"]
                if puestos:
                    med = np.median(puestos)
                    print("     por %-10s los controles caen en el percentil %.0f "
                          "(puesto %d de %d)"
                          % (nombre, 100 * (1 - med / len(vals)), med, len(vals)))

    # ---- 3. Ranking de candidatos, como SÍ se puede ordenar ----
    print()
    print("=" * 72)
    print("3. RANKING DE CANDIDATOS (por residual y por eficiencia, NO por dG)")
    print("=" * 72)
    for t in sorted({r["target"] for r in filas}):
        fond = [r for r in filas if r["target"] == t and r["tipo"] == "fondo"]
        todos = [r for r in filas if r["target"] == t]
        if len(fond) < 8:
            continue
        coef = np.polyfit([r["ha"] for r in fond], [r["dG"] for r in fond], 1)
        for r in todos:
            r["resid"] = r["dG"] - np.polyval(coef, r["ha"])
            r["efic"] = r["dG"] / r["ha"] if r["ha"] else 0.0
        print("\n--- %s ---" % t)
        print("  por RESIDUAL (mejor de lo que su tamaño predice):")
        for r in sorted(todos, key=lambda x: x["resid"])[:8]:
            print("    %-34s dG=%7.2f  %3d at.  residual=%+6.2f  [%s]"
                  % (r["ligand"][:34], r["dG"], r["ha"], r["resid"], r["tipo"]))
        print("  por EFICIENCIA DE LIGANDO (dG por atomo pesado):")
        for r in sorted(todos, key=lambda x: x["efic"])[:8]:
            print("    %-34s dG=%7.2f  %3d at.  efic=%+.4f  [%s]"
                  % (r["ligand"][:34], r["dG"], r["ha"], r["efic"], r["tipo"]))

    if errores:
        print()
        print("=" * 72)
        print("4. FILAS CON ERROR (no entraron en el análisis)")
        print("=" * 72)
        for r in errores:
            print("  %-10s %-28s %s" % (r["target"], r["ligand"], r.get("error", "")))
    return 0


if __name__ == "__main__":
    sys.exit(main())
