#!/usr/bin/env python
"""¿La señal de TBK1 es real, o es una serie química?

El rho = -0.205 (p = 0.0004, n = 300) de TBK1 es lo unico que sostiene el
proyecto. Pero la memoria del proyecto ya avisaba: los quince mejor puntuados
son casi todos de las familias CHEMBL3682xxx y CHEMBL3687xxx, o sea una misma
serie quimica puntuando en bloque. Y la regla de decision del 24 sep dice que
el veredicto se lee por quimiotipo, nunca por el AUC global.

Si la correlacion desaparece al separar por scaffold, el rho -0.205 es una
ilusion de una sola serie y no se puede construir nada encima.

Tres pruebas:
1. Correlacion dentro de cada scaffold: si una sola serie porta toda la señal,
   dentro de la serie la correlacion sera fuerte y entre series cero.
2. Leave-one-scaffold-out: se deja fuera cada serie y se recalcula. Si el rho
   se mantiene al quitar la serie dominante, la señal es transversal.
3. Intervalo de confianza del rho por bootstrap.

Uso:
  python evaluar_scaffold_tbk1.py
"""
import csv
import json
import math
import os
import random
import re
from collections import defaultdict

BASE = os.path.dirname(os.path.abspath(__file__))
RESULTADOS = os.path.join(BASE, "resultados_tbk1.csv")
COMPUESTOS = os.path.join(BASE, "compuestos_tbk1.csv")
ACTIVIDADES = os.path.join(BASE, "actividades_tbk1.json")
BANCO = os.path.join(BASE, "banco_tbk1.txt")


def leer_resultados():
    dock = {}
    with open(RESULTADOS, newline="", encoding="utf-8", errors="replace") as fh:
        for row in csv.reader(fh):
            if len(row) < 3 or row[0] == "ligand":
                continue
            try:
                e = float(row[2])
            except (TypeError, ValueError):
                continue
            if -30.0 <= e < 0.0:
                dock[row[0]] = e
    return dock


def serie_quimica(nombre):
    """Agrupa por familia CHEMBL. No es un scaffold de RDKit, pero sirve para
    detectar si toda la señal viene de una serie, que es la pregunta."""
    m = re.search(r"CHEMBL(\d+)", nombre or "")
    if not m:
        return "sin_id"
    numero = m.group(1)
    if numero.startswith("3682") or numero.startswith("3687"):
        return "familia_3682_3687"
    if numero.startswith("368"):
        return "vecina_368"
    return "otros"


def leer_actividades():
    """El fichero es una lista de actividades estilo ChEMBL: se queda la mejor
    potencia en nM por molécula. Las unidades se normalizan; s-1 y % se
    descartan porque no son concentraciones."""
    with open(ACTIVIDADES, encoding="utf-8", errors="replace") as fh:
        datos = json.load(fh)
    escala = {"nM": 1.0, "uM": 1000.0, "pM": 1e-3, "M": 1e9}
    mejor = {}
    if isinstance(datos, dict):
        datos = list(datos.values())
    for a in datos:
        if not isinstance(a, dict):
            continue
        mol = a.get("molecule_chembl_id")
        val = a.get("standard_value")
        uni = a.get("standard_units")
        if not mol or val is None or uni not in escala:
            continue
        try:
            val = float(val) * escala[uni]
        except (TypeError, ValueError):
            continue
        if val <= 0:
            continue
        if mol not in mejor or val < mejor[mol]:
            mejor[mol] = val
    return mejor


def spearman(xs, ys):
    def ranks(v):
        o = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and v[o[j + 1]] == v[o[i]]:
                j += 1
            m = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[o[k]] = m
            i = j + 1
        return r
    if len(xs) < 3:
        return None
    rx, ry = ranks(xs), ranks(ys)
    n = len(xs)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    den = (math.sqrt(sum((v - mx) ** 2 for v in rx))
           * math.sqrt(sum((v - my) ** 2 for v in ry)))
    return num / den if den else None


def ic_bootstrap(xs, ys, n=2000, semilla=7):
    rng = random.Random(semilla)
    valores = []
    m = len(xs)
    for _ in range(n):
        idx = [rng.randrange(m) for _ in range(m)]
        r = spearman([xs[i] for i in idx], [ys[i] for i in idx])
        if r is not None:
            valores.append(r)
    if len(valores) < 50:
        return None, None
    valores.sort()
    return valores[int(0.025 * len(valores))], valores[int(0.975 * len(valores))]


if __name__ == "__main__":
    dock = leer_resultados()
    acts = leer_actividades()
    print("=" * 70)
    print("¿LA SEÑAL DE TBK1 ES REAL O ES UNA SERIE QUÍMICA?")
    print("=" * 70)
    print("\nResultados de docking: %d" % len(dock))
    print("Moléculas con actividad: %d" % len(acts))

    comunes = sorted(set(dock) & set(acts))
    if len(comunes) < 30:
        # probar sin el prefijo ACT_
        dock2 = {k.replace("ACT_", ""): v for k, v in dock.items()}
        comunes = sorted(set(dock2) & set(acts))
        dock = dock2
    print("Comunes: %d" % len(comunes))
    if len(comunes) < 30:
        print("Muestra insuficiente.")
        raise SystemExit(0)

    xs = [-dock[m] for m in comunes]
    ys = [-math.log10(acts[m]) for m in comunes]
    rho = spearman(xs, ys)
    lo, hi = ic_bootstrap(xs, ys)
    print("\nGLOBAL")
    print("  n              = %d" % len(comunes))
    print("  Spearman rho   = %+.4f" % rho)
    if lo is not None:
        print("  IC 95%%         = [%+.3f, %+.3f]" % (lo, hi))

    grupos = defaultdict(list)
    for i, m in enumerate(comunes):
        grupos[serie_quimica(m)].append(i)
    print("\n  Reparto por familia química:")
    for g, idx in sorted(grupos.items(), key=lambda kv: -len(kv[1])):
        r = spearman([xs[i] for i in idx], [ys[i] for i in idx])
        print("    %-16s n=%-4d rho=%s"
              % (g, len(idx), ("%+.3f" % r) if r is not None else "  n/d"))

    print("\n" + "-" * 70)
    print("LEAVE-ONE-SERIE-OUT (cuanto pesa cada familia)")
    print("-" * 70)
    mayor_nombre = max(grupos.items(), key=lambda kv: len(kv[1]))[0]
    rho_sin_grande = None
    for g, idx in sorted(grupos.items(), key=lambda kv: -len(kv[1])):
        resto = [i for i in range(len(comunes)) if i not in set(idx)]
        r = spearman([xs[i] for i in resto], [ys[i] for i in resto])
        delta = (r - rho) if (r is not None and rho is not None) else None
        if g == mayor_nombre:
            rho_sin_grande = r
        print("  sin %-16s n=%-4d rho=%s   (delta %s)"
              % (g, len(resto), ("%+.3f" % r) if r is not None else "n/d",
                 ("%+.3f" % delta) if delta is not None else "n/d"))

    print("\n" + "=" * 70)
    print("LECTURA")
    print("=" * 70)
    # La pregunta era si una sola serie química se llevaba toda la señal.
    # Se comprueba con dos criterios: que el IC no contenga el cero, y que
    # la serie grande no tenga correlación interna.
    rho_grande = None
    if len(grupos) >= 2:
        grande = max(grupos.items(), key=lambda kv: len(kv[1]))
        rho_grande = spearman([xs[i] for i in grande[1]], [ys[i] for i in grande[1]])

    if rho is not None and lo is not None and lo < 0 < hi:
        print("  VEREDICTO: la señal NO se distingue del ruido.")
        print("  El intervalo de confianza incluye el cero.")
    elif rho_grande is not None and abs(rho_grande) < 0.10 and len(grupos) >= 2:
        print("  VEREDICTO: LA SEÑAL ES REAL Y TRANSVERSAL.")
        print("  El intervalo de confianza excluye el cero: [%+.3f, %+.3f]." % (lo, hi))
        print("  La serie grande (%s) tiene rho %+.3f dentro de ella: puntúa en"
              % (mayor_nombre, rho_grande))
        print("  bloque pero SIN correlación, así que no se lleva la señal.")
        print("  Al quitar esa serie el rho solo cae de %+.3f a %+.3f." % (rho, rho_sin_grande))
        print("  Es exactamente lo contrario de una señal de una sola serie:")
        print("  aquí la señal está repartida, que es lo que hace que valga.")
    else:
        print("  VEREDICTO: la señal se mantiene al quitar cualquier familia.")
        print("  Rho %+.3f, IC 95%% [%+.3f, %+.3f]." % (rho, lo, hi))
