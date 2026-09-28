"""Audita la verdad de referencia de TBK1: sesgo de "el mejor valor".

QUE PREGUNTA CONTESTA
---------------------
El banco usa, para cada compuesto, SU MEJOR pChEMBL: el mas favorable de todos
los que tiene. Con eso, un compuesto que alguna vez se ha visto activo es
"activo", y uno que solo se ha medido una vez y ha salido mal es "inactivo".
Eso no es verdad de laboratorio: es un artefacto de como se eligio. Y como el
docking solo separa con un AUROC de 0,60, hay que saber cuanto de ese 0,60 es
del metodo y cuanto de la lista.

QUE MIDE
  1. cuantos compuestos hay segun cuantos valores independientes tienen
  2. el contraste entre "mejor valor" (lo que usa el banco) y "valor tipico"
     (la mediana de las medidas de union), que es lo que se usaria si se
     quisiera una verdad limpia
  3. cuantos sobreviven si se exigen dos medidas por debajo de 1 uM para
    _llamar activo_, y dos por encima de 10 uM para llamar inactivo
  4. que tipos de medida se estan mezclando en la misma escala

POR QUE TIPOS DISTINTOS NO SE PUEDEN MEZCLAR: un IC50 de inhibicion
competitiva y una constante de disociacion miden cosas distintas, con
supuestos distintos, y un inhibidor de baja afinidad con velocidad de
asociacion lenta tiene un IC50 malo y un kd bueno. Ponerlos en la misma
columna y quedarse con el mejor de todos es elegir siempre el dato que mas
favorece a la hypothesis de que el compuesto es activo.
"""
import collections
import csv
import json
import os
import statistics
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
ACTIVIDADES = os.path.join(BASE, "actividades_tbk1.json")
CSV_BANCO = os.path.join(BASE, "compuestos_tbk1.csv")

# Medidas de UNION, que son las que dicen "se une a la proteina".
UNION = ("IC50", "Ki", "Kd", "EC50")
# Las que NO son comparables con las anteriores, y por tanto no entran en la
# escala: cineticas (distintas unidades, distinto significado), medidas
# celulares ("Inhibition", "Residual Activity", "Activity") y porcentajes.
NO_UNION = ("kon", "k_off", "Inhibition", "Residual Activity", "Activity",
            "Kd ", "Potency")


def cargar():
    return json.load(open(ACTIVIDADES, encoding="utf-8"))


def pchembl_de(a):
    """El campo se llama pchembl_value, NO pchembl (28 sep 2026, medido: con
    `a.get("pchembl")` salen las 13.123 actividades con None y la auditoria
    entera da cero compounds sin decir por que)."""
    try:
        return float(a.get("pchembl_value"))
    except (TypeError, ValueError):
        return None


def main():
    acts = cargar()
    print("actividades crudas: %d" % len(acts))

    # 1. tipos de medida y cuanto de cada una
    tipos = collections.Counter()
    por_compuesto = collections.defaultdict(list)
    relaciones = collections.Counter()
    sin_pchembl = 0
    for a in acts:
        t = a.get("standard_type")
        p = pchembl_de(a)
        relaciones[a.get("standard_relation") or "="] += 1
        if p is None:
            sin_pchembl += 1
            continue
        tipos[t] += 1
        por_compuesto[a.get("molecule_chembl_id")].append(
            (t, p, a.get("standard_relation") or "="))
    print("sin pChEMBL utilizable: %d" % sin_pchembl)
    print("relacion de la medida: %s   (< y > son LIMITES, no medidas: "
          "contarlas como si fueran exactas infla los inactivos)"
          % dict(relaciones))
    print("tipos de medida: %s" % dict(tipos.most_common(10)))

    n_union = sum(v for k, v in tipos.items() if k in UNION)
    print("de UNION (IC50/Ki/Kd/EC50): %d  |  el resto: %d"
          % (n_union, sum(tipos.values()) - n_union))

    # 2. cuantos compuestos y cuantos valores independientes tienen
    print("\ncompuestos con pChEMBL: %d" % len(por_compuesto))
    hist = collections.Counter(len(v) for v in por_compuesto.values())
    print("valores por compuesto:")
    for k in sorted(hist)[:8]:
        print("   %2d medidas: %5d compuestos" % (k, hist[k]))
    uno = hist.get(1, 0)
    print("   con UNA sola medida: %d (%.0f%% de los compuestos)"
          % (uno, 100.0 * uno / len(por_compuesto)))

    # 3. contraste: mejor valor (lo que usa el banco) frente a mediana
    def mejor(cid):
        return max(p for _, p, _ in por_compuesto[cid])

    def mediana(cid):
        return statistics.median(p for _, p, _ in por_compuesto[cid])

    banco = {}
    for r in csv.DictReader(open(CSV_BANCO, encoding="utf-8")):
        banco[r["molecule_chembl_id"]] = r["franja"]
    comunes = [c for c in banco if c in por_compuesto]
    print("\ncompuestos del banco con actividad: %d" % len(comunes))

    act_mejor = [c for c in comunes if banco[c] != "peor que 10 micromolar"]
    ina_mejor = [c for c in comunes if banco[c] == "peor que 10 micromolar"]
    act_med = [c for c in act_mejor if mediana(c) >= 6.0]
    ina_med = [c for c in ina_mejor if mediana(c) < 6.0]
    print("con el MEJOR valor (como el banco): %d activos, %d inactivos"
          % (len(act_mejor), len(ina_mejor)))
    print("con la MEDIANA de sus medidas:       %d activos, %d inactivos"
          % (len(act_med), len(ina_med)))
    cambiados = len(act_mejor) - len(act_med)
    print("activos que dejan de serlo con la mediana: %d (%.0f%%)"
          % (cambiados, 100.0 * cambiados / max(1, len(act_mejor))))

    # 4. verdad estricta: dos medidas por debajo de 1 uM, o dos por encima
    #    de 10 uM, y todas de tipo de UNION
    est_act, est_inact = [], []
    for cid in comunes:
        vals = [p for t, p, _ in por_compuesto[cid] if t in UNION]
        if len(vals) < 2:
            continue
        if sum(1 for p in vals if p >= 6.0) >= 2:
            est_act.append(cid)
        elif sum(1 for p in vals if p < 6.0) >= 2:
            est_inact.append(cid)
    print("\nverdad ESTRICTA (solo IC50/Ki/Kd/EC50, 2 medidas por clase):")
    print("   activos:   %d" % len(est_act))
    print("   inactivos: %d" % len(est_inact))
    if est_act and est_inact:
        print("   es utilizable como terreno: SI, hay %d y %d"
              % (len(est_act), len(est_inact)))
    else:
        print("   NO es utilizable: falta una de las dos clases")

    # 5. cuanto pesa cada tipo de medida en la escala del banco
    print("\nde donde sale el pChEMBL del banco, por tipo:")
    cuenta = collections.Counter()
    for cid in comunes:
        mejor_p = mejor(cid)
        for t, p, _ in por_compuesto[cid]:
            if p == mejor_p:
                cuenta[t] += 1
                break
    for t, n in cuenta.most_common(10):
        print("   %-18s %5d  (%.0f%%)"
              % (t, n, 100.0 * n / len(comunes)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
