"""Amplia la verdad de referencia de TBK1, que es lo que bloquea al proyecto.

EL PROBLEMA QUE RESUELVE
------------------------
Con la lista del banco (el mejor pChEMBL de cada compuesto) el docking de Vina
sale con un AUROC de 0,535: azar. Con una lista estricta (dos medidas de union
coherentes) sube a 0,672. O sea, buena parte de ese cero coma cinco era ruido
de nuestra propia lista, no del metodo. Pero la lista estricta deja solo 17
inactivos, y con 17 el margen de error es de seis centesimas: no sirve.

Aqui se sacan mas inactivos de verdad, de dos sitios que la lista del banco
tiraba a la basura sin mirar:

1. IC50 EXACTOS SIN pChEMBL. Hay 3.865 IC50, pero solo 2.642 tienen pChEMBL.
   Los 1.223 que se pierden se reparten en: 518 con un aviso de ChEMBL de
   "los valores parecen un orden de magnitud distintos de los publicados, las
   unidades pueden estar mal" (ESTOS SI HAY QUE TIRARLOS), 631 limites ("<" o
   ">", que no son medidas sino cotas) y **578 medidas exactas en nM sin
  NINGUN aviso**, que se pueden recuperar easy: pChEMBL = 9 - log10(IC50 en nM).
   Medido: 327 con "<" y 304 con ">" dentro de los que no tienen aviso.

2. Kd DERIVADO DE LA CINETICA. Hay 3.028 medidas de kon y 3.028 de k_off. Con
   las dos del mismo compuesto sale un Kd que es una medida de union tan
   valida como un IC50: Kd = koff / kon. Una kon de 10^4 M^-1 s^-1 es un
   inactivo, y asi se ven los compuestos que no unionan, que son justamente
   los que falta.

QUE SE DEVUELVE
  verdad_tbk1_v2.csv   id, smiles, n_medidas, pchembl_mediana, clase
  informe en la salida

Y el criterio de clase es el mismo para todos: activo si la mediana de sus
medidas da >= 6 (o sea 1 uM o mejor), inactivo si da < 6, y se exige un
minimo de medidas para no volver a caer en el problema de las medidas unicas.
"""
import collections
import csv
import json
import math
import os
import statistics
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
ACTIVIDADES = os.path.join(BASE, "actividades_tbk1.json")
SALIDA = os.path.join(BASE, "verdad_tbk1_v2.csv")

UNION = ("IC50", "Ki", "Kd", "EC50")
# ChEMBL avisa de estos: los valores son un orden de magnitud distintos de los
# publicados y las unidades pueden estar mal. No se pueden recuperar.
AVISO_UNIDADES = "order of magnitude different"
MIN_MEDIDAS = 2


def pchembl_ic50_nm(v):
    """pChEMBL de un IC50 en nM: 9 - log10(nM)."""
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    if v <= 0:
        return None
    return 9.0 - math.log10(v)


def medida_de(a):
    """Un valor de union utilizable de esta actividad, o None.

    Devuelve (pchembl, de_donde). Se acepta:
      - pChEMBL tal cual, si lo hay
      - IC50 exacto en nM sin aviso de unidades, calcularlo
      - Ki y Kd exactos en nM, igual
    Y se rechazan los limites (<, >) y los avisos de unidades.
    """
    rel = a.get("standard_relation") or "="
    if rel != "=":
        return None, None
    dv = a.get("data_validity_description") or ""
    if AVISO_UNIDADES in dv:
        return None, None
    p = a.get("pchembl_value")
    if p is not None:
        try:
            return float(p), a.get("standard_type")
        except (TypeError, ValueError):
            return None, None
    t = a.get("standard_type")
    u = a.get("standard_units")
    if t in UNION and u == "nM":
        v = pchembl_ic50_nm(a.get("standard_value"))
        if v is not None and 0.0 < v <= 15.0:
            return v, t + "(nM)"
    return None, None


def kd_derivado(por):
    """Kd = koff / kon, cuando el mismo compuesto tiene las dos."""
    kd = {}
    for cid, acts in por.items():
        kon = kfo = None
        for a in acts:
            if (a.get("standard_relation") or "=") != "=":
                continue
            if a.get("data_validity_description"):
                continue
            if a.get("standard_units") != "nM":
                continue
            v = a.get("standard_value")
            try:
                v = float(v)
            except (TypeError, ValueError):
                continue
            if v <= 0:
                continue
            if a.get("standard_type") == "kon":
                kon = max(kon or 0.0, v)          # 1/nM = 1e9/M
            elif a.get("standard_type") == "k_off":
                kfo = min(kfo or 1e30, v)         # nM/s
        if kon and kfo:
            # kon en 1/nM, koff en nM/s  ->  Kd (nM) = koff/kon
            valor = kfo / kon
            if 0.0 < valor <= 1e7:
                kd[cid] = -math.log10(valor * 1e-9)
    return kd


def main():
    acts = json.load(open(ACTIVIDADES, encoding="utf-8"))
    por = collections.defaultdict(list)
    for a in acts:
        por[a.get("molecule_chembl_id")].append(a)

    medidas, origen = collections.defaultdict(list), collections.Counter()
    for cid, lista in por.items():
        for a in lista:
            p, de = medida_de(a)
            if p is not None and 0.0 < p <= 15.0:
                medidas[cid].append(p)
                origen[de] += 1

    print("medidas utilizables por compuesto: %d" % sum(len(v) for v in medidas.values()))
    print("de donde (cuenta de medidas): %s" % dict(origen.most_common(8)))

    kd = kd_derivado(por)
    print("Kd derivados de la cinetica (kon + k_off): %d compuestos" % len(kd))
    # POR QUE SON CERO, y esta medido que no es un fallo del script: hay
    # 3.028 actividades de kon y 3.028 de k_off en el fichero, y las 6.056
    # tienen `value` y `standard_value` a None. ChEMBL guarda el registro del
    # ensayo (tipo, unidades) pero NO el numero. Son el 46% de las actividades
    # de la diana y no aportan ni una cifra. La cinematica se descarta entera.

    banco = {}
    for r in csv.DictReader(open(os.path.join(BASE, "compuestos_tbk1.csv"),
                                encoding="utf-8")):
        banco[r["molecule_chembl_id"]] = r

    filas = []
    for cid, vals in medidas.items():
        if cid not in banco:
            continue
        n = len(vals)
        med = statistics.median(vals)
        if n < MIN_MEDIDAS:
            continue
        filas.append((cid, banco[cid]["smiles"], n, med,
                      "activo" if med >= 6.0 else "inactivo"))
    for cid, p in kd.items():
        if cid not in banco:
            continue
        filas.append((cid + "_kin", banco[cid]["smiles"], 2, p,
                      "activo" if p >= 6.0 else "inactivo"))

    filas.sort(key=lambda f: -f[3])
    with open(SALIDA, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["molecule_chembl_id", "smiles", "n_medidas",
                    "pchembl_mediana", "clase"])
        for f in filas:
            w.writerow([f[0], f[1], f[2], "%.2f" % f[3], f[4]])

    act = sum(1 for f in filas if f[4] == "activo")
    ina = len(filas) - act
    print("\nVERDAD v2, con %d medidas o mas por compuesto:" % MIN_MEDIDAS)
    print("   activos:   %d" % act)
    print("   inactivos: %d" % ina)
    print("   (antes: 1552 activos y 141 inactivos con el mejor valor; "
          "346 y 17 con el criterio estricto de solo pChEMBL)")
    print("escrito %s" % SALIDA)
    return 0


if __name__ == "__main__":
    sys.exit(main())
