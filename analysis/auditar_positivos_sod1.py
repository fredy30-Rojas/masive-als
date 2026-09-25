# -*- coding: utf-8 -*-
"""¿El conjunto de positivos de SOD1 es una serie congenerica? (20 sep 2026).

El pilar del proyecto es SOD1 «validado» con AUC 0,815 frente a señuelos
emparejados en propiedades. Pero nadie había mirado si esos 20 activos son
20 compuestos independientes o un mismo esqueleto repetido, que es el defecto
que la literatura documenta como sesgo de análogos (Chen 2019, sobre DUD-E:
con activos analógicamente próximos y señuelos fáciles, cualquier modelo saca
enriquecimientos de fantasía).

Se mide, con los datos ya calculados (no se acopla nada):

  1. Esqueletos de Murcko repetidos dentro de los 20 activos.
  2. Similitud máxima de cada activo con otro activo.
  3. AUC con todos, AUC con la serie, y **AUC con sólo los que están fuera de
     la serie** — que es la prueba honesta.
  4. El puesto de cada activo dentro del conjunto, para ver quién sube.

Uso: python auditar_positivos_sod1.py
"""
import csv
import os
import statistics
import sys
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ANALISIS = os.path.join(r"C:\Users\Fredy\masive-als", "analysis")
ACTIVOS = os.path.join(ANALISIS, "activos_sod1.csv")
VALIDACION = os.path.join(ANALISIS, "_validacion_SOD1", "validacion_SOD1_trp32.csv")
# Prefijos de los identificadores que forman la serie congenerica (IDs seguidos
# = un mismo trabajo de quimica medicinal)
PREFIJOS_SERIE = ("CHEMBL2165", "CHEMBL1643", "CHEMBL1939")


def huella(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m else None


def auc(a, b):
    """AUC con empates a 0,5. Menos afinidad = mejor, asi que a<b es acierto."""
    if not a or not b:
        return None
    g = 0.0
    for x in a:
        for y in b:
            g += 0.5 if x == y else (1.0 if x < y else 0.0)
    return g / (len(a) * len(b))


def main():
    filas = list(csv.DictReader(open(VALIDACION, encoding="utf-8")))
    for r in filas:
        r["aff"] = float(r["affinity"])
        r["nombre"] = r["ligand"].replace("ACT_", "").replace("DEC_", "")
    act = [r for r in filas if r["rol"] == "activo"]
    dec = [r for r in filas if r["rol"] == "decoy"]
    print("En la validacion: %d activos, %d señuelos" % (len(act), len(dec)))

    # --- 1 y 2: ¿son una serie? ---
    rows = [r for r in csv.DictReader(open(ACTIVOS, encoding="utf-8")) if r.get("smiles")]
    fps = {r["name"]: huella(r["smiles"]) for r in rows}
    fps = {k: v for k, v in fps.items() if v}
    print("\n=== 1. ESQUELETOS DE MURCKO REPETIDOS DENTRO DE LOS ACTIVOS ===")
    esq = Counter()
    for r in rows:
        m = Chem.MolFromSmiles(r["smiles"])
        if m:
            esq[MurckoScaffold.MurckoScaffoldSmiles(mol=m)] += 1
    for k, v in esq.most_common(8):
        print("   %2dx  %s" % (v, k))
    mx = [max([DataStructs.TanimotoSimilarity(fps[a], fps[b]) for b in fps if b != a] or [0])
          for a in fps]
    print("\n   Similitud maxima de cada activo con otro activo:")
    print("   mediana %.2f | media %.2f | minimo %.2f" % (
        statistics.median(mx), statistics.mean(mx), min(mx)))
    print("   Activos que tienen otro activo a similitud >= 0,5: %d de %d"
          % (sum(1 for v in mx if v >= 0.5), len(mx)))

    # --- 3: el AUC, con y sin la serie ---
    serie = [r for r in act if r["nombre"].upper().startswith(PREFIJOS_SERIE)]
    otros = [r for r in act if r not in serie]
    df = [r["aff"] for r in dec]
    print("\n=== 2. EL AUC, CON LA SERIE Y SIN ELLA ===")
    print("   activos totales %d = serie %d + fuera de la serie %d"
          % (len(act), len(serie), len(otros)))
    print("   AUC con TODOS:                         %.3f" % auc([r["aff"] for r in act], df))
    print("   AUC con SOLO la serie:                 %.3f  (n=%d)"
          % (auc([r["aff"] for r in serie], df), len(serie)))
    if otros:
        print("   AUC con SOLO los de fuera de serie:    %.3f  (n=%d)   <-- la prueba honesta"
              % (auc([r["aff"] for r in otros], df), len(otros)))

    # --- 4: puestos ---
    print("\n=== 3. QUIEN SUBE ===")
    for r in sorted(act, key=lambda x: x["aff"]):
        puesto = 1 + sum(1 for d in dec if d["aff"] < r["aff"])
        pct = 100.0 * sum(1 for d in dec if d["aff"] < r["aff"]) / len(dec)
        fuera = "  <-- FUERA de la serie" if r in otros else ""
        print("   puesto %3d/%-3d  %-14s %8.3f  mejor que el %3.0f%% del fondo%s"
              % (puesto, len(filas), r["nombre"], r["aff"], pct, fuera))
    puesto = [1 + sum(1 for d in dec if d["aff"] < r["aff"]) for r in serie]
    print("\n   La serie ocupa los puestos: min %d, mediana %d, max %d de %d"
          % (min(puesto), int(statistics.median(puesto)), max(puesto), len(filas)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
