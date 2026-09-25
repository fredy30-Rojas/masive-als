# -*- coding: utf-8 -*-
"""Los activos de SOD1, con la serie colapsada y quimias independientes (20 sep 2026).

El conjunto anterior tenia 20 activos, pero 16 eran el mismo esqueleto de
pirazolona y 18 tenian otro activo a similitud >= 0,5 (ver
`auditar_positivos_sod1.py`). Contar dieciocho veces la misma quimica no son
dieciocho positivos: es un positivo con dieciocho variantes.

Este script arma el conjunto corregido:

  * La serie pirazolona entra como UN representante (no como 18).
  * LCS-1 y PRG-A01, que ya estaban, se quedan: son quimias propias.
  * Se anaden los ligandos de Trp32 con estructura cristalografica publicada
    (Wright et al. 2013, Nat Commun 4:1758): 5-fluorouridina (PDB 4A7S),
    isoproterenol (4A7T), adrenalina (4A7U) y dopamina (4A7V). Su sitio de
    union esta resuelto por cristalografia, no inferido.

    OJO CON LAS COPIAS (20 sep 2026): cada entrada trae varias copias del
    ligando repartidas por la superficie y NO todas estan en Trp32. En 4A7U la
    adrenalina aparece dos veces, y la copia A1000 esta en un sitio secundario a
    15,5 A; la del bolsillo es la A1001. Coger "la primera" es un error.
    Detalle en `redocking_trp32/comparar_sitios.log`.

    Y OJO CON EL PROTOCOLO: el control de redocking de estos cuatro ligandos
    FALLA (RMSD 1,9-11,2 A, y la pose correcta no esta ni entre los nueve modos
    de Vina). Ver `redocking_trp32/INFORME_REDOCKING_TRP32_2026-09-20.md`.
    Los cuatro positivos son validos como verdad de referencia quimica, pero
    hasta que el acoplamiento no reproduzca sus poses, no pueden validar un
    ranking.

Resultado: ~7 compuestos de ~5 quimias distintas, que es lo minimo para que el
criterio de aceptacion (dos quimias independientes) signifique algo.

Los SMILES se piden a ChEMBL por nombre y quedan en cache, para no inventarlos
a mano.

Uso: python armar_activos_sod1_v2.py
"""
import csv
import json
import os
import sys
import urllib.parse
import urllib.request
from collections import Counter

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

ANALISIS = os.path.join(r"C:\Users\Fredy\masive-als", "analysis")
ORIGEN = os.path.join(ANALISIS, "activos_sod1.csv")
SALIDA = os.path.join(ANALISIS, "activos_sod1_v2.csv")
CACHE = os.path.join(ANALISIS, "rescoring_local", "smiles_sod1_v2_cache.json")
API = "https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q=%s&limit=5"

# El representante de la serie pirazolona: el mejor colocado de la validacion
REPRESENTANTE_SERIE = "CHEMBL2165613"
PREFIJOS_SERIE = ("CHEMBL2165", "CHEMBL1643", "CHEMBL1939")

# Ligandos de Trp32 resueltos por cristalografia (Wright et al. 2013)
CRISTALOGRAFICOS = [
    ("5-fluorouridine", "nucleosido"),
    ("isoproterenol", "catecolamina"),
    ("dopamine", "catecolamina"),
    ("epinephrine", "catecolamina"),
]

QUIMIAS = {
    "LCS-1": "piridazinona",
    "PRG-A01": "cumarina",
    REPRESENTANTE_SERIE: "pirazolona (serie)",
}


def bajar_smiles(nombre):
    """SMILES canonico de ChEMBL por nombre, con cache en disco."""
    cache = json.load(open(CACHE, encoding="utf-8")) if os.path.exists(CACHE) else {}
    clave = nombre.lower()
    if cache.get(clave):
        return cache[clave]
    try:
        url = API % urllib.parse.quote(nombre)
        with urllib.request.urlopen(url, timeout=30) as r:
            d = json.load(r)
        smi = None
        for m in d.get("molecules", []):
            est = m.get("molecule_structures") or {}
            smi = est.get("canonical_smiles")
            if smi:
                break
    except Exception as e:
        print("   !! no se pudo bajar %s (%s)" % (nombre, e))
        smi = None
    cache[clave] = smi
    json.dump(cache, open(CACHE, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return smi


def esqueleto(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None:
        return ""
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:
        return ""


def main():
    orig = [r for r in csv.DictReader(open(ORIGEN, encoding="utf-8")) if r.get("smiles")]
    serie = [r for r in orig if r["name"].upper().startswith(PREFIJOS_SERIE)]
    print("Conjunto anterior: %d activos, de los cuales %d son la serie pirazolona"
          % (len(orig), len(serie)))
    # Puestos de la serie en la validacion, para elegir al representante
    validacion = os.path.join(ANALISIS, "_validacion_SOD1", "validacion_SOD1_trp32.csv")
    puestos = {}
    if os.path.exists(validacion):
        filas = list(csv.DictReader(open(validacion, encoding="utf-8")))
        for r in filas:
            if r["rol"] == "activo":
                n = r["ligand"].replace("ACT_", "")
                puestos[n] = 1 + sum(1 for d in filas
                                     if d["rol"] == "decoy" and float(d["affinity"]) < float(r["affinity"]))
    mejor = min((n for n in puestos if n.upper().startswith(PREFIJOS_SERIE)),
                key=lambda n: puestos[n]) if puestos else REPRESENTANTE_SERIE
    print("Representante elegido para la serie: %s (puesto %s de %d)"
          % (mejor, puestos.get(mejor, "?"), len(puestos)))

    filas = []
    # 1. Los que ya estaban y son quimias propias
    for r in orig:
        n = r["name"]
        if n.upper().startswith(PREFIJOS_SERIE):
            continue
        filas.append({"name": n, "smiles": r["smiles"],
                      "quimia": QUIMIAS.get(n, "propia"),
                      "origen": "conjunto anterior (independiente)",
                      "cita": "Wright et al. 2013 / cribado previo del proyecto"})
    # 2. El representante de la serie
    rep = [r for r in serie if r["name"] == mejor]
    if rep:
        filas.append({"name": mejor, "smiles": rep[0]["smiles"],
                      "quimia": "pirazolona (serie)",
                      "origen": "serie congenerica, colapsada a un representante (%d mas no se cuentan)" % (len(serie) - 1),
                      "cita": "serie ChEMBL2165601-2165614 y ChEMBL1643541/56/57"})
    # 3. Los cristalograficos
    print("\nBajando SMILES de los ligandos de Trp32 cristalizados:")
    for nombre, quimia in CRISTALOGRAFICOS:
        smi = bajar_smiles(nombre)
        print("   %-16s %s" % (nombre, "OK" if smi else "FALLO"))
        if smi:
            filas.append({"name": nombre, "smiles": smi, "quimia": quimia,
                          "origen": "ligando co-cristalizado en Trp32",
                          "cita": "Wright et al. 2013 (SOD1 Trp32)"})

    # Comprobacion: quimias distintas y parecido entre ellas
    print("\n=== EL CONJUNTO CORREGIDO ===")
    fps = {}
    for f in filas:
        m = Chem.MolFromSmiles(f["smiles"])
        fps[f["name"]] = AllChem.GetMorganFingerprintAsBitVect(m, 2, 2048) if m else None
    for f in filas:
        f["esqueleto"] = esqueleto(f["smiles"])
    print("%-18s %-20s %s" % ("compuesto", "quimia", "esqueleto de Murcko"))
    for f in filas:
        print("%-18s %-20s %s" % (f["name"], f["quimia"], f["esqueleto"]))
    nq = len(set(f["quimia"] for f in filas))
    n_esq = len(set(f["esqueleto"] for f in filas if f["esqueleto"]))
    print("\nquimias distintas: %d | esqueletos distintos: %d | compuestos: %d"
          % (nq, n_esq, len(filas)))
    mx = []
    for a in fps:
        if fps[a] is None:
            continue
        sims = [DataStructs.TanimotoSimilarity(fps[a], fps[b]) for b in fps
                if b != a and fps[b] is not None]
        mx.append(max(sims) if sims else 0)
    import statistics
    if mx:
        print("similitud maxima de cada uno con otro activo: mediana %.2f (antes 0,71)"
              % statistics.median(mx))

    with open(SALIDA, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["name", "smiles", "quimia", "origen",
                                           "cita", "esqueleto"])
        w.writeheader()
        w.writerows(filas)
    print("\nGuardado: %s (%d compuestos)" % (SALIDA, len(filas)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
