# -*- coding: utf-8 -*-
"""La lista focalizada: misma química que los activos, ya acoplada (17 sep 2026).

Lo que se midió antes, contra un fondo aleatorio de 200 moléculas de la librería:

    similitud con los activos conocidos        AUC 0,885
    media de las 3 similitudes más altas       AUC 0,902
    similitud con la familia plana             AUC 0,923
    afinidad de acoplamiento (lo de antes)     AUC 0,656

Es decir: para esta diana, la química pesa más que la energía. Tiene sentido y
tiene un motivo: el sitio de unión de ARN del TDP-43 es una superficie plana que
reconoce alcaloides planos y catiónicos (berberina, sanguinarina, coptisina,
berberrubina, nitidina, epiberberina son todos la misma familia de
bencilisoquinolinas). El acoplamiento prefiere al que más hueco rellena; la
similitud reconoce la familia.

Este script explota eso y devuelve algo usable: de toda la librería que tiene
SMILES, las moléculas de esa misma química, cruzadas con su afinidad ya
calculada. No gasta GPU: la afinidad ya está, de la corrida de la librería.

**Aviso honesto de lo que esto es y no es:** enriquece dentro de una familia
química conocida, así que sirve para ampliar la serie (análogos) y para rescatar
fármacos ya existentes de esa familia; NO descubre quimiotipos nuevos. Y la
afinidad que se muestra arrastra el sesgo de tamaño ya medido: sirve para
desempatar dentro de la familia, no como criterio único.

Uso: python lista_focalizada.py [umbral_similitud]
"""
import csv
import json
import os
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs, Descriptors

RDLogger.DisableLog("rdApp.*")

BASE = r"C:\Users\Fredy\masive-als"
ANALISIS = os.path.join(BASE, "analysis")
LIB_SMILES = os.path.join(ANALISIS, "_lib_consolidada.csv")
TOTAL = os.path.join(BASE, "gpu_dock", "resultados_libreria", "resultados_libreria_total.csv")
SALIDA = os.path.join(ANALISIS, "rescoring_local", "lista_focalizada_tdp43v2.csv")
CACHE_FICHA = os.path.join(ANALISIS, "rescoring_local", "fichas_chembl.json")
TARGET = "TDP43_v2"
# Como el acoplamiento, para desempatar dentro de la familia
UMBRAL = float(sys.argv[1]) if len(sys.argv) > 1 else 0.25
# La familia que reconoce el sitio: alcaloides planos y catiónicos
PLANOS = {"sanguinarine", "coptisine", "berberrubine", "BERBERINE", "nitidine", "epiberberine"}
EXCLUIDOS = {"cepharantina", "cepharanthene", "cepharanthine", "cepharanthine_limpio"}


def smiles_activos():
    out = {}
    for ruta, campos in ((os.path.join(ANALISIS, "controles_calibracion.csv"),
                          ("ligand_libreria", "control")),
                         (os.path.join(ANALISIS, "rescoring_local", "lista_controles_tdp43.csv"),
                          ("ligand",))):
        if not os.path.exists(ruta):
            continue
        for r in csv.DictReader(open(ruta, encoding="utf-8")):
            for c in campos:
                n = (r.get(c) or "").strip()
                s = (r.get("smiles") or "").strip()
                if n and s:
                    out.setdefault(n, s)
    return out


def huella(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m is not None else None


def fichas_chembl(nombres):
    """Nombre comercial y fase clínica, para saber qué es cada acierto."""
    cache = json.load(open(CACHE_FICHA, encoding="utf-8")) if os.path.exists(CACHE_FICHA) else {}
    faltan = [n for n in nombres if n not in cache]

    def uno(n):
        try:
            with urllib.request.urlopen(
                    "https://www.ebi.ac.uk/chembl/api/data/molecule/%s.json" % n, timeout=25) as r:
                d = json.load(r)
            return n, {"nombre": d.get("pref_name"), "fase": d.get("max_phase"),
                       "peso": d.get("molecule_properties", {}).get("full_mwt")}
        except Exception:
            return n, {}

    if faltan:
        with ThreadPoolExecutor(max_workers=8) as ex:
            for n, f in ex.map(uno, faltan):
                cache[n] = f
        json.dump(cache, open(CACHE_FICHA, "w", encoding="utf-8"), ensure_ascii=False)
    return cache


def main():
    activos = smiles_activos()
    ref = {n: huella(activos[n]) for n in PLANOS if n in activos}
    ref = {n: h for n, h in ref.items() if h is not None}
    print("Activos de la familia plana usados como referencia: %d -> %s" % (len(ref), ", ".join(ref)))
    if not ref:
        print("No encontré los SMILES de la familia plana.")
        return 1

    # Afinidad ya calculada (acoplamiento de la caja actual)
    afin = {}
    with open(TOTAL, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["target"] == TARGET:
                try:
                    afin[r["ligand"]] = float(r["affinity"])
                except ValueError:
                    pass
    print("Afinidades disponibles: %d" % len(afin))

    filas, descartados, sin_huella = [], 0, 0
    for r in csv.DictReader(open(LIB_SMILES, encoding="utf-8")):
        nombre, smi = (r.get("ligand") or "").strip(), (r.get("smiles") or "").strip()
        if not nombre or not smi or nombre in EXCLUIDOS:
            continue
        fp = huella(smi)
        if fp is None:
            sin_huella += 1
            continue
        sims = DataStructs.BulkTanimotoSimilarity(fp, [ref[n] for n in ref])
        sim = max(sims)
        if sim < UMBRAL:
            continue
        mol = Chem.MolFromSmiles(smi)
        # Esqueleto químico (primera parte del InChIKey): sirve para no contar dos
        # veces la misma molécula por venir en forma de sal distinta. En la
        # librería hay berberina, cloruro de berberina y sulfato de berberina: es
        # el mismo compuesto tres veces.
        esq = ""
        if mol is not None:
            try:
                esq = Chem.MolToInchiKey(mol).split("-")[0]
            except Exception:
                esq = ""
        filas.append({
            "ligand": nombre,
            "similitud_activos": round(sim, 3),
            "afinidad_tdp43v2": afin.get(nombre),
            "atomos_pesados": mol.GetNumHeavyAtoms() if mol else "",
            "peso_molecular": round(Descriptors.MolWt(mol), 1) if mol else "",
            "esqueleto": esq,
            "smiles": smi,
        })
        descartados += 1
    print("Moleculas que pasan el umbral: %d | sin huella (SMILES raro): %d" % (descartados, sin_huella))
    print("Por encima del umbral de similitud %.2f: %d" % (UMBRAL, len(filas)))

    # Orden: primero las más parecidas; y entre parecidas, la afinidad como desempate
    filas.sort(key=lambda x: (-x["similitud_activos"],
                              x["afinidad_tdp43v2"] if x["afinidad_tdp43v2"] is not None else 0))
    # Sin repetir el mismo compuesto vestido de sal: se queda la mejor de cada
    # esqueleto (la de mayor similitud, que es como está ordenada la lista).
    unicos, vistos = [], set()
    for f in filas:
        e = f["esqueleto"] or f["ligand"]
        if e in vistos:
            continue
        vistos.add(e)
        unicos.append(f)
    print("Esqueletos distintos: %d de %d filas (las sales y duplicados no cuentan dos veces)"
          % (len(unicos), len(filas)))

    cabecera = ["ligand", "similitud_activos", "afinidad_tdp43v2", "atomos_pesados",
                "peso_molecular", "esqueleto", "nombre_chembl", "fase", "smiles"]

    # Ficha de los primeros: qué es cada uno (fármaco conocido o no)
    primeros = [f["ligand"] for f in filas[:60] if f["ligand"].upper().startswith("CHEMBL")]
    fichas = fichas_chembl(primeros) if primeros else {}
    for f in filas:
        fic = fichas.get(f["ligand"], {})
        f["nombre_chembl"] = fic.get("nombre") or ""
        f["fase"] = fic.get("fase") or ""

    # La lista buena es la de esqueletos únicos; la completa se guarda aparte
    with open(SALIDA, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cabecera)
        w.writeheader()
        w.writerows(unicos)
    completa = SALIDA.replace(".csv", "_completa.csv")
    with open(completa, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cabecera)
        w.writeheader()
        w.writerows(filas)
    print("Guardado: %s (%d esqueletos distintos)" % (SALIDA, len(unicos)))
    print("Completa, con sales y duplicados: %s (%d filas)" % (completa, len(filas)))

    conocidos = [f for f in unicos if f["fase"] not in ("", None)]
    print("\nDe esos %d esqueletos, ya son fármacos con fase clínica: %d" % (len(unicos), len(conocidos)))
    print("\nLas 15 primeras (misma química que los activos, con afinidad como desempate):")
    for f in unicos[:15]:
        print("  %-14s sim %.2f | afinidad %-6s | %-8s %s" % (
            f["ligand"], f["similitud_activos"], f["afinidad_tdp43v2"],
            ("fase " + str(f["fase"])) if f["fase"] else "nuevo", f["nombre_chembl"] or ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
