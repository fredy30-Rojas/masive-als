# -*- coding: utf-8 -*-
"""La vía que no se había probado: buscar por QUÍMICA PARECIDA a los activos (17 sep 2026).

Estado de la investigación:
  - La energía (Vina) contra fondo aleatorio: AUC 0,66. Contra el fondo sesgado
    del recálculo parece un desastre, pero eso es la prueba, no el método.
  - MM-GBSA está dominado por el tamaño (−1,26 kcal/mol por átomo pesado).
  - Cambiar la caja al sitio de unión de ARN no cambia nada (0,657 → 0,655).
  - Contar contactos es medir tamaño otra vez; el apilamiento aromático informa
    algo, pero al combinarlo con la energía empeora (mete ruido de tamaño).
  - Y el detalle que da la pista: con la afinidad, los que suben son los
    alcaloides planos y catiónicos (sanguinarina percentil 88, coptisina 83,
    berberrubina 79) y el único que se cae es el ketoconazol, que no es de esa
    familia. O sea: **el sitio de ARN reconoce una química concreta**.

Hipótesis: en un caso así, con pocos activos conocidos y una química clara, la
búsqueda por SIMILITUD a los activos supera al acoplamiento, y además se combina
bien con él (son dos informaciones distintas: forma química y complementariedad).

Aquí se mide, contra el mismo fondo aleatorio de 200 moléculas. Los SMILES de
esas 200 se piden a ChEMBL por su identificador y quedan en caché.

Uso: python similitud_activos.py
"""
import csv
import json
import os
import random
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs

RDLogger.DisableLog("rdApp.*")   # RDKit escribe mucho por consola y aquí molesta

BASE = r"C:\Users\Fredy\masive-als"
ANALISIS = os.path.join(BASE, "analysis")
LIGLIB = os.path.join(BASE, "gpu_dock", "libreria_ligands")
TABLA = os.path.join(ANALISIS, "rescoring_local", "puntuacion_interaccion.csv")
CACHE = os.path.join(ANALISIS, "rescoring_local", "smiles_panel_chembl.json")
SALIDA = os.path.join(ANALISIS, "rescoring_local", "similitud_activos.json")
SEMILLA = 20260917
N_ALEATORIOS = 200
EXCLUIDOS = {"cepharantina", "cepharanthine", "cepharanthine_limpio"}
# Los activos de la familia plana y catiónica (sin el ketoconazol, que es de otra química)
PLANOS = {"sanguinarine", "coptisine", "berberrubine", "BERBERINE", "nitidine", "epiberberine"}


def smiles_locales():
    """SMILES de los controles, de las listas del proyecto."""
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


def bajar_smiles(nombres):
    """SMILES de ChEMBL por identificador, con caché en disco."""
    cache = {}
    if os.path.exists(CACHE):
        cache = json.load(open(CACHE, encoding="utf-8"))
    faltan = [n for n in nombres if not cache.get(n)]
    print("SMILES en caché: %d | por bajar: %d" % (len(cache) - len(faltan), len(faltan)))

    def uno(nombre):
        if not nombre.upper().startswith("CHEMBL"):
            return nombre, None
        url = "https://www.ebi.ac.uk/chembl/api/data/molecule/%s.json" % nombre
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                d = json.load(r)
            est = d.get("molecule_structures") or {}
            return nombre, (est.get("canonical_smiles") or est.get("standard_smiles"))
        except Exception:
            return nombre, None

    if faltan:
        with ThreadPoolExecutor(max_workers=8) as ex:
            for nombre, smi in ex.map(uno, faltan):
                cache[nombre] = smi
        json.dump(cache, open(CACHE, "w", encoding="utf-8"), ensure_ascii=False)
    return cache


def huella(smi):
    mol = Chem.MolFromSmiles(smi) if smi else None
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def main():
    filas = [r for r in csv.DictReader(open(TABLA, encoding="utf-8"))]
    controles = [r["ligand"] for r in filas if r["grupo"] == "control"]
    random.seed(SEMILLA)
    todos = sorted(f[:-6] for f in os.listdir(LIGLIB) if f.endswith(".pdbqt"))
    aleatorios = random.sample(todos, N_ALEATORIOS)

    locales = smiles_locales()
    bajar = [n for n in aleatorios + controles if n not in locales]
    chembl = bajar_smiles(bajar)

    panel = {}
    for n in aleatorios + controles:
        smi = locales.get(n) or chembl.get(n)
        h = huella(smi) if smi else None
        if h is not None:
            panel[n] = h
    print("Con huella química: %d de %d (controles: %d de %d)"
          % (len(panel), len(aleatorios) + len(controles),
             sum(1 for c in controles if c in panel), len(controles)))

    def similitud(fp, nombre):
        """Similitud con los activos conocidos, SIN contarse a sí mismo.

        Importante: si un control se compara consigo mismo sale 1,0 y el AUC se
        infla solo. Cada molécula se mide contra los activos menos ella.
        """
        todos_ref = [n for n in panel if n in controles and n != nombre]
        planos_ref = [n for n in todos_ref if n in PLANOS]
        s_todos = DataStructs.BulkTanimotoSimilarity(fp, [panel[n] for n in todos_ref]) or [0.0]
        s_plan = DataStructs.BulkTanimotoSimilarity(fp, [panel[n] for n in planos_ref]) or [0.0]
        return {
            "sim_max": max(s_todos),
            "sim_media_top3": sum(sorted(s_todos, reverse=True)[:3]) / min(3, len(s_todos)),
            "sim_planos": max(s_plan),
        }

    # Afinidad de la caja de ARN, de la tabla ya calculada
    afinidades = {}
    for r in filas:
        if r.get("arn_afinidad"):
            try:
                afinidades[r["ligand"]] = float(r["arn_afinidad"])
            except ValueError:
                pass

    datos = {}
    for n, fp in panel.items():
        d = dict(similitud(fp, n))
        d["grupo"] = "control" if n in controles else "fondo"
        d["afinidad"] = afinidades.get(n)
        datos[n] = d

    def auc(campo, grupo_a, grupo_b):
        a = [datos[n][campo] for n in grupo_a if n in datos and datos[n].get(campo) is not None]
        b = [datos[n][campo] for n in grupo_b if n in datos and datos[n].get(campo) is not None]
        if not a or not b:
            return None, 0, 0
        gana = 0.0
        for x in a:
            for y in b:
                gana += 0.5 if x == y else (1.0 if x > y else 0.0)
        return gana / (len(a) * len(b)), len(a), len(b)

    ct = [c for c in controles if c in datos]
    fo = [a for a in aleatorios if a in datos]
    print("\n============ SIMILITUD CONTRA EL FONDO ALEATORIO ============")
    resumen = {"controles": len(ct), "fondo": len(fo)}
    for campo, desc in (("sim_max", "similitud máxima con cualquier activo conocido"),
                        ("sim_media_top3", "media de las 3 similitudes más altas"),
                        ("sim_planos", "similitud máxima con la familia plana")):
        a, nc, nf = auc(campo, ct, fo)
        if a is None:
            continue
        resumen[campo] = round(a, 4)
        print("%-42s AUC %.3f  (n=%d controles, %d fondo)" % (desc, a, nc, nf))
    # Referencia: la afinidad sola, invertida porque menos es mejor
    af_ct = [-datos[c]["afinidad"] for c in ct if datos[c].get("afinidad") is not None]
    af_fo = [-datos[a]["afinidad"] for a in fo if datos[a].get("afinidad") is not None]
    if af_ct and af_fo:
        gana = sum(1.0 if x > y else (0.5 if x == y else 0.0) for x in af_ct for y in af_fo)
        print("%-42s AUC %.3f  (n=%d controles, %d fondo)" % (
            "afinidad de acoplamiento (referencia)", gana / (len(af_ct) * len(af_fo)),
            len(af_ct), len(af_fo)))
        resumen["afinidad"] = round(gana / (len(af_ct) * len(af_fo)), 4)

    # Dónde cae cada control por similitud
    vals = [datos[a]["sim_planos"] for a in fo]
    print("\nPercentil de cada control por similitud con la familia plana:")
    for c in sorted(ct, key=lambda n: -datos[n]["sim_planos"]):
        peores = sum(1 for v in vals if datos[c]["sim_planos"] > v)
        print("   %-16s percentil %3.0f  (sim %.2f, afinidad %s)" % (
            c, 100.0 * peores / len(vals), datos[c]["sim_planos"],
            datos[c]["afinidad"] if datos[c]["afinidad"] is not None else "n/d"))

    json.dump(resumen, open(SALIDA, "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    print("\nGuardado: %s" % SALIDA)
    return 0


if __name__ == "__main__":
    sys.exit(main())
