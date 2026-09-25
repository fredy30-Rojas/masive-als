# -*- coding: utf-8 -*-
"""Auditoría de la métrica de química (20 sep 2026).

La pregunta, y es la que decide si la lista de 179 vale algo:

    ¿La similitud ECFP4 con los activos mide que la molécula SE UNE a TDP-43,
    o mide sólo que la molécula es un alcaloide plano y catiónico de la
    familia de la berberina?

Motivo de la duda: de los 8 activos conocidos, SEIS son la misma familia
(berberrubina, berberina, sanguinarina, epiberberina, coptisina, nitidina).
Los otros dos no lo son: ketoconazol (azol antifúngico) y PE859 (inhibidor de
agregación, piridil-pirazol). Si la métrica funciona únicamente porque seis de
los ocho positivos son iguales entre sí, el AUC alto es una tautología y la
lista es "análogos de berberina", no candidatos.

Lo que se mide, todo contra el MISMO fondo aleatorio de 200 de antes
(semilla 20260917), sin gastar GPU:

  1. Percentil de cada control por cada medida.
  2. Descomposición: ¿cuánto aporta CADA activo por separado? Si un solo
     compuesto (p. ej. la berberrubina) da casi tanto AUC como los ocho
     juntos, la métrica es "¿es una berberina?".
  3. Reglas tontas de comparación: contar anillos aromáticos, mirar si hay
     catión permanente, medir el peso. Si una regla de dos líneas empata con
     la huella de 2048 bits, la huella no aporta nada que no estuviera en la
     etiqueta "familia".
  4. Prueba de generalidad cruzada: referencia SÓLO la familia plana y ver si
     encuentra al ketoconazol y al PE859 (y al revés).

Uso: python auditar_similitud_tdp43.py
"""
import csv
import json
import os
import random
import statistics
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, DataStructs, Descriptors
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

BASE = r"C:\Users\Fredy\masive-als"
ANALISIS = os.path.join(BASE, "analysis")
RL = os.path.join(ANALISIS, "rescoring_local")
LIGLIB = os.path.join(BASE, "gpu_dock", "libreria_ligands")
CACHE = os.path.join(RL, "smiles_panel_chembl.json")
CONTROLES = os.path.join(RL, "lista_controles_tdp43.csv")
TABLA = os.path.join(RL, "puntuacion_interaccion.csv")
FOCALIZADA = os.path.join(RL, "rescoring_focalizada_arn.csv")
SALIDA = os.path.join(RL, "auditoria_similitud_tdp43.json")

SEMILLA = 20260917
N_ALEATORIOS = 200
EXCLUIDOS = {"cepharanthine", "cepharanthine_limpio", "cepharantina"}

# La familia plana y catiónica: seis de los ocho activos
PLANOS = {"sanguinarine", "coptisine", "berberrubine", "BERBERINE",
          "nitidine", "epiberberine"}
OTROS = {"KETOCONAZOLE", "PE859"}


def huella(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=2048) if m else None


def auc(a, b):
    """AUC con empates a 0,5, como en similitud_activos.py (a = positivos)."""
    if not a or not b:
        return None
    gana = 0.0
    for x in a:
        for y in b:
            gana += 0.5 if x == y else (1.0 if x > y else 0.0)
    return gana / (len(a) * len(b))


def percentil(valor, fondo):
    peores = sum(1 for v in fondo if v is None or valor > v)
    return 100.0 * peores / len(fondo)


def rasgos_tontos(smi):
    """Lo mínimo que uno diría de esta familia química, sin huella ninguna."""
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None:
        return None
    n_arom = rdMolDescriptors.CalcNumAromaticRings(m)
    cation = 1 if any(a.GetFormalCharge() > 0 for a in m.GetAtoms()) else 0
    pesados = m.GetNumHeavyAtoms()
    sp2 = sum(1 for a in m.GetAtoms()
              if a.GetIsAromatic() or any(b.GetBondTypeAsDouble() == 2.0
                                          for b in a.GetBonds()))
    return {
        "anillos_aromaticos": n_arom,
        "cation": cation,
        # El paquete de dos líneas: plano (muchos anillos) y con carga
        "plano_cation": n_arom + 3 * cation,
        "atomos_pesados": pesados,
        "peso": round(Descriptors.MolWt(m), 1),
        "sp2_por_pesado": sp2 / pesados if pesados else 0.0,
    }


def main():
    cache = json.load(open(CACHE, encoding="utf-8")) if os.path.exists(CACHE) else {}

    # --- controles, con su SMILES y su apellido (familia o no) ---
    ctr = {}
    for r in csv.DictReader(open(CONTROLES, encoding="utf-8")):
        n, s = (r.get("ligand") or "").strip(), (r.get("smiles") or "").strip()
        if n and s and n not in EXCLUIDOS:
            ctr[n] = s
    if not ctr:
        print("No encontré los controles.")
        return 1

    # --- fondo aleatorio, el mismo de siempre ---
    random.seed(SEMILLA)
    todos = sorted(f[:-6] for f in os.listdir(LIGLIB) if f.endswith(".pdbqt"))
    aleatorios = random.sample(todos, N_ALEATORIOS)

    panel, sin_smiles = {}, []
    for n in list(ctr) + aleatorios:
        smi = ctr.get(n) or cache.get(n)
        fp = huella(smi)
        if fp is None:
            sin_smiles.append(n)
            continue
        panel[n] = {"fp": fp, "smi": smi, "rasgos": rasgos_tontos(smi)}
    fondo = [a for a in aleatorios if a in panel]
    conocidos = [c for c in ctr if c in panel]
    print("Con huella: %d controles, %d fondo (de %d aleatorios; sin SMILES: %d)"
          % (len(conocidos), len(fondo), len(aleatorios), len(sin_smiles)))
    familia = [c for c in conocidos if c in PLANOS]
    otros = [c for c in conocidos if c in OTROS]
    print("Familia plana: %d -> %s" % (len(familia), ", ".join(familia)))
    print("Otros: %d -> %s\n" % (len(otros), ", ".join(otros)))

    def sim(fp, refs):
        return DataStructs.BulkTanimotoSimilarity(fp, [panel[n]["fp"] for n in refs]) or [0.0]

    # ============================================================
    print("=" * 74)
    print("1. ¿CUÁNTO APORTA CADA ACTIVO POR SEPARADO?")
    print("=" * 74)
    medidas = {}
    medidas["sim_max_8"] = {n: max(sim(panel[n]["fp"], [c for c in conocidos if c != n]))
                            for n in panel}
    medidas["sim_max_planos"] = {n: max(sim(panel[n]["fp"], [c for c in familia if c != n]))
                                 for n in panel}
    for c in familia:
        medidas["sim_a_" + c.lower()] = {n: (sim(panel[n]["fp"], [c])[0])
                                         for n in panel}
    for c in otros:
        medidas["sim_a_" + c.lower()] = {n: (sim(panel[n]["fp"], [c])[0])
                                         for n in panel}
    # El mínimo esfuerzo posible: un solo compuesto de referencia
    print("%-46s %s" % ("medida", "AUC (positivos = los 8 activos)"))
    resumen = {}
    for nombre, vals in medidas.items():
        a = [vals[c] for c in conocidos]
        b = [vals[f] for f in fondo]
        v = auc(a, b)
        resumen[nombre] = round(v, 4) if v else None
        print("%-46s %.3f" % (nombre, v))

    print("\n  Desglose de esa AUC, por grupos de positivos:")
    for nombre in ("sim_max_8", "sim_max_planos"):
        vals = medidas[nombre]
        f_solo = auc([vals[c] for c in familia], [vals[f] for f in fondo])
        o_solo = auc([vals[c] for c in otros], [vals[f] for f in fondo])
        print("    %-16s -> familia sola %.3f | ketoconazol+PE859 solos %s"
              % (nombre, f_solo, ("%.3f" % o_solo) if o_solo else "n/d"))

    # ============================================================
    print("\n" + "=" * 74)
    print("2. REGLAS TONTAS: ¿la huella aporta algo que no esté en la etiqueta?")
    print("=" * 74)
    rasgos = {}
    for campo in ("anillos_aromaticos", "cation", "plano_cation", "peso",
                  "atomos_pesados", "sp2_por_pesado"):
        vals = {n: panel[n]["rasgos"][campo] for n in panel if panel[n]["rasgos"]}
        rasgos[campo] = vals
        a = [vals[c] for c in conocidos if c in vals]
        b = [vals[f] for f in fondo if f in vals]
        v = auc(a, b)
        print("%-46s %.3f" % (campo, v))
        resumen["regla_" + campo] = round(v, 4) if v else None

    # ============================================================
    print("\n" + "=" * 74)
    print("3. PERCENTIL DE CADA CONTROL EN EL FONDO ALEATORIO")
    print("=" * 74)
    print("%-15s %-7s %7s %7s %7s" % ("control", "grupo", "sim8", "sim_plan", "anillos+cat"))
    detalle = {}
    for c in sorted(conocidos, key=lambda n: -medidas["sim_max_8"][n]):
        r = panel[c]["rasgos"]
        detalle[c] = {
            "grupo": "familia" if c in PLANOS else "otro",
            "sim_max_8": round(medidas["sim_max_8"][c], 3),
            "pct_sim_max_8": round(percentil(medidas["sim_max_8"][c],
                                             [medidas["sim_max_8"][f] for f in fondo]), 1),
            "sim_planos": round(medidas["sim_max_planos"][c], 3),
            "pct_sim_planos": round(percentil(medidas["sim_max_planos"][c],
                                              [medidas["sim_max_planos"][f] for f in fondo]), 1),
            "pct_plano_cation": round(percentil(rasgos["plano_cation"][c],
                                                [rasgos["plano_cation"][f] for f in fondo]), 1),
        }
        print("%-15s %-7s %7.2f %7.2f %10.0f%%"
              % (c, detalle[c]["grupo"], detalle[c]["sim_max_8"],
                 detalle[c]["sim_planos"], detalle[c]["pct_plano_cation"]))

    # ============================================================
    print("\n" + "=" * 74)
    print("4. GENERALIDAD CRUZADA (¿sirve fuera de su familia?)")
    print("=" * 74)
    a = [max(sim(panel[c]["fp"], [x for x in familia if x != c])) for c in otros]
    b = [max(sim(panel[f]["fp"], familia)) for f in fondo]
    print("Referencia SÓLO familia -> ¿encuentra al ketoconazol y al PE859?  AUC %.3f"
          % auc(a, b))
    a2 = [sim(panel[c]["fp"], otros)[0] for c in familia]
    b2 = [max(sim(panel[f]["fp"], otros)) for f in fondo]
    print("Referencia SÓLO ketoconazol+PE859 -> ¿encuentra a los 6 planos?  AUC %.3f"
          % auc(a2, b2))

    # ============================================================
    print("\n" + "=" * 74)
    print("5. ¿Y EN LA LISTA DE 179, QUÉ ES LO QUE SUBE?")
    print("=" * 74)
    if os.path.exists(FOCALIZADA):
        foc = [r for r in csv.DictReader(open(FOCALIZADA, encoding="utf-8"))
               if r.get("mmgbsa_dG") not in ("", "None")]
        for r in foc:
            r["dg"] = float(r["mmgbsa_dG"])
        foc.sort(key=lambda r: r["dg"])
        n_f = len(foc)
        # ¿Cuántos de los 179 tendrían la etiqueta tonta de plano+catioónico?
        etiquetados = 0
        for r in foc:
            rr = rasgos_tontos(r.get("smiles") or "")
            if rr and rr["plano_cation"] >= 6:
                etiquetados += 1
        print("De los %d de la lista, con la regla tonta (>=3 anillos arom. y catión): %d (%.0f%%)"
              % (n_f, etiquetados, 100.0 * etiquetados / n_f))
        print("Top 5 por energía, con su etiqueta tonta:")
        for r in foc[:5]:
            rr = rasgos_tontos(r.get("smiles") or "")
            vals = {k: round(v, 2) for k, v in (rr or {}).items() if k != "peso"}
            print("   %-14s dG=%7.2f  %s" % (r["ligand"], r["dg"], vals))

    json.dump({"resumen": resumen, "controles": detalle}, open(SALIDA, "w", encoding="utf-8"),
              indent=2, ensure_ascii=False)
    print("\nGuardado: %s" % SALIDA)
    return 0


if __name__ == "__main__":
    sys.exit(main())
