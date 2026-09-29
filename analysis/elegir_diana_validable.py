#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""elegir_diana_validable.py — elegir la diana donde el embudo se pueda VALIDAR.

POR QUE ESTE SCRIPT
-------------------
El 29 de septiembre de 2026 se midio que las validaciones de TDP-43 y SOD1 no tienen
poder estadistico: con 7 y 11 positivos, el IC 95% del AUC tiene mas de 0,4 de ancho y
el criterio que decide sale con p de 0,18 a 0,81. Con esos numeros no se puede decir
nada, ni del embudo ni de la funcion de puntuacion.

La conclusion tecnica es que la diana, no el metodo, es el problema: TDP-43 no tiene
estructura con farmaco pequeno ni serie de potencia medida. Este script busca una diana
donde la verdad de referencia **si exista**, que es lo unico que convierte el proyecto
en falsable.

LOS DOS REQUISITOS, Y POR QUE LOS DOS
-------------------------------------
1. **Activos MEDIDOS.** Se cuentan las actividades con `pchembl_value` (potencia
   estandarizada, en -log10 de molar) y tipo IC50, Ki o Kd. El   `pchembl` importa: sin el, un numero sin unidades comparables no sirve, y el
   proyecto ya tiene ejemplos de referencias sesgadas (16 de 20 "activos" de SOD1 eran
   el mismo nucleo pirazolona).
2. **Un farmaco co-cristalizado.** Sin estructura con ligando no hay control de
   redocking, y sin control no se sabe si el receptor o la caja estan bien. Se
   comprueba de verdad: se leen las entradas PDB que el propio registro de ChEMBL asocia
   a la diana y se mira si el resumen de PDBe trae un ligando de polimero no
   proteico.

LO QUE NO HACE
--------------
No se fia de la memoria. Los nombres de la lista corta son HIPOTESIS de partida (los
que se nos ocurren como dianas caracterizadas) y TODO lo que se cuenta sale de la API:
si una no cumple, se cae. La lista corta son 12; hacer las 102 de DUD-E seria el
siguiente paso si esto no sale.

Uso:
    python elegir_diana_validable.py
"""
import csv
import json
import os
import time
import urllib.parse
import urllib.request

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
PDBE = "https://www.ebi.ac.uk/pdbe/api"

# Minimo de activos medidos para que el margen del AUC baje a algo utilizable. Sale de
# `bootstrap_auc.py`: con n+ = 40 el semiancho ronda 0,10, y con n+ = 7 es 0,23.
MIN_ACTIVOS = 50
# Cuantas entradas PDB se miran por diana para judgear si hay farmaco co-cristalizado.
MAX_ENTRADAS = 6

# Lista CORTA de partida. Son nombres, no hechos: el script comprueba cada uno.
CANDIDATOS = [
    "Carbonic anhydrase 2",
    "Acetylcholinesterase",
    "Cyclooxygenase-2",
    "Carbonic anhydrase 1",
    "Acetylcholinesterase (human)",
    "Chymotrypsin",
    "Trypsin",
    "Cyclooxygenase-1",
    "Acyl-protein thioesterase 1",
    "Beta-secretase",
    "Choline esterase",
    "Aldose reductase",
]

SALIDA_CSV = os.path.join(BASE, "dianas_validables.csv")


def log(m):
    print(m, flush=True)


def get(url, intentos=3):
    """GET con reintentos y sin cache. Devuelve dict o None."""
    for i in range(intentos):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json",
                                                       "User-Agent": "masive-als"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            if i == intentos - 1:
                log("      (fallo la consulta: %s)" % repr(e)[:90])
                return None
            time.sleep(2 * (i + 1))
    return None


def n_activos(target_id):
    """Actividades con pchembl (potencia estandarizada) de IC50/Ki/Kd."""
    q = urllib.parse.urlencode({
        "target_chembl_id": target_id,
        "pchembl_value__isnull": "false",
        "standard_type__in": "IC50,Ki,Kd",
        "limit": 1})
    d = get("%s/activity.json?%s" % (CHEMBL, q))
    return None if not d else d["page_meta"]["total_count"]


def dianas_por_nombre(nombre):
    """ChEMBL ids de dianas cuyo nombre contiene `nombre`, de humano y proteina unica."""
    q = urllib.parse.urlencode({"q": nombre, "limit": 20})
    d = get("%s/target/search.json?%s" % (CHEMBL, q))
    if not d:
        return []
    out = []
    for t in d.get("targets", []):
        if t.get("target_type") != "SINGLE PROTEIN":
            continue
        if (t.get("organism") or "") != "Homo sapiens":
            continue
        if t.get("pref_name") and nombre.lower() in t["pref_name"].lower():
            out.append((t["target_chembl_id"], t["pref_name"]))
    return out


def entradas_pdb(target_id):
    """Entradas PDB que el registro de la diana asocia (xrefs de PDBe)."""
    d = get("%s/target/%s.json" % (CHEMBL, target_id))  # ChEMBL SI admite .json
    if not d:
        return []
    pdb = set()
    for comp in d.get("target_components", []):
        for x in comp.get("target_component_xrefs", []):
            if x.get("xref_src_db") == "PDBe":
                pdb.add(x["xref_id"].lower())
    return sorted(pdb)


def ligandos_de_pdb(pdb_id):
    """(n_ligandos_cristalizados, sinonimos) de una entrada PDB, segun PDBe.

    OJO con la URL: el endpoint de PDBe NO admite el sufijo `.json` (con el devuelve
    404 y    el script concluye "no hay ligando" sin decir que fallo la consulta). La
    entidad que se busca es la de `molecule_type = "bound"`, que es el farmaco en el
    bolsillo, no el "non-polymer" del resumen, que tambien cuenta el glicol y los iones.
    """
    d = get("%s/pdb/entry/molecules/%s" % (PDBE, pdb_id))
    if not d:
        return None, []  # la consulta fallo: None, no cero, para no contar como "sin ligando"
    n, nombres = 0, []
    for e in d.values():
        for ent in e:
            if ent.get("molecule_type") == "bound":
                n += 1
                sin = ent.get("synonym") or (ent.get("molecule_name") or ["?"])[0]
                nombres.append(sin)
    return n, nombres


def perfil_activos(target_id, n=250):
    """Mediana de peso y flexibilidad de los activos MEDIDOS de una diana.

    POR QUE ESTE SEGUNDO PASO
    -------------------------
    Ser diana validable no basta: tambien tiene que ser una diana *ganable* por un
    cribado. El 29 de septiembre se midio que un fragmento de 11-15 atomos pesados no
    puede ganar a un fondo de unidores de ARN de 29-45, y que el problema es el rango.
    La misma logica se aplica al elegir la diana: si sus activos medidos son, de
    mediana, moleculas de 45 atomos con ocho enlaces rotatorios, el cribado va a
    ordenar por tamaño y no vamos a aprender nada, por bien medida que este la serie.

    Asi que para cada diana candidata se mide el perfil de SU serie: numero de atomos
    pesados, enlaces rotatorios y cLogP de una muestra de activos reales.
    """
    q = urllib.parse.urlencode({
        "target_chembl_id": target_id, "pchembl_value__isnull": "false",
        "standard_type__in": "IC50,Ki,Kd", "limit": n})
    d = get("%s/activity.json?%s" % (CHEMBL, q))
    if not d:
        return None
    smiles, vistos = [], set()
    for a in d.get("activities", []):
        s = a.get("canonical_smiles")
        mol_id = a.get("molecule_chembl_id")
        if not s or mol_id in vistos:
            continue
        vistos.add(mol_id)
        smiles.append(s)
    if not smiles:
        return None
    from rdkit import Chem, RDLogger
    from rdkit.Chem import Crippen, Descriptors
    RDLogger.DisableLog("rdApp.*")
    pesados, rot, logp = [], [], []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        pesados.append(m.GetNumHeavyAtoms())
        rot.append(Descriptors.NumRotatableBonds(m))
        logp.append(Crippen.MolLogP(m))
    if not pesados:
        return None
    return {"n": len(pesados),
            "pesados": float(np.median(pesados)),
            "rot": float(np.median(rot)),
            "logp": float(np.median(logp))}


def main():
    log("DIANA DONDE EL EMBUDO SE PUEDA VALIDAR")
    log("======================================")
    log("Criterio: >= %d activos con pchembl (IC50/Ki/Kd) Y un farmaco"
        " co-cristalizado verificado." % MIN_ACTIVOS)
    log("Todo se cuenta por API; la lista de nombres son solo hipótesis de partida.")
    log("")

    filas, ya_vistos = [], set()
    for nombre in CANDIDATOS:
        for tid, pref in dianas_por_nombre(nombre):
            if tid in ya_vistos:
                continue
            ya_vistos.add(tid)
            log("%-42s %s" % (pref, tid))
            n_act = n_activos(tid)
            if n_act is None:
                continue
            log("   activos medidos (pchembl, IC50/Ki/Kd): %d" % n_act)
            pdb = entradas_pdb(tid)
            con_ligando = []
            fallos_pdb = 0
            for pdb_id in pdb[:MAX_ENTRADAS]:
                n_lig, ids = ligandos_de_pdb(pdb_id)
                if n_lig:
                    con_ligando.append((pdb_id, ids))
                if n_lig is None:
                    fallos_pdb += 1
                time.sleep(0.2)
            log("   entradas PDB asociadas: %d | con farmaco en el bolsillo, de las"
                " %d revisadas: %d%s"
                % (len(pdb), min(MAX_ENTRADAS, len(pdb)), len(con_ligando),
                   "  (%d consultas fallidas)" % fallos_pdb if fallos_pdb else ""))
            for pdb_id, ids in con_ligando[:3]:
                log("      %s: %s" % (pdb_id, ", ".join(str(x)[:38] for x in ids[:3])))
            cumple_act = n_act >= MIN_ACTIVOS
            cumple_cristal = bool(con_ligando)
            log("   -> %s" % ("CUMPLE LOS DOS" if (cumple_act and cumple_cristal)
                              else "descartada (faltan %s)"
                              % ("activos" if not cumple_act else "cristal")))
            filas.append({"target_chembl_id": tid, "nombre": pref,
                          "activos_pchembl": n_act,
                          "entradas_pdb": len(pdb),
                          "con_ligando": len(con_ligando),
                          "cristal": "%s %s" % (con_ligando[0][0],
                                                ",".join(con_ligando[0][1][:4]))
                          if con_ligando else "",
                          "cumple": "SI" if (cumple_act and cumple_cristal) else "no"})
            time.sleep(0.3)
        log("")

    buenas = [f for f in filas if f["cumple"] == "SI"]
    buenas.sort(key=lambda f: -f["activos_pchembl"])

    log("")
    log("=" * 70)
    log("SEGUNDO FILTRO: ¿SON GANABLES? perfil de los activos medidos de cada una")
    log("=" * 70)
    log("Un cribado ordena por tamaño. Si la serie de la diana es grande y flexible,")
    log("el rango se come la señal (lo medido el 29 de septiembre con los fragmentos de")
    log("TDP-43). Se mide la mediana de una muestra de 250 activos reales de cada diana.")
    log("")
    log("%-30s %6s %9s %8s %8s" % ("diana", "muestra", "pesados", "rotat.", "cLogP"))
    for f in buenas:
        p = perfil_activos(f["target_chembl_id"])
        if not p:
            log("%-30s %6s  (no se pudo medir)" % (f["nombre"][:30], "-"))
            continue
        f.update({"muestra": p["n"], "pesados_med": p["pesados"],
                  "rot_med": p["rot"], "logp_med": p["logp"]})
        log("%-30s %6d %9.0f %8.0f %8.2f"
            % (f["nombre"][:30], p["n"], p["pesados"], p["rot"], p["logp"]))
        time.sleep(0.3)
    log("")
    log("Un cribado con Vina ordena de forma sistematica por numero de atomos")
    log("pesados. Las dianas cuya mediana esta por debajo de ~30 pesados y con menos de")
    log("~5 enlaces rotatorios son las que dan un resultado interpretable.")

    campos = ["target_chembl_id", "nombre", "activos_pchembl", "entradas_pdb",
              "con_ligando", "cristal", "cumple", "muestra", "pesados_med",
              "rot_med", "logp_med"]
    with open(SALIDA_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        for f in filas:
            w.writerow(f)

    log("=" * 70)
    log("DIANAS QUE CUMPLEN LOS DOS REQUISITOS")
    log("=" * 70)
    if not buenas:
        log("Ninguna de la lista corta. Habria que ampliar a las 102 de DUD-E, que es")
        log("el paso siguiente natural: todas tienen farmaco co-cristalizado, la")
        log("pregunta es solo cuales tienen serie medida.")
    for f in buenas:
        log("%-42s %s  %6d activos   cristal %s"
            % (f["nombre"], f["target_chembl_id"], f["activos_pchembl"], f["cristal"]))
    log("")
    log("csv: %s" % os.path.basename(SALIDA_CSV))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
