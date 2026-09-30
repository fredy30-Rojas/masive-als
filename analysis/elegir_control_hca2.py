#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Elige el co-cristal de CONTROL de hCA2, con el mismo criterio pericial que se uso
para `andr` (el unico control de redocking que ha pasado en este proyecto).

POR QUE
-------
El 30 de septiembre se audito la serie de hCA2 antes de acoplar: 1.200 activos
medidos, 281 esqueletos, 44 independientes, MCS de 4 atomos. La serie esta limpia y
el AUC va a poder interpretarse. Ahora falta lo otro: que el RECEPTOR y la CAJA esten
bien, y eso solo se demuestra con un co-cristal cuyo ligando se pueda volver a
acoplar en su sitio.

hCA2 tiene 1.135 entradas PDB, y NO se puede coger la primera. El fallo que mato al
control de CDK2 fue de pericia del cristal, no del motor: dos copias del ligando al
50 %, factores B de 47 y un brazo amonio que no tocaba ni proteina ni agua. Con 1.135
cristales a mano, escoger sin criterio es esperar a que salga un CDK2.

POR QUE EN DOS FASES
--------------------
La pericia del cristal (distancias, factores B, aguas) obliga a bajar el .cif.gz de
cada entrada, y son ~250 KB cada uno: bajar 1.135 son 280 MB para nada. Pero el
filtro QUIMICO (que el farmaco sea sulfonamida, rigido y de tamano normal) sale de la
API sin bajar un byte. Asi que primero se criba por quimica sobre las 1.135, y solo
bajan las coordenadas de las que pasan. Es el mismo orden que el que se uso con DUD-E:
primero la criba de rigidez sobre los 102 ligandos, despues el anclaje de los que
ganaron.

LOS CRITERIOS
-------------
FASE 1 (sin bajar coordenadas, con la API de RCSB y PDBe):
  1. **Que el farmaco sea sulfonamida.** Es lo que une al zinc, y sin interaccion con
     el metal ni la caja ni la puntuacion significan nada en hCA2. Ademas la
     auditoria dio 16,9 % de sulfonamidas en la serie: la subfamilia es pequena, y el
     control tiene que salir de ahi para que el AUC no premie solo "sulfonamida".
  2. **RIGIDEZ**: 0-2 enlaces rotatorios.
  3. **TAMANO**: 12-45 atomos pesados, para que quepa en la caja de 24 A.

FASE 2 (bajando el .cif de las que pasaron):
  4. **El zinc esta en el archivo y el N del sulfonamida lo coordina** (< 2,6 A). El
     mismo caso que el zinc de SOD1, que ya se resolvio: si el metal cae al limpiar
     el receptor, todo lo medido carece de sentido.
  5. **ANCLAJE**: cada N/O del ligando con un N/O de proteina a menos de 3,5 A. Un
     polar sin pareja es un brazo suelto; en CDK2 se iba a 6-9 A.
  6. **SIN AGUAS PUENTE**: pocas aguas a menos de 3,5 A del ligando. Los cristales de
     hCA2 vienen empapados (6t4o tiene 388 aguas en el archivo) y ese total no dice
     nada; lo que decide es si el ligando depende de un puente de agua, porque Vina
     no lo va a predecir. Jarabes: ninguna. Cristal con puente: el control falla por
     una razon que no es del motor.
  7. **FACTOR B bajo**: ligando bien ordenado. B alto es ocupacion parcial y pose
     ambigua, que es media mitad del problema de CDK2.
  8. **UN SOLO LIGANDO, UNA SOLA CONFORMACION**: una sola copia y ninguna
     conformacion alternativa del farmaco. Dos copias al 50 % mataron a CDK2.

QUE DESCARTA Y QUE NO
---------------------
Esto NO dice que el redocking vaya a pasar. Dice que el cristal es el adecuado para
comprobarlo. La prueba de fuego es acoplar, y esa va con la GPU.

Salida:
    analysis/calibracion_hca2/control_hca2_candidatos.csv
    analysis/calibracion_hca2/control_hca2.txt

Uso:
    python elegir_control_hca2.py                       # las 1.135, fase 1 + fase 2
    python elegir_control_hca2.py --fase 1              # solo la criba quimica
    python elegir_control_hca2.py --offset 0 --top 200   # por tandas
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import time
import urllib.request
from collections import defaultdict

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

try:
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict
except ImportError:  # pragma: no cover
    MMCIF2Dict = None

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(BASE, "calibracion_hca2")
CACHE = os.path.join(DIR, "_cache_hca2")
CSV_SALIDA = os.path.join(DIR, "control_hca2_candidatos.csv")
INFORME = os.path.join(DIR, "control_hca2.txt")
CACHE_JSON = os.path.join(DIR, "_cache_api.json")

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
PDBE = "https://www.ebi.ac.uk/pdbe/api"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/%s"
RCSB_NPE = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/%s/%s"
RCSB_COMP = "https://data.rcsb.org/rest/v1/core/chemcomp/%s"
PDB_CIF = "https://files.rcsb.org/download/%s.cif.gz"

TARGET = "CHEMBL205"
# El H explicito importa: `[NX3H]` es un nitrogeno con UN hidrogeno, y la
# sulfonamida primaria (la que mas une al Zn) lleva DOS. Con el patron primero, la
# acetazolamida y el 5-acetamido-1,3,4-tiadiazol-2-sulfonamida de hCA2 NO
# contaban como sulfonamidas y la criba devolvia cero candidatos. Se vio mirando
# el farmaco de 1a42 en la consola, no de memoria.
SULFONAMIDA = Chem.MolFromSmarts("[SX4](=O)(=O)[NX3;H1,H2]")

# El criterio de rigidez se abre a 4 a proposito, y el motivo esta medido: el
# farmaco de 1a42 (BZU) tiene 7 enlaces rotatorios y sin embargo es una
# sulfonamida CICLICA, con el nitrogeno y el azufre en el mismo anillo. Los brazos
# que se mueven estan en la corona y no en lo que toca el Zn. Adjudicar la rigidez
# aqui con el numero de rotatorios rechazaba justamente los ligandos que hacen la
# union al metal, que es lo que hay que medir. La rigidez que de verdad decide la
# decide la fase 2: ningun polar del ligando sin pareja en la proteina.
MAX_ROTATORIOS = 4
MIN_PESADOS, MAX_PESADOS = 12, 45
CORTE_ZN = 2.6          # N(sulfonamida) - Zn de coordinacion
CORTE_PAREJA = 3.5      # polar del ligando sin N/O de proteina cerca = suelto
MAX_AGUAS = 60
MAX_B = 45.0

_cache_api = {}


def get_json(url, intentos=3):
    """GET a la API, con cache en memoria y en disco.

    La cache en disco no es un lujo: son 1.135 entradas por 3 llamadas cada una, y sin
    cache una reejecucion se baja 3.400 peticiones otra vez para volver a lo mismo.
    """
    if url in _cache_api:
        return _cache_api[url]
    for i in range(intentos):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json",
                                                       "User-Agent": "masive-als"})
            with urllib.request.urlopen(req, timeout=60) as r:
                d = json.loads(r.read().decode("utf-8"))
            _cache_api[url] = d
            _cache_disco()
            return d
        except Exception:  # noqa: BLE001
            if i == intentos - 1:
                return None
            time.sleep(1.5 * (i + 1))
    return None


def _cache_disco():
    try:
        with open(CACHE_JSON, "w", encoding="utf-8") as f:
            json.dump(_cache_api, f)
    except Exception:  # noqa: BLE001
        pass


def cargar_cache():
    global _cache_api
    if os.path.exists(CACHE_JSON):
        try:
            with open(CACHE_JSON, encoding="utf-8") as f:
                _cache_api = json.load(f)
        except Exception:  # noqa: BLE001
            _cache_api = {}


def bajar(url, destino, pausa=0.0):
    if os.path.exists(destino) and os.path.getsize(destino) > 0:
        return open(destino, "rb").read()
    for intento in range(3):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "masive-als"})
            with urllib.request.urlopen(req, timeout=120) as r:
                d = r.read()
            os.makedirs(os.path.dirname(destino), exist_ok=True)
            with open(destino, "wb") as f:
                f.write(d)
            time.sleep(pausa)
            return d
        except Exception as e:  # noqa: BLE001
            if intento == 2:
                return b"__FALLO__" + str(e).encode()
            time.sleep(2 * (intento + 1))


def entradas_pdb_hca2():
    """Entradas PDB que ChEMBL asocia a hCA2 (las mismas 1.135 del 29 de septiembre)."""
    d = get_json("%s/target/%s.json" % (CHEMBL, TARGET))
    if not d:
        return []
    out = set()
    for comp in d.get("target_components", []):
        for x in comp.get("target_component_xrefs", []):
            if x.get("xref_src_db") == "PDBe":
                out.add(x["xref_id"].lower())
    return sorted(out)


def ligandos_de_entrada(pdb_id):
    """[(codigo, nombre)] de los farmacos co-cristalizados, por PDBe.

    Se usa `molecule_type = "bound"`, que es el farmaco en el bolsillo, y no el
    "non-polymer" del resumen, que tambien cuenta el zinc, el glicol y los iones.
    """
    d = get_json("%s/pdb/entry/molecules/%s" % (PDBE, pdb_id))
    if not d:
        return None
    out = []
    for e in d.values():
        for ent in e:
            if ent.get("molecule_type") == "bound":
                nombre = ent.get("molecule_name")
                if isinstance(nombre, list):
                    nombre = nombre[0] if nombre else "?"
                sin = ent.get("synonym") or nombre or "?"
                out.append((sin, nombre if isinstance(nombre, str) else "?"))
    return out


def smiles_de_comp(comp_id):
    """SMILES canonico de un componente quimico, por la API de RCSB."""
    d = get_json(RCSB_COMP % comp_id)
    if not d:
        return None, None
    nombre = (d.get("chem_comp") or {}).get("name")
    for desc in d.get("pdbx_chem_comp_descriptor", []) or []:
        if desc.get("type") in ("SMILES_CANONICAL", "SMILES"):
            return nombre, desc["descriptor"]
    return nombre, None


def comps_de_entrada(pdb_id):
    """Codigos de 3 letras de los no-polimeros de una entrada, por RCSB."""
    e = get_json(RCSB_ENTRY % pdb_id.upper())
    if not e:
        return None
    ids = (e.get("rcsb_entry_container_identifiers") or {}).get(
        "non_polymer_entity_ids") or []
    out = []
    for eid in ids:
        n = get_json(RCSB_NPE % (pdb_id.upper(), eid))
        if not n:
            continue
        np_ = n.get("pdbx_entity_nonpoly") or {}
        cid = np_.get("comp_id")
        if cid:
            out.append(cid)
    return out


# ------------------------------------------------------------------ fase 2

def medir_cristal(pdb_id, comp_id):
    """Pericia sobre las coordenadas reales del .cif.

    (n_zn, d_al_zn, sin_pareja, polares, dmin, expuesto, aguas_cerca, b_medio,
     n_copias, n_altloc)

    TRES COSAS QUE HAY QUE HACER BIEN Y QUE SALEN MIRANDO UN CRISTAL, NO DE
    MEMORIA:

    1. **Las conformaciones alternativas (altloc) se tiran.** 5jes tiene 141 atomos
       en conformacion A y 141 en B: el ligando esta modelado en DOS posiciones.
       Si se miden las dos a la vez, la distancia minima entre el ligando y la
       proteina sale 0,25 A, que parece un choque y en realidad son las dos
       conformaciones superpuestas. Se queda con la A. Cuantas alternativas hay se
       cuenta y se dice, porque un ligando con dos conformaciones a la vez es
       justamente lo que mato al control de CDK2 (dos copias al 50 %).
    2. **Las aguas se cuentan las que TOCAN el ligando, no todas.** Los cristales de
       hCA2 vienen empapados: 6t4o tiene 388 aguas en el archivo. Un numero de 388
       no dice nada, porque lo que importa es si el ligando depende de un puente de
       agua, y eso se ve contando las que estan a menos de 3,5 A del ligando.
       Jarabes: ninguna, cristal perfecto. Con agua puente: Vina no la va a predecir
       y el control falla por una razon que no es del motor.
    3. **Las copias se cuentan por cadena, no por conformacion.** Un dimer con dos
       copias del farmaco son dos copias reales, y eso tambien hay que verlo antes de
       prometer el redocking.
    """
    if MMCIF2Dict is None:
        return None, "biopython no esta instalado"
    ruta = os.path.join(CACHE, "%s.cif.gz" % pdb_id)
    datos = bajar(PDB_CIF % pdb_id.upper(), ruta)
    if datos.startswith(b"__FALLO__"):
        return None, "no se bajo el cif"
    try:
        with gzip.open(ruta, "rt", encoding="utf-8", errors="ignore") as h:
            d = MMCIF2Dict(h)
    except Exception as e:  # noqa: BLE001
        return None, "cif ilegible: %r" % e
    if "_atom_site.group_PDB" not in d:
        return None, "cif sin _atom_site"

    g = d["_atom_site.group_PDB"]
    comp = d["_atom_site.label_comp_id"]
    atomo = d["_atom_site.label_atom_id"]
    elem = d.get("_atom_site.type_symbol") or [a[:1] for a in atomo]
    try:
        xyz = np.stack([np.array([float(v) for v in d["_atom_site.Cartn_x"]]),
                        np.array([float(v) for v in d["_atom_site.Cartn_y"]]),
                        np.array([float(v) for v in d["_atom_site.Cartn_z"]])],
                       axis=1)
    except Exception as e:  # noqa: BLE001
        return None, "coordenadas: %r" % e
    bs = d.get("_atom_site.B_iso_or_equiv")
    bfac = np.array([float(v) if v not in (".", "?") else 0.0
                     for v in bs]) if bs else np.zeros(len(g))

    # Los atomos de hidrogeno no cuentan para nada de lo que se mide aqui: ni los
    # choques, ni el anclaje, ni el radio del ligando.
    heavy = np.array([str(e).upper() != "H" for e in elem])

    # 1) un solo modelo y una sola conformacion alternativa
    sel = np.ones(len(g), dtype=bool)
    mod = d.get("_atom_site.pdbx_PDB_model_num")
    if mod:
        sel &= np.array([m == mod[0] for m in mod])
    alt = d.get("_atom_site.label_alt_id")
    n_alt = 0
    if alt:
        lig_alt = sorted({a for a, ci, s in zip(alt, comp, sel)
                         if s and ci == comp_id and a not in (".", "?")})
        n_alt = len(lig_alt)
        # se queda la conformacion A si la hay, y si no la primera. Todas las
        # alternativas del farmaco se cuentan aparte (n_alt) para poder descartar el
        # cristal: un farmaco modelado en dos posiciones a la vez es el problema
        # exacto de CDK2, aqui con otro nombre.
        quedarse = lig_alt[0] if lig_alt else None
        sel &= np.array([(a in (".", "?")) or (quedarse is None or a == quedarse)
                         for a in alt])

    es_prot = np.array([gi == "ATOM" for gi in g]) & sel & heavy
    es_lig = np.array([ci == comp_id for ci in comp]) & sel & heavy
    es_zn = np.array([ci in ("ZN", "ZN2") for ci in comp]) & sel
    es_hoh = np.array([ci in ("HOH", "DOD", "WAT") for ci in comp]) & sel

    P, L, Zn, HOH = xyz[es_prot], xyz[es_lig], xyz[es_zn], xyz[es_hoh]
    if len(L) == 0:
        return None, "el farmaco no esta en el primer modelo"
    if len(P) == 0:
        return None, "sin proteina en el archivo"
    Le = [str(elem[i]).upper() for i in np.where(es_lig)[0]]
    b_medio = float(np.mean(bfac[es_lig]))
    copias = 1
    try:
        asym = d["_atom_site.label_asym_id"]
        copias = len({a for a, m in zip(asym, es_lig) if m})
    except Exception:  # noqa: BLE001
        pass

    dm = np.linalg.norm(L[:, None, :] - P[None, :, :], axis=2)
    dmin = float(dm.min())
    idx_pol = [i for i, e in enumerate(Le) if e in ("N", "O")]
    sin_pareja = sum(1 for i in idx_pol if dm[i].min() > CORTE_PAREJA)
    expuesto = float((dm.min(axis=1) > 4.5).mean())
    d_zn = None
    if len(Zn):
        for i, e in enumerate(Le):
            if e == "N":
                v = float(np.linalg.norm(Zn - L[i], axis=1).min())
                d_zn = v if d_zn is None else min(d_zn, v)
    # 2) aguas que TOCAN el ligando, no todas
    aguas_cerca = (len(HOH) and int((np.linalg.norm(
        HOH[:, None, :] - L[None, :, :], axis=2).min(axis=1) < 3.5).sum())) or 0
    return (len(Zn), d_zn, sin_pareja, len(idx_pol), dmin, expuesto, aguas_cerca,
            b_medio, copias, n_alt), None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fase", type=int, default=2)
    ap.add_argument("--top", type=int, default=0,
                    help="cuantas entradas se criban (0 = todas)")
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--fase2_max", type=int, default=60,
                    help="cuantas de la fase 1 bajan a la pericia de cristal")
    ap.add_argument("--csv", default=None,
                    help="nombre del CSV de salida. Hace falta para correr la criba "
                         "en varios procesos por tramos: si todos escribieran el "
                         "mismo fichero, el ultimo que terminara pisaria a los "
                         "demas y se perderian tres cuartos de las entradas sin "
                         "avisar.")
    ap.add_argument("--desde_csv", default=None,
                    help="se salta la fase 1 y lee los finalistas de este CSV. La "
                         "fase 2 no necesita volver a preguntar a la API lo mismo "
                         "que la fase 1, y con 1.135 entradas eso son casi dos horas "
                         "de reloj por no haber guardado el paso.")
    ap.add_argument("--cache", default=None,
                    help="fichero de cache de la API. Cada proceso por tramo lleva "
                         "el suyo: comparten uno y cada escritura lo pisa con la "
                         "copia que ese proceso tinha al arrancar, y se pierden las "
                         "consultas ya hechas. No da resultados falsos, pero tira "
                         "la cache a la basura.")
    args = ap.parse_args()

    os.makedirs(CACHE, exist_ok=True)
    global CSV_SALIDA, CACHE_JSON, INFORME
    if args.csv:
        CSV_SALIDA = os.path.join(DIR, args.csv)
    if args.cache:
        CACHE_JSON = os.path.join(DIR, args.cache)
    if args.csv:
        INFORME = os.path.join(DIR, args.csv.replace(".csv", ".txt"))
    cargar_cache()
    lineas = []

    def log(m):
        print(m, flush=True)
        lineas.append(m)

    log("ELEGIR EL CO-CRISTAL DE CONTROL DE hCA2")
    log("   %s" % time.strftime("%Y-%m-%d %H:%M"))
    log("   Con 1.135 entradas no se coge la primera: se miden los criterios que ya")
    log("   hicieron pasar a andr, para no repetir el fallo de CDK2.")

    entradas = entradas_pdb_hca2()
    log("   entradas PDB asociadas a hCA2: %d" % len(entradas))
    if not entradas:
        log("   no se pudieron listar (API caida). Se reintenta.")
        return 1
    if args.top:
        lote = entradas[args.offset:args.offset + args.top]
    else:
        lote = entradas[args.offset:]
    log("   cribando %d entradas (offset %d)" % (len(lote), args.offset))

    # ---------------- FASE 1: quimica, sin bajar coordenadas ----------------
    if args.desde_csv:
        log("")
        log("   (modo --desde_csv: se salta la criba y se leen los finalistas de %s)"
            % args.desde_csv)
        with open(args.desde_csv, newline="", encoding="utf-8") as f:
            fase1 = list(csv.DictReader(f))
        log("   %d finalistas desde el CSV" % len(fase1))
        with open(INFORME, "w", encoding="utf-8") as f:
            f.write("\n".join(lineas) + "\n")
        return _fase2(args, log, lineas, fase1)

    log("")
    log("=" * 78)
    log("FASE 1. EL FARMACO ES SULFONAMIDA, RIGIDO Y DE TAMANO NORMAL? (sin bajar")
    log("       ni un .cif: sale de la API de PDBe y de la de RCSB)")
    log("=" * 78)
    fase1, fallan = [], []
    n_sin_bound = n_sin_smiles = 0
    for i, pdb_id in enumerate(lote, 1):
        ligs = ligandos_de_entrada(pdb_id)
        if ligs is None:
            fallan.append((pdb_id, "PDBe no responde"))
            continue
        if not ligs:
            n_sin_bound += 1
            continue
        comps = comps_de_entrada(pdb_id)
        if not comps:
            fallan.append((pdb_id, "RCSB no responde"))
            continue
        mejor = None
        for sin, nombre in ligs:
            for cid in comps:
                if cid.upper() in ("ZN", "ZN2", "HOH", "DOD", "WAT", "GOL", "SO4",
                                   "PO4", "CL", "NA", "K", "MG", "CA", "EDO",
                                   "MES", "ACT", "DMS", "TRS", "IMD", "FMT",
                                   "PEG", "NO3", "CIT", "EPE", "AZI", "BME",
                                   "IPA", "BCT", "NH4", "SCN", "CO3", "FLC"):
                    continue
                cnom, smi = smiles_de_comp(cid)
                if not smi:
                    continue
                m = Chem.MolFromSmiles(smi)
                if m is None:
                    continue
                es_sulf = bool(m.HasSubstructMatch(SULFONAMIDA))
                rot = int(rdMolDescriptors.CalcNumRotatableBonds(m))
                np_ = m.GetNumHeavyAtoms()
                if not es_sulf:
                    continue
                if rot > MAX_ROTATORIOS:
                    continue
                if not (MIN_PESADOS <= np_ <= MAX_PESADOS):
                    continue
                cand = {"pdb": pdb_id, "comp": cid,
                        "nombre": (cnom or sin or "")[:52],
                        "pesados": np_, "rotatorios": rot, "sulfonamida": "SI",
                        "smiles": smi}
                if mejor is None or np_ < mejor["pesados"]:
                    mejor = cand
        if mejor is None:
            n_sin_smiles += 1
            continue
        fase1.append(mejor)
        if i % 40 == 0:
            log("   ...%d/%d  (%d pasan la quimica)" % (i, len(lote), len(fase1)))

    log("")
    log("   %d entradas con farmaco medido y %d con sulfonamida rigida de tamano"
        % (len(lote) - n_sin_bound - len(fallan), len(fase1)))
    log("   %d sin farmaco tipo 'bound' | %d sin sulfonamida que pase | %d fallos"
        % (n_sin_bound, n_sin_smiles, len(fallan)))
    if fallan:
        log("   fallos: %s" % ", ".join("%s (%s)" % f for f in fallan[:6]))
    log("")
    log("   las 25 primeras por tamano (menor = mas encajada en el bolsillo):")
    for r in sorted(fase1, key=lambda x: x["pesados"])[:25]:
        log("   %-6s %-4s %3d pes. %d rot | %s" % (r["pdb"], r["comp"], r["pesados"],
                                                  r["rotatorios"], r["nombre"][:44]))

    with open(CSV_SALIDA, "w", newline="", encoding="utf-8") as f:
        campos = ["pdb", "comp", "nombre", "pesados", "rotatorios", "sulfonamida",
                  "n_zn", "d_al_zn", "sin_pareja", "polares", "dmin", "expuesto",
                  "aguas", "b_medio", "copias", "puntos", "nota", "smiles"]
        w = csv.DictWriter(f, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        for r in fase1:
            w.writerow(r)

    if args.fase < 2 or not fase1:
        with open(INFORME, "w", encoding="utf-8") as f:
            f.write("\n".join(lineas) + "\n")
        log("")
        log("csv: %s" % os.path.basename(CSV_SALIDA))
        return 0

    return _fase2(args, log, lineas, fase1)


def _fase2(args, log, lineas, fase1):
    """Pericia de cristal sobre los finalistas de la fase 1.

    Va en su propia funcion para que el modo --desde_csv pueda llamarla sin repasar
    las 1.135 entradas de la API otra vez: la criba ya esta hecha y escrita en el CSV,
    y volver a preguntar lo mismo son casi dos horas de reloj.
    """
    # ---------------- FASE 2: pericia del cristal ----------------
    log("")
    log("=" * 78)
    log("FASE 2. PERICIA DEL CRISTAL (baja el .cif de las %d primeras)"
        % min(args.fase2_max, len(fase1)))
    log("=" * 78)
    log("%-6s %-4s %3s %6s %6s %7s %6s %6s %5s %5s %4s"
        % ("pdb", "comp", "pts", "nZn", "N-Zn", "sinPareja", "dmin", "agua<3.5",
           "B", "cop", "alt"))
    cribadas = []
    for r in sorted(fase1, key=lambda x: x["pesados"])[:args.fase2_max]:
        m, err = medir_cristal(r["pdb"], r["comp"])
        if m is None:
            r["nota"] = err
            log("%-6s %-4s  no se pudo medir: %s" % (r["pdb"], r["comp"], err))
            continue
        (n_zn, d_zn, sin_pareja, polares, dmin, expuesto, aguas, b_medio, copias,
         n_alt) = m
        puntos = 0
        if n_zn >= 1:
            puntos += 10
        if d_zn is not None and d_zn <= CORTE_ZN:
            puntos += 20
        puntos += max(0, 6 - sin_pareja) * 3
        if aguas <= MAX_AGUAS:
            puntos += 10
        if b_medio <= 30:
            puntos += 10
        elif b_medio <= MAX_B:
            puntos += 5
        if copias <= 1:
            puntos += 10
        if n_alt == 0:
            puntos += 10
        if dmin < 2.5:
            puntos -= 25          # receptor y ligando no son del mismo cristal
        r.update({"n_zn": n_zn,
                  "d_al_zn": "" if d_zn is None else "%.2f" % d_zn,
                  "sin_pareja": sin_pareja, "polares": polares,
                  "dmin": "%.2f" % dmin, "expuesto": "%.2f" % expuesto,
                  "aguas": aguas, "b_medio": "%.1f" % b_medio, "copias": copias,
                  "n_altloc": n_alt, "puntos": puntos, "nota": ""})
        cribadas.append(r)
        log("%-6s %-4s %3d %6d %6s %7s %6.2f %6d %5.0f %5d %4d"
            % (r["pdb"], r["comp"], puntos, n_zn,
               "-" if d_zn is None else "%.2f" % d_zn,
               "%d/%d" % (sin_pareja, polares), dmin, aguas, b_medio, copias,
               n_alt))

    log("")
    log("LAS QUE PASAN TODO")
    log("   sulfonamida + Zn coordinado (N a <= %.1f A) + rigido + ningun polar suelto"
        % CORTE_ZN)
    log("   + <= %d aguas puente + B <= %.0f + una sola copia y una sola conformacion"
        % (MAX_AGUAS, MAX_B))
    # El CSV se lee y vuelve como texto, y comparar "16" con 45 en Python da error.
    # No es un detalle: es la diferencia entre que el filtro funcione o reviente.
    def num(r, k, defecto):
        try:
            return float(r.get(k))
        except (TypeError, ValueError):
            return defecto

    ganan = [r for r in cribadas
             if num(r, "n_zn", 0) >= 1
             and num(r, "d_al_zn", 99) <= CORTE_ZN
             and num(r, "sin_pareja", 99) == 0
             and num(r, "aguas", 999) <= MAX_AGUAS
             and num(r, "b_medio", 999) <= MAX_B
             and num(r, "copias", 99) <= 1
             and num(r, "n_altloc", 0) == 0
             and num(r, "dmin", 99) >= 2.5]
    log("   %d de %d" % (len(ganan), len(cribadas)))
    for r in sorted(ganan, key=lambda x: -num(x, "puntos", 0))[:15]:
        log("   %-6s %3.0f pts | %-4s | %2.0f pes. %.0f rot | N-Zn %.2f A |"
            " %.0f agua puente | B %.0f | %.0f copias"
            % (r["pdb"], num(r, "puntos", 0), r["comp"], num(r, "pesados", 0),
               num(r, "rotatorios", 0), num(r, "d_al_zn", 0), num(r, "aguas", 0),
               num(r, "b_medio", 0), num(r, "copias", 0)))
    if not ganan:
        log("")
        log("   ninguna reune todo. Las 10 mejor puntuadas, con lo que les falla:")
        for r in sorted(cribadas, key=lambda x: -num(x, "puntos", 0))[:10]:
            fallas = []
            if num(r, "n_zn", 0) < 1:
                fallas.append("sin Zn")
            elif num(r, "d_al_zn", 99) > CORTE_ZN:
                fallas.append("N no llega al Zn (%s)" % r.get("d_al_zn", "-"))
            if num(r, "sin_pareja", 0):
                fallas.append("%.0f polar suelto" % num(r, "sin_pareja", 0))
            if num(r, "aguas", 0) > MAX_AGUAS:
                fallas.append("%.0f agua puente" % num(r, "aguas", 0))
            if num(r, "b_medio", 0) > MAX_B:
                fallas.append("B %.0f" % num(r, "b_medio", 0))
            if num(r, "copias", 1) > 1:
                fallas.append("%.0f copias" % num(r, "copias", 1))
            if num(r, "n_altloc", 0) > 0:
                fallas.append("%.0f conformaciones alternativas"
                              % num(r, "n_altloc", 0))
            if num(r, "dmin", 99) < 2.5:
                fallas.append("choque a %.2f A" % num(r, "dmin", 99))
            log("   %-6s %3.0f pts | %-28s | %s"
                % (r["pdb"], num(r, "puntos", 0), r["nombre"][:28],
                   ", ".join(fallas) if fallas else "CUMPLE"))

    with open(CSV_SALIDA, "w", newline="", encoding="utf-8") as f:
        campos = ["pdb", "comp", "nombre", "pesados", "rotatorios", "sulfonamida",
                  "n_zn", "d_al_zn", "sin_pareja", "polares", "dmin", "expuesto",
                  "aguas", "b_medio", "copias", "n_altloc", "puntos", "nota",
                  "smiles"]
        w = csv.DictWriter(f, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        for r in cribadas:
            w.writerow(r)
    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    log("")
    log("csv: %s" % os.path.basename(CSV_SALIDA))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
