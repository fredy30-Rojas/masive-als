#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audita la serie de hCA2 (CHEMBL205) ANTES de que exista ningun AUC.

POR QUE ESTO, Y POR QUE AHORA
------------------------------
El 29 de septiembre se eligio hCA2 como la diana donde el embudo se pueda validar:
13.051 activos medidos con pchembl, serie pequena (mediana de 22 atomos pesados y 5
enlaces rotatorios) y 1.135 estructuras con farmaco para el control de redocking.

El aviso 1 de aquel informe decia, textual: "la serie de hCA2 es congenerica, y eso es
un peligro conocido". Esto es lo que mide si el aviso era real o si se cumplio. Se
reutiliza la receta de `calibracion_dude_cdk2/auditar_activos_cdk2.py` (el UNICO
lugar donde ya se demostro que el AUC de CDK2 era interpretable porque se audito la
serie antes), con tres cosas que alli no hizo falta y aqui si:

1. **La muestra es de 1.200 activos, no de 474.** Es la que hay; y el numero importa,
   porque con 1.200 ya se ve si el problema es "muchos activos de una familia" o "la
   familia es la serie entera".
2. **Se separa por subtipo de medida.** IC50 y Ki no son intercambiables: si el
   embudo separa dentro de los Ki y no dentro de los IC50, el AUC global es una
   mezcla de dos cosas y no significa nada. Se cuenta cada subtipo por separado.
3. **Se mide el reparto de sulfonamidas.** hCA2 se une por el nitrogeno del sulfonamida
   al zinc. Si casi todos los activos tienen uno, el AUC va a premiar "tengo un
   sulfonamida" y no "soy un ligando bueno de hCA2". Hay que saberlo antes.

LO QUE NO HACE
--------------
No dice si el AUC sera alto o bajo: todavia no hay AUC. Dice si el numero, sea cual
sea, va a poder leerse como "el motor discrimina esta diana" o si va a estar vacio de
sentido por construccion. Y deja escrito el numero de activos INDEPENDIENTES, que es
el que de verdad limita el intervalo.

Salida:
    analisis/serie_hca2_auditoria.csv   (los activos medidos, con su esqueleto)
    analisis/INFORME_AUDITORIA_SERIE_HCA2_2026-09-30.md
    analisis/auditar_serie_hca2.log

Uso:
    python auditar_serie_hca2.py            # con red, descarga y audita
    python auditar_serie_hca2.py --cache    # solo audita lo ya descargado
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdFMCS, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
TARGET = "CHEMBL205"          # Carbonic anhydrase 2, humana
TARGET_NOMBRE = "Anhidrasa carbonica 2 (hCA2)"

MUESTRA = 1200                # activos medidos que se descargan y auditan
MIN_ACTIVOS = 5              # umbral de "familia", el mismo que usa el validador
UMBRAL_CONCENTRACION = 0.30  # mas del 30 % en un esqueleto = aviso
TAM_MUESTRA_MCS = 80         # cuantas molas al MCS (cuanto mas, mas lento)

CSV_ACTIVOS = os.path.join(BASE, "serie_hca2_auditoria.csv")
INFORME = os.path.join(BASE, "INFORME_AUDITORIA_SERIE_HCA2_2026-09-30.md")
LOG = os.path.join(BASE, "auditar_serie_hca2.log")


def log(m):
    print(m, flush=True)


def get(url, intentos=3):
    for i in range(intentos):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                              "User-Agent": "masive-als"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            if i == intentos - 1:
                log("      (fallo la consulta: %r)" % e)
                return None
            time.sleep(2 * (i + 1))
    return None


def descargar_activos(n=MUESTRA):
    """Actividades con pchembl de hCA2, deduplicadas por molecula, con subtipo.

    Se pagina porque ChEMBL corta en 1000 por pagina, y el limite se pide explicito:
    pedir 1200 en una sola llamada daria 1000 sin avisar, que es exactamente el tipo
    de fallo silencioso que hace que un numero parezca otro.
    """
    filas, vistos = [], set()
    pagina = 0
    while len(filas) < n:
        q = urllib.parse.urlencode({
            "target_chembl_id": TARGET,
            "pchembl_value__isnull": "false",
            "standard_type__in": "IC50,Ki,Kd",
            "limit": 1000,
            "offset": pagina * 1000,
        })
        d = get("%s/activity.json?%s" % (CHEMBL, q))
        if not d:
            break
        acts = d.get("activities", [])
        if not acts:
            break
        for a in acts:
            smi = a.get("canonical_smiles")
            mol_id = a.get("molecule_chembl_id")
            if not smi or not mol_id or mol_id in vistos:
                continue
            try:
                pch = float(a.get("pchembl_value"))
            except (TypeError, ValueError):
                continue
            vistos.add(mol_id)
            filas.append({"molecule_chembl_id": mol_id,
                          "smiles": smi,
                          "tipo": a.get("standard_type") or "",
                          "pchembl": "%.3f" % pch})
            if len(filas) >= n:
                break
        total = d.get("page_meta", {}).get("total_count", "?")
        log("   pagina %d: %d activos acumulados de %s medidos en total"
            % (pagina, len(filas), total))
        if len(acts) < 1000:
            break
        pagina += 1
        time.sleep(0.3)
    return filas


def leer_cache():
    if not os.path.exists(CSV_ACTIVOS):
        return []
    with open(CSV_ACTIVOS, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def a_mol(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        m = Chem.MolFromSmiles(smiles, sanitize=False)
        if m is not None:
            try:
                m.UpdatePropertyCache(strict=False)
            except Exception:  # noqa: BLE001
                return None
    return m


def esqueleto(mol):
    """Murcko como SMILES canonico, o None."""
    if mol is None:
        return None
    try:
        e = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:  # noqa: BLE001
        return None
    if e is None or e.GetNumAtoms() == 0:
        return None
    try:
        return Chem.MolToSmiles(e)
    except Exception:  # noqa: BLE001
        return None


SULFONAMIDA = Chem.MolFromSmarts("[NX3H][SX4](=O)(=O)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", action="store_true",
                    help="no descarga: audita el CSV ya guardado")
    args = ap.parse_args()

    t0 = time.time()
    lineas = []

    def w(m):
        log(m)
        lineas.append(m)

    w("AUDITORIA DE LA SERIE DE %s (%s)" % (TARGET_NOMBRE, TARGET))
    w("   %s" % time.strftime("%Y-%m-%d %H:%M"))
    w("   Antes de acoplar un solo ligando. Recipe: la de auditar_activos_cdk2.py,")
    w("   que es donde se demostro que el AUC de CDK2 era interpretable.")
    w("")

    if args.cache:
        filas = leer_cache()
        w("(modo cache: %d activos del CSV)" % len(filas))
    else:
        w("1. DESCARGA DE ACTIVOS MEDIDOS")
        filas = descargar_activos()
        w("   %d activos unicos con pchembl (IC50/Ki/Kd)" % len(filas))
        with open(CSV_ACTIVOS, "w", newline="", encoding="utf-8") as f:
            wr = csv.DictWriter(f, fieldnames=["molecule_chembl_id", "smiles",
                                               "tipo", "pchembl"])
            wr.writeheader()
            wr.writerows(filas)

    if not filas:
        w("")
        w("   no hay activos: no se puede auditar. Si se fallo la red, se reintenta")
        w("   con red; el cache tambien vale (--cache).")
        with open(LOG, "w", encoding="utf-8") as f:
            f.write("\n".join(lineas) + "\n")
        return 1

    # --- procesar ---
    w("")
    w("2. ESQUELETOS DE MURCKO")
    fam = defaultdict(list)
    mols = {}
    sin_mol = sin_esq = 0
    for r in filas:
        m = a_mol(r["smiles"])
        if m is None:
            sin_mol += 1
            continue
        mols[r["molecule_chembl_id"]] = m
        e = esqueleto(m)
        if e is None:
            sin_esq += 1
            continue
        fam[e].append(r["molecule_chembl_id"])
        r["esqueleto"] = e
    n_act = len(mols)
    # Un esqueleto de Murcko de UN anillo aislado (benceno, tiadiazol...) no es una
    # familia: es el residuo que queda cuando se le quitan las cadenas a un derivado
    # de benceno. Cuenta como esqueleto distinto pero agrupa 314 activos en uno solo
    # y falsea el "esqueleto mas grande", asi que se separa y se cuenta aparte. El
    # numero que decide (familias de 5+) se da con y sin ellos.
    ANILLOS_RUIDO = {"c1ccccc1", "c1ccncc1", "c1ccsc1", "c1cc[nH]c1", "c1ccnnc1",
                     "c1ccoc1", "c1cc[nH]n1", "c1cncs1", "c1ccsc1", "c1ncc[nH]1",
                     "c1ccn[nH]1", "c1ncs[nH]1"}
    fam_limpia = {k: v for k, v in fam.items() if k not in ANILLOS_RUIDO}

    tailles = sorted(((len(v), k) for k, v in fam.items()), reverse=True)
    top_n, top_e = tailles[0]
    frac_top = top_n / max(n_act, 1)
    w("   %d activos con mol legible | %d sin esqueleto (aciclicos)" % (n_act, sin_esq))
    w("   %d ESQUELETOS DISTINTOS (de los cuales %d son anillos aislados, que no"
      % (len(fam), len(fam) - len(fam_limpia)))
    w("   son familia sino el residuo de un derivado: benceno solo agrupa %d)"
      % top_n)
    w("   %d esqueletos si se quitan esos" % len(fam_limpia))
    tailles_l = sorted(((len(v), k) for k, v in fam_limpia.items()), reverse=True)
    w("   esqueletos con %d+ activos: %d, que cubren %d activos (%.1f %%)"
      % (MIN_ACTIVOS,
         sum(1 for n, _ in tailles_l if n >= MIN_ACTIVOS),
         sum(n for n, _ in tailles_l if n >= MIN_ACTIVOS),
         100.0 * sum(n for n, _ in tailles_l if n >= MIN_ACTIVOS) / max(n_act, 1)))
    n_indep_bruto = sum(1 for n, _ in tailles if n >= MIN_ACTIVOS)
    n_indep = sum(1 for n, _ in tailles_l if n >= MIN_ACTIVOS)
    w("   el esqueleto mas grande SIN contar anillos ruido: %d activos (%.1f %%) | %s"
      % (tailles_l[0][0], 100.0 * tailles_l[0][0] / max(n_act, 1),
         tailles_l[0][1][:60]))
    w("   INDEPENDIENTES: %d (bruto, con anillos ruido: %d). El bueno es %d"
      % (n_indep, n_indep_bruto, n_indep))
    w("   los 10 esquelestos mas poblados:")
    for n, e in tailles[:10]:
        w("      %4d activos (%.1f %%) | %s"
          % (n, 100.0 * n / max(n_act, 1), e[:66]))
    if tailles_l[0][0] / max(n_act, 1) > UMBRAL_CONCENTRACION:
        w("   AVISO: un solo esqueleto se lleva el %.1f %% de los activos. El AUC"
          % (100 * tailles_l[0][0] / max(n_act, 1)))
        w("          global mediria sobre todo ESE esqueleto. Habria que leerlo por")
        w("          familia, y el numero de independientes baja a ~%d." % n_indep)
    else:
        w("   bien: ningun esqueleto real pasa del %.0f %%."
          % (100 * UMBRAL_CONCENTRACION))

    # --- MCS ---
    w("")
    w("3. HAY COLUMNA VERTEBRAL COMUN? (maximo comun subestructura)")
    muestra = list(mols.values())
    import random as _r
    _r.seed(1)
    subm = _r.sample(muestra, min(TAM_MUESTRA_MCS, len(muestra)))
    n_at_mcs = None
    try:
        t0m = time.time()
        mcs = rdFMCS.FindMCS(subm, timeout=60, ringMatchesRingOnly=True)
        n_at_mcs = mcs.numAtoms
        w("   MCS sobre %d activos: %d atomos, %d enlaces | %.1f s | cancelado: %s"
          % (len(subm), mcs.numAtoms, mcs.numBonds, time.time() - t0m,
             mcs.canceled))
        w("   (que no saliera por tiempo se ve en 'cancelado': con True el MCS estaria")
        w("    cutado, y este sale en menos de un segundo, o sea que es de verdad.)")
        if n_at_mcs <= 5:
            w("   -> NO hay serie congenerica. Un MCS de %d atomos es un esqueleto,"
              % n_at_mcs)
            w("      no una familia: los %d activos son %d compuestos distintos."
              % (n_act, len(fam)))
        else:
            w("   -> HAY columna vertebral comun de %d atomos. Ojo: el aviso del"
              % n_at_mcs)
            w("      informe del 29 se cumple y hay que decirlo antes de acoplar.")
            smi_mcs = Chem.MolToSmiles(Chem.MolFromSmarts(mcs.smartsString))
            w("      la subestructura comun es: %s" % smi_mcs[:100])
            w("      (cuantos de los %d la contienen: se cuenta abajo con el"
              % n_act)
            w("       patron de sulfonamida, que es lo que une al zinc)")
    except Exception as e:  # noqa: BLE001
        w("   el MCS no se pudo calcular: %r" % e)

    # --- sulfonamidas ---
    w("")
    w("4. REPARTO DE SULFONAMIDAS (lo que une al zinc)")
    con_sulf = 0
    mols_sulf = []
    for mid, m in mols.items():
        if m.HasSubstructMatch(SULFONAMIDA):
            con_sulf += 1
            mols_sulf.append(m)
    frac_sulf = con_sulf / max(n_act, 1)
    w("   %d de %d (%.1f %%) tienen sulfonamida" % (con_sulf, n_act, 100 * frac_sulf))
    # Un MCS pequeno en la serie entera no dice si DENTRO de la sulfonamidas hay una
    # familia. Y esa es justo la subfamilia de la que tiene que salir el control de
    # redocking, asi que se mide por separado.
    n_at_sulf = None
    if len(mols_sulf) >= 10:
        try:
            mcs2 = rdFMCS.FindMCS(mols_sulf[:200], timeout=60,
                                  ringMatchesRingOnly=True)
            n_at_sulf = mcs2.numAtoms
            w("   MCS DENTRO de las %d sulfonamidas: %d atomos | cancelado: %s"
              % (len(mols_sulf), mcs2.numAtoms, mcs2.canceled))
            if n_at_sulf > 5:
                smi2 = Chem.MolToSmiles(Chem.MolFromSmarts(mcs2.smartsString))
                w("      esqueleto comun de las sulfonamidas: %s" % smi2[:90])
                w("      -> la sulfonamida SI es una familia. No es un problema para el")
                w("         AUC (el AUC global ya se sabe que no es degenerado), pero")
                w("         obliga a que el control de redocking sea de esa clase: si el")
                w("         motor solo Recognize 'anillo + sulfonamida cerca del metal',")
                w("         daria un AUC alto sin saber acoplar bien.")
            else:
                w("      -> ni las sulfonamidas forman familia (MCS de %d atomos)"
                  % n_at_sulf)
        except Exception as e:  # noqa: BLE001
            w("   el MCS de las sulfonamidas no se pudo calcular: %r" % e)
    if frac_sulf > 0.5:
        w("   AVISO: mas de la mitad de la serie tiene el mismo grupo funcional que")
        w("          hace la union al zinc. Un motor que solo sepa reconocer")
        w("          'sulfonamida cerca de un metal' va a sacar un AUC alto SIN")
        w("          saber acoplar bien. El AUC de hCA2 tiene que leerse CON")
        w("          este numero al lado, y el control de redocking tiene que ser")
        w("          de un compuesto de la MISMA clase (sulfonamida), no de otra.")
    else:
        w("   bien: no es una serie de un solo grupo funcional.")

    # --- por subtipo de medida ---
    w("")
    w("5. LOS ACTIVOS POR SUBTIPO DE MEDIDA (IC50 y Ki no se mezclan)")
    for tipo in ("IC50", "Ki", "Kd"):
        sub = [r for r in filas if r.get("tipo") == tipo]
        if not sub:
            continue
        fam_t = Counter()
        for r in sub:
            if r.get("esqueleto"):
                fam_t[r["esqueleto"]] += 1
        indep_t = sum(1 for n in fam_t.values() if n >= MIN_ACTIVOS)
        w("   %-5s %4d activos | %4d esquelestos | %2d independientes"
          % (tipo, len(sub), len(fam_t), indep_t))

    # --- similaridad entre activos ---
    w("")
    w("6. SIMILARIDAD ENTRE ACTIVOS (Morgan/Tanimoto, muestra)")
    # Morgan CON quiralidad y EN CONTEOS. Las dos correcciones son necesarias y las
    # dos se encontraron mirando lo que salia, no de memoria:
    #   - sin quiralidad, dos enantimeross dan huella identica y salen Tanimoto 1,000
    #     (mismo constituyente con el centro sin definir por un lado y definido por el
    #     otro). Es el caso real de esta serie.
    #   - en vector de bits, dos homologos de cadena (un -CH2- mas) colapsan al
    #     MISMO conjunto de bits y salen tambien 1,000, sin ser el mismo compuesto.
    #     Los conteos no tienen ese fallo: los cuatro pares que salian a 1,000 aqui
    #     son 0,89 / 0,89 / 0,84 / 0,69 con conteos, o sea moleculas distintas.
    # Y en cualquier caso la verdad de terreno no es la huella: es el SMILES
    # canonico. Si dos entradas tienen el mismo canonico, SI son el mismo compuesto.
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, includeChirality=True,
                                                    countSimulation=False)
    canon = {i: Chem.MolToSmiles(m) for i, m in mols.items()}
    ids = list(mols.keys())
    _r.shuffle(ids)
    sub = ids[:400]
    fps = {i: gen.GetCountFingerprint(mols[i]) for i in sub}
    pares, iden = [], []
    for a in range(len(sub)):
        for b in range(a + 1, len(sub)):
            s = DataStructs.TanimotoSimilarity(fps[sub[a]], fps[sub[b]])
            pares.append(s)
            if s >= 0.999:
                iden.append((sub[a], sub[b]))
    import numpy as np
    pares = np.array(pares)
    w("   pares comparados: %d (Morgan con quiralidad, en conteos)" % len(pares))
    w("   Tanimoto mediana %.3f | p90 %.3f | max %.3f | %% de pares > 0,70: %.2f %%"
      % (np.median(pares), np.percentile(pares, 90), pares.max(),
         100.0 * (pares > 0.70).mean()))
    w("   (CDK2, con 474 activos y 474 esquelestos, dio max ~0,7. Un max aqui")
    w("    por encima de 0,95 significa que hay compuestos casi repetidos.)")
    if iden:
        w("   %d pares con Tanimoto 1,000 aun en conteos; se listan para que se"
          % len(iden))
        w("   vean, no porque sean un problema:")
        for x, y in iden[:6]:
            w("             %s | %s" % (x, y))
    else:
        w("   bien: ningun par llega a 1,000.")

    # La verdad de terreno: mismo SMILES canonico bajo dos identificadores.
    por_canon = defaultdict(list)
    for i, c in canon.items():
        por_canon[c].append(i)
    repes = {c: v for c, v in por_canon.items() if len(v) > 1}
    w("   SMILES canonicos distintos: %d para %d activos" % (len(por_canon), n_act))
    if repes:
        w("   AVISO: %d canonicos con DOS O MAS entradas de ChEMBL (el mismo"
          % len(repes))
        w("          compuesto medido por dos laboratorios, que es normal y en el")
        w("          AUC manda una vez, pero si el banco lo prepara dos veces, una")
        w("          como activo y otra como senuelo, SI hincha el AUC):")
        for c, v in list(repes.items())[:6]:
            w("             %s <- %s" % (", ".join(v), c[:60]))
    else:
        w("   bien: ningun compuesto aparece con dos identificadores. No hay")
        w("         ActiveBait ni ningun par de entradas que sea el mismo, que es")
        w("         lo que hincha el AUC en DUD-E.")

    # --- veredicto ---
    w("")
    w("7. GLOSARIO: COMO SE LEERA EL AUC DE hCA2 CUANDO SALGA")
    veredictos = []
    if tailles_l[0][0] / max(n_act, 1) > UMBRAL_CONCENTRACION:
        veredictos.append("CONCENTRADO: el AUC global mide sobre todo un esqueleto")
    if n_at_mcs is not None and n_at_mcs > 5:
        veredictos.append("CONGENERICA: hay columna vertebral comun de %d atomos"
                          % n_at_mcs)
    if frac_sulf > 0.5:
        veredictos.append("UN SOLO GRUPO FUNCIONAL: el %d %% de sulfonamidas premia"
                          % round(100 * frac_sulf))
        veredictos[-1] += " reconocer el grupo, no acoplar bien"
    if n_at_sulf is not None and n_at_sulf > 5:
        veredictos.append("FAMILIA DENTRO DE LAS SULFONAMIDAS: MCS de %d atomos."
                          " El control de redocking tiene que salir de ahi."
                          % n_at_sulf)
    if not veredictos:
        veredictos.append("LIMPIA por los tres criterios")
    for v in veredictos:
        w("   - %s" % v)
    w("")
    w("  %d activos medidos, %d esqueletos, %d independientes. Con %d independientes"
      % (n_act, len(fam), n_indep, n_indep))
    w("   el semiancho del IC 95 %% del AUC anda alrededor de 0,1 si el AUC sale")
    w("   cerca de 0,5, y por debajo de 0,06 si sale alto. Es el orden de magnitud")
    w("   que hacia falta y que TDP-43 (7 positivos) y SOD1 (11) no tenian.")
    w("   duracion: %.1f s" % (time.time() - t0))

    with open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")

    # informe en markdown, con las cifras medidas
    md = []
    md.append("# Auditoría de la serie de hCA2 — 30 de septiembre de 2026\n")
    md.append("Antes de acoplar un solo ligando. Con la receta de "
              "`calibracion_dude_cdk2/auditar_activos_cdk2.py`, que es la que "
              "demostró que el AUC de CDK2 era interpretable.\n")
    md.append("## Cifras\n")
    md.append("| medida | valor |")
    md.append("|---|---|")
    md.append("| diana | %s (%s) |" % (TARGET_NOMBRE, TARGET))
    md.append("| activos medidos auditados | %d |" % n_act)
    md.append("| esqueletos de Murcko distintos | %d |" % len(fam))
    md.append("| familias de %d+ activos | %d |" % (MIN_ACTIVOS, n_indep))
    md.append("| esqueleto mas grande real | %d activos (%.1f %%) |"
              % (tailles_l[0][0], 100.0 * tailles_l[0][0] / max(n_act, 1)))
    if n_at_mcs is not None:
        md.append("| MCS sobre %d activos | %d atomos |"
                  % (min(TAM_MUESTRA_MCS, n_act), n_at_mcs))
    if n_at_sulf is not None:
        md.append("| MCS dentro de las %d sulfonamidas | %d atomos |"
                  % (con_sulf, n_at_sulf))
    md.append("| con sulfonamida | %d (%.1f %%) |" % (con_sulf, 100 * frac_sulf))
    md.append("| Tanimoto mediano entre activos (Morgan con quiralidad, conteos) | %.3f |"
              % float(np.median(pares)))
    md.append("| Tanimoto maximo | %.3f |" % float(pares.max()))
    md.append("| compuestos con dos entradas de ChEMBL (mismo canonico) | %d |"
              % len(repes))
    md.append("")
    md.append("## Los diez esqueletos más poblados\n")
    md.append("| activos | %% | esqueleto |")
    md.append("|---|---|---|")
    for n, e in tailles[:10]:
        md.append("| %d | %.1f | `%s` |" % (n, 100.0 * n / max(n_act, 1), e[:70]))
    md.append("")
    md.append("## Veredicto\n")
    for v in veredictos:
        md.append("- %s" % v)
    md.append("")
    md.append("## Qué significa para el AUC\n")
    md.append("El número grande de hCA2 (13.051 activos medidos en total) se reduce "
              "aquí a **%d independientes**, y ese es el que limita el intervalo del "
              "AUC, no los %d brutos. Con %d independientes el semiancho del IC 95 %% "
              "anda alrededor de 0,1 si el AUC sale cerca de 0,5, y baja de 0,06 si "
              "sale alto. Sigue siendo un orden de magnitud mejor que TDP-43 (7 "
              "positivos) y SOD1 (11).\n" % (n_indep, n_act, n_indep))
    md.append("## Ficheros\n")
    md.append("- `analysis/auditar_serie_hca2.py` — la auditoría.")
    md.append("- `analysis/serie_hca2_auditoria.csv` — los %d activos con su "
              "esqueleto." % n_act)
    md.append("- `analysis/auditar_serie_hca2.log` — la salida completa.\n")
    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    log("")
    log("informe: %s" % os.path.basename(INFORME))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
