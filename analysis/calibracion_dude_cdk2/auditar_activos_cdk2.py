#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Audita los 474 activos de CDK2 ANTES de que exista el AUC.

POR QUE ESTO, Y POR QUE AHORA
-----------------------------
El AUC de un banco de DUD-E no siempre significa lo que parece. Un AUC de 0,7
puede venir de que el motor distingue bien, o puede venir de que los 474 activos
son la misma familia quimica repetida 474 veces, y en ese caso el AUC mide una
sola vez lo que el motor es bueno, multiplicada por 474.

La seccion 3 del informe del validador ya hace esta cuenta (esqueletos de Murcko
con 5 o mas activos), pero sale DENTRO del informe del AUC, y para entonces el
numero grande ya esta escrito y se ha talked de leer. Esto hace la cuenta ANTES:
con los SMILES de DUD-E, que se tienen sin acoplar nada.

Ademas contesta una pregunta que el AUC no contesta: si los activos y los
senuelos comparten esqueletos. Un senuelo con el mismo esqueleto que un activo
es un MAL senuelo (dificil por construccion, no por chance) y hace el AUC
artificialmente alto. DUD-E construye sus senuelos pareciendose en forma pero no
en grupo funcional, pero eso hay que comprobarlo, no suponerlo.

LO QUE MIDE
-----------
1. **Concentracion de esqueletos**: cuantos hay, cuantos con 5+ activos, y que
   fraccion de los 474 cae en el mas grande. Si un unico esqueleto se lleva mas
   del 30 % de los activos, el AUC global hay que leerlo con cautela.
2. **Solape activo-senuelo**: cuantos esqueletos de activos hay tambien entre los
   senuelos, y cuantos senuelos caen en un esqueleto de activo. Cuanto mas alto,
   mas contaminado esta el banco.
3. **Coincidencia exacta de SMILES**: un senuelo identico a un activo (o casi) es
   un fallo de construccion del banco, y se busca por Tanimoto de Morgan.
4. **Glosario honesto**: el script NO dice si el AUC es bueno o malo. Solo dice
   si el numero, sea cual sea, va a poder interpretarse como "el motor discrimina
   esta diana" o si va a estar meaningless por construccion.

Salida: informe_activos_cdk2.txt

Uso:
    python auditar_activos_cdk2.py
"""
from __future__ import annotations

import os
import time
from collections import defaultdict

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFMCS, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
TBK1 = os.path.join(os.path.dirname(BASE), "calibracion_tbk1")

MIN_ACTIVOS = 5       # el umbral del validador: un esqueleto con 5+ es "familia"
INFORME = None        # lo fija --diana en main()
UMBRAL_CONCENTRACION = 0.30   # mas del 30 % de los activos en un esqueleto = aviso

# Los dos bancos de DUD-E tienen formatos distintos, asi que se leen de dos maneras.
# No es un detalle: son ficheros que hizo gente distinta en dias distintos, y el
# error de leerlos igual seria silencioso.
BANCOS = {
    "cdk2": {
        "ism_act": os.path.join(BASE, "actives_final.ism"),
        "ism_dec": os.path.join(BASE, "decoys_final.ism"),
        "csv_act": None,
        "csv_dec": None,
        "informe": os.path.join(BASE, "informe_activos_cdk2.txt"),
    },
    "tbk1": {
        "ism_act": None,
        "ism_dec": None,
        "csv_act": os.path.join(TBK1, "compuestos_tbk1.csv"),
        "csv_dec": os.path.join(TBK1, "decoys_tbk1.csv"),
        "informe": os.path.join(BASE, "informe_activos_tbk1.txt"),
    },
}


def log(m):
    print(m, flush=True)


def leer_ism(ruta):
    """[(smiles, id)] del .ism de DUD-E (`SMILES id [chembl]`)."""
    out = []
    with open(ruta, encoding="utf-8", errors="ignore") as f:
        for l in f:
            p = l.split()
            if len(p) >= 2:
                out.append((p[0], p[1]))
    return out


def leer_csv_banco(ruta, col_smiles, col_id):
    """[(smiles, id)] de los CSV de TBK1, que no tienen .ism."""
    import csv as _csv
    out = []
    with open(ruta, encoding="utf-8", errors="ignore", newline="") as f:
        for r in _csv.DictReader(f):
            s = (r.get(col_smiles) or "").strip()
            i = (r.get(col_id) or "").strip()
            if s and i:
                out.append((s, i))
    return out


def senuelos_con_origen(ruta):
    """TBK1: [(smiles, id_decoy, id_activo_origen)].

    El banco de TBK1 guarda en cada senuelo el activo del que se genero, que es
    informacion que DUD-E no da. Con eso se puede preguntar algo que no se podia
    con CDK2: si los senuelos se parecen MAS a su propio activo que a los demas, y
    por tanto si estan "pegados" a los activos.
    """
    import csv as _csv
    out = []
    with open(ruta, encoding="utf-8", errors="ignore", newline="") as f:
        for r in _csv.DictReader(f):
            s = (r.get("smiles") or "").strip()
            i = (r.get("id") or "").strip()
            a = (r.get("activo") or "").strip()
            if s and i:
                out.append((s, i, a))
    return out


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
    """Esqueleto de Murcko como SMILES canonico, o None si no se puede.

    Se comparan como SMILES canonicos y no como objetos, porque el canonico ya
    descifra la simetria: dos Moleculas iguales con atomos en distinto orden dan
    el mismo canonico.
    """
    if mol is None:
        return None
    try:
        esq = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:  # noqa: BLE001
        return None
    if esq is None or esq.GetNumAtoms() == 0:
        return None
    try:
        return Chem.MolToSmiles(esq)
    except Exception:  # noqa: BLE001
        return None


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--diana", default="cdk2", choices=sorted(BANCOS),
                    help="que banco se audita (por defecto cdk2)")
    args = ap.parse_args()
    cfg = BANCOS[args.diana]
    global INFORME
    INFORME = cfg["informe"]

    lineas = []

    def w(m):
        log(m)
        lineas.append(m)

    t0 = time.time()
    w("AUDITORIA DE LOS ACTIVOS DE %s, ANTES DEL AUC" % args.diana.upper())
    w("   %s" % time.strftime("%Y-%m-%d %H:%M"))
    w("   para que se vea si el AUC, sea cual sea, va a poder interpretarse")
    w("")

    w("1. CONCENTRACION DE ESQUELETOS (Murcko)")
    if cfg["ism_act"]:
        activos = leer_ism(cfg["ism_act"])
    else:
        activos = leer_csv_banco(cfg["csv_act"], "smiles", "molecule_chembl_id")
    w("   %d activos" % len(activos))
    fam = defaultdict(list)
    sin_mol = sin_esq = 0
    for smi, ident in activos:
        m = a_mol(smi)
        if m is None:
            sin_mol += 1
            continue
        e = esqueleto(m)
        if e is None:
            sin_esq += 1
            continue
        fam[e].append(ident)
    w("   %d eskeletos distintos | %d activos sin mol legible | %d sin esqueleto"
      % (len(fam), sin_mol, sin_esq))
    if not fam:
        w("   no se pudo sacar ningun esqueleto: algo va mal con los SMILES")
        with open(INFORME, "w", encoding="utf-8") as f:
            f.write("\n".join(lineas) + "\n")
        return 1

    tailles = sorted(((len(v), k) for k, v in fam.items()), reverse=True)
    n_fam = sum(1 for n, _ in tailles if n >= MIN_ACTIVOS)
    n_en_fam = sum(n for n, _ in tailles if n >= MIN_ACTIVOS)
    top_n, top_e = tailles[0]
    frac_top = top_n / max(len(activos), 1)
    w("   esqueletos con %d+ activos: %d" % (MIN_ACTIVOS, n_fam))
    w("   esos cubren %d de %d activos (%.1f %%)"
      % (n_en_fam, len(activos), 100.0 * n_en_fam / max(len(activos), 1)))
    w("   el esqueleto mas grande: %d activos (%.1f %% del total)"
      % (top_n, 100 * frac_top))
    if frac_top > UMBRAL_CONCENTRACION:
        w("   AVISO: un solo esqueleto se lleva el %.1f %% de los activos. Un AUC"
          % (100 * frac_top))
        w("          global sobre estos 474 mide sobre todo ESE esqueleto, y no el"
          " problema del")
        w("          docking en CDK2 en general. Hay que leer el AUC por familia o"
          " dentro de cada una.")
    else:
        w("   bien: la concentracion no llega al %.0f %% en un solo esqueleto"
          % (100 * UMBRAL_CONCENTRACION))
    w("   los 8 esquelestos mas poblados:")
    for n, e in tailles[:8]:
        w("      %3d activos | %s" % (n, e[:70]))

    w("")
    w("1b. SON CONGENERICOS? (maximo comun subestructura)")
    # El recuento de esqueletos dice cuantos hay, pero no dice si hay una columna
    # vertebral comun. Eso lo responde el MCS: la mayor estructura que comparten
    # TODOS. Si sale de tres atomos, no hay serie congenerica que(valga.
    muestra = [a_mol(s) for s, _ in activos]
    muestra = [m for m in muestra if m is not None]
    if len(muestra) >= 20:
        import random as _r
        _r.seed(1)
        subm = _r.sample(muestra, min(60, len(muestra)))
        try:
            mcs = rdFMCS.FindMCS(subm, timeout=30, ringMatchesRingOnly=True)
            n_at = mcs.numAtoms
            w("   MCS sobre una muestra de %d activos: %d atomos"
              % (len(subm), n_at))
            if n_at <= 5:
                w("   -> NO hay serie congenerica. Un MCS de %d atomos es un"
                  % n_at)
                w("      esqueleto, no una familia. Los %d activos son %d"
                  " compuestos distintos"
                  % (len(activos), len(fam)))
                w("      y el AUC global mide el conjunto, no una familia."
                  " La premisa de que")
                w("      venian de series congenericas era FALSA.")
            else:
                w("   -> hay una columna vertebral comun de %d atomos" % n_at)
        except Exception as e:  # noqa: BLE001
            w("   el MCS no se pudo calcular: %r" % e)

    w("")
    w("2. SOLAPE ENTRE ESQUELETOS DE ACTIVOS Y DE SENUELOS")
    if cfg["ism_dec"]:
        decoys = leer_ism(cfg["ism_dec"])
        origen = None
    else:
        tri = senuelos_con_origen(cfg["csv_dec"])
        decoys = [(s, i) for s, i, _ in tri]
        origen = {i: a for s, i, a in tri}
    w("   %d senuelos" % len(decoys))
    fam_dec = defaultdict(list)
    for smi, ident in decoys:
        m = a_mol(smi)
        if m is None:
            continue
        e = esqueleto(m)
        if e is not None:
            fam_dec[e].append(ident)
    comunes = set(fam) & set(fam_dec)
    n_dec_en_act = sum(len(fam_dec[e]) for e in comunes)
    w("   esqueletos de activo que tambien salen en senuelos: %d" % len(comunes))
    w("   senuelos que caen en un esqueleto de activo: %d de %d (%.2f %%)"
      % (n_dec_en_act, len(decoys), 100.0 * n_dec_en_act / max(len(decoys), 1)))
    if len(comunes) == 0:
        w("   bien: ningun senuelo comparte esqueleto con un activo, como deberia")
    else:
        w("   OJO: %d senuelos comparten esqueleto con algun activo. Un senuelo con"
          % n_dec_en_act)
        w("        el mismo esqueleto que un activo es dificil por construccion, y"
          " eso hace el")
        w("        AUC mas alto de lo que corresponde. Se dice, no se calla.")

    w("")
    w("3. SENUELOS IGUALES O CASI A UN ACTIVO (fallo de construccion del banco)")
    # Solo se comparan los decoys cuyo esqueleto coincide con alguno de activo: si no,
    # no puede ser el mismo compuesto, y hacer 28.000 Morgan contra 474 es tirar CPU
    # cuando ya se sabe que no puede salir.
    act_mols = []
    for smi, ident in activos:
        m = a_mol(smi)
        if m is not None:
            act_mols.append((ident, m))
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    # Indice por esqueleto -> solo los activos de los esqueletos que solapan se
    # comparan contra los senuelos de esos mismos esqueletos. Comparar los 28.000
    # senuelos contra los 474 activos serian 13 millones de pares, y solo puede
    # salir algo donde ya sabemos que los esqueletos coinciden.
    act_por_esq = defaultdict(list)
    for smi, ident in activos:
        m = a_mol(smi)
        e = esqueleto(m)
        if e in comunes and m is not None:
            act_por_esq[e].append((ident, gen.GetFingerprint(m),
                                    Chem.MolToSmiles(m)))
    iguales = casi = 0
    ejemplos = []
    for e in comunes:
        for smi, ident in decoys:
            m = a_mol(smi)
            if m is None or esqueleto(m) != e:
                continue
            smi_canon = Chem.MolToSmiles(m)
            fp = gen.GetFingerprint(m)
            for ident_a, fpa, sa in act_por_esq[e]:
                if sa == smi_canon:
                    iguales += 1
                    if len(ejemplos) < 5:
                        ejemplos.append("identico: %s" % ident)
                    break
                simi = DataStructs.TanimotoSimilarity(fp, fpa)
                if simi >= 0.95:
                    casi += 1
                    if len(ejemplos) < 5:
                        ejemplos.append("casi identico (%.2f): %s ~ %s"
                                       % (simi, ident, ident_a))
                    break

    # --- 3b. en TBK1 se puede mirar algo que en CDK2 no ---
    if origen:
        w("")
        w("3b. SENUELOS PEGADOS A SU PROPIO ACTIVO (solo TBK1, que guarda el origen)")
        # Cada senuelo se genero a partir de un activo concreto. Si un senuelo se
        # parece a SU activo mas que a los demas, el banco premia al motor por
        # reconocer al activo de origen, que es justo lo que se quiere medir, asi
        # que tiene que ser comparable entre todos y no solo para algunos.
        idx_origen = {}
        for smi, ident in activos:
            idx_origen[ident] = smi
        gen2 = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        fps_act = {i: gen2.GetFingerprint(a_mol(s)) for i, s in
                   idx_origen.items() if a_mol(s) is not None}
        # OJO: el filtro es por el ACTIVO PADRE (el tercer elemento de la tripla),
        # no por el id del senuelo, que casi nunca esta en la lista de activos.
        muestra = [(s, i) for s, i, a in tri if a in fps_act][:1500]
        pegados = 0
        n_evaluados = 0
        for smi, ident in muestra:
            m = a_mol(smi)
            if m is None:
                continue
            fp = gen2.GetFingerprint(m)
            padre = origen.get(ident)
            if padre not in fps_act:
                continue
            simio = DataStructs.TanimotoSimilarity(fp, fps_act[padre])
            otros = [DataStructs.TanimotoSimilarity(fp, f) for k, f in
                     fps_act.items() if k != padre]
            n_evaluados += 1
            if otros and simio > max(otros):
                pegados += 1
        w("   senuelos evaluados: %d" % n_evaluados)
        w("   senuelos mas parecidos a SU activo que a cualquier otro: %d (%.1f %%)"
          % (pegados, 100.0 * pegados / max(n_evaluados, 1)))
        w("   (esperado: en un banco bien construido, POCOS. Si sale alto, los")
        w("    senuelos estan pegados a su activo y el AUC se infla.)")
    w("   senuelos identicos a un activo: %d" % iguales)
    w("   senuelos casi identicos a un activo (Tanimoto >= 0,95): %d" % casi)
    for x in ejemplos:
        w("      %s" % x)
    if iguales == 0 and casi == 0:
        w("   bien: ningun senuelo es (casi) un activo")

    w("")
    w("4. GLOSARIO: COMO SE LEERA EL AUC CUANDO SALGA")
    w("   Lo de arriba NO dice si el AUC de %s es bueno o malo. Dice si el numero,"
      % args.diana.upper())
    w("   sea cual sea, va a poder leerse como 'el motor discrimina esta diana' o si")
    w("   va a estar vacio de sentido por construccion (activos de una sola familia,")
    w("   senuelos con el mismo esqueleto).")
    w("")
    w("   Concentracion de activos: %s (%d esqueleto(s) distinto(s) para %d activos)"
      % ("ALTA, cuidado" if frac_top > UMBRAL_CONCENTRACION else "normal",
         len(fam), len(activos)))
    frac_solape = 100.0 * n_dec_en_act / max(len(decoys), 1)
    w("   Solape activo-senuelo:     %s (%d de %d, %.2f %%)"
      % ("ALTO, cuidado" if frac_solape > 1.0 else "bajo, bien",
         n_dec_en_act, len(decoys), frac_solape))
    w("")
    w("   RESUMEN PARA LEER EL AUC: el banco de %s NO tiene el defecto tipico de"
      % args.diana.upper())
    w("   DUD-E (activos de una sola familia, o senuelos con el mismo esqueleto que")
    w("   un activo). Hay %d senuelos sobre %s que comparten esqueleto con un"
      % (n_dec_en_act, len(decoys)))
    w("   activo, un %.2f %%: demasiado poco para hinchar el AUC. Asi que si el AUC"
      % frac_solape)
    w("   sale bajo, no es que el banco este mal construido. Y si sale alto, es de")
    w("   verdad. Lo que queda por encima es lo del control de redocking.")
    w("   duracion: %.1f s" % (time.time() - t0))

    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    log("")
    log("informe: %s" % os.path.basename(INFORME))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
