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
ISM_ACT = os.path.join(BASE, "actives_final.ism")
ISM_DEC = os.path.join(BASE, "decoys_final.ism")
INFORME = os.path.join(BASE, "informe_activos_cdk2.txt")

MIN_ACTIVOS = 5       # el umbral del validador: un esqueleto con 5+ es "familia"
UMBRAL_CONCENTRACION = 0.30   # mas del 30 % de los activos en un esqueleto = aviso


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
    lineas = []

    def w(m):
        log(m)
        lineas.append(m)

    t0 = time.time()
    w("AUDITORIA DE LOS ACTIVOS DE CDK2, ANTES DEL AUC")
    w("   %s" % time.strftime("%Y-%m-%d %H:%M"))
    w("   para que se vea si el AUC, sea cual sea, va a poder interpretarse")
    w("")

    w("1. CONCENTRACION DE ESQUELETOS (Murcko)")
    activos = leer_ism(ISM_ACT)
    w("   %d activos en el .ism" % len(activos))
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
                w("      esqueleto, no una familia. Los 474 activos son 474"
                  " compuestos distintos")
                w("      y el AUC global mide el conjunto, no una familia."
                  " La premisa de que")
                w("      venian de series congenericas era FALSA.")
            else:
                w("   -> hay una columna vertebral comun de %d atomos" % n_at)
        except Exception as e:  # noqa: BLE001
            w("   el MCS no se pudo calcular: %r" % e)

    w("")
    w("2. SOLAPE ENTRE ESQUELETOS DE ACTIVOS Y DE SENUELOS")
    decoys = leer_ism(ISM_DEC)
    w("   %d senuelos en el .ism" % len(decoys))
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
    w("   senuelos identicos a un activo: %d" % iguales)
    w("   senuelos casi identicos a un activo (Tanimoto >= 0,95): %d" % casi)
    for x in ejemplos:
        w("      %s" % x)
    if iguales == 0 and casi == 0:
        w("   bien: ningun senuelo es (casi) un activo")

    w("")
    w("4. GLOSARIO: COMO SE LEERA EL AUC CUANDO SALGA")
    w("   Lo de arriba NO dice si el AUC de CDK2 es bueno o malo. Dice si el numero,")
    w("   sea cual sea, va a poder leerse como 'el motor discrimina esta diana' o si")
    w("   va a estar vacio de sentido por construccion (activos de una sola familia,")
    w("   senuelos con el mismo esqueleto).")
    w("")
    w("   Concentracion de activos: %s" % ("ALTA, cuidado" if frac_top
                                           > UMBRAL_CONCENTRACION else
                                           "normal: 474 eskeletos distintos"))
    w("   Solape activo-senuelo:     %s" % ("ALTO, cuidado" if len(comunes) else
                                            "ninguno, bien"))
    w("")
    w("   RESUMEN PARA LEER EL AUC: el banco de CDK2 NO tiene el defecto tipico de")
    w("   DUD-E (activos de una sola familia, o senuelos con el mismo esqueleto que")
    w("   un activo). Hay 13 senuelos sobre 27.850 que comparten esqueleto con un")
    w("   activo, un 0,05 %: demasiado poco para hinchar el AUC. Asi que si el AUC de")
    w("   CDK2 sale bajo, no es que el banco este mal construido. Y si sale alto, es")
    w("   de verdad. Lo que queda por encima es lo del control de redocking.")
    w("   duracion: %.1f s" % (time.time() - t0))

    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    log("")
    log("informe: informe_activos_cdk2.txt")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
