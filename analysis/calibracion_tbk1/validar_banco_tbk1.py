#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Que dice el banco de TBK1 sobre el embudo del proyecto.

QUE MIDE
--------
Las mismas medidas que el proyecto ya tiene escritas, para que los numeros se
puedan poner al lado de los de TDP-43 sin traducir nada:

  * **AUC crudo** de los activos contra los decoys, con los empates a medias;
  * **AUC por atomo pesado**, que es la medida con la que el proyecto lee todas
    sus listas desde el 24 de septiembre, porque el crudo paga el tamano;
  * **EF1% y EF5%**, que es lo que de verdad se pregunta un investigador: cuantos
    activos entran en una lista corta;
  * **BEDROC (alpha 20)**, la medida que pesa mas lo que esta arriba. Se usa la
    implementacion de RDKit y se acompaña de su **linea de base al azar**, medida
    por permutacion: sin eso, un BEDROC de 0,4 no se sabe si es bueno o malo;
  * **el intervalo bootstrap del 95 %** del AUC (2000 remuestreos), porque la
    regla de decision del proyecto solo declara un bloque cuando el limite
    inferior pasa de 0,5;
  * **el reparto por quimiotipo**, contando cuantos esqueletos distintos de
    activos baten al azar. Es obligatorio mirarlo aqui: la quimica de TBK1 es muy
    congenerica y un AUC global inflado por una sola serie no vale nada.

Uso:
    python validar_banco_tbk1.py
    python validar_banco_tbk1.py --min-quimiotipo 5

Salida en `analysis/calibracion_tbk1/`:
    validar_tbk1.csv      ligando, papel, energia, atomos pesados y las dos
                          puntuaciones, para rehacer cualquier cuenta
    INFORME_CALIBRACION_TBK1.md   el informe
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.ML.Scoring import Scoring

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(BASE))
LIGANDS = os.path.join(BASE, "ligands")
RESULTADOS = os.path.join(BASE, "resultados_tbk1.csv")
ACTIVOS = os.path.join(BASE, "compuestos_tbk1.csv")
DECOYS = os.path.join(BASE, "decoys_tbk1.csv")
CSV_SALIDA = os.path.join(BASE, "validar_tbk1.csv")
INFORME = os.path.join(BASE, "INFORME_CALIBRACION_TBK1.md")

sys.path.insert(0, RAIZ)
import preparar_ligando as PL  # noqa: E402

B = 2000          # remuestreos del bootstrap
SEMILLA = 20260927
ALPHA = 20.0      # el alpha de BEDROC de Truchon y Bayly
N_AZAR = 200      # permutaciones para la linea de base de BEDROC


def auc_np(act, dec):
    """AUC con empates a medias, sin matriz cartesiana de activos x decoys."""
    act = np.asarray(act, dtype=float)
    dec = np.asarray(dec, dtype=float)
    # Las puntuaciones se definen como -energia/atomos: mas alto es mejor.
    # Para cada activo, searchsorted cuenta decoys menores y empatados. Esto
    # conserva exactamente P(activo > decoy) + 0.5*P(empate), usando O(n) memoria.
    ordenados = np.sort(dec)
    menores = np.searchsorted(ordenados, act, side="left")
    menores_o_iguales = np.searchsorted(ordenados, act, side="right")
    ganancias = menores + 0.5 * (menores_o_iguales - menores)
    return float(ganancias.sum()) / (act.size * dec.size)


def auc_rdkit(act, dec):
    """AUC de RDKit; con empates entre clases su orden es dependiente de lista."""
    lista = [(float(v), 1) for v in act] + [(float(v), 0) for v in dec]
    # CalcAUC espera los positivos primero en la lista descendente; el mismo
    # sentido que auc_np(): una puntuacion mayor debe favorecer al activo.
    lista.sort(key=lambda t: -t[0])
    return float(Scoring.CalcAUC(lista, 1))


def bedroc(act, dec, alpha=ALPHA):
    lista = [(float(v), 1) for v in act] + [(float(v), 0) for v in dec]
    lista.sort(key=lambda t: -t[0])
    return float(Scoring.CalcBEDROC(lista, 1, alpha))


def factor_enriquecimiento(act, dec, pct):
    """EF a pct; comparte fraccionalmente el puesto de corte si hay empates."""
    a = np.asarray(act, dtype=float)
    d = np.asarray(dec, dtype=float)
    n_a, total = a.size, a.size + d.size
    corte = max(1, int(round(total * pct / 100.0)))
    puntuaciones = np.concatenate((a, d))
    umbral = np.sort(puntuaciones)[-corte]
    activos_encima = int(np.count_nonzero(a > umbral))
    total_encima = int(np.count_nonzero(puntuaciones > umbral))
    empatados = int(np.count_nonzero(puntuaciones == umbral))
    activos_empatados = int(np.count_nonzero(a == umbral))
    cupos_empate = corte - total_encima
    n_en_corte = activos_encima + cupos_empate * activos_empatados / empatados
    esperado = corte * (n_a / total)
    return (n_en_corte / esperado) if esperado else 0.0


def bootstrap_auc(act, dec, remuestreos=B, semilla=SEMILLA):
    rng = np.random.default_rng(semilla)
    a = np.asarray(act, dtype=float)
    d = np.asarray(dec, dtype=float)
    vals = []
    for _ in range(remuestreos):
        av = a[rng.integers(0, a.size, a.size)]
        dv = d[rng.integers(0, d.size, d.size)]
        vals.append(auc_np(av, dv))
    return float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))


def esqueleto(smiles):
    m = Chem.MolFromSmiles(smiles or "")
    if m is None:
        return None
    try:
        return MurckoScaffold.MurckoScaffoldSmiles(mol=m)
    except Exception:  # noqa: BLE001
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-quimiotipo", type=int, default=5,
                    help="activos minimos para que un esqueleto cuente")
    args = ap.parse_args()

    if not os.path.exists(RESULTADOS):
        raise SystemExit("falta %s: primero hay que acoplar el banco"
                         % os.path.relpath(RESULTADOS, RAIZ))

    # --- los datos ---
    act = {r["molecule_chembl_id"]: r for r in
           csv.DictReader(open(ACTIVOS, encoding="utf-8"))}
    dec = {r["id"]: r for r in csv.DictReader(open(DECOYS, encoding="utf-8"))}
    energia = {}
    for r in csv.DictReader(open(RESULTADOS, encoding="utf-8")):
        energia[r["ligand"]] = float(r["energy"])

    filas = []
    sin_pose_act = sin_pose_dec = 0
    for ident, grupo, etiqueta in ((act, "activos", 1), (dec, "decoys", 0)):
        for k, v in ident.items():
            nombre = ("ACT_" if etiqueta else "DEC_") + k
            if nombre not in energia:
                if etiqueta:
                    sin_pose_act += 1
                else:
                    sin_pose_dec += 1
                continue
            n, glue = PL.contar_fichero(os.path.join(LIGANDS, nombre + ".pdbqt"))
            if not n:
                continue
            e = energia[nombre]
            filas.append({"ligando": nombre, "papel": grupo, "chembl": k,
                          "smiles": v.get("smiles", ""),
                          "energia": e, "pesados": n,
                          "crudo": -e, "por_atomo": -e / n,
                          "pchembl": v.get("pchembl", "")})
    print("con pose: %d activos (%d sin pose) y %d decoys (%d sin pose)"
          % (sum(1 for f in filas if f["papel"] == "activos"), sin_pose_act,
             sum(1 for f in filas if f["papel"] == "decoys"), sin_pose_dec))

    A = [f for f in filas if f["papel"] == "activos"]
    D = [f for f in filas if f["papel"] == "decoys"]
    ratio = len(D) / max(len(A), 1)

    lineas = []
    lineas.append("INFORME DE CALIBRACION DE TBK1 — %s" % time.strftime("%Y-%m-%d"))
    lineas.append("")
    lineas.append("banco: %d activos con pose y %d decoys con pose (%.1f por activo)"
                  % (len(A), len(D), ratio))
    lineas.append("activos sin pose: %d | decoys sin pose: %d"
                  % (sin_pose_act, sin_pose_dec))

    resultados = {}
    for nombre, campo in (("crudo", "crudo"), ("por atomo pesado", "por_atomo")):
        a = [f[campo] for f in A]
        d = [f[campo] for f in D]
        auc = auc_np(a, d)
        # RDKit CalcAUC no asigna automáticamente medio acierto a empates entre
        # clases; compara solo cuando las puntuaciones no se solapan.
        if set(a).isdisjoint(d):
            rd = auc_rdkit(a, d)
            if abs(auc - rd) > 1e-6:
                lineas.append("AVISO: el AUC propio (%.6f) y el de RDKit (%.6f) no"
                              " coinciden" % (auc, rd))
        else:
            lineas.append("AUC RDKit no usado como contraste: hay puntuaciones"
                          " empatadas entre activos y decoys; auc_np las cuenta a medias.")
        lo, hi = bootstrap_auc(a, d)
        b = bedroc(a, d)
        ef1 = factor_enriquecimiento(a, d, 1.0)
        ef5 = factor_enriquecimiento(a, d, 5.0)
        resultados[nombre] = {"auc": auc, "lo": lo, "hi": hi, "bedroc": b,
                              "ef1": ef1, "ef5": ef5}
        lineas.append("")
        lineas.append("* %s:" % nombre)
        lineas.append("  AUC %.3f (bootstrap 95 %%: %.3f - %.3f) -> %s"
                      % (auc, lo, hi, "PASA" if lo > 0.5 else "NO PASA"))
        lineas.append("  EF1 %% = %.2f | EF5 %% = %.2f | BEDROC(20) = %.3f"
                      % (ef1, ef5, b))

    # --- linea de base de BEDROC, medida y no supuesta ---
    rng = np.random.default_rng(SEMILLA)
    a_at = np.array([f["por_atomo"] for f in A])
    d_at = np.array([f["por_atomo"] for f in D])
    azar = []
    for _ in range(N_AZAR):
        cat = np.concatenate([a_at, d_at])
        rng.shuffle(cat)
        azar.append(bedroc(cat[:len(a_at)], cat[len(a_at):]))
    base = float(np.mean(azar))
    lineas.append("")
    lineas.append("BEDROC al azar con este banco (%d permutaciones): %.3f"
                  % (N_AZAR, base))
    lineas.append("   el BEDROC de arriba solo vale si esta claramente por encima"
                  " de esa base, no por encima de cero.")

    # --- quimiotipos ---
    grupos = {}
    for f in A:
        esc = esqueleto(f["smiles"])
        if esc:
            grupos.setdefault(esc, []).append(f["por_atomo"])
    validos = [(e, v) for e, v in grupos.items() if len(v) >= args.min_quimiotipo]
    pasan = 0
    detalle = []
    for esc, vals in sorted(validos, key=lambda t: -len(t[1])):
        auc_esc = auc_np(vals, [f["por_atomo"] for f in D])
        if auc_esc > 0.5:
            pasan += 1
        detalle.append((len(vals), esc, auc_esc))
    lineas.append("")
    lineas.append("quimiotipos (esqueletos de Murcko con >= %d activos): %d"
                  % (args.min_quimiotipo, len(validos)))
    lineas.append("   de ellos, con AUC por atomo > 0,5: %d" % pasan)
    lineas.append("   es la cuenta que decide: la quimica de TBK1 es muy"
                  " congenerica y un AUC global alto puede venir de una sola serie.")
    for n, esc, auc_esc in detalle[:12]:
        lineas.append("     %2d activos  AUC %.3f  %s"
                      % (n, auc_esc, esc[:70]))

    lineas.append("")
    lineas.append("LO QUE HAY QUE LEER, DECIDIDO ANTES DE VERLO (plan de calibracion):")
    lineas.append("   * si TBK1 sale bien (activos arriba, EF1 % mayor que 1 y varios")
    lineas.append("     quimiotipos pasando), el embudo funciona y el problema de")
    lineas.append("     TDP-43 es su bolsillo, no la herramienta;")
    lineas.append("   * si sale mal, lo que hay que rehacer es el motor (Vinardo, o el")
    lineas.append("     rescoring por interacciones) contra este banco;")
    lineas.append("   * si sale a medias, el acoplamiento solo vale como prefiltro.")

    texto = "\n".join(lineas)
    print("\n" + texto)
    with open(INFORME, "w", encoding="utf-8") as f:
        f.write(texto + "\n")

    with open(CSV_SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ligando", "papel", "chembl", "energia",
                                          "pesados", "crudo", "por_atomo", "pchembl"])
        w.writeheader()
        w.writerows([{k: v for k, v in r.items() if k in w.fieldnames} for r in filas])
    print()
    print("CSV: %s" % os.path.relpath(CSV_SALIDA, RAIZ))
    print("informe: %s" % os.path.relpath(INFORME, RAIZ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
