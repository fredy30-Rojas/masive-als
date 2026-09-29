#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Que dice el banco CDK2 de DUD-E sobre el embudo del proyecto.

POR QUE CDK2 Y NO TBK1
----------------------
La calibracion de TBK1 dejo un techo claro: AUC crudo 0,602, por atomo pesado
0,617, 0,637 con la media de las poses, y el MM-GBSA en 0,498. Con esos numeros
solos NO se sabe si el techo es del embudo o de la diana: para distinguirlo hace
falta una diana donde la respuesta sea publica y el banco lo hayan hecho otros.

Eso es CDK2 de DUD-E: 474 activos y 27.850 senuelos emparejados por fisicoquimica
y por forma, publicados con la estructura 1h00, y con un numero de referencia
encima de la mesa —el DOCK 3.6 de los propios autores da AUC 0,791 (Mysinger y
col., "Directory of Useful Decoys, Enhanced", J. Med. Chem. 2012).

ENTONCES LA LECTURA ESTA DECIDIDA DE ANTEMANO
--------------------------------------------
  * CDK2 en ~0,79 -> el embudo esta bien y el techo de 0,62 es de TBK1;
  * CDK2 en ~0,60 -> el techo es del MOTOR (Vina-GPU con search_depth 20), no de
    la diana, y lo que hay que cambiar es el motor;
  * CDK2 por debajo de 0,55 -> el embudo no separa ni donde otros separan, y
    todo lo que se apoya encima (rescoring, MM-GBSA) se apoya en nada.

QUE MIDE
--------
Las MISMAS medidas que `calibracion_tbk1/validar_banco_tbk1.py`, para que los
numeros se puedan poner al lado sin traducir nada:

  * AUC crudo y AUC por atomo pesado (la medida con la que el proyecto lee sus
    listas desde el 24 de septiembre, porque el crudo paga el tamano);
  * EF1 % y EF5 %, que es lo que de verdad se pregunta un investigador;
  * BEDROC (alpha 20) de RDKit, acompanado de su base al azar MEDIDA por
    permutacion: sin eso un BEDROC de 0,4 no se sabe si es bueno o malo;
  * el intervalo bootstrap del 95 % (2000 remuestreos, semilla fija), porque la
    regla del proyecto solo declara un bloque si el limite inferior pasa de 0,5;
  * el reparto por quimiotipo (esqueletos de Murcko), que en DUD-E es
    obligatorio: los 474 activos de CDK2 vienen de series congenericas y un AUC
    global alto puede venir de una sola familia;
  * la comparacion directa con TBK1, recalculada aqui con el mismo codigo desde
    `resultados_tbk1.csv`, que es la unica forma honesta de ponerlas juntas.

Uso:
    python validar_banco_cdk2.py
    python validar_banco_cdk2.py --sin-tbk1

Salida en `analysis/calibracion_dude_cdk2/`:
    validar_cdk2.csv              ligando, papel, energia, pesados y las dos
                                  puntuaciones, para rehacer cualquier cuenta
    INFORME_CALIBRACION_CDK2.md   el informe
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
TBK1 = os.path.join(os.path.dirname(BASE), "calibracion_tbk1")

LIGANDS = os.path.join(BASE, "ligands")
BANCO = os.path.join(BASE, "banco_cdk2.csv")
RESULTADOS = os.path.join(BASE, "resultados_cdk2.csv")
ISM_ACT = os.path.join(BASE, "actives_final.ism")
ISM_DEC = os.path.join(BASE, "decoys_final.ism")
CSV_SALIDA = os.path.join(BASE, "validar_cdk2.csv")
INFORME = os.path.join(BASE, "INFORME_CALIBRACION_CDK2.md")

sys.path.insert(0, RAIZ)
import preparar_ligando as PL  # noqa: E402

B = 2000          # remuestreos del bootstrap
SEMILLA = 20260929
ALPHA = 20.0      # el alpha de BEDROC de Truchon y Bayly
N_AZAR = 200      # permutaciones para la linea de base de BEDROC

# El numero contra el que se lee todo esto. No es una medida nuestra: es la de
# los autores de DUD-E con DOCK 3.6, y por eso sirve de referencia externa.
AUC_PUBLICADO = 0.791
AUC_PUBLICADO_FUENTE = "Mysinger y col. 2012 (DUD-E), DOCK 3.6 sobre CDK2"

# El techo que dejo TBK1 el 28 de septiembre de 2026, que es lo que hay que
# explicar: o es de la diana o es del motor.
TBK1_TECHO = (0.602, 0.617, 0.637)   # crudo, por atomo, media de 9 poses


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


def leer_ism(ruta):
    """{id: (smiles, chembl)} de un `*_final.ism` de DUD-E.

    El formato es `SMILES id [chembl_id]`, separado por espacios. El id es lo que
    da nombre al PDBQT (`ACT_301074`), asi que es la llave que une el banco con
    el SMILES. Los senuelos vienen sin tercera columna.
    """
    salida = {}
    with open(ruta, encoding="utf-8", errors="ignore") as f:
        for l in f:
            p = l.split()
            if len(p) >= 2:
                salida[p[1]] = (p[0], p[2] if len(p) > 2 else "")
    return salida


def leer_tbk1():
    """(por_atomo_activos, por_atomo_decoys) de TBK1, recalculado aqui.

    Se recalcula en vez de copiar el numero del informe a proposito: es la unica
    forma de estar seguro de que los dos AUC/atomo salen del MISMO codigo, y la
    leccion del 24 de septiembre fue justo esa (una metrica no puede tener dos
    formulas). Si falta algo, devuelve None y el informe lo dice.
    """
    try:
        act = {r["molecule_chembl_id"] for r in
               csv.DictReader(open(os.path.join(TBK1, "compuestos_tbk1.csv"),
                                   encoding="utf-8"))}
        dec = {r["id"] for r in
               csv.DictReader(open(os.path.join(TBK1, "decoys_tbk1.csv"),
                                   encoding="utf-8"))}
        res = {r["ligand"]: float(r["energy"]) for r in
               csv.DictReader(open(os.path.join(TBK1, "resultados_tbk1.csv"),
                                   encoding="utf-8"))}
    except OSError:
        return None, None
    ligdir = os.path.join(TBK1, "ligands")
    valores = {1: [], 0: []}
    for etiqueta, grupo in (("ACT_", act), ("DEC_", dec)):
        for k in grupo:
            nombre = etiqueta + k
            e = res.get(nombre)
            if e is None:
                continue
            n, _ = PL.contar_fichero(os.path.join(ligdir, nombre + ".pdbqt"))
            if n:
                valores[1 if etiqueta == "ACT_" else 0].append(-e / n)
    # Un AUC sin fondo no es un AUC. El lanzador de TBK1 solo reescribe su CSV
    # al terminar, asi que mientras corre lo que hay en disco es la corrida
    # ANTERIOR (que en el banco v2 traia solo activos). Sin las dos clases no se
    # recalcula nada: se devuelve vacio y el informe cita el numero documentado,
    # que es mas honesto que imprimir el AUC de 340 activos contra cero senuelos.
    if not valores[1] or not valores[0]:
        return None, None
    return valores[1], valores[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-quimiotipo", type=int, default=5,
                    help="activos minimos para que un esqueleto cuente")
    ap.add_argument("--sin-tbk1", action="store_true",
                    help="no recalcular la comparacion con TBK1")
    args = ap.parse_args()

    if not os.path.exists(RESULTADOS):
        raise SystemExit("falta %s: primero hay que acoplar el banco"
                         % os.path.relpath(RESULTADOS, RAIZ))
    if not os.path.exists(BANCO):
        raise SystemExit("falta %s: lo escribe preparar_banco_cdk2.py"
                         % os.path.relpath(BANCO, RAIZ))

    # --- los datos ---
    # El manifiesto dice quien es activo y quien senuelo, y el .ism da el
    # SMILES. La energia sale del CSV que deja el lanzador.
    smi_act, smi_dec = leer_ism(ISM_ACT), leer_ism(ISM_DEC)
    energia = {}
    for r in csv.DictReader(open(RESULTADOS, encoding="utf-8")):
        energia[r["ligand"]] = float(r["energy"])

    filas = []
    en_banco = {"ACT": 0, "DEC": 0}
    sin_pose = {"ACT": 0, "DEC": 0}
    sin_smiles = 0
    for r in csv.DictReader(open(BANCO, encoding="utf-8")):
        nombre, etiqueta, ident = r["nombre"], r["etiqueta"], r["id"]
        en_banco[etiqueta] = en_banco.get(etiqueta, 0) + 1
        if nombre not in energia:
            sin_pose[etiqueta] = sin_pose.get(etiqueta, 0) + 1
            continue
        smiles, chembl = (smi_act if etiqueta == "ACT" else smi_dec).get(
            ident, ("", ""))
        if not smiles:
            sin_smiles += 1
        n, _ = PL.contar_fichero(os.path.join(LIGANDS, nombre + ".pdbqt"))
        if not n:
            continue
        e = energia[nombre]
        filas.append({"ligando": nombre,
                      "papel": "activos" if etiqueta == "ACT" else "decoys",
                      "ident": ident, "chembl": chembl, "smiles": smiles,
                      "energia": e, "pesados": n,
                      "crudo": -e, "por_atomo": -e / n})
    print("con pose: %d activos (%d sin pose) y %d decoys (%d sin pose)"
          % (sum(1 for f in filas if f["papel"] == "activos"), sin_pose["ACT"],
             sum(1 for f in filas if f["papel"] == "decoys"), sin_pose["DEC"]))
    if sin_smiles:
        print("sin SMILES en el .ism: %d (no entran en el reparto por quimiotipo)"
              % sin_smiles)

    A = [f for f in filas if f["papel"] == "activos"]
    D = [f for f in filas if f["papel"] == "decoys"]
    ratio = len(D) / max(len(A), 1)
    # Si al banco le falta gente, es que se acoplo la submuestra 1:10 y el
    # intervalo de abajo ya lo paga. Se dice, no se calla.
    completo = (en_banco["ACT"] and len(A) == en_banco["ACT"]
                and len(D) == en_banco["DEC"])

    lineas = []
    lineas.append("INFORME DE CALIBRACION DE CDK2 (DUD-E) — %s"
                  % time.strftime("%Y-%m-%d"))
    lineas.append("")
    lineas.append("banco: %d activos con pose y %d decoys con pose (%.1f por activo)"
                  % (len(A), len(D), ratio))
    lineas.append("activos sin pose: %d | decoys sin pose: %d | sin SMILES: %d"
                  % (sin_pose["ACT"], sin_pose["DEC"], sin_smiles))
    if not completo:
        lineas.append("ATENCION: es la SUBMUESTRA de senuelos (%d de %d), no el"
                      " banco entero." % (len(D), en_banco["DEC"]))
        lineas.append("   El AUC no depende de la proporcion (es un rango), asi"
                      " que la submuestra no")
        lineas.append("   cambia el numero esperado, solo su precision: el"
                      " intervalo de abajo ya lo dice.")

    resultados = {}
    for nombre, campo in (("crudo", "crudo"), ("por atomo pesado", "por_atomo")):
        a = [f[campo] for f in A]
        d = [f[campo] for f in D]
        auc = auc_np(a, d)
        # RDKit CalcAUC no asigna automaticamente medio acierto a empates entre
        # clases; compara solo cuando las puntuaciones no se solapan.
        if set(a).isdisjoint(d):
            rd = auc_rdkit(a, d)
            if abs(auc - rd) > 1e-6:
                lineas.append("AVISO: el AUC propio (%.6f) y el de RDKit (%.6f)"
                              " no coinciden" % (auc, rd))
        else:
            lineas.append("AUC RDKit no usado como contraste: hay puntuaciones"
                          " empatadas entre activos y decoys; auc_np las cuenta"
                          " a medias.")
        lo, hi = bootstrap_auc(a, d)
        resultados[nombre] = {"auc": auc, "lo": lo, "hi": hi,
                              "bedroc": bedroc(a, d),
                              "ef1": factor_enriquecimiento(a, d, 1.0),
                              "ef5": factor_enriquecimiento(a, d, 5.0)}
        lineas.append("")
        lineas.append("* %s:" % nombre)
        lineas.append("  AUC %.3f (bootstrap 95 %%: %.3f - %.3f) -> %s"
                      % (auc, lo, hi, "PASA" if lo > 0.5 else "NO PASA"))
        lineas.append("  EF1 %% = %.2f | EF5 %% = %.2f | BEDROC(20) = %.3f"
                      % (resultados[nombre]["ef1"], resultados[nombre]["ef5"],
                         resultados[nombre]["bedroc"]))

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
    lineas.append("   el BEDROC de arriba solo vale si esta claramente por"
                  " encima de esa base, no por encima de cero.")

    # --- quimiotipos ---
    grupos = {}
    for f in A:
        esc = esqueleto(f["smiles"])
        if esc:
            grupos.setdefault(esc, []).append(f["por_atomo"])
    validos = [(e, v) for e, v in grupos.items()
               if len(v) >= args.min_quimiotipo]
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
    lineas.append("   es la cuenta que decide: los 474 activos de CDK2 vienen de"
                  " series congenericas y un")
    lineas.append("   AUC global alto puede venir de una sola familia.")
    for n, esc, auc_esc in detalle[:12]:
        lineas.append("     %2d activos  AUC %.3f  %s" % (n, auc_esc, esc[:70]))

    # --- TBK1, recalculado con este mismo codigo ---
    if not args.sin_tbk1:
        a_t, d_t = leer_tbk1()
        lineas.append("")
        lineas.append("COMPARACION CON TBK1 (la pregunta del dia)")
        if a_t and d_t:
            lineas.append("   TBK1, recalculado aqui con ESTE mismo codigo desde"
                          " resultados_tbk1.csv:")
            lineas.append("     %d activos y %d decoys con pose | AUC/atomo %.3f"
                          % (len(a_t), len(d_t), auc_np(a_t, d_t)))
            lineas.append("     (si el CSV es el parcial del 28 de septiembre, ese"
                          " numero NO es el")
            lineas.append("      definitivo: el banco de TBK1 sigue acoplandose.")
            lineas.append("   TBK1 documentado del 28 de septiembre (banco de"
                          " 35.207): crudo %.3f,")
            lineas.append("     por atomo %.3f, media de 9 poses %.3f."
                          % TBK1_TECHO)
            lineas.append("   CDK2 medido hoy: AUC/atomo %.3f (crudo %.3f)"
                          % (resultados["por atomo pesado"]["auc"],
                             resultados["crudo"]["auc"]))
            lineas.append("   Referencia publica externa: AUC %.3f (%s)"
                          % (AUC_PUBLICADO, AUC_PUBLICADO_FUENTE))
        else:
            lineas.append("   TBK1 NO se recalcula aqui: su CSV de resultados no"
                          " trae las dos clases")
            lineas.append("   (el lanzador solo lo reescribe al terminar, y el"
                          " banco sigue acoplandose).")
            lineas.append("   Se cita lo documentado del 28 de septiembre, banco"
                          " de 35.207 poses:")
            lineas.append("     crudo %.3f | por atomo pesado %.3f | media de 9"
                          " poses %.3f." % TBK1_TECHO)
        auc_cdk2 = resultados["por atomo pesado"]["auc"]
        lo = resultados["por atomo pesado"]["lo"]
        lineas.append("")
        lineas.append("LO QUE HAY QUE LEER, DECIDIDO ANTES DE VER ESTO:")
        if lo > 0.5 and auc_cdk2 >= 0.70:
            lineas.append("   CDK2 sale bien: el embudo separa donde otros"
                          " separan, y el techo de 0,62")
            lineas.append("   es de TBK1, no de la herramienta.")
        elif lo > 0.5 and auc_cdk2 >= 0.55:
            lineas.append("   CDK2 sale a medias: el embudo ordena por encima del"
                          " azar, pero lejos del")
            lineas.append("   0,79 publicado. El techo es del MOTOR (Vina-GPU con"
                          " search_depth 20 y")
            lineas.append("   num_modes 3), no de la diana: lo que hay que"
                          " cambiar es el motor.")
        else:
            lineas.append("   CDK2 NO pasa: el embudo no separa ni donde el banco"
                          " ya esta hecho y otros")
            lineas.append("   separan. Todo lo que se apoya encima (rescoring,"
                          " MM-GBSA) se apoya en nada.")

    texto = "\n".join(lineas)
    print("\n" + texto)
    with open(INFORME, "w", encoding="utf-8") as f:
        f.write(texto + "\n")

    with open(CSV_SALIDA, "w", newline="", encoding="utf-8") as f:
        campos = ["ligando", "papel", "ident", "chembl", "smiles", "energia",
                  "pesados", "crudo", "por_atomo"]
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows([{k: v for k, v in r.items() if k in campos} for r in filas])
    print()
    print("CSV: %s" % os.path.relpath(CSV_SALIDA, RAIZ))
    print("informe: %s" % os.path.relpath(INFORME, RAIZ))
    return 0


if __name__ == "__main__":
    sys.exit(main())
