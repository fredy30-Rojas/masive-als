#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""rescoring_estrategia_estratos.py — ¿quitar el sesgo de tamaño cambia el veredicto?

LA PREGUNTA
-----------
La validacion de SOD1 contra el receptor limpio dejo esto medido (`INFORME_VALIDACION_SOD1_LIMPIA_2026-09-21.md`):

  * AUC crudo de Vina **0,439**: la lista ordenada por energia NO discrimina;
  * pero **emparejando cada positivo con señuelos de su mismo tamaño, 5 de 8
    quimiotipos ganan**.

O sea: el embudo reconoce quimica cuando se le pregunta en igualdad de tamaño, y lo
que no sabe es ORDENAR, porque su energia es casi una funcion del numero de atomos
(dentro del fondo, correlacion -0,720; -0,106 kcal/mol por carbono pesado).

Aqui se prueba si otras funciones de puntuacion tienen el mismo defecto, aplicadas
a las MISMAS poses:

  * `vina`    — la que se uso para acoplar (control: debe reproducir la energia del
                acoplamiento; si no, algo va mal en el montaje);
  * `vinardo` — la otra funcion incluida en AutoDock Vina 1.2.3;
  * `mmgbsa`  — el rescoring fisico (openmm + GAFF + OBC2) que corre en Oracle, sobre
                una muestra emparejada por tamaño (`preparar_mmgbsa_estratos.py`).

Y para cada una se miden las cuatro cosas, con el mismo criterio de siempre:

  1. AUC crudo;
  2. correlacion con el tamaño DENTRO del fondo;
  3. AUC residual (descontada la recta del tamaño);
  4. **el veredicto**: cuantos de los 8 quimiotipos baten a los señuelos DE SU MISMO
     TAMAÑO (corte +-2 atomos pesados). Es la medida que decide.

Uso:
  python rescoring_estrategia_estratos.py              # vina + vinardo (completo)
  python rescoring_estrategia_estratos.py --sin-recalcular   # usa lo ya calculado
Salida: rescoring_estrategia_estratos.log y .csv
"""
import argparse
import csv
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
GPU = os.path.join(RAIZ, "gpu_dock")
VINA = os.path.join(RAIZ, "tools", "vina.exe")
RECEPTOR = os.path.join(GPU, "SOD1_limpio.pdbqt")
# La caja Trp32 de la validacion. `--score_only` NO se puede omitir aunque la
# documentacion diga que si: sin ella responde "ERROR: Grid box dimensions must be
# greater than 0 Angstrom" y no puntua nada.
CENTRO = (46.5, 80.0, 73.3)
TAMANO = 22
POSES = os.path.join(BASE, "validacion_SOD1_v5", "out")
LIGS = os.path.join(BASE, "validacion_SOD1_v5", "ligands")
MMGBSA = os.path.join(BASE, "mmgbsa_estratos", "resultado_oracle.csv")
SALIDA = os.path.join(BASE, "rescoring_estrategia_estratos")

HILOS = 8
MARGEN_TAMANO = 2      # atomos pesados de tolerancia al emparejar
BINS = [(0, 12), (13, 17), (18, 22), (23, 27), (28, 99)]


def log(m):
    print(m, flush=True)
    with open(SALIDA + ".log", "a", encoding="utf-8") as f:
        f.write(m + "\n")


def pesados_smiles(smi):
    m = Chem.MolFromSmiles(smi) if smi else None
    return None if m is None else m.GetNumHeavyAtoms()


def smi_de_texto(txt):
    for l in txt.splitlines():
        if l.startswith("REMARK SMILES "):
            return l[len("REMARK SMILES "):].strip()
        if l.startswith(("ATOM", "HETATM")):
            break
    return None


def afinidad_vina(ruta):
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            try:
                return float(l.split()[3])
            except (IndexError, ValueError):
                return None
    return None


# ------------------------------------------------------------------ score_only
def _score_una(t):
    """Primera pose de un fichero, puntuada con `--score_only`. Devuelve (nombre, E)."""
    vina, receptor, ruta, scoring, nombre = t
    txt = open(ruta, encoding="utf-8", errors="ignore").read()
    # `--score_only` exige un fichero SIN etiquetas MODEL: la pose trae tres y con
    # solo recortar por ENDMDL sigue respondiendo "PDBQT parsing error: Unexpected
    # multi-MODEL tag found in flex residue or ligand PDBQT file". Hay que quedarse
    # con la primera pose y quitarle tambien las lineas MODEL/ENDMDL.
    if "ENDMDL" in txt:
        txt = txt.split("ENDMDL")[0]
    txt = "\n".join(l for l in txt.splitlines()
                    if not l.startswith(("MODEL", "ENDMDL"))) + "\n"
    fd, tmp = tempfile.mkstemp(suffix=".pdbqt")
    os.close(fd)
    try:
        with open(tmp, "w", encoding="utf-8") as f:
            f.write(txt)
        p = subprocess.run([vina, "--receptor", receptor, "--ligand", tmp,
                            "--center_x", str(CENTRO[0]),
                            "--center_y", str(CENTRO[1]),
                            "--center_z", str(CENTRO[2]),
                            "--size_x", str(TAMANO), "--size_y", str(TAMANO),
                            "--size_z", str(TAMANO),
                            "--score_only", "--scoring", scoring],
                           capture_output=True, text=True, timeout=600)
        for l in p.stdout.splitlines():
            if l.startswith("Estimated Free Energy of Binding") or \
               l.strip().startswith("Affinity:"):
                try:
                    return nombre, float(l.split(":")[1].split()[0])
                except (IndexError, ValueError):
                    pass
        return nombre, None
    except Exception:
        return nombre, None
    finally:
        try:
            os.remove(tmp)
        except OSError:
            pass


def score_only(scoring, nombres_ficheros):
    tareas = [(VINA, RECEPTOR, ruta, scoring, nombre)
              for nombre, ruta in nombres_ficheros.items()]
    out = {}
    with ProcessPoolExecutor(max_workers=HILOS) as ex:
        for nombre, e in ex.map(_score_una, tareas):
            if e is not None:
                out[nombre] = e
    return out


# ------------------------------------------------------------------- metricas
def auc(act, dec):
    if not act or not dec:
        return None
    g = 0.0
    for a in act:
        for d in dec:
            g += 0.5 if a == d else (1.0 if a < d else 0.0)
    return g / (len(act) * len(dec))


def ef(act, dec, pct):
    todo = sorted([(s, 1) for s in act] + [(s, 0) for s in dec])
    k = max(1, int(round(len(todo) * pct / 100.0)))
    n_top = sum(1 for _, lab in todo[:k] if lab == 1)
    return (n_top / len(act)) / (pct / 100.0) if act else None


def quimios_recuperadas(scores, positivos, decoys, pesados, quimias):
    """Cuantos quimiotipos baten a los señuelos de su mismo tamaño (la medida justa)."""
    ganan = {}
    for q, miembros in quimias.items():
        for n in miembros:
            if n not in scores:
                continue
            tam = pesados[n]
            par = [d for d in decoys if abs(pesados[d] - tam) <= MARGEN_TAMANO]
            if len(par) < 5:
                continue
            a = auc([scores[n]], [scores[d] for d in par])
            if a is not None and a > 0.5:
                ganan.setdefault(q, []).append(n)
    return ganan


def analizar(etiqueta, scores, positivos, decoys, pesados, quimias):
    a = [scores[k] for k in positivos if k in scores]
    d = [scores[k] for k in decoys if k in scores]
    if not a or not d:
        return None
    xs = np.array([pesados[k] for k in decoys if k in scores], dtype=float)
    ys = np.array([scores[k] for k in decoys if k in scores], dtype=float)
    corr = float(np.corrcoef(xs, ys)[0, 1])
    b, a0 = np.polyfit(xs, ys, 1)
    resid = {k: scores[k] - (a0 + b * pesados[k]) for k in scores}
    auc_crudo = auc(a, d)
    auc_resid = auc([resid[k] for k in positivos if k in scores],
                    [resid[k] for k in decoys if k in scores])
    ganan = quimios_recuperadas(scores, positivos, decoys, pesados, quimias)

    log("")
    log("=" * 78)
    log("%s" % etiqueta)
    log("=" * 78)
    log("   AUC crudo %.3f | residual %.3f | EF5%% %.2f"
        % (auc_crudo, auc_resid, ef(a, d, 5.0)))
    log("   correlacion con el tamaño DENTRO del fondo: %+.3f"
        " (%+.3f kcal/mol por carbono)" % (corr, b))
    log("   quimiotipos que baten a los señuelos de su mismo tamaño: %d de %d -> %s"
        % (len(ganan), len(quimias), ", ".join(sorted(ganan)) or "ninguno"))

    # --- por estratos de tamaño, que es lo que se pidio ver ---
    log("")
    log("   por estratos de tamaño (AUC dentro de cada estrato):")
    log("      %-12s %5s %6s %8s" % ("estrato", "pos", "señuel", "AUC"))
    pesos, valores = [], []
    for lo, hi in BINS:
        pa = [scores[k] for k in positivos if k in scores and lo <= pesados[k] <= hi]
        de = [scores[k] for k in decoys if k in scores and lo <= pesados[k] <= hi]
        v = auc(pa, de)
        log("      %-12s %5d %6d %8s"
            % ("%d-%d atomos" % (lo, hi), len(pa), len(de),
               "-" if v is None else "%.3f" % v))
        if v is not None and pa and de:
            pesos.append(len(pa) + len(de))
            valores.append(v)
    estratificado = (np.average(valores, weights=pesos)
                     if valores else None)
    log("      %-12s %5s %6s %8s" % ("ESTRATIFICADO", "", "",
                                     "-" if estratificado is None else "%.3f"
                                     % estratificado))

    return {"etiqueta": etiqueta, "auc_crudo": auc_crudo, "auc_residual": auc_resid,
            "ef5": ef(a, d, 5.0), "correlacion_tamano": corr,
            "pendiente_kcal_por_atomo": b,
            "quimios_recuperadas": len(ganan), "quimias": len(quimias),
            "auc_estratificado": estratificado,
            "detalle_quimias": ";".join(sorted(ganan))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sin-recalcular", action="store_true")
    args = ap.parse_args()

    open(SALIDA + ".log", "w", encoding="utf-8").close()
    log("SESGO DE TAMAÑO Y RESCORING POR ESTRATOS   %s"
        % __import__("time").strftime("%Y-%m-%d %H:%M"))
    log("receptor %s" % os.path.basename(RECEPTOR))

    # ------------------------------------------------------------------- ligandos
    positivos_meta = {}
    for r in csv.DictReader(open(os.path.join(BASE, "verdad_de_referencia.csv"),
                                 encoding="utf-8")):
        if r["target"] == "SOD1" and r["apto"] == "si":
            positivos_meta["ACT_" + r["ligand"]] = r["quimia"]
    quimias = {}
    for k, q in positivos_meta.items():
        quimias.setdefault(q, []).append(k)

    ficheros, pesados, positivos, decoys = {}, {}, [], []
    for p in sorted(os.listdir(POSES)):
        if not p.endswith("_out.pdbqt"):
            continue
        nombre = p.replace("_out.pdbqt", "")
        ruta = os.path.join(POSES, p)
        lig = os.path.join(LIGS, nombre + ".pdbqt")
        if not os.path.exists(lig):
            continue
        n = pesados_smiles(smi_de_texto(open(lig, encoding="utf-8").read()))
        if n is None:
            continue
        ficheros[nombre] = ruta
        pesados[nombre] = n
        (positivos if nombre.startswith("ACT_") else decoys).append(nombre)
    positivos = [k for k in positivos if k in positivos_meta]
    log("positivos: %d en %d quimias | fondo: %d"
        % (len(positivos), len(quimias), len(decoys)))

    resultados = []
    series = []          # (etiqueta, puntuaciones) para la comparacion justa
    cache = SALIDA + "_cache.csv"

    # --------------------------------------------------- funcs de AutoDock Vina
    for scoring in ("vina", "vinardo"):
        scores = None
        if args.sin_recalcular and os.path.exists(cache):
            scores = {}
            for r in csv.DictReader(open(cache, encoding="utf-8")):
                v = r.get(scoring, "")
                if v:
                    scores[r["ligand"]] = float(v)
        if not scores:
            log("")
            log("puntuando %d poses con --score_only --scoring %s..."
                % (len(ficheros), scoring))
            scores = score_only(scoring, ficheros)
            log("   puntuadas: %d" % len(scores))
            # se guarda en la cache compartida
            todas = {}
            for r in csv.DictReader(open(cache, encoding="utf-8")) \
                    if os.path.exists(cache) else []:
                todas[r["ligand"]] = r
            for k, v in scores.items():
                todas.setdefault(k, {"ligand": k})[scoring] = "%.4f" % v
            with open(cache, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=["ligand", "vina", "vinardo"])
                w.writeheader()
                for k, r in sorted(todas.items()):
                    w.writerow({"ligand": k, "vina": r.get("vina", ""),
                                "vinardo": r.get("vinardo", "")})

        # el acoplado original, para comparar
        acoplado = {k: afinidad_vina(ficheros[k]) for k in ficheros}
        acoplado = {k: v for k, v in acoplado.items() if v is not None}
        if scoring == "vina":
            resultados.append(analizar(
                "CONTROL: energia del ACOPLAMIENTO (toda la lista, N=%d)" % len(acoplado),
                acoplado, positivos, decoys, pesados, quimias))
            series.append(("ACOPLAMIENTO", acoplado))
        resultados.append(analizar(
            "RESCORING %s (--score_only sobre las mismas poses, N=%d)"
            % (scoring.upper(), len(scores)),
            scores, positivos, decoys, pesados, quimias))
        series.append(("RESCORING %s" % scoring.upper(), scores))

    # ------------------------------------------------------------------- MM-GBSA
    if os.path.exists(MMGBSA):
        filas = list(csv.DictReader(open(MMGBSA, encoding="utf-8")))
        mm = {}
        for r in filas:
            try:
                mm[r["ligand"]] = float(r["mmgbsa_dG"])
            except (KeyError, TypeError, ValueError):
                continue
        if mm:
            # la muestra tiene su propio fondo: se analiza sobre ella
            decoys_mm = [k for k in mm if k in decoys]
            pos_mm = [k for k in mm if k in positivos]
            resultados.append(analizar(
                "RESCORING MM-GBSA (muestra emparejada por tamaño: %d positivos,"
                " %d señuelos)" % (len(pos_mm), len(decoys_mm)),
                mm, positivos, decoys_mm, pesados, quimias))
            series.append(("RESCORING MM-GBSA", mm))

            # No vale comparar un AUC medido sobre 82 ligandos emparejados con otro
            # medido sobre 513: la muestra emparejada es mas facil para cualquier
            # funcion. Aqui se pasan las CUATRO funciones por los MISMOS ligandos.
            muestra = set(mm)
            justas = []
            log("")
            log("=" * 78)
            log("COMPARACION JUSTA: las mismas funciones sobre los MISMOS %d ligandos"
                % len(muestra))
            log("=" * 78)
            for etq, sc in series:
                filtrado = {k: v for k, v in sc.items() if k in muestra}
                r = analizar("EN LA MUESTRA DEL MM-GBSA (%s)" % etq, filtrado,
                             positivos, decoys_mm, pesados, quimias)
                if r:
                    r["etiqueta"] = etq
                    justas.append(r)
            if justas:
                with open(SALIDA + "_muestra.csv", "w", newline="",
                          encoding="utf-8") as f:
                    w = csv.DictWriter(f, fieldnames=list(justas[0].keys()))
                    w.writeheader()
                    for r in justas:
                        w.writerow(r)
                log("")
                log("guardado: %s_muestra.csv" % SALIDA)
        else:
            log("")
            log("MM-GBSA: %s existe pero sin valores utilizables" % MMGBSA)
    else:
        log("")
        log("MM-GBSA: aun no hay resultados en %s" % os.path.basename(MMGBSA))

    # ------------------------------------------------------------------ resumen
    validos = [r for r in resultados if r]
    if validos:
        log("")
        log("=" * 78)
        log("RESUMEN: ninguno de los numeros clasicos decide; decide la ultima columna")
        log("=" * 78)
        log("   %-34s %8s %9s %8s %6s %10s"
            % ("funcion", "AUC", "residual", "corr-tam", "EF5%", "quimias"))
        for r in validos:
            log("   %-34s %8.3f %9.3f %+8.3f %6.2f %6d/%d"
                % (r["etiqueta"][:34], r["auc_crudo"], r["auc_residual"],
                   r["correlacion_tamano"], r["ef5"], r["quimios_recuperadas"],
                   r["quimias"]))
        with open(SALIDA + ".csv", "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(validos[0].keys()))
            w.writeheader()
            for r in validos:
                w.writerow(r)
        log("")
        log("guardado: %s.csv" % SALIDA)
    return 0


if __name__ == "__main__":
    sys.exit(main())
