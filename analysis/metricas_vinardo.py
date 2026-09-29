#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""metricas_vinardo.py — las mismas metricas de la validacion, con Vinardo.

QUE PREGUNTA CONTESTA
---------------------
El 29 sep 2026 se vio que los tres fragmentos de Nshogoza caen al final del
ranking de TDP-43 y se sospecho que era la FUNCION DE PUNTUACION: Vina premia la
superficie enterrada, y un fragmento de 11-15 atomos pesados no puede competir
con un farmaco de 25-34. Vinardo no premia la superficie, asi que si el problema
fosse la funcion, los fragmentos recuperarian puesto.

Este script NO re-acopla nada: `_repunuar_vinardo.py` ya dejo las 576 poses de
TDP-43 puntuadas con las dos funciones en `vinardo_tdp43.csv` (la afinidad es la
misma pose, solo cambia el termino). Aqui se leen esas puntuaciones y se calculan
las metricas con `evaluar()` de `validar_sod1_limpia.py`, que es la unica copia
de AUC, EF, residual y del criterio deodetiquia: si las cuentas se copiaran,
este proyecto volveria a tener dos recetas y a discrepar.

ASI SE SABE SI EL PROBLEMA ES LA FUNCION O EL RANGO
---------------------------------------------------
Si con Vinardo los fragmentos suben, el fallo era de puntuacion y basta cambiar
la funcion. Si se quedan igual, el fallo es del RANGO: no existe ninguna funcion
que haga ganara un fragmento pequeno, porque la pregunta del cribado (separar
TDP-43 de un unidor generico de ARN) la ganan las moleculas grandes. Ese caso se
documenta y se deja de reintentar.

Uso:
    python metricas_vinardo.py --csv vinardo_tdp43.csv --target TDP43
"""
import argparse
import csv
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)

CORTE_PCT = [1, 2, 5, 10, 20]


def log(m=""):
    print(m)
    sys.stdout.flush()


def cargar_csv(ruta):
    """Lee el CSV del re-puntado -> (vina, vinardo, papel)."""
    vina, vinar, papel = {}, {}, {}
    with open(ruta, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            vina[r["ligand"]] = float(r["vina"])
            vinar[r["ligand"]] = float(r["vinardo"])
            papel[r["ligand"]] = r["papel"]
    return vina, vinar, papel


def resolver(V, nombre, dirs_ligs):
    """Traduce el nombre del CSV al fichero preparado y devuelve sus atomos pesados.

    El CSV del re-puntado llama a los positivos por su nombre desnudo
    (`fragmento_1`) y el fichero preparado lleva el prefijo `ACT_`, porque los
    positivos se prepararon desde `verdad_de_referencia.csv`. Se prueban los dos
    nombres y gana el que exista; si no hay ninguno, no se inventa el numero: el
    ligando queda fuera y se dice en el log.
    """
    for cand in (nombre, "ACT_" + nombre, nombre.replace("ACT_", "")):
        for d in dirs_ligs:
            p = os.path.join(d, cand + ".pdbqt")
            if os.path.exists(p):
                n = V.pesados_smiles(V.smi_del_fichero(p))
                if n:
                    return n
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=os.path.join(BASE, "vinardo_tdp43.csv"))
    ap.add_argument("--target", default="TDP43", choices=["TDP43", "SOD1"])
    ap.add_argument("--salida", default=None)
    args = ap.parse_args()

    import validar_sod1_limpia as V
    import validar_diana_limpia as D
    cfg = D.DIANAS[args.target]

    vina, vinar, papel = cargar_csv(args.csv)
    log("lectura de %s: %d ligandos (vina %d, vinardo %d)"
        % (os.path.basename(args.csv), len(papel), len(vina), len(vinar)))
    faltan = [k for k in papel if k not in vina or k not in vinar]
    if faltan:
        log("!! %d ligandos sin alguna de las dos puntuaciones: %s"
            % (len(faltan), ", ".join(sorted(faltan)[:8])))
        return 2

    # --- positively known: las filas de verdad_de_referencia de esta diana --------
    filas = [r for r in csv.DictReader(
        open(os.path.join(BASE, "verdad_de_referencia.csv"), encoding="utf-8"))
        if r["target"] == args.target and r["apto"] == "si"]

    # El nombre del positivo en el CSV del re-puntado puede venir con o sin `ACT_`.
    # Se usa el que exista; si no existe ninguno, el positivo no cuenta y se dice.
    def nombre_pos(r):
        n = r["ligand"]
        return "ACT_" + n if "ACT_" + n in papel else n

    quimias = {}
    for r in filas:
        quimias.setdefault(r["quimia"], []).append(nombre_pos(r))
    log("positivos de union medida: %d en %d quimiotipos"
        % (len(filas), len(quimias)))

    positivos = [nombre_pos(r) for r in filas]
    perdidos = [r["ligand"] for r in filas if nombre_pos(r) not in papel]
    if perdidos:
        log("!! positivos ausentes del CSV (no cuentan): %s" % ", ".join(perdidos))
    positivos = [p for p in positivos if p in papel]

    fondos = {"fondo blando": [], "fondo duro": []}
    for k, p in papel.items():
        if p in fondos:
            fondos[p].append(k)
    for k, v in fondos.items():
        log("%-12s %d ligandos" % (k, len(v)))

    # --- atomos pesados ----------------------------------------------------------
    # Tres carpetas, en el mismo orden que usa el validador: la del fondo blando, la
    # del fondo duro y la de la salida de la corrida (alli acaban los positivos que
    # hubo que re-acoplar, como los tres fragmentos de Nshogoza). Sin la tercera,
    # esos tres se quedarian fuera y el veredicto seria de un subconjunto distinto.
    dirs = [cfg["ligands"], cfg["ligands2"], os.path.join(cfg["salida"], "ligands")]
    pesados = {}
    sin_pesados = []
    for k in papel:
        n = resolver(V, k, dirs)
        if n:
            pesados[k] = n
        else:
            sin_pesados.append(k)
    if sin_pesados:
        log("!! %d ligandos sin atomos pesados legibles (no entran): %s"
            % (len(sin_pesados), ", ".join(sorted(sin_pesados)[:8])))
    log("atomos pesados leidos: %d de %d" % (len(pesados), len(papel)))

    pos_f = [p for p in positivos if p in pesados]
    fondos_f = {k: [x for x in v if x in pesados] for k, v in fondos.items()}

    # --- metricas ----------------------------------------------------------------
    log("")
    log("=" * 78)
    log("METRICAS CON VINA Y CON VINARDO (mismas poses, misma caja)")
    log("=" * 78)
    log("%-11s %-7s %7s %7s %9s %7s %7s %9s"
        % ("fondo", "func", "AUC", "AUC/at", "residual", "EF5%",
           "quim 5%", "veredicto"))
    filas_out = []
    for nombre_f, decoys in (("fondo blando", fondos_f["fondo blando"]),
                             ("fondo duro", fondos_f["fondo duro"]),
                             ("los dos", fondos_f["fondo blando"] + fondos_f["fondo duro"])):
        for func, sc in (("vina", vina), ("vinardo", vinar)):
            res = V.evaluar(pos_f, decoys, sc, pesados, quimias,
                            "%s (%s)" % (nombre_f, func), detalle=True)
            if not res:
                log("%-11s %-7s  (sin datos)" % (nombre_f, func))
                continue
            filas_out.append((nombre_f, func, res))
            log("%-11s %-7s %7.3f %7.3f %9.3f %7.2f %7d %9s"
                % (nombre_f, func, res["auc_crudo"], res["auc_atomo"],
                   res["auc_residual"], res["ef5"],
                   len(res["quimias_corte_5pct"]), res["veredicto"]))

    # --- el veredicto de Fredy: suben los fragmentos? ----------------------------
    log("")
    log("=" * 78)
    log("LOS TRES FRAGMENTOS DE RRM2: PUESTO CONTRA CADA FONDO, VINA -> VINARDO")
    log("=" * 78)
    log("%-24s %5s %-9s %-14s %-14s %6s"
        % ("ligando", "pesad", "quimia", "puesto con vina", "puesto con vinardo",
           "delta"))
    fila = {}
    for p in positivos:
        if "fragmento" not in p.lower() or p not in pesados:
            continue
        for nombre_f, decoys in (("blando", fondos_f["fondo blando"]),
                                 ("duro", fondos_f["fondo duro"])):
            S = {k: sc[k] for k in pos_f + decoys}
            orden = sorted(S.items(), key=lambda kv: kv[1])
            puesto = {k: i + 1 for i, (k, _) in enumerate(orden)}
            for func, sc in (("vina", vina), ("vinardo", vinar)):
                S = {k: sc[k] for k in pos_f + decoys}
                orden = sorted(S.items(), key=lambda kv: kv[1])
                puesto = {k: i + 1 for i, (k, _) in enumerate(orden)}
                q = next((q for q, m in quimias.items() if p in m), "?")
                # se imprime una vez por par (vina, vinardo) con los dos puestos
                if func == "vina":
                    pv = puesto[p]
                    fila[p, nombre_f] = pv
                else:
                    pv = fila[p, nombre_f]
                    signo = "+" if puesto[p] > pv else ""
                    log("%-24s %5d %-9s %-14s %-14s %6s"
                        % (p.replace("ACT_", "") + " / " + nombre_f,
                           pesados.get(p, 0), q[:9],
                           "%d de %d" % (pv, len(S)),
                           "%d de %d" % (puesto[p], len(S)),
                           "%s%d" % (signo, puesto[p] - pv)))

    # --- donde caen los fragmentos por tamano ------------------------------------
    log("")
    log("=" * 78)
    log("PUESTO DENTRO DE SU MISMO TAMANO (fondo duro) — la pregunta que decide")
    log("=" * 78)
    log("El puesto global mezcla dos cosas: si la quimia se reconoce y si la")
    log("molecula es grande. El puesto dentro del grupo de ligandos de tamano")
    log("parecido quita lo segundo: si el fragmento pierde tambien ahi, ningun")
    log("cambio de funcion lo arregla, porque el problema es el rango de la diana.")
    decoys = fondos_f["fondo duro"]
    for func, sc in (("vina", vina), ("vinardo", vinar)):
        log("")
        log("%s:" % func)
        log("   %-16s %6s %-16s %-18s %s"
            % ("ligando", "pesad", "puesto global", "puesto en su tamano",
               "quimias que ganan asi"))
        for p in pos_f:            # Grupo de tamano: el propio positivo y los del fondo con menos de
            # 4 atomos pesados de diferencia, la misma tolerancia que usa
            # `pareado_por_tamano` en el validador.
            grupo = [k for k in pos_f + decoys if abs(pesados[k] - pesados[p]) <= 4]
            Sg = {k: sc[k] for k in grupo}
            orden = sorted(Sg.items(), key=lambda kv: kv[1])
            puesto_g = orden.index((p, sc[p])) + 1
            S = {k: sc[k] for k in pos_f + decoys}
            puesto = {k: i + 1 for i, (k, _) in enumerate(sorted(S.items(), key=lambda kv: kv[1]))}
            # Quimias que ganan: su mejor miembro cae en la mitad superior de SU
            # propio grupo de tamano. El corte es la mediana de ese grupo y no un
            # porcentaje de la lista entera, que con grupos de 11 a 17 no
            # significaria nada.
            ganadoras = 0
            total_q = 0
            mediana_g = float(np.median(list(Sg.values())))
            for q, miembros in quimias.items():
                ms = [m for m in miembros if m in Sg]
                if not ms:
                    continue
                total_q += 1
                if min(Sg[m] for m in ms) <= mediana_g:
                    ganadoras += 1
            log("   %-16s %6d %-16s %-18s %d de %d"
                % (p.replace("ACT_", ""), pesados[p],
                   "%d de %d" % (puesto[p], len(S)),
                   "%d de %d" % (puesto_g, len(grupo)),
                   ganadoras, total_q))

    if args.salida:
        with open(args.salida, "w", encoding="utf-8") as fh:
            fh.write("fondo,func,auc_crudo,auc_atomo,auc_residual,ef5,"
                     "quimias_5pct,veredicto\n")
            for nombre_f, func, res in filas_out:
                fh.write("%s,%s,%.4f,%.4f,%.4f,%.4f,%d,%s\n"
                         % (nombre_f, func, res["auc_crudo"], res["auc_atomo"],
                            res["auc_residual"], res["ef5"],
                            len(res["quimias_corte_5pct"]), res["veredicto"]))
        log("")
        log("escrito %s" % args.salida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
