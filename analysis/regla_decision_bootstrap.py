#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""regla_decision_bootstrap.py — la regla que decide, con su error medido.

POR QUE
-------
El corte por quimiotipos no aguanta: con el mismo motor y el mismo conjunto da 2 de 5
a exhaustividad 8 y 16 y 3 de 5 a 32 (una quimia cruza por 0,59 kcal/mol), y al cambiar
de motor da 1 de 5 (INFORME_BARRIDO_EXHAUSTIVIDAD_TDP43_2026-09-24.md §1). El corte mide
donde cae la cola de un fondo de 122 ligandos, y esa cola se mueve con el esfuerzo.

Esta es la sustituta. La medida es el **AUC por atomo pesado**, que es el unico numero
que quedo por encima del azar en las CUATRO corridas (0,666-0,738) y el unico estable.
Lo que le faltaba era el error, y el error se mide remuestreando.

LA REGLA, FIJADA POR ADELANTADO
-------------------------------
    Medida: AUC por atomo pesado contra el fondo DURO (verdaderos unidores de ARN),
    que es el unico fondo que puede separar quimia especifica de quimia generica.

    Se declara "reconoce quimia" solo si el limite INFERIOR del intervalo del 95%
    es mayor que 0,5 en las DOS formas de remuestreo:

      (a) remuestreando el FONDO (2000 veces, con reemplazo) — mide la cola del
          fondo, que es exactamente lo que hacia bailar el corte por quimiotipos;
      (b) remuestreando el fondo Y LOS POSITIVOS (2000) — mide ademas el muestreo de
          positivos, que con 7 ligandos es lo que de verdad aprieta.

    Si solo lo supera (a): NO se declara. Queda "sin evidencia" y lo que hay que
    ampliar es el numero de positivos, no la exhaustividad.

Y una segunda condicion, que es la que hace util la regla: **el veredicto tiene que ser
el mismo en todas las corridas**. Una regla que cambia de veredicto al cambiar el motor
o el esfuerzo no sirve, aunque su intervalo sea bonito. Eso lo comprueba `comparar()`.

Los tres veredictos posibles son, y no hay mas:

  * PASA           — el limite inferior de (a) y de (b) supera 0,5. Se puede reportar.
  * SIN EVIDENCIA  — (a) si, (b) no: el orden no es un capricho del fondo, pero con
                     esos positivos no se puede afirmar que reconozca quimia.
  * NO PASA        — ni (a). El embudo no ordena mejor que el azar.

Se importan `auc` y `residual` de `validar_sod1_limpia`, que es la unica copia de las
metricas: la regla no puede calcular su numero con una formula distinta a la del informe
que puntua (leccion del 24 sep con `auditar_tipos_cribado`).

Uso:
    python regla_decision_bootstrap.py
    python regla_decision_bootstrap.py --nombre FUS "exh 8=_barrido/validar_exh8.csv"
Salida: regla_decision.log y regla_decision.csv
"""
import argparse
import csv
import os
import sys
import time

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)

from validar_sod1_limpia import auc, residual        # noqa: E402  (unica copia)

# --------------------------------------------------------------- fijado por adelantado
B = 2000                 # remuestreos
NIVEL = 95.0             # nivel del intervalo
UMBRAL = 0.5             # por debajo de esto el embudo no ordena
SEMILLA = 20260924       # para que el intervalo sea reproducible

SALIDA = os.path.join(BASE, "regla_decision")

# Las corridas que tienen que dar el MISMO veredicto, que es la prueba de la regla.
TDP43 = [
    ("exh 8 (CPU)", os.path.join(BASE, "_barrido_TDP43", "validar_exh8.csv")),
    ("exh 16 (CPU)", os.path.join(BASE, "validar_tdp43_limpia.csv")),
    ("exh 32 (CPU)", os.path.join(BASE, "_barrido_TDP43", "validar_exh32.csv")),
    ("GPU (search_depth 20)", os.path.join(BASE, "_control_gpu_TDP43", "validar_gpu.csv")),
]
SOD1 = [("exh 8 (CPU)", os.path.join(BASE, "validar_sod1_limpia.csv"))]

BLOQUES = [("blando", "fondo", "señuelos emparejados por tamaño"),
           ("duro", "fondo2", "unidores de ARN de R-BIND 2.0"),
           ("juntos", ("fondo", "fondo2"), "los dos fondos juntos")]


def log(m):
    print(m, flush=True)
    os.makedirs(SALIDA, exist_ok=True)
    with open(os.path.join(SALIDA, "regla_decision.log"), "a",
              encoding="utf-8") as f:
        f.write(m + "\n")


def leer_corrida(ruta):
    """Devuelve (positivos, bloques) de un CSV ya puntuado.

    Se lee del CSV de la validacion y no de las poses: asi la regla se aplica al
    MISMO numero que ya esta en el informe, y no a un recalculo paralelo.
    """
    scores, pesados, papeles = {}, {}, {}
    with open(ruta, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            k = r["ligand"]
            scores[k] = float(r["afinidad"])
            pesados[k] = int(r["pesados"])
            papeles[k] = r["papel"]
    positivos = [k for k in scores if papeles[k] == "positivo"]
    hay = {p: any(v == p for v in papeles.values()) for p in ("fondo", "fondo2")}
    bloques = {}
    for nombre, papeles_bloque, _ in BLOQUES:
        # Un bloque solo existe si SU fondo existe: sin fondo duro no hay bloque duro,
        # y sin los dos no hay "juntos" (si no, SOD1 saldria dos veces con el mismo
        # fondo y pareceria que se midio dos cosas distintas).
        if isinstance(papeles_bloque, str):
            papeles_bloque = (papeles_bloque,)
        if any(not hay.get(p, False) for p in papeles_bloque):
            continue
        fondo = [k for k in scores if papeles[k] in papeles_bloque]
        if fondo:
            bloques[nombre] = fondo
    return positivos, scores, pesados, bloques


def auc_np(act, dec):
    """El MISMO AUC de validar_sod1_limpia (empates a medias), vectorizado.

    `comprobar_auc()` verifica que coincide al decimal con el importado antes de
    usarlo en los 2000 remuestreos.
    """
    act = np.asarray(act, dtype=float)
    dec = np.asarray(dec, dtype=float)
    if act.size == 0 or dec.size == 0:
        return None
    g = (act[:, None] < dec[None, :]).sum() + 0.5 * (act[:, None] == dec[None, :]).sum()
    return float(g) / (act.size * dec.size)


def comprobar_auc(positivos, fondo, scores):
    a = [scores[k] for k in positivos]
    d = [scores[k] for k in fondo]
    lento, rapido = auc(a, d), auc_np(a, d)
    if abs(lento - rapido) > 1e-12:
        raise SystemExit("el AUC vectorizado no coincide con el importado: %r vs %r"
                         % (lento, rapido))


def metricas(positivos, fondo, scores, pesados, rng, remuestrear_positivos):
    """Un remuestreo del bloque, con las tres medidas.

    `remuestrear_positivos=False` -> solo se mueve el fondo (la cola del fondo es lo
    que se quiere medir). Con True se mueven los dos, que es el intervalo honesto.
    """
    n_f, n_p = len(fondo), len(positivos)

    def toma(pool, n):
        idx = rng.integers(0, n, size=n)
        return np.array(pool, dtype=float)[idx]

    S = {k: scores[k] for k in positivos + fondo}
    P = {k: pesados[k] for k in positivos + fondo}
    if remuestrear_positivos:
        pos = [positivos[i] for i in rng.integers(0, n_p, size=n_p)]
    else:
        pos = positivos
    dec = [fondo[i] for i in rng.integers(0, n_f, size=n_f)]

    # el residual se reajusta al fondo remuestreado, como manda el metodo
    resid, _, _ = residual(S, P, dec)
    return (auc_np([S[k] for k in pos], [S[k] for k in dec]),
            auc_np([S[k] / P[k] for k in pos], [S[k] / P[k] for k in dec]),
            auc_np([resid[k] for k in pos], [resid[k] for k in dec]))


def intervalo(positivos, fondo, scores, pesados, solo_fondo):
    rng = np.random.default_rng(SEMILLA)
    crudo, atomo, res = [], [], []
    for _ in range(B):
        c, a, r = metricas(positivos, fondo, scores, pesados, rng,
                           remuestrear_positivos=not solo_fondo)
        crudo.append(c), atomo.append(a), res.append(r)
    lo = (100.0 - NIVEL) / 2.0
    hi = 100.0 - lo
    return {n: (float(np.percentile(v, lo)), float(np.percentile(v, hi)))
            for n, v in (("crudo", crudo), ("atomo", atomo), ("residual", res))}


def evaluar_bloque(positivos, fondo, scores, pesados):
    comprobar_auc(positivos, fondo, scores)
    S = {k: scores[k] for k in positivos + fondo}
    P = {k: pesados[k] for k in positivos + fondo}
    ic_fondo = intervalo(positivos, fondo, scores, pesados, solo_fondo=True)
    ic_todo = intervalo(positivos, fondo, scores, pesados, solo_fondo=False)
    resid, a, b = residual(S, P, fondo)
    return {
        "n_pos": len(positivos), "n_fondo": len(fondo),
        "auc_crudo": auc([S[k] for k in positivos], [S[k] for k in fondo]),
        "auc_atomo": auc([S[k] / P[k] for k in positivos],
                         [S[k] / P[k] for k in fondo]),
        "auc_residual": auc([resid[k] for k in positivos],
                            [resid[k] for k in fondo]),
        "ic_fondo": ic_fondo, "ic_todo": ic_todo, "pendiente": b,
    }


def veredicto(r):
    """La regla, sin margen de interpretacion."""
    if r["ic_todo"]["atomo"][0] > UMBRAL:
        return "PASA"
    if r["ic_fondo"]["atomo"][0] > UMBRAL:
        return "SIN EVIDENCIA"
    return "NO PASA"


def positivos_necesarios(r):
    """Cuantos positivos harian falta para que el limite inferior tocara el umbral.

    El error del AUC/atomo cuando se remuestrean los positivos baja como 1/sqrt(n), y
    con estos numeros esa parte es la que manda (la del fondo es cinco veces menor).
    Es una estimacion, no una promesa: supone que el embudo conserva el mismo orden
    al añadir positivos, y eso no se sabe hasta tenerlos. Devuelve None si el AUC/atomo
    ni siquiera esta por encima del umbral: entonces no hay positivos que basten.
    """
    A = r["auc_atomo"]
    lo, hi = r["ic_todo"]["atomo"]
    if A <= UMBRAL or hi <= lo:
        return None
    se = (hi - lo) / (2 * 1.96)          # error observado con n_pos positivos
    se_obj = (A - UMBRAL) / 1.96         # error que haria falta para tocar el umbral
    if se <= se_obj:
        return r["n_pos"]
    return int(np.ceil(r["n_pos"] * (se / se_obj) ** 2))


def pintar(etiqueta, bloque, r):
    log("")
    log("-" * 78)
    log("%s | bloque %s (%s)" % (etiqueta, bloque,
                                 BLOQUES[[b[0] for b in BLOQUES].index(bloque)][2]))
    log("   %d positivos contra %d del fondo" % (r["n_pos"], r["n_fondo"]))
    log("   %-16s %-8s %-22s %-22s" % ("medida", "punto", "IC 95% (solo fondo)",
                                       "IC 95% (fondo y positivos)"))
    for nombre, etiqueta_m in (("crudo", "AUC crudo"), ("atomo", "AUC/atomo"),
                               ("residual", "AUC residual")):
        p = r["auc_" + nombre]
        f1 = r["ic_fondo"][nombre]
        f2 = r["ic_todo"][nombre]
        log("   %-16s %-8.3f [%.3f, %.3f]%-8s [%.3f, %.3f]%s"
            % (etiqueta_m, p, f1[0], f1[1],
               "" if f1[0] > UMBRAL else "  <", f2[0], f2[1],
               "" if f2[0] > UMBRAL else "  <"))
    log("   pendiente de tamaño %+.4f kcal/mol por atomo pesado" % r["pendiente"])
    log("   VEREDICTO DEL BLOQUE: %s" % veredicto(r))
    falta = positivos_necesarios(r)
    if veredicto(r) == "SIN EVIDENCIA":
        log("      El AUC/atomo aguanta el remuestreo del FONDO, pero no el de los")
        log("      positivos: con %d ligandos no se puede afirmar que el embudo" % r["n_pos"])
        log("      reconozca quimia. Lo que falta son positivos, no exhaustividad.")
        if falta:
            log("      Para que el limite inferior tocara el %s harian falta ~%d"
                " positivos (hoy %d)." % (UMBRAL, falta, r["n_pos"]))
        else:
            log("      Y con este AUC/atomo (%.3f, en el azar o por debajo) ningun"
                " numero de positivos lo arregla." % r["auc_atomo"])


def comparar(nombre, corridas):
    """El veredicto tiene que ser el MISMO en todas las corridas. Eso es la prueba."""
    log("")
    log("=" * 78)
    log("PRUEBA DE LA REGLA: %s, el mismo conjunto corrido de %d %s"
        % (nombre, len(corridas), "manera" if len(corridas) == 1 else "maneras"))
    log("=" * 78)
    log("   %-34s %-9s %-9s %-9s %-9s %s"
        % ("corrida (bloque)", "AUC/at.", "IC f. inf", "IC t. inf", "residual", "veredicto"))
    filas, veredictos = [], {}
    for etiqueta, ruta in corridas:
        if not os.path.exists(ruta):
            log("   %-24s FALTA %s" % (etiqueta, ruta))
            continue
        positivos, scores, pesados, bloques = leer_corrida(ruta)
        for nombre_bloque, fondo in bloques.items():
            r = evaluar_bloque(positivos, fondo, scores, pesados)
            v = veredicto(r)
            filas.append({"target": nombre, "corrida": etiqueta, "bloque": nombre_bloque,
                          "n_pos": r["n_pos"], "n_fondo": r["n_fondo"],
                          "auc_crudo": round(r["auc_crudo"], 4),
                          "auc_atomo": round(r["auc_atomo"], 4),
                          "atomo_ic_fondo_lo": round(r["ic_fondo"]["atomo"][0], 4),
                          "atomo_ic_fondo_hi": round(r["ic_fondo"]["atomo"][1], 4),
                          "atomo_ic_todo_lo": round(r["ic_todo"]["atomo"][0], 4),
                          "atomo_ic_todo_hi": round(r["ic_todo"]["atomo"][1], 4),
                          "auc_residual": round(r["auc_residual"], 4),
                          "residual_ic_fondo_lo": round(r["ic_fondo"]["residual"][0], 4),
                          "residual_ic_todo_lo": round(r["ic_todo"]["residual"][0], 4),
                          "positivos_que_harian_falta": (
                              "" if v == "PASA" else (positivos_necesarios(r) or "-")),
                          "veredicto": v})
            if nombre_bloque == "duro":
                veredictos.setdefault(etiqueta, v)
            log("   %-34s %-9.3f %-9.3f %-9.3f %-9.3f %s"
                % ("%s (%s)" % (etiqueta, nombre_bloque), r["auc_atomo"],
                   r["ic_fondo"]["atomo"][0], r["ic_todo"]["atomo"][0],
                   r["auc_residual"], v))
            pintar(etiqueta, nombre_bloque, r)
    if veredictos:
        unicos = set(veredictos.values())
        log("")
        if len(veredictos) == 1:
            log("   UNA SOLA CORRIDA: no hay nada que comparar; estabilidad sin probar."
                " Su veredicto sobre el fondo duro es %s." % unicos.pop())
        elif len(unicos) == 1:
            log("   LA REGLA AGUANTA: las %d corridas dan el mismo veredicto sobre el"
                " fondo duro -> %s" % (len(veredictos), unicos.pop()))
        else:
            log("   LA REGLA NO AGUANTA: los veredictos cambian segun la corrida -> %s"
                % ", ".join("%s=%s" % kv for kv in sorted(veredictos.items())))
    return filas


def main():
    global B
    ap = argparse.ArgumentParser()
    ap.add_argument("corridas", nargs="*", metavar="ETIQUETA=RUTA.CSV",
                    help="corridas a comparar; por defecto, las de TDP-43 y SOD1")
    ap.add_argument("--nombre", default="corrida",
                    help="nombre de la diana cuando se pasan corridas a mano")
    ap.add_argument("--remuestreos", type=int, default=B)
    args = ap.parse_args()
    B = args.remuestreos

    os.makedirs(SALIDA, exist_ok=True)
    open(os.path.join(SALIDA, "regla_decision.log"), "w", encoding="utf-8").close()
    log("REGLA DE DECISION CON ERROR MEDIDO (bootstrap)   %s"
        % time.strftime("%Y-%m-%d %H:%M"))
    log("medida: AUC por atomo pesado | umbral %s | %d remuestreos | nivel %s%%"
        % (UMBRAL, B, NIVEL))
    log("semilla %d (reproducible)" % SEMILLA)
    log("")
    log("   (a) IC solo fondo      = mueve la cola del fondo, que es lo que hacia")
    log("                            bailar el corte por quimiotipos")
    log("   (b) IC fondo+positivos = mide ademas el muestreo de los positivos")
    log("   PASA si el limite INFERIOR de (a) y de (b) supera %s." % UMBRAL)

    if args.corridas:
        pares = []
        for c in args.corridas:
            etiqueta, ruta = c.split("=", 1) if "=" in c else (c, c)
            pares.append((etiqueta, ruta))
        filas = comparar(args.nombre, pares)
    else:
        filas = comparar("TDP-43", TDP43) + comparar("SOD1", SOD1)

    # El AUC importado y el vectorizado tienen que dar lo mismo: se deja escrito.
    log("")
    log("   (comprobado en cada bloque: el AUC vectorizado coincide al decimal con el")
    log("    de validar_sod1_limpia, que es la unica copia de la metrica)")

    ruta = os.path.join(SALIDA, "regla_decision.csv")
    with open(ruta, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0]))
        w.writeheader()
        w.writerows(filas)
    log("")
    log("guardado: %s" % ruta)
    return 0


if __name__ == "__main__":
    sys.exit(main())
