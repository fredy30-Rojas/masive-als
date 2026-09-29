#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""bootstrap_auc.py — ¿se puede distinguir el AUC del azar con estos positivos?

LA PREGUNTA QUE ESTA NUNCA SE HA HECHO
--------------------------------------
En todo el proyecto se han reportado AUC crudos (0,301 / 0,357 / 0,332 en TDP-43 y
0,440 / 0,313 / 0,409 en SOD1) y se ha concluding de ellos que el embudo NO reconoce
quimia por encima del azar. Pero un AUC con 7 positivos no es un numero: es un numero
con un margen, y si el margen se come la mitad, la conclusion "no pasa" no es una
medida de que el embudo falle, es una medida de que **no hay poder estadistico**.

Esto no es un matiz de Statistical. Cambia lo que se puede pedir al laboratorio y lo
que se puede escribir en el paper. Con 7 positivos, ninguna diferencia entre dos
metodos de puntuacion es medible, por muy grande que parezca el numero.

EL METODO
---------
Se re-muestrea **solo el conjunto de positivos** (con reemplazo) 5.000 veces contra
el fondo fijo, y se lee el percentil 2,5 y el 97,5 del AUC. El fondo no se re-muestrea
porque son cientos de ligandos y su error es despreciable al lado del de siete.

Se lee tambien la probabilidad de que el AUC real sea mayor que 0,5, que es la pregunta
directa: ¿esta el embudo por encima o por debajo del azar, con los datos que hay?

Uso:
    python bootstrap_auc.py
"""
import csv
import os
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)

N_BOOT = 5000
SEMILLA = 20260929


def log(m):
    print(m)
    sys.stdout.flush()


def auc_de(pos, neg):
    """AUC de Mann-Whitney con la CONVENCION DEL PROYECTO: menos energia es mejor.

    Cuenta los pares en los que el positivo tiene afinidad MAS NEGATIVA que el fondo,
    que es lo que hace `auc()` en `validar_sod1_limpia.py`. La version clasica
    (proporcion de pares en los que el positivo puntua mas) sale invertida, y aqui se
    vio: daba 0,699 donde el validador da 0,301, que es su complemento exacto. Esa
    comprobacion es la que asegura que las dos copias cuentan lo mismo.
    """
    pos = np.asarray(pos, dtype=float).ravel()
    neg = np.asarray(neg, dtype=float).ravel()
    if pos.size == 0 or neg.size == 0:
        return None
    comparar = neg[None, :] - pos[:, None]      # positivo mejor = mas negativo
    return float((np.sum(comparar > 0) + 0.5 * np.sum(comparar == 0))
                 / comparar.size)


def carga_afinidad(ruta):
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            try:
                return float(l.split()[3])
            except (IndexError, ValueError):
                return None
    return None


def ligandos_de(d, modo="todos"):
    """{nombre: afinidad} de las poses de una carpeta.

    `modo` decide que se queda, y son tres casos distintos que NO se pueden mezclar:
      * `solo_act`: solo las que llevan prefijo `ACT_`, que es donde esta el positivo;
      * `sin_act`:  solo las que NO lo llevan, que es el fondo. Los sedantes con
        prefijo `ACT_` que no tienen union medida se dejan fuera de todo, que es lo
        que hace el validador;
      * `todos`:   todo, para diagnostico.
    """
    out = {}
    if not d or not os.path.isdir(d):
        return out
    for p in sorted(os.listdir(d)):
        if not p.endswith("_out.pdbqt"):
            continue
        crudo = p[: -len("_out.pdbqt")]
        es_act = crudo.startswith("ACT_")
        if modo == "solo_act" and not es_act:
            continue
        if modo == "sin_act" and es_act:
            continue
        aff = carga_afinidad(os.path.join(d, p))
        if aff is not None:
            out[crudo[4:] if es_act else crudo] = aff
    return out


def main():
    import validar_sod1_limpia as V
    import validar_diana_limpia as D

    log("INTERVALO DE CONFIANZA DEL AUC  (%d re-muestreos, semilla %d)"
        % (N_BOOT, SEMILLA))
    log("Se re-muestrea SOLO el conjunto de positivos: con 7, el error es suyo,")
    log("y el del fondo (cientos de ligandos) es despreciable al lado.")
    log("")
    log("%-7s %-12s %5s %7s %18s %9s"
        % ("diana", "fondo", "n+", "AUC", "IC 95% del AUC", "P(AUC>0,5)"))
    log("-" * 74)

    rng = np.random.default_rng(SEMILLA)
    resumen = []
    for diana in ("TDP43", "SOD1"):
        cfg = D.DIANAS[diana]
        out_salida = os.path.join(cfg["salida"], "out")
        poses2 = cfg.get("poses2") or os.path.join(
            os.path.dirname(V.LIGS_FONDO2), "out")

        truths = {r["ligand"] for r in csv.DictReader(
            open(os.path.join(BASE, "verdad_de_referencia.csv"), encoding="utf-8"))
            if r["target"] == diana and r["apto"] == "si"}
        positivos = {}
        for d in (out_salida, cfg["poses"]):
            for n, a in ligandos_de(d, modo="solo_act").items():
                if n in truths:
                    positivos[n] = a
        blando = ligandos_de(cfg["poses"], modo="sin_act")
        duro = ligandos_de(poses2, modo="sin_act")

        for etiqueta, fondo in (("blando", blando), ("duro", duro),
                                ("los dos", {**blando, **duro})):
            if not positivos or not fondo:
                continue
            p = [positivos[k] for k in sorted(positivos)]
            n = [fondo[k] for k in sorted(fondo)]
            punto = auc_de(p, n)
            p = np.asarray(p, dtype=float)
            idx = rng.integers(0, len(p), size=(N_BOOT, len(p)))
            muestras = np.array([auc_de(p[i], n) for i in idx])
            lo, hi = np.percentile(muestras, [2.5, 97.5])
            prob = float(np.mean(muestras > 0.5))
            resumen.append((diana, etiqueta, len(p), punto, lo, hi, prob))
            log("%-7s %-12s %5d %7.3f    [%.3f , %.3f] %9.3f"
                % (diana, etiqueta, len(p), punto, lo, hi, prob))

    log("")
    log("COMPROBACION CONTRA EL VALIDADOR: el punto tiene que salir igual que las")
    log("cifras publicadas. Con menos energia es mejor (misma convencion que")
    log("`auc()` del validador): TDP-43 blando 0,301 | duro 0,357 | los dos 0,332;")
    log("SOD1 blando 0,440 | duro 0,313 | los dos 0,409.")
    log("")
    log("LECTURA")
    log("-------")
    for diana, etiqueta, n, punto, lo, hi, prob in resumen:
        ancho = hi - lo
        if lo > 0.5:
            veredicto = "por encima del azar, y el margen NO toca el azar"
        elif hi < 0.5:
            veredicto = "por debajo del azar, y el margen NO toca el azar"
        else:
            veredicto = "NO DISTINGUIBLE del azar (el margen cruza 0,5)"
        log("%-7s %-12s n+ = %2d | AUC %.3f | IC [%.3f, %.3f] (ancho %.3f) | %s"
            % (diana, etiqueta, n, punto, lo, hi, ancho, veredicto))
    log("")
    log("CUANTOS POSITIVOS HACEN FALTA PARA PODER DECIR ALGO")
    log("---------------------------------------------------")
    log("El margen de un AUC.bootstrap escala como 1/sqrt(n). Para estrecharlo hay que")
    log("multiplicar el numero de positivos por (margen actual / margen deseado)^2.")
    log("Con los numeros de arriba:")
    for diana, etiqueta, n, punto, lo, hi, prob in resumen:
        semiancho = (hi - lo) / 2.0
        for objetivo, etiqueta_obj in ((0.10, "IC de +-0,10"), (0.05, "IC de +-0,05")):
            n_need = int(round(n * (semiancho / objetivo) ** 2))
            log("   %-7s %-12s semiancho %.3f -> para %s harian falta ~%d positivos"
                " (hoy %d)" % (diana, etiqueta, semiancho, etiqueta_obj,
                               n_need, n))
    log("")
    log("   Leidas al reves: para poder afirmar que el embudo separa TDP-43 del azar")
    log("   con una confianza del 95% haria falta del orden de 150 positivos MEDIDOS de")
    log("   TDP-43, y hay 7. El cribado no puede responder esa pregunta, y ningun")
    log("   cambio de funcion de puntuacion la responde.")
    log("")
    log("CONSECUENCIA PRACTICA")
    log("---------------------")
    log("Un intervalo de este anchura no deja afirmar NADA sobre la funcion de")
    log("puntuacion: por eso las doce celdas con Vinardo no se podian comparar entre si,")
    log("y por eso el 'NO PASA' repetido no demuestra que el embudo falle, demuestra que")
    log("el numero de positivos medidos es demasiado pequeno para decidirlo.")
    log("Aumentar funciones de puntuacion sin aumentar positivos medidos es mas de lo")
    log("mismo. Lo que falta es medicion.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
