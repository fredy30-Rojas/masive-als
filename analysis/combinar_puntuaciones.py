# -*- coding: utf-8 -*-
"""¿Se puede arreglar el ranking combinando energía y contactos? (17 sep 2026)

Parte de la tabla que dejó puntuacion_interaccion.py, que midió esto contra un
fondo aleatorio de 200 moléculas (no contra el top-100 por Vina, que está
sesgado):

    afinidad 0,657 | apilamiento 0,605 | apilamiento+puentes 0,598
    contactos 0,425 | contactos por átomo 0,379  ← peor que el azar

Dos lecturas de ahí:
  - Contar contactos es medir tamaño otra vez (por eso sale por debajo de 0,5).
  - Lo que sí informa es la QUÍMICA del contacto: apilarse con las aromáticas
    que apilan las bases del ARN.

Ahora se prueba si combinar ayuda, y con una regla de disciplina: **las
combinaciones no se ajustan a los controles**. Pesos iguales sobre puntuaciones
tipificadas (z), nada de buscar el peso que mejor queda con 8 moléculas. Si algo
mejora, mejora con la regla puesta de antemano.

Y se corrige el tamaño de la afinidad con un ajuste hecho sobre el FONDO (no
sobre los controles): residuo = afinidad − (a + b·átomos pesados).

Uso: python combinar_puntuaciones.py
"""
import csv
import json
import os

ANALISIS = r"C:\Users\Fredy\masive-als\analysis"
TABLA = os.path.join(ANALISIS, "rescoring_local", "puntuacion_interaccion.csv")
SALIDA = os.path.join(ANALISIS, "rescoring_local", "combinaciones_auc.json")


def media(v):
    return sum(v) / len(v) if v else 0.0


def desv(v):
    if len(v) < 2:
        return 1.0
    m = media(v)
    return (sum((x - m) ** 2 for x in v) / (len(v) - 1)) ** 0.5 or 1.0


def regresion(xs, ys):
    """Ajuste lineal simple y = a + b·x."""
    mx, my = media(xs), media(ys)
    sxx = sum((x - mx) ** 2 for x in xs) or 1e-9
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    return my - b * mx, b


def auc_positivo(ct, fo):
    """AUC con "más alto = mejor"."""
    if not ct or not fo:
        return None
    gana = 0.0
    for x in ct:
        for y in fo:
            gana += 0.5 if x == y else (1.0 if x > y else 0.0)
    return gana / (len(ct) * len(fo))


def main():
    filas = list(csv.DictReader(open(TABLA, encoding="utf-8")))
    controles = [r for r in filas if r["grupo"] == "control"]
    fondo = [r for r in filas if r["grupo"] == "fondo"]
    print("Controles: %d | fondo: %d" % (len(controles), len(fondo)))

    resumen = {}
    for caja, pre in (("caja_actual", "actual_"), ("caja_arn", "arn_")):
        # Solo las filas con todos los datos que hacen falta
        def num(r, c):
            try:
                return float(r[pre + c])
            except (KeyError, TypeError, ValueError):
                return None
        def valido(r):
            return all(num(r, c) is not None for c in ("afinidad", "atomos", "apilamiento", "puentes"))
        ct = [r for r in controles if valido(r)]
        fo = [r for r in fondo if valido(r)]
        if not ct or not fo:
            continue

        # Corrección de tamaño aprendida SOLO del fondo
        a, b = regresion([num(r, "atomos") for r in fo], [num(r, "afinidad") for r in fo])
        print("\n%s: afinidad = %.2f %+.4f·atomos (ajustado sobre el fondo, n=%d)" % (
            caja, a, b, len(fo)))

        def residual(r):
            return num(r, "afinidad") - (a + b * num(r, "atomos"))

        # Puntuaciones tipificadas con media y desviación del fondo
        def z(campo, fn):
            vals = [fn(r) for r in fo]
            m, s = media(vals), desv(vals)
            return lambda r: (fn(r) - m) / s

        z_afin = z("afinidad", lambda r: -num(r, "afinidad"))      # menos es mejor → se invierte
        z_resid = z("residual", lambda r: -residual(r))
        z_apila = z("apilamiento", lambda r: num(r, "apilamiento"))
        z_puente = z("puentes", lambda r: num(r, "puentes"))

        combinaciones = {
            "afinidad": z_afin,
            "afinidad_corregida_tamano": z_resid,
            "apilamiento": z_apila,
            "afinidad + apilamiento": lambda r: z_afin(r) + z_apila(r),
            "afinidad_corregida + apilamiento": lambda r: z_resid(r) + z_apila(r),
            "afinidad + apilamiento + puentes": lambda r: z_afin(r) + z_apila(r) + z_puente(r),
            "afinidad_corregida + apilamiento + puentes": lambda r: z_resid(r) + z_apila(r) + z_puente(r),
        }

        print("%-44s %s" % ("puntuación (pesos iguales, nada ajustado)", "AUC"))
        print("-" * 56)
        for nombre, fn in combinaciones.items():
            a_ = auc_positivo([fn(r) for r in ct], [fn(r) for r in fo])
            resumen["%s|%s" % (caja, nombre)] = round(a_, 4)
            print("%-44s %.3f" % (nombre, a_))

        # Dónde cae cada control con la mejor de las combinaciones simples
        mejor = max(combinaciones.items(), key=lambda kv: auc_positivo([kv[1](r) for r in ct],
                                                                      [kv[1](r) for r in fo]))
        print("\nCon «%s» (la mejor de esta caja), puesto de cada control entre %d moléculas:"
              % (mejor[0], len(fo) + len(ct)))
        valores_fondo = [mejor[1](r) for r in fo]
        for r in sorted(ct, key=mejor[1], reverse=True):
            mejor_que = sum(1 for v in valores_fondo if mejor[1](r) > v)
            print("   %-16s percentil %3.0f  (afinidad %s)" % (
                r["ligand"], 100.0 * mejor_que / len(valores_fondo), r[pre + "afinidad"]))

    with open(SALIDA, "w", encoding="utf-8") as f:
        json.dump(resumen, f, indent=2, ensure_ascii=False)
    print("\nGuardado: %s" % SALIDA)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
