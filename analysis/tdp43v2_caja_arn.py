# -*- coding: utf-8 -*-
"""Elige la caja nueva de TDP43_v2 mirando el sitio de unión de ARN, y la compara
con la que usa el cribado en igualdad de condiciones.

Lo que se aprendió en tdp43v2_sitio_arn.py:
  - El receptor del cribado y la estructura 4BS2 (TDP-43 + ARN UG-rico) tienen
    EXACTAMENTE el mismo sistema de coordenadas (desviación media del CA 0,00 Å).
    Así que el sitio de unión de ARN se puede leer de la estructura sin alinear
    nada.
  - 39 residuos tocan al ARN, repartidos por 35,6 x 42,0 x 31,6 Å (el ARN cruza
    los dos dominios RRM).
  - La caja del cribado NO está en otro bolsillo: está a 7,6 Å del centro de esa
    interfaz y contiene el 49 % de sus átomos. El problema no es que mire a otro
    sitio, es que mira a medias.

Aquí se prueban varios centros posibles, todos con el MISMO tamaño que la caja
vieja (26 Å), para poder comparar enterramiento sin trampa: cuanto más hueco
encierra una caja, más premia el acoplamiento al que rellena cavidad en vez de
al que se apoya en una superficie.

Salida: por pantalla y un JSON con la caja elegida.
"""
import json
import math
import os

from tdp43v2_sitio_arn import CAJA_VIEJA, CORTE_CONTACTO, dist, leer_modelo_1, leer_receptor

SALIDA = os.path.join(r"C:\Users\Fredy\masive-als\analysis\rescoring_local", "caja_tdp43v2_arn.json")
LADO = 26.0     # el mismo que la caja del cribado
# Plataformas aromáticas que apilan las bases del ARN (RNP1 y vecinas).
AROMATICOS = [147, 149, 194, 221, 229, 231]
# Zona clásica de reconocimiento del ARN en RRM1: RNP2 y RNP1.
RRM1_RNP = [109, 113, 147, 149]


def centro_de(atomos, residuos=None):
    sel = [a for a in atomos if residuos is None or a[0] in residuos]
    if not sel:
        return None
    return tuple(sum(a[i] for a in sel) / len(sel) for i in (3, 4, 5))


def dentro(p, centro, lado):
    return all(abs(p[i + 3] - centro[i]) <= lado / 2 for i in range(3))


def main():
    prot, rna = leer_modelo_1()
    recept = leer_receptor()

    # Interfaz del ARN (mismos números que el script anterior)
    tocan = set()
    for p in prot:
        for r in rna:
            if abs(p[3] - r[3]) > CORTE_CONTACTO:
                continue
            if dist((p[3], p[4], p[5]), (r[3], r[4], r[5])) <= CORTE_CONTACTO:
                tocan.add((p[0], p[1]))
                break
    interfaz = [p for p in prot if (p[0], p[1]) in tocan]

    candidatos = {
        "interfaz_arn_completa": centro_de(interfaz),
        "plataformas_aromaticas": centro_de(prot, AROMATICOS),
        "rrm1_rnp": centro_de(prot, RRM1_RNP),
        "caja_actual_del_cribado": CAJA_VIEJA[0],
    }

    print("Todas las cajas con lado %.0f Å, para poder comparar sin trampa.\n" % LADO)
    print("%-26s %-24s %8s %10s %10s" % ("candidato", "centro", "interfaz", "atomos", "dist. caja"))
    print("%-26s %-24s %8s %10s %10s" % ("", "", "de ARN", "receptor", "actual"))
    print("-" * 84)

    filas = {}
    for nombre, c in candidatos.items():
        if c is None:
            continue
        n_if = sum(1 for a in interfaz if dentro(a, c, LADO))
        n_rec = sum(1 for a in recept if dentro(a, c, LADO))
        d = dist(c, CAJA_VIEJA[0])
        filas[nombre] = {"centro": [round(x, 3) for x in c], "interfaz_dentro": n_if,
                         "atomos_receptor": n_rec, "dist_a_caja_actual": round(d, 2)}
        print("%-26s %-24s %6d/%-3d %10d %9.1f" % (
            nombre, "(%.1f, %.1f, %.1f)" % c, n_if, len(interfaz), n_rec, d))

    total_if = len(interfaz)
    print("\nLectura: el candidato que cubre más interfaz de ARN encerrando MENOS")
    print("receptor es el que de verdad mide unión y no relleno de cavidad.")

    # Elección: de las cajas que encierran MENOS receptor que la actual (o sea,
    # menos relleno de cavidad), la que más interfaz de ARN cubre. Se exige lo
    # primero porque el enterramiento es justo el artefacto que se quiere
    # esquivar; y entre las que cumplen, se premia cubrir el sitio real.
    aforo_actual = filas["caja_actual_del_cribado"]["atomos_receptor"]
    validos = {k: v for k, v in filas.items()
               if k != "caja_actual_del_cribado" and v["atomos_receptor"] <= aforo_actual}
    elegido = max(validos, key=lambda k: validos[k]["interfaz_dentro"]) if validos else None
    if elegido:
        print("\nElegida: %s  centro (%.2f, %.2f, %.2f)  lado %.0f Å" % (
            elegido, *filas[elegido]["centro"], LADO))

    datos = {
        "origen": "4BS2 modelo 1 (TDP-43 tandem RRM + ARN UG-rico); mismo sistema de coordenadas que el receptor",
        "lado": LADO,
        "candidatos": filas,
        "elegida": elegido,
        "caja_actual": {"centro": list(CAJA_VIEJA[0]), "lado": CAJA_VIEJA[1]},
        "criterio": ("De las cajas que encierran menos receptor que la actual (menos relleno de "
                     "cavidad), la que mas atomos del sitio de union de ARN cubre"),
    }
    if elegido:
        datos["centro_elegido"] = filas[elegido]["centro"]
    with open(SALIDA, "w", encoding="utf-8") as f:
        json.dump(datos, f, indent=2, ensure_ascii=False)
    print("\nGuardado: %s" % SALIDA)
    return 0 if elegido else 1


if __name__ == "__main__":
    raise SystemExit(main())
