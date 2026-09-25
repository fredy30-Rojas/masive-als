#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Los ocho de la tabla del CR, acoplados en el CR: distingue el sitio a los activos?

LA PREGUNTA
-----------
La pareja XL20/XL23 dejo claro que el acoplamiento no ordena dos compuestos cuando se
llevan 8 atomos pesados de diferencia. Queda la pregunta de al lado, que es mas util:
**los compuestos de la misma tabla que SI tienen actividad funcional (XL21 y XL23
inhiben la agregacion del LCD in vitro a 100 uM) se colocan de otra manera en el sitio
del Trp334 que los que no tienen nada reportado (XL22, XL25, XL26)?** Y de paso, donde
caen el unico con union medida (XL20) y el que empeora la muerte neuronal (XL24).

Por que la pregunta tiene sentido: los tres grupos no se separan por tamano. Los dos
inhibidores tienen 23 y 34 atomos pesados; los tres sin nada reportado, 22, 29 y 22. Si
el resultado ordenara por tamano, esta prueba lo ensenaria.

QUE HACE
--------
1. Prepara los SIETE companeros de XL20 con la receta canonica (`preparar_ligando.py`,
   la de todos) leyendo sus SMILES de `suplementario_tdp43_xl21_27.csv`. **No se
   reinventa nada**: la lectura de poses, los contactos a 4 A, la distancia al anillo
   del Trp334 y la llamada a Vina se importan del script de la pareja.
2. Los acopla con **los mismos parametros y la misma caja**: 6 modelos de RMN del 2N2C,
   3 semillas (42/2026/777), exhaustividad 8, 9 modos. Las poses de XL20 y XL23 que ya
   estan calculadas **se reutilizan del mismo directorio**, no se vuelven a correr: asi
   las ocho estan hechas exactamente con lo mismo.
3. Compara por grupos de evidencia funcional y **mira si las distribuciones se solapan
   o no**, en vez de comparar dos numeros sueltos.

QUE NO SE PUEDE CONCLUIR CON ESTO, Y VA ESCRITO EN LA SALIDA
------------------------------------------------------------
Con dos compuestos que inhiben y tres que no, **no hay test estadistico posible**: los
grupos son de dos y de tres. Lo que se puede ver es si las distribuciones se solapan o
estan separadas, y si la separacion es mayor que el ruido del propio metodo (el reparto
de un mismo compuesto entre modelos y semillas). Y no hay fondo de señuelos en el CR, asi
que una afinidad absoluta aqui no significa nada por si sola: solo cuenta la comparacion
entre estos ocho.

Uso:
    python analysis/acoplar_familia_cr.py [--exhaustividad 8] [--semillas 42,2026,777]

Salida en analysis/_cr_receptor/:
    ligands/XL2*.pdbqt              los ligandos preparados
    out/<ligando>_modelo<n>_s<s>.pdbqt   las poses (carpeta compartida con la pareja)
    familia_cr.csv                  una fila por pose, con su grupo de evidencia
    familia_cr.txt                  el resumen legible
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import statistics
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)
sys.path.insert(0, os.path.dirname(BASE))
import acoplar_xl20_xl23_cr as P  # noqa: E402  (la maquinaria de la pareja, sin copiarla)
from preparar_ligando import escribir as escribir_ligando  # noqa: E402

SALIDA = P.SALIDA
MODELOS = P.MODELOS
LIGANDS = P.LIGANDS
OUT = P.OUT
CSV_PAREJA = os.path.join(BASE, "suplementario_tdp43_xl21_27.csv")

# Los ocho de la tabla del articulo, con de donde sale su SMILES y su evidencia
# funcional. La evidencia es LITERAL del CSV: no se resume ni se interpreta aqui.
# (el texto del grupo sale de la columna `tipo_ensayo` del propio fichero)
FUENTES = {
    "XL20": ("controles_tdp43_xl20.csv", "union medida"),
    "XL21": (CSV_PAREJA, "inhibe la agregacion"),
    "XL23": (CSV_PAREJA, "inhibe la agregacion"),
    "XL22": (CSV_PAREJA, "sin nada reportado"),
    "XL25": (CSV_PAREJA, "sin nada reportado"),
    "XL26": (CSV_PAREJA, "sin nada reportado"),
    "XL24": (CSV_PAREJA, "empeora"),
    "XL27": (CSV_PAREJA, "neuroprotege solo a 100 uM"),
}

# Los grupos que se comparan. Se declaran aqui y no se deducen del texto para que se
# vea de un vistazo quien esta en cada uno y se pueda discutir sin mirar el codigo.
GRUPOS = {
    "union medida": ["XL20"],
    "inhibe la agregacion (100 uM)": ["XL21", "XL23"],
    "sin nada reportado": ["XL22", "XL25", "XL26"],
    "otros": ["XL24", "XL27"],
}


def log(m):
    print(m, flush=True)


def leer(archivo, ligando):
    """(smiles, evidencia) de un ligando, de su CSV."""
    with open(os.path.join(BASE, archivo), encoding="utf-8-sig", newline="") as f:
        for fila in csv.DictReader(f):
            if (fila.get("ligand") or "").strip().upper() == ligando:
                return (fila.get("smiles") or "").strip(), (fila.get("tipo_ensayo") or "")
    raise SystemExit("no encontrado %s en %s" % (ligando, archivo))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustividad", type=int, default=8)
    ap.add_argument("--semillas", default="42,2026,777")
    args = ap.parse_args()
    semillas = [int(s) for s in args.semillas.split(",") if s.strip()]

    os.makedirs(LIGANDS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(SALIDA, "caja.json"), encoding="utf-8") as f:
        caja = json.load(f)
    modelos = sorted(caja["modelos"], key=int)

    # --- los ocho ligandos, con la receta canonica ---
    evidencia, pesados, pdbqt = {}, {}, {}
    for nombre, (archivo, grupo) in FUENTES.items():
        smi, texto = leer(archivo, nombre)
        ruta = escribir_ligando(nombre, smi, LIGANDS)
        if ruta is None:
            log("%s NO se pudo preparar" % nombre)
            return 1
        pdbqt[nombre] = ruta
        evidencia[nombre] = (grupo, texto)
        pesados[nombre] = len(P.atomos(ruta))
    log("Los %d compuestos de la tabla del CR, preparados con la receta canonica:"
        % len(pdbqt))
    for nombre in FUENTES:
        log("  %-5s %-30s %2d atomos pesados | %s"
            % (nombre, FUENTES[nombre][1], pesados[nombre], evidencia[nombre][1][:70]))
    log("")

    # cuantas tandas hay que correr de verdad (las de XL20 y XL23 ya estan)
    pendientes = 0
    for nombre in pdbqt:
        for clave in modelos:
            for s in semillas:
                f = os.path.join(OUT, "%s_modelo%s_s%d.pdbqt" % (nombre, clave, s))
                if not (os.path.exists(f) and os.path.getsize(f) > 100):
                    pendientes += 1
    log("tandas de Vina por correr: %d de %d (las de la pareja ya estaban hechas)"
        % (pendientes, len(pdbqt) * len(modelos) * len(semillas)))
    log("")

    filas = []
    for nombre, ligando in pdbqt.items():
        for clave in modelos:
            info = caja["modelos"][clave]
            receptor = os.path.join(BASE, info["pdbqt"])
            rec = P.atomos(os.path.join(MODELOS, "cr_modelo%s.pdb" % clave))
            anillo = np.array([xyz for (res, nom, el, xyz) in rec
                               if res == P.TRP and nom in P.ANILLO_TRP])
            centro = anillo.mean(axis=0) if len(anillo) else None
            rec_xyz = np.array([a[3] for a in rec])
            for semilla in semillas:
                out = os.path.join(OUT, "%s_modelo%s_s%d.pdbqt" % (nombre, clave, semilla))
                if not (os.path.exists(out) and os.path.getsize(out) > 100):
                    r = P.acoplar(ligando, receptor, info["centro"][0], info["centro"][1],
                                  info["centro"][2], info["tamano"], args.exhaustividad,
                                  semilla, out)
                    if not (os.path.exists(out) and os.path.getsize(out) > 100):
                        log("  %s modelo %s semilla %d: FALLO Vina: %s"
                            % (nombre, clave, semilla, (r.stderr or r.stdout)[-160:]))
                        continue
                    log("  %-5s modelo %-2s semilla %-5d ok" % (nombre, clave, semilla))
                for modo, afinidad, lig in P.poses_del_pdbqt(out):
                    L = np.array([a[1] for a in lig])
                    d = np.sqrt(((L[:, None, :] - rec_xyz[None, :, :]) ** 2).sum(-1))
                    residuos = sorted({rec[i][0] for i in np.nonzero(d.min(axis=0) < P.CORTE)[0]})
                    en_contacto = int((d.min(axis=1) < P.CORTE).sum())
                    dist_anillo = (float(np.linalg.norm(
                        L[:, None, :] - anillo[None, :, :], axis=-1).min())
                        if len(anillo) else None)
                    dist_centro = (float(np.linalg.norm(L - centro, axis=1).min())
                                   if centro is not None else None)
                    filas.append({
                        "ligando": nombre, "grupo": evidencia[nombre][0],
                        "evidencia_funcional": evidencia[nombre][1],
                        "atomos_pesados": pesados[nombre],
                        "modelo": clave, "semilla": semilla, "modo": modo,
                        "afinidad": afinidad,
                        "contactos": len(residuos),
                        "atomos_del_ligando_en_contacto": en_contacto,
                        "fraccion_ligando_en_contacto": round(en_contacto / float(len(L)), 2),
                        "residuos": " ".join(str(r) for r in residuos),
                        "toca_trp334": P.TRP in residuos,
                        "dist_min_al_anillo_trp334": (round(dist_centro, 2)
                                                      if dist_centro is not None else None),
                        "dist_min_a_un_atomo_del_anillo": (round(dist_anillo, 2)
                                                           if dist_anillo is not None else None),
                    })

    campos = list(filas[0].keys())
    with open(os.path.join(SALIDA, "familia_cr.csv"), "w", encoding="utf-8",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(filas)

    cabecera = [
        "Los ocho compuestos de la tabla del CR, acoplados en el CR",
        "receptor %s (%s), %d modelos de RMN, caja de %d A sobre el anillo del Trp334"
        " de cada modelo" % (caja["pdb"], caja["residuos"], len(modelos),
                             caja["modelos"][modelos[0]]["tamano"]),
        "exhaustividad %d, 9 modos, semillas %s | contactos a %.1f A"
        % (args.exhaustividad, ", ".join(str(s) for s in semillas), P.CORTE),
        "",
    ]
    texto = informe(filas, modelos, semillas, pesados, cabecera)
    with open(os.path.join(SALIDA, "familia_cr.txt"), "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    log("")
    log(texto)
    return 0


def afinidad(f):
    """La afinidad de una fila, como numero y no como texto.

    OJO CON ESTO, que ya costo un susto: las filas que vienen de un CSV traen la
    afinidad como TEXTO, y `min` sobre textos no compara numeros, compara letras: entre
    "-3.901" y "-4.865" elige "-4.865" (porque '4' < '5'), o sea la PEOR. Con `filas`
    hechas en memoria (float) el error no aparece, y por eso se cuela solo cuando el
    mismo codigo se usa para leer resultados ya guardados. Aqui se convierte siempre.
    """
    try:
        return float(f["afinidad"])
    except (TypeError, ValueError):
        return float("inf")


def mediana_por_compuesto(filas, modelos, semillas):
    """{compuesto: [afinidad de la mejor pose de cada (modelo, semilla)], ordenadas}."""
    out = {}
    for nombre in {f["ligando"] for f in filas}:
        v = []
        for clave in modelos:
            for semilla in semillas:
                suyas = [f for f in filas if f["ligando"] == nombre
                         and f["modelo"] == clave and str(f["semilla"]) == str(semilla)
                         and f["afinidad"] not in ("", None)]
                if suyas:
                    v.append(afinidad(min(suyas, key=afinidad)))
        out[nombre] = sorted(v)
    return out


def solape(a, b):
    """(se solapan, peor de a, mejor de b) de dos listas de afinidades."""
    return (not (max(a) < min(b) or max(b) < min(a))), max(a), min(b)


def residuos_de_tamano(medianas, pesados):
    """(pendiente, {compuesto: residuo}) tras quitar el efecto del numero de atomos.

    Vina premia los contactos, asi que parte de la afinidad es el tamano. Con ocho
    compuestos esto es descriptivo y no un modelo; el informe lo dice asi.
    """
    xs = np.array([pesados[n] for n in medianas], dtype=float)
    ys = np.array([medianas[n] for n in medianas])
    pendiente, corte = np.polyfit(xs, ys, 1)
    return float(pendiente), {n: float(ys[i] - (pendiente * xs[i] + corte))
                              for i, n in enumerate(medianas)}


def informe(filas, modelos, semillas, pesados, cabecera):
    """El analisis de la familia, en UN solo sitio: lo usan el rigido y el flexible.

    `cabecera` son las lineas de arriba (que receptor, con que caja, cuantos modelos y
    que semillas), que el script del receptor flexible cambia por las suyas para que el
    informe diga de donde sale cada numero.
    """
    lineas = list(cabecera)
    valores = mediana_por_compuesto(filas, modelos, semillas)
    lineas.append("  %-5s %-31s %-4s %-8s %-8s %-8s %-8s %-10s"
                  % ("lig", "evidencia funcional", "at.", "mejor", "mediana", "peor",
                     "kcal/at.", "sobre Trp334"))
    for nombre in FUENTES:
        v = valores[nombre]
        sobre = sum(1 for f in filas if f["ligando"] == nombre
                    and (f["dist_min_a_un_atomo_del_anillo"] or 99) < P.CORTE)
        total = sum(1 for f in filas if f["ligando"] == nombre)
        lineas.append("  %-5s %-31s %-4d %-8.2f %-8.2f %-8.2f %-8.3f %-10s"
                      % (nombre, FUENTES[nombre][1], pesados[nombre], v[0],
                         statistics.median(v), v[-1],
                         statistics.median(v) / pesados[nombre],
                         "%d de %d poses" % (sobre, total)))

    lineas.append("")
    lineas.append("  POR GRUPOS (mejor y mediana de las %d combinaciones de modelo y"
                  " semilla de cada compuesto)" % (len(modelos) * len(semillas)))
    grupos = {}
    for etiqueta, miembros in GRUPOS.items():
        todos = [x for n in miembros for x in valores[n]]
        grupos[etiqueta] = todos
        lineas.append("    %-32s n=%d | de %.2f a %.2f | mediana %.2f"
                      % (etiqueta, len(miembros), min(todos), max(todos),
                         statistics.median(todos)))
    lineas.append("")
    lineas.append("  SE SOLAPAN O ESTAN SEPARADAS (lo que de verdad se puede leer):")
    a = grupos["inhibe la agregacion (100 uM)"]
    b = grupos["sin nada reportado"]
    solapan, peor_a, mejor_b = solape(a, b)
    lineas.append("    inhibidores (XL21, XL23) frente a sin nada reportado (XL22, XL25, XL26):")
    lineas.append("      peor de los inhibidores: %.2f | mejor de los otros: %.2f"
                  % (peor_a, mejor_b))
    lineas.append("      -> %s" % ("SE SOLAPAN: el sitio no los separa"
                                   if solapan else "sin solape: el sitio los separa"))
    lineas.append("      (XL20, el unico con union medida: de %.2f a %.2f; XL24, que empeora la"
                  " muerte neuronal: de %.2f a %.2f)"
                  % (min(valores["XL20"]), max(valores["XL20"]),
                     min(valores["XL24"]), max(valores["XL24"])))
    lineas.append("")
    # --- el tamano, apartado: quien queda mejor de lo que le toca por tamano ---
    # Vina premia los contactos, asi que parte de la afinidad es el numero de atomos.
    # Se ajusta una recta (afinidad mediana frente a atomos pesados) con los ocho y se
    # miran los residuos: quien esta por debajo de la recta se coloca mejor de lo que su
    # tamano predice. Con ocho compuestos esto es descriptivo, no un modelo, y asi se dice.
    pendiente, residuos = residuos_de_tamano(
        {n: statistics.median(valores[n]) for n in valores}, pesados)
    lineas.append("")
    lineas.append("  QUITANDO EL EFECTO DEL TAMANO (recta de afinidad mediana frente a"
                  " atomos pesados:")
    lineas.append("  %.3f kcal/mol por atomo). El residuo es lo que le sobra o le falta a cada"
                  % pendiente)
    lineas.append("  uno respecto a lo que su tamano predice; negativo = mejor de lo que le toca.")
    lineas.append("  %-5s %-31s %-4s %-9s %s"
                  % ("lig", "evidencia funcional", "at.", "residuo", "lectura"))
    for nombre in sorted(residuos, key=lambda n: residuos[n]):
        lineas.append("  %-5s %-31s %-4d %-9.3f %s"
                      % (nombre, FUENTES[nombre][1], pesados[nombre],
                         residuos[nombre], "mejor de lo que le toca"
                         if residuos[nombre] < 0 else "peor de lo que le toca"))
    arriba = [n for n in sorted(residuos, key=lambda n: residuos[n])[:3]]
    abajo = [n for n in sorted(residuos, key=lambda n: residuos[n])[-3:]]
    lineas.append("    los tres mejores por residuo: %s" % ", ".join(arriba))
    lineas.append("    los tres peores por residuo:  %s" % ", ".join(abajo))
    # Cuanto de todo esto es ruido del propio metodo: el reparto de UN MISMO compuesto
    # entre modelos y semillas. Es contra esto contra lo que hay que comparar cualquier
    # diferencia de la tabla de arriba, no contra cero.
    desviaciones = []
    for n in valores:
        v = valores[n]
        if len(v) > 1:
            desviaciones.append(statistics.stdev(v))
    ruido = statistics.median(desviaciones)
    lineas.append("    (de cuanto: el reparto de un mismo compuesto entre modelos y semillas"
                  " tiene %.2f kcal/mol" % ruido)
    lineas.append("     de desviacion tipica; el residuo va de %.2f a %.2f, del mismo orden)"
                  % (min(residuos.values()), max(residuos.values())))
    lineas.append("")
    lineas.append("  LO QUE ESTO NO PUEDE DECIR:")
    lineas.append("    Dos compuestos que inhiben y tres sin nada reportado no son un test")
    lineas.append("    estadistico: los grupos son de dos y de tres, y el ruido del propio")
    lineas.append("    metodo (el reparto de un mismo compuesto entre modelos y semillas)")
    lineas.append("    es de decimas de kcal/mol, igual que las diferencias que se ven.")
    lineas.append("    Y no hay fondo de señuelos en el CR: aqui solo cuenta la comparacion")
    lineas.append("    entre estos ocho, nunca una afinidad absoluta.")

    return "\n".join(lineas)


if __name__ == "__main__":
    sys.exit(main())
