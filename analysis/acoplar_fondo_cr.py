#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fondo de señuelos emparejados en propiedades en el CR: ¿destaca alguno de los ocho?

LA PREGUNTA
-----------
Todo lo que se ha dicho del CR hasta ahora compara los ocho compuestos ENTRE ELLOS. Eso
sirve para ver quien va delante de quien, pero no para saber si alguno se une bien: sin
un fondo con el que comparar, una afinidad de -6 kcal/mol no significa nada por si sola.
Este script pone ese fondo: **señuelos emparejados en propiedades con cada uno de los
ocho** (mismo tamaño, misma hidrofobia, misma polaridad, mismos anillos, misma
flexibilidad) y **distintos en el esqueleto**, que es la convencion que ya usa el
proyecto para los señuelos emparejados.

Si ningún compuesto de los ocho destaca sobre su propio fondo, entonces el acoplamiento
en este sitio no mide union: mide cuánto abulta la molécula. Y eso zanja la discusion.

QUE HACE
--------
1. **Elige los señuelos** de las librerías que ya estan en el repositorio (no se baja
   nada): `compounds/decoys_library.smi` (6.461), `compounds/libreria_extra_chembl.csv`,
   `compounds/full_fda_library.csv` y `compounds/batch2_chembl.csv`. Para cada uno de los
   ocho busca moléculas con los descriptores emparejados y **Tanimoto ECFP4 < 0,35**
   (mismo perfil, distinta quimica), una por esqueleto de Murcko. Si con las ventanas
   estrechas no hay bastantes, se ensanchan **en pasos declarados** y se guarda que
   ventana hizo falta para cada compuesto: eso tambien es un resultado.
2. **Acopla los señuelos Y vuelve a acoplar los ocho** con `--cpu 1` en el mismo proceso
   y con los mismos ajustes (6 modelos, semilla 42, exhaustividad 8, 9 modos). Se
   reacoplan los ocho a proposito: el numero de nucleos de Vina cambia como se reparte
   la busqueda, asi que comparar contra sus corridas de antes (que usaron los 20 nucleos)
   no seria limpio. Mismo proceso para todos, o no vale.
3. **Compara**: posicion de cada uno de los ocho dentro del fondo (percentil), cuantos
   señuelos lo superan, y el ranking completo de la tabla.

LO QUE ESTE FONDO NO ES
-----------------------
No es un fondo de inactivos: son moléculas parecidas en propiedades que **no se sabe** si
se unen al CR (de hecho, de ninguna de ellas hay dato). Por eso lo que se lee no es «se
une o no se une», sino **si el sitio distingue a unos de otros**; y si los ocho caen
dentro del fondo, la respuesta es que no distingue.

Uso:
    python analysis/acoplar_fondo_cr.py [--por-ligando 8] [--hilos 6]

Salida en analysis/_cr_receptor/:
    fondo_moleculas.csv    los señuelos elegidos, con su compuesto de referencia y la ventana usada
    fondo_out/             las poses (los ocho y los señuelos, mismo proceso)
    fondo_cr.csv           una fila por compuesto: mejor, mediana, peor y percentil
    fondo_cr.txt           el informe con percentiles y ranking
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import statistics
import sys
from concurrent.futures import ProcessPoolExecutor

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
sys.path.insert(0, BASE)
sys.path.insert(0, RAIZ)
import acoplar_xl20_xl23_cr as P  # noqa: E402
import acoplar_familia_cr as F  # noqa: E402
from preparar_ligando import escribir as escribir_ligando  # noqa: E402

SALIDA = P.SALIDA
FONDO = os.path.join(SALIDA, "fondo_out")
FONDO_LIG = os.path.join(SALIDA, "fondo_ligands")
POOL = [
    os.path.join(RAIZ, "compounds", "decoys_library.smi"),
    os.path.join(RAIZ, "compounds", "libreria_extra_chembl.csv"),
    os.path.join(RAIZ, "compounds", "full_fda_library.csv"),
    os.path.join(RAIZ, "compounds", "batch2_chembl.csv"),
    os.path.join(BASE, "_rbind", "rbind_fondo_sm.csv"),
]

# Ventanas de emparejamiento. La primera es la estrecha; las siguientes se aplican solo
# si con la anterior no hay señuelos suficientes, y queda escrito cual hizo falta.
VENTANAS = [
    dict(pesados=2, logp=0.7, tpsa=15.0, donantes=1, aceptores=2, rotables=2, anillos=1),
    dict(pesados=3, logp=1.0, tpsa=25.0, donantes=2, aceptores=3, rotables=3, anillos=2),
    dict(pesados=4, logp=1.5, tpsa=40.0, donantes=3, aceptores=4, rotables=4, anillos=2),
]
TANIMOTO_MAX = 0.35        # parecido maximo al compuesto de referencia
SEMILLA = 42               # la misma semilla para todos: los ocho y el fondo
EXHAUSTIVIDAD = 8
MODELOS = [str(i) for i in range(1, 7)]


def log(m):
    print(m, flush=True)


def limpiar(nombre):
    """Nombre de fichero sin espacios ni simbolos (el fondo trae nombres con parentesis)."""
    return re.sub(r"[^A-Za-z0-9_.-]", "_", (nombre or "senuelo").strip())


def propiedades(mol):
    return {
        "pesados": mol.GetNumHeavyAtoms(),
        "logp": Descriptors.MolLogP(mol),
        "tpsa": Descriptors.TPSA(mol),
        "donantes": Descriptors.NumHDonors(mol),
        "aceptores": Descriptors.NumHAcceptors(mol),
        "rotables": Descriptors.NumRotatableBonds(mol),
        "anillos": Descriptors.RingCount(mol),
        "murcko": MurckoScaffold.MurckoScaffoldSmiles(mol=mol),
    }


def leer_pool():
    """(nombre, smiles) de todas las librerias locales."""
    moleculas = {}
    for ruta in POOL:
        if not os.path.exists(ruta):
            continue
        if ruta.endswith(".smi"):
            for l in open(ruta, encoding="utf-8", errors="ignore"):
                partes = l.split()
                if len(partes) >= 2:
                    moleculas.setdefault(partes[1], partes[0])
            continue
        with open(ruta, encoding="utf-8-sig", newline="") as f:
            for fila in csv.DictReader(f):
                nombre = (fila.get("chembl_id") or fila.get("name")
                          or fila.get("nombre") or fila.get("rbind_id") or "")
                smi = (fila.get("smiles") or "").strip()
                if smi:
                    moleculas.setdefault(nombre or smi, smi)
    return moleculas


def cabe(props, ref, ventana):
    return (abs(props["pesados"] - ref["pesados"]) <= ventana["pesados"]
            and abs(props["logp"] - ref["logp"]) <= ventana["logp"]
            and abs(props["tpsa"] - ref["tpsa"]) <= ventana["tpsa"]
            and abs(props["donantes"] - ref["donantes"]) <= ventana["donantes"]
            and abs(props["aceptores"] - ref["aceptores"]) <= ventana["aceptores"]
            and abs(props["rotables"] - ref["rotables"]) <= ventana["rotables"]
            and abs(props["anillos"] - ref["anillos"]) <= ventana["anillos"])


def elegir(objetivo, pool, cuantos, usados):
    """(señuelos, ventana_usada) emparejados con un compuesto y distintos de quimica."""
    ref_mol = Chem.MolFromSmiles(objetivo["smiles"])
    if ref_mol is None:
        return [], None
    ref = propiedades(ref_mol)
    fp_ref = AllChem.GetMorganFingerprintAsBitVect(ref_mol, 2, 2048)
    for ventana in VENTANAS:
        elegidos, esqueletos = [], set()
        for nombre, smi in pool.items():
            if nombre in usados:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None or mol.GetNumHeavyAtoms() < 12 or mol.GetNumHeavyAtoms() > 60:
                continue
            # Una sola molecula: las librerias traen sales, y si se elige una sal, los
            # descriptores emparejados serian los del conjunto y no los del ligando que
            # de verdad se acopla. Se descartan aqui y no se "quitan sales" despues, para
            # que el emparejamiento sea con la especie que se acopla.
            if len(Chem.GetMolFrags(mol)) != 1:
                continue
            props = propiedades(mol)
            if props["murcko"] in esqueletos or not props["murcko"]:
                continue
            if not cabe(props, ref, ventana):
                continue
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
            if DataStructs.TanimotoSimilarity(fp_ref, fp) > TANIMOTO_MAX:
                continue
            elegidos.append({"nombre": nombre, "smiles": smi, "props": props,
                             "tanimoto": round(DataStructs.TanimotoSimilarity(fp_ref, fp), 3)})
            esqueletos.add(props["murcko"])
            if len(elegidos) >= cuantos:
                break
        if len(elegidos) >= cuantos:
            return elegidos, ventana
    return elegidos, VENTANAS[-1]


def acoplar_uno(tarea):
    """Una tanda de Vina: (ligando.pdbqt, receptor.pdbqt, caja, salida)."""
    ligando, receptor, caja, salida = tarea
    if os.path.exists(salida) and os.path.getsize(salida) > 100:
        return os.path.basename(salida), "ya estaba"
    import subprocess
    cmd = [P.R.VINA, "--receptor", receptor, "--ligand", ligando,
           "--center_x", "%.3f" % caja[0], "--center_y", "%.3f" % caja[1],
           "--center_z", "%.3f" % caja[2],
           "--size_x", "22", "--size_y", "22", "--size_z", "22",
           "--exhaustiveness", str(EXHAUSTIVIDAD), "--num_modes", "9",
           "--cpu", "1", "--seed", str(SEMILLA), "--out", salida]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    if not (os.path.exists(salida) and os.path.getsize(salida) > 100):
        return os.path.basename(salida), ("FALLO: " + (r.stderr or r.stdout)[-120:])
    return os.path.basename(salida), "ok"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--por-ligando", type=int, default=8,
                    help="señuelos emparejados por cada uno de los ocho")
    ap.add_argument("--hilos", type=int, default=6)
    ap.add_argument("--solo-elegir", action="store_true",
                    help="elegir los señuelos y parar (para mirarlos antes de acoplar)")
    args = ap.parse_args()

    os.makedirs(FONDO, exist_ok=True)
    os.makedirs(FONDO_LIG, exist_ok=True)
    pool = leer_pool()
    log("Librerias locales: %d moleculas con SMILES" % len(pool))

    # --- 1. los ocho y su fondo ---
    objetivos = {}
    for nombre, (archivo, grupo) in F.FUENTES.items():
        smi, texto = F.leer(archivo, nombre)
        objetivos[nombre] = {"smiles": smi, "grupo": grupo, "evidencia": texto}

    usados, filas = set(), []
    for nombre, obj in objetivos.items():
        elegidos, ventana = elegir(obj, pool, args.por_ligando, usados)
        for e in elegidos:
            usados.add(e["nombre"])
            filas.append({"objetivo": nombre, "señuelo": e["nombre"], "smiles": e["smiles"],
                          "pesados": e["props"]["pesados"],
                          "logp": round(e["props"]["logp"], 2),
                          "tpsa": round(e["props"]["tpsa"], 1),
                          "rotables": e["props"]["rotables"],
                          "anillos": e["props"]["anillos"],
                          "murcko": e["props"]["murcko"],
                          "tanimoto_ecfp4": e["tanimoto"],
                          "ventana": VENTANAS.index(ventana) + 1,
                          "ventana_detalle": "pesados+-%d logp+-%.1f tpsa+-%.0f"
                          % (ventana["pesados"], ventana["logp"], ventana["tpsa"])})
        ref = propiedades(Chem.MolFromSmiles(obj["smiles"]))
        log("  %-5s %2d señuelos (ventana %d: %s) | referencia: %d atomos, logP %.2f,"
            " TPSA %.0f"
            % (nombre, len(elegidos), VENTANAS.index(ventana) + 1,
               filas[-1]["ventana_detalle"] if filas else "-", ref["pesados"],
               ref["logp"], ref["tpsa"]))
    campos = list(filas[0].keys())
    with open(os.path.join(SALIDA, "fondo_moleculas.csv"), "w", encoding="utf-8",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(filas)
    log("  total de señuelos: %d (%s)" % (len(filas), "analysis/_cr_receptor/fondo_moleculas.csv"))
    if args.solo_elegir:
        return 0

    # --- 2. preparar y acoplar: los ocho Y su fondo, en el mismo proceso ---
    with open(os.path.join(SALIDA, "caja.json"), encoding="utf-8") as f:
        caja = json.load(f)

    etiquetas = {}   # etiqueta -> (referencia, es_señuelo, ventana)
    tareas = []
    for nombre, obj in objetivos.items():
        etiqueta = "CR_" + nombre          # los ocho, otra vez y con el mismo proceso
        ruta = escribir_ligando(etiqueta, obj["smiles"], FONDO_LIG)
        if ruta:
            etiquetas[etiqueta] = (nombre, False, None)
            for clave in MODELOS:
                info = caja["modelos"][clave]
                tareas.append((ruta, os.path.join(BASE, info["pdbqt"]), info["centro"],
                               os.path.join(FONDO, "%s_modelo%s.pdbqt" % (etiqueta, clave))))
    for fila in filas:
        etiqueta = "D_" + limpiar(fila["señuelo"])
        ruta = escribir_ligando(etiqueta, fila["smiles"], FONDO_LIG)
        if not ruta:
            continue
        etiquetas[etiqueta] = (fila["objetivo"], True, fila["ventana"])
        for clave in MODELOS:
            info = caja["modelos"][clave]
            tareas.append((ruta, os.path.join(BASE, info["pdbqt"]), info["centro"],
                           os.path.join(FONDO, "%s_modelo%s.pdbqt" % (etiqueta, clave))))

    log("")
    log("tandas de Vina: %d (los %d compuestos y %d señuelos, %d modelos cada uno)"
        % (len(tareas), len(objetivos), len(filas), len(MODELOS)))
    hechas = 0
    with ProcessPoolExecutor(max_workers=args.hilos) as ex:
        for salida, estado in ex.map(acoplar_uno, tareas):
            hechas += 1
            if estado != "ya estaba" and hechas % 25 == 0:
                log("  %d de %d" % (hechas, len(tareas)))
            if estado.startswith("FALLO"):
                log("  %s %s" % (salida, estado))

    # --- 3. leer y comparar ---
    valores = {}   # etiqueta -> [afinidad por modelo]
    for etiqueta, (referencia, es_senuelo, ventana) in etiquetas.items():
        v = []
        for clave in MODELOS:
            ruta = os.path.join(FONDO, "%s_modelo%s.pdbqt" % (etiqueta, clave))
            if not os.path.exists(ruta):
                continue
            poses = P.poses_del_pdbqt(ruta)
            afs = [a for (_, a, _) in poses if a is not None]
            if afs:
                v.append(min(afs))
        if v:
            valores[etiqueta] = v

    fondo_vals = sorted(a for et in valores if et.startswith("D_")
                        for a in valores[et])

    lineas = []
    lineas.append("Fondo de señuelos emparejados en propiedades en el CR")
    lineas.append("receptor 2N2C (%s), 6 modelos de RMN, caja de 22 A, semilla %d,"
                  " exhaustividad %d, --cpu 1 en todos"
                  % (caja["residuos"], SEMILLA, EXHAUSTIVIDAD))
    lineas.append("señuelos: %d emparejados con los ocho en tamaño, logP, TPSA, donantes,"
                  " aceptores, rotables y anillos," % len(filas))
    lineas.append("         con Tanimoto ECFP4 < %.2f al de referencia y uno por esqueleto"
                  " de Murcko" % TANIMOTO_MAX)
    lineas.append("")

    def percentil(v):
        """% de valores del fondo PEORES que el mio (mas negativo = mejor)."""
        if not fondo_vals:
            return None
        return 100.0 * sum(1 for x in fondo_vals if x > v) / len(fondo_vals)

    lineas.append("  LOS OCHO DENTRO DE SU PROPIO FONDO (mejor y mediana de los 6 modelos)")
    lineas.append("  %-5s %-8s %-8s %-10s %-10s %-9s %s"
                  % ("lig", "mejor", "mediana", "percentil", "percentil", "señuelos", "lectura"))
    lineas.append("  %-5s %-8s %-8s %-10s %-10s %-9s %s"
                  % ("", "", "", "(mejor)", "(mediana)", "mejores", ""))
    for nombre in F.FUENTES:
        etiqueta = "CR_" + nombre
        if etiqueta not in valores:
            continue
        v = valores[etiqueta]
        pm, pmed = percentil(min(v)), percentil(statistics.median(v))
        mejores_mejor = sum(1 for et in valores if et.startswith("D_")
                            and min(valores[et]) < min(v))
        mejores_med = sum(1 for et in valores if et.startswith("D_")
                          and statistics.median(valores[et]) < statistics.median(v))
        # La lectura se hace por la MEDIANA y no por el mejor modo: quedarse con el
        # mejor de 9 modos x 6 modelos es quedarse con el mas generoso de 54 intentos, y
        # eso sube a cualquiera. El veredicto de "destaca o no" lo da el bloque del
        # margen, mas abajo, que es el unico que lo compara con el ruido.
        if pmed >= 90:
            lectura = "destaca dentro del fondo"
        elif pmed >= 60:
            lectura = "algo por encima del fondo"
        elif pmed >= 40:
            lectura = "dentro del fondo"
        else:
            lectura = "por debajo del fondo"
        lineas.append("  %-5s %-8.2f %-8.2f %-10.1f %-10.1f %-9s %s"
                      % (nombre, min(v), statistics.median(v), pm, pmed,
                         "%d de %d" % (mejores_mejor, len(filas)), lectura))
    lineas.append("")
    if fondo_vals:
        lineas.append("  EL FONDO, para comparar: mejor %.2f | mediana %.2f | peor %.2f |"
                      " cuartil 25%% %.2f"
                      % (min(fondo_vals), statistics.median(fondo_vals),
                         max(fondo_vals), sorted(fondo_vals)[len(fondo_vals) // 4]))

    # --- el margen, medido contra el ruido del propio metodo ---
    # "Destacar" sobre un fondo no es quedar el primero de la lista: es quedar por
    # delante por MAS de lo que cambia un mismo compuesto al mover el receptor y la
    # semilla. Sin esta comparacion, una diferencia de decimas parece un resultado.
    ruidos = [statistics.stdev(valores[et]) for et in valores if len(valores[et]) > 1]
    ruido = statistics.median(ruidos) if ruidos else 0.0
    medias = {et: statistics.mean(valores[et]) for et in valores}
    mejor_fondo = min((medias[et], et) for et in valores if et.startswith("D_"))
    mejor_ocho = min((medias[et], et) for et in valores if et.startswith("CR_"))
    if mejor_fondo and mejor_ocho:
        margen = mejor_fondo[0] - mejor_ocho[0]
        lineas.append("")
        lineas.append("  EL MARGEN DEL PRIMERO (media de los 6 modelos, que es lo comparable)")
        lineas.append("    mejor de los ocho:   %s %.2f" % (mejor_ocho[1][3:], mejor_ocho[0]))
        lineas.append("    mejor del fondo:     %s %.2f" % (mejor_fondo[1][2:], mejor_fondo[0]))
        lineas.append("    margen: %.2f kcal/mol" % margen)
        lineas.append("    ruido del metodo aqui (desviacion tipica mediana entre los 6 modelos):"
                      " %.2f kcal/mol" % ruido)
        lineas.append("    -> %s"
                      % (("DESTACA: el margen es mayor que el ruido" if margen < -ruido else
                          "EMPATE: el margen es del tamano del ruido, no es un destacar")))
    lineas.append("")
    lineas.append("  LOS OCHO EN EL RANKING COMPLETO (por mediana de los 6 modelos, 1 = el mejor)")
    ranking = sorted(valores, key=lambda et: statistics.median(valores[et]))
    for i, et in enumerate(ranking, 1):
        if et.startswith("CR_"):
            nombre = et[3:]
            lineas.append("    puesto %3d de %d | %-5s | %.2f | %s"
                          % (i, len(valores), nombre, statistics.median(valores[et]),
                             F.FUENTES[nombre][1]))
    lineas.append("")
    lineas.append("  LO QUE ESTO DICE Y LO QUE NO")
    lineas.append("    Los señuelos son parecidos en propiedades, no inactivos conocidos: de")
    lineas.append("    ninguno hay dato de union al CR. Asi que esto no lee «se une o no se")
    lineas.append("    une», lee si el sitio distingue a unos de otros. Si los ocho caen dentro")
    lineas.append("    del fondo, y algun señuelo los supera, la respuesta es que no distingue:")
    lineas.append("    la afinidad que da Vina aqui es funcion del tamaño y la forma, no de la")
    lineas.append("    quimica que separa a estos compuestos en el laboratorio.")
    # resumen por compuesto, en CSV, para poder ordenarlo sin releer el texto
    with open(os.path.join(SALIDA, "fondo_cr.csv"), "w", encoding="utf-8",
              newline="") as f:
        w = csv.writer(f)
        w.writerow(["etiqueta", "tipo", "objetivo", "mejor", "mediana", "peor",
                    "n_modelos", "percentil_mejor", "percentil_mediana", "ventana"])
        for etiqueta in sorted(valores, key=lambda e: statistics.median(valores[e])):
            referencia, es_senuelo, ventana = etiquetas[etiqueta]
            v = valores[etiqueta]
            w.writerow([etiqueta, "señuelo" if es_senuelo else "compuesto", referencia,
                        round(min(v), 2), round(statistics.median(v), 2), round(max(v), 2),
                        len(v), ("%.1f" % percentil(min(v))) if fondo_vals else "",
                        ("%.1f" % percentil(statistics.median(v))) if fondo_vals else "",
                        ventana or ""])
    texto = "\n".join(lineas)
    with open(os.path.join(SALIDA, "fondo_cr.txt"), "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    log("")
    log(texto)
    return 0


if __name__ == "__main__":
    sys.exit(main())
