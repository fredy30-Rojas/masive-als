# -*- coding: utf-8 -*-
"""¿Y si en vez de energía hay que mirar los contactos? Investigación del 17 sep 2026.

El camino hasta aquí, en cuatro frases:
  1. La dG bruta de MM-GBSA ordena al revés, y está medido que es un contador de
     tamaño disfrazado (−1,26 kcal/mol por átomo pesado en SOD1).
  2. Contra un fondo aleatorio, el acoplamiento SÍ distingue a los ligandos
     conocidos (AUC 0,66). Contra el top-100 por Vina, no distingue nada, pero
     eso es porque ese fondo está elegido por puntuar bien.
  3. Mover la caja al sitio de unión de ARN no cambia nada (0,657 → 0,655).
  4. Queda el sitio de verdad donde no se ha mirado: **la energía no dice si la
     pose toca lo que tiene que tocar**. Los ligandos conocidos de TDP-43 son
     alcaloides planos que se apilan con las plataformas aromáticas del ARN
     (berberina, sanguinarina, nitidina...). Eso no se mide con una energía: se
     mide contando contactos.

Hipótesis a probar: una puntuación de INTERACCIÓN calculada sobre la pose
(apilamiento con las aromáticas que unen ARN, puentes de hidrógeno con los
residuos polares de la interfaz, y contactos con el sitio de ARN) separa a los
conocidos del fondo aleatorio mejor que la energía, y **sin gastar GPU**, porque
se calcula sobre las poses que ya están acopladas.

Todo se mide con el mismo panel: los 8 controles y las 200 moléculas al azar
(semilla 20260917), comparados contra sus propias poses. AUC 0,5 = azar.

Uso: python puntuacion_interaccion.py [caja_arn|caja_actual|las_dos]
"""
import csv
import json
import math
import os
import random
import sys

from tdp43v2_sitio_arn import CORTE_CONTACTO, leer_modelo_1, leer_receptor

BASE = r"C:\Users\Fredy\masive-als"
ANALISIS = os.path.join(BASE, "analysis")
LIGLIB = os.path.join(BASE, "gpu_dock", "libreria_ligands")
POSES_VIEJAS = os.path.join(BASE, "gpu_dock", "resultados_libreria", "results_TDP43_v2")
POSES_ARN = os.path.join(BASE, "gpu_dock", "resultados_caja_arn", "results_TDP43_v2_arn_aleatorio")
TOTAL = os.path.join(BASE, "gpu_dock", "resultados_libreria", "resultados_libreria_total.csv")
SALIDA = os.path.join(ANALISIS, "rescoring_local", "puntuacion_interaccion.csv")
SEMILLA = 20260917
N_ALEATORIOS = 200
EXCLUIDOS = {"cepharantina", "cepharanthine", "cepharanthine_limpio"}

D_CONTACTO = 4.5      # Å para decir que dos átomos se tocan
D_APILA = 4.5         # Å entre aromaticos para apilamiento (pi-pi / T)
D_PUENTE = 3.5        # Å entre donador y aceptor de puente de hidrogeno
AROMATICAS = {"PHE", "TYR", "TRP", "HIS"}


def tipo_de(linea):
    """Tipo de atomo de AutoDock: la ultima columna de la linea."""
    return linea.rsplit(None, 1)[-1] if len(linea.split()) > 1 else ""


def es_aromatico(t):
    return t == "A"


def es_polar(t):
    return t in ("N", "NA", "OA", "O", "SA") or t in ("NS",)


def leer_ligando(ruta):
    """Atomos del PRIMER modelo de la pose: (x, y, z, tipo)."""
    atomos = []
    for ln in open(ruta, encoding="utf-8", errors="replace"):
        if ln.startswith("MODEL"):
            if atomos:      # ya tenemos el primer modelo completo
                break
            continue
        if ln.startswith(("ATOM", "HETATM")):
            try:
                atomos.append((float(ln[30:38]), float(ln[38:46]), float(ln[46:54]), tipo_de(ln)))
            except ValueError:
                continue
    return atomos


def afinidad_pose(ruta):
    for ln in open(ruta, encoding="utf-8", errors="replace"):
        if "VINA RESULT" in ln:
            return float(ln.split()[3])
    return None


def main():
    # ── Residuos que tocan al ARN, de la estructura ────────────────
    prot, rna = leer_modelo_1()
    tocan = set()
    for p in prot:
        for r in rna:
            if abs(p[3] - r[3]) > CORTE_CONTACTO:
                continue
            if math.dist((p[3], p[4], p[5]), (r[3], r[4], r[5])) <= CORTE_CONTACTO:
                tocan.add(p[0])
                break
    print("Residuos del sitio de ARN: %d" % len(tocan))

    # Receptor: átomos con su posición, su tipo de AutoDock, su residuo y su nombre
    rec_at = []
    for ln in open(os.path.join(BASE, "gpu_dock", "TDP43_v2.pdbqt"), encoding="utf-8", errors="replace"):
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        try:
            x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
        except ValueError:
            continue
        rec_at.append((x, y, z, tipo_de(ln), int(ln[22:26]), ln[17:20].strip()))
    print("Atomos del receptor: %d" % len(rec_at))

    # ── El panel: 8 controles + 200 al azar (mismo sorteo que la prueba anterior) ──
    filas = [r for r in csv.DictReader(open(os.path.join(ANALISIS, "rescoring_local",
                                                         "lista_recalculo.csv"), encoding="utf-8"))
             if r["target"] == "TDP43_v2"]
    controles = sorted({r["ligand"] for r in filas
                        if r["tipo"].startswith("control") and r["ligand"] not in EXCLUIDOS})
    todos = sorted(f[:-6] for f in os.listdir(LIGLIB) if f.endswith(".pdbqt"))
    random.seed(SEMILLA)
    aleatorios = random.sample(todos, N_ALEATORIOS)
    print("Controles: %d | fondo aleatorio: %d" % (len(controles), len(aleatorios)))

    def puntuar(ruta):
        """Saca de una pose: afinidad, tamaño y las medidas de interacción."""
        lig = leer_ligando(ruta)
        if not lig:
            return None
        pesados = [a for a in lig if a[3] != "HD" and a[3] not in ("H",)]
        afinidad = afinidad_pose(ruta)
        n_contacto = n_arn = n_apila = n_puente = 0
        for (lx, ly, lz, lt) in pesados:
            for (rx, ry, rz, rt, rnum, rname) in rec_at:
                d2 = (lx - rx) ** 2 + (ly - ry) ** 2 + (lz - rz) ** 2
                if d2 > D_CONTACTO ** 2:
                    continue
                n_contacto += 1
                if rnum in tocan:
                    n_arn += 1
                if es_aromatico(lt) and rname in AROMATICAS and d2 <= D_APILA ** 2:
                    n_apila += 1
                if es_polar(lt) and es_polar(rt) and d2 <= D_PUENTE ** 2:
                    n_puente += 1
        n_heavy = len(pesados)
        return {
            "afinidad": afinidad,
            "atomos": n_heavy,
            "contactos": n_contacto,
            "contactos_arn": n_arn,
            "apilamiento": n_apila,
            "puentes": n_puente,
            "contactos_por_atomo": n_contacto / n_heavy if n_heavy else 0,
            "arn_por_atomo": n_arn / n_heavy if n_heavy else 0,
            "apila_por_atomo": n_apila / n_heavy if n_heavy else 0,
            "interaccion": n_apila + n_puente,          # lo que define el reconocimiento del ARN
            "interaccion_por_atomo": (n_apila + n_puente) / n_heavy if n_heavy else 0,
        }

    datos = {}
    for nombre in aleatorios + controles:
        grupo = "control" if nombre in controles else "fondo"
        fila = {"ligand": nombre, "grupo": grupo}
        p_vieja = os.path.join(POSES_VIEJAS, nombre + "_out.pdbqt")
        if os.path.exists(p_vieja):
            for k, v in (puntuar(p_vieja) or {}).items():
                fila["actual_" + k] = v
        p_arn = os.path.join(POSES_ARN, nombre + "_out.pdbqt")
        if os.path.exists(p_arn):
            for k, v in (puntuar(p_arn) or {}).items():
                fila["arn_" + k] = v
        datos[nombre] = fila

    # ── AUC de cada medida, en cada caja ──────────────────────────
    def auc_campo(campo, mejor_es_menor):
        ct = [datos[c][campo] for c in controles if datos[c].get(campo) is not None]
        fo = [datos[a][campo] for a in aleatorios if datos[a].get(campo) is not None]
        if not ct or not fo:
            return None, 0, 0
        gana = 0.0
        for x in ct:
            for y in fo:
                if x == y:
                    gana += 0.5
                elif (x < y) if mejor_es_menor else (x > y):
                    gana += 1.0
        return gana / (len(ct) * len(fo)), len(ct), len(fo)

    medidas = [
        ("afinidad", "afinidad", True, "energía de acoplamiento (la de siempre)"),
        ("contactos", "contactos", False, "nº de contactos ligando-receptor"),
        ("contactos_por_atomo", "contactos_por_atomo", False, "contactos por átomo pesado"),
        ("contactos_arn", "contactos_arn", False, "contactos con los residuos que unen ARN"),
        ("arn_por_atomo", "arn_por_atomo", False, "esos contactos por átomo pesado"),
        ("apilamiento", "apilamiento", False, "apilamiento con aromáticas (bases del ARN)"),
        ("puentes", "puentes", False, "puentes de hidrógeno"),
        ("interaccion", "interaccion", False, "apilamiento + puentes"),
        ("interaccion_por_atomo", "interaccion_por_atomo", False, "apilamiento + puentes por átomo"),
    ]

    resumen = {}
    for caja, prefijo in (("caja actual", "actual_"), ("caja de ARN", "arn_")):
        print("\n================ %s ================" % caja.upper())
        print("%-24s %-7s %s" % ("medida", "AUC", "qué mide"))
        for etiqueta, campo, menor, desc in medidas:
            a, nc, nf = auc_campo(prefijo + campo, menor)
            if a is None:
                print("%-24s %-7s %s" % (etiqueta, "n/d", desc))
                continue
            resumen["%s|%s" % (caja, etiqueta)] = round(a, 4)
            print("%-24s %-7.3f %s  (n=%d controles, %d fondo)" % (etiqueta, a, desc, nc, nf))

    # ── Guardar la tabla completa para poder seguir tirando del hilo ──
    campos = ["ligand", "grupo"]
    for k in sorted({k for f in datos.values() for k in f} - {"ligand", "grupo"}):
        campos.append(k)
    with open(SALIDA, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for nombre in aleatorios + controles:
            w.writerow(datos[nombre])
    with open(SALIDA.replace(".csv", "_auc.json"), "w", encoding="utf-8") as f:
        json.dump(resumen, f, indent=2, ensure_ascii=False)
    print("\nTabla: %s" % SALIDA)
    print("AUC por medida: %s" % SALIDA.replace(".csv", "_auc.json"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
