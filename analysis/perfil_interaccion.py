# -*- coding: utf-8 -*-
"""La huella de interaccion 3D, leida de las poses ya calculadas (20 sep 2026).

POR QUE ESTO Y NO OTRA COSA

Las dos medidas que el proyecto tenia estan rotas, y por motivos distintos:

  * La ENERGIA (Vina, MM-GBSA) esta dominada por el tamano: -0,59 a -1,26
    kcal/mol por atomo pesado. Ordena al mas grande, no al mejor.
  * La SIMILITUD 2D (ECFP4) es un detector de familia: da AUC 1,000 dentro de
    su propia familia y 0,538 fuera de ella.

Ninguna de las dos mira lo unico que de verdad describe una union: DONDE y
COMO se apoya la molecula en la proteina. Eso si esta en las poses, que ya
estan calculadas y guardadas.

Aqui se construye la HUELLA DE INTERACCION (IFP, interaction fingerprint):
para cada pose, la lista de atomos del receptor que el ligando toca, con su
residuo y su tipo. La huella no depende del tamano del ligando ni de si se
parece a otro en 2D: depende de con que parte de la proteina habla.

Puntuacion: parecido entre la huella del ligando y la huella de consenso de
los activos (a cada activo se le mide contra el consenso SIN el mismo, para no
inflarlo). Se calculan tres, para ver cual aguanta:
    recall   = fraccion del consenso que el ligando satisface (premia cubrir)
    jaccard  = interseccion / union (castiga los contactos de relleno)
    coseno   = parecido de vectores (intermedio)

Uso:
    python perfil_interaccion.py                 # SOD1, conjunto de validacion
    python perfil_interaccion.py --receptor X --poses CARPETA --tabla Y.csv
"""
import argparse
import csv
import math
import os
import statistics
import sys

import numpy as np

BASE = r"C:\Users\Fredy\masive-als"
ANALISIS = os.path.join(BASE, "analysis")

# Corte por tipo de atomo: los polares se tocan desde mas lejos que los carbonos
CORTE = {"polar": 3.9, "carbon": 4.5, "metal": 3.0, "halogeno": 4.0}
TIPOS_POLARES = {"N", "NA", "OA", "O", "SA", "NS"}
TIPOS_CARBONO = {"C", "A"}
TIPOS_METAL = {"Cu", "Zn", "Ca", "Fe", "Mn", "Mg"}
TIPOS_HALOGENO = {"Cl", "Br", "F", "I"}


def categoria(tipo):
    if tipo in TIPOS_POLARES:
        return "polar"
    if tipo in TIPOS_METAL:
        return "metal"
    if tipo in TIPOS_HALOGENO:
        return "halogeno"
    return "carbon"


def leer_receptor(ruta):
    """Atomos del receptor: coordenadas, residuo y tipo."""
    nombres, resnames, resnums, tipos, coords = [], [], [], [], []
    with open(ruta, "r", errors="ignore") as f:
        for linea in f:
            if not linea.startswith(("ATOM", "HETATM")):
                continue
            try:
                x = float(linea[30:38]); y = float(linea[38:46]); z = float(linea[46:54])
            except ValueError:
                continue
            partes = linea.split()
            tipo = partes[-1] if partes else "C"
            # El tipo de AutoDock es la ultima columna; si viene corrida, cae en
            # el nombre del elemento, que es suficientemente bueno.
            if len(tipo) > 2:
                tipo = linea[76:79].strip() or "C"
            nombres.append(linea[12:16].strip())
            resnames.append(linea[17:20].strip())
            resnums.append(linea[22:26].strip())
            tipos.append(tipo)
            coords.append((x, y, z))
    return {"nombre": nombres, "resname": resnames, "resnum": resnums,
            "tipo": tipos, "xyz": np.array(coords, dtype=np.float64)}


def leer_poses(ruta):
    """Todas las poses (todos los MODEL) de un fichero de Vina."""
    poses, actual = [], None
    with open(ruta, "r", errors="ignore") as f:
        for linea in f:
            if linea.startswith("MODEL"):
                actual = []
            elif linea.startswith("ENDMDL"):
                if actual:
                    poses.append(actual)
                actual = None
            elif linea.startswith(("ATOM", "HETATM")) and actual is not None:
                try:
                    x = float(linea[30:38]); y = float(linea[38:46]); z = float(linea[46:54])
                except ValueError:
                    continue
                partes = linea.split()
                tipo = partes[-1] if partes else "C"
                if len(tipo) > 2:
                    tipo = linea[76:79].strip() or "C"
                actual.append((x, y, z, tipo))
    return poses


def huella(pose, receptor, residuos_diana=None):
    """Atomos del receptor que toca esta pose, como conjunto de etiquetas."""
    if not pose:
        return {}
    lig = np.array([[p[0], p[1], p[2]] for p in pose], dtype=np.float64)
    tipos_lig = [p[3] for p in pose]
    xyz = receptor["xyz"]
    # Distancias de cada atomo del ligando a todos los del receptor (vectorizado)
    d = np.sqrt(((lig[:, None, :] - xyz[None, :, :]) ** 2).sum(axis=2))
    cuenta = {}
    for i, tipo_l in enumerate(tipos_lig):
        cat_l = categoria(tipo_l)
        for j in np.nonzero(d[i] <= CORTE["carbon"])[0]:
            tipo_r = receptor["tipo"][j]
            cat_r = categoria(tipo_r)
            if residuos_diana and receptor["resnum"][j] not in residuos_diana:
                continue
            # El corte fino depende de los dos extremos
            corte = CORTE["metal"] if cat_r == "metal" else (
                CORTE["polar"] if (cat_l == "polar" and cat_r == "polar") else
                CORTE["halogeno"] if cat_l == "halogeno" else CORTE["carbon"])
            if d[i][j] > corte:
                continue
            etiqueta = "%s:%s:%s|%s-%s" % (receptor["resnum"][j], receptor["resname"][j],
                                           receptor["nombre"][j], cat_l, cat_r)
            cuenta[etiqueta] = cuenta.get(etiqueta, 0) + 1
    return cuenta


def vector(cuenta, vocabulario):
    v = np.zeros(len(vocabulario), dtype=np.float64)
    for i, etq in enumerate(vocabulario):
        if etq in cuenta:
            v[i] = 1.0
    return v


def puntuar(v, consenso, metodo):
    inter = float((v * consenso).sum())
    if metodo == "recall":
        base = consenso.sum()
        return inter / base if base else 0.0
    if metodo == "jaccard":
        union = float(((v + consenso) > 0).sum())
        return inter / union if union else 0.0
    # coseno
    nv = math.sqrt(float((v * v).sum()))
    nc = math.sqrt(float((consenso * consenso).sum()))
    return inter / (nv * nc) if nv and nc else 0.0


def auc(a, b):
    if not a or not b:
        return None
    g = 0.0
    for x in a:
        for y in b:
            g += 0.5 if x == y else (1.0 if x > y else 0.0)
    return g / (len(a) * len(b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", default=os.path.join(BASE, "gpu_dock", "SOD1.pdbqt"))
    ap.add_argument("--poses", default=os.path.join(ANALISIS, "_validacion_SOD1", "out"))
    ap.add_argument("--tabla", default=os.path.join(ANALISIS, "_validacion_SOD1",
                                                    "validacion_SOD1_trp32.csv"),
                    help="CSV con ligand,rol,affinity. Si no existe, se deduce de "
                         "los nombres de las poses (ACT_ = activo, resto = señuelo)")
    ap.add_argument("--salida", default=os.path.join(ANALISIS, "rescoring_local",
                                                     "perfil_interaccion_sod1.csv"))
    ap.add_argument("--umbral-consenso", type=float, default=0.5,
                    help="fraccion de activos que debe tener una etiqueta para entrar en el consenso")
    ap.add_argument("--referencia", default="",
                    help="activos que definen el consenso, separados por comas. "
                         "Por defecto, todos los de la tabla.")
    args = ap.parse_args()

    print("Leyendo receptor: %s" % args.receptor)
    receptor = leer_receptor(args.receptor)
    print("   atomos del receptor: %d" % len(receptor["tipo"]))

    if os.path.exists(args.tabla):
        filas = list(csv.DictReader(open(args.tabla, encoding="utf-8")))
        print("Ligandos en la tabla de validacion: %d" % len(filas))
    else:
        # Sin tabla: el papel sale del prefijo del nombre y la afinidad de la pose.
        filas = []
        for p in sorted(os.listdir(args.poses)):
            if not p.endswith("_out.pdbqt"):
                continue
            nombre = p.replace("_out.pdbqt", "")
            casa = leer_poses(os.path.join(args.poses, p))
            aff = None
            with open(os.path.join(args.poses, p), encoding="utf-8", errors="ignore") as fh:
                for l in fh:
                    if l.startswith("REMARK VINA RESULT:"):
                        try:
                            aff = float(l.split()[3])
                        except (IndexError, ValueError):
                            pass
                        break
            if not casa or aff is None:
                continue
            filas.append({"ligand": nombre,
                          "rol": "activo" if nombre.startswith("ACT_") else "decoy",
                          "affinity": aff})
        print("Falta %s: se deducen %d ligandos de las poses"
              % (args.tabla, len(filas)))

    # Huellas de todos los ligandos
    print("Calculando la huella de interaccion de cada pose...")
    datos = {}
    for k, r in enumerate(filas, 1):
        nombre = r["ligand"]
        ruta = os.path.join(args.poses, nombre + "_out.pdbqt")
        if not os.path.exists(ruta):
            continue
        poses = leer_poses(ruta)
        if not poses:
            continue
        # Se queda la mejor pose (la primera: Vina las escribe ordenadas)
        h = huella(poses[0], receptor)
        if not h:
            continue
        datos[nombre] = {"rol": r["rol"], "vina": float(r["affinity"]), "huella": h}
        if k % 50 == 0:
            print("   %d/%d" % (k, len(filas)))
    print("Con huella: %d de %d" % (len(datos), len(filas)))

    act_nombres = [n for n, d in datos.items() if d["rol"] == "activo"]
    dec_nombres = [n for n, d in datos.items() if d["rol"] == "decoy"]
    print("Activos: %d | señuelos: %d" % (len(act_nombres), len(dec_nombres)))

    # Quienes definen el consenso: por defecto todos los activos; con
    # --referencia, solo los que se indiquen (para probar quimias sueltas).
    ref_nombres = act_nombres
    if args.referencia.strip():
        pedidos = [x.strip() for x in args.referencia.split(",") if x.strip()]
        ref_nombres = []
        for p in pedidos:
            cand = [n for n in datos if n == p or n.replace("ACT_", "") == p]
            if not cand:
                print("   !! no esta en los datos: %s" % p)
            ref_nombres.extend(cand)
        print("Consenso definido SOLO por: %s" % ", ".join(ref_nombres))

    # Vocabulario: todas las etiquetas que aparecen en cualquier huella
    vocab = sorted({e for d in datos.values() for e in d["huella"]})
    print("Tamano del vocabulario de interacciones: %d" % len(vocab))
    idx = {e: i for i, e in enumerate(vocab)}
    for d in datos.values():
        d["vec"] = vector(d["huella"], vocab)

    # --- Vina, como referencia ---
    print("\n" + "=" * 74)
    print("REFERENCIA: VINA (menos energia = mejor)")
    print("=" * 74)
    a_v = [-datos[n]["vina"] for n in act_nombres]
    d_v = [-datos[n]["vina"] for n in dec_nombres]
    print("AUC con los %d activos: %.3f" % (len(act_nombres), auc(a_v, d_v)))

    # --- La huella de consenso, sin contar al propio activo ---
    print("\n" + "=" * 74)
    print("LA HUELLA DE INTERACCION 3D")
    print("=" * 74)
    resultados = {}
    for metodo in ("recall", "jaccard", "coseno"):
        for n in datos:
            otros = [m for m in ref_nombres if m != n] if n in ref_nombres else ref_nombres
            if not otros:
                continue
            if not otros:
                continue
            # Consenso: etiquetas presentes en al menos el umbral de los activos
            M = np.vstack([datos[m]["vec"] for m in otros])
            consenso = (M.mean(axis=0) >= args.umbral_consenso).astype(np.float64)
            if metodo == "jaccard":
                consenso = (M.mean(axis=0) > 0).astype(np.float64)
            datos[n][metodo] = puntuar(datos[n]["vec"], consenso, metodo)
        a_m = [datos[n][metodo] for n in ref_nombres if metodo in datos[n]]
        d_m = [datos[n][metodo] for n in dec_nombres if metodo in datos[n]]
        v = auc(a_m, d_m)
        resultados[metodo] = v
        print("%-10s AUC %.3f  (activos de referencia %d, señuelos %d)"
          % (metodo, v, len(a_m), len(d_m)))
        # Puestos de los activos de referencia y del resto de activos
        ordenados = sorted(datos.items(), key=lambda x: -(x[1].get(metodo) or 0))
        for n in ref_nombres:
            if n in datos:
                print("     %-18s puesto %3d/%d   huella=%d" % (
                    n.replace("ACT_", ""), ordenados.index((n, datos[n])) + 1,
                    len(ordenados), len(datos[n]["huella"])))

    # --- Por quimias: la prueba que importa ---
    print("\n" + "=" * 74)
    print("LA PRUEBA QUE IMPORTA: por quimias, no por el conjunto entero")
    print("=" * 74)
    v2 = os.path.join(ANALISIS, "activos_sod1_v2.csv")
    if os.path.exists(v2):
        quimias = {r["name"].upper(): r["quimia"]
                   for r in csv.DictReader(open(v2, encoding="utf-8"))}
        grupos = {}
        for n in act_nombres:
            base = n.replace("ACT_", "").upper()
            q = quimias.get(base, "?")
            grupos.setdefault(q, []).append(n)
        print("Quimias representadas en los activos de la validacion:")
        for q, ns in sorted(grupos.items(), key=lambda x: -len(x[1])):
            print("   %-22s %d compuestos" % (q, len(ns)))
        # AUC de cada quimia por separado, con la mejor medida
        mejor = max(resultados, key=lambda k: resultados[k] or 0)
        print("\nAUC por quimia (medida: %s):" % mejor)
        for q, ns in sorted(grupos.items()):
            a_q = [datos[n][mejor] for n in ns if mejor in datos[n]]
            if not a_q:
                continue
            print("   %-22s AUC %.3f  (n=%d)  vina %.3f"
                  % (q, auc(a_q, d_v and [datos[n][mejor] for n in dec_nombres]),
                     len(a_q),
                     auc([-datos[n]["vina"] for n in ns], d_v)))
        # Y el caso concreto: LCS-1
        for especial in ("ACT_LCS-1", "ACT_PRG-A01"):
            if especial in datos:
                puesto = 1 + sum(1 for n in dec_nombres if datos[n][mejor] > datos[especial][mejor])
                print("\n   %s: puesto %d de %d por %s  |  por Vina: puesto %d"
                      % (especial.replace("ACT_", ""), puesto, len(datos), mejor,
                         1 + sum(1 for n in dec_nombres if datos[n]["vina"] < datos[especial]["vina"])))
    else:
        print("(Falta %s: sin el no se puede agrupar por quimias)" % v2)

    # --- Guardar ---
    os.makedirs(os.path.dirname(args.salida), exist_ok=True)
    campos = ["ligand", "rol", "vina", "recall", "jaccard", "coseno", "n_contactos"]
    with open(args.salida, "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(campos)
        for n, d in sorted(datos.items(), key=lambda x: -(x[1].get("recall") or 0)):
            w.writerow([n, d["rol"], d["vina"],
                        round(d.get("recall", 0), 4), round(d.get("jaccard", 0), 4),
                        round(d.get("coseno", 0), 4), len(d["huella"])])
    print("\nGuardado: %s" % args.salida)

    # Correlacion con el tamano, que es lo que hundia a la energia
    tam = np.array([len(d["huella"]) for d in datos.values()], dtype=float)
    rec = np.array([d.get("recall", 0) for d in datos.values()], dtype=float)
    vina = np.array([d["vina"] for d in datos.values()], dtype=float)
    if len(tam) > 2:
        print("\nCorrelacion con el numero de contactos (proxy del tamano):")
        print("   Vina   : %+.3f   <- cuanto mas grande, mejor energia" % np.corrcoef(tam, vina)[0, 1])
        print("   huella : %+.3f   <- cuanto mas grande, mas huella" % np.corrcoef(tam, rec)[0, 1])
    return 0


if __name__ == "__main__":
    sys.exit(main())
