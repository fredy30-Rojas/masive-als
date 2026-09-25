#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""reparar_validaciones_glue.py — arregla el sesgo «glue» en las validaciones.

POR QUE
-------
Meeko, con los ajustes por defecto, no cierra los anillos de 7 eslabones en
adelante: los parte en dos ramas y pega los extremos con pseudo-atomos de
pegamento (tipos CG0 y G0), DUPLICANDO dos atomos por anillo abierto.

Esos ficheros envenenan la validacion de dos maneras:
  * Vina ve el anillo abierto y su termino de energia interna del ligando LIBRE
    (REMARK UNBOUND) sale a +8 o +15 kcal/mol en vez de ~0. Como el score que
    reporta Vina es inter + intra - unbound, el compuesto recibe un regalo de
    varios kcal/mol que no existen.
  * El ligando tiene 2 atomos de mas, asi que cualquier normalizacion por
    tamano queda tambien mal.

El proyecto ya lo sabia para la LIBRERIA (`gpu_dock/reparar_libreria_glue.py`,
3.803 afectados) y se re-acoplaron solo contra TDP43_v2. Los conjuntos de
VALIDACION se quedaron con los ficheros rotos, y en todos ellos el glue ocupa la
cabeza del ranking:

    TDP43 v6 (4BS2)   puestos 1-15 de 131   (14 senuelos + la activa cefarantina)
    FUS v2            puestos  1-3  de  89
    SOD1 (v1)         puestos  1-2  de 219
    SOD1 v3           11 afectados

Eso invalida el diagnostico escrito en `REVISION_VALIDACION_BOLSILLOS`: los
decoys top de TDP43 NO son «compuestos grandes que puntuan mejor por artefacto
de tamano de Vina». Son moleculas normales (22-54 atomos) cuyo score estaba
inflado por el anillo abierto.

QUE HACE
--------
1. Detecta los ligandos con tipos CG0/G0 de cada conjunto.
2. Guarda aparte la pose vieja (en `<out>/_antes_glue/`) y el ligando viejo
   (en `<ligands>/_roto_glue/`), no borra nada.
3. Vuelve a preparar el ligando con `rigid_macrocycles=True` y comprueba que el
   numero de atomos coincide con su SMILES.
4. Re-acopla SOLO esos ligandos con el MISMO receptor, caja y exhaustividad del
   conjunto original, para que el antes y el despues sean comparables.
5. Recalcula las metricas (AUC crudo, por atomo pesado y por raiz de atomos)
   antes y despues de la reparacion.

Uso:
  python reparar_validaciones_glue.py            # todos los conjuntos
  python reparar_validaciones_glue.py --solo tdp43_v6
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

BASE = os.path.dirname(os.path.abspath(__file__))
GPU = os.path.abspath(os.path.join(BASE, "..", "gpu_dock"))
VINA_CPU = os.path.abspath(os.path.join(BASE, "..", "tools", "vina.exe"))

sys.path.insert(0, BASE)
from reparar_ligandos_glue import contar_fichero, reparar   # noqa: E402

LOG = os.path.join(BASE, "reparar_validaciones_glue.log")
HILOS = 6

# receptor, centro, tamano, exhaustividad: los del lanzamiento ORIGINAL de cada
# conjunto. La exhaustividad de TDP43 v6 (16) esta en el informe del bolsillo
# 4BS2; la de FUS (8) se dedujo reproduciendo el score de un ligando intacto.
CONJUNTOS = {
    "tdp43_v6": {
        "nombre": "TDP43 v6 (4BS2, interfaz RRM1-RRM2)",
        "ligandos": os.path.join(BASE, "_validacion_TDP43", "ligands"),
        "poses": os.path.join(BASE, "_validacion_TDP43", "out"),
        "csv": os.path.join(BASE, "validacion_TDP43_v6_4BS2.csv"),
        "receptor": os.path.join(BASE, "_tdp43_bolsillo_v2", "4BS2_ph74.pdbqt"),
        "centro": (24.23, 16.89, -15.87),
        "tamano": 26,
        "exhaustividad": 16,
    },
    "fus_v2": {
        "nombre": "FUS v2 (modelo 1)",
        "ligandos": os.path.join(BASE, "_validacion_FUS", "ligands"),
        "poses": os.path.join(BASE, "_validacion_FUS", "out"),
        "csv": os.path.join(BASE, "validacion_FUS_v2.csv"),
        "receptor": os.path.join(GPU, "FUS.pdbqt"),
        "centro": (-14.5, 15.1, -7.8),
        "tamano": 18,
        "exhaustividad": 8,
    },
    "sod1_v4": {
        # La v4 se quedo sin reparar cuando se reparo la v5, y sus dos ficheros
        # con pseudo-atomos son justo los que dieron el «positivo» falso
        # (4MQ -10,07 y ZO0 -9,83). Los ligands de la v5 son los de la v4, asi que
        # este conjunto es el original de la ronda que produjo aquel numero.
        "nombre": "SOD1 v4 (controles de la ronda del 20 sep)",
        "ligandos": os.path.join(BASE, "validacion_SOD1_v4", "ligands"),
        "poses": os.path.join(BASE, "validacion_SOD1_v4", "out"),
        "csv": os.path.join(BASE, "validacion_SOD1_v4", "no_existe.csv"),
        "receptor": os.path.join(GPU, "SOD1.pdbqt"),
        "centro": (46.5, 80.0, 73.3),
        "tamano": 22,
        "exhaustividad": 8,
    },
    "sod1_trp32": {
        # El conjunto historico del bolsillo Trp32 (219 ligandos, 2 activos).
        # Se conserva a proposito el receptor ORIGINAL gpu_dock/SOD1.pdbqt (el
        # que lleva aguas y las 18 copias del cristal) para que el «antes» y el
        # «despues» sean comparables: aqui solo se mide el efecto del glue. La
        # re-validacion con el receptor limpio es la v5.
        "nombre": "SOD1 trp32 (v1 historica)",
        "ligandos": os.path.join(BASE, "_validacion_SOD1", "ligands"),
        "poses": os.path.join(BASE, "_validacion_SOD1", "out"),
        "csv": os.path.join(BASE, "_validacion_SOD1", "validacion_SOD1_trp32.csv"),
        "receptor": os.path.join(GPU, "SOD1.pdbqt"),
        "centro": (46.5, 80.0, 73.3),
        "tamano": 22,
        "exhaustividad": 8,
    },
}


def norm(n):
    """Nombre comparable entre el CSV y el fichero.

    El CSV de TDP43 lista `CHEMBL607833` y el fichero es `DEC_CHEMBL607833`;
    el de FUS lista `DEC_Sertraline` y el fichero es igual. Se quitan prefijos.
    """
    n = os.path.basename(n).replace(".pdbqt", "")
    for p in ("ACT_", "DEC_", "DECH_", "DECM_"):
        if n.startswith(p):
            return n[len(p):]
    return n


def log(m):
    linea = "[%s] %s" % (time.strftime("%H:%M:%S"), m)
    print(linea, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(linea + "\n")


def leer_afinidad(p, modelo=1):
    """Afinidad del modelo pedido de una pose PDBQT de Vina."""
    visto = 0
    for l in open(p, encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            visto += 1
        if l.startswith("REMARK VINA RESULT:") and visto == modelo:
            try:
                return float(l.split()[3])
            except Exception:
                return None
    return None


def leer_unbound(p):
    for l in open(p, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK UNBOUND:"):
            try:
                return float(l.split(":")[1].split()[0])
            except Exception:
                return None
    return None


def _dock(t):
    vina, rec, lig, out, cx, cy, cz, tam, exh = t
    cmd = [vina, "--receptor", rec, "--ligand", lig,
           "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", str(tam), "--size_y", str(tam), "--size_z", str(tam),
           "--exhaustiveness", str(exh), "--num_modes", "3",
           "--out", out, "--cpu", "1"]
    subprocess.run(cmd, capture_output=True, timeout=3600)
    return out


def auc(act, dec):
    if not act or not dec:
        return None
    n = 0.0
    for a in act:
        n += sum(1 for d in dec if d >= a)
    return n / (len(act) * len(dec))


def tres_auc(scores, rol, peso):
    """AUC crudo, por atomo pesado y por raiz de atomos de un diccionario de scores.

    `peso` puede faltar (los conjuntos viejos no lo traen); en ese caso solo se
    devuelve el crudo.
    """
    act = [(k, v) for k, v in scores.items() if rol.get(k) == "activo"]
    dec = [(k, v) for k, v in scores.items() if rol.get(k) != "activo"]
    crudo = auc([v for _, v in act], [v for _, v in dec])
    if not peso or not all(k in peso for k, _ in act + dec):
        return crudo, None, None
    por_at = auc([v / peso[k] for k, v in act], [v / peso[k] for k, v in dec])
    por_rz = auc([v / (peso[k] ** 0.5) for k, v in act],
                 [v / (peso[k] ** 0.5) for k, v in dec])
    return crudo, por_at, por_rz


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solo", default=None, help="clave de CONJUNTOS")
    args = ap.parse_args()

    for clave, c in CONJUNTOS.items():
        if args.solo and clave != args.solo:
            continue
        if os.path.exists(c["csv"]):
            filas_csv = list(csv.DictReader(open(c["csv"], encoding="utf-8",
                                               errors="ignore")))
            col = "affinity" if "affinity" in filas_csv[0] else "energy"
        else:
            # Sin tabla: el papel sale del prefijo y la afinidad de la propia pose.
            filas_csv = []
            for p in sorted(os.listdir(c["poses"])):
                if not p.endswith("_out.pdbqt"):
                    continue
                nombre = p.replace("_out.pdbqt", "")
                aff = leer_afinidad(os.path.join(c["poses"], p))
                if aff is None:
                    continue
                filas_csv.append({
                    "ligand": nombre,
                    "rol": "activo" if nombre.startswith("ACT_") else "decoy",
                    "affinity": aff})
            col = "affinity"
            log("    (falta la tabla %s: se deducen %d ligandos de las poses)"
                % (os.path.basename(c["csv"]), len(filas_csv)))
        rol = {norm(r["ligand"]): r["rol"] for r in filas_csv}
        peso = {}
        for r in filas_csv:
            for k in ("heavy_atoms", "atomos"):
                if r.get(k):
                    try:
                        peso[norm(r["ligand"])] = float(r[k])
                    except ValueError:
                        pass
        afin_antes = {}
        for r in filas_csv:
            try:
                afin_antes[norm(r["ligand"])] = float(r[col])
            except (ValueError, TypeError):
                pass

        log("")
        log("=" * 78)
        log("=== %s" % c["nombre"])
        log("    ligando   %s" % c["ligandos"])
        log("    receptor  %s" % os.path.basename(c["receptor"]))
        log("    caja      %s tamano %d exhaustividad %d"
            % (c["centro"], c["tamano"], c["exhaustividad"]))

        rotos = []
        for p in sorted(os.listdir(c["ligandos"])):
            if not p.endswith(".pdbqt"):
                continue
            f = os.path.join(c["ligandos"], p)
            n, glue = contar_fichero(f)
            if glue:
                rotos.append((os.path.splitext(p)[0], n, glue))

        if not rotos:
            log("    ningun ligando con pseudo-atomos. Nada que hacer.")
            continue
        log("    ligandos con pseudo-atomos: %d" % len(rotos))

        os.makedirs(os.path.join(c["poses"], "_antes_glue"), exist_ok=True)
        os.makedirs(os.path.join(c["ligandos"], "_roto_glue"), exist_ok=True)

        tareas, cambios = [], []
        for fichero, n, glue in rotos:
            clave = norm(fichero)
            p = os.path.join(c["ligandos"], fichero + ".pdbqt")
            pose = os.path.join(c["poses"], fichero + "_out.pdbqt")
            unbound_antes = leer_unbound(pose) if os.path.exists(pose) else None
            if os.path.exists(pose):
                shutil.copy(pose, os.path.join(c["poses"], "_antes_glue", os.path.basename(pose)))
            txt, n_ok, motivo = reparar(p)
            if txt is None:
                log("    %-26s NO reparado: %s" % (fichero, motivo))
                continue
            shutil.copy(p, os.path.join(c["ligandos"], "_roto_glue", fichero + ".pdbqt"))
            with open(p, "w", encoding="utf-8") as fh:
                fh.write(txt)
            if os.path.exists(pose):
                os.remove(pose)
            tareas.append((VINA_CPU, c["receptor"], p, pose,
                           c["centro"][0], c["centro"][1], c["centro"][2],
                           c["tamano"], c["exhaustividad"]))
            cambios.append((fichero, clave, n, n_ok, glue, unbound_antes))

        log("    re-acoplando %d ligandos (Vina CPU, exhaustividad %d)..."
            % (len(tareas), c["exhaustividad"]))
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=HILOS) as ex:
            list(ex.map(_dock, tareas))
        log("    acoplados en %.0f s" % (time.time() - t0))

        log("")
        log("    %-26s %-7s %9s %9s %9s %10s"
            % ("ligando", "rol", "antes", "despues", "cambio", "UNBOUND"))
        log("    " + "-" * 74)
        afin_despues = dict(afin_antes)
        for fichero, clave, n, n_ok, glue, unb in cambios:
            pose = os.path.join(c["poses"], fichero + "_out.pdbqt")
            nuevo = leer_afinidad(pose)
            viejo = afin_antes.get(clave)
            if nuevo is None:
                log("    %-26s %-7s %9s %9s   (sin pose)"
                    % (fichero, rol.get(clave, "?"),
                       "%.2f" % viejo if viejo is not None else "-", "-"))
                continue
            afin_despues[clave] = nuevo
            log("    %-26s %-7s %9.2f %9.2f %+9.2f %10s"
                % (fichero, rol.get(clave, "?"),
                   viejo if viejo is not None else float("nan"), nuevo,
                   (nuevo - viejo) if viejo is not None else float("nan"),
                   "%.2f" % unb if unb is not None else "-"))

        # --- metricas antes / despues -------------------------------------
        log("")
        log("    metricas (activos=%d, decoys=%d)"
            % (sum(1 for v in rol.values() if v == "activo"),
               sum(1 for v in rol.values() if v != "activo")))
        for etiqueta, dic in (("ANTES", afin_antes), ("DESPUES", afin_despues)):
            cruda, por_at, por_rz = tres_auc(dic, rol, peso)
            if cruda is None:
                # sin fondo no hay AUC que calcular (conjuntos de solo controles)
                log("    %-9s AUC n/d (no hay señuelos en este conjunto)" % etiqueta)
                continue
            extra = ""
            if por_at is not None:
                extra = " | por atomo %.3f | por raiz %.3f" % (por_at, por_rz)
            log("    %-9s AUC crudo %.3f%s" % (etiqueta, cruda, extra))

    log("")
    log("informe de consola replicable; originales guardados en _roto_glue/ y _antes_glue/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
