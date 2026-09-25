#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""donde_caen_positivos_sod1.py — ¿dónde se coloca de verdad cada positivo de SOD1?

POR QUE ESTO AHORA
------------------
El resultado que hundio la metrica de huella de interaccion fue este: los dos
positivos independientes de SOD1, LCS-1 y PRG-A01, comparten 0,038 de los
residuos que contactan, MENOS que dos señuelos al azar (0,297). Si dos ligandos
que deberian unirse al mismo sitio no comparten colocacion, ninguna metrica
basada en las poses puede funcionar. El problema es anterior a la metrica.

Antes de tocar nada mas hay que responder una pregunta de hecho:

    ¿dónde se coloca cada positivo cuando el acoplamiento lo busca?

Y hay un patron de medida que no admite discusion: los cuatro ligandos
co-cristalizados en Trp32 de la serie de Wright 2013, porque su sitio esta
resuelto por cristalografia, no inferido. Si ELLOS caen en Trp32 al acoplarlos y
LCS-1 no, el problema es la colocacion de LCS-1. Si no cae ninguno, la caja
Trp32 no es el sitio de estos ligandos y hay que replantearla — el mismo camino
que llevo a corregir la caja de TDP-43 cuando se descubrio que estaba a 17 A del
sitio real.

QUE HACE
--------
1. Lee la referencia de la cristalografia: para cada uno de los cuatro
   co-cristalizados (5UD 5-fluorouridina, 5FW isoprenalina, ALE epinefrina,
   LDP dopamina) saca **que residuos de Trp32 toca** en su propia estructura.
   Con ellos se construye el conjunto consenso del sitio.
2. Prepara y acopla los siete positivos de `activos_sod1_v2.csv` con la misma
   caja Trp32 de la validacion historica (centro 46.5/80.0/73.3, 22 A,
   exhaustividad 8), contra DOS receptores: el historico `gpu_dock/SOD1.pdbqt`
   (con aguas y las 18 copias del cristal) y el corregido `SOD1_limpio.pdbqt`.
3. Para cada pose: afinidad, numero de contactos con el bolsillo, si toca
   Trp32 (residuo 32) y parecido (Jaccard) con el conjunto consenso del cristal.
4. Veredicto explicito por ligando.

Los residuos del bolsillo son los del inspector anterior
(`pose_nativa_en_cribado.py`): 19-23, 28-33, 96-101 y 135.

Uso:
    python donde_caen_positivos_sod1.py            # prepara, acopla y analiza
    python donde_caen_positivos_sod1.py --solo-analizar
Salida: donde_caen_positivos_sod1.log y .csv
"""
import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
GPU = os.path.join(RAIZ, "gpu_dock")
CRISTAL = os.path.join(BASE, "trp32_cristal")
VINA_CPU = os.path.join(RAIZ, "tools", "vina.exe")
SALIDA = os.path.join(BASE, "donde_caen_positivos_sod1")

sys.path.insert(0, RAIZ)
from preparar_ligando import contar_fichero, escribir, construir   # noqa: E402

# La caja Trp32 de la validacion historica (validar_sod1_v3.py lineas 72-75).
CENTRO = (46.5, 80.0, 73.3)
TAMANO = 22
EXHAUSTIVIDAD = 8
HILOS = 6

RESIDUOS_BOLSILLO = {str(n) for n in
                     list(range(19, 24)) + list(range(28, 34)) +
                     list(range(96, 102)) + [135]}
# el residuo que da nombre al sitio
RESIDUO_TRP32 = "32"

RECEPTORES = {
    "limpio": os.path.join(GPU, "SOD1_limpio.pdbqt"),
    "historico": os.path.join(GPU, "SOD1.pdbqt"),
}

# Los cuatro co-cristalizados en Trp32 (serie Wright 2013) y su estructura.
CRISTALOGRAFICOS = {
    "5UD": ("4a7s", "5-FLUOROURIDINE", "5-fluorouridine"),
    "5FW": ("4a7t", "ISOPRENALINE", "isoproterenol"),
    "ALE": ("4a7u", "L-EPINEPHRINE", "epinephrine"),
    "LDP": ("4a7v", "L-DOPAMINE", "dopamine"),
}

CORTE_CONTACTO = 4.5


def log(m):
    linea = m
    print(linea, flush=True)
    with open(SALIDA + ".log", "a", encoding="utf-8") as f:
        f.write(linea + "\n")


# ---------------------------------------------------------------- cristalografia
def leer_pdb(ruta):
    """(proteinas, ligandos) de un PDB: listas de (resname, resnum, elem, xyz)."""
    prot, lig = [], []
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        resname = l[17:20].strip()
        resnum = l[22:26].strip()
        elem = l[76:78].strip().upper() or l[12:16].strip()[:1].upper()
        try:
            xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        except ValueError:
            continue
        (lig if l.startswith("HETATM") else prot).append((resname, resnum, elem, xyz))
    return prot, lig


def referencia_cristal():
    """Residuos de Trp32 que toca cada co-cristalizado, en su estructura."""
    log("=" * 78)
    log("REFERENCIA DE LA CRISTALOGRAFIA: que residuos toca cada co-cristalizado")
    log("=" * 78)
    refs = {}
    for cod, (pdb, resname, nombre) in CRISTALOGRAFICOS.items():
        ruta = os.path.join(CRISTAL, "pdb", pdb + ".pdb")
        if not os.path.exists(ruta):
            log("   falta %s -> se salta %s" % (ruta, cod))
            continue
        prot, lig = leer_pdb(ruta)
        # todas las copias del ligando en esa estructura
        copias = {}
        for rn, num, elem, xyz in lig:
            if rn == cod:
                copias.setdefault(num, []).append(xyz)
        if not copias:
            log("   %s: no aparece el ligando %s en %s" % (nombre, cod, pdb))
            continue
        # la copia que MAS residuos del bolsillo toca
        mejor, mejor_res = None, set()
        for num, xyzs in copias.items():
            ligx = np.array(xyzs, dtype=float)
            res = set()
            for rn, rnum, elem, xyz in prot:
                d = np.sqrt(((ligx - np.array(xyz)) ** 2).sum(axis=1)).min()
                if d <= CORTE_CONTACTO and rnum in RESIDUOS_BOLSILLO:
                    res.add(rnum)
            if len(res) > len(mejor_res):
                mejor, mejor_res = num, res
        centro = np.array(copias[mejor], dtype=float).mean(axis=0)
        refs[cod] = {"nombre": nombre, "pdb": pdb, "res": mejor_res,
                     "centro": centro, "n_copias": len(copias)}
        log("   %-16s (%s, %s) toca %d residuos del bolsillo: %s"
            % (nombre, pdb.upper(), cod, len(mejor_res),
               " ".join(sorted(mejor_res, key=int)) or "-"))

    # consenso: residuos tocados por al menos 2 de los co-cristalizados
    cuenta = {}
    for d in refs.values():
        for r in d["res"]:
            cuenta[r] = cuenta.get(r, 0) + 1
    consenso = {r for r, n in cuenta.items() if n >= 2}
    log("")
    log("   consenso (>= 2 de los 4): %s"
        % (" ".join(sorted(consenso, key=int)) or "NINGUNO"))
    log("   residuos que toca cada uno: %s"
        % ", ".join("%s:%d" % (d["nombre"], len(d["res"]))
                    for d in refs.values()))
    return refs, consenso


# ---------------------------------------------------------------------- docking
def _dock(t):
    vina, rec, lig, out, cx, cy, cz, tam, exh = t
    if os.path.exists(out) and os.path.getsize(out) > 100:
        return out
    subprocess.run([vina, "--receptor", rec, "--ligand", lig,
                    "--center_x", str(cx), "--center_y", str(cy),
                    "--center_z", str(cz),
                    "--size_x", str(tam), "--size_y", str(tam), "--size_z", str(tam),
                    "--exhaustiveness", str(exh), "--num_modes", "3",
                    "--out", out, "--cpu", "1"],
                   capture_output=True, timeout=3600)
    return out


def leer_poses(ruta):
    poses, actual = [], None
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            actual = []
        elif l.startswith("ENDMDL"):
            if actual:
                poses.append(actual)
            actual = None
        elif l.startswith(("ATOM", "HETATM")) and actual is not None:
            try:
                xyz = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
            except ValueError:
                continue
            actual.append(xyz)
    return poses


def leer_receptor(ruta):
    resnum, xyz = [], []
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        try:
            xyz.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
        except ValueError:
            continue
        resnum.append(l[22:26].strip())
    return np.array(xyz, dtype=float), resnum


def contactos(pose, rec_xyz, rec_res):
    """Residuos del bolsillo que toca una pose, y distancia al centro del bolsillo."""
    ligx = np.array(pose, dtype=float)
    d = np.sqrt(((ligx[:, None, :] - rec_xyz[None, :, :]) ** 2).sum(axis=2))
    cerca = np.nonzero((d <= CORTE_CONTACTO).any(axis=0))[0]
    res = {rec_res[j] for j in cerca if rec_res[j] in RESIDUOS_BOLSILLO}
    return res


def jaccard(a, b):
    if not a and not b:
        return None
    return len(a & b) / len(set(a) | set(b))


# ------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solo-analizar", action="store_true")
    args = ap.parse_args()

    open(SALIDA + ".log", "w", encoding="utf-8").close()
    log("¿DONDE CAEN LOS POSITIVOS DE SOD1?  %s" % time.strftime("%Y-%m-%d %H:%M"))
    log("")

    refs, consenso = referencia_cristal()

    # --- los siete positivos -------------------------------------------------
    activos = list(csv.DictReader(open(os.path.join(BASE, "activos_sod1_v2.csv"),
                                       encoding="utf-8")))
    log("")
    log("=" * 78)
    log("LOS SIETE POSITIVOS INDEPENDIENTES")
    log("=" * 78)
    for a in activos:
        log("   %-16s %-14s %s" % (a["name"], a["quimia"], a["origen"]))

    ligdir = os.path.join(SALIDA, "ligands")
    os.makedirs(ligdir, exist_ok=True)
    tareas, fallos = [], []
    for a in activos:
        # quitar_sales: las tres catecolaminas vienen con su sal en el SMILES
        # (clorhidrato, tartrato) y meeko rechaza moleculas de 2 fragmentos.
        p = escribir("ACT_" + a["name"], a["smiles"], ligdir, quitar_sales=True)
        if p is None:
            fallos.append(a["name"])
            continue
        for etq, rec in RECEPTORES.items():
            outdir = os.path.join(SALIDA, "out_" + etq)
            os.makedirs(outdir, exist_ok=True)
            out = os.path.join(outdir, "ACT_%s_out.pdbqt" % a["name"])
            tareas.append((VINA_CPU, rec, p, out, CENTRO[0], CENTRO[1], CENTRO[2],
                           TAMANO, EXHAUSTIVIDAD))
    if fallos:
        log("   NO preparados (fuera del analisis): %s" % ", ".join(fallos))

    if not args.solo_analizar:
        log("")
        log("acoplando %d ligandos x receptor (exhaustividad %d)..."
            % (len(activos) - len(fallos), EXHAUSTIVIDAD))
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=HILOS) as ex:
            list(ex.map(_dock, tareas))
        log("acoplados en %.0f s" % (time.time() - t0))

    # --- analisis ------------------------------------------------------------
    resultados = []
    for etq, rec in RECEPTORES.items():
        if not os.path.exists(rec):
            log("   falta el receptor %s" % rec)
            continue
        rec_xyz, rec_res = leer_receptor(rec)
        log("")
        log("=" * 78)
        log("RECEPTOR %s  (%s, %d atomos)"
            % (etq, os.path.basename(rec), len(rec_res)))
        log("=" * 78)
        log("   %-16s %7s %6s %8s %7s %8s  %s"
            % ("ligando", "vina", "contac", "Trp32?", "Jaccard", "cristal?", "veredicto"))
        log("   " + "-" * 92)
        for a in activos:
            nombre = a["name"]
            out = os.path.join(SALIDA, "out_" + etq, "ACT_%s_out.pdbqt" % nombre)
            if not os.path.exists(out):
                log("   %-16s   (sin pose)" % nombre)
                continue
            poses = leer_poses(out)
            if not poses:
                log("   %-16s   (pose vacia)" % nombre)
                continue
            mejor = poses[0]
            res = contactos(mejor, rec_xyz, rec_res)
            jac = jaccard(res, consenso) if consenso else None
            cod = next((c for c, d in CRISTALOGRAFICOS.items()
                        if d[2] == nombre), None)
            # veredicto
            if RESIDUO_TRP32 in res:
                veredicto = "EN Trp32"
            elif res:
                veredicto = "en el bolsillo, pero sin tocar Trp32"
            else:
                veredicto = "FUERA del bolsillo"
            if cod and cod in refs and cod in [c for c in (cod,)]:
                veredicto += "  (co-cristalizado en %s, toca %d res.)" % (
                    refs[cod]["pdb"].upper(), len(refs[cod]["res"]))
            log("   %-16s %7.2f %6d %8s %7s %8s  %s"
                % (nombre, _afinidad(out), len(res),
                   "si" if RESIDUO_TRP32 in res else "no",
                   "-" if jac is None else "%.3f" % jac,
                   cod or "-", veredicto))
            resultados.append({
                "receptor": etq, "ligando": nombre, "quimia": a["quimia"],
                "origen": a["origen"], "vina": _afinidad(out),
                "n_contactos": len(res), "toca_trp32": RESIDUO_TRP32 in res,
                "residuos": " ".join(sorted(res, key=int)),
                "jaccard_consenso": jac,
                "co_cristalizado": cod or "",
            })

    # --- resumen -------------------------------------------------------------
    log("")
    log("=" * 78)
    log("RESUMEN")
    log("=" * 78)
    for etq in RECEPTORES:
        filas = [r for r in resultados if r["receptor"] == etq]
        if not filas:
            continue
        dentro = [r for r in filas if r["toca_trp32"]]
        crist = [r for r in filas if r["co_cristalizado"]]
        crist_ok = [r for r in crist if r["toca_trp32"]]
        log("   receptor %-10s: %d de %d positivos tocan Trp32"
            % (etq, len(dentro), len(filas)))
        log("      de los 4 co-cristalizados: %d de %d" % (len(crist_ok), len(crist)))
        for r in sorted(filas, key=lambda x: x["vina"]):
            log("      %-16s %7.2f  %2d contactos  Trp32=%-3s"
                % (r["ligando"], r["vina"], r["n_contactos"],
                   "si" if r["toca_trp32"] else "no"))

    campos = ["receptor", "ligando", "quimia", "origen", "vina", "n_contactos",
              "toca_trp32", "residuos", "jaccard_consenso", "co_cristalizado"]
    with open(SALIDA + ".csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for r in resultados:
            w.writerow(r)
    log("")
    log("guardado: %s.csv" % SALIDA)
    return 0


def _afinidad(ruta):
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            try:
                return float(l.split()[3])
            except (IndexError, ValueError):
                return float("nan")
    return float("nan")


if __name__ == "__main__":
    sys.exit(main())
