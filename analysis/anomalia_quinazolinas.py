#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""La anomalia de las quinazolinas: de donde salen los -10 kcal/mol (21 sep 2026).

EL DATO QUE HAY QUE EXPLICAR
----------------------------
En la validacion v4 (INFORME_VALIDACION_SOD1_V4_2026-09-20.md) tres quinazolinas
practicamente identicas se separan CINCO kcal/mol dentro del mismo receptor y la
misma caja:

    4MQ  (diazepano, 7 anillos)   -10,02
    ZO0  (diazepano + CF3)         -9,77
    12I  (piperazina, 6 anillos)   -4,97

Un CH2 no vale 5 kcal/mol. Y ese resultado es el unico «positivo» que el
proyecto tenia hasta hoy: en la lista ordenada, 4MQ y ZO0 entran en el 1 %
superior y son dos quimias independientes por delante del fondo.

DOS COSAS QUE SE MIDE AQUI
--------------------------
1. REPRODUCIBILIDAD. La validacion v4 acoplo sin `--seed`, o sea con la semilla
   aleatoria de Vina. Si el -10 sale de una sola corrida estocastica, no es un
   dato: es un numero. Se repiten los tres ligandos con tres semillas fijas
   (42, 2026, 777) y los mismos parametros (caja 46,5 / 80,0 / 73,3; 22 A; ex 8).

2. CONTACTOS. Se comparan las poses ya guardadas de los tres (carpeta de la v4):
   que residuos toca cada una, cuantos pares de atomos hay a menos de 4 A, y
   cuanto queda enterrado cada atomo del ligando. Si los dos diazepanos comparten
   un contacto que la piperazina no hace, ahi esta la respuesta.

3. EL RECEPTOR. Y aqui aparece lo que puede ser la causa de todo: SOD1.pdbqt
   NO es una proteina limpia. Lleva dentro las 1.763 AGUAS CRISTALOGRAFICAS de
   1HL5 (tipo OA en el PDBQT) y los iones de zinc, y ademas 18 copias de la
   proteina (todo el ensamblaje del fichero depositado). Cinco de esas aguas
   caen DENTRO del bolsillo de Trp32, a menos de 10 A del centro de la caja.
   Vina no puede desplazar una molecula de agua: la trata como parte rigida del
   receptor. Un ligando que se aprieta ENTRE la proteina y las aguas congeladas
   cobra terminos favorables que en la realidad no existirian (en la realidad
   esas aguas se van). Se repite el acoplamiento con el mismo receptor pero sin
   aguas para ver si los -10 sobreviven.

Uso: python anomalia_quinazolinas.py [--sin-acoplar]
Salida: anomalia_quinazolinas.log y anomalia_quinazolinas.csv
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
GPU = os.path.join(BASE, "..", "gpu_dock")
RECEPTOR = os.path.abspath(os.path.join(GPU, "SOD1.pdbqt"))
VINA = r"C:\Users\Fredy\masive-als\tools\vina.exe"
CENTRO = (46.5, 80.0, 73.3)
TAMANO = 22
EXHAUSTIVIDAD = 8
SEMILLAS = (42, 2026, 777)

LIGANDOS = {
    "4MQ": ("diazepanoquinazolina", "validacion_SOD1_v4/out/ACT_diazepanoquinazolina_4MQ_out.pdbqt"),
    "ZO0": ("cf3-diazepanoquinazolina", "validacion_SOD1_v4/out/ACT_cf3quinazolina_ZO0_out.pdbqt"),
    "12I": ("piperazinaquinazolina", "validacion_SOD1_v4/out/ACT_quinazolina_12I_out.pdbqt"),
}
ENTRADAS = {
    "4MQ": "validacion_SOD1_v4/ligands/ACT_diazepanoquinazolina_4MQ.pdbqt",
    "ZO0": "validacion_SOD1_v4/ligands/ACT_cf3quinazolina_ZO0.pdbqt",
    "12I": "validacion_SOD1_v4/ligands/ACT_quinazolina_12I.pdbqt",
}
NUCLEO = "quinazolina"


def log(m):
    print(m, flush=True)


def atomos_receptor():
    """Atomos pesados del receptor: (nombre_residuo, cadena, resSeq, xyz)."""
    out = []
    for l in open(RECEPTOR, encoding="utf-8", errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        tipo = l.rsplit(None, 1)[-1].upper()
        if tipo in ("H", "HD", "HS", "D", "DD"):
            continue
        out.append((l[17:20].strip() + l[22:27].strip(), l[21], l[22:27].strip(),
                    np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


def primera_pose(path):
    """Coordenadas de los atomos pesados del primer modo de un PDBQT."""
    out = []
    for l in open(path, encoding="utf-8", errors="ignore"):
        if out and l.startswith("ENDMDL"):
            break
        if not l.startswith(("ATOM", "HETATM")):
            continue
        tipo = l.rsplit(None, 1)[-1].upper()
        if tipo in ("H", "HD", "HS", "D", "DD"):
            continue
        out.append(np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
    return np.array(out)


def afinidad(path):
    for l in open(path, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            return float(re.search(r"([-+]?\d+\.\d+)", l).group(1))
    return None


def contactos(pose, coords, nombres, corte_res=4.5, corte_par=4.0):
    """Residuos tocados, pares de atomos a menos de 4 A y atomos enterrados."""
    por_res = {}
    pares = 0
    for x in pose:
        d = np.sqrt(((coords - x) ** 2).sum(-1))
        for i in np.where(d <= corte_res)[0]:
            clave = nombres[i]
            r = por_res.setdefault(clave, {"n_atomos": 0, "dmin": 99.0})
            r["n_atomos"] += 1
            r["dmin"] = min(r["dmin"], float(d[i]))
        pares += int((d <= corte_par).sum())
    enterrados = sum(1 for x in pose
                     if (np.sqrt(((coords - x) ** 2).sum(-1)) <= corte_res).sum() >= 3)
    return por_res, pares, enterrados


def receptor_sin_disolvente(origen, destino):
    """Copia el receptor quitando las aguas cristalograficas.

    Los iones (ZN, CU) se conservan: son parte de la estructura de la proteina.
    Lo que se quita es el disolvente, que en un acoplamiento no se puede
    desplazar y por tanto no es un companero de union legitimo.
    """
    n = 0
    with open(destino, "w", encoding="utf-8") as f:
        for l in open(origen, encoding="utf-8", errors="ignore"):
            if l.startswith(("ATOM", "HETATM")):
                if l[17:20].strip() in ("HOH", "DOD", "WAT"):
                    continue
                f.write(l)
                n += 1
            elif l.startswith("REMARK"):
                f.write(l)
    return n


def acoplar(ligando, receptor=RECEPTOR, etiqueta=""):
    """Repite el acoplamiento con semillas fijas (la v4 no uso semilla)."""
    ent = os.path.join(BASE, ENTRADAS[ligando])
    if not os.path.exists(ent):
        log("  falta el ligando de entrada %s" % ent)
        return []
    os.makedirs(os.path.join(BASE, "out"), exist_ok=True)
    salidas = []
    for semilla in SEMILLAS:
        out = os.path.join(BASE, "out", "anomalia_%s%s_s%d.pdbqt"
                           % (ligando, etiqueta, semilla))
        if not (os.path.exists(out) and os.path.getsize(out) > 100):
            r = subprocess.run([VINA, "--receptor", receptor, "--ligand", ent,
                            "--center_x", str(CENTRO[0]), "--center_y", str(CENTRO[1]),
                            "--center_z", str(CENTRO[2]),
                            "--size_x", str(TAMANO), "--size_y", str(TAMANO),
                            "--size_z", str(TAMANO),
                            "--exhaustiveness", str(EXHAUSTIVIDAD),
                            "--num_modes", "9", "--seed", str(semilla),
                            "--out", out],
                           capture_output=True, text=True, timeout=3600)
        if not (os.path.exists(out) and os.path.getsize(out) > 100):
            log("  FALLO acoplando %s (semilla %d): %s"
                % (ligando, semilla, str(getattr(r, "stderr", ""))[-200:]))
            continue
        salidas.append((semilla, out))
    return salidas


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sin-acoplar", action="store_true",
                    help="no repetir el acoplamiento; solo analizar lo guardado")
    args = ap.parse_args()

    info = atomos_receptor()
    coords = np.array([a[3] for a in info], dtype=float)
    nombres = [a[0] for a in info]
    log("receptor: %d atomos pesados de %s" % (len(info), os.path.basename(RECEPTOR)))

    # receptor sin aguas: la prueba de si los -10 son de la proteina o del disolvente
    sin_agua = os.path.join(BASE, "receptores", "SOD1_sin_agua.pdbqt")
    os.makedirs(os.path.join(BASE, "receptores"), exist_ok=True)
    if not os.path.exists(sin_agua):
        n = receptor_sin_disolvente(RECEPTOR, sin_agua)
        log("receptor sin aguas escrito: %d atomos (de %d) -> %s"
            % (n, len(info), sin_agua))

    lineas, filas = [], []
    for lig, (nombre, ruta_v4) in LIGANDOS.items():
        log("\n=== %s (%s) ===  caja %.1f,%.1f,%.1f  %d A  ex %d"
            % (lig, nombre, *CENTRO, TAMANO, EXHAUSTIVIDAD))
        if not args.sin_acoplar:
            salidas = acoplar(lig)
            afs = [afinidad(o) for _, o in salidas]
            log("  con semilla fija: %s"
                % ", ".join("s%d %.2f" % (s, a) for (s, _), a in zip(salidas, afs)))
            if afs:
                log("  -> mejor %.2f | dispersion %.2f kcal/mol entre semillas"
                    % (min(afs), max(afs) - min(afs)))

        if not args.sin_acoplar:
            afs_sa = [afinidad(o) for _, o in acoplar(lig, sin_agua, "_sinagua")]
            if afs_sa:
                log("  SIN AGUAS en el receptor: %s | mejor %.2f"
                    % (", ".join("%.2f" % a for a in afs_sa), min(afs_sa)))

        for etiqueta, ruta in (("v4 (sin semilla)", os.path.join(BASE, ruta_v4)),
                               ("semilla 42", os.path.join(BASE, "out",
                                                           "anomalia_%s_s42.pdbqt" % lig)),
                               ("sin aguas", os.path.join(BASE, "out",
                                                          "anomalia_%s_sinagua_s42.pdbqt" % lig))):
            if not os.path.exists(ruta):
                continue
            pose = primera_pose(ruta)
            if not len(pose):
                continue
            por_res, pares, enterrados = contactos(pose, coords, nombres)
            log("  [%s] afinidad %.2f | %d atomos pesados | %d pares <4 A | "
                "%d residuos a <4,5 A | %d/%d atomos enterrados"
                % (etiqueta, afinidad(ruta), len(pose), pares, len(por_res),
                   enterrados, len(pose)))
            detalle = sorted(por_res.items(), key=lambda kv: -kv[1]["n_atomos"])
            log("      %s" % "  ".join("%s(%d)" % (k, v["n_atomos"]) for k, v in detalle[:12]))
            filas.append({"ligando": lig, "nombre": nombre, "origen": etiqueta,
                          "afinidad": afinidad(ruta), "pares_4A": pares,
                          "residuos": len(por_res), "enterrados": enterrados,
                          "pesados": len(pose),
                          "residuos_detalle": " ".join("%s:%d" % (k, v["n_atomos"])
                                                       for k, v in detalle)})
            lineas.append("%-4s %-22s %-16s afinidad %7.2f | pares<4A %4d | "
                          "residuos %2d | enterrados %2d/%2d"
                          % (lig, nombre, etiqueta, afinidad(ruta), pares,
                             len(por_res), enterrados, len(pose)))

    texto = "\n".join(lineas)
    with open(os.path.join(BASE, "anomalia_quinazolinas.log"), "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    campos = ["ligando", "nombre", "origen", "afinidad", "pares_4A", "residuos",
              "enterrados", "pesados", "residuos_detalle"]
    with open(os.path.join(BASE, "anomalia_quinazolinas.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for fila in filas:
            w.writerow({**{c: "" for c in campos}, **fila})
    log("\nguardado: anomalia_quinazolinas.csv / .log")
    return 0


if __name__ == "__main__":
    sys.exit(main())
