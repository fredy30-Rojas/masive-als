#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Los ocho en el CR con el Trp334 FLEXIBLE: aparece la senal que el rigido no da?

LA PREGUNTA
-----------
El acoplamiento rigido (`acoplar_familia_cr.py`) no distingue a los dos que inhiben la
agregacion (XL21, XL23) de los tres que no tienen nada reportado (XL22, XL25, XL26):
las distribuciones se solapan y el residuo tras quitar el tamano va al reves. Una
explicacion posible es que **el receptor este de mas**: la cara de la helice se calcula
con el anillo del Trp334 clavado en su sitio, y el propio articulo dice que ese residuo
importa (mutarlo baja la union). Si el indol puede girar para acomodar al ligando, puede
aparecer la diferencia que con el receptor quieto no se ve.

Este script hace exactamente eso y nada mas: mismo receptor, misma caja, mismas tres
semillas, mismos parametros; **lo unico que cambia es que la cadena lateral del Trp334
se mueve** (meeko `-f A:334`). Los resultados del rigido se leen de `familia_cr.csv`,
no se recalculan, para que la comparacion sea con lo que ya esta publicado.

QUE MIDE, Y POR QUE DE ESA MANERA
---------------------------------
1. Acoplamiento flexible ligando + cadena lateral con Vina (`--receptor <rigid>
   --flex <flex>`), que es el mismo esquema que uso el proyecto para el control de
   receptor flexible del Trp32 de SOD1 (`redocking_trp32/redock_flexible.py`).
2. **Los contactos se miden contra la posicion MOVIDA del indol**, que es la que Vina
   devuelve en cada modo: si el residuo se ha movido, la distancia del ligando a su
   anillo no puede calcularse con la de la estructura de partida.
3. La misma comparacion por grupos y el mismo ajuste por tamano que el rigido, con las
   funciones del script de la familia (una sola copia del analisis).

LO QUE HAY QUE TENER PRESENTE AL LEERLO
---------------------------------------
Soltar la cadena lateral **sube todas las afinidades** (hay mas libertad para
acomodarse), asi que las puntuaciones flexibles no se comparan con las rigidas una a
una: lo que se compara es **la forma de la tabla** (quien va delante de quien, y si los
grupos se separan). Y un residuo flexible tampoco arregla lo de fondo: si el CR es una
helice de una micela, mover el Trp334 sigue siendo mover un residuo de esa helice.

Uso:
    python analysis/acoplar_familia_cr_flexible.py [--exhaustividad 8] [--lado 22]
                                                   [--semillas 42,2026,777]

Salida en analysis/_cr_receptor/:
    modelos/flex334_modelo<n>_rigid.pdbqt, _flex.pdbqt   receptores por modelo
    out/<ligando>_modelo<n>_s<s>_flex334.pdbqt            las poses
    familia_cr_flexible.csv                               una fila por pose
    familia_cr_flexible.txt                               el informe, con el contraste
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import statistics
import subprocess
import sys
import time

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)
sys.path.insert(0, os.path.dirname(BASE))
import acoplar_xl20_xl23_cr as P  # noqa: E402
import acoplar_familia_cr as F  # noqa: E402
from preparar_ligando import escribir as escribir_ligando  # noqa: E402

SALIDA = P.SALIDA
MODELOS = P.MODELOS
LIGANDS = P.LIGANDS
OUT = P.OUT
RESIDUO_FLEXIBLE = P.TRP          # el Trp334: el residuo al que el articulo senala


def log(m):
    print(m, flush=True)


def preparar_flexible(clave, centro, lado):
    """(rigid.pdbqt, flex.pdbqt) del modelo, con la cadena del Trp334 suelta."""
    crudo = os.path.join(MODELOS, "cr_modelo%s.pdb" % clave)
    base = os.path.join(MODELOS, "flex%d_modelo%s" % (RESIDUO_FLEXIBLE, clave))
    rigid, flex = base + "_rigid.pdbqt", base + "_flex.pdbqt"
    if os.path.exists(rigid) and os.path.exists(flex):
        return rigid, flex
    cmd = [P.R.PYTHON, P.R.MK_RECEPTOR, "--read_pdb", crudo, "-o", base,
           "-p", "-j", "-f", "A:%d" % RESIDUO_FLEXIBLE, "--default_altloc", "A",
           "--box_center", "%.3f" % centro[0], "%.3f" % centro[1], "%.3f" % centro[2],
           "--box_size", str(lado), str(lado), str(lado)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    if not (os.path.exists(rigid) and os.path.exists(flex)):
        log("  modelo %s: FALLO el receptor flexible: %s"
            % (clave, ((r.stdout or "") + (r.stderr or ""))[-300:]))
        return None, None
    n_flex = sum(1 for l in open(flex, encoding="utf-8", errors="ignore")
                 if l.startswith("ATOM"))
    log("  modelo %s: receptor flexible listo (%d atomos sueltos del Trp%d)"
        % (clave, n_flex, RESIDUO_FLEXIBLE))
    return rigid, flex


def poses_flexibles(ruta):
    """[(modo, afinidad, [atomo del ligando], [atomo del residuo movido])].

    En la salida de un acoplamiento flexible, cada MODEL trae DOS cosas: el ligando
    (residuo UNL) y la conformacion del residuo que se ha soltado (TRP A 334). Hay que
    separarlas, porque si se mezclan los contactos salen mal: el ligando "tocaria" sus
    propios atomos y el residuo movido contaria como receptor dos veces.
    """
    poses = []
    actual = None
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            actual = {"af": None, "lig": [], "flex": []}
            continue
        if l.startswith("ENDMDL"):
            if actual is not None:
                poses.append((len(poses) + 1, actual["af"], actual["lig"], actual["flex"]))
            actual = None
            continue
        if actual is None:
            continue
        if l.startswith("REMARK VINA RESULT:"):
            m = re.search(r"(-?\d+\.\d+)", l)
            actual["af"] = float(m.group(1)) if m else None
            continue
        if not l.startswith(("ATOM", "HETATM")):
            continue
        nombre = l[12:16].strip()
        elemento = (l[76:78].strip() or nombre[0]).upper()
        if elemento == "H" or nombre.startswith("H"):
            continue
        xyz = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        if l[17:20].strip() == "TRP":
            actual["flex"].append((nombre, xyz))
        else:
            actual["lig"].append((nombre, xyz))
    return poses


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exhaustividad", type=int, default=8)
    ap.add_argument("--semillas", default="42,2026,777")
    ap.add_argument("--lado", type=int, default=22)
    args = ap.parse_args()
    semillas = [int(s) for s in args.semillas.split(",") if s.strip()]

    os.makedirs(LIGANDS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(SALIDA, "caja.json"), encoding="utf-8") as f:
        caja = json.load(f)
    modelos = sorted(caja["modelos"], key=int)

    evidencia, pesados, pdbqt = {}, {}, {}
    for nombre, (archivo, grupo) in F.FUENTES.items():
        smi, texto = F.leer(archivo, nombre)
        ruta = escribir_ligando(nombre, smi, LIGANDS)
        if ruta is None:
            log("%s NO se pudo preparar" % nombre)
            return 1
        pdbqt[nombre] = ruta
        evidencia[nombre] = (grupo, texto)
        pesados[nombre] = len(P.atomos(ruta))

    log("Los %d compuestos con el Trp%d flexible (mismo receptor, caja y semillas que"
        " el rigido):" % (len(pdbqt), RESIDUO_FLEXIBLE))
    log("")

    receptores = {}
    for clave in modelos:
        info = caja["modelos"][clave]
        rigid, flex = preparar_flexible(clave, info["centro"], args.lado)
        if rigid:
            receptores[clave] = (rigid, flex, info)

    filas = []
    t0 = time.time()
    for nombre, ligando in pdbqt.items():
        for clave in modelos:
            if clave not in receptores:
                continue
            rigid, flex, info = receptores[clave]
            info_flex = info
            # el receptor quieto de este modelo: sus atomos pesados y sus residuos
            rec = P.atomos(rigid)
            rec_xyz = np.array([a[3] for a in rec])
            cx, cy, cz = info_flex["centro"]
            for semilla in semillas:
                out = os.path.join(OUT, "%s_modelo%s_s%d_flex%d.pdbqt"
                                   % (nombre, clave, semilla, RESIDUO_FLEXIBLE))
                if not (os.path.exists(out) and os.path.getsize(out) > 100):
                    r = subprocess.run(
                        [P.R.VINA, "--receptor", rigid, "--flex", flex,
                         "--ligand", ligando,
                         "--center_x", "%.3f" % cx, "--center_y", "%.3f" % cy,
                         "--center_z", "%.3f" % cz,
                         "--size_x", str(args.lado), "--size_y", str(args.lado),
                         "--size_z", str(args.lado),
                         "--exhaustiveness", str(args.exhaustividad),
                         "--num_modes", "9", "--seed", str(semilla), "--out", out],
                        capture_output=True, text=True, timeout=3600)
                    if not (os.path.exists(out) and os.path.getsize(out) > 100):
                        log("  %s modelo %s semilla %d: FALLO Vina: %s"
                            % (nombre, clave, semilla, (r.stderr or r.stdout)[-160:]))
                        continue
                    log("  %-5s modelo %-2s semilla %-5d ok (%.0f s acumulados)"
                        % (nombre, clave, semilla, time.time() - t0))
                for modo, afinidad, lig, flex_xyz in poses_flexibles(out):
                    L = np.array([a[1] for a in lig])
                    # contactos: la parte quieta del receptor MAS el indol en la
                    # posicion que Vina le ha dado en esta pose.
                    todos = rec_xyz
                    if flex_xyz:
                        todos = np.vstack([rec_xyz, np.array([a[1] for a in flex_xyz])])
                    d = np.sqrt(((L[:, None, :] - todos[None, :, :]) ** 2).sum(-1))
                    residuos = sorted({rec[i][0] for i
                                       in np.nonzero(d.min(axis=0) < P.CORTE)[0]
                                       if i < len(rec)})
                    hay_flex = bool(d[:, len(rec):].size
                                    and (d[:, len(rec):].min(axis=0) < P.CORTE).any())
                    if hay_flex:
                        residuos = sorted(set(residuos) | {RESIDUO_FLEXIBLE})
                    en_contacto = int((d.min(axis=1) < P.CORTE).sum())
                    anillo = np.array([xyz for (nom, xyz) in flex_xyz
                                       if nom in P.ANILLO_TRP])
                    dist_atomo = (float(np.linalg.norm(
                        L[:, None, :] - anillo[None, :, :], axis=-1).min())
                        if len(anillo) else None)
                    dist_centro = (float(np.linalg.norm(L - anillo.mean(axis=0), axis=1).min())
                                   if len(anillo) else None)
                    filas.append({
                        "ligando": nombre, "grupo": evidencia[nombre][0],
                        "evidencia_funcional": evidencia[nombre][1],
                        "atomos_pesados": pesados[nombre],
                        "residuo_flexible": RESIDUO_FLEXIBLE,
                        "modelo": clave, "semilla": semilla, "modo": modo,
                        "afinidad": afinidad,
                        "contactos": len(residuos),
                        "atomos_del_ligando_en_contacto": en_contacto,
                        "fraccion_ligando_en_contacto": round(en_contacto / float(len(L)), 2),
                        "residuos": " ".join(str(r) for r in residuos),
                        "toca_trp334": RESIDUO_FLEXIBLE in residuos,
                        "dist_min_al_anillo_trp334": (round(dist_centro, 2)
                                                      if dist_centro is not None else None),
                        "dist_min_a_un_atomo_del_anillo": (round(dist_atomo, 2)
                                                           if dist_atomo is not None else None),
                    })

    campos = list(filas[0].keys())
    with open(os.path.join(SALIDA, "familia_cr_flexible.csv"), "w", encoding="utf-8",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(filas)

    cabecera = [
        "Los ocho compuestos de la tabla del CR, acoplados con el Trp%d FLEXIBLE"
        % RESIDUO_FLEXIBLE,
        "receptor %s (%s), %d modelos de RMN, caja de %d A; cadena lateral del Trp%d"
        " suelta (meeko -f)" % (caja["pdb"], caja["residuos"], len(modelos), args.lado,
                                RESIDUO_FLEXIBLE),
        "exhaustividad %d, 9 modos, semillas %s | contactos a %.1f A medidos contra la"
        " posicion MOVIDA del indol" % (args.exhaustividad,
                                        ", ".join(str(s) for s in semillas), P.CORTE),
        "",
    ]
    texto = F.informe(filas, modelos, semillas, pesados, cabecera)

    # --- el contraste con el rigido, que es a lo que se ha venido ---
    rigido_csv = os.path.join(SALIDA, "familia_cr.csv")
    contraste = []
    if os.path.exists(rigido_csv):
        with open(rigido_csv, encoding="utf-8", newline="") as f:
            rigidas = list(csv.DictReader(f))
        rig_med = F.mediana_por_compuesto(rigidas, modelos, semillas)
        flex_med = F.mediana_por_compuesto(filas, modelos, semillas)
        rsol, rpeor, rmejor = F.solape(rig_med["XL21"] + rig_med["XL23"],
                                       rig_med["XL22"] + rig_med["XL25"] + rig_med["XL26"])
        fsol, fpeor, fmejor = F.solape(flex_med["XL21"] + flex_med["XL23"],
                                       flex_med["XL22"] + flex_med["XL25"] + flex_med["XL26"])
        _, rres = F.residuos_de_tamano({n: statistics.median(v)
                                        for n, v in rig_med.items()}, pesados)
        _, fres = F.residuos_de_tamano({n: statistics.median(v)
                                        for n, v in flex_med.items()}, pesados)

        def mediana(v):
            return statistics.median(v)

        contraste.append("")
        contraste.append("  RIGIDO FRENTE A FLEXIBLE (mediana de las %d combinaciones de"
                         " cada compuesto)" % (len(modelos) * len(semillas)))
        contraste.append("  %-5s %-9s %-9s %-8s %-13s %-13s"
                         % ("lig", "rigido", "flexible", "cambio", "Trp334 rigido",
                            "Trp334 flexible"))
        for nombre in F.FUENTES:
            s_rig = sum(1 for r in rigidas if r["ligando"] == nombre
                        and r["dist_min_a_un_atomo_del_anillo"] not in ("", None)
                        and float(r["dist_min_a_un_atomo_del_anillo"]) < P.CORTE)
            t_rig = sum(1 for r in rigidas if r["ligando"] == nombre)
            s_flex = sum(1 for r in filas if r["ligando"] == nombre
                         and (r["dist_min_a_un_atomo_del_anillo"] or 99) < P.CORTE)
            t_flex = sum(1 for r in filas if r["ligando"] == nombre)
            contraste.append("  %-5s %-9.2f %-9.2f %-8.2f %-13s %-13s"
                             % (nombre, mediana(rig_med[nombre]), mediana(flex_med[nombre]),
                                mediana(flex_med[nombre]) - mediana(rig_med[nombre]),
                                "%d de %d" % (s_rig, t_rig),
                                "%d de %d" % (s_flex, t_flex)))
        contraste.append("")
        contraste.append("  LO QUE CAMBIA CON LA CADENA SUELTA:")
        contraste.append("    solape rigido:   %s (peor inhibidor %.2f, mejor sin nada %.2f)"
                         % ("SI, se solapan" if rsol else "NO, se separan", rpeor, rmejor))
        contraste.append("    solape flexible: %s (peor inhibidor %.2f, mejor sin nada %.2f)"
                         % ("SI, se solapan" if fsol else "NO, se separan", fpeor, fmejor))
        contraste.append("    por residuo (tras quitar el tamano) -- los tres mejores:")
        contraste.append("       rigido:   %s"
                         % ", ".join(sorted(rres, key=lambda n: rres[n])[:3]))
        contraste.append("       flexible: %s"
                         % ", ".join(sorted(fres, key=lambda n: fres[n])[:3]))
        contraste.append("    (hay que leerlo con esto en la mano: soltar la cadena sube"
                         " TODAS las")
        contraste.append("     afinidades, asi que lo que se compara es la forma de la"
                         " tabla, no los")
        contraste.append("     numeros uno a uno; y el Trp334 sigue siendo un residuo de"
                         " esta helice)")
    texto = texto + "\n" + "\n".join(contraste)
    with open(os.path.join(SALIDA, "familia_cr_flexible.txt"), "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    log("")
    log("\n".join(contraste) if contraste else "(no hay corrida rigida con la que comparar)")
    log("")
    log("guardado: analysis/_cr_receptor/familia_cr_flexible.txt")
    return 0


if __name__ == "__main__":
    sys.exit(main())
