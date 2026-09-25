#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tercer control de redocking: receptor FLEXIBLE en el bolsillo Trp32 (20 sep 2026).

POR QUE EXISTE
--------------
El segundo control (`protocolo_alternativo.py`) partio el fallo en dos:

- isoproterenol y adrenalina: con el ligando RIGIDO el acoplamiento reproduce la
  pose cristalina a 0,33-0,65 A, y con el ligando flexible se va a 3,55-3,98 A
  por solo 0,4-0,5 kcal/mol. El problema medido es de torsionos del ligando.
- dopamina y 5-fluorouridina: ni rigidas se reproducen (2,97 y 9,8-10,0 A).

Quedaba sin probar la otra mitad del muestreo: la CADENA LATERAL del receptor.
Todos los controles anteriores usaron receptor rigido. Aqui se sueltan los
residuos que definen el bolsillo (los mismos que el primer informe identifico
como contactos: Glu21, Gln22, Lys23, Lys30 y Glu100; se puede anadir Trp32 con
--residuos) y se repite el mismo control, con las mismas cajas, las mismas
semillas y el mismo criterio de RMSD.

Pasos:
  1. mk_prepare_receptor -f  ->  <base>_rigid.pdbqt + <base>_flex.pdbqt
  2. acoplamiento flexible (ligando y residuos) con vina y/o vinardo
  3. RMSD sin superponer contra la pose del cristal (criterio <= 2,0 A)

Uso:
    python redock_flexible.py --entradas 4A7T --semillas 42      # prueba de tiempo
    python redock_flexible.py                                    # las cuatro, 3 semillas
"""
import argparse
import json
import os
import re
import subprocess
import sys
import time

import redock_trp32 as R
import protocolo_alternativo as P

BASE = os.path.dirname(os.path.abspath(__file__))
PYDIR = R.PYDIR
PYTHON = R.PYTHON
MK_RECEPTOR = R.MK_RECEPTOR
VINA = R.VINA
SISTEMAS = R.SISTEMAS
CAJA = 24

# residuos del bolsillo de Trp32 que se dejan flexibles (cadena A de I113T)
RESIDUOS_FLEX = "21,22,23,30,100"


def log(m):
    print(m, flush=True)


def rmsd_alineado(cry_pdb, out_pdbqt):
    """RMSD del mejor modo DESPUES de superponerlo al cristal (Kabsch).

    No sustituye al RMSD en el sitio (el criterio del control es el de siempre:
    ambas poses ya estan en el sistema del receptor). Sirve para caracterizar el
    fallo: con receptor flexible aparecen desplazamientos de 12 A en el sitio
    pero de solo 4,5 A de centroide, y hay que poder decir si el ligando se fue
    del bolsillo o si se quedo cerca pero girado.
    """
    import numpy as np
    from rdkit import Chem
    cry = Chem.MolFromPDBBlock(open(cry_pdb).read(), removeHs=False,
                               proximityBonding=True)
    pose = R.molecula_dockeada(out_pdbqt)
    if cry is None or pose is None:
        return None
    p = R.aplanar(Chem.RemoveHs(Chem.Mol(pose)))
    c = R.aplanar(Chem.RemoveHs(cry))
    if p.GetNumAtoms() != c.GetNumAtoms():
        return None
    P3 = p.GetConformer().GetPositions()
    C3 = c.GetConformer().GetPositions()
    Cc = C3 - C3.mean(0)
    mejor = None
    for match in p.GetSubstructMatches(c, uniquify=False, maxMatches=500,
                                       useChirality=False):
        A = P3[list(match)]
        Ac = A - A.mean(0)
        U, S, Vt = np.linalg.svd(Ac.T @ Cc)
        D = np.diag([1.0, 1.0, float(np.sign(np.linalg.det(Vt.T @ U.T)))])
        Ar = Ac @ (Vt.T @ D @ U.T).T
        v = float(np.sqrt(((Ar - Cc) ** 2).sum() / len(A)))
        if mejor is None or v < mejor:
            mejor = v
    return mejor


def dist_centroides(a, b):
    """Distancia entre los centroides de dos ficheros de ligando (atomo pesado).

    Sirve para distinguir «acoplo mal dentro del bolsillo» de «el ligando se
    marcho del bolsillo»: con 12 A de RMSD la diferencia importa.
    """
    A, B = P.coord_pesados(a), P.coord_pesados(b)
    if not A or not B:
        return None
    ca = [sum(c[i] for c in A) / len(A) for i in range(3)]
    cb = [sum(c[i] for c in B) / len(B) for i in range(3)]
    return sum((ca[i] - cb[i]) ** 2 for i in range(3)) ** 0.5


def preparar_receptor_flex(entrada, residuos, caja, cx, cy, cz):
    """Devuelve (rigid.pdbqt, flex.pdbqt) con los residuos pedidos flexibles."""
    rec_pdb = os.path.join(BASE, "receptores", "%s_proteina.pdb" % entrada)
    if not os.path.exists(rec_pdb):
        return None, None
    etiqueta = residuos.replace(",", "-")
    base = os.path.join(BASE, "receptores", "%s_flex_%s" % (entrada, etiqueta))
    rigid, flex = base + "_rigid.pdbqt", base + "_flex.pdbqt"
    if os.path.exists(rigid) and os.path.exists(flex):
        return rigid, flex
    # --default_altloc A: estas estructuras traen atomos en dos posiciones
    # alternativas y meeko se niega a adivinar cual usar (el primer control ya
    # usaba la misma opcion, por eso alli si funciono).
    cmd = [PYTHON, MK_RECEPTOR, "--read_pdb", rec_pdb, "-o", base,
           "-p", "-j", "-f", "A:%s" % residuos, "--default_altloc", "A",
           "--box_center", "%.3f" % cx, "%.3f" % cy, "%.3f" % cz,
           "--box_size", str(caja), str(caja), str(caja)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
    if not (os.path.exists(rigid) and os.path.exists(flex)):
        log("  FALLO preparando el receptor flexible: %s"
            % ((r.stdout or "") + (r.stderr or ""))[-400:])
        return None, None
    n_flex = sum(1 for l in open(flex, encoding="utf-8", errors="ignore")
                 if l.startswith("ATOM"))
    log("  receptor flexible listo: %d atomos en la parte flexible" % n_flex)
    return rigid, flex


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entradas", default=",".join(SISTEMAS))
    ap.add_argument("--puntuaciones", default="vina")
    ap.add_argument("--residuos", default=RESIDUOS_FLEX,
                    help="residuos flexibles de la cadena A, separados por comas")
    ap.add_argument("--exhaustividad", type=int, default=8)
    ap.add_argument("--caja", type=int, default=CAJA)
    ap.add_argument("--semillas", default="42,2026,777")
    args = ap.parse_args()

    semillas = [int(s) for s in args.semillas.split(",") if s.strip()]
    resultados = []
    for entrada in [e.strip() for e in args.entradas.split(",") if e.strip()]:
        codigo, cadena, nombre = SISTEMAS[entrada]
        cry_pdb = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
        sdf = os.path.join(BASE, "ligands", "%s_ideal.sdf" % codigo)
        cry_pdbqt = os.path.join(BASE, "ligands", "%s_cristal_flexp.pdbqt" % entrada)
        if not os.path.exists(cry_pdbqt) and not P.pdbqt_del_cristal(
                cry_pdb, sdf, cry_pdbqt):
            log("=== %s: no se pudo preparar la pose del cristal" % entrada)
            continue
        cx, cy, cz = R.centroide(cry_pdb)
        log("\n=== %s (%s) ===" % (entrada, nombre))
        rigid, flex = preparar_receptor_flex(entrada, args.residuos, args.caja,
                                             cx, cy, cz)
        if not rigid:
            continue

        for puntuacion in [p.strip() for p in args.puntuaciones.split(",") if p.strip()]:
            # 1. puntuacion de la pose del cristal, con los residuos en su sitio
            txt = P.correr_vina(["--receptor", rigid, "--flex", flex,
                                 "--ligand", cry_pdbqt, "--scoring", puntuacion,
                                 "--autobox", "--score_only"], timeout=3600)
            e_cristal = P.afinidad_stdout(txt)
            log("  [%s] E(cristal, residuos cristalinos) = %s kcal/mol"
                % (puntuacion,
                   "%.2f" % e_cristal if e_cristal is not None else "no leida"))
            if e_cristal is None:
                log("     ultimas lineas: %s" % str(txt.strip().splitlines()[-3:]))

            por_semilla = []
            for semilla in semillas:
                out = os.path.join(BASE, "out", "%s_flex%s_%s_caja%d_e%d_s%d.pdbqt"
                                   % (entrada, args.residuos.replace(",", "-"),
                                      puntuacion, args.caja,
                                      args.exhaustividad, semilla))
                t0 = time.time()
                r = subprocess.run(
                    [VINA, "--receptor", rigid, "--flex", flex,
                     "--ligand", cry_pdbqt,
                     "--center_x", "%.3f" % cx, "--center_y", "%.3f" % cy,
                     "--center_z", "%.3f" % cz,
                     "--size_x", str(args.caja), "--size_y", str(args.caja),
                     "--size_z", str(args.caja),
                     "--scoring", puntuacion,
                     "--exhaustiveness", str(args.exhaustividad),
                     "--num_modes", "9", "--seed", str(semilla), "--out", out],
                    capture_output=True, text=True, timeout=21600)
                dt = time.time() - t0
                af = P.afinidad_de_salida(out)
                val, dist = None, None
                if os.path.exists(out):
                    try:
                        val = R.rmsd_en_sitio(cry_pdb, R.molecula_dockeada(out))
                    except Exception as e:
                        log("  AVISO: no se pudo medir el RMSD: %s" % e)
                    dist = dist_centroides(cry_pdbqt, out)
                log("  [%s] semilla %-5d afinidad %s | RMSD vs cristal %s | "
                    "centroide a %s | %.0f s"
                    % (puntuacion, semilla,
                       "%.2f" % af if af is not None else "no leida",
                       "%.2f A" % val if val is not None else "no calculable",
                       "%.2f A" % dist if dist is not None else "n/d", dt))
                por_semilla.append({"semilla": semilla, "afinidad": af,
                                    "rmsd": round(val, 2) if val is not None else None,
                                    "dist_centroide": round(dist, 2) if dist is not None else None,
                                    "segundos": round(dt, 1)})

            rmsds = [x["rmsd"] for x in por_semilla if x["rmsd"] is not None]
            mejor = min(rmsds) if rmsds else None
            mediana = sorted(rmsds)[len(rmsds) // 2] if rmsds else None
            mejor_af = min([x["afinidad"] for x in por_semilla
                            if x["afinidad"] is not None], default=None)
            de = (e_cristal - mejor_af) if (e_cristal is not None
                                            and mejor_af is not None) else None
            log("  [%s] -> mejor %s | mediana %s | pasa %d/%d | dE %s"
                % (puntuacion,
                   "%.2f A" % mejor if mejor is not None else "n/d",
                   "%.2f A" % mediana if mediana is not None else "n/d",
                   sum(1 for v in rmsds if v <= 2.0), len(rmsds),
                   "%+.2f kcal/mol" % de if de is not None else "n/d"))
            resultados.append({
                "entrada": entrada, "ligando": codigo, "nombre": nombre,
                "puntuacion": puntuacion, "residuos_flexibles": args.residuos,
                "caja": args.caja, "exhaustividad": args.exhaustividad,
                "E_cristal_quieta": e_cristal,
                "E_mejor_modo": mejor_af, "delta_E": round(de, 2) if de is not None else None,
                "por_semilla": por_semilla,
                "rmsd_mejor": mejor, "rmsd_mediana": mediana,
                "semillas_que_pasan": sum(1 for v in rmsds if v <= 2.0),
                "criterio": "PASA" if (mediana is not None and mediana <= 2.0) else "FALLA",
            })

    destino = os.path.join(BASE, "resultado_redocking_flexible_%s_%s.json"
                           % (args.residuos.replace(",", "-"),
                              args.entradas.replace(",", "-")))
    with open(destino, "w", encoding="utf-8") as f:
        json.dump(resultados, f, indent=2, ensure_ascii=False)
    log("\n=== RESUMEN (residuos flexibles A:%s) ===" % args.residuos)
    for r in resultados:
        log("  %-5s %-8s E(cristal) %-7s E(mejor) %-7s dE %-7s mejor %-7s mediana %-7s %s"
            % (r["entrada"], r["puntuacion"], r["E_cristal_quieta"],
               r["E_mejor_modo"], r["delta_E"], r["rmsd_mejor"], r["rmsd_mediana"],
               r["criterio"]))
    log("guardado: %s" % destino)
    return 0


if __name__ == "__main__":
    sys.exit(main())
