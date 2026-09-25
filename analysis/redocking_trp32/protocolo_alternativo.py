#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Segundo control de redocking: otras funciones de puntuacion (20 sep 2026).

POR QUE EXISTE
--------------
`redock_trp32.py` habia medido el control con un medidor de RMSD que emparejaba
los atomos al reves e inflaba el resultado. Con el medidor corregido ese control
pasa en dos de los cuatro ligandos: isoproterenol 0,48-1,28 A y adrenalina
0,66-0,74 A; sigue fallando en dopamina (3,13-3,15 A) y en 5-fluorouridina
(11,1-11,3 A). La pregunta que este script responde es si ese resto de fallo es
**el programa y su funcion de puntuacion** o algo mas profundo (receptor
rigido, tipado de atomos, o el propio sistema). Ver
`INFORME_CONTROLES_CORREGIDOS_2026-09-20.md`.

Este script mide tres cosas, de la mas barata a la mas costosa, sobre las
MISMAS cuatro estructuras y las MISMAS cajas:

1. `--score_only`: puntuar la pose del cristal tal cual, sin moverla. Comparada
   con la afinidad del mejor modo acoplado da el numero decisivo:
       dE = E(cristal) - E(mejor modo)
   - dE > 0  -> la funcion prefiere una pose que NO es la del cristal: el fallo
               esta en el paisaje de puntuacion (la busqueda no puede ganar).
   - dE < 0  -> la pose del cristal puntua mejor pero no aparece: el fallo esta
               en la busqueda.
2. `--local_only`: minimizacion local partiendo del cristal. Mide si la pose
   cristalina es un minimo local real de esa funcion (deriva < 1 A) o si el
   paisaje la repele activamente (deriva > 2 A).
3. Acoplamiento completo con cada funcion (`vina`, `vinardo`, `ad4`), caja de
   24 A centrada en el centroide del cristal, exhaustividad 8, tres semillas,
   midiendo RMSD sin superponer con el mismo criterio que el primer control.

NOTA SOBRE `tools/smina.exe`: esta corrupto (9 bytes, contiene el texto
«Not Found»), asi que no se pudo usar. Vina 1.2.3 ya trae `vinardo` y `ad4`,
que es lo que se usa aqui en su lugar.

Uso:
    python protocolo_alternativo.py --puntuaciones vina,vinardo,ad4
    python protocolo_alternativo.py --solo-score     # pasos 1 y 2 nada mas
"""
import argparse
import json
import os
import re
import subprocess
import sys

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

import redock_trp32 as R

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
VINA = R.VINA
SISTEMAS = R.SISTEMAS
CAJA = 24
SEMILLAS = [42, 2026, 777]


def log(m):
    print(m, flush=True)


def pdbqt_del_cristal(cry_pdb, sdf, destino):
    """Convierte la pose del cristal a PDBQT conservando sus coordenadas.

    El PDB del cristal no trae ordenes de enlace, y RDKit los percibe mal (el
    anillo del catecol lo lee como ciclohexano saturado, el uracilo pierde los
    dobles enlaces). Eso arruinaria el tipado de atomos del PDBQT y con el la
    puntuacion. Se arregla copiando los ordenes de enlace desde el ligando
    ideal de RCSB (`AssignBondOrdersFromTemplate`) SIN tocar las coordenadas.
    """
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    ideal = Chem.MolFromMolFile(sdf)          # solo atomos pesados
    if ideal is None:
        return False
    cry = Chem.MolFromPDBBlock(open(cry_pdb).read(), removeHs=False,
                               proximityBonding=True)
    if cry is None:
        return False
    cry = Chem.RemoveHs(cry)
    if cry.GetNumAtoms() != ideal.GetNumAtoms():
        log("  AVISO: el cristal tiene %d atomos y el ideal %d; no se puede copiar "
            "el tipado" % (cry.GetNumAtoms(), ideal.GetNumAtoms()))
        return False
    try:
        m = AllChem.AssignBondOrdersFromTemplate(ideal, cry)
    except Exception as e:
        log("  AVISO: fallo copiando ordenes de enlace: %s" % e)
        return False
    m = Chem.AddHs(m, addCoords=True)
    prep = MoleculePreparation(rigid_macrocycles=True)
    try:
        setups = prep.prepare(m)
        pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    except Exception as e:
        log("  AVISO: fallo preparando el ligando del cristal: %s" % e)
        return False
    if not ok:
        return False
    open(destino, "w", encoding="utf-8").write(pdbqt)
    return True


def afinidad_stdout(txt):
    """Energia libre estimada que Vina imprime en --score_only y --local_only.

    Ojo con la convencion: es la MISMA que la afinidad que Vina reporta en un
    acoplamiento normal (incluye el termino torsional (3)), asi que se puede
    comparar directamente con la afinidad del mejor modo.
    """
    if not txt:
        return None
    vals = re.findall(r"Estimated Free Energy of Binding\s*:\s*([-+]?\d+\.\d+)", txt)
    return float(vals[-1]) if vals else None


def afinidad_de_salida(f):
    """Afinidad del primer modo de un PDBQT de salida de un acoplamiento."""
    if not os.path.exists(f):
        return None
    for l in open(f, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            return float(re.search(r"([-+]?\d+\.\d+)", l).group(1))
    return None


def hacer_rigido(entrada_pdbqt, salida):
    """Convierte un ligando flexible en uno RIGIDO en su conformacion actual.

    Poniendo todos los atomos dentro del ROOT y TORSDOF 0, el acoplamiento deja
    de muestrear torsiones y solo busca traslacion y rotacion. Sirve para
    separar las dos causas de un redocking fallido: si con el ligando rigido la
    pose del cristal SI se reproduce, el problema estaba en el muestreo de
    torsiones; si tampoco, el problema es la funcion de puntuacion.
    Las coordenadas de los atomos de los BRANCH ya estan en el sistema del
    receptor, asi que se pueden meter tal cual dentro del ROOT.
    """
    atomos = [l for l in open(entrada_pdbqt, encoding="utf-8", errors="ignore")
              if l.startswith(("ATOM  ", "HETATM"))]
    if not atomos:
        return False
    with open(salida, "w", encoding="utf-8") as f:
        f.write("ROOT\n")
        f.writelines(atomos)
        f.write("ENDROOT\nTORSDOF 0\n")
    return True


def coord_pesados(pdbqt, solo_primer_modelo=True):
    """Coordenadas de los atomos pesados de un PDBQT, en el orden del fichero.

    Las salidas de Vina traen los 9 modos concatenados dentro del mismo fichero
    (un MODEL por modo). Para comparar con el ligando de entrada hay que leer
    solo el primer modelo; si no, el recuento de atomos no cuadra.
    """
    out = []
    for l in open(pdbqt, encoding="utf-8", errors="ignore"):
        if solo_primer_modelo and out and l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM  ", "HETATM")):
            tipo = l[77:79].strip() if len(l) > 78 else ""
            if tipo.startswith("H"):
                continue
            out.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


def rmsd_posicional(a, b):
    """RMSD por orden de atomo (vale para el ligando rigido).

    Con el ligando rigido Vina conserva el orden de atomos del fichero de
    entrada, asi que la correspondencia es directa. Se usa solo cuando no hay
    ambiguedad de simetria, que es el caso de las tres catecolaminas y de la
    5-fluorouridina (sustituyentes distintos en cada posicion del anillo).
    """
    A, B = coord_pesados(a), coord_pesados(b)
    if not A or len(A) != len(B):
        return None
    s = 0.0
    for i in range(len(A)):
        for j in range(3):
            s += (A[i][j] - B[i][j]) ** 2
    return (s / len(A)) ** 0.5


def correr_vina(args, timeout=3600):
    r = subprocess.run([VINA] + args, capture_output=True, text=True,
                       timeout=timeout)
    return (r.stdout or "") + (r.stderr or "")


def evaluar(entrada, puntuacion, rec_pdbqt, caja, exhaustividad, semillas,
            solo_score=False, rigido=False):
    codigo, cadena, nombre = SISTEMAS[entrada]
    cry_pdb = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
    sdf = os.path.join(BASE, "ligands", "%s_ideal.sdf" % codigo)
    cry_pdbqt = os.path.join(BASE, "ligands", "%s_cristal_%s.pdbqt"
                             % (entrada, puntuacion))
    if not os.path.exists(cry_pdbqt):
        if not pdbqt_del_cristal(cry_pdb, sdf, cry_pdbqt):
            log("  FALLO: no se pudo preparar la pose del cristal")
            return None
    etiqueta = puntuacion + ("-rigido" if rigido else "")
    if rigido:
        rig = os.path.join(BASE, "ligands", "%s_cristal_rigido.pdbqt" % entrada)
        if not os.path.exists(rig):
            if not hacer_rigido(cry_pdbqt, rig):
                log("  FALLO: no se pudo hacer rigido el ligando")
                return None
        cry_pdbqt = rig

    cx, cy, cz = R.centroide(cry_pdb)
    base_args = ["--receptor", rec_pdbqt, "--ligand", cry_pdbqt,
                 "--scoring", puntuacion, "--autobox"]
    # las etiquetas de las trazas llevan «-rigido» para no confundir corridas

    # 1. puntuacion de la pose del cristal, sin moverla
    txt = correr_vina(base_args + ["--score_only"], timeout=1200)
    e_cristal = afinidad_stdout(txt)
    log("  [%s] E(cristal, quieta)      = %s kcal/mol"
        % (etiqueta, "%.2f" % e_cristal if e_cristal is not None else "no leida"))
    if e_cristal is None:
        log("     salida de score_only (ultimas lineas): %s"
            % str(txt.strip().splitlines()[-3:]))

    # 2. minimizacion local desde el cristal
    out_min = os.path.join(BASE, "out", "%s_min_%s.pdbqt" % (entrada, etiqueta))
    txt = correr_vina(base_args + ["--local_only", "--out", out_min],
                      timeout=1800)
    e_min = afinidad_stdout(txt)
    deriva = None
    if os.path.exists(out_min) and os.path.getsize(out_min) > 100:
        if rigido:
            deriva = rmsd_posicional(cry_pdbqt, out_min)
        else:
            try:
                deriva = R.rmsd_en_sitio(cry_pdb, R.molecula_dockeada(out_min))
            except Exception as e:
                log("  AVISO: no se pudo medir la deriva: %s" % e)
    log("  [%s] E(cristal minimizado) = %s kcal/mol | deriva %.2f A"
        % (etiqueta,
           "%.2f" % e_min if e_min is not None else "no leida",
           deriva if deriva is not None else float("nan")))

    if solo_score:
        return {"entrada": entrada, "ligando": codigo, "nombre": nombre,
                "puntuacion": etiqueta, "E_cristal_quieta": e_cristal,
                "E_cristal_minimizado": e_min,
                "deriva_minimizacion": round(deriva, 2) if deriva is not None else None}

    # 3. acoplamiento completo
    por_semilla = []
    for semilla in semillas:
        out = os.path.join(BASE, "out", "%s_%s_caja%d_e%d_s%d.pdbqt"
                           % (entrada, etiqueta, caja, exhaustividad, semilla))
        args = ["--receptor", rec_pdbqt, "--ligand", cry_pdbqt,
                "--center_x", "%.3f" % cx, "--center_y", "%.3f" % cy,
                "--center_z", "%.3f" % cz,
                "--size_x", str(caja), "--size_y", str(caja), "--size_z", str(caja),
                "--scoring", puntuacion, "--exhaustiveness", str(exhaustividad),
                "--num_modes", "9", "--seed", str(semilla), "--out", out]
        correr_vina(args)
        af = afinidad_de_salida(out)
        val = None
        if os.path.exists(out):
            if rigido:
                val = rmsd_posicional(cry_pdbqt, out)
            else:
                try:
                    val = R.rmsd_en_sitio(cry_pdb, R.molecula_dockeada(out))
                except Exception as e:
                    log("  AVISO: no se pudo medir el RMSD: %s" % e)
        log("  [%s] semilla %-5d afinidad %s | RMSD vs cristal %s"
            % (etiqueta, semilla,
               "%.2f" % af if af is not None else "no leida",
               "%.2f A" % val if val is not None else "no calculable"))
        por_semilla.append({"semilla": semilla, "afinidad": af,
                            "rmsd": round(val, 2) if val is not None else None})

    rmsds = [x["rmsd"] for x in por_semilla if x["rmsd"] is not None]
    mejor = min(rmsds) if rmsds else None
    mediana = sorted(rmsds)[len(rmsds) // 2] if rmsds else None
    mejor_af = min([x["afinidad"] for x in por_semilla
                    if x["afinidad"] is not None], default=None)
    de = (e_cristal - mejor_af) if (e_cristal is not None and mejor_af is not None) else None
    log("  [%s] -> mejor RMSD %s | mediana %s | pasa %d/%d | dE(cristal-mejor) %s"
        % (etiqueta,
           "%.2f A" % mejor if mejor is not None else "n/d",
           "%.2f A" % mediana if mediana is not None else "n/d",
           sum(1 for v in rmsds if v <= 2.0), len(rmsds),
           "%+.2f kcal/mol" % de if de is not None else "n/d"))
    return {
        "entrada": entrada, "ligando": codigo, "nombre": nombre,
        "puntuacion": etiqueta, "caja": caja,
        "E_cristal_quieta": e_cristal,
        "E_cristal_minimizado": e_min,
        "deriva_minimizacion": round(deriva, 2) if deriva is not None else None,
        "E_mejor_modo": mejor_af,
        "delta_E": round(de, 2) if de is not None else None,
        "por_semilla": por_semilla,
        "rmsd_mejor": mejor, "rmsd_mediana": mediana,
        "semillas_que_pasan": sum(1 for v in rmsds if v <= 2.0),
        "criterio": "PASA" if (mediana is not None and mediana <= 2.0) else "FALLA",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entradas", default=",".join(SISTEMAS))
    ap.add_argument("--puntuaciones", default="vina,vinardo,ad4")
    ap.add_argument("--exhaustividad", type=int, default=8)
    ap.add_argument("--caja", type=int, default=CAJA)
    ap.add_argument("--semillas", default=",".join(str(s) for s in SEMILLAS))
    ap.add_argument("--solo-score", action="store_true",
                    help="solo pasos 1 y 2, sin acoplamiento")
    ap.add_argument("--rigido", action="store_true",
                    help="acoplar el ligando RIGIDO en su conformacion cristalina: "
                         "separa el problema de muestreo del de puntuacion")
    args = ap.parse_args()

    semillas = [int(s) for s in args.semillas.split(",") if s.strip()]
    resultados = []
    for entrada in [e.strip() for e in args.entradas.split(",") if e.strip()]:
        rec_pdbqt = os.path.join(BASE, "receptores", entrada + ".pdbqt")
        if not os.path.exists(rec_pdbqt):
            log("=== %s: falta el receptor %s" % (entrada, rec_pdbqt))
            continue
        log("\n=== %s (%s) ===" % (entrada, SISTEMAS[entrada][2]))
        for p in [x.strip() for x in args.puntuaciones.split(",") if x.strip()]:
            r = evaluar(entrada, p, rec_pdbqt, args.caja, args.exhaustividad,
                        semillas, solo_score=args.solo_score, rigido=args.rigido)
            if r:
                resultados.append(r)

    sufijo = "score" if args.solo_score else ("rigido" if args.rigido else "completo")
    nombre_json = "resultado_protocolo_alternativo_%s_%s.json" % (
        sufijo, args.puntuaciones.replace(",", "_"))
    destino = os.path.join(BASE, nombre_json)
    with open(destino, "w", encoding="utf-8") as f:
        json.dump(resultados, f, indent=2, ensure_ascii=False)
    log("\n=== RESUMEN ===")
    for r in resultados:
        log("  %-5s %-8s E(cristal) %-7s E(mejor) %-7s dE %-7s deriva %-6s mejor %-7s mediana %-7s %s"
            % (r["entrada"], r["puntuacion"],
               r.get("E_cristal_quieta"), r.get("E_mejor_modo"), r.get("delta_E"),
               r.get("deriva_minimizacion"), r.get("rmsd_mejor"),
               r.get("rmsd_mediana"), r.get("criterio", "")))
    log("guardado: %s" % destino)
    return 0


if __name__ == "__main__":
    sys.exit(main())
