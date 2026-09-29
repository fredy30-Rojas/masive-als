#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Elige una diana de DUD-E para la que el control de redocking valga de verdad.

POR QUE
-------
El control de CDK2 no sirve: su ligando cristalino se modela en DOS conformaciones
al 50 %, el brazo del amonio no toca ni proteina ni agua, y ni Vina ni el CNN de GNINA
prefieren la geometria cristalografica. Un control que falla asi no dice nada del
embudo, asi que hace falta OTRO: una diana cuyo ligando sea rigido y este anclado de
verdad, y en la que la pose del cristal sea la que cualquier motor deberia encontrar.

LOS CRITERIOS, Y POR QUE SE PUEDEN MEDIR SIN ACOPLAR
-----------------------------------------------------
Un redocking se acierta cuando la pose del cristal es un minimo LOCAL OBVIO del
paisaje energetico, y eso se puede prever sin docking con dos medidas:

  1. **RIGIDEZ** — cuantos enlaces rotatorios tiene el ligando. Un ligando con
     brazos flexibles tiene muchas conformaciones casi equivalentes, y ninguna es
     "la" del cristal. Se quiere 0 o 1 enlace rotatorio.
  2. **ANCLAJE** — si cada polar del ligando (N y O) tiene una pareja en la
     proteina a menos de 3,5 A. Un polar sin pareja es un brazo suelto: lo que en
     CDK2 se va de 6 a 9 A. Se cuenta cuantos hay SIN pareja.
  3. **SIN CHOQUES** — la distancia minima heavy-heavy no puede ser de contacto
     covalentico (DUD-E quita el ligando del receptor, asi que deberia salir
     >= 2,5 A). Si sale menos, el receptor y el ligando no son del mismo cristal.
  4. **TAMAÑO** — que quepa en la caja de 24 A del proyecto sin que el motor tenga
     que recortarlo.
  5. **Y el cristal, no el mol2** — para las ganadoras se baja el PDB original y se
     mira cuantas copias del ligando hay y cuales son sus factores B. Eso es
     exactamente lo queTiempo de 京赛车女ㅁ desacREDITÓ en CDK2 (dos copias al 50 %,
     B ~47), asi que hay que comprobarlo antes de prometer nada.

QUE DESCARTA Y QUE NO
---------------------
Esto NO dice que el redocking vaya a pasar. Dice que la diana es la adecuada para
comprobarlo: descarta las que no pueden servir y ordena las que si. La prueba de
fuego es acoplar, y esa va con la GPU.

Uso:
    python elegir_control_dude.py            # fase 1 (102 ligandos) y fase 2
    python elegir_control_dude.py --fase 1   # solo el cribado de rigidez
Salida: control_dude_candidatos.csv, control_dude.txt
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import time
import urllib.request

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "_cache_dude")
CSV_SALIDA = os.path.join(BASE, "control_dude_candidatos.csv")
INFORME = os.path.join(BASE, "control_dude.txt")
BASE_URL = "https://dude.docking.org/targets"
PDB_URL = "https://files.rcsb.org/download/%s.pdb"

# Criterios de la criba.
MAX_ROTATORIOS = 1      # 0 o 1: si tiene brazos, el cristal no puede decidir
MIN_PESADOS = 12
MAX_PESADOS = 45
CORTE_PAREJA = 3.5      # un polar del ligando sin N/O de proteina a esta distancia
                        # esta suelto
CORTE_EXPUESTO = 4.5    # un atomo del ligando sin nada de proteina a esta distancia
                        # esta suelto en el espacio


def bajar(url, destino, pausa=0.3):
    """Descarga a cache. Si ya esta, no vuelve a pedirla."""
    if os.path.exists(destino) and os.path.getsize(destino) > 0:
        return open(destino, "rb").read()
    for intento in range(3):
        try:
            d = urllib.request.urlopen(url, timeout=60).read()
            with open(destino, "wb") as f:
                f.write(d)
            time.sleep(pausa)
            return d
        except Exception as e:  # noqa: BLE001
            if intento == 2:
                return b"__FALLO__%s" % str(e).encode()
            time.sleep(2)


def leer_mol2(datos, ruta_cache):
    """Molecula del ligando cristalino: sin H, con sanitizado si se puede."""
    with open(ruta_cache, "wb") as f:
        f.write(datos)
    m = Chem.MolFromMol2File(ruta_cache, removeHs=True, sanitize=True)
    if m is None:
        m = Chem.MolFromMol2File(ruta_cache, removeHs=True, sanitize=False)
        if m is None:
            return None
        m.UpdatePropertyCache(strict=False)
        try:
            Chem.SanitizeMol(m, Chem.SanitizeFlags.SANITIZE_ALL
                             ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                             ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
        except Exception:  # noqa: BLE001
            return None
    return m


def proteina_heavy(ruta):
    """Coords y (elemento, residuo) de los atomos pesados del receptor.

    Los `receptor.pdb` de DUD-E son de formato viejo (54 columnas, sin columna de
    elemento), asi que el elemento se saca del NOMBRE del atomo, como se hizo con
    CDK2. Un "H" al principio del nombre es un hidrogeno; "CA" es el alfa-carbono, no
    calcio, porque los iones serian HETATM y DUD-E no trae ninguno.
    """
    coords, elem, res = [], [], []
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if not l.startswith("ATOM"):
            continue
        nombre = l[12:16].strip()
        if not nombre or nombre[0] == "H":
            continue
        e = nombre[0].upper()
        if e not in ("C", "N", "O", "S", "P", "F", "I", "B"):
            continue
        coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        elem.append(e)
        res.append((l[17:20].strip().upper(), int(l[22:26])))
    return np.array(coords), elem, res


def medir_anclaje(mol, ruta_rec):
    """(sin_pareja, distancia_minima, fraccion_expuesta, residuos_tocados)."""
    P, elem, res = proteina_heavy(ruta_rec)
    if len(P) == 0:
        return None
    conf = mol.GetConformer()
    L = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                   conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    d = np.linalg.norm(L[:, None, :] - P[None, :, :], axis=2)
    dmin = float(d.min())
    # polar del ligando sin N/O de proteina cerca = brazo suelto
    idx_pol = np.array([a.GetIdx() for a in mol.GetAtoms()
                        if a.GetSymbol() in ("N", "O")], dtype=int)
    idx_prot_pol = np.array([i for i, e in enumerate(elem) if e in ("N", "O")],
                            dtype=int)
    if len(idx_pol) and len(idx_prot_pol):
        dpol = d[np.ix_(idx_pol, idx_prot_pol)]
        sin_pareja = int(np.count_nonzero(dpol.min(axis=1) > CORTE_PAREJA))
        n_polares = len(idx_pol)
    else:
        sin_pareja, n_polares = 0, 0
    expuesto = float((d.min(axis=1) > CORTE_EXPUESTO).mean())
    tocados = len({r for r, dm in zip(res, d.min(axis=1)) if dm < 4.0})
    return sin_pareja, n_polares, dmin, expuesto, tocados


def id_pdb_del_mol2(datos):
    """El codigo PDB esta en el nombre del mol2: `FAP_A1300_1H00` -> `1H00`."""
    m = re.search(rb"@<TRIPOS>MOLECULE\s*\n\s*(\S+)", datos)
    if not m:
        return None
    partes = m.group(1).decode("utf-8", "ignore").split("_")
    for p in reversed(partes):
        if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", p):
            return p.upper()
    return None


def fase1(lineas, log):
    """Criba de rigidez sobre las 102 dianas, sin bajar ningun receptor."""
    dianas = [l.strip() for l in
              open(os.path.join(BASE, "_dianas_dude.txt"), encoding="utf-8")
              if l.strip()]
    os.makedirs(CACHE, exist_ok=True)
    filas, fallan = [], []
    log("")
    log("FASE 1. RIGIDEZ Y TAMANO DEL LIGANDO CRISTALINO (%d dianas)" % len(dianas))
    for i, d in enumerate(dianas, 1):
        datos = bajar("%s/%s/crystal_ligand.mol2" % (BASE_URL, d),
                      os.path.join(CACHE, "%s_crystal_ligand.mol2" % d))
        if datos.startswith(b"__FALLO__"):
            fallan.append(d)
            continue
        mol = leer_mol2(datos, os.path.join(CACHE, "%s_mol2.mol2" % d))
        if mol is None:
            fallan.append(d)
            continue
        rot = int(rdMolDescriptors.CalcNumRotatableBonds(mol))
        n = mol.GetNumAtoms()
        n_pol = sum(1 for a in mol.GetAtoms() if a.GetSymbol() in ("N", "O"))
        anillos = rdMolDescriptors.CalcNumRings(mol)
        filas.append({"diana": d, "pesados": n, "rotatorios": rot,
                      "anillos": anillos, "polares": n_pol,
                      "pdb": id_pdb_del_mol2(datos) or ""})
        if i % 20 == 0:
            log("   %d/%d" % (i, len(dianas)))
    filas.sort(key=lambda r: (r["rotatorios"], r["pesados"]))
    log("   %d ligandos leidos, %d que no se pudieron abrir" % (len(filas),
                                                               len(fallan)))
    if fallan:
        log("   no se leyeron: %s" % ", ".join(fallan))
    log("")
    log("   los 20 mas rigidos (0 o 1 enlace rotatorio, tamano razonable):")
    log("   %-10s %-7s %-9s %-8s %-8s %s"
        % ("diana", "pesados", "rotatorios", "anillos", "polares", "pdb"))
    preseleccion = [r for r in filas
                    if r["rotatorios"] <= MAX_ROTATORIOS
                    and MIN_PESADOS <= r["pesados"] <= MAX_PESADOS]
    for r in preseleccion[:20]:
        log("   %-10s %-7d %-9d %-8d %-8d %s"
            % (r["diana"], r["pesados"], r["rotatorios"], r["anillos"],
               r["polares"], r["pdb"]))
    log("")
    log("   con 0 enlaces rotatorios en total: %d"
        % sum(1 for r in filas if r["rotatorios"] == 0))
    return filas, preseleccion


def fase2(lineas, log, filas, preseleccion, cuantos):
    """Baja el receptor de las mejor preseleccionadas y mide el anclaje."""
    os.makedirs(CACHE, exist_ok=True)
    anchos = preseleccion[:cuantos]
    log("")
    log("FASE 2. ANCLAJE (baja el receptor de las %d mejor preseleccionadas)"
        % len(anchos))
    log("   %-10s %-9s %-11s %-9s %-9s %s"
        % ("diana", "rotatorios", "sin pareja", "d min", "expuesto", "residuos"))
    out = []
    for r in anchos:
        datos = bajar("%s/%s/receptor.pdb" % (BASE_URL, r["diana"]),
                      os.path.join(CACHE, "%s_receptor.pdb" % r["diana"]))
        if datos.startswith(b"__FALLO__"):
            log("   %-10s no se pudo bajar el receptor" % r["diana"])
            continue
        mol = Chem.MolFromMol2File(os.path.join(CACHE, "%s_mol2.mol2" % r["diana"]),
                                   removeHs=True, sanitize=True)
        if mol is None:
            continue
        m = medir_anclaje(mol, os.path.join(CACHE, "%s_receptor.pdb" % r["diana"]))
        if m is None:
            continue
        sin_pareja, n_pol, dmin, expuesto, tocados = m
        r.update({"sin_pareja": sin_pareja, "polares_reales": n_pol,
                  "dmin": dmin, "expuesto": expuesto, "residuos": tocados})
        out.append(r)
        log("   %-10s %-9d %-11s %-9.2f %-9.2f %d"
            % (r["diana"], r["rotatorios"], "%d de %d" % (sin_pareja, n_pol),
               dmin, expuesto, tocados))
    return out


def fase3(lineas, log, anchos, cuantos):
    """El PDB original de las ganadoras: cuantas copias del ligando y que factores B."""
    log("")
    log("FASE 3. EL CRISTAL DE VERDAD (esto es lo que mato a CDK2)")
    out = []
    for r in anchos[:cuantos]:
        if not r.get("pdb"):
            continue
        datos = bajar(PDB_URL % r["pdb"],
                      os.path.join(CACHE, "%s.pdb" % r["pdb"]))
        if datos.startswith(b"__FALLO__"):
            log("   %-10s %s: no se pudo bajar el PDB original" % (r["diana"],
                                                                   r["pdb"]))
            continue
        texto = datos.decode("utf-8", "ignore")
        nombres = {}
        for l in texto.splitlines():
            if not l.startswith("HETATM"):
                continue
            het = l[17:20].strip()
            if het in ("HOH", "WAT"):
                continue
            if het in ("ACE", "NME", "DOD"):
                continue
            nombres.setdefault(het, []).append(l)
        bmed = []
        for het, lineas_h in nombres.items():
            bs = [float(l[60:66]) for l in lineas_h
                  if l[60:66].strip()]
            bmed.append((het, len(lineas_h), float(np.mean(bs)) if bs else 0.0))
        n_copias = len(nombres)
        peor_b = max((b for _, _, b in bmed), default=0.0)
        r.update({"hetero": len(nombres), "b_medio": peor_b})
        out.append(r)
        log("   %-10s %s | %d grupo(s) HETATM: %s"
            % (r["diana"], r["pdb"], n_copias,
               ", ".join("%s (%d at., B %.1f)" % x for x in bmed[:4])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fase", type=int, default=3,
                    help="hasta que fase llegar (1 = solo rigidez)")
    ap.add_argument("--corting", type=int, default=14,
                    help="cuantas dianas pasan a la fase 2")
    ap.add_argument("--cristal", type=int, default=6,
                    help="cuantas pasan a la fase 3")
    ap.add_argument("--rotatorios", type=int, default=None,
                    help="cuantos enlaces rotatorios se admiten (por defecto el "
                         "de MAX_ROTATORIOS; poner 2 abre la criba)")
    args = ap.parse_args()
    if args.rotatorios is not None:
        # Se cambia la variable del modulo, no la de la funcion: `fase1` la lee
        # desde el ambito del modulo.
        globals()["MAX_ROTATORIOS"] = args.rotatorios

    lineas = []

    def log(m):
        print(m, flush=True)
        lineas.append(m)

    log("ELEGIR UNA DIANA DE DUD-E PARA EL QUE EL CONTROL VALGA")
    log("   %s | criterios: rotatorios <= %d, %d-%d atomos, sin choques"
        % (time.strftime("%Y-%m-%d %H:%M"), MAX_ROTATORIOS, MIN_PESADOS,
           MAX_PESADOS))
    filas, preseleccion = fase1(lineas, log)
    anchos = []
    if args.fase >= 2:
        anchos = fase2(lineas, log, filas, preseleccion, args.corting)
    if args.fase >= 3 and anchos:
        anchos = fase3(lineas, log, anchos, args.cristal)

    with open(CSV_SALIDA, "w", newline="", encoding="utf-8") as f:
        campos = ["diana", "pesados", "rotatorios", "anillos", "polares", "pdb",
                  "sin_pareja", "polares_reales", "dmin", "expuesto", "residuos",
                  "hetero", "b_medio"]
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for r in filas:
            w.writerow({k: r.get(k, "") for k in campos})
    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    print("")
    print("csv: %s | informe: %s" % (CSV_SALIDA, INFORME))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
