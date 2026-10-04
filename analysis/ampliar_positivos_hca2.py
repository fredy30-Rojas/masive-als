#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Amplia el banco de positivos de hCA2 de 60 a ~300, sin tocar lo ya acoplado.

POR QUE
-------
El 4 de octubre de 2026 el banco de hCA2 dejo su veredicto con la ventana de 9
modos: AUC crudo 0,446 con IC95 [0,360, 0,527]. El intervalo **cruza el azar**, y
el propio informe lo dice bien: no es prueba de que el embudo falle, es prueba de
que **este banco no alcanza**. Con 60 positivos no hay poder estadistico.

Y hay un hallazgo aparte que no se debe perder: el docking **si encuentra el
mecanismo** (coordinan el zinc 33 de 60 positivos frente a 265 de 2.331 senuelos,
diferencia con el intervalo entero por encima de cero). Lo que no ordena es la
kcal/mol. Eso es un resultado; lo que falta es un numero con intervalo que no
cruce el azar.

ChEMBL tiene **13.051 activos medidos** de hCA2 y solo se habian descargado 1.200.
Con mas activos hay mas familias de Murcko distintas, y por tanto mas positivos
independientes que se pueden tomar sin repetir familia.

QUE HACE, Y QUE NO
------------------
Amplia los POSITIVOS. No toca los 2.403 senuelos ya preparados ni el acople ya
hecho, y no reescribe `banco_hca2.txt` (es el registro de lo que se hizo el 30 de
septiembre y no se pisa).

  1. Descarga mas activos medidos de CHEMBL205 (`serie_hca2_ampliada.csv`),
     con la misma consulta que `auditar_serie_hca2.py`.
  2. Mantiene **los 60 positivos de siempre tal cual** (para que lo ya acoplado
     siga siendo comparable) y añade positivos nuevos, **uno por familia de
     Murcko nueva**, por pchembl descendente, hasta el objetivo.
  3. Prepara los nuevos con la receta canonica (`preparar_ligando.py`), al lado de
     los que ya estan, sin volver a preparar los 60.
  4. Audita: esqueletos, sulfonamidas y **cuantos de los 2.403 senuelos actuales
     comparten ahora esqueleto con un positivo nuevo** — esos inflarian el AUC por
     construccion y hay que quitarlos antes de medir.

LO QUE NO ARREGLA
-----------------
Ampliar positivos **no sube** el AUC: lo vuelve legible. Si hCA2 ordena de verdad
a 0,45, con 300 positivos el intervalo se estrecha alrededor de 0,45 y el veredicto
pasa de «no alcanza» a «no ordena, con confianza». Eso es justo lo falsable, y hay
que decirlo asi antes de ver el numero.

Uso:
    python ampliar_positivos_hca2.py --solo-listar
    python ampliar_positivos_hca2.py --objetivo 300 --descargar 6000
    python ampliar_positivos_hca2.py --cache      # sin red, con lo descargado
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from collections import defaultdict

from rdkit import Chem, RDLogger
from rdkit.Chem import rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
DIR = os.path.join(BASE, "calibracion_hca2")

sys.path.insert(0, BASE)
sys.path.insert(0, RAIZ)
# Se reutilizan las funciones del banco original a proposito: dos copias de una
# receta es el fallo mas caro de este proyecto (los pseudo-atomos de pegamento).
import preparar_banco_hca2 as PB  # noqa: E402
import preparar_ligando as PL     # noqa: E402

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
TARGET = "CHEMBL205"

SERIE_VIEJA = os.path.join(BASE, "serie_hca2_auditoria.csv")
SERIE_NUEVA = os.path.join(BASE, "serie_hca2_ampliada.csv")
ACTIVOS_VIEJOS = os.path.join(DIR, "activos_hca2.csv")
ACTIVOS_NUEVOS = os.path.join(DIR, "activos_hca2_ampliado.csv")
DECOYS = os.path.join(DIR, "decoys_hca2.csv")
LIGANDS = os.path.join(DIR, "ligands")
INFORME = os.path.join(BASE, "INFORME_AMPLIACION_HCA2.md")
LOG = os.path.join(BASE, "ampliar_positivos_hca2.log")

lineas = []


def log(m):
    print(m, flush=True)
    lineas.append(m)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")


def get(url, intentos=3):
    for i in range(intentos):
        try:
            req = urllib.request.Request(
                url, headers={"Accept": "application/json",
                              "User-Agent": "masive-als"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode("utf-8"))
        except Exception as e:  # noqa: BLE001
            if i == intentos - 1:
                log("      (fallo la consulta: %r)" % e)
                return None
            time.sleep(2 * (i + 1))
    return None


def descargar(n, ya_vistos):
    """Actividades con pchembl, deducplicadas por molecula, con subtipo.

    Misma consulta que `auditar_serie_hca2.py`: paginada, porque ChEMBL corta en
    1000 por pagina y pedir mas en una llamada devuelve 1000 sin avisar.
    """
    filas, vistos = [], set(ya_vistos)
    pagina = 0
    while len(filas) < n:
        q = urllib.parse.urlencode({
            "target_chembl_id": TARGET,
            "pchembl_value__isnull": "false",
            "standard_type__in": "IC50,Ki,Kd",
            "limit": 1000,
            "offset": pagina * 1000,
        })
        d = get("%s/activity.json?%s" % (CHEMBL, q))
        if not d:
            break
        acts = d.get("activities", [])
        if not acts:
            break
        for a in acts:
            smi = a.get("canonical_smiles")
            mid = a.get("molecule_chembl_id")
            if not smi or not mid or mid in vistos:
                continue
            try:
                pch = float(a.get("pchembl_value"))
            except (TypeError, ValueError):
                continue
            vistos.add(mid)
            filas.append({"molecule_chembl_id": mid, "smiles": smi,
                          "tipo": a.get("standard_type") or "",
                          "pchembl": "%.3f" % pch})
            if len(filas) >= n:
                break
        total = d.get("page_meta", {}).get("total_count", "?")
        log("   pagina %d: %d nuevos acumulados (de %s medidos en total)"
            % (pagina, len(filas), total))
        if len(acts) < 1000:
            break
        pagina += 1
        time.sleep(0.3)
    return filas


def leer_csv(ruta):
    if not os.path.exists(ruta):
        return []
    with open(ruta, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def esqueleto_de(smiles):
    m = Chem.MolFromSmiles(smiles or "")
    return PB.esqueleto(m) if m is not None else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--objetivo", type=int, default=300,
                    help="positivos totales a los que se quiere llegar")
    ap.add_argument("--descargar", type=int, default=6000,
                    help="activos medidos nuevos a bajar de ChEMBL")
    ap.add_argument("--cache", action="store_true",
                    help="no descargar; usar solo lo que ya hay en serie_hca2_ampliada.csv")
    ap.add_argument("--solo-listar", action="store_true",
                    help="elegir y auditar, pero no preparar PDBQT")
    args = ap.parse_args()

    log("AMPLIACION DE LOS POSITIVOS DE hCA2 — %s"
        % time.strftime("%Y-%m-%d %H:%M"))
    log("   objetivo: %d positivos" % args.objetivo)

    # ---- 1. la serie medida ----
    vieja = leer_csv(SERIE_VIEJA)
    log("")
    log("SERIE MEDIDA DE hCA2")
    log("   ya descargada: %d activos (%s)"
        % (len(vieja), os.path.basename(SERIE_VIEJA)))
    if os.path.exists(SERIE_NUEVA):
        nueva = leer_csv(SERIE_NUEVA)
        log("   ampliada en disco: %d activos" % len(nueva))
    else:
        nueva = []
    if not args.cache:
        vistos = {r["molecule_chembl_id"] for r in vieja + nueva}
        log("   bajando hasta %d activos nuevos de ChEMBL (CHEMBL205)..."
            % args.descargar)
        bajados = descargar(args.descargar, vistos)
        log("   bajados: %d" % len(bajados))
        if bajados:
            todas = vieja + nueva + bajados
            with open(SERIE_NUEVA, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=["molecule_chembl_id", "smiles",
                                                  "tipo", "pchembl"])
                w.writeheader()
                for r in todas:
                    w.writerow({k: r.get(k, "") for k in w.fieldnames})
            log("   serie escrita: %s (%d activos)"
                % (os.path.basename(SERIE_NUEVA), len(todas)))
    serie = leer_csv(SERIE_NUEVA) or vieja
    log("   serie de trabajo: %d activos medidos" % len(serie))

    # ---- 2. elegir los positivos: los 60 de siempre + familias nuevas ----
    viejos = leer_csv(ACTIVOS_VIEJOS)
    ids_viejos = [r["id"] for r in viejos]
    log("")
    log("POSITIVOS")
    log("   se conservan los %d de siempre (lo ya acoplado sigue siendo"
        " comparable)" % len(ids_viejos))

    mols, fam = {}, defaultdict(list)
    for r in serie:
        mid, smi = r.get("molecule_chembl_id"), r.get("smiles")
        if not mid or not smi:
            continue
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        try:
            pch = float(r.get("pchembl") or "nan")
        except ValueError:
            continue
        mols[mid] = (m, r.get("tipo") or "", pch)
        e = PB.esqueleto(m)
        if e:
            fam[e].append(mid)
    log("   activos con mol legible: %d" % len(mols))
    log("   familias de Murcko de 5+ activos: %d"
        % sum(1 for v in fam.values() if len(v) >= 5))
    log("   familias distintas en total: %d" % len(fam))

    fams_viejas = set()
    for r in viejos:
        e = esqueleto_de(r.get("smiles"))
        if e:
            fams_viejas.add(e)
    log("   familias ya cubiertas por los %d positivos: %d"
        % (len(ids_viejos), len(fams_viejas)))

    faltan = max(0, args.objetivo - len(ids_viejos))
    orden = sorted(mols, key=lambda c: -mols[c][2])
    fams, usados = set(fams_viejas), set(ids_viejos)

    def tirar(n):
        """Hasta n candidatos nuevos, uno por familia de Murcko nueva."""
        out = []
        for cid in orden:
            if len(out) >= n:
                break
            if cid in usados:
                continue
            e = PB.esqueleto(mols[cid][0])
            if e and e in fams:
                continue                      # familia ya representada
            usados.add(cid)
            if e:
                fams.add(e)
            out.append(cid)
        return out

    def preparar(lista):
        """(buenos, fallos). Idempotente: si el PDBQT ya esta, no lo rehace."""
        buenos, malos, sales = [], [], 0
        os.makedirs(LIGANDS, exist_ok=True)
        t0 = time.time()
        for k, cid in enumerate(lista, 1):
            nombre = "ACT_%s" % cid
            if os.path.exists(os.path.join(LIGANDS, nombre + ".pdbqt")):
                buenos.append(cid)
                continue
            smi = Chem.MolToSmiles(mols[cid][0])
            ruta = PL.escribir(nombre, smi, LIGANDS)
            if ruta is None or not os.path.exists(ruta):
                ruta = PL.escribir(nombre, smi, LIGANDS, quitar_sales=True)
                if ruta is None or not os.path.exists(ruta):
                    malos.append(cid)
                    continue
                sales += 1
            buenos.append(cid)
            if k % 25 == 0:
                log("   %d/%d (%.0f s)" % (k, len(lista), time.time() - t0))
        log("   preparados %d de %d en %.0f s"
            % (len(buenos), len(lista), time.time() - t0))
        if sales:
            log("   de los cuales %d salieron tras quitar la sal" % sales)
        return buenos, malos

    anadidos = []
    if args.solo_listar:
        anadidos = tirar(faltan)
    else:
        log("")
        log("PREPARANDO LOS NUEVOS (receta canonica del proyecto)")
        pendientes = tirar(faltan)
        while pendientes and len(anadidos) < faltan:
            buenos, malos = preparar(pendientes)
            anadidos += buenos
            if malos:
                log("   DESCARTADOS %d (siguen dando mas de un fragmento): %s"
                    % (len(malos), ", ".join(malos[:10])))
                log("   Un positivo que se cae sin avisar baja el AUC sin que se"
                    " note, asi que se repone con el siguiente de la lista.")
            pendientes = (tirar(faltan - len(anadidos))
                          if len(anadidos) < faltan else [])

    log("   anadidos: %d (objetivo total %d)"
        % (len(anadidos), len(ids_viejos) + len(anadidos)))
    if not anadidos:
        log("   nada que anadir; la serie descargada no aporta familias nuevas.")
    elif len(ids_viejos) + len(anadidos) < args.objetivo:
        log("   AVISO: se pidieron %d y solo se llega a %d. La serie descargada"
            % (args.objetivo, len(ids_viejos) + len(anadidos)))
        log("          no da mas familias de Murcko distintas. Descargar mas"
            " activos, no relajar la regla.")
    if anadidos:
        ps = [mols[c][2] for c in anadidos]
        log("   pchembl de los anadidos: max %.2f | min %.2f"
            % (max(ps), min(ps)))

    # ---- la sulfonamida: lo que se puede medir y lo que no ----
    if anadidos and not args.solo_listar:
        con_sulf, con_hd = 0, 0
        for cid in anadidos:
            ruta = os.path.join(LIGANDS, "ACT_%s.pdbqt" % cid)
            if not os.path.exists(ruta):
                continue
            if not mols[cid][0].HasSubstructMatch(PB.SULFONAMIDA):
                continue
            con_sulf += 1
            if PB.re_hd(open(ruta, encoding="utf-8", errors="ignore").read()):
                con_hd += 1
        log("")
        log("SULFONAMIDA (esto es lo que se puede decir, y lo que no)")
        log("   sulfonamidas entre los nuevos: %d de %d"
            % (con_sulf, len(anadidos)))
        log("   ficheros con algun hidrogeno polar (HD/HS): %d de %d"
            % (con_hd, con_sulf))
        log("   OJO CON ESE NUMERO, y por eso no se lee como un aviso:")
        log("     `re_hd` mira si el PDBQT trae ALGUN hidrogeno polar, no si el")
        log("     hidrogeno esta EN el nitrogeno del sulfonamida. Con cualquier")
        log("     otro donador en la molecula ya da que si, asi que no decide")
        log("     nada. El banco de 60 lo leyo como '60 de 60 mal' y ese aviso")
        log("     era un artefacto de la comprobacion, no una medida.")
        log("     Mapear el N del sulfonamida por indice tampoco vale: meeko")
        log("     reordena los atomos al escribir el PDBQT (fue el hallazgo del")
        log("     control de CDK2 el 29 de septiembre).")
        log("   Queda pendiente una comprobacion de verdad, con correspondencia de")
        log("   atomos deducida como en `diagnostico_orden_poses.py`. No se")
        log("   inventa aqui.")

    # ---- 4. escribir el CSV de positivos y auditar ----
    filas = list(viejos)
    for cid in anadidos:
        m = mols[cid][0]
        filas.append({"id": cid, "smiles": Chem.MolToSmiles(m),
                      "tipo": mols[cid][1], "pchembl": "%.3f" % mols[cid][2],
                      "sulfonamida": "SI" if m.HasSubstructMatch(PB.SULFONAMIDA)
                                     else "no",
                      "pesados": m.GetNumHeavyAtoms(),
                      "rotatorios": rdMolDescriptors.CalcNumRotatableBonds(m)})
    campos = ["id", "smiles", "tipo", "pchembl", "sulfonamida", "pesados",
              "rotatorios"]
    with open(ACTIVOS_NUEVOS, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for r in filas:
            w.writerow({k: r.get(k, "") for k in campos})
    log("")
    log("positivos escritos: %s (%d)" % (os.path.basename(ACTIVOS_NUEVOS),
                                         len(filas)))

    fams_pos = set()
    for r in filas:
        e = esqueleto_de(r.get("smiles"))
        if e:
            fams_pos.add(e)
    decoys = leer_csv(DECOYS)
    chocan, limpios = [], []
    for r in decoys:
        e = esqueleto_de(r.get("smiles"))
        (chocan if (e and e in fams_pos) else limpios).append(r)
    # El banco original NO se toca: se escribe al lado el que se usaria al medir.
    # Un senuelo con el mismo esqueleto que un positivo es casi indistinguible de
    # el y haria el AUC mas alto por construccion, no por el motor; es justo lo
    # que la auditoria de CDK2 enseño a mirar.
    ruta_limpios = os.path.join(DIR, "decoys_hca2_ampliado.csv")
    with open(ruta_limpios, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "smiles", "activo"])
        w.writeheader()
        for r in limpios:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})
    ruta_chocan = os.path.join(DIR, "decoys_hca2_excluidos.csv")
    with open(ruta_chocan, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "smiles", "activo"])
        w.writeheader()
        for r in chocan:
            w.writerow({k: r.get(k, "") for k in w.fieldnames})
    log("")
    log("AUDITORIA DEL BANCO AMPLIADO")
    log("   esqueletos de los %d positivos: %d" % (len(filas), len(fams_pos)))
    log("   senuelos actuales: %d" % len(decoys))
    log("   senuelos que comparten esqueleto con ALGUN positivo: %d"
        % len(chocan))
    log("   senuelos limpios (los que se usarian al medir): %d -> %s"
        % (len(limpios), os.path.basename(ruta_limpios)))
    log("   excluidos: %d -> %s" % (len(chocan), os.path.basename(ruta_chocan)))
    if chocan:
        log("   ATENCION: esos %d inflarian el AUC por construccion (son casi"
            % len(chocan))
        log("   indistinguibles de un positivo). Se dejan fuera al medir; el"
            " banco original NO se toca.")

    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("INFORME DE AMPLIACION DE LOS POSITIVOS DE hCA2 — %s\n\n"
                % time.strftime("%Y-%m-%d"))
        f.write("\n".join(lineas) + "\n")
    log("")
    log("informe: %s" % os.path.relpath(INFORME, RAIZ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
