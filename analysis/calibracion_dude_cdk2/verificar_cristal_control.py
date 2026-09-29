#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pericia forense del cristal de las candidatas a diana de control.

POR QUE ESTE SCRIPT
-------------------
La criba por rigidez y anclaje dice que el ligando esta bien sujeto. Pero en CDK2
el ligando tambien estaba bien sujeto en el modelo, y aun asi el control no valia, por
tres motivos que solo se ven mirando el PDB ORIGINAL:

  1. el ligando venia DOS veces modelado (FAP y FCP, altloc A y B, ocupacion 0,50),
     o sea que el cristal no dice cual de las dos geometrias es la buena;
  2. los factores B de los 47 de media dicen que la densidad es debilisima;
  3. habia un brazo entero (el del amonio) flotando sin tocar nada.

Este script comprueba las tres cosas en las candidatas, ADEMAS de una cuarta que
no se comprueba en ningun sitio y es la mas peligrosa: que el mol2 de DUD-E sea de verdad
el mismo grupo que esta en el PDB. Si el receptor que bajo DUD-E y el mol2 vinieran
de cristales distintos, todas las medidas de anclaje serian mentira.

COMO SE COMPRUEBA LA IDENTIDAD
------------------------------
No se puede alinear por el orden de los atomos porque nadie lo garantiza. Se emparejan
con LAP por elemento: si el multiset de elementos del mol2 y el del HETATM coincide y
el emparejamiento optimo da 0,00 A, el mol2 ES ese grupo de ese cristal.

Uso:
    python verificar_cristal_control.py andr akt1 aces mapk2
"""
from __future__ import annotations

import os
import re
import sys
import time
import urllib.request

import numpy as np
from rdkit import Chem, RDLogger
from scipy.optimize import linear_sum_assignment

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(BASE, "_cache_dude")
PDB_URL = "https://files.rcsb.org/download/%s.pdb"


def bajar(url, destino, pausa=0.3):
    """Descarga a cache; si ya esta, no vuelve a pedirla."""
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
    return b"__FALLO__sin intentos"


def id_pdb_del_mol2(datos):
    """El codigo PDB esta en el nombre del mol2: `TES_A1300_1H00` -> `1H00`."""
    m = re.search(rb"@<TRIPOS>MOLECULE\s*\n\s*(\S+)", datos)
    if not m:
        return None
    for p in reversed(m.group(1).decode("utf-8", "ignore").split("_")):
        if re.fullmatch(r"[0-9][A-Za-z0-9]{3}", p):
            return p.upper()
    return None

# Los ligandos de DUD-E se llaman asi en el PDB. No se adivinan: se deduce del
# residuo que aparece en el nombre del mol2.
DESCARTE_HET = {"HOH", "WAT", "ACE", "NME", "DOD", "GOL", "EDO", "DMS", "PEG"}


def coords_mol2(ruta):
    """(coords Nx3, simbolos) del mol2 SOLO con atomos pesados.

    Los mol2 de DUD-E traen los hidrogenos explicitos, y `removeHs` no siempre los
    quita si la lectura no llega a sanitizar. Como el emparejamiento con el PDB se
    hace por numero de atomos, dejar los H dentro descuadra el emparejamiento y hace
    que se elija el grupo HETATM equivocado. Se filtran a mano.
    """
    m = Chem.MolFromMol2File(ruta, removeHs=True, sanitize=True)
    if m is None:
        m = Chem.MolFromMol2File(ruta, removeHs=True, sanitize=False)
    if m is None:
        return None
    c = m.GetConformer()
    L, e = [], []
    for a in m.GetAtoms():
        if a.GetSymbol() in ("H", "D"):
            continue
        p = c.GetAtomPosition(a.GetIdx())
        L.append([p.x, p.y, p.z])
        e.append(a.GetSymbol())
    return np.array(L), e


def coords_pdb_het(texto, resname):
    """(coords Nx3, elementos, altloc, ocupacion, B) del HETATM pedido.

    El elemento sale del nombre del atomo: los PDB antiguos no tienen columna de
    elemento. Un "H" al principio del nombre es hidrogeno, no helio.
    """
    L, e, alt, occ, b = [], [], [], [], []
    for l in texto.splitlines():
        if not l.startswith("HETATM"):
            continue
        if l[17:20].strip().upper() != resname:
            continue
        nombre = l[12:16].strip()
        if not nombre or nombre[0] == "H":
            continue
        L.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        e.append(elemento_de(nombre))
        alt.append(l[16].strip() or ".")
        try:
            occ.append(float(l[54:60]))
        except ValueError:
            occ.append(1.0)
        try:
            b.append(float(l[60:66]))
        except ValueError:
            b.append(0.0)
    if not L:
        return None
    return (np.array(L), e, alt, np.array(occ), np.array(b))


def elemento_de(nombre):
    """Elemento a partir del nombre del atomo, sin columna de elemento.

    Hay que mirar DOS letras antes que una: `CL1` es cloro y no carbono. Con una sola
    letra se leia como C y el emparejamiento con el mol2 fallaba.
    """
    dos = nombre[:2].upper()
    if dos in ("CL", "BR", "SE", "NA", "MG", "ZN", "FE", "MN", "CA", "CU", "NI",
               "CO", "CD", "HG", "PT", "AU", "AG", "LI", "SI", "AL", "TI", "CR",
               "PD", "RU", "IR", "SN", "SB", "BA", "CS", "RB", "SR", "TE", "AS"):
        return dos.capitalize() if dos in ("CL", "BR", "SE") else dos.capitalize()
    return nombre[0].upper()


def proteina(texto):
    """Coords y elementos de los ATOM (no HETATM): la parte proteica."""
    L, e, res = [], [], []
    for l in texto.splitlines():
        if not l.startswith("ATOM"):
            continue
        nombre = l[12:16].strip()
        if not nombre or nombre[0] == "H":
            continue
        e0 = elemento_de(nombre)
        if e0 not in ("C", "N", "O", "S", "P", "F", "I", "B", "CL", "BR"):
            continue
        L.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        e.append(e0)
        res.append((l[17:20].strip().upper(), int(l[22:26])))
    return np.array(L), e, res


def emparejar(A, ea, B, eb):
    """Emparejamiento optimo que exige el mismo elemento. Devuelve el RMSD."""
    if len(A) != len(B):
        return None
    D = np.linalg.norm(A[:, None, :] - B[None, :, :], axis=2)
    grande = 1e6
    for i, x in enumerate(ea):
        for j, y in enumerate(eb):
            if x != y:
                D[i, j] = grande
    r, c = linear_sum_assignment(D)
    if D[r, c].max() >= grande:
        return None
    d = D[r, c]
    return float(np.sqrt((d ** 2).mean()))


def analizar(diana, log):
    mol2 = os.path.join(CACHE, "%s_mol2.mol2" % diana)
    rec = os.path.join(CACHE, "%s_receptor.pdb" % diana)
    if not os.path.exists(mol2) or not os.path.exists(rec):
        log("   %-8s FALTA la cache de esta diana" % diana)
        return None
    cm = coords_mol2(mol2)
    if cm is None:
        log("   %-8s no se pudo abrir el mol2" % diana)
        return None
    Lm, em = cm

    # El nombre del residuo del ligando: el ultimo grupo de 3-4 letras del nombre
    # del mol2, o el nombre que coincida por numero de atomos con un HETATM.
    texto = open(rec, encoding="utf-8", errors="ignore").read()
    # OJO: el receptor.pdb de DUD-E viene SIN ningun HETATM (DUD-E le quita el
    # ligando). Los HETATM, y por tanto el ligando cristalino, solo estan en el PDB
    # original de RCSB. Por eso hacen falta los dos ficheros.
    pdb_id = id_pdb_del_mol2(open(mol2, "rb").read())
    if not pdb_id:
        log("   %-8s no se deduce el codigo PDB del mol2" % diana)
        return None
    pdb_txt = os.path.join(CACHE, "%s.pdb" % pdb_id)
    if not os.path.exists(pdb_txt) or os.path.getsize(pdb_txt) == 0:
        ptexto = bajar(PDB_URL % pdb_id, pdb_txt)
    else:
        ptexto = open(pdb_txt, "rb").read()
    if ptexto.startswith(b"__FALLO__"):
        log("   %-8s %s: no se pudo bajar el PDB original" % (diana, pdb_id))
        return None
    ptexto = ptexto.decode("utf-8", "ignore")

    # Grupo HETATM que tiene el MISMO numero de atomos que el mol2: esa es la
    # forma de identificar el ligando sin fiarnos de los nombres.
    cand = {}
    for l in ptexto.splitlines():
        if not l.startswith("HETATM"):
            continue
        het = l[17:20].strip().upper()
        if het in DESCARTE_HET:
            continue
        cand[het] = cand.get(het, 0) + 1
    con_pdb = []
    for het, n in sorted(cand.items(), key=lambda kv: abs(kv[1] - len(Lm))):
        p = coords_pdb_het(ptexto, het)
        if p is not None:
            con_pdb.append((het, p, abs(n - len(Lm))))

    elegido = None
    for het, (Lp, ep, alt, occ, b), dif in con_pdb:
        if dif == 0:
            elegido = (het, Lp, ep, alt, occ, b)
            break
    if elegido is None and con_pdb:
        het, (Lp, ep, alt, occ, b), dif = con_pdb[0]
        elegido = (het, Lp, ep, alt, occ, b)
    if elegido is None:
        log("   %-8s el receptor no tiene ningun HETATM con el tamano del ligando"
            % diana)
        return None
    het, Lp, ep, alt, occ, b = elegido

    log("")
    log("   %s | ligando %s | %d atomos en el mol2, %d en el PDB"
        % (diana, het, len(Lm), len(Lp)))
    rmsd = emparejar(Lp, ep, Lm, em)
    if rmsd is None:
        log("      IDENTIDAD: no se puede emparejar (elementos distintos)")
    else:
        log("      IDENTIDAD mol2 <-> PDB: RMSD %.4f A" % rmsd)
        if rmsd > 0.05:
            log("      AVISO: el mol2 y el PDB no son el mismo grupo")
        else:
            log("      -> el mol2 de DUD-E ES el HETATM de este receptor")

    # COPIAS, ALTLOCS Y OCUPACION
    alts = sorted(set(alt))
    n_res = set()
    for l in ptexto.splitlines():
        if l.startswith("HETATM") and l[17:20].strip().upper() == het:
            n_res.add((int(l[22:26]), l[16].strip() or "."))
    log("      copias en el cristal: %d (residuo+altloc: %s)"
        % (len(n_res), ", ".join("%d%s" % t for t in sorted(n_res))))
    log("      altloc presentes: %s" % (", ".join(alts) if len(alts) > 1
                                        else "ninguno"))
    log("      ocupacion: min %.2f max %.2f" % (occ.min(), occ.max()))
    log("      factor B: medio %.1f, min %.1f, max %.1f"
        % (b.mean(), b.min(), b.max()))

    # ANCLAJE CONTRA LA PROTEINA (recalculado aqui, sin depender de la criba)
    P, ep_, resp = proteina(texto)
    d = np.linalg.norm(Lp[:, None, :] - P[None, :, :], axis=2)
    dmin = d.min()
    dm_lig = d.min(axis=1)
    log("      d min ligando-proteina: %.2f A" % dmin)
    log("      atomos del ligando a >4.5 A de todo: %d de %d"
        % (int((dm_lig > 4.5).sum()), len(Lp)))
    vecinos = {}
    for i in range(len(Lp)):
        for j in np.where(d[i] < 4.0)[0]:
            vecinos.setdefault(resp[j], 0)
            vecinos[resp[j]] += 1
    top = sorted(vecinos.items(), key=lambda kv: -kv[1])[:6]
    log("      residuos a <4 A: %s"
        % ", ".join("%s%d (%d)" % (r, n, c) for (r, n), c in top))

    # CADA POLAR, SU PAREJA
    idx_pol = [i for i, x in enumerate(ep) if x in ("N", "O")]
    idx_pp = [j for j, x in enumerate(ep_) if x in ("N", "O")]
    sueltos = 0
    if idx_pol and idx_pp:
        for i in idx_pol:
            dp = d[i, idx_pp].min()
            if dp > 3.5:
                sueltos += 1
                log("      polar %d (%s) SIN pareja: el N/O mas cercano esta a %.2f A"
                    % (i + 1, ep[i], dp))
    log("      polares sin pareja: %d de %d" % (sueltos, len(idx_pol)))

    # AGUAS CERCANAS (la hipotesis que se refuto en CDK2)
    if ptexto:
        aguas = []
        for l in ptexto.splitlines():
            if l.startswith("HETATM") and l[17:20].strip().upper() in ("HOH", "WAT"):
                try:
                    aguas.append([float(l[30:38]), float(l[38:46]),
                                  float(l[46:54])])
                except ValueError:
                    pass
        if aguas:
            A = np.array(aguas)
            dw = np.linalg.norm(Lp[:, None, :] - A[None, :, :], axis=2).min()
            log("      agua mas cercana: %.2f A (%d aguas en el fichero)"
                % (dw, len(A)))
            idx_pol_w = [i for i in idx_pol]
            if idx_pol_w:
                dpol = np.linalg.norm(
                    Lp[idx_pol_w][:, None, :] - A[None, :, :], axis=2).min(axis=1)
                log("      polares con agua a <4.5 A: %d de %d"
                    % (int((dpol < 4.5).sum()), len(idx_pol_w)))

    return {"diana": diana, "het": het, "rmsd": rmsd, "copias": len(n_res),
            "altloc": len(alts) > 1, "b": float(b.mean()), "dmin": float(dmin),
            "sueltos": sueltos, "n_pol": len(idx_pol)}


def main():
    dianas = sys.argv[1:] or ["andr", "akt1", "aces", "mapk2"]
    lineas = []

    def log(m):
        print(m, flush=True)
        lineas.append(m)

    log("PERICIA FORENSE DE LAS CANDIDATAS A DIANA DE CONTROL")
    res = []
    for d in dianas:
        r = analizar(d, log)
        if r is not None:
            res.append(r)

    log("")
    log("RESUMEN")
    log("   %-8s %-8s %-9s %-7s %-7s %-7s %s"
        % ("diana", "ligando", "identidad", "copias", "altloc", "B medio",
           "sin pareja"))
    for r in res:
        log("   %-8s %-8s %-9s %-7d %-7s %-7.1f %d de %d"
            % (r["diana"], r["het"],
               "%.3f A" % r["rmsd"] if r["rmsd"] is not None else "n/d",
               r["copias"], "SI" if r["altloc"] else "no", r["b"],
               r["sueltos"], r["n_pol"]))
    with open(os.path.join(BASE, "control_dude_pericia.txt"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
