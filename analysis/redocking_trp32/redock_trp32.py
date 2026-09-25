#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Control de redocking del bolsillo Trp32 de SOD1 (20 sep 2026).

POR QUE EXISTE ESTE SCRIPT
--------------------------
Hasta hoy el proyecto daba por bueno que el bolsillo Trp32 (caja centrada en
(46.5, 80.0, 73.3) sobre el receptor 1HL5) es el sitio de union de los ligandos
de SOD1, y ordenaba la libreria con esa caja. Pero los cuatro ligandos que se
usan como «patron de medida» (5-fluorouridina, isoproterenol, dopamina y
adrenalina) no salen de 1HL5: salen de las estructuras cristalograficas de
Wright 2013, que son del mutante I113T (4A7S, 4A7T, 4A7V, 4A7U).

Antes de usar esas poses como referencia hay que responder una pregunta basica
que el proyecto nunca se hizo: **el protocolo de acoplamiento reproduce la pose
cristalografica?** Eso es un control de redocking, y es el estandar minimo de
cualquier trabajo de docking (RMSD <= 2,0 A del modo mejor puntuado respecto al
cristal). Si falla, ni la caja ni los parametros valen, y todo ranking que
salga de ahi es papel mojado.

QUE HACE
--------
Para cada estructura: extrae el ligando nativo, prepara el receptor (meeko),
acopla el ligando nativo en su propia caja (Vina, exhaustividad 8, los mismos
parametros que usa el proyecto) y mide el RMSD con simetria corregida
(GetBestRMS) entre el mejor modo y la pose del cristal.

Uso:
    python redock_trp32.py                 # las cuatro estructuras I113T
    python redock_trp32.py --entradas 4A7V # solo una
"""
import argparse
import json
import os
import re
import subprocess
import sys

import numpy as np

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
PYDIR = r"C:\Users\Fredy\AppData\Local\Python\pythoncore-3.14-64"
MK_RECEPTOR = os.path.join(PYDIR, "Scripts", "mk_prepare_receptor.exe")
PYTHON = os.path.join(PYDIR, "python.exe")
VINA = r"C:\Users\Fredy\masive-als\tools\vina.exe"
PDB = os.path.join(BASE, "pdb")

# entrada -> (codigo del ligando, cadena cuya proteina se usa, nombre)
#
# OJO CON LAS COPIAS: cada entrada trae VARIAS copias del ligando en sitios
# distintos de la superficie (no todas en Trp32). En 4A7U, por ejemplo, la
# adrenalina aparece dos veces: A1000 en un sitio secundario a 15,5 A y A1001 en
# el bolsillo de Trp32. Elegir «la primera copia» a dedo fue un error de la
# primera version de este script y llevo a concluir en falso que la adrenalina
# no estaba en Trp32. La copia se elige ahora por criterio explicito: la que mas
# residuos del bolsillo de Trp32 toca.
#
# SEGUNDA RONDA (21 sep 2026). Cuando el barrido del PDB (ligandos_cristal_trp32.py)
# encontro siete quimias NUEVAS en Trp32, el control habia que repetirlo con ellas:
# son justo las que sostienen el unico resultado positivo del proyecto y ninguna
# se habia comparado aun con su pose cristalina. Esas estructuras NO son el mutante
# I113T de Wright (8GSQ, 6A9O y 5YTO son otras construcciones y otra numeracion),
# asi que para ellas la cadena se deja en None y la copia se elige sola por
# cercania al anillo de Trp32 (ver elegir_copia_trp32).
SISTEMAS = {
    # --- ronda 1: catecolaminas y nucleosido (Wright 2013, I113T)
    "4A7S": ("5UD", "A", "5-fluorouridina"),
    "4A7T": ("5FW", "A", "isoproterenol"),
    "4A7U": ("ALE", "A", "adrenalina"),
    "4A7V": ("LDP", "A", "dopamina"),
    # --- ronda 2: quimias nuevas del barrido del PDB (20 sep 2026, validacion v4)
    "4A7Q": ("4MQ", None, "diazepanoquinazolina"),
    "4A7G": ("12I", None, "piperazinaquinazolina"),
    "2WZ6": ("ZO0", None, "cf3-diazepanoquinazolina"),
    "2WZ0": ("ZZT", None, "anilina (2-metoxi-5-metilanilina)"),
    "8GSQ": ("K4I", None, "paliperidona"),
    "6A9O": ("6B3", None, "Lig9 (fenantridinona)"),
    "5YTO": ("946", None, "aminoalcohol naftalenico"),
}

# residuos que definen el bolsillo de Trp32 en la cadena A (Wright 2013)
RESIDUOS_TRP32 = {"21", "22", "23", "28", "29", "30", "32", "100"}

# anillo indol de la Trp32: los atomos que definen el bolsillo
ANILLO_TRP32 = ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2")


def _buscar_pdb(entrada):
    """Localiza el PDB de una entrada en las carpetas que usa el proyecto.

    La ronda 1 vive en redocking_trp32/pdb; las estructuras del barrido estan en
    la cache de trp32_cristal/pdb (descargadas por ligandos_cristal_trp32.py).
    """
    for carpeta in (os.path.join(BASE, "pdb"),
                    os.path.join(os.path.dirname(BASE), "trp32_cristal", "pdb")):
        for nombre in (entrada + ".pdb", entrada.lower() + ".pdb"):
            p = os.path.join(carpeta, nombre)
            if os.path.exists(p):
                return p
    return None


def elegir_copia_trp32(pdb, codigo):
    """Elige la copia del ligando mas pegada al anillo de Trp32. Devuelve (cadena, resSeq).

    Criterio geometrico puro, sin listas de residuos: distancia minima entre
    cualquier atomo del ligando y cualquier atomo del ANILLO INDOL de la Trp32 de
    SU MISMA CADENA. Sirve para cualquier estructura de SOD1, tenga la numeracion
    que tenga (la ronda 2 incluye construcciones que no son el mutante I113T).
    """
    cadenas = {}
    ligandos = {}
    for l in open(pdb, encoding="utf-8", errors="ignore"):
        if l.startswith("ATOM") and l[17:20].strip() == "TRP" and l[22:26].strip() == "32" \
                and l[12:16].strip() in ANILLO_TRP32:
            cadenas.setdefault(l[21], []).append(
                (float(l[30:38]), float(l[38:46]), float(l[46:54])))
        elif l.startswith("HETATM") and l[17:20].strip() == codigo:
            ligandos.setdefault((l[21], l[22:27].strip()), []).append(
                (float(l[30:38]), float(l[38:46]), float(l[46:54])))
    mejor, mejor_d = None, None
    for clave, xyz in sorted(ligandos.items()):
        anillo = cadenas.get(clave[0])
        if not anillo:
            continue
        d = min((x[0] - a[0]) ** 2 + (x[1] - a[1]) ** 2 + (x[2] - a[2]) ** 2
                for x in xyz for a in anillo) ** 0.5
        log("     copia %s%s: %.2f A al anillo de Trp32%s"
            % (clave[0], clave[1], d, "   <-- elegida" if mejor_d is None or d < mejor_d else ""))
        if mejor_d is None or d < mejor_d:
            mejor, mejor_d = clave, d
    return mejor if mejor else (None, None)

TAMANO_CAJA = 24  # igual que la caja de Trp32 del proyecto (22-24)


def log(m):
    print(m, flush=True)


def extraer_ligando(pdb, codigo, cadena, resseq, destino):
    """Escribe un PDB con los atomos de UNA copia del ligando nativo."""
    lineas = []
    for l in open(pdb, encoding="utf-8", errors="ignore"):
        if l.startswith("HETATM") and l[17:20].strip() == codigo \
                and l[21] == cadena and l[22:27].strip() == resseq:
            # marcar como ATOM para que RDKit lo lea como molecula
            lineas.append("ATOM  " + l[6:])
    if not lineas:
        return None
    lineas.append("END\n")
    with open(destino, "w", encoding="utf-8") as f:
        f.writelines(lineas)
    return destino


# Residuos que hay que descartar antes de preparar el receptor porque meeko no
# sabe construir la cadena peptidica con ellos.
#   2WZ6: el par Glu132-Glu133 de las dos cadenas tiene las dos conformaciones a
#   media ocupacion (0,50 y 0,30) y meeko detecta DOS enlaces entre los dos
#   residuos: «Expected 4 paddings for (A:132, A:133) with bonds [(0,4),(13,14)],
#   but got 2». Estan a 40 A de la caja (centro 15,5 / -1,6 / -10,1) y ni se ven
#   desde el bolsillo, asi que se quitan: sin ellos el receptor se prepara.
RESIDUOS_DESCARTADOS = {
    "2WZ6": ["A:132", "A:133", "F:132", "F:133"],
}


def extraer_receptor(pdb, destino, cadenas=None, descartar=()):
    """Proteina a solo, con UNA sola conformacion por atomo.

    Quita aguas, iones y hetero-atomos. Ademas resuelve las conformaciones
    alternas, que es donde meeko se atasca con estas entradas: en 2WZ6 se cae con
    «Requested altlocs not found for: F:102, F:132» (residuos que existen con
    altloc B y no con A) y en 8GSQ con «Template matching failed». Aqui se decide
    atomo por atomo: si el atomo tiene un registro SIN altloc, ese manda; si solo
    existe con altloc, se toma el de MAYOR OCUPACION (y 'A' como desempate). El
    altloc se borra del fichero, de modo que meeko ya no tiene nada que elegir.
    Es mas fiel que `--default_altloc A`, que reconstruye cadenas laterales del
    bolsillo en vez de conservar las coordenadas depositadas (ver
    INFORME_CONTROLES_CORREGIDOS_2026-09-20.md, seccion 4).
    """
    descartar = set(descartar)
    lineas = []
    for l in open(pdb, encoding="utf-8", errors="ignore"):
        if l.startswith("ATOM") or (l.startswith("HETATM") and l[17:20].strip() in
                                    ("MSE", "CYX")):
            if cadenas and l[21] not in cadenas:
                continue
            if "%s:%s" % (l[21], l[22:27].strip()) in descartar:
                continue
            lineas.append(l)

    elegida = {}
    for l in lineas:
        clave = (l[21], l[22:27], l[17:20], l[12:16])
        try:
            occ = float(l[54:60])
        except ValueError:
            occ = 0.0
        # menor es mejor: sin altloc primero, luego mayor ocupacion
        valor = (0 if l[16] == " " else 1, -occ, l[16])
        if clave not in elegida or valor < elegida[clave][0]:
            elegida[clave] = (valor, l)

    n = 0
    with open(destino, "w", encoding="utf-8") as f:
        for l in lineas:
            clave = (l[21], l[22:27], l[17:20], l[12:16])
            if elegida[clave][1] is not l:
                continue
            f.write(l[:16] + " " + l[17:])   # altloc fuera
            n += 1
        f.write("END\n")
    return n


def atomos_receptor(entrada, pdb_path=None):
    """Atomos de protein del PDB de una entrada: (cadena, residuo, nombre, xyz)."""
    out = []
    for l in open(pdb_path or os.path.join(PDB, entrada + ".pdb"), encoding="utf-8",
                  errors="ignore"):
        if l.startswith("ATOM"):
            out.append((l[21], l[22:27].strip(), l[17:20].strip(),
                        np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


def elegir_copia(pdb, codigo, cadena):
    """Devuelve el resSeq de la copia del ligando que ocupa el bolsillo Trp32.

    Criterio: la copia de esa cadena que toca mas residuos del conjunto
    {Glu21, Gln22, Lys23, Pro28, Val29, Lys30, Trp32, Glu100} a 4,5 A. Se
    reporta cual se eligio y cuantas copias habia, para que se vea a ojo.
    """
    atomos_rec = {}
    copias = {}
    for l in open(pdb, encoding="utf-8", errors="ignore"):
        if l.startswith("ATOM") and l[21] == cadena:
            atomos_rec.setdefault(l[22:27].strip(), []).append(
                (float(l[30:38]), float(l[38:46]), float(l[46:54])))
        if l.startswith("HETATM") and l[17:20].strip() == codigo and l[21] == cadena:
            copias.setdefault(l[22:27].strip(), []).append(
                (float(l[30:38]), float(l[38:46]), float(l[46:54])))
    mejor, mejor_n, detalle = None, -1, []
    for rs, xyz in copias.items():
        tocados = set()
        for r2 in RESIDUOS_TRP32:
            if r2 not in atomos_rec:
                continue
            if any(min((x[0] - a[0]) ** 2 + (x[1] - a[1]) ** 2 + (x[2] - a[2]) ** 2
                       for a in atomos_rec[r2]) <= 4.5 ** 2 for x in xyz):
                tocados.add(r2)
        detalle.append((rs, len(tocados), sorted(tocados, key=int)))
        if len(tocados) > mejor_n:
            mejor, mejor_n = rs, len(tocados)
    log("  copias de %s en la cadena %s: %d" % (codigo, cadena, len(copias)))
    for rs, n, toc in detalle:
        log("     resSeq %-5s toca %d residuos de Trp32 %s%s"
            % (rs, n, toc, "   <-- elegida" if rs == mejor else ""))
    return mejor


def centroide(pdb_lig):
    xs, ys, zs = [], [], []
    for l in open(pdb_lig):
        if l.startswith(("ATOM", "HETATM")):
            xs.append(float(l[30:38]))
            ys.append(float(l[38:46]))
            zs.append(float(l[46:54]))
    return sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs)


def aplanar(mol):
    """Deja la molecula en conectividad pura (todo enlace simple, sin aromas).

    Es necesario porque el PDB del cristal no trae ordenes de enlace y RDKit los
    percibe mal (el anillo de la dopamina lo lee como ciclohexano saturado). Para
    EMPAREJAR atomos basta la conectividad; el RMSD de redocking se mide sin
    superponer (las dos poses ya estan en el sistema del receptor), asi que los
    ordenes de enlace no entran en la cuenta.
    """
    m = Chem.RWMol(mol)
    for a in m.GetAtoms():
        a.SetIsAromatic(False)
        a.SetNoImplicit(True)
    for b in m.GetBonds():
        b.SetBondType(Chem.BondType.SINGLE)
        b.SetIsAromatic(False)
    m = m.GetMol()
    m.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(m)
    return m


def molecula_dockeada(out_pdbqt):
    from meeko import PDBQTMolecule, RDKitMolCreate
    pmol = PDBQTMolecule.from_file(out_pdbqt, skip_typing=True)
    mols = RDKitMolCreate.from_pdbqt_mol(pmol)
    return mols[0] if mols else None


def rmsd_en_sitio(pdb_lig, pose, conf_id=0):
    """RMSD de una pose frente al cristal, sin superponer poses.

    Se prueban todos los emparejamientos posibles (la simetria del anillo de
    catecol y de los metilos da varios) y se toma el minimo, que es el criterio
    estandar. conf_id permite pedir un modo concreto cuando el fichero trae
    varios.

    OJO — ERROR CORREGIDO EL 20 DE SEPTIEMBRE DE 2026. Hasta hoy esta funcion
    emparejaba los atomos al reves: hacia `cp[i]` contra `cc[match[i]]`. RDKit
    devuelve, en `match[i]`, el indice del atomo del OBJETIVO (la pose) que
    corresponde al atomo i de la CONSULTA (el cristal), asi que el
    emparejamiento correcto es `cp[match[i]]` contra `cc[i]`. Las dos formas
    solo coinciden si el emparejamiento es involutivo, cosa que casi nunca pasa
    cuando el orden de atomos del PDBQT (raiz primero) no es el del PDB del
    cristal. El efecto era INFLAR el RMSD: el isoproterenol de 4A7T pasaba de
    1,88 A a 0,46 A y la adrenalina de 4A7U de 3,74 A a 0,66 A. Los numeros de
    los controles anteriores a esa fecha estan recalculados en
    `recalcular_rmsd.py`; la comprobacion independiente es `rdMolAlign.CalcRMS`.
    """
    cry = Chem.MolFromPDBBlock(open(pdb_lig).read(), removeHs=False,
                               proximityBonding=True)
    if cry is None or pose is None:
        return None
    p = aplanar(Chem.RemoveHs(Chem.Mol(pose)))
    c = aplanar(Chem.RemoveHs(cry))
    if p.GetNumAtoms() != c.GetNumAtoms():
        return None
    cp = p.GetConformer(conf_id).GetPositions()
    cc = c.GetConformer().GetPositions()
    mejor = None
    for match in p.GetSubstructMatches(c, uniquify=False, maxMatches=500,
                                       useChirality=False):
        d2 = 0.0
        for i, j in enumerate(match):
            # i: atomo del cristal (consulta); match[i]: atomo homólogo de la pose
            dx = cp[j][0] - cc[i][0]
            dy = cp[j][1] - cc[i][1]
            dz = cp[j][2] - cc[i][2]
            d2 += dx * dx + dy * dy + dz * dz
        val = (d2 / len(match)) ** 0.5
        if mejor is None or val < mejor:
            mejor = val
    return mejor


def preparar_ligando(sdf, destino_pdbqt):
    """Convierte el ligando ideal (SDF de RCSB) a PDBQT, como el proyecto."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    m = Chem.MolFromMolFile(sdf, removeHs=False)
    if m is None:
        return False
    m = Chem.AddHs(m, addCoords=True)
    try:
        AllChem.MMFFOptimizeMolecule(m, maxIters=200)
    except Exception:
        pass
    prep = MoleculePreparation(rigid_macrocycles=True)
    setups = prep.prepare(m)
    pdbqt, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        return False
    open(destino_pdbqt, "w", encoding="utf-8").write(pdbqt)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--entradas", default=",".join(SISTEMAS))
    ap.add_argument("--exhaustividad", type=int, default=8)
    ap.add_argument("--caja", type=int, default=TAMANO_CAJA,
                    help="lado de la caja en amstrong (el proyecto usa 22-24)")
    ap.add_argument("--rehacer-receptor", action="store_true",
                    help="vuelve a preparar el receptor aunque ya exista "
                         "(obligatorio si cambio extraer_receptor)")
    ap.add_argument("--salida", default=None,
                    help="nombre del JSON de resultados (por defecto, uno por caja)")
    ap.add_argument("--semillas", default="42,2026,777",
                    help="semillas del generador de Vina (el proyecto usa 42/2026/777). "
                         "Sin semilla fija el resultado no es reproducible: el RMSD del "
                         "isoproterenol salio 1,89 A en una corrida y 2,55 A en la siguiente.")
    args = ap.parse_args()

    resultados = []
    for entrada in [e.strip() for e in args.entradas.split(",") if e.strip()]:
        codigo, cadena, nombre = SISTEMAS[entrada]
        log("\n=== %s  (%s, ligando %s, cadena %s) ===" % (entrada, nombre, codigo, cadena))
        pdb = _buscar_pdb(entrada)
        if pdb is None:
            log("  falta el PDB de %s" % entrada)
            continue
        if cadena is None:
            cadena, resseq = elegir_copia_trp32(pdb, codigo)
        else:
            resseq = elegir_copia(pdb, codigo, cadena)
        if not resseq:
            log("  no se encontro ninguna copia de %s" % codigo)
            continue
        lig_pdb = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
        rec_pdb = os.path.join(BASE, "receptores", "%s_proteina.pdb" % entrada)
        rec_base = os.path.join(BASE, "receptores", entrada)
        sdf = os.path.join(BASE, "ligands", "%s_ideal.sdf" % codigo)
        lig_pdbqt = os.path.join(BASE, "ligands", "%s_dock.pdbqt" % entrada)

        if not os.path.exists(sdf):
            cmd = ["curl", "-s", "--max-time", "40", "-o", sdf,
                   "https://files.rcsb.org/ligands/download/%s_ideal.sdf" % codigo]
            subprocess.run(cmd, capture_output=True)
            if not os.path.exists(sdf) or os.path.getsize(sdf) < 50:
                log("  no se pudo bajar el ligando ideal %s" % codigo)
                continue

        extraer_ligando(pdb, codigo, cadena, resseq, lig_pdb)
        n = extraer_receptor(pdb, rec_pdb,
                             descartar=RESIDUOS_DESCARTADOS.get(entrada, ()))
        log("  receptor: %d atomos | ligando: %d atomos" %
            (n, sum(1 for l in open(lig_pdb) if l.startswith("ATOM"))))

        cx, cy, cz = centroide(lig_pdb)
        log("  centro de caja (centroide del cristal): %.1f, %.1f, %.1f" % (cx, cy, cz))

        # receptor PDBQT (meeko escribe <base>.pdbqt sin residuos flexibles)
        rec_pdbqt = rec_base + ".pdbqt"
        if args.rehacer_receptor or not os.path.exists(rec_pdbqt):
            # -a (allow_bad_res): descarta residuos que meeko no sabe emparejar con
            # una plantilla. Sin esto fallan 8GSQ, 6A9O y 5YTO, que traen restos
            # terminales («residuo 0») y una serina incompleta (5YTO C:142SER).
            rr = subprocess.run([PYTHON, MK_RECEPTOR, "--read_pdb", rec_pdb,
                                 "-o", rec_base, "-p", "-j", "-a",
                                 "--default_altloc", "A"],
                                capture_output=True, text=True)
            if not os.path.exists(rec_pdbqt):
                log("  meeko: %s" % ((rr.stdout or "") + (rr.stderr or ""))[-400:])
        if not os.path.exists(rec_pdbqt):
            log("  FALLO: no se preparo el receptor")
            continue

        if not os.path.exists(lig_pdbqt):
            if not preparar_ligando(sdf, lig_pdbqt):
                log("  FALLO: no se preparo el ligando")
                continue

        por_semilla = []
        for semilla in [int(s) for s in args.semillas.split(",") if s.strip()]:
            out = os.path.join(BASE, "out", "%s_caja%d_e%d_s%d.pdbqt"
                               % (entrada, args.caja, args.exhaustividad, semilla))
            cmd = [VINA, "--receptor", rec_pdbqt, "--ligand", lig_pdbqt,
                   "--center_x", "%.3f" % cx, "--center_y", "%.3f" % cy,
                   "--center_z", "%.3f" % cz,
                   "--size_x", str(args.caja), "--size_y", str(args.caja),
                   "--size_z", str(args.caja),
                   "--exhaustiveness", str(args.exhaustividad),
                   "--num_modes", "9", "--seed", str(semilla), "--out", out]
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            # Vina devuelve la salida por stdout y algunos fallos solo se ven ahi.
            if "Parse error" in (r.stdout or "") or "Parse error" in (r.stderr or ""):
                log("  FALLO en Vina (semilla %d): error de parseo del receptor/ligando" % semilla)
            if not os.path.exists(out) or os.path.getsize(out) < 100:
                log("  FALLO en Vina (semilla %d): %s"
                    % (semilla, (r.stderr or r.stdout)[-200:]))
                continue
            af = None
            for l in open(out):
                if l.startswith("REMARK VINA RESULT:"):
                    af = float(re.search(r"([-+]?\d+\.\d+)", l).group(1))
                    break
            pose = molecula_dockeada(out)
            val = rmsd_en_sitio(lig_pdb, pose)
            log("  semilla %-5d afinidad %.2f kcal/mol | RMSD vs cristal: %s"
                % (semilla, af if af is not None else float("nan"),
                   ("%.2f A" % val) if val is not None else "no calculable"))
            por_semilla.append({"semilla": semilla, "afinidad": af,
                                "rmsd": round(val, 2) if val is not None else None})

        rmsds = [x["rmsd"] for x in por_semilla if x["rmsd"] is not None]
        mejor = min(rmsds) if rmsds else None
        mediana = sorted(rmsds)[len(rmsds) // 2] if rmsds else None
        pases = sum(1 for v in rmsds if v <= 2.0)
        log("  -> mejor %.2f A | mediana %.2f A | pasa %d de %d semillas" %
            (mejor if mejor is not None else float("nan"),
             mediana if mediana is not None else float("nan"), pases, len(rmsds)))
        resultados.append({
            "entrada": entrada, "ligando": codigo, "nombre": nombre,
            "caja": args.caja, "cadena": cadena, "resSeq_elegido": resseq,
            "centro": [round(cx, 2), round(cy, 2), round(cz, 2)],
            "por_semilla": por_semilla,
            "rmsd_mejor": mejor, "rmsd_mediana": mediana,
            "semillas_que_pasan": pases,
            "criterio": "PASA" if (mediana is not None and mediana <= 2.0) else "FALLA",
        })

    destino = os.path.join(BASE, args.salida or ("resultado_redocking_caja%d.json" % args.caja))
    with open(destino, "w", encoding="utf-8") as f:
        json.dump(resultados, f, indent=2, ensure_ascii=False)
    log("\n=== RESUMEN (caja de %d A) ===" % args.caja)
    for r in resultados:
        log("  %s (%s)  mejor %.2f A | mediana %.2f A | %d/%d semillas pasan -> %s"
            % (r["entrada"], r["nombre"],
               r["rmsd_mejor"] if r["rmsd_mejor"] is not None else float("nan"),
               r["rmsd_mediana"] if r["rmsd_mediana"] is not None else float("nan"),
               r["semillas_que_pasan"], len(r["por_semilla"]), r["criterio"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
