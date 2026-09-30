#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara el receptor de hCA2 del cristal 6T4P CON SU ZINC, la caja y el ligando de
referencia del cristal. Con esto empieza el redocking de control.

POR QUE hCA2 ES DISTINTA Y HAY QUE TENER CUIDADO CON EL ZINC
----------------------------------------------------------
hCA2 es una enzima de zinc, y el farmaco (la sulfonamida) se une al METAL: no al
carbono, al zinc. Tres cosas que en las otras dianas no importan y aqui son el
proyecto entero:

1. **El zinc tiene que estar en el receptor.** Si al limpiar el receptor se cae el
   HETATM del Zn, todo lo medido carece de sentido: estariamos puntuando un sitio que
   no existe. Es el mismo caso que el zinc de SOD1, que ya se resolvio en este
   proyecto. Aqui se comprueba DESPUES de preparar, sobre el pdbqt final, que hay un
   Zn y que su distancia al nitrogeno del farmaco es la de coordinacion.
2. **His94, His96, His119 y Thr199.** Son los cuatro que sujetan el zinc. Los tres
   histidinas donan un nitrogeno de la cadena lateral y la treonina dona el
   oxigeno. Si el protonado sale mal, el Zn queda con la carga que no es y Vina
   coloca el farmaco en otro sitio. Por eso el script escribe el pdbqt del receptor
   y ADEMAS el pdb sin protonar, para poder comparar.
3. **La carga de Zn en PDBQT.** Zn(II) va como ZN con carga +2. Se comprueba en el
   pdbqt, no se supone.

LA CAJA
-------
24 A de lado, centrada en el centroide de los atomos PESADOS del farmaco del
cristal. Es la convencion de DUD-E y la del proyecto, y con un farmaco de 14 atomos
pesados y un radio de 3,6 A el ligando esta holgadamente dentro: no lo va a recortar
el borde, que era una de las tres medidas que abrio el paso anterior.

OJO: NO SE USA EL MOL2 DE NADIE
-------------------------------
El farmaco sale del .cif del propio cristal 6T4P, con sus coordenadas reales. Un
mol2 de una base de datos tendria una conformacion inventada, y el RMSD se compararia
contra una pose que nunca existio.

Salida (en `analysis/calibracion_hca2/`):
    receptor_hca2.pdb      proteina normalizada, sin H y con columna de elemento
    receptor_hca2.pdbqt    preparado con meeko, con el zinc dentro
    ligando_cristal.pdb    los atomos del farmaco tal cual estan en el cristal
    ligando_cristal.pdbqt  el farmaco preparado para acoplar (se carga al docking)
    caja_hca2.json         centro y tamano
    preparar_receptor_hca2.log

Uso:
    python preparar_receptor_hca2.py
"""
from __future__ import annotations

import gzip
import json
import os
import subprocess

import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
DIR = os.path.join(BASE, "calibracion_hca2")
CACHE = os.path.join(DIR, "_cache_hca2")

PDB_ORIG = "6T4P"
COMP = "MHK"
NOMBRE = "naphthalene-1-sulfonamide"
SMILES = "NS(=O)(=O)c1cccc2ccccc12"

PDB_LIMPIO = os.path.join(DIR, "receptor_hca2.pdb")
PDBQT = os.path.join(DIR, "receptor_hca2.pdbqt")
LIG_PDB = os.path.join(DIR, "ligando_cristal.pdb")
LIG_PDBQT = os.path.join(DIR, "ligando_cristal.pdbqt")
CAJA = os.path.join(DIR, "caja_hca2.json")
LOG = os.path.join(DIR, "preparar_receptor_hca2.log")

TAMANO_CAJA = 24
RADIO_BOLSILLO = 8.0
CORTE_ZN = 2.6

PYDIR = os.path.join(os.environ["LOCALAPPDATA"], "Python", "pythoncore-3.14-64")
MK_RECEPTOR = os.path.join(PYDIR, "Scripts", "mk_prepare_receptor.exe")
MK_LIGAND = os.path.join(PYDIR, "Scripts", "mk_prepare_ligand.exe")
OBABEL = os.path.join(PYDIR, "Scripts", "obabel.exe")

# Las 20 aminoacidas estandar mas las variantes de histidina que meeko trae en sus
# plantillas. Meeko se para con cualquier otra y hay que saberlo antes, no despues.
ESTANDAR = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "CYX", "CYM", "ASH", "GLH", "LYN",
}

lineas = []


def log(m):
    print(m, flush=True)
    lineas.append(m)



def linea_pdb(registro, serial, nombre, resname, chain, resseq,
              x, y, z, b, elemento, altloc=" "):
    """Escribe UNA linea PDB colococating cada campo en su columna de verdad.

    POR QUE ESTO Y NO UN FORMATO CON %

    Con "%8.3f%8.3f%8.3f" hay que contar caracteres de memoria para saber en que
    columna cae el elemento (77-78), y si la cadena se queda corta el elemento se
    sale del registro. Entonces meeko lo lee como vacio y para con
    "Element '' not found", un error que no senala que el problema es el formato.
    Con columnas explicitas no hay nada que contar: el elemento va SIEMPRE en
    77-78 porque para ahi.

        1-6   registro      13-16  nombre      43-54  x y z
        7-11  serial        17-20  resname     55-60  occupancy
        12    altloc        22-26  resseq      61-66  factor B
        (13 es el hueco del elemento de dos letras)   77-78  elemento
    """
    nombre = nombre.strip()
    # Un elemento de dos letras tiene que empezar en la columna 13; por eso el
    # nombre se centra en 13-16 con un hueco delante.
    nombre = (" " + nombre).rjust(4)[:4] if len(nombre) >= 3 else nombre.rjust(4)
    # El hueco entre el nombre del atomo (13-16) y el del residuo (17-20) es
    # OBLIGATORIO. Sin el, "N" y "HIS" salen pegados como "NHIS" y meeko lee
    # un residuo llamado NHIS, luego no encuentra su plantilla y acaba
    # interpreteando nombres de atomos de otros sitios (por eso aparecian
    # HA, SN, SP y LA como si fueran residuos).
    s = (
        "%-6s%5d %-4s %-3s %s%4s    %8.3f%8.3f%8.3f%6.2f%6.2f"
        % (registro, serial, nombre, resname, chain, resseq, x, y, z,
           min(b, 99.99), 1.00)
    )
    s = s.ljust(76)
    return s + elemento.rjust(2)


def cargar_cif(pdb_id):
    ruta = os.path.join(CACHE, "%s.cif.gz" % pdb_id)
    if not os.path.exists(ruta):
        raise SystemExit("falta el .cif en cache: corre elegir_control_hca2.py antes")
    with gzip.open(ruta, "rt", encoding="utf-8", errors="ignore") as h:
        return MMCIF2Dict(h)


def main():
    log("PREPARAR EL RECEPTOR DE hCA2 DEL CRISTAL %s" % PDB_ORIG)
    log("   %s" % __import__("time").strftime("%Y-%m-%d %H:%M"))
    os.makedirs(DIR, exist_ok=True)
    d = cargar_cif(PDB_ORIG)

    g = d["_atom_site.group_PDB"]
    comp = d["_atom_site.label_comp_id"]
    atomo = d["_atom_site.label_atom_id"]
    elem = d.get("_atom_site.type_symbol") or [a[:1] for a in atomo]
    residuo = d.get("_atom_site.label_seq_id") or d.get("_atom_site.auth_seq_id") or \
        ["?"] * len(g)
    xyz = np.stack([np.array([float(v) for v in d["_atom_site.Cartn_x"]]),
                    np.array([float(v) for v in d["_atom_site.Cartn_y"]]),
                    np.array([float(v) for v in d["_atom_site.Cartn_z"]])], axis=1)
    alt = d.get("_atom_site.label_alt_id") or ["."] * len(g)
    bs = d.get("_atom_site.B_iso_or_equiv") or ["0"] * len(g)
    mod = d.get("_atom_site.pdbx_PDB_model_num") or ["1"] * len(g)

    def Atom(i):
        # Todo lo que sale del .cif es texto. Se normaliza aqui una vez, porque si no
        # cada sitio que lo imprima con %d revienta el script, y el error que sale no
        # senala la causa: senala la linea de un formateo.
        num = str(residuo[i]).strip()
        if num in (".", "?", ""):
            num = "1"
        b = str(bs[i]).strip()
        if b in (".", "?", ""):
            b = 0.0
        else:
            try:
                b = float(b)
            except ValueError:
                b = 0.0
        return {"x": float(xyz[i][0]), "y": float(xyz[i][1]),
                "z": float(xyz[i][2]),
                "elem": str(elem[i]).strip().upper(), "name": atomo[i].strip(),
                "res": comp[i].strip(), "num": num, "alt": str(alt[i]),
                "b": b}

    # ---- el farmaco del cristal: coordenadas de verdad, un solo modelo/altloc ----
    idx_lig = [i for i in range(len(g))
               if g[i] == "HETATM" and comp[i] == COMP
               and mod[i] == mod[0] and alt[i] in (".", "?")
               and str(elem[i]).upper() != "H"]
    if not idx_lig:
        log("no encontre el farmaco %s en el cristal" % COMP)
        return 1
    L = xyz[idx_lig]
    centro = L.mean(axis=0)
    log("")
    log("FARMACO DEL CRISTAL")
    log("   %s (%s): %d atomos pesados" % (COMP, NOMBRE, len(idx_lig)))
    log("   centroide %.2f, %.2f, %.2f" % tuple(centro))
    log("   radio alrededor del centroide %.2f A"
        % float(np.linalg.norm(L - centro, axis=1).max()))

    # el PDB del farmaco, tal cual, para tener la pose de referencia del RMSD
    with open(LIG_PDB, "w", encoding="utf-8") as f:
        for k, i in enumerate(idx_lig, 1):
            a = Atom(i)
            nombre = a["name"][:4].ljust(4)
            f.write(linea_pdb("HETATM", k, nombre, a["res"], "A", "1",
                              a["x"], a["y"], a["z"], a["b"], a["elem"]) + "\n")
        f.write("END\n")
    log("   pose de referencia: ligando_cristal.pdb")

    # ---- el zinc: esta es la parte que hay que mirar dos veces ----
    idx_zn = [i for i in range(len(g))
              if comp[i] in ("ZN", "ZN2") and mod[i] == mod[0]
              and alt[i] in (".", "?")]
    if not idx_zn:
        log("")
        log("ATENCION: no hay ZN en el archivo. El farmaco se une al metal, y sin el")
        log("           metal este control no mide nada. Se para aqui.")
        return 1
    Zn = xyz[idx_zn]
    log("")
    log("EL ZINC")
    log("   %d atomo(s) de Zn en el archivo" % len(Zn))
    zn_cerca = Zn[np.linalg.norm(Zn - centro, axis=1) < 6.0]
    log("   %d a menos de 6 A del farmaco" % len(zn_cerca))
    if len(zn_cerca) == 0:
        log("   el Zn mas cercano esta a %.2f A del farmaco: no es el sitio"
            % float(np.linalg.norm(Zn - centro, axis=1).min()))
        return 1
    # distancia del N del sulfonamida al Zn: tiene que ser de coordinacion
    d_nzn = None
    for i in idx_lig:
        if str(elem[i]).upper() == "N":
            v = float(np.linalg.norm(zn_cerca - xyz[i], axis=1).min())
            d_nzn = v if d_nzn is None else min(d_nzn, v)
    log("   N del sulfonamida a %.2f A del Zn %s"
        % (d_nzn, "(coordinacion, correcto)" if d_nzn <= CORTE_ZN
           else "DEMASIADO LEJOS: no esta coordinando"))
    if d_nzn > CORTE_ZN:
        log("   el farmaco de este cristal no toca el Zn: revisa el control antes de")
        log("   acoplar. Puede que la sulfonamida este en otro bolsillo.")

    # ---- los cuatro que sujetan el zinc ----
    log("")
    log("LOS RESIDUOS QUE SUJETAN EL ZINC (deben estar a 2,0-2,4 A)")
    proteina = [i for i in range(len(g)) if g[i] == "ATOM"
                and mod[i] == mod[0] and alt[i] in (".", "?", "A")
                and str(elem[i]).upper() != "H"]
    P = xyz[proteina]
    objetivo = zn_cerca[0]
    for i in proteina:
        d = float(np.linalg.norm(xyz[i] - objetivo))
        if d < 2.6:
            a = Atom(i)
            log("   %s%s %-4s  %-4s  %.2f A"
                % (a["res"], a["num"], a["name"], a["elem"], d))

    # ---- receptor: proteina sin H, con columna de elemento, mas el Zn ----
    log("")
    n_prot = n_zn_out = 0
    raros = {}
    with open(PDB_LIMPIO, "w", encoding="utf-8") as f:
        for i in proteina:
            a = Atom(i)
            res = a["res"].upper()
            if res not in ESTANDAR:
                raros[res] = raros.get(res, 0) + 1
            nombre = a["name"]
            if not nombre[:1].isalpha() or len(nombre.strip()) > 4:
                nombre = nombre[-4:]
            nombre = nombre[:4].rjust(4)
            f.write(linea_pdb("ATOM", n_prot + 1, nombre, res, "A", a["num"],
                              a["x"], a["y"], a["z"], a["b"], a["elem"]) + "\n")
            n_prot += 1
        # El zinc entra como HETATM con numero de residuo propio, y se COMENTA que
        # es un metal: no es una proteina y va con carga +2, que se comprueba abajo
        # en el pdbqt.
        for k, i in enumerate(idx_zn, 1):
            a = Atom(i)
            # El numero de residuo se imprime como texto a proposito: en un .cif
            # puede venir "901" como cadena y con %d revienta el script entero.
            # Los marcadores y los argumentos tienen que casar uno a uno, y el
            # numero de residuo sale del dato, no de un contador.
            num_zn, bzn = a["num"], a["b"]
            x, y, z = a["x"], a["y"], a["z"]
            # La linea tiene que ocupar EXACTAMENTE 78 columnas: la del elemento
            # va en 77-78 y si la cadena se queda corta, el elemento acaba fuera
            # del registro, meeko lo lee como vacio y para con "Element '' not
            # found". Se le anaden los huecos que falten en vez de confiar en que
            # salen bien.
            f.write(linea_pdb("HETATM", 900 + k, "ZN", "ZN", "A", num_zn,
                              x, y, z, bzn, a["elem"]) + "\n")
            n_zn_out += 1
        f.write("END\n")
    log("RECEPTOR NORMALIZADO")
    log("   receptor_hca2.pdb: %d atomos de proteina + %d Zn"
        % (n_prot, n_zn_out))
    if raros:
        log("   AVISO: residuos fuera de la tabla de meeko: %s"
            % ", ".join("%s (%d at.)" % (k, v) for k, v in sorted(raros.items())))
    else:
        log("   residuos: todos en la tabla de meeko")

    # residuos del bolsillo, para dejarlo escrito
    dist = np.linalg.norm(P - centro, axis=1)
    residuos = {}
    for k in np.where(dist < RADIO_BOLSILLO)[0]:
        a = Atom(proteina[k])
        etiqueta = "%s%s" % (a["res"].capitalize(), a["num"])
        residuos[etiqueta] = min(residuos.get(etiqueta, 1e9), dist[k])
    log("   residuos a menos de %.1f A del farmaco: %d"
        % (RADIO_BOLSILLO, len(residuos)))
    for etiqueta in sorted(residuos, key=lambda x: residuos[x])[:10]:
        log("      %-8s %.2f A" % (etiqueta, residuos[etiqueta]))

    # ---- meeko: el receptor, CON caja ----
    log("")
    log("PREPARANDO EL RECEPTOR CON MEEKO")
    for p in (MK_RECEPTOR,):
        if not os.path.exists(p):
            log("FALTA %s" % p)
            return 1
    if os.path.exists(PDBQT):
        os.remove(PDBQT)
    cmd = [MK_RECEPTOR, "--read_pdb", PDB_LIMPIO, "-o", PDBQT[:-6], "-p", "-j",
           "-a", "--default_altloc", "A",
           "--box_center", "%.3f" % centro[0], "%.3f" % centro[1],
           "%.3f" % centro[2],
           "--box_size", str(TAMANO_CAJA), str(TAMANO_CAJA), str(TAMANO_CAJA)]
    log("   %s" % " ".join(os.path.basename(c) for c in cmd[:3]) + " ...")
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except subprocess.TimeoutExpired:
        log("   meeko se paso de tiempo")
        return 1
    if not os.path.exists(PDBQT):
        log("   FALLO preparando el receptor:\n%s"
            % ((r.stdout or "") + (r.stderr or ""))[-1800:])
        return 1
    log("   meeko: %s" % ((r.stdout or "").strip()[-300:] or "sin salida"))

    # ---- LA COMPROBACION QUE IMPORTA: el Zn sigue en el pdbqt ----
    log("")
    log("COMPROBANDO QUE EL ZINC SOBREVIVIO A LA PREPARACION")
    lineas_rec = [l for l in open(PDBQT, encoding="utf-8", errors="ignore")
                  if l.startswith(("ATOM", "HETATM"))]
    zn_rec = [l for l in lineas_rec if l[77:79].strip().upper() in ("ZN", "D", "Q")
              or l[12:16].strip().upper().startswith("ZN")]
    log("   receptor_hca2.pdbqt: %d atomos" % len(lineas_rec))
    log("   atomos de Zn en el pdbqt: %d" % len(zn_rec))
    if len(zn_rec) == 0:
        log("   FALLO: EL ZINC SE HA PERDIDO. Meeko se lo comio. Todo lo que se")
        log("   acople contra este receptor seria sobre un sitio que no existe.")
        log("   Se para aqui en vez de devolver un resultado bonito y falso.")
        return 1
    coord = []
    for l in zn_rec:
        coord.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
    log("   Zn en el pdbqt: %s" % ", ".join("(%.2f, %.2f, %.2f)" % c
                                            for c in coord))
    # Lo que hay que mirar NO es la distancia del Zn al CENTROIDE del farmaco, sino
    # al NITROGENO de la sulfonamida. El zinc se coordina a ese nitrogeno concreto,
    # y el centroide cae en el naftaleno, que esta en el otro extremo: dan 4,88 A
    # uno y 2,00 A el otro. Comparar con el centroide daba un aviso falso que
    # celebrations aparecio el Zn intacto.
    zn_pdb = np.array(coord[0])
    idx_n = [i for i in idx_lig if str(elem[i]).upper() == "N"]
    if idx_n:
        d_n_pdb = float(min(np.linalg.norm(zn_pdb - xyz[i]) for i in idx_n))
        log("   Zn (pdbqt) - N del sulfonamida (cristal): %.2f A %s"
            % (d_n_pdb, "COORDINACION, correcto" if d_n_pdb <= CORTE_ZN
               else "DEMASIADO LEJOS"))
        if d_n_pdb > CORTE_ZN:
            log("   AVISO: el Zn del pdbqt no queda a distancia de coordinacion del")
            log("          nitrogeno. Habria que revisar el receptor antes de acoplar.")
    # ademas, que el Zn no se haya movido al preparar
    zn_cif = Zn[0]
    mov = float(np.linalg.norm(zn_pdb - zn_cif))
    log("   desplazamiento del Zn al preparar: %.3f A %s"
        % (mov, "(intacto)" if mov < 0.01 else "(MOVIDO)"))
    if mov >= 0.01:
        log("   AVISO: meeko movio el zinc. Las coordenadas de referencia del RMSD")
        log("          y las del docking estarian en sitios distintos.")
    # carga del Zn
    try:
        mk = json.load(open(PDBQT[:-6] + ".json", encoding="utf-8"))
        log("   meeko escribio el json de indices: %d entradas" % len(mk))
    except Exception:  # noqa: BLE001
        pass

    # ---- el ligando del cristal, preparado para acoplar ----
    log("")
    log("PREPARANDO EL LIGANDO DEL CRISTAL")
    mol = Chem.MolFromSmiles(SMILES)
    if mol is None:
        log("   el SMILES del control no parsea: %s" % SMILES)
        return 1
    # Meeko exige HIDROGENOS EXPLICITOS: con H implicitos dice "RDKit molecule has
    # implicit Hs. Need explicit Hs" y no escribe el pdbqt. Y no es un capricho del
    # script: sin los H, el acoplamiento no sabe calcular bien los puentes de
    # hidrogeno, que es justamente de lo que aqui depende la union al zinc.
    mol = Chem.AddHs(mol)
    Chem.MolToMolFile(mol, os.path.join(DIR, "_lig_tmp.mol"))
    if os.path.exists(LIG_PDBQT):
        os.remove(LIG_PDBQT)
    cmd = [MK_LIGAND, "-i", os.path.join(DIR, "_lig_tmp.mol"), "-o", LIG_PDBQT]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if not os.path.exists(LIG_PDBQT):
        # meeko a veces no lee .mol: se pasa por sdf
        mol_b = Chem.Mol(mol)
        mol_b.SetProp("_Name", "MHK")
        w = Chem.SDWriter(os.path.join(DIR, "_lig_tmp.sdf"))
        w.write(mol_b)
        w.close()
        cmd = [MK_LIGAND, "-i", os.path.join(DIR, "_lig_tmp.sdf"), "-o", LIG_PDBQT]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
    if not os.path.exists(LIG_PDBQT):
        log("   FALLO preparando el ligando:\n%s"
            % ((r.stdout or "") + (r.stderr or ""))[-1200:])
        return 1
    log("   meeko: %s" % ((r.stdout or "").strip()[-300:] or "sin salida"))
    n_lig = sum(1 for l in open(LIG_PDBQT, encoding="utf-8", errors="ignore")
                if l.startswith(("ATOM", "HETATM")))
    log("   ligando_cristal.pdbqt: %d atomos" % n_lig)
    for aux in ("_lig_tmp.mol", "_lig_tmp.sdf"):
        q = os.path.join(DIR, aux)
        if os.path.exists(q):
            os.remove(q)

    # ---- la caja ----
    json.dump({"centro_caja": [round(float(x), 3) for x in centro],
               "tamano_caja": TAMANO_CAJA,
               "diana": "hCA2 (anhidrasa carbonica 2)",
               "pdb": PDB_ORIG,
               "ligando_cristal": "%s (%s, %d atomos pesados, Ki p6.88)"
                                  % (COMP, NOMBRE, len(idx_lig)),
               "smiles": SMILES,
               "receptor_pdbqt": "receptor_hca2.pdbqt",
               "receptor_pdb": "receptor_hca2.pdb",
               "ligando_pdb": "ligando_cristal.pdb",
               "ligando_pdbqt": "ligando_cristal.pdbqt",
               "zn_en_pdbqt": len(zn_rec),
               "d_n_zn": None if d_nzn is None else round(d_nzn, 2),
               "residuos_bolsillo": sorted(residuos)},
              open(CAJA, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    log("")
    log("caja escrita: caja_hca2.json (%.1f A, centro %.2f %.2f %.2f)"
        % (TAMANO_CAJA, *centro))
    log("receptor con zinc: receptor_hca2.pdbqt")
    log("todo listo para el redocking")

    with open(LOG, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
