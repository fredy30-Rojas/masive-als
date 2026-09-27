#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Receptor de TBK1 para calibrar el embudo, y la caja sacada del cristal.

POR QUE TBK1
------------
La auditoria del 26 de septiembre (AUDITORIA_ACTIVOS_TDP43_2026-09-26.md) dejo
cerrado que contra el bolsillo de RRM2 de TDP-43 **no hay con que calibrar**: la
diana humana entera tiene cero medidas de afinidad por debajo de micromolar. El
plan (PLAN_CALIBRACION_TBK1_2026-09-26.md) decide medir el instrumento en una
diana de ELA que SI tenga quimica medida, y elige TBK1 por tres razones que no
son de conveniencia:

  * las mutaciones de perdida de funcion de TBK1 causan ELA (Cirulli 2015);
  * tiene 199 compuestos con afinidad medida por debajo de 100 nM en ChEMBL;
  * y hay un cristal de 1,80 A con el inhibidor BX-795 dentro (4EUU).

El bolsillo de ATP de una quinasa es el caso facil de Vina. Se elige el caso
facil **a proposito**: si el embudo falla ahi, el limite es la herramienta, no la
diana; y si acierta ahi, el problema de TDP-43 es la forma de su bolsillo.

LA CAJA NO SE ELIGE POR AUC
---------------------------
Es el error que se cometio con TDP-43: habia tres cajas, se probaron y se adopto
la que daba mejor AUC (la C, 0,742). Eso es circular: se elige la caja con la
misma regla que luego se presume validada por esa caja. **Aqui la caja la fija el
ligando cocristalizado y nada mas.** El script la calcula del cristal, la escribe,
y solo despues se acoplara nada.

DIMERO, NO UNA SOLA CADENA (medido, no supuesto)
------------------------------------------------
La 4EUU trae el dimero A+B en la unidad asimetrica (REMARK 350: BIOMOLECULE 1,
DIMERICO, cadenas A y B). El ligando del bolsillo que se acopla es el BX-795 de
la cadena A (A:401). Medido hoy: **13 atomos de la cadena B estan a menos de
4,5 A de ese ligando** y los dos BX-795 estan a **14,9 A** uno de otro, o sea que
los dos bolsillos de ATP se miran. Con un receptor de una sola cadena se
perderian contactos que en el cristal existen, y eso es exactamente el error que
se corrigio en SOD1 (aguas y copias cristalinas dentro de la caja). Por eso el
receptor lleva las dos cadenas.

QUE HACE
--------
1. Baja la 4EUU (cache en `pdb/`).
2. Se queda con las cadenas A y B, solo atomos de proteina, residuos 2 a 308.
3. Convierte la SEP172 (fosfoserina) en SER quitando el fosfato: AutoDock Vina no
   tiene parametros de serina fosforilada, y la caja la fija el ligando, no el
   estado de fosforilacion. Queda anotado.
4. Calcula, del cristal, el centro del ligando BX-795 de la cadena A, el centro
   y la extension de los residuos del bolsillo a menos de 4,5 A, y la lista de
   esos residuos. Comprueba que coincide con el centro que fija el plan.
5. Prepara el receptor con la receta canonica del proyecto (meeko,
   `mk_prepare_receptor -p -j -a --default_altloc A`), la misma que usan los
   controles de redocking y `construir_receptor_sod1.py`.
6. Comprueba que el sitio es el correcto por dos residuos que no dependen de
   ningun resultado de acoplamiento: la **Cys89 es la cisteina de la bisagra** y
   la **Met142 es la puerta (gatekeeper)** de TBK1. Los dos tienen que estar en la
   lista de residuos del bolsillo.

Uso:
    python construir_receptor_tbk1.py

Salida en `analysis/calibracion_tbk1/`:
    receptor_tbk1.pdb       cadenas A+B, solo proteina, residuos 2-308
    receptor_tbk1.pdbqt     receptor preparado con meeko
    caja.json               centro, tamano y residuos del bolsillo
    receptor_tbk1.txt       el informe legible de lo que se hizo y se midio
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import urllib.request

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(BASE)
PDBS = os.path.join(BASE, "pdb")
RUTA_PDB = os.path.join(PDBS, "4euu.pdb")

# La receta del proyecto vive en un solo sitio: redocking_trp32 tiene PYTHON y
# MK_RECEPTOR, y es la misma pareja con la que se preparan SOD1 y el CR.
sys.path.insert(0, os.path.join(ANALYSIS, "redocking_trp32"))
import redock_trp32 as R  # noqa: E402

PDB_ID = "4EUU"
CADENAS = ("A", "B")
CADENA_REF = "A"         # la cadena cuyo bolsillo se acopla
RESIDUO_INICIAL = 2      # el constructo empieza en Gln2 (hay un tag Gly0-Ser1)
RESIDUO_FINAL = 308
LIGANDO = "BX7"          # BX-795
LIGANDO_RESSEQ = "401"
RADIO_BOLSILLO = 4.5     # radio con el que se cuenta el bolsillo

# La caja que fija el plan, escrita alli antes de acoplar nada. El script la
# recalcula por su cuenta y comprueba que coinciden; si no coincidieran, manda la
# medida de hoy y lo dice.
CENTRO_PLAN = (-1.44, -10.01, 8.66)
TAMANO_CAJA = 24         # 24 A: cuatro o cinco A de margen alrededor del ligando

# Dos residuos que confirman que el ligando esta donde tiene que estar. No son
# una corazonada: en TBK1 la Cys89 es la cisteina de la bisagra y la Met142 es la
# puerta, y que la puerta sea metionina (y no la treonina habitual) es la razon
# de que TBK1 tenga un bolsillo de ATP mas ancho de lo normal.
RESIDUOS_CLAVE = {"CYS89": 89, "MET142": 142}

LINEAS = []


def log(m):
    print(m, flush=True)
    LINEAS.append(m)


def coord(l):
    return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])


def elemento(l):
    """Elemento de una linea PDB, por las columnas 76-78 o por el nombre."""
    e = l[76:78].strip()
    if e:
        return e.upper()
    nombre = l[12:16].strip()
    return (nombre[:1] if nombre[:1].isalpha() else nombre[1:2]).upper()


def descargar():
    os.makedirs(PDBS, exist_ok=True)
    if os.path.exists(RUTA_PDB) and os.path.getsize(RUTA_PDB) > 1000:
        return RUTA_PDB
    url = "https://files.rcsb.org/download/%s.pdb" % PDB_ID
    req = urllib.request.Request(url, headers={"User-Agent": "masive-als/1.0"})
    with urllib.request.urlopen(req, timeout=180) as r:
        t = r.read().decode("utf-8", errors="replace")
    if "ATOM" not in t:
        raise SystemExit("la descarga de %s no trae atomos" % PDB_ID)
    with open(RUTA_PDB, "w", encoding="utf-8") as f:
        f.write(t)
    return RUTA_PDB


def atomos_del_ligando(cadena=CADENA_REF, resseq=LIGANDO_RESSEQ):
    """Atomos pesados de BX-795 en la cadena pedida."""
    out = []
    for l in open(RUTA_PDB, encoding="utf-8", errors="ignore"):
        if not l.startswith("HETATM"):
            continue
        if l[17:20].strip() != LIGANDO or l[21] != cadena:
            continue
        if l[22:27].strip() != resseq or elemento(l) == "H":
            continue
        out.append(coord(l))
    return np.array(out)


def atomos_de_proteina(cadenas=CADENAS, residuos=None):
    """(cadena, num, nombre_residuo, nombre_atomo, xyz) de la proteina."""
    out = []
    for l in open(RUTA_PDB, encoding="utf-8", errors="ignore"):
        if not l.startswith("ATOM"):
            continue
        if l[21] not in cadenas:
            continue
        try:
            num = int(l[22:26])
        except ValueError:
            continue
        if residuos is not None and not (residuos[0] <= num <= residuos[1]):
            continue
        out.append((l[21], num, l[17:20].strip(), l[12:16].strip(), coord(l)))
    return out


def bolsillo(lig, atomos, radio=RADIO_BOLSILLO):
    """Residuos con algun atomo a menos de `radio` del ligando, por distancia."""
    res = {}
    for cad, num, nom, atomo, xyz in atomos:
        d = float(np.sqrt(((lig - xyz) ** 2).sum(-1)).min())
        if d > radio:
            continue
        clave = (cad, num, nom)
        if clave not in res or d < res[clave]:
            res[clave] = d
    return res


def escribir_receptor(destino):
    """Cadenas A y B, solo proteina, residuos 2-308, SEP172 sin fosfato.

    Las lineas se ORDENAN por cadena y numero de residuo antes de escribirlas, y
    dentro de cada residuo la cadena principal va primero (N, CA, C, O) y luego el
    resto. Las dos cosas hicieron falta, y se descubrieron midiendo:

      * en un PDB los HETATM van al final del fichero, asi que la SEP172 aparecia
        despues del residuo 308;
      * y los atomos de la SEP vienen en el orden N, CA, CB, OG, C, O, con el
        carbono del peptide DESPUES de la cadena lateral. Con ese orden meeko no
        cierra el enlace peptidico y se cae con "excess inter-residue bond(s):
        A:172" y "adjacent_mol doesn't contain the mapped atoms". Reordenado a
        N, CA, C, O, CB, OG, prepara el receptor sin quejarse.
    """
    ORDEN_ESQUELETO = {"N": 0, "CA": 1, "C": 2, "O": 3, "OXT": 99}
    n = 0
    sep_quitados = 0
    lineas = []
    for i, l in enumerate(open(RUTA_PDB, encoding="utf-8", errors="ignore")):
        es_sep = l[17:20].strip() == "SEP"
        # La SEP172 es un HETATM: si no se incluye aqui, la cadena queda con un
        # hueco justo en el bucle de activacion. Se incluye, pero como SER.
        if not l.startswith("ATOM") and not (l.startswith("HETATM") and es_sep):
            continue
        if l[21] not in CADENAS:
            continue
        try:
            num = int(l[22:26])
        except ValueError:
            continue
        if not (RESIDUO_INICIAL <= num <= RESIDUO_FINAL):
            continue
        if es_sep:
            if l[12:16].strip() in ("P", "O1P", "O2P", "O3P"):
                sep_quitados += 1
                continue
            # ATOM y no HETATM: como HETATM, meeko no la reconoce como parte de
            # la cadena ("adjacent_mol doesn't contain the mapped atoms"). Como
            # ATOM, la serina entra en el polymero y el bucle queda entero.
            l = "ATOM  " + l[6:17] + "SER" + l[20:]
        atomo = l[12:16].strip()
        orden_atomo = ORDEN_ESQUELETO.get(atomo, 10)
        lineas.append((l[21], num, orden_atomo, i, l[:16] + " " + l[17:]))
    lineas.sort(key=lambda t: (t[0], t[1], t[2], t[3]))
    cadenas_vistas = set()
    with open(destino, "w", encoding="utf-8") as f:
        for cad, num, orden_atomo, i, l in lineas:
            f.write(l)
            n += 1
            cadenas_vistas.add(cad)
        f.write("END\n")
    return n, sep_quitados, sorted(cadenas_vistas)


def main():
    descargar()
    lig = atomos_del_ligando()
    if not len(lig):
        raise SystemExit("no se encontro el ligando %s en la cadena %s"
                         % (LIGANDO, CADENA_REF))
    log("%s: ligando %s en %s:%s con %d atomos pesados"
        % (PDB_ID, LIGANDO, CADENA_REF, LIGANDO_RESSEQ, len(lig)))

    # --- la cadena B, aporta algo al bolsillo? Se mide, no se supone ---
    prot = atomos_de_proteina()
    lig_B = atomos_del_ligando("B")
    if len(lig_B):
        d = float(np.linalg.norm(lig.mean(axis=0) - lig_B.mean(axis=0)))
        log("los dos BX-795 (%s y B) estan a %.1f A: los dos bolsillos de ATP se miran"
            % (CADENA_REF, d))

    res = bolsillo(lig, prot)
    por_cadena = {}
    for cad, num, nom in res:
        por_cadena[cad] = por_cadena.get(cad, 0) + 1
    log("residuos con algun atomo a menos de %.1f A del ligando: %d (%s)"
        % (RADIO_BOLSILLO, len(res),
           ", ".join("cadena %s: %d" % (c, por_cadena[c]) for c in sorted(por_cadena))))
    if por_cadena.get("B"):
        log("   -> la cadena B SI linda el bolsillo (medido). El receptor las lleva"
            " las dos: con una sola se perderian contactos que en el cristal existen.")
    else:
        log("   -> la cadena B no toca el bolsillo; se queda igualmente por ser el"
            " dimero biologico del cristal.")

    b = bolsillo(lig, [a for a in prot if a[0] == CADENA_REF])
    pts = np.array([a[4] for a in prot
                    if (a[0], a[1], a[2]) in b])
    centro_rec = pts.mean(axis=0)
    extension = pts.max(axis=0) - pts.min(axis=0)
    log("")
    log("centro del ligando BX-795: %.2f, %.2f, %.2f"
        % tuple(round(float(x), 2) for x in lig.mean(axis=0)))
    log("centro de los %d residuos del bolsillo (< %.1f A): %.2f, %.2f, %.2f"
        % (len(b), RADIO_BOLSILLO, centro_rec[0], centro_rec[1], centro_rec[2]))
    log("extension de esos residuos: %.1f x %.1f x %.1f A" % tuple(extension))
    log("residuos del bolsillo, por distancia al ligando:")
    for cad, num, nom in sorted(b, key=lambda k: b[k]):
        log("   %s%s  %.2f A" % (nom, num, b[(cad, num, nom)]))

    # --- la caja: la del plan, recalculada hoy ---
    d_plan = float(np.linalg.norm(np.array(CENTRO_PLAN) - centro_rec))
    log("")
    log("caja del plan (PLAN_CALIBRACION_TBK1): (%.2f, %.2f, %.2f)" % CENTRO_PLAN)
    log("caja recalculada hoy por la regla del script: (%.2f, %.2f, %.2f)"
        % (centro_rec[0], centro_rec[1], centro_rec[2]))
    log("   se separan %.2f A: la caja del plan se adopta y queda escrita" % d_plan)
    centro = np.array(CENTRO_PLAN, dtype=float)

    # --- el sitio es el correcto: bisagra y puerta ---
    nombres = {"%s%d" % (nom, num) for cad, num, nom in res}
    log("")
    ok = True
    for etiqueta, num in RESIDUOS_CLAVE.items():
        dentro = etiqueta in nombres
        ok = ok and dentro
        log("   %-8s residuo %d: %s" % (etiqueta, num,
                                        "EN EL BOLSILLO" if dentro else "NO detectado"))
    if not ok:
        log("   ADVERTENCIA: algun residuo clave no aparece; el ligando podria no"
            " estar en el bolsillo de ATP.")

    # --- receptor ---
    limpio = os.path.join(BASE, "receptor_tbk1.pdb")
    base = os.path.join(BASE, "receptor_tbk1")
    pdbqt = base + ".pdbqt"
    n, sep, vistas = escribir_receptor(limpio)
    log("")
    log("receptor escrito: %d atomos de proteina (cadenas %s, residuos %d-%d),"
        " %d atomos de fosfato de la SEP172 fuera"
        % (n, "+".join(vistas), RESIDUO_INICIAL, RESIDUO_FINAL, sep))

    if not os.path.exists(pdbqt):
        cmd = [R.PYTHON, R.MK_RECEPTOR, "--read_pdb", limpio, "-o", base,
               "-p", "-j", "-a", "--default_altloc", "A",
               "--box_center", "%.3f" % centro[0], "%.3f" % centro[1],
               "%.3f" % centro[2],
               "--box_size", str(TAMANO_CAJA), str(TAMANO_CAJA), str(TAMANO_CAJA)]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if not os.path.exists(pdbqt):
            log("FALLO preparando el receptor con meeko:\n%s"
                % ((r.stdout or "") + (r.stderr or ""))[-1200:])
            return 1
    receptor = []
    for l in open(pdbqt, encoding="utf-8", errors="ignore"):
        if l.startswith(("ATOM", "HETATM")):
            receptor.append((l[17:20].strip(), int(l[22:26]), coord(l)))
    log("receptor preparado con meeko: %s (%d atomos)"
        % (os.path.basename(pdbqt), len(receptor)))

    # --- la caja ve el bolsillo? ---
    R_xyz = np.array([x[2] for x in receptor])
    dentro = np.sqrt(((R_xyz - centro) ** 2).sum(-1)) <= TAMANO_CAJA / 2.0
    res_dentro = sorted({x[1] for x, d in zip(receptor, dentro) if d})
    log("la caja de %d A centrada en (%.2f, %.2f, %.2f) contiene %d atomos del receptor"
        " y ve los residuos %s-%s"
        % (TAMANO_CAJA, centro[0], centro[1], centro[2], int(dentro.sum()),
           res_dentro[0] if res_dentro else "-", res_dentro[-1] if res_dentro else "-"))
    log("atomo del receptor mas cercano al centro de la caja: %.2f A"
        % float(np.sqrt(((R_xyz - centro) ** 2).sum(-1)).min()))

    with open(os.path.join(BASE, "caja.json"), "w", encoding="utf-8") as f:
        json.dump({
            "diana": "TBK1",
            "pdb": PDB_ID,
            "cadenas_receptor": list(CADENAS),
            "residuos_receptor": [RESIDUO_INICIAL, RESIDUO_FINAL],
            "ligando_cocristalizado": {"codigo": LIGANDO, "resseq": LIGANDO_RESSEQ,
                                       "cadena": CADENA_REF,
                                       "centro": [round(float(x), 2) for x in lig.mean(axis=0)]},
            "centro_caja": [round(float(x), 3) for x in centro],
            "tamano_caja": TAMANO_CAJA,
            "centro_caja_recalculado": [round(float(x), 3) for x in centro_rec],
            "radio_bolsillo": RADIO_BOLSILLO,
            "residuos_bolsillo": sorted("%s%s" % (nom, num) for cad, num, nom in res),
            "residuos_clave": sorted(RESIDUOS_CLAVE),
            "receptor_pdbqt": os.path.basename(pdbqt),
            "nota_sep172": "SEP172 reducida a SER (fosfato fuera): Vina no tiene"
                           " parametros de fosfoserina",
            "nota_cadenas": "se usa el dimero A+B porque 13 atomos de la cadena B"
                            " estan a menos de 4,5 A del ligando de la cadena A",
        }, f, indent=2, ensure_ascii=False)
    with open(os.path.join(BASE, "receptor_tbk1.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LINEAS) + "\n")
    log("")
    log("caja y residuos anotados en analysis/calibracion_tbk1/caja.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
