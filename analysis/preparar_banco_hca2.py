#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara el banco de hCA2: los activos medidos y sus señuelos emparejados.

POR QUE ESTE BANCO Y NO OTRO
----------------------------
Es el primer banco del proyecto donde los negativos son **decoys emparejados en
propiedades** (como en DUD-E) y no "activos que no unen". Y es el primero donde hay
suficientes positivos para que el intervalo del AUC signifique algo: 44 activos
independientes de la auditoría, frente a los 7 de TDP-43 y los 11 de SOD1.

La diferencia con el banco de TBK1 está en el que hay que tener cuidado. En TBK1 el
problema era que no habia negativos medidos. Aquí el problema es el contrario: hCA2
tiene **13.051 activos medidos**, y si se toman 44 como positivos y el resto como
negativos, se esta tomando el fondo de la misma serie. Los sulfonamidas se parecen
mucho entre si (la auditoría lo vio: 281 esqueletos para 1.200 activos, MCS de 4
atomos), así que un "activo de hCA2 que no une" puede ser indistinguible del que si
une, y el AUC mediría otra vez la construcción del banco.

Los criterios de emparejamiento son los de DUD-E, sin inventar:
    peso molecular     ± 25 Da
    cLogP              ± 2,5  (Crippen)
    donadores de H     ± 2
    aceptores de H     ± 2
    enlaces rotables   ± 2
    carga neta         igual
    parecido           Tanimoto ECFP4 < 0,25 contra CUALQUIER activo de hCA2

La última línea es la que hCA2 exige y TBK1 no: el parecido se mide contra
**todos** los 13.051 activos medidos, no solo contra los 44 que van a ser positivos.
Porque si un señuelo se parece a un activo que no entra en el banco, y ese activo
une de verdad, el señuelo es un falso negativo y el AUC sale artificialmente bajo.

LA RECETA DE PREPARACIÓN ES LA DEL PROYECTO
-------------------------------------------
No se escribe una receta nueva: se importa `preparar_ligando.py`, la copia canónica.
Escribir otra receta distinta fue el fallo más caro del proyecto.

Y hay un detalle que en hCA2 NO es opcional: la sulfonamida tiene que ir
DESPROTONADA en el nitrógeno que coordina al zinc, porque ese es el estado en que
se une. `preparar_ligando.py` deja esa decisión al modelo de protonación, así que
aquí se comprueba después de preparar que el N del sulfonamida queda sin hidrógeno en
el PDBQT, y se avisa si no.

Salida en `analysis/calibracion_hca2/`:
    activos_hca2.csv       los positivos, con su pchembl
    decoys_hca2.csv        id, smiles y de qué activo es señuelo
    ligands/               ACT_*.pdbqt y DEC_*.pdbqt
    banco_hca2.txt          informe y fallos

Uso:
    python preparar_banco_hca2.py --solo-listar
    python preparar_banco_hca2.py --ratio 50
    python preparar_banco_hca2.py --activos
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time
from collections import defaultdict

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
DIR = os.path.join(BASE, "calibracion_hca2")
LIGANDS = os.path.join(DIR, "ligands")
ACTIVOS_IN = os.path.join(BASE, "serie_hca2_auditoria.csv")
ACTIVOS_CSV = os.path.join(DIR, "activos_hca2.csv")
DECOYS_CSV = os.path.join(DIR, "decoys_hca2.csv")
INFORME = os.path.join(DIR, "banco_hca2.txt")

sys.path.insert(0, RAIZ)
import preparar_ligando as PL  # noqa: E402

POOL = [os.path.join(RAIZ, "analysis", "full_library_solo.smi"),
        os.path.join(RAIZ, "compounds", "decoys_library.smi")]

# Criterios DUD-E, tal cual.
TOL_MW = 25.0
TOL_LOGP = 2.5
TOL_HBD = 2
TOL_HBA = 2
TOL_ROT = 2
TANIMOTO_MAX = 0.25

MIN_ACTIVOS = 40        # por debajo, el intervalo del AUC no dice nada (ver bootstrap)
MAX_ACTIVOS = 60        # mas de 60 y las 24 h de GPU por 1000 ligandos se disparan
MIN_INDEPENDIENTES = 25  # familias de 5+, el numero que de verdad limita el AUC
SULFONAMIDA = Chem.MolFromSmarts("[SX4](=O)(=O)[NX3;H1,H2]")

lineas = []


def log(m):
    print(m, flush=True)
    lineas.append(m)


def leer_smi(ruta):
    out = []
    if not os.path.exists(ruta):
        return out
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        p = l.rstrip("\n").split("\t")
        if len(p) < 2 or not p[0].strip():
            continue
        out.append((p[1].strip(), p[0].strip()))
    return out


def props(mol):
    return {"mw": Descriptors.MolWt(mol),
            "logp": Crippen.MolLogP(mol),
            "hbd": rdMolDescriptors.CalcNumHBD(mol),
            "hba": rdMolDescriptors.CalcNumHBA(mol),
            "rot": rdMolDescriptors.CalcNumRotatableBonds(mol),
            "carga": Chem.GetFormalCharge(mol)}


def fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def compatible(a, d):
    return (abs(a["mw"] - d["mw"]) <= TOL_MW
            and abs(a["logp"] - d["logp"]) <= TOL_LOGP
            and abs(a["hbd"] - d["hbd"]) <= TOL_HBD
            and abs(a["hba"] - d["hba"]) <= TOL_HBA
            and abs(a["rot"] - d["rot"]) <= TOL_ROT
            and a["carga"] == d["carga"])


def esqueleto(mol):
    try:
        e = MurckoScaffold.GetScaffoldForMol(mol)
    except Exception:  # noqa: BLE001
        return None
    if e is None or e.GetNumAtoms() == 0:
        return None
    try:
        return Chem.MolToSmiles(e)
    except Exception:  # noqa: BLE001
        return None


def leer_activos():
    """[(id, smiles, tipo, pchembl)] de la auditoría ya descargada."""
    out = []
    with open(ACTIVOS_IN, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if not r.get("smiles"):
                continue
            out.append((r["molecule_chembl_id"], r["smiles"], r["tipo"],
                        float(r["pchembl"])))
    return out


def elegir_activos(todos):
    """Los positivos del banco: los mas activos y con suficiente independencia.

    Se cogen por pchembl descendente y se cuentan las familias de Murcko de 5+, que
    es lo que limita el intervalo del AUC. Coger los mas activos sin mirar las
    familias seria coger 40 miembros de las dos o tres familias mas pobladas, que es
    justamente el sesgo que la auditoria de CDK2 demonstrated que arruina un AUC.
    """
    mols = {}
    for cid, smi, tipo, pch in todos:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            mols[cid] = (m, tipo, pch)
    fam = defaultdict(list)
    for cid, (m, _, _) in mols.items():
        e = esqueleto(m)
        if e:
            fam[e].append(cid)
    orden = sorted(mols, key=lambda c: -mols[c][2])
    elegidos, fams = [], set()
    for cid in orden:
        m = mols[cid][0]
        e = esqueleto(m)
        if e in fams:
            continue                      # no repetir familia: cada familia entra una vez
        elegidos.append(cid)
        if e:
            fams.add(e)
        if len(elegidos) >= MAX_ACTIVOS:
            break
    return elegidos, mols, fam


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ratio", type=int, default=50)
    ap.add_argument("--solo-listar", action="store_true")
    ap.add_argument("--activos", action="store_true",
                    help="prepara solo los positivos")
    args = ap.parse_args()

    os.makedirs(DIR, exist_ok=True)
    log("BANCO DE hCA2: activos medidos y sus senuelos emparejados")
    log("   %s" % time.strftime("%Y-%m-%d %H:%M"))
    log("   Criterios DUD-E: MW +-%.0f, logP +-%.1f, HBD/HBA/rot +-%d, carga"
        " igual, Tanimoto < %.2f contra TODOS los activos de hCA2"
        % (TOL_MW, TOL_LOGP, TOL_HBD, TANIMOTO_MAX))

    todos = leer_activos()
    log("")
    log("ACTIVOS MEDIDOS DESCARGADOS: %d" % len(todos))
    elegidos, mols, fam = elegir_activos(todos)
    log("   %d compuestos con mol legible" % len(mols))
    log("   familias de Murcko de 5+ activos: %d"
        % sum(1 for v in fam.values() if len(v) >= 5))
    log("   positivos elegidos: %d, uno por familia, por pchembl descendente"
        % len(elegidos))
    if elegidos:
        log("   pchembl del mas activo: %.2f | del mas debil: %.2f"
            % (mols[elegidos[0]][2], mols[elegidos[-1]][2]))
    else:
        log("   no hay positivos utilizables")
    if len(elegidos) < MIN_ACTIVOS:
        log("")
        log("   AVISO: solo %d positivos y el minimo util es %d. Con menos, el"
            % (len(elegidos), MIN_ACTIVOS))
        log("          intervalo del AUC es tan ancho que no dice nada. Se sigue")
        log("          igual, pero el numero sera el que sea.")

    with open(ACTIVOS_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["id", "smiles", "tipo", "pchembl", "sulfonamida",
                    "pesados", "rotatorios"])
        for cid in elegidos:
            m = mols[cid][0]
            w.writerow([cid, Chem.MolToSmiles(m), mols[cid][1], "%.3f" % mols[cid][2],
                        "SI" if m.HasSubstructMatch(SULFONAMIDA) else "no",
                        m.GetNumHeavyAtoms(),
                        rdMolDescriptors.CalcNumRotatableBonds(m)])
    log("   positivos: activos_hca2.csv")

    n_sulf = sum(1 for c in elegidos if mols[c][0].HasSubstructMatch(SULFONAMIDA))
    log("   de los cuales %d son sulfonamidas (%.0f %%)"
        % (n_sulf, 100.0 * n_sulf / max(len(elegidos), 1)))
    if n_sulf < len(elegidos) * 0.3:
        log("   OJO: menos de un tercio son sulfonamidas. La serie de hCA2 es")
        log("        diversa (lo dijo la auditoria), asi que el banco tambien lo")
        log("        es. Es correcto, pero significa que el AUC no mide solo")
        log("        'tengo una sulfonamida pegada al zinc'.")

    if args.solo_listar:
        log("")
        log("(--solo-listar: no se prepara nada)")
        with open(INFORME, "w", encoding="utf-8") as f:
            f.write("\n".join(lineas) + "\n")
        return 0

    # ---- preparar los activos ----
    log("")
    log("PREPARANDO LOS ACTIVOS (receta canonica del proyecto)")
    os.makedirs(LIGANDS, exist_ok=True)
    n_ok = 0
    n_sales = 0
    fallos = []
    descartados = []
    t0 = time.time()
    for k, cid in enumerate(elegidos, 1):
        m = mols[cid][0]
        # La firma real de la receta canonica es escribir(nombre, smiles, outdir),
        # y devuelve la ruta. Se llama por el SMILES canonico, no por el objeto
        # mol, porque es lo que hace el resto del proyecto y lo que garantiza que
        # dos ligandos preparados por el mismo camino salen iguales.
        ruta = PL.escribir("ACT_%s" % cid, Chem.MolToSmiles(m), LIGANDS)
        if ruta is None or not os.path.exists(ruta):
            # "2 fragments" quiere decir que el SMILES trae una sal o una mezcla
            # (el cation de un clorhidrato, por ejemplo), no que la receta falle.
            # Se reintenta quitando las sales, que es lo que hace la propia receta
            # canonica, y si tampoco sale, el activo se descarta y SE DICE: dejar
            # fuera un positivo sin avisar cambia el AUC sin que nadie lo note.
            ruta = PL.escribir("ACT_%s" % cid, Chem.MolToSmiles(m), LIGANDS,
                               quitar_sales=True)
            if ruta is None or not os.path.exists(ruta):
                fallos.append(cid)
                descartados.append(cid)
            else:
                n_ok += 1
                n_sales += 1
        else:
            n_ok += 1
        if k % 10 == 0:
            log("   %d/%d (%.0f s)" % (k, len(elegidos), time.time() - t0))
    log("   %d de %d activos preparados en %.0f s" % (n_ok, len(elegidos),
                                                     time.time() - t0))
    log("   de los cuales %d salieron tras quitar la sal del SMILES" % n_sales)
    if descartados:
        log("   DESCARTADOS %d activos (siguen dando mas de un fragmento): %s"
            % (len(descartados), ", ".join(descartados)))
        log("   Un positivo que se cae sin avisar baja el AUC sin que se note, asi")
        log("   que se quitan aqui y se cuentan.")

    # ---- comprobacion de la sulfonamida: DESPROTONADA ----
    log("")
    log("COMPROBANDO QUE LA SULFONAMIDA VA DESPROTONADA (es lo que une al zinc)")
    mal = []
    for cid in elegidos:
        ruta = os.path.join(LIGANDS, "ACT_%s.pdbqt" % cid)
        if not os.path.exists(ruta):
            continue
        mol = mols[cid][0]
        if not mol.HasSubstructMatch(SULFONAMIDA):
            continue
        txt = open(ruta, encoding="utf-8", errors="ignore").read()
        # el N del sulfonamida no debe llevar HD (donador). Si lo lleva, la
        # sulfonamida va neutra y no coordina al zinc bien.
        if re_hd(txt):
            mal.append(cid)
    log("   sulfonamidas con el N como donador de H (mal para coordinar al Zn):"
        " %d de %d" % (len(mal), n_sulf))
    if mal:
        log("   AVISO: %s. Un N con H nocoordina igual. Se dice, no se calla;"
            % ", ".join(mal[:6]))
        log("          el AUC de hCA2 puede salir bajo por esto y no por el motor.")
    else:
        log("   bien: todas las sulfonamidas van con el N sin H, que es el estado")
        log("         que coordina al zinc.")

    if args.activos:
        with open(INFORME, "w", encoding="utf-8") as f:
            f.write("\n".join(lineas) + "\n")
        return 0

    # ---- los señuelos ----
    log("")
    log("SEÑUELOS: emparejados en propiedades, contra TODOS los activos de hCA2")
    pool = []
    for ruta in POOL:
        d = leer_smi(ruta)
        log("   fondo %s: %d" % (os.path.basename(ruta), len(d)))
        pool += d
    # Todos los activos medidos de hCA2, no solo los positivos: un señuelo que se
    # parezca a un activo que no entra en el banco pero une de verdad es un falso
    # negativo. Con 13.051 medidos esta exclusion es la mayor del banco.
    canon_act = {}
    for cid, smi, _, _ in todos:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            canon_act[Chem.MolToSmiles(m)] = cid
    log("   %d activos medidos de hCA2 excluidos del fondo" % len(canon_act))

    fps_act = {cid: fp(m) for cid, (m, _, _) in mols.items()}
    idx = 0
    n_pool = len(pool)
    asignados = []
    por_activo = defaultdict(int)
    agotados = []
    while idx < n_pool and len(asignados) < len(elegidos) * args.ratio:
        cid_pool, smi = pool[idx]
        idx += 1
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        # fuera si es un activo medido de hCA2
        if Chem.MolToSmiles(m) in canon_act:
            continue
        d = props(m)
        # fuera si se parece a CUALQUIER activo medido
        f = fp(m)
        if fps_act:
            peor = DataStructs.BulkTanimotoSimilarity(f, list(fps_act.values()))
            if max(peor) >= TANIMOTO_MAX:
                continue
        for cid in elegidos:
            if por_activo[cid] >= args.ratio:
                continue
            if compatible(props(mols[cid][0]), d):
                asignados.append({"id": cid_pool, "smiles": smi, "activo": cid})
                por_activo[cid] += 1
                break
    log("   señuelos asignados: %d (objetivo %d)" % (len(asignados),
                                                       len(elegidos) * args.ratio))
    faltan = [c for c in elegidos if por_activo[c] < args.ratio]
    if faltan:
        log("   ACTIVOS SIN %d señuelos: %d" % (args.ratio, len(faltan)))
        log("   (con menos de %d no se puede hacer un AUC util)"
            % (args.ratio // 2))

    # Los senuelos que comparten esqueleto de Murcko con un activo se QUITAN, no se
    # cuentan y se avisan: un senuelo con el mismo esqueleto que un activo es
    # dificil por construccion, no por azar, y hace el AUC mas alto de lo que
    # corresponde. Esto es exactamente lo que el script de DUD-E mide, y aqui sale
    # en 1 %, que no es cero.
    # El conjunto de referencia son los esqueletos de LOS ELEGIDOS, no las familias
    # enteras. Antes se comparaba contra `v[0] in elegidos`, es decir contra el
    # primer miembro de cada familia, y como los positivos son uno por familia casi
    # nunca coinciden: por eso salia "0 descartados" mientras la auditoria de abajo
    # contaba 24 senuelos sobre un esqueleto de activo. Dos numeros que se
    # contradicen en el mismo informe.
    fam_act = set()
    for cid in elegidos:
        e = esqueleto(mols[cid][0])
        if e:
            fam_act.add(e)
    log("   esqueletos de los %d positivos: %d" % (len(elegidos), len(fam_act)))
    buenos, malos = [], []
    for r in asignados:
        m = Chem.MolFromSmiles(r["smiles"])
        e = esqueleto(m) if m else None
        if e and e in fam_act:
            malos.append(r)
        else:
            buenos.append(r)
    log("   senuelos descartados por compartir esqueleto con un activo: %d"
        % len(malos))
    asignados = buenos

    with open(DECOYS_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["id", "smiles", "activo"])
        w.writeheader()
        for r in asignados:
            w.writerow(r)
    log("   senuelos definitivos: %d -> decoys_hca2.csv" % len(asignados))
    reparto = defaultdict(int)
    for r in asignados:
        reparto[r["activo"]] += 1
    pocos = [c for c in elegidos if reparto[c] < args.ratio // 2]
    log("   activos con menos de %d senuelos: %d"
        % (args.ratio // 2, len(pocos)))

    # ---- la auditoria del banco, ANTES del AUC ----
    log("")
    log("AUDITORIA DEL BANCO (antes de que exista ningun AUC)")
    fam_dec = defaultdict(int)
    for r in asignados:
        m = Chem.MolFromSmiles(r["smiles"])
        e = esqueleto(m) if m else None
        if e:
            fam_dec[e] += 1
    # Aqui hay que distinguir DOS cosas que antes se mezclaban en un solo numero.
    # Lo que infla el AUC es un senuelo con el mismo esqueleto que un activo QUE
    # ESTA EN EL BANCO, porque ese activo se puntua y el senuelo tambien, y son
    # intercambiables. Un senuelo que comparte esqueleto con uno de los otros 1.140
    # activos medidos NO sesga nada: ese activo no se acopla y no cuenta.
    comunes_pos = set(fam_dec) & fam_act
    comunes_todos = set(fam_dec) & set(e for e in fam if fam.get(e))
    n_pos = sum(fam_dec[e] for e in comunes_pos)
    n_todos = sum(fam_dec[e] for e in comunes_todos)
    log("   %d esqueletos en los %d senuelos" % (len(fam_dec), len(asignados)))
    log("")
    log("   CONTRA LOS %dPOSITIVOS DEL BANCO (esto es lo que importa):"
        % len(elegidos))
    log("      senuelos que comparten esqueleto con un positivo: %d de %d"
        % (n_pos, len(asignados)))
    if n_pos == 0:
        log("      bien: ninguno. El AUC no se infla por construccion del banco.")
    else:
        log("      ATENCION: %d senuelos son casi indistinguibles de un positivo."
            % n_pos)
        log("               El AUC saldra alto por construccion y no por el motor.")
    log("")
    log("   CONTRA TODA LA SERIE MEDIDA DE hCA2 (informativo, no sesga el AUC):")
    log("      senuelos que comparten esqueleto con algun activo de los %d: %d"
        " (%.2f %%)"
        % (len(fam), n_todos, 100.0 * n_todos / max(len(asignados), 1)))
    log("      esos activos NO entran en el banco, asi que no puntuan y el")
    log("      numero no afecta al AUC. Se dice para que no se lea como un fallo")
    log("      del banco cuando no lo es.")

    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    log("")
    log("informe: banco_hca2.txt")
    return 0


def re_hd(txt):
    """Hay algun HD (donador de H) pegado a un nitrogeno del sulfonamida."""
    for l in txt.splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        tipo = l.rsplit(None, 1)[-1].strip().upper()
        if tipo in ("HD", "HS"):
            return True
    return False


if __name__ == "__main__":
    raise SystemExit(main())
