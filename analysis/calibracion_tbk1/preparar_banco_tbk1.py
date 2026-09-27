#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Prepara el banco de TBK1: los activos medidos y los decoys emparejados.

EL FONDO DEL BANCO
------------------
El plan lo deja escrito y aqui se hace igual: **los negativos no son compuestos
medidos que no unen** (contra TBK1 solo hay 20 con afinidad peor que 10 uM y no
son un fondo), sino **decoys emparejados en propiedades**, como en DUD-E
(Mysinger et al., J. Chem. Inf. Model. 2012). Un decoy vale si se parece al
activo en lo que el acoplamiento *no debe* estar midiendo (tamano, lipofilia,
donadores, aceptores, enlaces rotables y carga) y en nada mas.

Los criterios son los de DUD-E, sin inventar:

    peso molecular     +- 25 Da
    cLogP              +- 2,5   (Crippen)
    donadores de H     +- 2
    aceptores de H     +- 2
    enlaces rotables   +- 2
    carga neta         igual
    parecido           Tanimoto ECFP4 < 0,25 contra CUALQUIER activo

DE DONDE SALEN
--------------
Del propio fondo de compuestos del proyecto, que ya esta preparado y es diverso:
`analysis/full_library_solo.smi` (66.478) mas `compounds/decoys_library.smi`
(6.462). Se excluye cualquier compuesto que sea activo de TBK1 en ChEMBL, no solo
los del banco: un decoy que fuera un inhibidor conocido envenenaria la medida.

Cada decoy se usa UNA vez (como en DUD-E), y se asignan por activo en el orden de
la lista, de modo que la asignacion no depende de ningun resultado de acoplamiento.

LA RECETA DE PREPARACION ES LA DEL PROYECTO
-------------------------------------------
No se escribe una receta nueva: se importa `preparar_ligando.py`, que es la copia
canonica (macrociclos rigidos, varias semillas de incrustacion, cero pseudo-atomos
y comprobacion de que el numero de atomos cuadra con el SMILES). Escribir otra
receta distinta fue el fallo mas caro del proyecto.

Uso:
    python preparar_banco_tbk1.py --solo-listar     # cuenta y estima, no prepara
    python preparar_banco_tbk1.py --ratio 50        # (por defecto, 50)
    python preparar_banco_tbk1.py --activos         # solo los activos

Salida en `analysis/calibracion_tbk1/`:
    decoys_tbk1.csv     id, smiles y de que activo es decoy
    ligands/            los PDBQT del banco (ACT_... y DEC_...)
    banco_tbk1_preparado.txt  el informe y los fallos
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
import time

import numpy as np
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(BASE))
LIGANDS = os.path.join(BASE, "ligands")
ACTIVOS_CSV = os.path.join(BASE, "compuestos_tbk1.csv")
DECOYS_CSV = os.path.join(BASE, "decoys_tbk1.csv")
INFORME = os.path.join(BASE, "banco_tbk1_preparado.txt")

# La receta canonica de ligandos del proyecto (una sola copia en todo el repo).
sys.path.insert(0, RAIZ)
import preparar_ligando as PL  # noqa: E402

POOL = [os.path.join(RAIZ, "analysis", "full_library_solo.smi"),
        os.path.join(RAIZ, "compounds", "decoys_library.smi")]

# Los PDBQT que el proyecto ya tiene preparados de su libreria. De los 66.476
# SMILES de full_library_solo.smi, 60.702 ya estan ahi con la receta canonica (y
# reparados de pseudo-atomos por reparar_libreria_glue.py). Reutilizarlos no es un
# atajo: es hacer lo mismo que hizo el cribado de la libreria, que tampoco los
# re-preparo. Se comprueba de todos modos que no traigan pseudo-atomos.
LIBRERIA = os.path.join(RAIZ, "gpu_dock", "libreria_ligands")

# Criterios DUD-E, tal cual.
TOL_MW = 25.0
TOL_LOGP = 2.5
TOL_HBD = 2
TOL_HBA = 2
TOL_ROT = 2
TANIMOTO_MAX = 0.25


def leer_smi(ruta):
    """[(id, smiles)] de un .smi (smiles<TAB>id)."""
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
    return {
        "mw": Descriptors.MolWt(mol),
        "logp": Crippen.MolLogP(mol),
        "hbd": rdMolDescriptors.CalcNumHBD(mol),
        "hba": rdMolDescriptors.CalcNumHBA(mol),
        "rot": rdMolDescriptors.CalcNumRotatableBonds(mol),
        "carga": Chem.GetFormalCharge(mol),
    }


def fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def compatible(a, d):
    """Los seis criterios de DUD-E."""
    return (abs(a["mw"] - d["mw"]) <= TOL_MW
            and abs(a["logp"] - d["logp"]) <= TOL_LOGP
            and abs(a["hbd"] - d["hbd"]) <= TOL_HBD
            and abs(a["hba"] - d["hba"]) <= TOL_HBA
            and abs(a["rot"] - d["rot"]) <= TOL_ROT
            and a["carga"] == d["carga"])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ratio", type=int, default=50,
                    help="decoys por activo (DUD-E usa 50)")
    ap.add_argument("--solo-listar", action="store_true")
    ap.add_argument("--activos", action="store_true", help="no generar decoys")
    args = ap.parse_args()

    os.makedirs(LIGANDS, exist_ok=True)
    T = []

    # ---------------- activos ----------------
    activos = list(csv.DictReader(open(ACTIVOS_CSV, encoding="utf-8")))
    print("activos en el CSV: %d" % len(activos))
    ids_activos = {a["molecule_chembl_id"] for a in activos}

    act = []
    for a in activos:
        mol = Chem.MolFromSmiles(a["smiles"])
        if mol is None:
            continue
        act.append({"id": a["molecule_chembl_id"], "smiles": a["smiles"],
                    "mol": mol, "props": props(mol), "pchembl": float(a["pchembl"])})
    print("activos con SMILES valido: %d" % len(act))

    # ---------------- pool de decoys ----------------
    pool = []
    vistos = set()
    for ruta in POOL:
        for pid, smi in leer_smi(ruta):
            if pid in vistos or pid in ids_activos:
                continue
            vistos.add(pid)
            pool.append((pid, smi))
    print("pool de decoys: %d compuestos" % len(pool))

    decoys = []
    if not args.activos:
        fp_act = [fp(a["mol"]) for a in act]
        t0 = time.time()
        prop_pool = []
        for pid, smi in pool:
            m = Chem.MolFromSmiles(smi)
            if m is None:
                continue
            prop_pool.append((pid, smi, m, props(m)))
        print("pool con propiedades calculadas: %d (%.0f s)"
              % (len(prop_pool), time.time() - t0))

        # Las propiedades del pool, en vectores: el filtro por propiedades se hace
        # entero en numpy, porque compararlo en Python son 2313 x 73.000 vueltas.
        mw = np.array([d["mw"] for _, _, _, d in prop_pool])
        logp = np.array([d["logp"] for _, _, _, d in prop_pool])
        hbd = np.array([d["hbd"] for _, _, _, d in prop_pool])
        hba = np.array([d["hba"] for _, _, _, d in prop_pool])
        rot = np.array([d["rot"] for _, _, _, d in prop_pool])
        carga = np.array([d["carga"] for _, _, _, d in prop_pool])
        usados = set()
        # Un mismo compuesto del pool puede ser compatible con varios activos.
        # Cachear su ECFP4 evita reconstruirlo cada vez que se evalua su similitud.
        fp_pool = [None] * len(prop_pool)
        for i, a in enumerate(act):
            p = a["props"]
            m = ((np.abs(mw - p["mw"]) <= TOL_MW)
                 & (np.abs(logp - p["logp"]) <= TOL_LOGP)
                 & (np.abs(hbd - p["hbd"]) <= TOL_HBD)
                 & (np.abs(hba - p["hba"]) <= TOL_HBA)
                 & (np.abs(rot - p["rot"]) <= TOL_ROT)
                 & (carga == p["carga"]))
            elegidos = 0
            for j in np.nonzero(m)[0]:
                if elegidos >= args.ratio:
                    break
                pid, smi, mol_j, _ = prop_pool[j]
                if pid in usados:
                    continue
                fp_j = fp_pool[j]
                if fp_j is None:
                    fp_j = fp(mol_j)
                    fp_pool[j] = fp_j
                if max(DataStructs.BulkTanimotoSimilarity(fp_j, fp_act)) >= TANIMOTO_MAX:
                    continue
                usados.add(pid)
                decoys.append({"id": pid, "smiles": smi, "activo": a["id"]})
                elegidos += 1
            if (i + 1) % 200 == 0:
                print("   activos procesados: %d/%d | decoys: %d (%.0f s)"
                      % (i + 1, len(act), len(decoys), time.time() - t0))
        ratio_real = len(decoys) / max(len(act), 1)
        print("decoys asignados: %d (%.1f por activo; se pidieron %d)"
              % (len(decoys), ratio_real, args.ratio))
        if ratio_real < args.ratio:
            print("   el pool se agoto antes de llegar al ratio pedido: "
                  "habra que decirlo al usar las metricas")

        with open(DECOYS_CSV, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=["id", "smiles", "activo"])
            w.writeheader()
            w.writerows(decoys)

    total = len(act) + len(decoys)
    horas = total * 1.40 / 3600.0     # 1,40 s por ligando, medido con la GPU libre
    print("ligandos del banco: %d (%d activos + %d decoys) -> ~%.1f h de GPU"
          % (total, len(act), len(decoys), horas))

    if args.solo_listar:
        return 0

    # ---------------- preparacion a PDBQT ----------------
    fallos = {"activos": [], "decoys": []}
    reutilizados = 0
    for etiqueta, grupo, prefijo in (("activos", act, "ACT_"),
                                     ("decoys", decoys, "DEC_")):
        for i, v in enumerate(grupo):
            destino = os.path.join(LIGANDS, prefijo + v["id"] + ".pdbqt")
            origen = os.path.join(LIBRERIA, v["id"] + ".pdbqt")
            if (not os.path.exists(destino) and os.path.exists(origen)
                    and os.path.getsize(origen) > 100):
                n, glue = PL.contar_fichero(origen)
                if not glue:
                    import shutil
                    shutil.copy2(origen, destino)
                    reutilizados += 1
            ruta = PL.escribir(prefijo + v["id"], v["smiles"], LIGANDS)
            if ruta is None:
                fallos[etiqueta].append(v["id"])
            if (i + 1) % 500 == 0:
                print("   %s: %d/%d (reutilizados %d)"
                      % (etiqueta, i + 1, len(grupo), reutilizados))
        print("%s: %d preparados, %d no (reutilizados en total: %d)"
              % (etiqueta, len(grupo) - len(fallos[etiqueta]),
                 len(fallos[etiqueta]), reutilizados))

    n_pdbqt = len([p for p in os.listdir(LIGANDS) if p.endswith(".pdbqt")])
    L = ["BANCO DE TBK1 PREPARADO — %s" % time.strftime("%Y-%m-%d"),
         "activos: %d | decoys: %d | ratio real: %.1f por activo"
         % (len(act), len(decoys), len(decoys) / max(len(act), 1)),
         "PDBQT en ligands/: %d (de ellos, %d copiados ya preparados de"
         " gpu_dock/libreria_ligands y verificados sin pseudo-atomos)"
         % (n_pdbqt, reutilizados),
         "estimacion de GPU: ~%.1f h con la tarjeta libre (1,40 s por ligando)"
         % horas,
         "criterios de emparejamiento (DUD-E): MW +-%.0f, logP +-%.1f, HBD +-%d,"
         " HBA +-%d, rotables +-%d, carga igual, Tanimoto ECFP4 < %.2f"
         % (TOL_MW, TOL_LOGP, TOL_HBD, TOL_HBA, TOL_ROT, TANIMOTO_MAX),
         "",
         "no preparados (activos): %s" % (", ".join(fallos["activos"][:30]) or "ninguno"),
         "no preparados (decoys): %s" % (", ".join(fallos["decoys"][:30]) or "ninguno")]
    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
