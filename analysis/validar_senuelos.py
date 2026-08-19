#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""validar_senuelos.py — Protocolo de validación con señuelos (decoy validation).

Mide si el ranking de docking discrimina ligandos activos conocidos (controles
positivos) de señuelos inactivos con propiedades similares (property-matched
decoys). Métricas: ROC-AUC y factor de enriquecimiento EF1%/EF5%.
Un AUC > 0.7-0.8 apoya la utilidad del protocolo de cribado.

Flujo:
  1. Leer activos conocidos (CSV: name, smiles, target).
  2. Generar decoys desde la librería fuente: MW/logP/rotables similares y
     similitud Tanimoto (Morgan r=2) < 0.35 frente a cada activo.
  3. Convertir activos + decoys a PDBQT (RDKit ETKDGv3 + MMFF + meeko,
     el mismo pipeline que la conversión ZINC 2M).
  4. Dockear con Vina-GPU usando la misma caja del cribado.
  5. Calcular ROC-AUC y EF1%/EF5% por target.

Ejemplo:
  python validar_senuelos.py --activos activos_conocidos.csv \
      --libreria ../compounds/full_library.smi --target SOD1 \
      --receptor ../gpu_dock/SOD1.pdbqt --centro 27.9,111.8,64.4 \
      --decoys-por-activo 30 --salida validacion_SOD1.csv
"""
import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys
import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, DataStructs
from rdkit import RDLogger
from meeko import MoleculePreparation, PDBQTWriterLegacy

RDLogger.DisableLog("rdApp.*")

VINA = r"C:\Users\Fredy\masive-als\gpu_dock\Vina-GPU-2.1-win.exe"
VINA_CPU = r"C:\Users\Fredy\masive-als\tools\vina.exe"
GBASE = r"C:\Users\Fredy\masive-als\gpu_dock"  # opencl_binary_path
SIZE = 25
NUM_MODES = 3
THREAD = 8000


LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "validacion_senuelos.log")


def log(m):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), m)
    print(line, flush=True)
    try:
        with open(LOG_FILE, "a", encoding="utf-8") as f:
            f.write(line + "\n")
    except Exception:
        pass


# ---------------- lectura ----------------
def leer_activos(ruta):
    act = []
    with open(ruta, encoding="utf-8") as f:
        rd = csv.DictReader(f)
        for r in rd:
            name = (r.get("name") or "").strip()
            smi = (r.get("smiles") or "").strip()
            tgt = (r.get("target") or "").strip()
            if name and smi and Chem.MolFromSmiles(smi):
                act.append({"name": name, "smiles": smi, "target": tgt})
    return act


def leer_libreria(ruta, limite=200000):
    """Lee SMILES de un .smi (smiles<TAB>name) o .csv (columna smiles)."""
    out = []
    if ruta.endswith(".smi"):
        with open(ruta, encoding="utf-8", errors="replace") as f:
            for line in f:
                p = line.rstrip("\n").split("\t")
                if p and p[0].strip() and Chem.MolFromSmiles(p[0].strip()):
                    name = p[1].strip() if len(p) > 1 else p[0].strip()
                    out.append((name, p[0].strip()))
                if len(out) >= limite:
                    break
    else:
        with open(ruta, encoding="utf-8", errors="replace") as f:
            rd = csv.DictReader(f)
            for r in rd:
                s = (r.get("smiles") or "").strip()
                n = (r.get("name") or r.get("chembl_id") or s).strip()
                if s and Chem.MolFromSmiles(s):
                    out.append((n, s))
                if len(out) >= limite:
                    break
    return out


# ---------------- decoys ----------------
def _fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)


def generar_decoys(activo_smi, libreria, n=30):
    a = Chem.MolFromSmiles(activo_smi)
    if a is None:
        return []
    mw_a = Descriptors.MolWt(a)
    logp_a = Descriptors.MolLogP(a)
    rot_a = Descriptors.NumRotatableBonds(a)
    fa = _fp(a)
    cands = []
    for name, s in libreria:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        if abs(Descriptors.MolWt(m) - mw_a) > 30:
            continue
        if abs(Descriptors.MolLogP(m) - logp_a) > 1.0:
            continue
        if abs(Descriptors.NumRotatableBonds(m) - rot_a) > 2:
            continue
        if DataStructs.TanimotoSimilarity(fa, _fp(m)) > 0.35:
            continue
        cands.append((name, s))
    cands.sort(key=lambda x: abs(Descriptors.MolWt(Chem.MolFromSmiles(x[1])) - mw_a))
    return cands[:n]


# ---------------- conversión SMILES -> PDBQT ----------------
def convertir_pdbqt(name, smi, outdir):
    outp = os.path.join(outdir, name + ".pdbqt")
    if os.path.exists(outp) and os.path.getsize(outp) > 100:
        return outp
    try:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            return None
        mol = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            return None
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
        except Exception:
            pass
        prep = MoleculePreparation()
        setups = prep.prepare(mol)
        if not setups:
            return None
        pdbqt, ok, _ = PDBQTWriterLegacy.write_string(setups[0])
        if not ok:
            return None
        with open(outp, "w", encoding="utf-8") as f:
            f.write(pdbqt)
        return outp if os.path.getsize(outp) > 100 else None
    except Exception:
        return None


# ---------------- docking ----------------
def parse_affinity(pdbqt_path):
    try:
        with open(pdbqt_path, encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("REMARK VINA RESULT:"):
                    m = re.search(r"([-+]?\d+\.\d+)", line)
                    if m:
                        return float(m.group(1))
    except Exception:
        pass
    return None


def dockear(receptor, ligand_dir, out_dir, cx, cy, cz):
    tmp = os.path.join(out_dir, "_tmp")
    os.makedirs(os.path.join(tmp, "ligands"), exist_ok=True)
    os.makedirs(os.path.join(tmp, "out"), exist_ok=True)
    for p in glob.glob(os.path.join(ligand_dir, "*.pdbqt")):
        shutil.copy(p, os.path.join(tmp, "ligands", os.path.basename(p)))
    cfg = os.path.join(tmp, "config.txt")
    with open(cfg, "w", encoding="utf-8") as f:
        f.write("receptor = %s\n" % receptor.replace("\\", "/"))
        f.write("ligand_directory = %s\n" % os.path.join(tmp, "ligands").replace("\\", "/"))
        f.write("output_directory = %s\n" % os.path.join(tmp, "out").replace("\\", "/"))
        f.write("opencl_binary_path = %s\n" % GBASE.replace("\\", "/"))
        f.write("center_x = %s\ncenter_y = %s\ncenter_z = %s\n" % (cx, cy, cz))
        f.write("size_x = %d\nsize_y = %d\nsize_z = %d\n" % (SIZE, SIZE, SIZE))
        f.write("num_modes = %d\nthread = %d\n" % (NUM_MODES, THREAD))
    subprocess.run([VINA, "--config", cfg], cwd=GBASE, capture_output=True,
                   timeout=4 * 3600)
    res = {}
    for p in glob.glob(os.path.join(tmp, "out", "*_out.pdbqt")):
        lig = os.path.basename(p).replace("_out.pdbqt", "")
        aff = parse_affinity(p)
        if aff is not None:
            res[lig] = aff
    shutil.rmtree(tmp, ignore_errors=True)
    return res


# ---------------- docking CPU (AutoDock Vina 1.2.3) ----------------
def _dock_one_cpu(args):
    vina_exe, receptor, lig_pdbqt, out_pdbqt, cx, cy, cz, exhaust = args
    cmd = [vina_exe, "--receptor", receptor, "--ligand", lig_pdbqt,
           "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", str(SIZE), "--size_y", str(SIZE), "--size_z", str(SIZE),
           "--exhaustiveness", str(exhaust), "--num_modes", "3",
           "--out", out_pdbqt, "--cpu", "1"]
    subprocess.run(cmd, capture_output=True, timeout=900)
    return out_pdbqt


def dockear_cpu(receptor, ligand_dir, out_dir, cx, cy, cz, workers=6, exhaustividad=4):
    from concurrent.futures import ProcessPoolExecutor
    ligs = glob.glob(os.path.join(ligand_dir, "*.pdbqt"))
    os.makedirs(out_dir, exist_ok=True)
    tasks = []
    for p in ligs:
        name = os.path.basename(p).replace(".pdbqt", "")
        outp = os.path.join(out_dir, name + "_out.pdbqt")
        tasks.append((VINA_CPU, receptor, p, outp, cx, cy, cz, exhaustividad))
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for _ in ex.map(_dock_one_cpu, tasks):
            pass
    res = {}
    for p in glob.glob(os.path.join(out_dir, "*_out.pdbqt")):
        lig = os.path.basename(p).replace("_out.pdbqt", "")
        aff = parse_affinity(p)
        if aff is not None:
            res[lig] = aff
    return res


# ---------------- métricas ----------------
def roc_auc(act, dec):
    if not act or not dec:
        return None
    act = np.array(act, dtype=float)
    dec = np.array(dec, dtype=float)
    auc = 0.0
    for a in act:
        auc += float(np.sum(a <= dec))
    return auc / (len(act) * len(dec))


def enrichment_factor(act, dec, pct=1.0):
    todo = [(s, 1) for s in act] + [(s, 0) for s in dec]
    todo.sort(key=lambda x: x[0])
    k = max(1, int(len(todo) * pct / 100.0))
    n_act_top = sum(1 for _, lab in todo[:k] if lab == 1)
    return (n_act_top / max(1, len(act))) / (pct / 100.0)


# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--activos", required=True)
    ap.add_argument("--libreria", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--centro", required=True, help="x,y,z")
    ap.add_argument("--decoys-por-activo", type=int, default=30)
    ap.add_argument("--cpu", action="store_true", help="usar Vina CPU (sin GPU)")
    ap.add_argument("--exhaustividad", type=int, default=4,
                    help="exhaustividad de Vina (default 4)")
    ap.add_argument("--salida", default="validacion_senuelos.csv")
    args = ap.parse_args()

    cx, cy, cz = [float(x) for x in args.centro.split(",")]
    activos = [a for a in leer_activos(args.activos) if a["target"] == args.target]
    log("Activos para %s: %d" % (args.target, len(activos)))
    if not activos:
        print("ERROR: no hay activos para el target %s" % args.target, file=sys.stderr)
        sys.exit(1)

    libreria = leer_libreria(args.libreria)
    log("Librería fuente: %d compuestos" % len(libreria))

    work = os.path.join(os.path.dirname(os.path.abspath(args.salida)),
                        "_validacion_" + args.target)
    ligdir = os.path.join(work, "ligands")
    outdir = os.path.join(work, "out")
    os.makedirs(ligdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)

    usados = set()
    for a in activos:
        convertir_pdbqt("ACT_" + a["name"], a["smiles"], ligdir)
        dec = generar_decoys(a["smiles"], libreria, args.decoys_por_activo)
        n_dec = 0
        for dn, ds in dec:
            if dn in usados:
                continue
            usados.add(dn)
            convertir_pdbqt("DEC_" + dn, ds, ligdir)
            n_dec += 1
        log("  %s: %d decoys generados" % (a["name"], n_dec))

    n_lig = len(glob.glob(os.path.join(ligdir, "*.pdbqt")))
    log("Ligandos preparados (activos+decoys): %d" % n_lig)

    if args.cpu:
        log("Dockeando con Vina CPU contra %s..." % args.target)
        affs = dockear_cpu(args.receptor, ligdir, outdir, cx, cy, cz,
                           exhaustividad=args.exhaustividad)
    else:
        log("Dockeando con Vina-GPU contra %s..." % args.target)
        affs = dockear(args.receptor, ligdir, outdir, cx, cy, cz)
    log("Docking completado: %d resultados" % len(affs))

    act_scores = [v for k, v in affs.items() if k.startswith("ACT_")]
    dec_scores = [v for k, v in affs.items() if k.startswith("DEC_")]

    auc = roc_auc(act_scores, dec_scores)
    ef1 = enrichment_factor(act_scores, dec_scores, 1.0)
    ef5 = enrichment_factor(act_scores, dec_scores, 5.0)

    log("RESULTADO %s:" % args.target)
    log("  activos: %d | decoys: %d" % (len(act_scores), len(dec_scores)))
    log("  media energía activos: %.2f | decoys: %.2f" % (
        np.mean(act_scores) if act_scores else 0.0,
        np.mean(dec_scores) if dec_scores else 0.0))
    log("  ROC-AUC: %.3f" % (auc if auc is not None else -1))
    log("  EF1%%: %.2f | EF5%%: %.2f" % (ef1, ef5))

    with open(args.salida, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["ligand", "rol", "affinity"])
        for k, v in sorted(affs.items(), key=lambda kv: kv[1]):
            rol = "activo" if k.startswith("ACT_") else "decoy"
            w.writerow([k, rol, v])
    log("CSV de detalle: %s" % args.salida)


if __name__ == "__main__":
    main()
