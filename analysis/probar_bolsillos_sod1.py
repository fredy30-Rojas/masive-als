#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""probar_bolsillos_sod1.py — Comparar bolsillos de docking para SOD1.

Reutiliza los ligandos ya preparados (activos + decoys) en
`_validacion_SOD1/ligands/` y los dockea (Vina CPU, sin tocar la GPU) contra
varios centros de caja, calculando ROC-AUC y EF1%/EF5% para cada uno.

Bolsa de bolsillos candidatos (coordenadas del receptor SOD1.pdbqt, 1hl5):
  1. trp32     : bolsillo de Trp32 (Wright 2013, 4A7S/4A7T). Residuos Trp32,
                 Glu21, Gln22, Pro28, Lys30, Ser98, Glu100. Centroide ~(46,80,73).
  2. metal     : sitio activo Cu/Zn (canal catalítico). Cu A154 ~(40.6,99.0,79.0),
                 Zn A155 ~(46.6,100.4,77.7). Centroide ~(43.6,99.7,78.3).
  3. dimer     : interfaz del dímero (cavidad hidrofóbica Val7-Gly147-Val148).
                 Centroide ~(35.7,87.6,84.4).
  4. actual    : caja actual del cribado (interfaz cristalina E/N). (27.9,111.8,64.4)

Uso:
  python probar_bolsillos_sod1.py --ligandos _validacion_SOD1/ligands \
      --receptor ../gpu_dock/SOD1.pdbqt --exhaustividad 8 --cpu-workers 6
"""
import argparse
import csv
import glob
import os
import re
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

VINA_CPU = r"C:\Users\Fredy\masive-als\tools\vina.exe"
LOG_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "probar_bolsillos.log")


def log(m):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), m)
    print(line, flush=True)
    with open(LOG_FILE, "a", encoding="utf-8") as f:
        f.write(line + "\n")


BOLSILLOS = {
    "trp32":  (46.5, 80.0, 73.3, 22),
    "metal":  (43.6, 99.7, 78.3, 22),
    "dimer":  (35.7, 87.6, 84.4, 22),
    "actual": (27.9, 111.8, 64.4, 25),
}


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


def _dock_one(args):
    receptor, lig, out, cx, cy, cz, size, exh = args
    cmd = [VINA_CPU, "--receptor", receptor, "--ligand", lig,
           "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
           "--size_x", str(size), "--size_y", str(size), "--size_z", str(size),
           "--exhaustiveness", str(exh), "--num_modes", "3",
           "--out", out, "--cpu", "1"]
    subprocess.run(cmd, capture_output=True, timeout=1200)
    return out


def dock_pocket(receptor, ligand_dir, out_dir, cx, cy, cz, size, exh, workers):
    os.makedirs(out_dir, exist_ok=True)
    ligs = sorted(glob.glob(os.path.join(ligand_dir, "*.pdbqt")))
    tasks = []
    for p in ligs:
        name = os.path.basename(p).replace(".pdbqt", "")
        outp = os.path.join(out_dir, name + "_out.pdbqt")
        if os.path.exists(outp) and os.path.getsize(outp) > 100:
            continue  # ya hecho, resumible
        tasks.append((receptor, p, outp, cx, cy, cz, size, exh))
    log("  %d ligandos a dockear (ya hechos: %d)" %
        (len(tasks), len(ligs) - len(tasks)))
    if tasks:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for _ in ex.map(_dock_one, tasks):
                pass
    res = {}
    for p in glob.glob(os.path.join(out_dir, "*_out.pdbqt")):
        lig = os.path.basename(p).replace("_out.pdbqt", "")
        aff = parse_affinity(p)
        if aff is not None:
            res[lig] = aff
    return res


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligandos", required=True)
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--exhaustividad", type=int, default=8)
    ap.add_argument("--cpu-workers", type=int, default=6)
    ap.add_argument("--solo", default="", help="nombre de bolsillo (coma sep)")
    args = ap.parse_args()

    receptor = os.path.abspath(args.receptor)
    base = os.path.dirname(os.path.abspath(args.ligandos))

    nombres = [b for b in BOLSILLOS] if not args.solo else \
        [b.strip() for b in args.solo.split(",") if b.strip()]

    resumen = []
    for nombre in nombres:
        cx, cy, cz, size = BOLSILLOS[nombre]
        out_dir = os.path.join(base, "out_" + nombre)
        log("=== Bolsillo %s: centro (%.1f,%.1f,%.1f) size=%d ===" %
            (nombre, cx, cy, cz, size))
        affs = dock_pocket(receptor, args.ligandos, out_dir, cx, cy, cz, size,
                           args.exhaustividad, args.cpu_workers)
        log("  docking completado: %d resultados" % len(affs))

        act = [v for k, v in affs.items() if k.startswith("ACT_")]
        dec = [v for k, v in affs.items() if k.startswith("DEC_")]
        auc = roc_auc(act, dec)
        ef1 = enrichment_factor(act, dec, 1.0)
        ef5 = enrichment_factor(act, dec, 5.0)
        mean_act = np.mean(act) if act else 0.0
        mean_dec = np.mean(dec) if dec else 0.0
        log("  activos=%d decoys=%d | media act=%.2f dec=%.2f | AUC=%.3f EF1%%=%.2f EF5%%=%.2f" %
            (len(act), len(dec), mean_act, mean_dec, auc if auc is not None else -1, ef1, ef5))
        resumen.append({
            "bolsillo": nombre, "centro": "(%.1f,%.1f,%.1f)" % (cx, cy, cz),
            "size": size, "activos": len(act), "decoys": len(dec),
            "media_activos": round(mean_act, 2), "media_decoys": round(mean_dec, 2),
            "AUC": round(auc, 3) if auc is not None else -1,
            "EF1": round(ef1, 2), "EF5": round(ef5, 2),
        })

        # detalle por bolsillo
        det = os.path.join(base, "validacion_SOD1_%s.csv" % nombre)
        with open(det, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            w.writerow(["ligand", "rol", "affinity"])
            for k, v in sorted(affs.items(), key=lambda kv: kv[1]):
                rol = "activo" if k.startswith("ACT_") else "decoy"
                w.writerow([k, rol, v])

    # resumen global
    resumen_path = os.path.join(base, "resumen_bolsillos_sod1.csv")
    with open(resumen_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(resumen[0].keys()))
        w.writeheader()
        for r in resumen:
            w.writerow(r)
    log("Resumen guardado: %s" % resumen_path)
    log("=== FIN ===")


if __name__ == "__main__":
    main()
