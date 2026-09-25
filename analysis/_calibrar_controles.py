# -*- coding: utf-8 -*-
"""
Calibracion del corte del ranking MASIVE-ALS con CONTROLES POSITIVOS.

Cruce por SMILES canonico (RDKit) entre:
  - controles positivos conocidos (activos_*.csv: TDP43, SOD1, FUS)
  - la libreria cribada en el GPU (gpu_dock/libreria_ligands/*.pdbqt con su
    nombre, y los resultados results_<target>/<nombre>_out.pdbqt)

Salida:
  analysis/controles_calibracion.csv  -> por control: en_libreria?, afinidad Vina real
  analysis/CALIBRACION_CONTROLES.md   -> informe con el corte recomendado
"""
import glob
import os
import csv
import sys

from rdkit import Chem
from rdkit.Chem import RDKFingerprint
from rdkit import DataStructs

BASE = r"C:\Users\Fredy\masive-als"
LIGDIR = os.path.join(BASE, "gpu_dock", "libreria_ligands")
RESDIR = os.path.join(BASE, "gpu_dock", "resultados_libreria")
ANALYSIS = os.path.join(BASE, "analysis")

CONTROLES = {
    "TDP43": os.path.join(ANALYSIS, "activos_tdp43_v2.csv"),
    "SOD1": os.path.join(ANALYSIS, "activos_sod1.csv"),
    "FUS": os.path.join(ANALYSIS, "activos_fus.csv"),
}


def can(smi):
    try:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            return None
        return Chem.MolToSmiles(m)
    except Exception:
        return None


def leer_afinidad(pdbqt_path):
    """Primer REMARK VINA RESULT (kcal/mol)."""
    try:
        with open(pdbqt_path, encoding="utf-8", errors="replace") as f:
            for l in f:
                if l.startswith("REMARK VINA RESULT"):
                    return float(l.split()[3])
    except Exception:
        pass
    return None


def main():
    # 1) indice nombre -> smiles canonico de la libreria cribada
    lib_idx = {}   # nombre -> can
    lib_can = {}   # can -> nombre
    for p in glob.glob(os.path.join(LIGDIR, "*.pdbqt")):
        nombre = os.path.basename(p).replace(".pdbqt", "")
        # el SMILES no esta en el pdbqt; lo buscamos en full_library_solo.smi
        # (construido el 20/08; puede no cubrir todo). Mejor: nombre directo.
        lib_idx[nombre] = None
    print("libreria pdbqt:", len(lib_idx))

    # 2) indice can -> nombre desde full_library_solo.smi
    smi_path = os.path.join(ANALYSIS, "full_library_solo.smi")
    n_smi = 0
    if os.path.exists(smi_path):
        with open(smi_path, encoding="utf-8", errors="replace") as f:
            for ln in f:
                ln = ln.strip()
                if not ln or "\t" not in ln:
                    continue
                smi, nombre = ln.rsplit("\t", 1)
                c = can(smi)
                if c and nombre in lib_idx:
                    lib_idx[nombre] = c
                    lib_can.setdefault(c, nombre)
                    n_smi += 1
    print("libreria con SMILES resuelto:", n_smi)

    filas = []
    for target, path in CONTROLES.items():
        with open(path, encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        print("\n[%s] controles: %d" % (target, len(rows)))
        resdir = os.path.join(RESDIR, "results_" + target)
        for r in rows:
            nombre = r["name"].strip()
            smi = r["smiles"].strip()
            c = can(smi)
            # buscar en libreria por nombre directo o por canonico
            en_libreria = "no"
            lig_nombre = None
            if c and c in lib_can:
                en_libreria = "si"
                lig_nombre = lib_can[c]
            elif nombre in lib_idx:
                en_libreria = "si"
                lig_nombre = nombre
            # afinidad real si existe output
            aff = None
            if lig_nombre:
                out = os.path.join(resdir, lig_nombre + "_out.pdbqt")
                if os.path.exists(out):
                    aff = leer_afinidad(out)
            filas.append({
                "target": target, "control": nombre, "smiles": smi,
                "en_libreria": en_libreria, "ligand_libreria": lig_nombre or "",
                "afinidad_vina": "" if aff is None else "%.3f" % aff,
            })
            estado = ("aff=%.3f" % aff) if aff is not None else ("en_libreria pero sin output" if en_libreria == "si" else "NO en libreria")
            print("  %-22s %s" % (nombre, estado))

    out = os.path.join(ANALYSIS, "controles_calibracion.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
        w.writeheader()
        w.writerows(filas)
    print("\nCSV:", out)


if __name__ == "__main__":
    main()