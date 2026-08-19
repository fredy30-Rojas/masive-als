# -*- coding: utf-8 -*-
"""Anade la columna pose_pdbqt a los CSVs de candidatos, mapeando cada
ligando al PDBQT de su pose de docking en las tandas locales."""
import csv
import os

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def buscar_pose(lig, tgt):
    for tanda in ["z001", "z002", "z003"]:
        d = os.path.join(BASE, "gpu_dock", "tanda_%s" % tanda, "results_%s" % tgt)
        p = os.path.join(d, lig + "_out.pdbqt")
        if os.path.exists(p):
            return p.replace("\\", "/")
    return ""


def main():
    for f in ["analysis/candidatos_filtrados.csv",
              "analysis/candidatos_total_cns.csv"]:
        ruta = os.path.join(BASE, f)
        if not os.path.exists(ruta):
            print(f, "no existe")
            continue
        rows = []
        with open(ruta, encoding="utf-8", errors="replace") as fh:
            rdr = csv.DictReader(fh)
            cols = list(rdr.fieldnames)
            for row in rdr:
                row["pose_pdbqt"] = buscar_pose(row["ligand"], row["target"])
                rows.append(row)
        outcols = cols + (["pose_pdbqt"] if "pose_pdbqt" not in cols else [])
        with open(ruta, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=outcols)
            w.writeheader()
            for r in rows:
                w.writerow(r)
        n_pose = sum(1 for r in rows if r["pose_pdbqt"])
        print("%s: %d filas, %d con pose encontrada" % (f, len(rows), n_pose))


if __name__ == "__main__":
    main()
