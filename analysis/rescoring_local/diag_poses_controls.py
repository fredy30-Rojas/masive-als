# -*- coding: utf-8 -*-
"""¿Los controles vienen acoplados al mismo sitio que el fondo?

Sospecha: `extract_receptor_pdb` construye el bolsillo alrededor del centroide
de la pose DE CADA LIGANDO (residuos a menos de 12 A de ese centroide). Si un
control quedó acoplado en otro sitio del receptor, se le fabrica un bolsillo
local propio y el MM-GBSA no lo castiga por estar mal colocado: sólo mide la
energía de interacción local de la pose que le dieron.

Si eso pasa, los controles no puntúan mal por el scoring, sino porque Vina los
dejó fuera del sitio bueno. Y entonces el arreglo está en el acoplamiento, no
en la métrica.

Mide, por diana: centroide de cada pose, distancia entre el centro del grupo de
controles y el del fondo, dispersión de cada grupo, y distancia de cada pose al
átomo de receptor más cercano (para descartar poses en el vacío).
"""
import csv
import os
import sys

import numpy as np

BASE = r"C:\Users\Fredy\masive-als"
if not os.path.isdir(BASE):
    BASE = "/mnt/c/Users/Fredy/masive-als"
LOCAL = os.path.join(BASE, "analysis", "rescoring_local")
POSES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
CSV = os.path.join(LOCAL, "rescoring_recalculo_validado.csv")
VAL = os.path.join(LOCAL, "validacion_controles.csv")


def atomos(pose_path):
    """Coordenadas de los atomos pesados de la primera pose del PDBQT."""
    pts = []
    for line in open(pose_path, encoding="utf-8", errors="replace"):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            ad = line[77:80].strip().upper()
            if ad.startswith("H"):
                continue
            pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except (ValueError, IndexError):
            continue
    return np.array(pts) if pts else None


def atomos_receptor(path):
    pts = []
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            try:
                pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except (ValueError, IndexError):
                continue
    return np.array(pts) if pts else None


def main():
    tipos = {(r["target"], r["ligand"]): r["tipo"]
             for r in csv.DictReader(open(VAL, encoding="utf-8"))}
    filas = [r for r in csv.DictReader(open(CSV, encoding="utf-8"))
             if r["mmgbsa_dG"]]

    for t in ("TDP43_v2", "SOD1", "FUS"):
        rec = os.path.join(LOCAL, "receptores_fijos", "snapshot_20260911",
                           "%s_fijo.pdb" % t)
        rec_pts = atomos_receptor(rec) if os.path.exists(rec) else None
        grupos = {"control": [], "fondo": []}
        faltan = 0
        for r in filas:
            if r["target"] != t:
                continue
            p = os.path.join(POSES, "results_" + t, r["ligand"] + "_out.pdbqt")
            if not os.path.exists(p):
                faltan += 1
                continue
            pts = atomos(p)
            if pts is None:
                faltan += 1
                continue
            cen = pts.mean(axis=0)
            dmin = float(np.linalg.norm(rec_pts - cen, axis=1).min()) if rec_pts is not None else float("nan")
            grupos[tipos.get((t, r["ligand"]), "fondo")].append((r["ligand"], cen, dmin, len(pts)))

        print("=== %s ===" % t)
        print("  poses leidas: control %d | fondo %d | sin pose %d"
              % (len(grupos["control"]), len(grupos["fondo"]), faltan))
        if not grupos["fondo"]:
            print("  (sin poses de fondo; nada que comparar)\n")
            continue
        ctrl_c = np.array([c for _, c, _, _ in grupos["control"]]) if grupos["control"] else None
        fond_c = np.array([c for _, c, _, _ in grupos["fondo"]])
        if ctrl_c is not None and len(ctrl_c):
            # distancia entre el centro de cada grupo y su dispersion interna
            dif = np.linalg.norm(ctrl_c.mean(axis=0) - fond_c.mean(axis=0))
            dis_c = float(np.linalg.norm(ctrl_c - ctrl_c.mean(axis=0), axis=1).mean())
            dis_f = float(np.linalg.norm(fond_c - fond_c.mean(axis=0), axis=1).mean())
            print("  centro controles vs centro fondo: %.1f A" % dif)
            print("  dispersion interna: controles %.1f A | fondo %.1f A" % (dis_c, dis_f))
        else:
            print("  (sin controles con pose)")
        # cuantos de cada grupo caen lejos del atomo de receptor mas cercano
        for etq in ("control", "fondo"):
            g = grupos[etq]
            if not g:
                continue
            d = np.array([x[2] for x in g])
            print("  %-8s distancia al receptor: mediana %.1f A  (%.1f a %.1f)"
                  % (etq, np.median(d), d.min(), d.max()))
        print()

    return 0


if __name__ == "__main__":
    sys.exit(main())
