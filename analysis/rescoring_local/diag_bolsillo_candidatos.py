# -*- coding: utf-8 -*-
"""¿El candidato estrella esta en el mismo bolsillo que los farmacos conocidos?

Sospecha concreta: CHEMBL520254 sale primero en SOD1 por las tres metricas
(dG -72,01, residual -37,6, mas de 13 kcal/mol por delante del siguiente),
pero su nota de Vina es solo -6,3, que es mediocre. Esa mezcla -- Vina
indiferente y MM-GBSA extraordinario -- suele significar que la pose se metio
en una cavidad que no es el sitio de union, donde se entierra muchisimo y el
MM-GBSA lo confunde con afinidad.

Se mide, para cada pose de la diana:
  - donde esta su centro
  - cuantos atomos de receptor tiene a menos de 4,5 A (enterramiento)
  - cual es su distancia minima al receptor (menos de 2,2 A = choque)
  - cuantos atomos tiene

Y se compara el candidato contra los controles y contra el fondo.

Uso:  python diag_bolsillo_candidatos.py --target SOD1 [--candidato CHEMBL520254]
"""
import argparse
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
RADIO_CONTACTO = 4.5


def leer_pdbqt(path, solo_pesados=True):
    pts = []
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            if solo_pesados and line[77:80].strip().upper().startswith("H"):
                continue
            pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
        except (ValueError, IndexError):
            continue
    return np.array(pts) if pts else None


def leer_receptor(path):
    pts = []
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith(("ATOM", "HETATM")):
            try:
                pts.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
            except (ValueError, IndexError):
                continue
    return np.array(pts) if pts else None


def num(txt, por_defecto=float("nan")):
    """Convierte a numero sin reventar: hay filas con la casilla vacia."""
    try:
        return float(txt)
    except (TypeError, ValueError):
        return por_defecto


def medir(pose, receptor):
    """Centro, contactos, distancia minima y radio de giro de una pose."""
    cen = pose.mean(axis=0)
    # Distancias de cada atomo de ligando a cada atomo de receptor.
    d = np.linalg.norm(pose[:, None, :] - receptor[None, :, :], axis=2)
    dmin = float(d.min())
    tocados = int((d.min(axis=0) < RADIO_CONTACTO).sum())  # atomos de receptor
    rg = float(np.sqrt(((pose - cen) ** 2).sum(axis=1).mean()))
    return cen, dmin, tocados, rg


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="SOD1")
    ap.add_argument("--candidato", default="CHEMBL520254")
    a = ap.parse_args()

    receptor = leer_receptor(os.path.join(
        LOCAL, "receptores_fijos", "snapshot_20260911", "%s_fijo.pdb" % a.target))
    if receptor is None:
        print("no encuentro el receptor fijo de %s" % a.target)
        return 1

    tipos = {}
    if os.path.exists(VAL):
        # OJO: este CSV trae controles Y fondo, con una columna `tipo` que lo
        # dice. Hay que leer el VALOR, no conformarse con que el ligando este
        # en el fichero: si no, todo el fondo sale etiquetado de control.
        for r in csv.DictReader(open(VAL, encoding="utf-8")):
            tipos[(r["target"], r["ligand"])] = (r.get("tipo") or "fondo").strip()

    filas = [r for r in csv.DictReader(open(CSV, encoding="utf-8"))
             if r["mmgbsa_dG"] and r["target"] == a.target]

    datos = []
    for r in filas:
        p = os.path.join(POSES, "results_" + a.target, r["ligand"] + "_out.pdbqt")
        if not os.path.exists(p):
            continue
        pose = leer_pdbqt(p)
        if pose is None:
            continue
        cen, dmin, tocados, rg = medir(pose, receptor)
        datos.append({
            "lig": r["ligand"], "cen": cen, "dmin": dmin, "tocados": tocados,
            "rg": rg, "nat": len(pose), "dg": num(r["mmgbsa_dG"]),
            "vina": num(r["vina_affinity"]),
            "tipo": tipos.get((a.target, r["ligand"]), "fondo"),
        })

    print("=== %s: %d poses medidas contra %d atomos de receptor ==="
          % (a.target, len(datos), len(receptor)))
    print()

    # Centro de referencia: donde se apuntó el acoplamiento (centro del fondo).
    fondos = [d for d in datos if d["tipo"] == "fondo"]
    ctrl = [d for d in datos if d["tipo"] == "control"]
    ref = np.array([d["cen"] for d in fondos]).mean(axis=0) if fondos else None

    def linea(d):
        dist = float(np.linalg.norm(d["cen"] - ref)) if ref is not None else 0.0
        return ("  %-22s %-7s nat=%3d  centro a %5.1f A  contactos=%4d  "
                "min=%4.2f A  dG=%+7.2f  vina=%+.1f"
                % (d["lig"][:22], d["tipo"], d["nat"], dist, d["tocados"],
                   d["dmin"], d["dg"], d["vina"]))

    cand = [d for d in datos if d["lig"] == a.candidato]
    if cand:
        print("EL CANDIDATO:")
        for d in cand:
            print(linea(d))
        # Comparar con el fondo en presupuesto de atomos parecido.
        nat = cand[0]["nat"]
        pares = [d for d in fondos if abs(d["nat"] - nat) <= 3]
        print()
        print("  con un tamano parecido (%d-%d atomos) hay %d del fondo:"
              % (nat - 3, nat + 3, len(pares)))
        if pares:
            tc = np.array([d["tocados"] for d in pares])
            print("     contactos: mediana %d  (min %d, max %d)"
                  % (np.median(tc), tc.min(), tc.max()))
            puesto = int((tc > cand[0]["tocados"]).sum()) + 1
            print("     el candidato tiene %d contactos: puesto %d de %d"
                  % (cand[0]["tocados"], puesto, len(pares) + 1))
        print()

    print("LOS 5 CON MAS ENTERRAMIENTO DE LA DIANA:")
    for d in sorted(datos, key=lambda x: -x["tocados"])[:5]:
        print(linea(d))
    print()
    print("LOS 5 CONTROLES:")
    for d in sorted(ctrl, key=lambda x: -x["tocados"])[:5]:
        print(linea(d))
    print()

    # ¿Hay choques? Una distancia minima muy corta delata una pose metida a
    # la fuerza dentro del receptor.
    chocan = [d for d in datos if d["dmin"] < 2.2]
    print("poses con choque esterico (menos de 2,2 A): %d de %d"
          % (len(chocan), len(datos)))
    for d in sorted(chocan, key=lambda x: x["dmin"])[:6]:
        print(linea(d))
    return 0


if __name__ == "__main__":
    sys.exit(main())
