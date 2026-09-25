# -*- coding: utf-8 -*-
"""Construye la lista de VALIDACIÓN del rescoring: controles + fondo.

Por qué (auditoría del 11/09/2026, hallazgo 7): los 20 controles positivos
conocidos de SOD1 y TDP-43 quedaron TODOS fuera del corte top-5 % del
cribado, así que ninguno entró en la lista corta ni fue reescalado. Sin un
solo activo conocido dentro del conjunto, no hay forma de medir si el
rescoring ordena bien: no se puede calcular enriquecimiento.

Qué se mide con esta lista: si los controles positivos conocidos reciben
un dG más negativo que un fondo aleatorio de compuestos de la propia lista
corta (todos acoplados, todos pasando los mismos filtros). Se reporta el
AUC (probabilidad de que un control puntúe mejor que un fondo) y, como el
dG sin desolvación premia el tamaño molecular, la correlación con el peso
para no confundir discriminación real con sesgo de tamaño.

No es un ensayo de enriquecimiento completo (harían falta señuelos
físicamente apareados, tipo DUD-E), sino la comprobación mínima e
imprescindible antes de ordenar candidatos: ¿el instrumento distingue lo
que sabemos que es activo de lo que no sabemos nada?

Uso:
    python preparar_validacion.py --fondo 60
"""
import argparse
import csv
import os
import random

BASE = r"C:\Users\Fredy\masive-als"
if not os.path.isdir(BASE):
    BASE = "/mnt/c/Users/Fredy/masive-als"
ANALYSIS = os.path.join(BASE, "analysis")
LOCAL = os.path.join(ANALYSIS, "rescoring_local")
RES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
CONTROLES = os.path.join(ANALYSIS, "controles_calibracion.csv")
CORTOS = os.path.join(ANALYSIS, "lista_corta_candidatos.csv")
SALIDA = os.path.join(LOCAL, "validacion_controles.csv")


def pose(t, lig):
    return os.path.join(RES, "results_" + t, lig + "_out.pdbqt")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fondo", type=int, default=60,
                    help="compuestos de fondo por diana (de la lista corta)")
    ap.add_argument("--semilla", type=int, default=42)
    a = ap.parse_args()
    rnd = random.Random(a.semilla)

    filas, usados = [], set()

    # 1) Controles positivos con pose acoplada disponible
    sin_pose = []
    with open(CONTROLES, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r.get("en_libreria") != "si":
                continue
            t, lig = r["target"], r["ligand_libreria"]
            # TDP-43: la diana del rescoring es TDP43_v2 (misma pose, caja v2)
            t_resc = "TDP43_v2" if t in ("TDP43", "TDP43_v2") else t
            if not os.path.exists(pose(t_resc, lig)):
                sin_pose.append((t, lig))
                continue
            filas.append({"target": t_resc, "ligand": lig,
                          "afinidad": r.get("afinidad_vina", ""),
                          "smiles": r.get("smiles", ""), "tipo": "control"})
            usados.add((t_resc, lig))

    # 2) Fondo: compuestos de la lista corta (misma diana), al azar con semilla
    cortos = {}
    with open(CORTOS, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            cortos.setdefault(r["target"], []).append(r)
    for t, rows in sorted(cortos.items()):
        libres = [r for r in rows if (t, r["ligand"]) not in usados
                  and os.path.exists(pose(t, r["ligand"]))]
        for r in rnd.sample(libres, min(a.fondo, len(libres))):
            filas.append({"target": t, "ligand": r["ligand"],
                          "afinidad": r["afinidad"], "smiles": r["smiles"],
                          "tipo": "fondo"})

    with open(SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["target", "ligand", "afinidad",
                                          "smiles", "tipo"])
        w.writeheader()
        w.writerows(filas)

    from collections import Counter
    c = Counter((r["target"], r["tipo"]) for r in filas)
    print("escrito: %s (%d filas)" % (SALIDA, len(filas)))
    for k in sorted(c):
        print("  %-12s %-8s %d" % (k[0], k[1], c[k]))
    if sin_pose:
        print("controles SIN pose (excluidos): %s" % (sin_pose,))


if __name__ == "__main__":
    main()
