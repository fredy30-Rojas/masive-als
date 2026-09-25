# -*- coding: utf-8 -*-
"""¿Cuántas moléculas de la librería se acoplaron SIN ser ellas mismas?

De dónde viene esto: la auditoría del 14 sep encontró que tres compuestos con
boro (CHEMBL15632, CHEMBL157117, CHEMBL4553125) tienen la pose con un CARBONO
donde debería estar el boro. El PDBQT trae el átomo NOMBRADO "B" pero TIPADO
como "C" en la última columna, que es la que manda en AutoDock:

    ATOM     47  B   UNL     1      52.451  82.543  66.028  ...  C

Así que esos compuestos se acoplaron y puntuaron como su análogo de carbono:
son inhibidores de proteasoma cuyo boro ES el grupo que se une al sitio activo.

Pero aquella auditoría sólo pudo comparar 60.750 de las 430.440 poses, porque
para las otras 369.652 no encontró el SMILES en ningún CSV. El "sólo hay tres"
no estaba demostrado: era "sólo se ven tres".

Este script cierra ese hueco por el lado que importa. El boro es el ÚNICO
elemento que puede sufrir esta conversión silenciosa, así que en vez de
auditar las 430.440 poses se buscan en la librería TODAS las moléculas que
contienen boro y se comprueba cada una contra su pose. Es un conjunto pequeño
y es el único donde este error es posible.

Cuenta también el caso hermano: las poses que son sólo un contraión (Br, Cl),
que son sales cuyo fragmento grande se perdió al preparar.

Uso:  python diag_boro.py            -> informe por pantalla
      python diag_boro.py --csv F    -> además, detalle a un CSV
"""
import argparse
import csv
import os
import re
import sys
from collections import defaultdict

BASE = r"C:\Users\Fredy\masive-als"
if not os.path.isdir(BASE):
    BASE = "/mnt/c/Users/Fredy/masive-als"

POSES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
CSVS_SMILES = [
    os.path.join(BASE, "analysis", "full_library_solo.smi"),
    os.path.join(BASE, "compounds", "decoys_library.smi"),
    os.path.join(BASE, "compounds", "full_library.smi"),
    os.path.join(BASE, "compounds", "libreria_extra_chembl.csv"),
    os.path.join(BASE, "compounds", "full_fda_library.csv"),
    os.path.join(BASE, "rescoring_datos", "chembl_smiles_extra.csv"),
]
DIANAS = ["SOD1", "TDP43", "TDP43_v2", "FUS"]

# Un boro de verdad: B mayúscula que no empieza otro símbolo (Br, Bi, Be) y que
# no va seguida de minúscula. [B-] con carga también cuenta.
BORO = re.compile(r"\[B|(?<![A-Za-z])B(?![a-z])")

# Tipos de AutoDock -> elemento. Aquí está la clave del asunto: el tipo manda.
TIPO_A_ELEMENTO = {
    "C": "C", "A": "C", "N": "N", "NA": "N", "O": "O", "OA": "O",
    "S": "S", "SA": "S", "P": "P", "F": "F", "CL": "CL", "BR": "BR",
    "I": "I", "B": "B", "HD": "H", "H": "H", "FE": "FE", "MG": "MG",
    "ZN": "ZN", "MN": "MN", "CA": "CA", "CU": "CU", "SI": "SI",
}
# Una pose de ligando de verdad ocupa bastante más que esto. Un contraión
# suelto, con sus diez modelos, no llega ni de lejos.
BYTES_MINIMOS = 4000


def leer_smiles():
    """SMILES de toda la librería, venga del formato que venga."""
    tabla = {}
    for ruta in CSVS_SMILES:
        if not os.path.exists(ruta):
            continue
        try:
            with open(ruta, encoding="utf-8", errors="replace") as f:
                primera = f.readline()
                f.seek(0)
                if "," in primera and "smiles" in primera.lower():
                    for d in csv.DictReader(f):
                        smi = _campo(d, ("smiles", "SMILES", "canonical_smiles"))
                        nom = _campo(d, ("chembl_id", "id", "name", "ID", "cid"))
                        if smi and nom:
                            tabla.setdefault(nom.strip(), smi.strip())
                else:
                    for linea in f:
                        partes = linea.rstrip("\r\n").split("\t")
                        if len(partes) >= 2 and partes[0].strip():
                            tabla.setdefault(partes[1].strip(), partes[0].strip())
        except OSError:
            continue
    return tabla


def _campo(d, nombres):
    for n in nombres:
        if d.get(n):
            return d[n]
    return ""


def composicion_smiles(smiles):
    """Átomos pesados del fragmento MAYOR, por elemento. Los H no cuentan."""
    from rdkit import Chem

    m = Chem.MolFromSmiles(smiles) if smiles else None
    if m is None:
        return None
    trozos = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    if not trozos:
        return None
    mayor = max(trozos, key=lambda t: t.GetNumHeavyAtoms())
    c = defaultdict(int)
    for a in mayor.GetAtoms():
        c[a.GetSymbol().upper()] += 1
    return dict(c)


def composicion_pose(ruta):
    """Átomos pesados del PRIMER modelo de la pose, contados por el TIPO de
    AutoDock (la última columna), que es la que decide cómo se puntúa."""
    c = defaultdict(int)
    if not os.path.exists(ruta):
        return None
    with open(ruta, encoding="utf-8", errors="replace") as f:
        for linea in f:
            if linea.startswith("ENDMDL"):
                break
            if not linea.startswith(("ATOM", "HETATM")):
                continue
            tipo = linea[77:79].strip().upper()
            if not tipo:
                continue
            ele = TIPO_A_ELEMENTO.get(tipo)
            if ele is None:
                # Tipo desconocido: se cuenta por el nombre como último recurso.
                ele = re.sub(r"[^A-Za-z]", "", linea[12:16]).upper()[:2] or "?"
                if ele not in ("BR", "CL"):
                    ele = ele[0]
            if ele == "H":
                continue
            c[ele] += 1
    return dict(c) if c else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="")
    a = ap.parse_args()

    print("Leyendo la librería de SMILES...")
    libreria = leer_smiles()
    print("  %d moléculas con SMILES conocido" % len(libreria))
    con_boro = {k: v for k, v in libreria.items() if BORO.search(v)}
    print("  de ellas, %d contienen boro  <-- el único grupo donde este error "
          "es posible" % len(con_boro))

    filas = []
    for diana in DIANAS:
        carpeta = os.path.join(POSES, "results_" + diana)
        if not os.path.isdir(carpeta):
            continue
        for lig, smi in con_boro.items():
            ruta = os.path.join(carpeta, lig + "_out.pdbqt")
            if not os.path.exists(ruta):
                continue
            quiero = composicion_smiles(smi)
            tengo = composicion_pose(ruta)
            if quiero is None or tengo is None:
                continue
            boro_smi = quiero.get("B", 0)
            boro_pose = tengo.get("B", 0)
            filas.append({
                "diana": diana, "ligando": lig, "boro_smiles": boro_smi,
                "boro_pose": boro_pose,
                "perdido": boro_smi > 0 and boro_pose < boro_smi,
                "carbonos_de_mas": tengo.get("C", 0) - quiero.get("C", 0),
                "atomos_pose": sum(tengo.values()),
                "smiles": smi,
            })

    print()
    print("=" * 74)
    print("MOLÉCULAS CON BORO QUE TIENEN POSE")
    print("=" * 74)
    if not filas:
        print("  ninguna: o no hay poses de estos compuestos, o no se "
              "emparejaron los nombres.")
        return 0

    perdidas = [f for f in filas if f["perdido"]]
    conservadas = [f for f in filas if not f["perdido"]]
    print("  poses encontradas      : %d" % len(filas))
    print("  con el boro PERDIDO    : %d" % len(perdidas))
    print("  con el boro conservado : %d" % len(conservadas))
    if conservadas:
        print()
        print("  (las conservadas, por si acaso:)")
        for f in conservadas[:10]:
            print("     %-12s %-18s boro_smiles=%d boro_pose=%d"
                  % (f["diana"], f["ligando"], f["boro_smiles"], f["boro_pose"]))
    print()
    por_lig = defaultdict(set)
    for f in perdidas:
        por_lig[f["ligando"]].add(f["diana"])
    print("  son %d COMPUESTOS distintos afectados:" % len(por_lig))
    for lig in sorted(por_lig):
        extra = [f["carbonos_de_mas"] for f in perdidas if f["ligando"] == lig]
        print("     %-18s en %-28s (carbonos de más en la pose: %s)"
              % (lig, ",".join(sorted(por_lig[lig])), extra[0]))

    # Caso hermano: poses que son sólo un contraión. Se filtra por tamaño de
    # fichero antes de abrir nada, porque son 430.000 ficheros.
    print()
    print("=" * 74)
    print("POSES SOSPECHOSAMENTE PEQUEÑAS (posible contraión suelto)")
    print("=" * 74)
    for diana in DIANAS:
        carpeta = os.path.join(POSES, "results_" + diana)
        if not os.path.isdir(carpeta):
            continue
        pequenas = []
        for nombre in os.listdir(carpeta):
            if not nombre.endswith("_out.pdbqt"):
                continue
            ruta = os.path.join(carpeta, nombre)
            try:
                if os.path.getsize(ruta) < BYTES_MINIMOS:
                    pequenas.append(nombre[:-len("_out.pdbqt")])
            except OSError:
                continue
        print("  %-10s %4d poses de menos de %d bytes" % (diana, len(pequenas),
                                                          BYTES_MINIMOS))
        if pequenas:
            print("             %s" % ", ".join(sorted(pequenas)[:8]))

    if a.csv:
        with open(a.csv, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(filas[0].keys()))
            w.writeheader()
            w.writerows(filas)
        print()
        print("detalle escrito en %s" % a.csv)
    return 0


if __name__ == "__main__":
    sys.exit(main())
