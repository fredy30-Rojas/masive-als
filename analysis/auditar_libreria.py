# -*- coding: utf-8 -*-
"""¿Hay más entradas de la librería mal etiquetadas, como la cefarantina?

El 14 sep 2026 se descubrió que la pose guardada como `cepharanthine` tiene 47
átomos pesados y 39 carbonos, cuando la cefarantina (C37H38N2O6) tiene 45
átomos y 37 carbonos. La pose es de OTRA molécula.

Si eso pasa una vez, puede pasar más veces, y entonces el cribado entero se
apoya en poses que no son de quien dicen ser. Este script audita TODA la
librería: para cada pose compara su composición de elementos pesados contra la
del SMILES de su ligando.

No hace falta ningún cálculo caro: sólo se leen las cabeceras de los PDBQT y la
fórmula del SMILES.

Uso:  python auditar_libreria.py [--dir <raiz masive-als>]
"""
import argparse
import collections
import csv
import os
import sys

from rdkit import Chem

# Tipos de átomo de AutoDock -> elemento
AD2ELEM = {
    "C": "C", "A": "C",
    "N": "N", "NA": "N", "NS": "N", "NX": "N",
    "O": "O", "OA": "O", "OS": "O",
    "S": "S", "SA": "S",
    "P": "P", "F": "F", "CL": "Cl", "BR": "Br", "I": "I",
    "H": "H", "HD": "H", "HS": "H",
}
METALES = {"CU", "ZN", "CA", "FE", "MG", "MN", "CO", "NI", "CD", "NA", "K"}


def composicion_pose(path):
    """Elementos pesados de la PRIMERA pose del PDBQT -> Counter.

    Se normaliza el simbolo a la grafia de RDKit (CL->Cl, SE->Se) y se quitan
    los metales, igual que se hace en el SMILES, para que la comparacion sea
    simetrica.
    """
    c = collections.Counter()
    for line in open(path, encoding="utf-8", errors="replace"):
        if line.startswith("ENDMDL"):
            break
        if not line.startswith(("ATOM", "HETATM")):
            continue
        t = line[77:80].strip().upper()
        e = AD2ELEM.get(t, t)
        e = e.capitalize() if len(e) > 1 else e
        if e == "H" or e.upper() in METALES:
            continue
        c[e] += 1
    return c


def composicion_smiles(smi):
    """Elementos pesados del fragmento mayor del SMILES -> Counter o None."""
    m = Chem.MolFromSmiles(smi) if smi else None
    if m is None:
        return None
    trozos = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    if not trozos:
        return None
    mayor = max(trozos, key=lambda x: x.GetNumHeavyAtoms())
    c = collections.Counter()
    for a in mayor.GetAtoms():
        n = a.GetAtomicNum()
        if n > 1 and a.GetSymbol().upper() not in METALES:
            c[a.GetSymbol()] += 1
    return c


def clasificar(ref, pose):
    """Qué tipo de desajuste es. No todos significan lo mismo."""
    nref, npose = sum(ref.values()), sum(pose.values())
    if npose <= 3:
        return "pose degenerada (casi vacia)"
    if npose < 0.6 * nref:
        return "pose fragmentada (le falta medio ligando)"
    if not (pose - ref):
        return "pose es un trozo del ligando"
    return "MOLECULA DISTINTA"



def buscar_smiles(raiz):
    """Ligando -> SMILES, de todos los CSV que tengan esas columnas."""
    smi = {}
    origen = {}
    for base, dirs, ficheros in os.walk(raiz):
        dirs[:] = [d for d in dirs
                   if d not in (".git", "__pycache__", "node_modules",
                                "rescoring_env", ".venv")]
        for fn in ficheros:
            if not fn.lower().endswith(".csv"):
                continue
            p = os.path.join(base, fn)
            try:
                if os.path.getsize(p) > 80 * 1024 * 1024:
                    continue
                with open(p, encoding="utf-8", errors="replace") as f:
                    cab = f.readline()
                    cols = [c.strip().strip('"').lower() for c in cab.split(",")]
                    if "smiles" not in cols:
                        continue
                    i_smi = cols.index("smiles")
                    i_lig = None
                    for cand in ("ligand", "ligando", "nombre", "name", "id",
                                 "chembl_id", "compound"):
                        if cand in cols:
                            i_lig = cols.index(cand)
                            break
                    if i_lig is None:
                        continue
                    for fila in csv.reader(f):
                        if len(fila) <= max(i_smi, i_lig):
                            continue
                        lig, s = fila[i_lig].strip(), fila[i_smi].strip()
                        if lig and s and lig not in smi:
                            smi[lig] = s
                            origen[lig] = os.path.relpath(p, raiz)
            except (OSError, UnicodeDecodeError, csv.Error):
                continue
    return smi, origen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default=r"C:\Users\Fredy\masive-als")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    raiz = a.dir
    lib = os.path.join(raiz, "gpu_dock", "resultados_libreria")
    if not os.path.isdir(lib):
        print("no encuentro la libreria en %s" % lib)
        return 1

    print("Buscando SMILES en %s ..." % raiz, flush=True)
    smi, origen = buscar_smiles(raiz)
    print("ligandos con SMILES conocido: %d" % len(smi), flush=True)

    malos, sin_smi, revisados, ok = [], [], 0, 0
    dianas = sorted(d for d in os.listdir(lib)
                    if d.startswith("results_") and
                    os.path.isdir(os.path.join(lib, d)))
    for d in dianas:
        diana = d[len("results_"):]
        ficheros = sorted(f for f in os.listdir(os.path.join(lib, d))
                          if f.endswith("_out.pdbqt"))
        print("  %-12s %d poses" % (diana, len(ficheros)), flush=True)
        for f in ficheros:
            lig = f[:-len("_out.pdbqt")]
            revisados += 1
            if lig not in smi:
                sin_smi.append((diana, lig))
                continue
            pose = composicion_pose(os.path.join(lib, d, f))
            ref = composicion_smiles(smi[lig])
            if ref is None:
                sin_smi.append((diana, lig))
                continue
            if pose == ref:
                ok += 1
            else:
                malos.append((diana, lig, ref, pose, smi[lig],
                              clasificar(ref, pose)))

    def fmt(c):
        return " ".join("%s%d" % (k, v) for k, v in sorted(c.items()))

    print()
    print("=" * 74)
    print("poses revisadas      : %d" % revisados)
    print("  coinciden          : %d" % ok)
    print("  NO COINCIDEN       : %d" % len(malos))
    print("  sin SMILES conocido: %d" % len(sin_smi))
    print("=" * 74)

    if malos:
        cuenta = collections.Counter(c for *_, c in malos)
        print()
        print("TIPOS DE DESAJUSTE:")
        for tipo, n in cuenta.most_common():
            print("  %-42s %6d" % (tipo, n))
        graves = [m for m in malos if m[5] == "MOLECULA DISTINTA"]
        print()
        print("MOLECULA DISTINTA (la pose es de otro compuesto de tamano "
              "parecido): %d" % len(graves))
        for diana, lig, ref, pose, s, _ in graves[:25]:
            print("  %-12s %-24s SMILES[%s]  pose[%s]"
                  % (diana, lig[:24], fmt(ref), fmt(pose)))
        if len(graves) > 25:
            print("  ... y %d mas" % (len(graves) - 25))

    if a.out:
        with open(a.out, "w", encoding="utf-8") as fh:
            fh.write("diana,ligando,tipo,smiles,composicion_smiles,"
                     "composicion_pose\n")
            for diana, lig, ref, pose, s, tipo in malos:
                fh.write("%s,%s,%s,%s,%s,%s\n"
                         % (diana, lig, tipo, s, fmt(ref), fmt(pose)))
        print()
        print("detalle escrito en %s" % a.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
