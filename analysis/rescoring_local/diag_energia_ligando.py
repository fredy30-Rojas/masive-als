# -*- coding: utf-8 -*-
"""¿Cuánto del dG de MM-GBSA es afinidad y cuánto es el ligando retorciéndose?

En el recalculo salio un candidato sospechoso, CHEMBL520254 en SOD1, con
dG = -72,01 cuando el siguiente esta a mas de 13 kcal/mol. Mirando la fila
cruda aparece esto:

    e_complex = -3814,35   e_receptor = -3852,38   e_ligand = +110,04
    dG = e_complex - e_receptor - e_ligand = -72,01

El termino `e_ligand` es la energia interna del ligando SOLO, en la geometria
que le dejo el acoplamiento. Un valor positivo y grande significa que la
conformacion esta forzada: la molecula esta retorcida, con los grupos
empujandose entre si. Eso no es afinidad, es tension.

Este script mide si eso es un caso suelto o un problema de toda la tabla, y
si contamina el ranking. Si `e_ligand` varia mucho de un ligando a otro
entonces el dG no se puede comparar entre moleculas distintas.

Uso:  python diag_energia_ligando.py [--csv <ruta>]
"""
import argparse
import csv
import os
import statistics as st

ARCHIVO = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "rescoring_recalculo_validado.csv")


def cargar(ruta):
    filas = []
    with open(ruta, encoding="utf-8", errors="replace") as f:
        for d in csv.DictReader(f):
            if (d.get("error") or "").strip():
                continue
            try:
                filas.append({
                    "target": d["target"],
                    "ligand": d["ligand"],
                    "smiles": d.get("smiles", ""),
                    "vina": float(d["vina_affinity"]),
                    "dg": float(d["mmgbsa_dG"]),
                    "sd": float(d["dG_sd"]) if d.get("dG_sd") else 0.0,
                    "e_c": float(d["e_complex"]),
                    "e_r": float(d["e_receptor"]),
                    "e_l": float(d["e_ligand"]),
                })
            except (KeyError, ValueError):
                continue
    return filas


def atomos(smiles):
    try:
        from rdkit import Chem
    except ImportError:
        return None
    m = Chem.MolFromSmiles(smiles) if smiles else None
    if m is None:
        return None
    trozos = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    return max(t.GetNumHeavyAtoms() for t in trozos) if trozos else None


def spearman(xs, ys):
    """Correlacion de rangos, sin depender de scipy."""
    if len(xs) < 3:
        return 0.0

    def rangos(v):
        orden = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(orden):
            j = i
            while j + 1 < len(orden) and v[orden[j + 1]] == v[orden[i]]:
                j += 1
            medio = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[orden[k]] = medio
            i = j + 1
        return r

    rx, ry = rangos(xs), rangos(ys)
    mx, my = st.mean(rx), st.mean(ry)
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else 0.0


def mediana(v):
    return st.median(v) if v else float("nan")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=ARCHIVO)
    a = ap.parse_args()
    filas = cargar(a.csv)
    if not filas:
        print("no he podido leer ninguna fila de %s" % a.csv)
        return 1
    print("filas con resultado: %d" % len(filas))

    for diana in sorted({f["target"] for f in filas}):
        grupo = [f for f in filas if f["target"] == diana]
        els = [f["e_l"] for f in grupo]
        dgs = [f["dg"] for f in grupo]
        print()
        print("=" * 72)
        print("%s  (%d moleculas)" % (diana, len(grupo)))
        print("=" * 72)
        print("  e_ligand (energia interna del ligando solo):")
        print("     minimo  %+9.2f" % min(els))
        print("     mediana %+9.2f" % mediana(els))
        print("     maximo  %+9.2f" % max(els))
        print("     rango   %9.2f kcal/mol  <-- si es grande, el dG no es "
              "comparable entre moleculas" % (max(els) - min(els)))
        altos = sorted(grupo, key=lambda f: -f["e_l"])[:5]
        print("  los 5 ligandos mas tensos (e_ligand mas alto):")
        for f in altos:
            print("     %-22s e_lig=%+8.2f  dG=%+8.2f  vina=%+.1f"
                  % (f["ligand"][:22], f["e_l"], f["dg"], f["vina"]))

        # Si el dG solo siguiera al ligando, esta correlacion seria fuerte.
        print("  correlacion (rangos) e_ligand vs dG : %+.3f" % spearman(els, dgs))
        print("  correlacion (rangos) e_ligand vs vina: %+.3f"
              % spearman(els, [f["vina"] for f in grupo]))

    # ¿Cambia el podio si se le quita al dG lo que es tension del ligando?
    print()
    print("=" * 72)
    print("EL CASO CHEMBL520254")
    print("=" * 72)
    todos = {f["ligand"]: f for f in filas if f["target"] == "SOD1"}
    if "CHEMBL520254" in todos:
        f = todos["CHEMBL520254"]
        print("  dG          %+9.2f" % f["dg"])
        print("  e_ligand    %+9.2f   <-- lo que aporta la tension sola" % f["e_l"])
        print("  sin tension %+9.2f   (dG + e_ligand, es decir e_complex - e_receptor)"
              % (f["dg"] + f["e_l"]))
        print("  vina        %+9.2f   (un -6,3 es mediocre, no extraordinario)"
              % f["vina"])

    print()
    print("=" * 72)
    print("PODIO DE SOD1, CON Y SIN LA TENSION DEL LIGANDO")
    print("=" * 72)
    sod1 = [f for f in filas if f["target"] == "SOD1"]
    for f in sod1:
        f["sin_tension"] = f["dg"] + f["e_l"]
        f["atomos"] = atomos(f["smiles"])
    limpio = [f for f in sod1 if f["atomos"]]
    if limpio:
        # Residual = lo que queda tras quitar el sesgo de tamaño.
        xs = [f["atomos"] for f in limpio]
        ys = [f["dg"] for f in limpio]
        pend = _pendiente(xs, ys)
        for f in limpio:
            f["resid"] = f["dg"] - (ys[0] + pend * (f["atomos"] - xs[0]))
            f["resid_st"] = f["sin_tension"] - (ys[0] + pend * (f["atomos"] - xs[0]))
        print("  %-22s %8s %9s %10s %10s" % ("ligando", "atomos", "dG", "residual",
                                             "sin_tension"))
        for f in sorted(limpio, key=lambda g: g["resid"])[:8]:
            print("  %-22s %8d %+9.2f %+10.2f %+10.2f"
                  % (f["ligand"][:22], f["atomos"], f["dg"], f["resid"],
                     f["sin_tension"]))
        print()
        print("  los mismos, ordenados quitando la tension:")
        for f in sorted(limpio, key=lambda g: g["sin_tension"])[:8]:
            print("  %-22s %8d %+9.2f %+10.2f %+10.2f"
                  % (f["ligand"][:22], f["atomos"], f["dg"], f["resid"],
                     f["sin_tension"]))
    return 0


def _pendiente(xs, ys):
    mx, my = st.mean(xs), st.mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    return num / den if den else 0.0


if __name__ == "__main__":
    raise SystemExit(main())
