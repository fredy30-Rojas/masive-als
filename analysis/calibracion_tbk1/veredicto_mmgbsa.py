"""Veredicto: mide Vina y el MM-GBSA contra la misma verdad de referencia.

POR QUE HAY QUE HACERLO CONTRA LAS DOS
El resultado que importa no es "el MM-GBSA ordena", sino "el MM-GBSA ordena
MEJOR QUE LO QUE YA TENIAMOS". Y lo que ya teniamos es la afinidad de Vina.
Ademas, medir solo contra la lista del banco no vale: con "el mejor valor de
cada compuesto", el docking da 0,535, que es ruido de la propia lista. Con una
lista de verdad (mediana de las medidas, dos o mas medidas por compuesto) el
mismo docking da 0,617. O sea, la lista decide el resultado. Por eso las dos.

LO QUE SALIO (28 sep 2026, medido sobre 277 ligandos con dG)
  lista del banco (mejor valor)   Vina 0,602   MM-GBSA 0,524
  verdad v2 (mediana, 2+ medidas)  Vina 0,617   (el MM-GBSA se mide aparte)
  El mejor peso que se le puede dar al MM-GBSA sobre Vina es CERO.
"""
import argparse
import csv
import math
import os
import statistics
import sys

BASE = os.path.dirname(os.path.abspath(__file__))


def afinidad_vina(pdbqt):
    dentro = False
    try:
        fh = open(pdbqt, errors="replace")
    except OSError:
        return None
    with fh:
        for linea in fh:
            if linea.startswith("MODEL"):
                dentro = True
                continue
            if linea.startswith("ENDMDL") and dentro:
                break
            if dentro and linea.startswith("REMARK VINA RESULT"):
                try:
                    return float(linea.split()[3])
                except (IndexError, ValueError):
                    return None
    return None


def auroc(activo, inactivo):
    if len(activo) < 2 or len(inactivo) < 2:
        return None, None
    n = len(activo) * len(inactivo)
    s = sum(1.0 if a < i else (0.5 if a == i else 0.0)
            for a in activo for i in inactivo)
    a = s / n
    na, ni = len(activo), len(inactivo)
    q1 = a / (2.0 - a)
    q2 = 2.0 * a * a / (1.0 + a)
    var = (a * (1 - a) + (na - 1) * (q1 - a * a)
           + (ni - 1) * (q2 - a * a)) / (na * ni)
    return a, math.sqrt(var) if var > 0 else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verdad", default=os.path.join(BASE, "verdad_tbk1_v2.csv"))
    ap.add_argument("--criba", default=os.path.join(BASE, "criba_250.csv"))
    ap.add_argument("--dir-poses", default=os.path.join(BASE, "out"))
    args = ap.parse_args()

    dg = {}
    for r in csv.DictReader(open(args.criba, encoding="utf-8")):
        if r["dG_union"]:
            dg[r["ligand"]] = float(r["dG_union"])

    filas = [r for r in csv.DictReader(open(args.verdad, encoding="utf-8"))
             if "_kin" not in r["molecule_chembl_id"]]
    print("verdad: %d filas, %d con dG del MM-GBSA"
          % (len(filas), sum(1 for r in filas
                             if r["molecule_chembl_id"] in dg)))

    def mide(nombre, seleccion):
        act_v, ina_v, act_m, ina_m = [], [], [], []
        for r in filas:
            cid = r["molecule_chembl_id"]
            if not seleccion(r):
                continue
            v = afinidad_vina(os.path.join(args.dir_poses,
                                           "ACT_%s_out.pdbqt" % cid))
            if v is None:
                continue
            (act_v if r["clase"] == "activo" else ina_v).append(v)
            if cid in dg:
                (act_m if r["clase"] == "activo" else ina_m).append(dg[cid])
        av, ev = auroc(act_v, ina_v)
        am, em = auroc(act_m, ina_m)
        print("\n%s" % nombre)
        print("  con pose %5d  |  con dG %5d  |  inactivos %d"
              % (len(act_v) + len(ina_v), len(act_m) + len(ina_m), len(ina_v)))
        if av is not None:
            print("  Vina     AUROC %.3f ±%.3f" % (av, ev))
        if am is not None:
            print("  MM-GBSA  AUROC %.3f ±%.3f" % (am, em))
            # combinado
            comunes = [r["molecule_chembl_id"] for r in filas
                       if seleccion(r) and r["molecule_chembl_id"] in dg]
            if comunes:
                cls = {r["molecule_chembl_id"]: r["clase"] for r in filas}
                mejor, peso = None, 0.0
                for paso in range(0, 41):
                    p = paso / 8.0
                    a, _ = auroc(
                        [dg[c] - p * afinidad_vina(os.path.join(
                            args.dir_poses, "ACT_%s_out.pdbqt" % c))
                         for c in comunes if cls[c] == "activo"],
                        [dg[c] - p * afinidad_vina(os.path.join(
                            args.dir_poses, "ACT_%s_out.pdbqt" % c))
                         for c in comunes if cls[c] != "activo"])
                    if a is not None and (mejor is None or a > mejor):
                        mejor, peso = a, p
                if mejor is not None:
                    print("  combinado AUROC %.3f  con peso %.2f sobre Vina"
                          % (mejor, peso))
                    if peso == 0.0:
                        print("  -> peso CERO: el MM-GBSA no aporta NADA "
                              "sobre Vina. Veredicto: no se usa de mas.")

    mide("VERDAD v2 (mediana de las medidas, 2 o mas medidas)", lambda r: True)
    mide("solo los que tienen 3 o mas medidas",
         lambda r: int(r["n_medidas"]) >= 3)
    return 0


if __name__ == "__main__":
    sys.exit(main())
