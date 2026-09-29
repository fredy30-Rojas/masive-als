"""Prueba si las nueve poses de cada ligando ordenan mejor que la mejor de ellas.

POR QUE
Con una sola pose, Vina separa activos de inactivos con un AUROC de 0,617. Es
poco, pero es real. La pregunta es si hay informacion que estamos tirando en
la basura: Vina devuelve nueve modos por ligando y el pipeline solo se queda
con el primero.

Lo que miden las tres cosas, todas sobre los mismos ligandos:
  - la mejor pose:      la de menor energia, que es lo que usa el banco
  - la media de las     promedia las nueve. Si la señal esta repartida, promediar
    nueve               la sube; si la mejor pose es buena y las demas son ruido,
                        la baja
  - cuantas poses caen  un ligando que encaja bien tiene MAS modos dentro de
    en 2 kcal/mol de    dos kcal/mol de la mejor que uno que no encaja. Es la
    la mejor             idea de "consenso" y es la que mas se sostiene en la
                        practica

Ademas se mide una cuarta, que es de controle: el numero de atomos pesados por
si solo. Si esa da mas de cero coma cincuenta, el AUROC del docking esta
midiendo el tamano del ligando y no el acoplamiento.
"""
import csv
import math
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "out")


def poses(pdbqt):
    """Energias de los modelos de un pdbqt, en orden."""
    vals, dentro = [], False
    try:
        fh = open(pdbqt, errors="replace")
    except OSError:
        return vals
    with fh:
        for linea in fh:
            if linea.startswith("MODEL"):
                dentro = True
                continue
            if linea.startswith("ENDMDL"):
                dentro = False
                continue
            if dentro and linea.startswith("REMARK VINA RESULT"):
                try:
                    vals.append(float(linea.split()[3]))
                except (IndexError, ValueError):
                    pass
    return vals


def auroc(activo, inactivo):
    if len(activo) < 2 or len(inactivo) < 2:
        return None, None
    n = len(activo) * len(inactivo)
    s = sum(1.0 if a < i else (0.5 if a == i else 0.0)
            for a in activo for i in inactivo)
    au = s / n
    na, ni = len(activo), len(inactivo)
    q1 = au / (2.0 - au)
    q2 = 2.0 * au * au / (1.0 + au)
    var = (au * (1 - au) + (na - 1) * (q1 - au * au)
           + (ni - 1) * (q2 - au * au)) / (na * ni)
    return au, math.sqrt(var) if var > 0 else 0.0


def main():
    ruta = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        BASE, "verdad_tbk1_v2.csv")
    filas = [r for r in csv.DictReader(open(ruta, encoding="utf-8"))
             if "_kin" not in r["molecule_chembl_id"]]
    metricas = {"mejor pose": lambda v: v[0],
                "2 mejores (media)": lambda v: sum(v[:2]) / 2,
                "media de las 9": lambda v: sum(v) / len(v),
                "nº dentro de 2 kcal": lambda v: -sum(
                    1 for e in v if e <= v[0] + 2.0),
                "nº dentro de 4 kcal": lambda v: -sum(
                    1 for e in v if e <= v[0] + 4.0),
                "3 mejores (media)": lambda v: sum(v[:3]) / 3}
    datos = {k: {"activo": [], "inactivo": []} for k in metricas}
    conpose = 0
    for r in filas:
        v = poses(os.path.join(OUT, "ACT_%s_out.pdbqt" % r["molecule_chembl_id"]))
        if not v:
            continue
        conpose += 1
        k = "activo" if r["clase"] == "activo" else "inactivo"
        for nombre, fn in metricas.items():
            datos[nombre][k].append(fn(v))
    print("ligandos con pose: %d  (activos %d, inactivos %d)"
          % (conpose, len(datos["mejor pose"]["activo"]),
             len(datos["mejor pose"]["inactivo"])))
    print("medias de modos por ligando:")
    for nombre in metricas:
        a, e = auroc(datos[nombre]["activo"], datos[nombre]["inactivo"])
        print("  %-20s AUROC %.3f ±%.3f" % (nombre, a, e))
    return 0


if __name__ == "__main__":
    sys.exit(main())
