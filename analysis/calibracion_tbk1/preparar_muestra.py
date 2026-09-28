"""Prepara la muestra para medir si el MM-GBSA ordena.

Toma los ligandos YA dockeados que ademas tienen actividad medida en ChEMBL, y
los parte en dos clases: activos y de la franja de inactivos. La comparacion
tiene que ser entre las dos clases, porque con una sola el area bajo la curva
no se puede calcular.

Reparto de los activos por franja: 50 de cada una de las tres franjas de
"activo de verdad" (>=10 nM, 10-100 nM y 0,1-1 uM). Repartirlos asi y no
cogerlos al azar es a proposito: si se cogen al azar, la mayoria caen en la
franja mejor y no hay rango de afinidad con el que trabajar.

Escribe:
  muestra_250.csv   ligand,clase,pchembl (aunque salgan 300, el nombre del
                    fichero dice 250: se cambio el reparto de franjas despues
                    de escribirlo y no vale la pena renombrarlo, se lee del
                    csv)
  poses250.tgz      las poses de esa muestra
"""
import csv
import os
import random
import tarfile

BASE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(BASE, "out")
VERDAD = os.path.join(BASE, "compuestos_tbk1.csv")

POR_FRANJA = 50
N_INACTIVOS = 100
SEMILLA = 7


def main():
    ref = {}
    for r in csv.DictReader(open(VERDAD, encoding="utf-8")):
        try:
            ref[r["molecule_chembl_id"]] = (float(r["pchembl"]), r["franja"])
        except (ValueError, KeyError):
            continue
    hechos = {f[:-10][4:] for f in os.listdir(OUT)
              if f.endswith("_out.pdbqt") and f.startswith("ACT_")}

    por_franja = {}
    inactivos = []
    for c in hechos:
        if c not in ref:
            continue
        pch, franja = ref[c]
        if "peor" in franja:
            inactivos.append(c)
        else:
            por_franja.setdefault(franja, []).append(c)

    azar = random.Random(SEMILLA)
    activos = []
    for franja, lista in sorted(por_franja.items()):
        elegidas = azar.sample(lista, min(POR_FRANJA, len(lista)))
        print("  %-32s %5d disponibles, %3d elegidas"
              % (franja, len(lista), len(elegidas)))
        activos += elegidas
    elegidos_inactivos = azar.sample(inactivos,
                                      min(N_INACTIVOS, len(inactivos)))
    print("  inactivos: %5d disponibles, %3d elegidas"
          % (len(inactivos), len(elegidos_inactivos)))

    muestra = ([(c, "activo") for c in activos]
               + [(c, "inactivo") for c in elegidos_inactivos])
    ruta_csv = os.path.join(BASE, "muestra_250.csv")
    with open(ruta_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["ligand", "clase", "pchembl"])
        for c, clase in muestra:
            w.writerow([c, clase, ref[c][0]])

    ruta_tgz = os.path.join(BASE, "poses250.tgz")
    faltan = []
    with tarfile.open(ruta_tgz, "w:gz") as tf:
        for c, _ in muestra:
            p = os.path.join(OUT, "ACT_%s_out.pdbqt" % c)
            if os.path.exists(p):
                tf.add(p, arcname="poses/ACT_%s_out.pdbqt" % c)
            else:
                faltan.append(c)
    print("\nmuestra: %d  poses empaquetadas: %d  sin pose: %d"
          % (len(muestra), len(muestra) - len(faltan), len(faltan)))
    if faltan:
        print("  sin pose:", ", ".join(faltan[:10]))
    print("escrito %s" % ruta_csv)
    print("escrito %s (%.1f MB)" % (ruta_tgz, os.path.getsize(ruta_tgz) / 1e6))


if __name__ == "__main__":
    main()
