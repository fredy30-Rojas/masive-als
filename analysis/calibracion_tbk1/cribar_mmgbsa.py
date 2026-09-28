"""Criba de MM-GBSA sobre una lista de ligandos, en paralelo. Corre en Oracle.

QUE HACE
Lee `lista.csv` (ligand,clase,pchembl), busca la pose de cada uno en
`--dir-poses` y lanza `rescoring_mmgbsa.py` en un proceso por nucleo.

POR QUE EN PARALELO (medido el 28 sep 2026 en Oracle): un sistema son 12 s, o
unos 40 s por ligando contando los tres (complejo, receptor, ligando). En
serie salen 40 s por ligando; con los 4 nucleos de la maquina, unos 10 s. Con
mil ligandos la diferencia son once horas contra hora y media, y no cuesta
nada: solo abrir un proceso por nucleo.

POR QUE NO KAGGLE: el MM-GBSA con Generalized Born lo AmberTools lo calcula en
el procesador, no en la tarjeta, asi que Kaggle no daria mas potencia, daria
nucleos, que es lo mismo. Y su cuota son 30 horas de GPU por semana que no se
recuperan. Se guarda para cuando se sepa que el numero ordena.
"""
import argparse
import csv
import os
import subprocess
import sys
import time
from multiprocessing import Pool

ESCORING = "/home/ubuntu/prueba_mmgbsa/rescoring_mmgbsa.py"
TRABAJO = "/home/ubuntu/criba_mmgbsa"


def un_ligando(args):
    """Rescoring de UN ligando. Devuelve (ligand, dG, nota)."""
    lig, receptor, poses = args
    pose = os.path.join(poses, "ACT_%s_out.pdbqt" % lig)
    if not os.path.exists(pose):
        return lig, None, "no esta la pose"
    cmd = [sys.executable, ESCORING, "--receptor", receptor, "--pose", pose]
    t0 = time.time()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
    except subprocess.TimeoutExpired:
        return lig, None, "se paso del tiempo"
    dG = None
    motivo = "sin dG (rc=%d)" % r.returncode
    for linea in r.stdout.splitlines():
        if "dG_union MM-GBSA" in linea:
            try:
                dG = float(linea.split("=")[1].split()[0])
            except (IndexError, ValueError):
                pass
        elif "FALLO" in linea or "solapamiento" in linea:
            motivo = linea.strip()[:130]
    if dG is None:
        return lig, None, "%s (%.0f s)" % (motivo, time.time() - t0)
    return lig, dG, "ok (%.0f s)" % (time.time() - t0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lista", required=True, help="CSV con ligand,clase,pchembl")
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--dir-poses", required=True)
    ap.add_argument("--nucleos", type=int, default=4)
    ap.add_argument("--salida", default=os.path.join(TRABAJO, "criba.csv"))
    args = ap.parse_args()

    with open(args.lista, encoding="utf-8") as fh:
        filas = [r for r in csv.DictReader(fh)]
    print("lista: %d ligandos" % len(filas))
    clases = {}
    for r in filas:
        clases[r["clase"]] = clases.get(r["clase"], 0) + 1
    print("por clase:", clases)

    trabajo = [(r["ligand"], args.receptor, args.dir_poses) for r in filas]
    t0 = time.time()
    print("lanzando %d rescores con %d procesos..." % (len(trabajo),
                                                       args.nucleos))
    with Pool(args.nucleos) as pool:
        res = pool.map(un_ligando, trabajo)

    con_dg = [r for r in res if r[1] is not None]
    print("\n%d de %d con dG, en %.0f s (%.1f s por ligando con %d nucleos)"
          % (len(con_dg), len(res), time.time() - t0,
             (time.time() - t0) * args.nucleos / max(1, len(res)),
             args.nucleos))

    info = {r["ligand"]: r for r in filas}
    with open(args.salida, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["ligand", "clase", "pchembl", "dG_union", "nota"])
        for lig, dG, nota in res:
            r = info[lig]
            w.writerow([lig, r["clase"], r["pchembl"],
                        "" if dG is None else "%.2f" % dG, nota])
    print("escrito en %s" % args.salida)

    for lig, dG, nota in res:
        marca = "%8.2f" % dG if dG is not None else "    FALLO"
        print("  %-16s %-9s %s  %s" % (lig, info[lig]["clase"], marca, nota))
    return 0 if len(con_dg) == len(res) else 2


if __name__ == "__main__":
    sys.exit(main())
