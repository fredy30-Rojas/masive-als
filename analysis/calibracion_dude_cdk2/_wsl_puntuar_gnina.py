#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GNINA (score_only) sobre las poses de CDK2, dentro de WSL.

POR QUE EN WSL Y NO EN KAGGLE
-----------------------------
La ruta de Kaggle era la buena (es la que corrio el CNN de TBK1) pero la cuota
semanal de GPU de la cuenta esta agotada: 30 horas, gastadas en la corrida de anoche.
Asi que se hace aqui, en el Ubuntu que ya estaba instalado en WSL, con el binario
estatico de GNINA 1.3.3 y las librerias CUDA bajadas de PyPI a mano (ver
`_wsl_bajar_libs.py`; el Ubuntu no trae pip ni ensurepip).

QUE RESPONDE
------------
Si el CNN de GNINA, que esta entrenado con poses cristalograficas, pone las DOS
conformaciones del cristal de 1h00 por delante de las 54 poses que devolvio Vina-GPU,
entonces hay una herramienta que si reconoce el modo de union nativo en esta diana y
el rescoring tiene sentido. Si las deja al final, no lo tiene.

Uso (dentro de WSL):
    python3 /mnt/c/Users/Fredy/masive-als/analysis/calibracion_dude_cdk2/_wsl_puntuar_gnina.py
"""
import csv
import glob
import os
import re
import subprocess
import sys
import time

BASE = ("/mnt/c/Users/Fredy/masive-als/analysis/calibracion_dude_cdk2/"
        "_gnina_poses")
GNINA = os.path.expanduser("~/gnina/gnina")
REC = os.path.join(BASE, "receptor_cdk2.pdbqt")
SALIDA = os.path.expanduser("~/gnina/gnina_cdk2.csv")


def preparar_librerias():
    """Pone LD_LIBRARY_PATH con las librerias CUDA bajadas a ~/gnina/libs.

    Se hace aqui y no en el shell a proposito: los subprocess heredan el entorno,
    asi que con esto ya no hace falta exportar nada fuera, y de paso se evita el
    enredo de comillas que fallo la primera vez.
    """
    libs = sorted(glob.glob(os.path.expanduser("~/gnina/libs/nvidia/*/lib")))
    if not libs:
        raise SystemExit("no hay librerias en ~/gnina/libs: ejecuta antes "
                         "_wsl_bajar_libs.py")
    os.environ["LD_LIBRARY_PATH"] = (":".join(libs) + ":"
                                     + os.environ.get("LD_LIBRARY_PATH", ""))
    return libs


def puntuar(pose):
    """(afinidad, cnnscore, cnnaffinity, salida cruda) de una pose, sin acoplar."""
    r = subprocess.run([GNINA, "-r", REC, "-l", pose, "--score_only"],
                       capture_output=True, text=True, timeout=1800)
    t = r.stdout + r.stderr
    def primero(patron):
        m = re.search(patron, t)
        return float(m.group(1)) if m else None
    return (primero(r"Affinity:\s*(-?[0-9]+\.[0-9]+)"),
            primero(r"CNNscore:\s*([0-9]+\.[0-9]+)"),
            primero(r"CNNaffinity:\s*(-?[0-9]+\.[0-9]+)"),
            t)


def main():
    libs = preparar_librerias()
    print("librerias: %d rutas" % len(libs))
    print("gnina:", subprocess.run([GNINA, "--version"], capture_output=True,
                                   text=True).stdout.strip().splitlines()[-1])
    poses = sorted(glob.glob(os.path.join(BASE, "out", "*.pdbqt")))
    print("receptor:", os.path.basename(REC), "| poses:", len(poses))
    if not poses:
        raise SystemExit("no hay poses en %s/out" % BASE)

    # Prueba con la primera: si el formato de salida no fuera el esperado, se ve ya.
    a, s, ca, t = puntuar(poses[0])
    open(os.path.expanduser("~/gnina/ejemplo_salida.txt"), "w").write(t)
    print("")
    print("=== salida cruda de la primera pose ===")
    print(t[-900:])
    print("=== valores: afinidad=%s cnnscore=%s cnnaffinity=%s ===" % (a, s, ca))

    filas = []
    t0 = time.time()
    for i, p in enumerate(poses, 1):
        nombre = os.path.basename(p)[:-len(".pdbqt")]
        try:
            a, s, ca, _ = puntuar(p)
        except Exception as e:  # noqa: BLE001
            a = s = ca = None
            print("FALLO", nombre, repr(e)[:90])
        filas.append({"fichero": nombre, "afinidad": a, "cnnscore": s,
                      "cnnaffinity": ca})
        if i % 5 == 0 or i == len(poses):
            print("%d/%d | %.1f s" % (i, len(poses), time.time() - t0), flush=True)
            with open(SALIDA, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["fichero", "afinidad",
                                                  "cnnscore", "cnnaffinity"])
                w.writeheader()
                w.writerows(filas)

    with open(SALIDA, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["fichero", "afinidad", "cnnscore",
                                          "cnnaffinity"])
        w.writeheader()
        w.writerows(filas)
    print("")
    print("csv escrito: %s (%d filas)" % (SALIDA, len(filas)))

    # --- donde queda el cristal ---
    print("")
    print("DONDE QUEDA EL CRISTAL DE ENTRE LAS %d POSES" % len(filas))
    for campo, mayor_es_mejor in (("cnnaffinity", True), ("cnnscore", True),
                                  ("afinidad", False)):
        validas = [f for f in filas if f[campo] is not None]
        if len(validas) < 5:
            print("   %-12s sin datos suficientes (%d)" % (campo, len(validas)))
            continue
        validas.sort(key=lambda f: f[campo], reverse=mayor_es_mejor)
        puestos = {f["fichero"]: i + 1 for i, f in enumerate(validas)}
        cristal = sorted((k, v) for k, v in puestos.items()
                         if k.startswith("cristal_"))
        print("   %-12s | puesto del cristal: %s"
              % (campo, " | ".join("%s = %d de %d" % (k, v, len(validas))
                                   for k, v in cristal)))
        print("      las 5 primeras: %s"
              % ", ".join("%s (%.3f)" % (f["fichero"], f[campo])
                          for f in validas[:5]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
