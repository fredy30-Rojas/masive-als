#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""GNINA (score_only) sobre las poses de ANDR, dentro de WSL.

REPLICA DE CDK2, CON LA PREGUNTA INVERTIDA
-------------------------------------------
Con CDK2 el CNN tampoco puso el cristal por delante (puestos 29 y 55 de 56, y las
ultimas de 56 en afinidad clasica), y eso se leyo como "la pose es discutible".
Aqui se pregunta al reves, y la pregunta si tiene respuesta buena: el ligando TES
es un esteroide con 0 enlaces rotatorios, UNA sola copia en el cristal (B medio
19,2) y los dos polares emparejados a menos de 3,5 A. Si el CNN pone ese cristal
por delante de las poses del motor, hay una segunda herramienta, INDEPENDIENTE de
Vina, que confirma que la pose cristalografica es la buena. Eso convierte el
control en una validacion y no en una discordancia.

POR QUE EN WSL Y NO EN KAGGLE
-----------------------------
La cuota semanal de GPU de Kaggle esta agotada (30 h, gastadas la noche del 28 al
29). Se usa el Ubuntu que ya estaba instalado, con el binario estatico de GNINA
1.3.3 y las librerias CUDA bajadas de PyPI a mano (`_wsl_bajar_libs.py`; ese Ubuntu
no trae pip ni ensurepip). Todo en `~/gnina`, se borra con `rm -rf ~/gnina`.

QUE RESPONDE
------------
Donde queda la pose del cristal entre todas las puntuadas, en los tres criterios
(afinidad clasica, CNNscore y CNNaffinity). El cristal esta en cabeza de las tres
si el rescorig con CNN tiene sentido en esta diana.

OJO SI SOLO HAY UNA POSE
------------------------
Si el control de `andr` aun no ha corrido, en la carpeta solo esta el cristal. En
ese caso GNINA no puede decir si el CNN lo prefiere, porque no hay con que
compararlo: puntuacion alta de una sola pose no es segunda opinion, es un numero
suelto. El script lo avisa en vez de dar un veredicto que no se sostiene.

Uso (dentro de WSL):
    python3 /mnt/c/Users/Fredy/masive-als/analysis/calibracion_dude_cdk2/_wsl_puntuar_gnina_andr.py
"""
import csv
import glob
import os
import re
import subprocess
import sys
import time

BASE = ("/mnt/c/Users/Fredy/masive-als/analysis/calibracion_dude_cdk2/"
        "_gnina_poses_andr")
GNINA = os.path.expanduser("~/gnina/gnina")
REC = os.path.join(BASE, "receptor_andr.pdbqt")
SALIDA = os.path.expanduser("~/gnina/gnina_andr.csv")
SALIDA_WINDOWS = ("/mnt/c/Users/Fredy/masive-als/analysis/"
                  "calibracion_dude_cdk2/gnina_andr.csv")

LISTON = 2.0   # el mismo liston del proyecto, para la coherencia de las cifras


def preparar_librerias():
    """LD_LIBRARY_PATH con las librerias CUDA de ~/gnina/libs.

    Se pone aqui y no en el shell: los subprocess heredan el entorno y asi no hace
    falta exportar nada fuera, que es el enredo de comillas que fallo la primera vez.
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
    v = subprocess.run([GNINA, "--version"], capture_output=True, text=True)
    print("gnina:", (v.stdout or "").strip().splitlines()[-1] if v.stdout else "?")

    poses = sorted(glob.glob(os.path.join(BASE, "out", "*.pdbqt")))
    print("receptor: %s | poses: %d" % (os.path.basename(REC), len(poses)))
    if not poses:
        raise SystemExit("no hay poses en %s/out" % BASE)

    a, s, ca, t = puntuar(poses[0])
    open(os.path.expanduser("~/gnina/ejemplo_andr.txt"), "w").write(t)
    print("")
    print("=== salida cruda de la primera pose ===")
    print(t[-700:])
    print("=== afinidad=%s cnnscore=%s cnnaffinity=%s ===" % (a, s, ca))

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
            escribir(filas)

    escribir(filas)
    print("")
    print("csv escrito: gnina_andr.csv (%d filas)" % len(filas))

    # --- donde queda el cristal ---
    n_cristal = sum(1 for f in filas if f["fichero"].startswith("cristal_"))
    if n_cristal == 0:
        print("\nNO HAY NINGUNA POSE DEL CRISTAL: no se puede concluir nada.")
        return 1
    if len(filas) < 5:
        print("\nSOLO HAY %d POSES. El control de andr aun no ha corrido, asi que"
              % len(filas))
        print("no hay poses del motor con las que comparar. Una puntuacion alta de")
        print("una sola pose NO es una segunda opinion: es un numero suelto. Hay que")
        print("volver a ejecutar esto cuando existan las carpetas _control_andr/*.")
        return 2

    print("")
    print("DONDE QUEDA EL CRISTAL DE ENTRE LAS %d POSES" % len(filas))
    veredictos = {}
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
        for k, v in cristal:
            veredictos[campo] = v
    return 0 if veredictos else 1


def escribir(filas):
    """El CSV en WSL y tambien en Windows, para poder leerlo sin entrar en WSL."""
    for ruta in (SALIDA, SALIDA_WINDOWS):
        try:
            with open(ruta, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=["fichero", "afinidad",
                                                  "cnnscore", "cnnaffinity"])
                w.writeheader()
                w.writerows(filas)
        except OSError:
            pass


if __name__ == "__main__":
    sys.exit(main())
