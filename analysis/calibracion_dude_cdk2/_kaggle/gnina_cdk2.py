"""GNINA (score_only) sobre las poses de CDK2: distingue el CNN la pose del cristal?

LA PREGUNTA
-----------
Vina-GPU no reproduce la pose del cristal de CDK2: el nucleo anclado si lo acierta
en 2 de 6 variantes (1,95 y 2,00 A) pero el brazo del amonio se va de 6 a 9 A. Y el
brazo resulta que en el cristal esta suelto (ni proteina a menos de 4 A ni agua a
menos de 4,5 A) y modelado en dos conformaciones al 50 %.

Si el CNN de GNINA, que esta entrenado con poses cristalograficas, pone las dos
conformaciones del cristal por delante de las 54 poses del motor, entonces hay una
herramienta que si distingue el modo de union nativo en esta diana, y el rescoring
tiene sentido. Si las pone al final, no lo hay.

DE DONDE SALE ESTO
------------------
Es el kernel de TBK1 (`yograbotodo/masive-als-gnina-score`) adaptado: mismo binario
estatico, mismas librerias CUDA por pip, y el mismo cuidado con las lineas
MODEL/ENDMDL. Aqui las poses ya vienen de una en una, asi que no hace falta partirlas.
"""
import csv
import glob
import os
import re
import subprocess

SALIDA = "/kaggle/working/gnina_cdk2.csv"

# --- 1. el binario y sus librerias (lo mismo que funciono en TBK1) ---
print("cuda:", end=" ")
try:
    import torch
    print(torch.cuda.is_available())
except Exception as e:  # noqa: BLE001
    print("no hay torch:", e)

subprocess.run("wget -q -O /tmp/gnina https://github.com/gnina/gnina/releases/"
               "download/v1.3.3/gnina.cuda12.8.static", shell=True, check=True)
subprocess.run("chmod +x /tmp/gnina", shell=True, check=True)
subprocess.run("pip -q install nvidia-cudnn-cu12 nvidia-cuda-runtime-cu12 "
               "nvidia-cublas-cu12 nvidia-cufft-cu12 nvidia-curand-cu12 "
               "nvidia-cusparse-cu12 nvidia-cusolver-cu12 nvidia-cuda-nvrtc-cu12",
               shell=True)
nvdirs = (glob.glob("/usr/local/lib/python*/dist-packages/nvidia/*/lib")
          + glob.glob("/opt/conda/lib/python*/site-packages/nvidia/*/lib")
          + glob.glob("/usr/lib/python3*/dist-packages/nvidia/*/lib"))
assert nvdirs, "no encontre las librerias nvidia del pip"
os.environ["LD_LIBRARY_PATH"] = (":".join(nvdirs) + ":"
                                 + os.environ.get("LD_LIBRARY_PATH", ""))
print(subprocess.run(["/tmp/gnina", "--version"], capture_output=True,
                     text=True).stdout[-200:])

# --- 2. los datos ---
rec = (glob.glob("/kaggle/input/*/receptor_cdk2.pdbqt")
       or glob.glob("/kaggle/input/**/receptor_cdk2.pdbqt", recursive=True))[0]
# Se buscan las poses EN TODO el dataset y no en una carpeta concreta: al subir un
# dataset, el cliente de Kaggle se salta las carpetas salvo que se le pida
# `--dir-mode`, y de esta forma da igual si acaban en `out/`, en la raiz o dentro de
# un zip que Kaggle descomprime solo.
poses = sorted(p for p in glob.glob("/kaggle/input/**/*.pdbqt", recursive=True)
               if os.path.basename(p) != "receptor_cdk2.pdbqt")
print("receptor:", rec, "| poses:", len(poses))
assert poses, "no hay poses en %s" % dir_poses

GNINA = "/tmp/gnina"


def puntuar(pose):
    """(afinidad, cnnscore, cnnaffinity, salida cruda) de una pose, sin acoplar."""
    r = subprocess.run([GNINA, "-r", rec, "-l", pose, "--score_only"],
                       capture_output=True, text=True, timeout=900)
    t = r.stdout + r.stderr
    # GNINA saca varias cosas; se recogen por nombre y no por posicion.
    def primero(patron, casteo=float):
        m = re.search(patron, t)
        return casteo(m.group(1)) if m else None
    return (primero(r"Affinity:\s*(-?[0-9]+\.[0-9]+)"),
            primero(r"CNNscore:\s*([0-9]+\.[0-9]+)"),
            primero(r"CNNaffinity:\s*(-?[0-9]+\.[0-9]+)"),
            t)


# diagnostico: la primera pose, cruda, para ver el formato exacto de salida
a, s, ca, t = puntuar(poses[0])
open("/kaggle/working/ejemplo_salida.txt", "w").write(t)
print("CRUDA:\n", t[:700])
print("valores:", a, s, ca)

filas = []
for i, p in enumerate(poses, 1):
    nombre = os.path.basename(p)[:-len(".pdbqt")]
    try:
        a, s, ca, _ = puntuar(p)
    except Exception as e:  # noqa: BLE001
        a = s = ca = None
        print("FALLO", nombre, repr(e)[:90])
    filas.append({"fichero": nombre, "afinidad": a, "cnnscore": s,
                  "cnnaffinity": ca})
    if i % 10 == 0:
        print("%d/%d" % (i, len(poses)))
        with open(SALIDA, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["fichero", "afinidad", "cnnscore",
                                              "cnnaffinity"])
            w.writeheader()
            w.writerows(filas)

with open(SALIDA, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["fichero", "afinidad", "cnnscore",
                                      "cnnaffinity"])
    w.writeheader()
    w.writerows(filas)
print("csv escrito:", SALIDA, len(filas), "filas")

# --- 3. donde queda el cristal ---
print("")
print("DONDE QUEDA EL CRISTAL")
for campo, mayor_es_mejor in (("cnnaffinity", True), ("cnnscore", True),
                              ("afinidad", False)):
    validas = [f for f in filas if f[campo] is not None]
    if len(validas) < 5:
        print("   %-12s sin datos suficientes (%d)" % (campo, len(validas)))
        continue
    validas.sort(key=lambda f: f[campo], reverse=mayor_es_mejor)
    puestos = {f["fichero"]: i + 1 for i, f in enumerate(validas)}
    cristal = {k: v for k, v in puestos.items() if k.startswith("cristal_")}
    print("   %-12s de %d poses | puesto del cristal: %s"
          % (campo, len(validas),
             " | ".join("%s = %d" % kv for kv in sorted(cristal.items()))))
    print("      las 3 primeras: %s"
          % ", ".join("%s (%.3f)" % (f["fichero"], f[campo]) for f in validas[:3]))

with open("/kaggle/working/veredicto.txt", "w") as f:
    f.write("listo\n")
print("hecho")
