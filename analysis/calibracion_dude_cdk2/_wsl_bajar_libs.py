#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Baja las librerias CUDA que pide GNINA, sin pip y sin apt.

POR QUE ASI
-----------
El binario estatico de GNINA 1.3.3 pide `libcudnn.so.9` y el Ubuntu de WSL no la
tiene. La forma normal es `pip install nvidia-cudnn-cu12`, pero este Ubuntu no trae
pip ni ensurepip, y meterlos pide `apt`. Asi que se baja el fichero .whl de PyPI
directamente (un .whl no es mas que un zip) y se descomprime con el propio python3,
que si esta. Todo queda en `~/gnina/libs`, que se borra con un `rm -rf` y no toca
nada del sistema.

Uso (dentro de WSL):
    python3 _wsl_bajar_libs.py nvidia-cudnn-cu12
"""
import json
import os
import sys
import urllib.request
import zipfile

DEST = os.path.expanduser("~/gnina/libs")


def elegir_wheel(pak):
    """El .whl mas nuevo para linux x86_64, prefiriendo el manylinux mas moderno."""
    d = json.load(urllib.request.urlopen(
        "https://pypi.org/pypi/%s/json" % pak, timeout=120))
    def clave(f):
        n = f["filename"]
        if not n.endswith(".whl") or "x86_64" not in n or "manylinux" not in n:
            return None
        if "manylinux_2_28" in n:
            return 4
        if "manylinux_2_27" in n:
            return 3
        if "manylinux_2_24" in n:
            return 2
        if "manylinux2014" in n:
            return 1
        return 0
    cand = [(clave(f), f) for f in d["urls"]]
    cand = [c for c in cand if c[0] is not None]
    if not cand:
        raise SystemExit("no hay wheel de linux x86_64 para %s" % pak)
    cand.sort(key=lambda c: c[0], reverse=True)
    return cand[0][1]


def main():
    paks = sys.argv[1:] or ["nvidia-cudnn-cu12"]
    os.makedirs(DEST, exist_ok=True)
    for pak in paks:
        f = elegir_wheel(pak)
        ruta = os.path.join(DEST, f["filename"])
        if os.path.exists(ruta):
            print("ya estaba:", f["filename"])
        else:
            print("bajando %s (%.0f MB)" % (f["filename"], f["size"] / 1e6),
                  flush=True)
            urllib.request.urlretrieve(f["url"], ruta)
        print("descomprimiendo", f["filename"], flush=True)
        with zipfile.ZipFile(ruta) as z:
            z.extractall(DEST)
    # Las librerias quedan en ~/gnina/libs/nvidia/<paquete>/lib
    libs = []
    for raiz, _, ficheros in os.walk(DEST):
        if os.path.basename(raiz) != "lib":
            continue
        if any(x.endswith(".so.9") or ".so." in x for x in ficheros):
            libs.append(raiz)
    print("")
    print("rutas con librerias:")
    for l in sorted(set(libs)):
        print("   ", l)
    print("")
    print("LD_LIBRARY_PATH=%s" % ":".join(sorted(set(libs))))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
