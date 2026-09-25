#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regenera PDBQT 3D para los 916 ligandos con SMILES corregido.

Replica EXACTAMENTE el preparador de ligandos del pipeline local
(convertir_tandas_v3.py / convertir_fda_full.py):

  1) Leer SMILES (formato smi) con OpenBabel
  2) AddHydrogens()
  3) OBBuilder.Build(mol)  -> genera coordenadas 3D
  4) UFF: SteepestDescent(150) + ConjugateGradients(50)
  5) Escribir PDBQT (tipos AutoDock de OpenBabel)

Resumible: salta archivos ya generados con tamaño valido.
Salida: analysis/_redock_corregidos/pdbqt/<ligand>.pdbqt
"""
import math
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE = r"C:\Users\Fredy\masive-als"
SMI = os.path.join(BASE, r"analysis\ligandos_a_reacoplar.smi")
OUT = os.path.join(BASE, r"analysis\_redock_corregidos\pdbqt")
LOG = os.path.join(BASE, r"analysis\_redock_corregidos\log.txt")
WORKERS = 10
MIN_SIZE = 300   # un PDBQT con atomos pesa > 300 bytes


def log(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def nombre_archivo(nombre):
    n = (nombre or "").strip()
    n = n.replace("/", "_").replace("\\", "_").replace(" ", "_").replace(":", "_")
    n = "".join(ch for ch in n if ch.isalnum() or ch in "._-")
    return n[:120] or "LIG"


def convert_one(args):
    """Un ligando: devuelve (nombre, ok, detalle)."""
    nombre, smiles = args
    from openbabel import openbabel as ob
    out_p = os.path.join(OUT, nombre + ".pdbqt")
    if os.path.exists(out_p) and os.path.getsize(out_p) > MIN_SIZE:
        return nombre, True, "ya_existia"
    try:
        conv = ob.OBConversion()
        conv.SetInFormat("smi")
        mol = ob.OBMol()
        if not conv.ReadString(mol, smiles):
            return nombre, False, "no_lee_smiles"
        mol.AddHydrogens()
        b = ob.OBBuilder()
        if not b.Build(mol):
            return nombre, False, "no_genera_3d"
        ff = ob.OBForceField.FindForceField("UFF")
        if ff:
            ff.Setup(mol)
            ff.SteepestDescent(150)
            ff.ConjugateGradients(50)
            ff.GetCoordinates(mol)
        conv.SetOutFormat("pdbqt")
        out_s = conv.WriteString(mol)
        if not out_s or len(out_s) < MIN_SIZE:
            return nombre, False, "pdbqt_vacio"
        with open(out_p, "w", encoding="utf-8") as f:
            f.write(out_s)
        if os.path.getsize(out_p) <= MIN_SIZE:
            os.remove(out_p)
            return nombre, False, "pdbqt_demasiado_pequeno"
        return nombre, True, "ok"
    except Exception as exc:
        return nombre, False, "excepcion: %s" % str(exc)[:100]


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    items = []
    with open(SMI, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2 and p[0].strip() and p[1].strip():
                items.append((nombre_archivo(p[1]), p[0].strip()))
    log("SMILES a convertir: %d" % len(items))

    ok = fallos = 0
    fallos_lista = []
    n = len(items)
    csize = max(1, math.ceil(n / WORKERS))
    subchunks = [items[j:j + csize] for j in range(0, n, csize)]

    with ProcessPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(convert_one, sc_item) for sc in subchunks for sc_item in sc]
        done = 0
        for fut in as_completed(futs):
            nombre, res, detalle = fut.result()
            done += 1
            if res:
                ok += 1
            else:
                fallos += 1
                fallos_lista.append((nombre, detalle))
            if done % 100 == 0:
                log("progreso %d/%d ok=%d fallos=%d" % (done, n, ok, fallos))

    log("FIN ok=%d fallos=%d tiempo=%.0fs" % (ok, fallos, time.time() - t0))
    if fallos_lista:
        log("--- fallos (%d) ---" % len(fallos_lista))
        for nombre, detalle in fallos_lista[:40]:
            log("  FALLO %s | %s" % (nombre, detalle))
    print("OUT:", OUT)


if __name__ == "__main__":
    main()