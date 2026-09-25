#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Rescata los ligandos que OBBuilder no pudo generar en 3D.

Estrategia: generar 3D con RDKit (ETKDG + MMFF94), escribir SDF y convertirlo a
PDBQT con OpenBabel (conservando las coordenadas 3D y anadiendo tipos AD).

Esto aplica SOLO a los fallos 'no_genera_3d' del conversor principal, para no
desviarnos del preparador estandar del pipeline (OpenBabel) en el resto.
"""
import csv
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom

BASE = r"C:\Users\Fredy\masive-als"
CLASIF = os.path.join(BASE, r"analysis\ligandos_fallidos_clasificados.csv")
OUT = os.path.join(BASE, r"analysis\_redock_corregidos\pdbqt")
FALLOS = ["CHEMBL113882", "CHEMBL1186247", "CHEMBL1186304", "CHEMBL1186306",
          "CHEMBL1186312", "CHEMBL1186318", "CHEMBL1186333", "CHEMBL1186344",
          "CHEMBL1186345", "CHEMBL1186395", "CHEMBL1187666", "CHEMBL1187708",
          "CHEMBL1187732", "CHEMBL1189155", "CHEMBL1196224"]
MIN_SIZE = 300


def rescatar(nombre, smiles):
    out_p = os.path.join(OUT, nombre + ".pdbqt")
    if os.path.exists(out_p) and os.path.getsize(out_p) > MIN_SIZE:
        return nombre, True, "ya_existia"
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return nombre, False, "rdkit_no_parsea"
        mol = Chem.AddHs(mol)
        embed_ok = AllChem.EmbedMolecule(mol, randomSeed=42) == 0
        if not embed_ok:
            # fallback: ETKDGv3 con torsiones de anillos y sin quiralidad forzada
            ps = rdDistGeom.ETKDGv3()
            ps.randomSeed = 2026
            ps.useSmallRingTorsions = True
            ps.useMacrocycleTorsions = True
            ps.enforceChirality = False
            embed_ok = AllChem.EmbedMolecule(mol, ps) == 0
        if not embed_ok:
            return nombre, False, "rdkit_no_embed"
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except Exception:
            pass
        sdf_tmp = os.path.join(OUT, "_tmp_%s.sdf" % nombre)
        w = Chem.SDWriter(sdf_tmp)
        w.write(mol)
        w.close()

        # leer el SDF como string y borrarlo ANTES de que OpenBabel lo abra
        # (evita WinError 32 por manejador de archivo retenido en Windows)
        with open(sdf_tmp, "r", encoding="utf-8", errors="replace") as fh:
            sdf_text = fh.read()
        os.remove(sdf_tmp)

        from openbabel import openbabel as ob
        conv = ob.OBConversion()
        conv.SetInFormat("sdf")
        omol = ob.OBMol()
        if not conv.ReadString(omol, sdf_text) or not omol.Has3D():
            return nombre, False, "ob_no_lee_sdf"
        conv.SetOutFormat("pdbqt")
        out_s = conv.WriteString(omol)
        if not out_s or len(out_s) < MIN_SIZE:
            return nombre, False, "pdbqt_vacio"
        with open(out_p, "w", encoding="utf-8") as f:
            f.write(out_s)
        if os.path.getsize(out_p) <= MIN_SIZE:
            os.remove(out_p)
            return nombre, False, "pdbqt_demasiado_pequeno"
        return nombre, True, "ok_rdkit"
    except Exception as exc:
        return nombre, False, "excepcion: %s" % str(exc)[:120]


def main():
    t0 = time.time()
    os.makedirs(OUT, exist_ok=True)
    with open(CLASIF, newline="", encoding="utf-8", errors="replace") as h:
        rows = list(csv.DictReader(h))
    por_lig = {r["ligand"]: r for r in rows}
    items = []
    for nombre in FALLOS:
        r = por_lig.get(nombre)
        if not r:
            print("sin registro:", nombre)
            continue
        smi = (r.get("smiles_corregido") or "").strip()
        if not smi:
            smi = (r.get("smiles_libreria") or "").strip()
        items.append((nombre, smi))

    ok = fallos = 0
    with ProcessPoolExecutor(max_workers=8) as ex:
        futs = [ex.submit(rescatar, n, s) for n, s in items]
        for fut in as_completed(futs):
            nombre, res, detalle = fut.result()
            if res:
                ok += 1
                print("OK", nombre)
            else:
                fallos += 1
                print("FALLO", nombre, "|", detalle)
    print("rescatados_ok=%d fallos=%d tiempo=%.0fs" % (ok, fallos, time.time() - t0))


if __name__ == "__main__":
    main()