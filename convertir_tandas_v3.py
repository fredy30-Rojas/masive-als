#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convierte seleccion_30000.tsv -> 100 tandas (l300..l399) de PDBQT 3D.
Corre en la PC (Windows) con Open Babel 3.1.0. Resumible por tanda empaquetada.
Uso: python convertir_tandas_v3.py
"""
import math, os, sys, tarfile, time
from concurrent.futures import ProcessPoolExecutor, as_completed

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
os.environ.setdefault("BABEL_DATADIR", r"C:\Program Files\OpenBabel-3.1.1\data")

BASE = r"C:\Users\Fredy\masive-als\compounds"
TSV = os.path.join(BASE, "seleccion_30000_v2.tsv")
OUT = os.path.join(BASE, "lotes_v4_pdbqt")
os.makedirs(OUT, exist_ok=True)

N_TANDAS = 100
PER_TANDA = 300
START = 400          # l400 .. l499
WORKERS = 10
LOG = os.path.join(BASE, "fabrica_v4_log.txt")


def log(msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    print(line, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def convert_batch(args):
    outdir, items = args
    from openbabel import openbabel as ob
    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    conv.SetOutFormat("pdbqt")
    ok = 0
    for cid, smi in items:
        outp = os.path.join(outdir, cid + ".pdbqt")
        if os.path.exists(outp) and os.path.getsize(outp) > 100:
            ok += 1
            continue
        mol = ob.OBMol()
        try:
            if not conv.ReadString(mol, smi):
                continue
            mol.AddHydrogens()
            b = ob.OBBuilder()
            if not b.Build(mol):
                continue
            ff = ob.OBForceField.FindForceField("UFF")
            if ff:
                ff.Setup(mol)
                ff.SteepestDescent(150)
                ff.ConjugateGradients(50)
                ff.GetCoordinates(mol)
            out_s = conv.WriteString(mol)
            if out_s and len(out_s) > 100:
                with open(outp, "w", encoding="utf-8") as f:
                    f.write(out_s)
                ok += 1
        except Exception:
            continue
    return ok


def main():
    t0 = time.time()
    items = []
    with open(TSV, encoding="utf-8") as f:
        for line in f:
            p = line.rstrip("\n").split("\t", 1)
            if len(p) == 2 and p[0].strip() and p[1].strip():
                items.append((p[0].strip(), p[1].strip()))
    log("compuestos leidos: %d" % len(items))

    tandas = []
    for i in range(N_TANDAS):
        chunk = items[i * PER_TANDA:(i + 1) * PER_TANDA]
        if not chunk:
            break
        tandas.append(("l%03d" % (START + i), chunk))
    log("tandas: %d" % len(tandas))

    manifest_new = []
    for tanda_id, chunk in tandas:
        tgz = os.path.join(BASE, "lote_%s_plano.tar.gz" % tanda_id)
        if os.path.exists(tgz) and os.path.getsize(tgz) > 1000:
            manifest_new.append("lote_%s_plano.tar.gz" % tanda_id)
            log("%s: ya empaquetada" % tanda_id)
            continue

        outdir = os.path.join(OUT, tanda_id)
        os.makedirs(outdir, exist_ok=True)

        n = len(chunk)
        csize = max(1, math.ceil(n / WORKERS))
        subchunks = [chunk[j:j + csize] for j in range(0, n, csize)]

        total_ok = 0
        with ProcessPoolExecutor(max_workers=WORKERS) as ex:
            futs = [ex.submit(convert_batch, (outdir, sc)) for sc in subchunks]
            for fut in as_completed(futs):
                total_ok += fut.result()

        pdbs = sorted(f for f in os.listdir(outdir) if f.endswith(".pdbqt"))
        if not pdbs:
            log("%s: 0 pdbqt, no empaqueto" % tanda_id)
            continue
        with tarfile.open(tgz, "w:gz") as t:
            for f in pdbs:
                t.add(os.path.join(outdir, f), arcname=f)
        log("%s: %d pdbqt -> lote_%s_plano.tar.gz (%d KB)"
            % (tanda_id, len(pdbs), tanda_id, os.path.getsize(tgz) // 1024))
        manifest_new.append("lote_%s_plano.tar.gz" % tanda_id)

    log("DONE: %d tandas en %.0f min" % (len(manifest_new), (time.time() - t0) / 60))


if __name__ == "__main__":
    main()
