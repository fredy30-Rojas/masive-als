#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convierte las tandas t3/t6/t7/t10/t12 de SMILES a PDBQT 3D plano y las
empaqueta como ligandos_batchN_plano.tar.gz (para GitHub / docking)."""
import csv, os, sys, tarfile, time

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
os.environ.setdefault('BABEL_DATADIR', r'C:\Program Files\OpenBabel-3.1.1\data')

from openbabel import openbabel as ob

BASE = r'C:\Users\Fredy\masive-als\compounds'
SRC = os.path.join(BASE, 'next_batches')
OUT = os.path.join(BASE, 'next_batches_pdbqt')
os.makedirs(OUT, exist_ok=True)

TANDAS = [('t3', 'batch3'), ('t6', 'batch6'), ('t7', 'batch7'),
          ('t10', 'batch10'), ('t12', 'batch12')]

def convertir(smiles_csv, outdir):
    """Devuelve (ok, vacios) por filas del CSV."""
    os.makedirs(outdir, exist_ok=True)
    rows = list(csv.DictReader(open(smiles_csv, encoding='utf-8')))
    conv = ob.OBConversion()
    conv.SetInFormat('smi')
    conv.SetOutFormat('pdbqt')
    ok, vacios = 0, 0
    for i, r in enumerate(rows):
        mid = (r.get('molecule_chembl_id') or '').strip()
        smi = (r.get('canonical_smiles') or r.get('smiles') or '').strip()
        if not mid or not smi:
            continue
        nombre = mid.replace('/', '_')
        salida = os.path.join(outdir, nombre + '.pdbqt')
        if os.path.exists(salida) and os.path.getsize(salida) > 100:
            ok += 1
            continue
        mol = ob.OBMol()
        try:
            if not conv.ReadString(mol, smi):
                vacios += 1
                continue
            mol.AddHydrogens()
            b = ob.OBBuilder()
            if not b.Build(mol):
                vacios += 1
                continue
            ff = ob.OBForceField.FindForceField('UFF')
            if ff:
                ff.Setup(mol)
                ff.SteepestDescent(150)
                ff.ConjugateGradients(50)
                ff.GetCoordinates(mol)
            outp = conv.WriteString(mol)
            if outp and len(outp) > 100:
                open(salida, 'w', encoding='utf-8').write(outp)
                ok += 1
            else:
                vacios += 1
        except Exception:
            vacios += 1
        if (i + 1) % 1500 == 0:
            print('  %s: %d/%d (ok=%d, vacios=%d)' % (os.path.basename(outdir), i + 1, len(rows), ok, vacios), flush=True)
    return ok, vacios

def empaquetar(outdir, tanda):
    tgz = os.path.join(BASE, 'ligandos_%s_plano.tar.gz' % tanda)
    with tarfile.open(tgz, 'w:gz') as t:
        for f in sorted(os.listdir(outdir)):
            if f.endswith('.pdbqt'):
                t.add(os.path.join(outdir, f), arcname=f)
    n = len([f for f in os.listdir(outdir) if f.endswith('.pdbqt')])
    print('%s: %d pdbqt -> %s (%d KB)' % (tanda, n, os.path.basename(tgz), os.path.getsize(tgz) // 1024), flush=True)
    return tgz

for t, b in TANDAS:
    csv_file = os.path.join(SRC, '%s_chembl.csv' % t)
    if not os.path.exists(csv_file):
        print('%s: CSV no existe, salto' % t, flush=True)
        continue
    outdir = os.path.join(OUT, b)
    tgz = os.path.join(BASE, 'ligandos_%s_plano.tar.gz' % b)
    if os.path.exists(tgz) and len([f for f in os.listdir(outdir) if f.endswith('.pdbqt')]) > 10:
        print('%s: ya convertido' % b, flush=True)
        continue
    print('=== Convirtiendo %s ===' % b, flush=True)
    t0 = time.time()
    ok, vacios = convertir(csv_file, outdir)
    print('  %s: ok=%d, vacios=%d, %.0fs' % (b, ok, vacios, time.time() - t0), flush=True)
    if ok > 0:
        empaquetar(outdir, b)

print('DONE', flush=True)
