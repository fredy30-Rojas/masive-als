#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Convierte cada tanda de lotes_100 a PDBQT 3D y empaqueta lote_lXXX_plano.tar.gz."""
import csv, os, sys, tarfile, time

sys.stdout.reconfigure(encoding='utf-8', errors='replace')
os.environ.setdefault('BABEL_DATADIR', r'C:\Program Files\OpenBabel-3.1.1\data')

from openbabel import openbabel as ob

BASE = r'C:\Users\Fredy\masive-als\compounds'
SRC = os.path.join(BASE, 'lotes_100')
OUT = os.path.join(BASE, 'lotes_100_pdbqt')
os.makedirs(OUT, exist_ok=True)

def convertir_y_empaquetar(csv_file):
    tanda = os.path.basename(csv_file).replace('_chembl.csv', '')
    outdir = os.path.join(OUT, tanda)
    os.makedirs(outdir, exist_ok=True)
    tgz = os.path.join(BASE, 'lote_%s_plano.tar.gz' % tanda)
    if os.path.exists(tgz) and os.path.getsize(tgz) > 1000:
        print('%s: ya empaquetado' % tanda, flush=True)
        return True
    rows = list(csv.DictReader(open(csv_file, encoding='utf-8')))
    conv = ob.OBConversion()
    conv.SetInFormat('smi')
    conv.SetOutFormat('pdbqt')
    ok, vacios = 0, 0
    for r in rows:
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
    if ok == 0:
        print('%s: 0 pdbqt (vacios=%d), no empaqueto' % (tanda, vacios), flush=True)
        return False
    with tarfile.open(tgz, 'w:gz') as t:
        for f in sorted(os.listdir(outdir)):
            if f.endswith('.pdbqt'):
                t.add(os.path.join(outdir, f), arcname=f)
    print('%s: %d pdbqt -> lote_%s_plano.tar.gz (%d KB)' % (tanda, ok, tanda, os.path.getsize(tgz) // 1024), flush=True)
    return True

if __name__ == '__main__':
    t0 = time.time()
    hecho = 0
    for f in sorted(os.listdir(SRC)):
        if f.endswith('_chembl.csv'):
            if convertir_y_empaquetar(os.path.join(SRC, f)):
                hecho += 1
    print('=== %d tandas empaquetadas en %.0fs ===' % (hecho, time.time() - t0), flush=True)
    print('DONE', flush=True)
