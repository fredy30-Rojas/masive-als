#!/usr/bin/env python3
"""Convierte full_fda_library.csv (3311 SMILES aprobados) a PDBQT 3D.
Con resume (salta los ya convertidos) y log de progreso."""
import os, csv, sys, time
from openbabel import openbabel as ob

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

CSV = r'C:\Users\Fredy\masive-als\compounds\full_fda_library.csv'
OUT = r'C:\Users\Fredy\masive-als\compounds\fda_full_pdbqt'
LOG = r'C:\Users\Fredy\masive-als\compounds\fda_convert_log.txt'
os.makedirs(OUT, exist_ok=True)

def log(msg):
    with open(LOG, 'a', encoding='utf-8') as f:
        f.write(time.strftime('%H:%M:%S') + ' ' + msg + '\n')

def generate_3d(smiles, out_path):
    conv = ob.OBConversion()
    conv.SetInFormat('smi')
    mol = ob.OBMol()
    if not conv.ReadString(mol, smiles):
        return False, 'no lee smi'
    mol.AddHydrogens()
    b = ob.OBBuilder()
    if not b.Build(mol):
        return False, 'no gen3d'
    ff = ob.OBForceField.FindForceField('UFF')
    if ff:
        ff.Setup(mol)
        ff.SteepestDescent(150)
        ff.ConjugateGradients(50)
        ff.GetCoordinates(mol)
    conv.SetOutFormat('pdbqt')
    return conv.WriteFile(mol, out_path), ''

rows = []
with open(CSV, encoding='utf-8') as f:
    rows = list(csv.DictReader(f))

def nombre_archivo(name, cid):
    n = (name or '').strip().replace('/', '_').replace(' ', '_').replace('\\', '_')
    if not n:
        n = cid
    return n[:80]

ok = 0
fail = 0
t0 = time.time()
log(f'INICIO: {len(rows)} compuestos')
for i, c in enumerate(rows):
    name = nombre_archivo(c.get('name'), c.get('chembl_id'))
    smiles = (c.get('smiles') or '').strip()
    if not smiles:
        fail += 1
        continue
    out = os.path.join(OUT, name + '.pdbqt')
    if os.path.exists(out) and os.path.getsize(out) > 200:
        ok += 1
        continue
    res, err = generate_3d(smiles, out)
    if res and os.path.exists(out) and os.path.getsize(out) > 200:
        ok += 1
    else:
        fail += 1
        if os.path.exists(out):
            try: os.remove(out)
            except: pass
    if (i + 1) % 200 == 0:
        el = time.time() - t0
        log(f'{i+1}/{len(rows)} ok={ok} fail={fail} {el:.0f}s ({el/(i+1)*1000:.0f}ms/comp)')

el = time.time() - t0
log(f'FIN: ok={ok} fail={fail} tiempo={el:.0f}s')
print(f'Convertidos OK: {ok} | Fallos: {fail} | Tiempo: {el:.0f}s')
print('OUT:', OUT)
