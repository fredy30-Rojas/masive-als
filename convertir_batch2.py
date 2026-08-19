#!/usr/bin/env python3
"""Convierte batch2_chembl.csv (bioactivos TDP43/SOD1/FUS) a PDBQT plano 3D.
Empaqueta ligandos_batch2_plano.tar.gz."""
import os, csv, sys, time, tarfile
from openbabel import openbabel as ob

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

CSV = r'C:\Users\Fredy\masive-als\compounds\batch2_chembl.csv'
OUT = r'C:\Users\Fredy\masive-als\compounds\batch2_pdbqt'
TAR = r'C:\Users\Fredy\masive-als\compounds\ligandos_batch2_plano.tar.gz'
os.makedirs(OUT, exist_ok=True)

def generate_3d(smiles, out_path):
    conv = ob.OBConversion()
    conv.SetInFormat('smi')
    mol = ob.OBMol()
    if not conv.ReadString(mol, smiles):
        return False
    mol.AddHydrogens()
    b = ob.OBBuilder()
    if not b.Build(mol):
        return False
    ff = ob.OBForceField.FindForceField('UFF')
    if ff:
        ff.Setup(mol)
        ff.SteepestDescent(150)
        ff.ConjugateGradients(50)
        ff.GetCoordinates(mol)
    conv.SetOutFormat('pdbqt')
    return conv.WriteFile(mol, out_path)

rows = list(csv.DictReader(open(CSV, encoding='utf-8')))
ok = 0
fail = 0
t0 = time.time()

for i, r in enumerate(rows):
    name = (r.get('name') or r.get('chembl_id') or '').strip()
    smiles = (r.get('smiles') or '').strip()
    if not name or not smiles:
        fail += 1
        continue
    # nombre unico y seguro (usa chembl_id si el nombre esta vacio o raro)
    base = name.replace('/', '_').replace(' ', '_').replace('\\', '_').replace(',', '').replace("'", '')
    if not base or base.upper() == name.upper() and len(base) < 3:
        base = (r.get('chembl_id') or '').strip()
    # asegurar unicidad: si el nombre es un CHEMBL id, usarlo; si no, prefijo
    cid = (r.get('chembl_id') or '').strip()
    fname = base[:70]
    out = os.path.join(OUT, fname + '.pdbqt')
    # si ya existe, saltar
    if os.path.exists(out) and os.path.getsize(out) > 200:
        ok += 1
        continue
    res = generate_3d(smiles, out)
    if res and os.path.exists(out) and os.path.getsize(out) > 200:
        ok += 1
    else:
        fail += 1
        if os.path.exists(out):
            try: os.remove(out)
            except: pass

el = time.time() - t0
print(f'Convertidos OK: {ok} | Fallos: {fail} | Tiempo: {el:.0f}s')

# empaquetar
n = len([f for f in os.listdir(OUT) if f.endswith('.pdbqt')])
print('PDBQT en carpeta:', n)
if os.path.exists(TAR):
    os.remove(TAR)
with tarfile.open(TAR, 'w:gz') as t:
    for f in os.listdir(OUT):
        if f.endswith('.pdbqt'):
            t.add(os.path.join(OUT, f), arcname=f)
print('TAR:', TAR, os.path.getsize(TAR), 'bytes')
