#!/usr/bin/env python3
"""Convierte full_library.csv (SMILES) a PDBQT 3D para docking con OpenBabel."""
import os, csv, sys, time
from openbabel import openbabel as ob

BASE = r'C:\Users\Fredy\masive-als'
CSV = os.path.join(BASE, 'compounds', 'full_library.csv')
OUT = r'C:\Users\Fredy\Desktop\MASIVE-ALS-Colab\fda_pdbqt'
os.makedirs(OUT, exist_ok=True)

def generate_3d_structure(smiles, output_pdbqt):
    conv = ob.OBConversion()
    conv.SetInFormat("smi")
    mol = ob.OBMol()
    if not conv.ReadString(mol, smiles):
        return False, "no lee smi"
    mol.AddHydrogens()
    builder = ob.OBBuilder()
    builder.Build(mol)
    ff = ob.OBForceField.FindForceField("UFF")
    if ff:
        ff.Setup(mol)
        ff.SteepestDescent(200)
        ff.ConjugateGradients(100)
        ff.GetCoordinates(mol)
    conv.SetOutFormat("pdbqt")
    return conv.WriteFile(mol, output_pdbqt), ""

rows = []
with open(CSV, encoding='utf-8') as f:
    for r in csv.DictReader(f):
        rows.append(r)

ok = 0
fail = []
t0 = time.time()
for i, comp in enumerate(rows):
    name = (comp.get('name') or '').strip().replace('/', '_').replace(' ', '_')
    smiles = (comp.get('smiles') or '').strip()
    if not name or not smiles or smiles in ('N/A', ''):
        continue
    out = os.path.join(OUT, name + '.pdbqt')
    if os.path.exists(out) and os.path.getsize(out) > 200:
        ok += 1
        continue
    res, err = generate_3d_structure(smiles, out)
    if res and os.path.exists(out) and os.path.getsize(out) > 200:
        ok += 1
    else:
        fail.append((name, err or 'fallo'))

print('Compuestos CSV:', len(rows))
print('Convertidos OK:', ok)
print('Fallos:', len(fail))
for n, e in fail:
    print('  FALLO:', n, '-', e)
print('Tiempo: %.1f s' % (time.time() - t0))
print('OUT:', OUT)
