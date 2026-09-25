#!/usr/bin/env python3
"""Fix SMILES in candidatos_limpios.csv - validate with RDKit and canonicalize."""
import csv
import os
from rdkit import Chem

CSV = '/home/ubuntu/mmgbsa/candidatos_limpios.csv'
POSES_DIR = '/home/ubuntu/mmgbsa/poses'

rows = []
skipped = 0
with open(CSV) as f:
    for row in csv.DictReader(f):
        ligand = row['ligand']
        pose_path = os.path.join(POSES_DIR, ligand + '_out.pdbqt')
        
        # Extract SMILES from PDBQT REMARK
        smiles = ''
        if os.path.exists(pose_path):
            with open(pose_path, encoding='utf-8', errors='ignore') as pf:
                for line in pf:
                    if 'REMARK SMILES ' in line and 'IDX' not in line:
                        m = __import__('re').search(r'REMARK SMILES (.+)', line)
                        if m:
                            smiles = m.group(1).strip()
                            break
        
        if not smiles:
            skipped += 1
            continue
            
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            skipped += 1
            continue
        
        row['smiles'] = Chem.MolToSmiles(mol)
        rows.append(row)

with open(CSV, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['ligand', 'target', 'affinity', 'smiles', 'pose_pdbqt'])
    w.writeheader()
    w.writerows(rows)

print(f'Validos: {len(rows)}, Skipped: {skipped}')
