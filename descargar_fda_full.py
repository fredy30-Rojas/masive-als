#!/usr/bin/env python3
"""Descarga todos los farmacos aprobados (max_phase=4) de ChEMBL.
Filtra solo moleculas pequenas con SMILES canonico valido.
Guarda full_fda_library.csv (name, chembl_id, smiles)."""
import urllib.request, json, csv, time, sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE = 'https://www.ebi.ac.uk/chembl/api/data/molecule.json'
OUT = r'C:\Users\Fredy\masive-als\compounds\full_fda_library.csv'
LIMIT = 500

rows = []
offset = 0
total = None
while True:
    url = f'{BASE}?max_phase=4&limit={LIMIT}&offset={offset}'
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (MASIVE-ALS)'})
    try:
        d = json.load(urllib.request.urlopen(req, timeout=60))
    except Exception as e:
        print('ERROR pagina offset', offset, ':', e)
        time.sleep(3)
        continue
    meta = d.get('page_meta', {})
    total = meta.get('total_count', total)
    mols = d.get('molecules', [])
    if not mols:
        break
    for m in mols:
        mtype = (m.get('molecule_type') or '')
        smi = (m.get('molecule_structures') or {}).get('canonical_smiles')
        if mtype != 'Small molecule' or not smi:
            continue
        name = (m.get('pref_name') or m.get('molecule_chembl_id') or '').strip()
        cid = m.get('molecule_chembl_id', '')
        rows.append({'name': name, 'chembl_id': cid, 'smiles': smi})
    offset += len(mols)
    print(f'offset={offset}/{total}  validos={len(rows)}', flush=True)
    if offset >= total:
        break

with open(OUT, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['name', 'chembl_id', 'smiles'])
    w.writeheader()
    w.writerows(rows)

print('TOTAL aprobados:', total)
print('Moleculas pequenas con SMILES:', len(rows))
print('CSV guardado:', OUT)
