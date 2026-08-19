#!/usr/bin/env python3
"""Descarga compuestos bioactivos contra TDP-43, SOD1 y FUS desde ChEMBL.
Usa fields= para reducir la respuesta y evitar el 500 del servidor.
Desduplica por molecula y guarda batch2_chembl.csv."""
import urllib.request, json, csv, time, sys

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

TARGETS = [
    ('TDP43', 'CHEMBL2362981'),
    ('SOD1', 'CHEMBL2354'),
    ('FUS', 'CHEMBL5724679'),
]
OUT = r'C:\Users\Fredy\masive-als\compounds\batch2_chembl.csv'
LIMIT = 250
FIELDS = 'molecule_chembl_id,parent_molecule_chembl_id,canonical_smiles,molecule_pref_name,standard_type,standard_flag'

def get(url, retries=6):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (MASIVE-ALS)'})
            return json.load(urllib.request.urlopen(req, timeout=90))
        except Exception as e:
            if i == retries - 1:
                raise
            time.sleep(3 * (i + 1))

uniques = {}

for tname, tid in TARGETS:
    offset = 0
    total = None
    while True:
        url = (f'https://www.ebi.ac.uk/chembl/api/data/activity.json'
               f'?target_chembl_id={tid}&limit={LIMIT}&offset={offset}&fields={FIELDS}')
        try:
            d = get(url)
        except Exception as e:
            print(f'{tname} ERROR offset {offset}: {e}', flush=True)
            break
        meta = d.get('page_meta', {})
        total = meta.get('total_count', total)
        acts = d.get('activities', [])
        if not acts:
            break
        for a in acts:
            parent = a.get('parent_molecule_chembl_id') or a.get('molecule_chembl_id')
            smi = (a.get('canonical_smiles') or '').strip()
            if not parent or not smi:
                continue
            if parent not in uniques:
                name = (a.get('molecule_pref_name') or parent).strip()
                uniques[parent] = (name, parent, smi)
        offset += len(acts)
        if total and offset >= total:
            break
        if offset % 2000 < LIMIT:
            print(f'{tname}: {offset}/{total} unicos={len(uniques)}', flush=True)

print('TOTAL compuestos unicos:', len(uniques), flush=True)

rows = [{'name': v[0], 'chembl_id': v[1], 'smiles': v[2]} for v in uniques.values()]
rows.sort(key=lambda r: r['name'].lower())
with open(OUT, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['name', 'chembl_id', 'smiles'])
    w.writeheader()
    w.writerows(rows)
print('CSV guardado:', OUT, 'con', len(rows), 'compuestos')
