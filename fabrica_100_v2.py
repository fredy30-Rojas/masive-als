#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fabrica ~100 tandas de ChEMBL (v2 con filtros verificados).
Parte 1: dianas ALS/neuro encontradas por busqueda + actives pchembl>=6.
Parte 2: rodajas de moleculas (fases, anos, flags) verificadas.
Checkpoint por CSV. Desduplica contra todo lo previo.
"""
import urllib.request, json, csv, time, sys, os, hashlib

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE = r'C:\Users\Fredy\masive-als\compounds'
OUTDIR = os.path.join(BASE, 'lotes_100')
os.makedirs(OUTDIR, exist_ok=True)

# ===== SMILES ya usados =====
ya_usados = set()
def cargar_smiles(path):
    if not os.path.exists(path):
        return
    try:
        for r in csv.DictReader(open(path, encoding='utf-8', errors='replace')):
            s = (r.get('canonical_smiles') or r.get('smiles') or '').strip()
            if s:
                ya_usados.add(hashlib.md5(s.encode()).hexdigest())
    except Exception:
        pass

cargar_smiles(os.path.join(BASE, 'batch2_chembl.csv'))
cargar_smiles(os.path.join(BASE, 'full_fda_library.csv'))
cargar_smiles(os.path.join(BASE, 'full_library.csv'))
cargar_smiles(os.path.join(BASE, 'fda_subset.csv'))
for d in ['next_batches', 'lotes_100']:
    dd = os.path.join(BASE, d)
    if os.path.isdir(dd):
        for f in os.listdir(dd):
            if f.endswith('_chembl.csv'):
                cargar_smiles(os.path.join(dd, f))
print('SMILES ya usados:', len(ya_usados), flush=True)

UA = {'User-Agent': 'Mozilla/5.0 (MASIVE-ALS)'}

def get(url, retries=6):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            return json.load(urllib.request.urlopen(req, timeout=120))
        except Exception as e:
            if i == retries - 1:
                raise
            time.sleep(2 * (i + 1))

def buscar_targets(palabra, limit=10):
    url = (f'https://www.ebi.ac.uk/chembl/api/data/target.json'
           f'?pref_name__icontains={urllib.parse.quote(palabra)}&limit={limit}&fields=target_chembl_id,pref_name,organism')
    try:
        d = get(url)
        return [(t.get('target_chembl_id'), t.get('pref_name'), t.get('organism')) for t in d.get('targets', [])]
    except Exception:
        return []

def query_activity(tid, pmin=6, limit=250, max_pages=6):
    res = {}
    offset = 0
    total = None
    fields = 'molecule_chembl_id,parent_molecule_chembl_id,canonical_smiles,pchembl_value'
    while True:
        url = (f'https://www.ebi.ac.uk/chembl/api/data/activity.json'
               f'?target_chembl_id={tid}&limit={limit}&offset={offset}&fields={fields}'
               f'&pchembl_value__gte={pmin}')
        try:
            d = get(url)
        except Exception as e:
            print('  ERR', tid, str(e)[:50], flush=True)
            break
        total = d.get('page_meta', {}).get('total_count', total)
        acts = d.get('activities', [])
        if not acts:
            break
        for a in acts:
            pid = a.get('parent_molecule_chembl_id') or a.get('molecule_chembl_id')
            smi = (a.get('canonical_smiles') or '').strip()
            if pid and smi:
                res[pid] = smi
        offset += limit
        if total is not None and offset >= total:
            break
        if offset >= limit * max_pages:
            break
        time.sleep(0.25)
    return res

def query_molecules(filtro, limit=250, max_pages=8):
    res = {}
    offset = 0
    total = None
    fields = 'molecule_chembl_id,canonical_smiles'
    while True:
        url = (f'https://www.ebi.ac.uk/chembl/api/data/molecule.json'
               f'?{filtro}&limit={limit}&offset={offset}&fields={fields}')
        try:
            d = get(url)
        except Exception as e:
            print('  ERR', str(e)[:50], flush=True)
            break
        total = d.get('page_meta', {}).get('total_count', total)
        mols = d.get('molecules', [])
        if not mols:
            break
        for m in mols:
            mid = m.get('molecule_chembl_id')
            smi = ((m.get('molecule_structures') or {}).get('canonical_smiles') or '').strip()
            if mid and smi:
                res[mid] = smi
        offset += limit
        if total is not None and offset >= total:
            break
        if offset >= limit * max_pages:
            break
        time.sleep(0.25)
    return res

def guardar(tanda, res, desc):
    nuevos = []
    for mid, smi in res.items():
        if not smi:
            continue
        h = hashlib.md5(smi.encode()).hexdigest()
        if h in ya_usados:
            continue
        ya_usados.add(h)
        nuevos.append((mid, smi))
    if not nuevos:
        print('%s: 0 nuevos (%s)' % (tanda, desc), flush=True)
        return False
    path = os.path.join(OUTDIR, '%s_chembl.csv' % tanda)
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['molecule_chembl_id', 'canonical_smiles', 'fuente'])
        for mid, smi in nuevos:
            w.writerow([mid, smi, desc])
    print('%s: %d nuevos (%s)' % (tanda, len(nuevos), desc), flush=True)
    return True

def ya_lista(tanda):
    return os.path.exists(os.path.join(OUTDIR, '%s_chembl.csv' % tanda))

import urllib.parse
idx = 99  # l100 en adelante

# ===== Parte 1: dianas ALS/neuro =====
PALABRAS = ['sigma-1', 'sigma receptor', 'glutamate', 'acetylcholinesterase', 'monoamine oxidase',
            'glycogen synthase kinase', 'cyclin dependent kinase 5', 'rho kinase', 'histone deacetylase 6',
            'sirtuin 1', 'poly polymerase', 'mTOR', 'AMPK', 'proteasome', 'ubiquitin', 'autophagy',
            'NRF2', 'KEAP1', 'TRPV', 'KCNQ', 'TRPM', 'GABA', 'serotonin', 'dopamine', 'opioid',
            'cannabinoid', 'tau protein', 'alpha-synuclein', 'superoxide dismutase', 'TARDBP',
            'FUS RNA', 'hnRNP', 'TIA1', 'ataxin', 'valosin', 'p62', 'optineurin', 'TBK1',
            'MATR3', 'sphingosine kinase', 'lysophosphatidic acid', 'sigma receptor', 'matrix metalloproteinase',
            'cathepsin', 'caspase', 'parkin', 'PINK1', 'LRRK2', 'dynein', 'kinesin', 'tubulin']

targets_usados = set()
n_dianas = 0
for pal in PALABRAS:
    for tid, pref, org in buscar_targets(pal, limit=8):
        if not tid or tid in targets_usados:
            continue
        # solo humanos
        if org and 'Homo sapiens' not in org:
            continue
        targets_usados.add(tid)
        idx += 1
        tanda = 'l%03d' % idx
        if ya_lista(tanda):
            continue
        r = query_activity(tid, pmin=6)
        if len(r) < 5:
            r2 = query_activity(tid, pmin=5)
            if len(r2) > len(r):
                r = r2
        ok = guardar(tanda, r, 'diana:%s (%s)' % (pref[:40], tid))
        if ok:
            n_dianas += 1
        time.sleep(0.3)

print('=== dianas con batch:', n_dianas, '===', flush=True)

# ===== Parte 2: rodajas de moleculas (verificadas) =====
SLICES = [
    ('max_phase=0', 'phase0'),
    ('oral=1', 'oral'),
    ('parenteral=1', 'parenteral'),
    ('orphan=1', 'orphan'),
    ('first_in_class=1', 'first_in_class'),
    ('prodrug=1', 'prodrug'),
    ('withdrawn_flag=1', 'withdrawn'),
    ('first_approval__gte=2020', 'fa_2020'),
    ('first_approval__gte=2015&first_approval__lte=2019', 'fa_1519'),
    ('first_approval__gte=2010&first_approval__lte=2014', 'fa_1014'),
    ('first_approval__gte=2000&first_approval__lte=2009', 'fa_0009'),
    ('first_approval__gte=1990&first_approval__lte=1999', 'fa_9099'),
    ('availability_type=1', 'avail_1'),
    ('chirality=1', 'chiral_1'),
    ('chirality=2', 'chiral_2'),
    ('usan_year__gte=2010', 'usan_2010'),
    ('molecule_type=Small%20molecule&max_phase=3', 'small_ph3'),
    ('natural_product=1&max_phase=3', 'nat_ph3'),
    ('natural_product=1&max_phase=4', 'nat_ph4'),
    ('oral=1&natural_product=1', 'oral_nat'),
]
for filtro, nombre in SLICES:
    idx += 1
    tanda = 'l%03d' % idx
    if ya_lista(tanda):
        continue
    r = query_molecules(filtro)
    guardar(tanda, r, nombre)
    time.sleep(0.3)

total = sum(1 for f in os.listdir(OUTDIR) if f.endswith('_chembl.csv'))
print('=== TOTAL tandas generadas:', total, '===', flush=True)
print('DONE', flush=True)
