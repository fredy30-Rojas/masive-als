#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fabrica ~100 tandas de ChEMBL. Cada tanda = una query distinta (diana x
tipo de ensayo x umbral, o rodaja de moleculas por fase/año/ATC/flags).
Checkpoint: si el CSV de la tanda ya existe, no la repite.
"""
import urllib.request, json, csv, time, sys, os, hashlib

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE = r'C:\Users\Fredy\masive-als\compounds'
OUTDIR = os.path.join(BASE, 'lotes_100')
os.makedirs(OUTDIR, exist_ok=True)

TARGETS = [('TDP43', 'CHEMBL2362981'), ('SOD1', 'CHEMBL2354'), ('FUS', 'CHEMBL5724679')]

# ===== SMILES ya usados (todas las tandas previas) =====
ya_usados = set()
def cargar_smiles(path, cols=None):
    if not os.path.exists(path):
        return
    if cols is None:
        cols = ['canonical_smiles', 'smiles', 'SMILES', 'smi', 'structure']
    try:
        for r in csv.DictReader(open(path, encoding='utf-8', errors='replace')):
            s = ''
            for c in cols:
                if c in r and r[c]:
                    s = r[c].strip()
                    break
            if not s:
                for v in r.values():
                    v = (v or '').strip()
                    if v and not v.isdigit() and len(v) > 5:
                        s = v
                        break
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

def query_activity(target_id, pmin=None, stype=None, limit=250, max_pages=8):
    res = {}
    offset = 0
    total = None
    fields = 'molecule_chembl_id,parent_molecule_chembl_id,canonical_smiles,pchembl_value'
    while True:
        url = (f'https://www.ebi.ac.uk/chembl/api/data/activity.json'
               f'?target_chembl_id={target_id}&limit={limit}&offset={offset}&fields={fields}')
        if pmin is not None:
            url += f'&pchembl_value__gte={pmin}'
        if stype:
            url += f'&standard_type={stype}'
        try:
            d = get(url)
        except Exception as e:
            print('  ERR offset', offset, str(e)[:60], flush=True)
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
    fields = 'molecule_chembl_id,canonical_smiles,molecule_pref_name'
    while True:
        url = (f'https://www.ebi.ac.uk/chembl/api/data/molecule.json'
               f'?{filtro}&limit={limit}&offset={offset}&fields={fields}')
        try:
            d = get(url)
        except Exception as e:
            print('  ERR offset', offset, str(e)[:60], flush=True)
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
            print('  (cap paginas, %d mols)' % len(res), flush=True)
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
        print('%s: 0 nuevos (desc=%s)' % (tanda, desc), flush=True)
        return
    path = os.path.join(OUTDIR, '%s_chembl.csv' % tanda)
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['molecule_chembl_id', 'canonical_smiles', 'fuente'])
        for mid, smi in nuevos:
            w.writerow([mid, smi, desc])
    print('%s: %d nuevos (%s)' % (tanda, len(nuevos), desc), flush=True)

def ya_lista(tanda):
    return os.path.exists(os.path.join(OUTDIR, '%s_chembl.csv' % tanda))

# ============ DEFINICION DE ~100 TANDAS ============
# Parte 1: dianas x tipo x umbral (48)
TIPOS = [('IC50', 5), ('IC50', 6), ('IC50', 7), ('IC50', 8),
         ('Ki', 5), ('Ki', 6), ('Ki', 7), ('Ki', 8),
         ('EC50', 5), ('EC50', 6), ('EC50', 7),
         ('Kd', 5), ('Kd', 6), ('Kd', 7),
         (None, 7), (None, 8)]
idx = 20
n_parte1 = 0
for tname, tid in TARGETS:
    for stype, pmin in TIPOS:
        idx += 1
        tanda = 'l%03d' % idx
        if ya_lista(tanda):
            continue
        label = '%s %s>=%s' % (tname, stype or 'pchembl', pmin)
        r = query_activity(tid, pmin=pmin, stype=stype)
        guardar(tanda, r, label)
        n_parte1 += 1
        time.sleep(0.3)

# Parte 2: rodajas de moleculas (fases, anos, flags, ATC)
SLICES = [
    ('max_phase=0', 'phase0'),
    ('num_ro5_violations=0', 'ro5_0'),
    ('num_ro5_violations=1', 'ro5_1'),
    ('num_ro5_violations=2', 'ro5_2'),
    ('oral=1', 'oral'),
    ('parenteral=1', 'parenteral'),
    ('topical=1', 'topical'),
    ('orphan=1', 'orphan'),
    ('first_in_class=1', 'first_in_class'),
    ('prodrug=1', 'prodrug'),
    ('withdrawn_flag=1', 'withdrawn'),
    ('first_approval__gte=2020', 'fa_2020'),
    ('first_approval__gte=2015&first_approval__lte=2019', 'fa_1519'),
    ('first_approval__gte=2010&first_approval__lte=2014', 'fa_1014'),
    ('first_approval__gte=2005&first_approval__lte=2009', 'fa_0509'),
    ('first_approval__gte=2000&first_approval__lte=2004', 'fa_0004'),
    ('first_approval__gte=1990&first_approval__lte=1999', 'fa_9099'),
    ('first_approval__gte=1980&first_approval__lte=1989', 'fa_8089'),
    ('natural_product=1&max_phase=3', 'nat_ph3'),
    ('natural_product=1&max_phase=4', 'nat_ph4'),
    ('molecule_type=Small%20molecule&max_phase=3', 'small_ph3'),
    ('availability_type=1', 'avail_1'),
    ('availability_type=4', 'avail_4'),
    ('chirality=1', 'chiral_1'),
    ('chirality=2', 'chiral_2'),
    ('usan_year__gte=2010', 'usan_2010'),
    ('atc_classification__icontains=N', 'atc_N_ner'),
    ('atc_classification__icontains=A', 'atc_A'),
    ('atc_classification__icontains=B', 'atc_B'),
    ('atc_classification__icontains=C', 'atc_C'),
    ('atc_classification__icontains=J', 'atc_J'),
    ('atc_classification__icontains=L', 'atc_L'),
    ('atc_classification__icontains=M', 'atc_M'),
    ('atc_classification__icontains=R', 'atc_R'),
    ('atc_classification__icontains=S', 'atc_S'),
    ('atc_classification__icontains=G', 'atc_G'),
    ('atc_classification__icontains=H', 'atc_H'),
]
n_parte2 = 0
for filtro, nombre in SLICES:
    idx += 1
    tanda = 'l%03d' % idx
    if ya_lista(tanda):
        continue
    r = query_molecules(filtro)
    guardar(tanda, r, nombre)
    n_parte2 += 1
    time.sleep(0.3)

# Parte 3: 3 dianas pchembl>=9 y act_5 para completar 100
for tname, tid in TARGETS:
    for pmin, extra in [(9, 'p9'), (5, None)]:
        if extra is None:
            continue
        idx += 1
        tanda = 'l%03d' % idx
        if ya_lista(tanda):
            continue
        r = query_activity(tid, pmin=pmin)
        guardar(tanda, r, '%s pchembl>=%s' % (tname, pmin))
        n_parte1 += 1
        time.sleep(0.3)

total = sum(1 for f in os.listdir(OUTDIR) if f.endswith('_chembl.csv'))
print('=== TOTAL tandas generadas:', total, '===', flush=True)
print('DONE', flush=True)
