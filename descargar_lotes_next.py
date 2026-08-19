#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Genera las proximas 10 tandas (lotes 3-12) desde ChEMBL con fuentes distintas.
Cada tanda: descarga SMILES -> desduplica contra lotes previos -> guarda CSV.
Conversion a PDBQT y empaquetado se hace con convertir_lotes_next.py despues.
"""
import urllib.request, json, csv, time, sys, os, hashlib

sys.stdout.reconfigure(encoding='utf-8', errors='replace')

BASE = r'C:\Users\Fredy\masive-als\compounds'
OUTDIR = os.path.join(BASE, 'next_batches')
os.makedirs(OUTDIR, exist_ok=True)

TARGETS = [
    ('TDP43', 'CHEMBL2362981'),
    ('SOD1', 'CHEMBL2354'),
    ('FUS', 'CHEMBL5724679'),
]

# ===== Cargar SMILES ya usados (para no duplicar) =====
ya_usados = set()
def cargar_csv_smiles(path, cols=None):
    if not os.path.exists(path):
        return
    if cols is None:
        cols = ['canonical_smiles', 'smiles', 'SMILES', 'smi', 'structure', 'canonical_smilesx']
    try:
        reader = csv.DictReader(open(path, encoding='utf-8', errors='replace'))
        for r in reader:
            s = ''
            for c in cols:
                if c in r and r[c]:
                    s = r[c].strip()
                    break
            if not s and r:
                # ultimo recurso: primer valor que parezca SMILES
                for v in r.values():
                    v = (v or '').strip()
                    if v and not v.isdigit() and len(v) > 5:
                        s = v
                        break
            if s:
                ya_usados.add(hashlib.md5(s.encode()).hexdigest())
    except Exception as e:
        print('AVISO cargando', path, str(e)[:60])

cargar_csv_smiles(os.path.join(BASE, 'batch2_chembl.csv'))
cargar_csv_smiles(os.path.join(BASE, 'full_fda_library.csv'))
cargar_csv_smiles(os.path.join(BASE, 'full_library.csv'))
cargar_csv_smiles(os.path.join(BASE, 'fda_subset.csv'))
for f in os.listdir(OUTDIR):
    if f.endswith('_chembl.csv') and not f.startswith('t'):
        cargar_csv_smiles(os.path.join(OUTDIR, f))
print('SMILES ya usados en lotes previos:', len(ya_usados), flush=True)

UA = {'User-Agent': 'Mozilla/5.0 (MASIVE-ALS)'}

def get(url, retries=8):
    for i in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            return json.load(urllib.request.urlopen(req, timeout=120))
        except Exception as e:
            if i == retries - 1:
                raise
            time.sleep(2.5 * (i + 1))

def query_activity(target_id, pchem_min=None, stype=None, limit=250):
    """Pagina la API de activity de ChEMBL y devuelve set de (molecule_chembl_id, smiles)."""
    res = {}
    offset = 0
    total = None
    fields = 'molecule_chembl_id,parent_molecule_chembl_id,canonical_smiles,molecule_pref_name,standard_type,pchembl_value'
    while True:
        url = (f'https://www.ebi.ac.uk/chembl/api/data/activity.json'
               f'?target_chembl_id={target_id}&limit={limit}&offset={offset}&fields={fields}')
        if pchem_min is not None:
            url += f'&pchembl_value__gte={pchem_min}'
        if stype:
            url += f'&standard_type={stype}'
        try:
            d = get(url)
        except Exception as e:
            print('  ERROR offset', offset, str(e)[:70], flush=True)
            break
        meta = d.get('page_meta', {})
        total = meta.get('total_count', total)
        acts = d.get('activities', [])
        if not acts:
            break
        for a in acts:
            pid = a.get('parent_molecule_chembl_id') or a.get('molecule_chembl_id')
            smi = (a.get('canonical_smiles') or '').strip()
            pc = a.get('pchembl_value')
            if not pid or not smi:
                continue
            try:
                pcv = float(pc) if pc else 0.0
            except Exception:
                pcv = 0.0
            if pchem_min is not None and pcv < pchem_min:
                continue
            if pid not in res or pcv > res[pid][0]:
                res[pid] = (pcv, smi)
        offset += limit
        if total is not None and offset >= total:
            break
        if offset > 20000:
            break
        time.sleep(0.3)
    return res

def query_molecules(filtro_extra, limit=250, max_pages=60):
    """Pagina la API de molecule de ChEMBL con un filtro (ej: max_phase=3)."""
    res = {}
    offset = 0
    total = None
    fields = 'molecule_chembl_id,canonical_smiles,molecule_pref_name'
    while True:
        url = (f'https://www.ebi.ac.uk/chembl/api/data/molecule.json'
               f'?{filtro_extra}&limit={limit}&offset={offset}&fields={fields}')
        try:
            d = get(url)
        except Exception as e:
            print('  ERROR offset', offset, str(e)[:70], flush=True)
            break
        meta = d.get('page_meta', {})
        total = meta.get('total_count', total)
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
            print('  (cota de paginas alcanzada: %d moleculas)' % len(res), flush=True)
            break
        time.sleep(0.3)
    return res

def guardar(tanda, res, filtro):
    """Filtra no-usados, escribe CSV y devuelve el nombre."""
    nuevos = []
    for mid, val in res.items():
        smi = val if isinstance(val, str) else val[1]
        if not smi:
            continue
        h = hashlib.md5(smi.encode()).hexdigest()
        if h in ya_usados:
            continue
        ya_usados.add(h)
        nuevos.append((mid, smi))
    if not nuevos:
        print(f'{tanda}: 0 compuestos nuevos', flush=True)
        return None
    path = os.path.join(OUTDIR, f'{tanda}_chembl.csv')
    with open(path, 'w', newline='', encoding='utf-8') as f:
        w = csv.writer(f)
        w.writerow(['molecule_chembl_id', 'canonical_smiles', 'fuente'])
        for mid, smi in nuevos:
            w.writerow([mid, smi, filtro])
    print(f'{tanda}: {len(nuevos)} compuestos nuevos -> {os.path.basename(path)}', flush=True)
    return path

# ================= TANDAS =================
def ya_lista(tanda):
    p = os.path.join(OUTDIR, f'{tanda}_chembl.csv')
    return os.path.exists(p)

print('=== T3: TDP-43 pChEMBL>=6 ===', flush=True)
if not ya_lista('t3'): guardar('t3', query_activity('CHEMBL2362981', pchem_min=6), 'TDP43 pChEMBL>=6')
else: print('t3 ya existe', flush=True)

print('=== T4: SOD1 pChEMBL>=6 ===', flush=True)
if not ya_lista('t4'): guardar('t4', query_activity('CHEMBL2354', pchem_min=6), 'SOD1 pChEMBL>=6')
else: print('t4 ya existe', flush=True)

print('=== T5: FUS pChEMBL>=5 ===', flush=True)
if not ya_lista('t5'): guardar('t5', query_activity('CHEMBL5724679', pchem_min=5), 'FUS pChEMBL>=5')
else: print('t5 ya existe', flush=True)

print('=== T6: max_phase=3 (fase III) ===', flush=True)
if not ya_lista('t6'): guardar('t6', query_molecules('max_phase=3'), 'max_phase=3')
else: print('t6 ya existe', flush=True)

print('=== T7: max_phase=2 (fase II) ===', flush=True)
if not ya_lista('t7'): guardar('t7', query_molecules('max_phase=2'), 'max_phase=2')
else: print('t7 ya existe', flush=True)

print('=== T8: 3 dianas pChEMBL>=7 (potentes) ===', flush=True)
if not ya_lista('t8'):
    t8 = {}
    for tname, tid in TARGETS:
        r = query_activity(tid, pchem_min=7)
        print(f'  {tname}: {len(r)} act.', flush=True)
        for k, v in r.items():
            if k not in t8 or v[0] > t8[k][0]:
                t8[k] = v
    guardar('t8', t8, '3dianas pChEMBL>=7')
else: print('t8 ya existe', flush=True)

print('=== T9: 3 dianas standard_type=Ki pChEMBL>=5 ===', flush=True)
if not ya_lista('t9'):
    t9 = {}
    for tname, tid in TARGETS:
        r = query_activity(tid, pchem_min=5, stype='Ki')
        print(f'  {tname}: {len(r)} act.', flush=True)
        for k, v in r.items():
            if k not in t9 or v[0] > t9[k][0]:
                t9[k] = v
    guardar('t9', t9, '3dianas Ki>=5')
else: print('t9 ya existe', flush=True)

print('=== T10: max_phase=1 (fase I) ===', flush=True)
if not ya_lista('t10'): guardar('t10', query_molecules('max_phase=1'), 'max_phase=1')
else: print('t10 ya existe', flush=True)

print('=== T11: 3 dianas standard_type=IC50 pChEMBL>=5 ===', flush=True)
if not ya_lista('t11'):
    t11 = {}
    for tname, tid in TARGETS:
        r = query_activity(tid, pchem_min=5, stype='IC50')
        print(f'  {tname}: {len(r)} act.', flush=True)
        for k, v in r.items():
            if k not in t11 or v[0] > t11[k][0]:
                t11[k] = v
    guardar('t11', t11, '3dianas IC50>=5')
else: print('t11 ya existe', flush=True)

print('=== T12: productos naturales (natural_product=1) ===', flush=True)
if not ya_lista('t12'): guardar('t12', query_molecules('natural_product=1', max_pages=40), 'natural_product=1')
else: print('t12 ya existe', flush=True)

print('=== RESUMEN ===', flush=True)
for f in sorted(os.listdir(OUTDIR)):
    if f.endswith('_chembl.csv'):
        print(' ', f, '-', sum(1 for _ in open(os.path.join(OUTDIR, f), encoding='utf-8')) - 1, 'compuestos', flush=True)
print('DONE', flush=True)
