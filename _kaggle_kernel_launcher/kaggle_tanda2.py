# # MASIVE-ALS — Cribado Tanda 2: 1000 ligandos nuevos contra TDP43_v2
#
# Descarga ligandos_tanda2.tar.gz desde GitHub (1000 ligandos de dsA)
# Usa Vina-GPU-2.1 con checkpoint. 3 semillas (42, 2026, 777).

try:
    subprocess.run(['nvidia-smi'], check=False)
except Exception:
    pass
print('GPU OK')

# CELDA 2: Obtener binario Vina-GPU
import os, glob, subprocess, urllib.request, tarfile, shutil

r = subprocess.run(['apt-get', 'update', '-qq'], capture_output=True, text=True)
r = subprocess.run(['apt-get', 'install', '-y', '-qq', 'libboost-program-options1.74.0'], capture_output=True, text=True)
print('boost instalado rc=%d' % r.returncode)
if r.returncode != 0:
    subprocess.run(["pip", "install", "-q", "libboost"], capture_output=True, text=True)

VINA_GPU_URL = 'https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/vinagpu_linux.tar.gz'
BIN_DIR = '/kaggle/working/vinagpu_linux'
os.makedirs(BIN_DIR, exist_ok=True)

def buscar_binario(base):
    if not base or not os.path.exists(base):
        return None
    hits = glob.glob(base + '/**/AutoDock-Vina-GPU-2-1', recursive=True)
    return hits[0] if hits else None

VINA_GPU_BIN = buscar_binario(BIN_DIR)
if not VINA_GPU_BIN:
    pkg = '/tmp/vinagpu_linux.tar.gz'
    urllib.request.urlretrieve(VINA_GPU_URL, pkg)
    with tarfile.open(pkg) as t:
        t.extractall(BIN_DIR)
    VINA_GPU_BIN = buscar_binario(BIN_DIR)
if not VINA_GPU_BIN:
    raise SystemExit('ERROR: no se pudo obtener Vina-GPU')
BIN_DIR = os.path.dirname(VINA_GPU_BIN)
os.chmod(VINA_GPU_BIN, 0o755)
print('Vina-GPU listo:', VINA_GPU_BIN)

# CELDA 3: Descargar ligandos tanda 2 + receptor desde GitHub
WORK = '/kaggle/working/masive_als'
for sub in ['receptores', 'ligandos', 'resultados', 'checkpoint']:
    os.makedirs(WORK + '/' + sub, exist_ok=True)

TANDA2_URL = 'https://raw.githubusercontent.com/fredy30-Rojas/masive-als/master/_kaggle_kernel/ligandos_tanda2_clean.tar.gz'
pkg2 = '/tmp/ligandos_tanda2.tar.gz'
print('Descargando ligandos_tanda2_clean.tar.gz ...')
urllib.request.urlretrieve(TANDA2_URL, pkg2)
with tarfile.open(pkg2) as t:
    t.extractall(WORK)

print('Receptores:', len(glob.glob(WORK + '/receptores/*.pdbqt')))
print('Ligandos:', len(glob.glob(WORK + '/ligandos/*.pdbqt')))

# CELDA 4: Receptor TDP43_v2
RECEPTORES = {
    'TDP43_v2': {'archivo': WORK + '/receptores/TDP43_v2.pdbqt', 'centro': [24.23, 16.89, -15.87], 'tamano': [26, 26, 26]},
}
print('Receptor TDP43_v2:', os.path.exists(RECEPTORES['TDP43_v2']['archivo']))

# CELDA 5: Pipeline de docking
import csv, glob, os, time, subprocess

LIG_DIR = WORK + '/ligandos'
OUT_DIR = WORK + '/resultados'
CSV = OUT_DIR + '/resultados_tanda2.csv'
CKPT = OUT_DIR + '/hechos_tanda2.txt'

if not os.path.exists(CSV):
    with open(CSV, 'w', newline='') as f:
        csv.writer(f).writerow(['ligand', 'target', 'energy', 'timestamp'])

hechos = set()
if os.path.exists(CKPT):
    hechos = set(l.strip() for l in open(CKPT) if l.strip())
print('Hechos antes:', len(hechos))

ligandos = sorted(glob.glob(LIG_DIR + '/*.pdbqt'))
print('Ligandos totales:', len(ligandos))

SEEDS = [42, 2026, 777]

def acoplar(lig, target, seed, thread=8000):
    info = RECEPTORES[target]
    cfg = '/tmp/cfg_%s.txt' % os.path.basename(lig).replace('.pdbqt', '')[:30]
    with open(cfg, 'w') as f:
        f.write('receptor = %s\n' % info['archivo'])
        f.write('ligand = %s\n' % lig)
        f.write('center_x = %s\n' % info['centro'][0])
        f.write('center_y = %s\n' % info['centro'][1])
        f.write('center_z = %s\n' % info['centro'][2])
        f.write('size_x = %s\n' % info['tamano'][0])
        f.write('size_y = %s\n' % info['tamano'][1])
        f.write('size_z = %s\n' % info['tamano'][2])
        f.write('num_modes = 3\n')
        f.write('seed = %d\n' % seed)
        f.write('thread = %d\n' % thread)
    try:
        r = subprocess.run([VINA_GPU_BIN, '--config', cfg], capture_output=True, text=True,
                           timeout=1800, cwd=BIN_DIR)
        out = (r.stdout or '') + (r.stderr or '')
        if r.returncode != 0:
            return None, 'rc=%d %s' % (r.returncode, out[-150:])
        for ln in out.splitlines():
            s = ln.split()
            if len(s) >= 2 and s[0] == '1':
                try:
                    return round(float(s[1]), 4), None
                except ValueError:
                    pass
        return None, 'sin afinidad: ' + out[-150:]
    except Exception as ex:
        return None, str(ex)[:100]

def tiene_atomos(lig):
    try:
        with open(lig, errors='replace') as f:
            return any(l.startswith(('ATOM', 'HETATM')) for l in f)
    except Exception:
        return False

t0 = time.time()
n_ok = 0
n_err = 0
for lig in ligandos:
    nombre = os.path.basename(lig).replace('.pdbqt', '')
    if nombre in hechos:
        continue
    if not tiene_atomos(lig):
        print('LIGANDO_VACIO', nombre)
        hechos.add(nombre)
        continue
    for target in RECEPTORES:
        energias = []
        err = None
        for sd in SEEDS:
            e, er = acoplar(lig, target, sd)
            if e is not None:
                energias.append(e)
            else:
                err = er
        if energias:
            energia = min(energias)
            with open(CSV, 'a', newline='') as f:
                csv.writer(f).writerow([nombre, target, energia, time.strftime('%Y-%m-%d %H:%M:%S')])
            n_ok += 1
        else:
            n_err += 1
            print('ERROR', nombre, target, err)
    with open(CKPT, 'a') as f:
        f.write(nombre + '\n')
    hechos.add(nombre)
    if n_ok % 10 == 0 and n_ok > 0:
        print('[%d acoplados, %d errores] %.1f min' % (n_ok, n_err, (time.time() - t0) / 60), flush=True)

print()
print('=== TANDA 2 COMPLETADA ===')
print('Acoplados:', n_ok, '| Errores:', n_err)
print('CSV:', CSV)

# CELDA 6: Resumen
import csv
rows = list(csv.DictReader(open(CSV)))
print('Total resultados tanda 2:', len(rows))
if rows:
    best = sorted(rows, key=lambda x: float(x['energy']))[:10]
    print()
    print('Top 10 tanda 2:')
    for r in best:
        print('  ', r['ligand'], r['target'], r['energy'])
print()
print('Descargar: resultados_tanda2.csv')
