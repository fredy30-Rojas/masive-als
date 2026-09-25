import csv
import glob
import os
import re
import stat
import subprocess
import tarfile
import time
import urllib.request

DATASET_ROOT = "/kaggle/input/masive-als-sod1-smoke-data"
WORK = "/kaggle/working/masive_als_sod1_smoke"
RECEPTOR = os.path.join(DATASET_ROOT, "SOD1.pdbqt")
EXPECTED = ["CHEMBL3311449", "CHEMBL9532", "CHEMBL8905"]
SEEDS = [42, 2026, 777]
CENTER = (46.5, 80.0, 73.3)
BOX = (22, 22, 22)


def fail(message):
    raise RuntimeError("VALIDATION FAILED: " + message)


def find_dataset_root():
    candidates = [DATASET_ROOT, "/kaggle/src", "/kaggle/working"]
    candidates.extend(glob.glob("/kaggle/input/*"))
    for root in candidates:
        if os.path.isfile(os.path.join(root, "SOD1.pdbqt")):
            return root
    # Kaggle scripts may expose attached data under a nested working directory.
    for root, _, files in os.walk("/kaggle"):
        if "SOD1.pdbqt" in files:
            return root
    fail("SOD1.pdbqt no encontrado en ninguna ruta de Kaggle")


def parse_affinity(text):
    patterns = [
        r"^\s*1\s+(-?\d+(?:\.\d+)?)\s+",
        r"REMARK VINA RESULT:\s+(-?\d+(?:\.\d+)?)",
    ]
    for line in text.splitlines():
        for pattern in patterns:
            match = re.search(pattern, line)
            if match:
                return float(match.group(1))
    return None


def ensure_runtime_dependencies():
    result = subprocess.run(
        ["apt-get", "update", "-qq"],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    if result.returncode != 0:
        fail("no se pudo actualizar el índice de paquetes: " + (result.stderr or "")[-300:])
    result = subprocess.run(
        ["apt-get", "install", "-y", "-qq", "libboost-program-options1.74.0"],
        capture_output=True,
        text=True,
        timeout=300,
        check=False,
    )
    if result.returncode != 0:
        fail("no se pudo instalar libboost-program-options1.74.0: " + (result.stderr or "")[-300:])


def locate_vina():
    existing = glob.glob("/kaggle/working/**/AutoDock-Vina-GPU-2-1", recursive=True)
    if existing:
        return existing[0]
    url = "https://raw.githubusercontent.com/fredy30-Rojas/masive-als-data/main/vinagpu_linux.tar.gz"
    archive = "/tmp/vinagpu_linux.tar.gz"
    urllib.request.urlretrieve(url, archive)
    target = "/kaggle/working/vinagpu_linux"
    os.makedirs(target, exist_ok=True)
    with tarfile.open(archive, "r:gz") as archive_file:
        archive_file.extractall(target)
    existing = glob.glob(target + "/**/AutoDock-Vina-GPU-2-1", recursive=True)
    if not existing:
        fail("no se pudo localizar AutoDock-Vina-GPU-2-1")
    binary = existing[0]
    os.chmod(binary, os.stat(binary).st_mode | stat.S_IEXEC)
    return binary


def run_one(binary, receptor, ligand, seed, output_dir):
    name = os.path.splitext(os.path.basename(ligand))[0]
    config = os.path.join(output_dir, "%s_seed_%s.conf" % (name, seed))
    with open(config, "w", encoding="utf-8") as handle:
        handle.write("receptor = %s\n" % receptor)
        handle.write("ligand = %s\n" % ligand)
        handle.write("center_x = %s\n" % CENTER[0])
        handle.write("center_y = %s\n" % CENTER[1])
        handle.write("center_z = %s\n" % CENTER[2])
        handle.write("size_x = %s\n" % BOX[0])
        handle.write("size_y = %s\n" % BOX[1])
        handle.write("size_z = %s\n" % BOX[2])
        handle.write("num_modes = 3\n")
        handle.write("seed = %s\n" % seed)
        handle.write("thread = 8000\n")
    result = subprocess.run(
        [binary, "--config", config],
        cwd=os.path.dirname(binary),
        capture_output=True,
        text=True,
        timeout=1800,
        check=False,
    )
    text = (result.stdout or "") + "\n" + (result.stderr or "")
    with open(os.path.join(output_dir, "%s_seed_%s.log" % (name, seed)), "w", encoding="utf-8") as handle:
        handle.write(text)
    if result.returncode != 0:
        return None, "returncode=%s" % result.returncode
    affinity = parse_affinity(text)
    if affinity is None:
        return None, "afinidad no encontrada"
    return affinity, ""


root = find_dataset_root()
source_receptor = os.path.join(root, "SOD1.pdbqt")
source_ligands = []
for compound in EXPECTED:
    path = os.path.join(root, compound + ".pdbqt")
    if not os.path.isfile(path):
        fail("falta el ligando " + compound)
    if os.path.getsize(path) == 0:
        fail("ligando vacío " + compound)
    source_ligands.append(path)

all_pdbqt = [p for p in glob.glob(os.path.join(root, "*.pdbqt")) if os.path.basename(p) != "SOD1.pdbqt"]
if sorted(os.path.splitext(os.path.basename(p))[0] for p in all_pdbqt) != sorted(EXPECTED):
    fail("el dataset no contiene exactamente los tres ligandos esperados")

# Vina writes *_out.pdbqt next to the ligand, so work from writable storage.
os.makedirs(WORK, exist_ok=True)
input_root = os.path.join(WORK, "inputs")
os.makedirs(input_root, exist_ok=True)
receptor = os.path.join(input_root, "SOD1.pdbqt")
os.replace(source_receptor, receptor) if source_receptor == receptor else __import__("shutil").copy2(source_receptor, receptor)
ligands = []
for source in source_ligands:
    destination = os.path.join(input_root, os.path.basename(source))
    __import__("shutil").copy2(source, destination)
    ligands.append(destination)

ensure_runtime_dependencies()
binary = locate_vina()
rows = []
for ligand in ligands:
    compound = os.path.splitext(os.path.basename(ligand))[0]
    energies = []
    errors = []
    for seed in SEEDS:
        energy, error = run_one(binary, receptor, ligand, seed, WORK)
        if energy is not None:
            energies.append(energy)
        else:
            errors.append(error)
    if not energies:
        fail("sin afinidad válida para " + compound + ": " + "; ".join(errors))
    rows.append({
        "ligand": compound,
        "target": "SOD1",
        "n_seeds": len(energies),
        "vina_affinity_min": min(energies),
        "vina_affinity_median": sorted(energies)[len(energies) // 2],
        "seed_affinities": ";".join("%.4f" % value for value in energies),
        "status": "ok",
    })

if len(rows) != 3 or any(row["status"] != "ok" for row in rows):
    fail("el CSV final no contiene tres resultados válidos")

csv_path = os.path.join(WORK, "resultado_smoke_sod1.csv")
with open(csv_path, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

print("=== PRUEBA SOD1 COMPLETADA ===")
print("Receptor:", receptor)
print("Ligandos validados: %d/3" % len(rows))
for row in rows:
    print("%s SOD1 min=%.4f median=%.4f semillas=%s" % (
        row["ligand"], row["vina_affinity_min"], row["vina_affinity_median"], row["seed_affinities"]
    ))
print("CSV:", csv_path)
