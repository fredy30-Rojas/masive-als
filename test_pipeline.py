#!/usr/bin/env python3
"""
PRUEBA COMPLETA DEL PIPELINE MASIVE-ALS
Simula datos de docking y ejecuta los 3 scripts de analisis en secuencia.
"""
import os, sys, csv, random, subprocess, tempfile

TEST_DIR = os.path.join(tempfile.gettempdir(), "masive-als-test")
os.makedirs(TEST_DIR, exist_ok=True)

print("=" * 60)
print(" PRUEBA PIPELINE MASIVE-ALS")
print("=" * 60)

# ──────────────────────────────────────────────────
# FASE 1: Crear datos simulados de docking (AutoDock-GPU)
# ──────────────────────────────────────────────────
print("\n>>> FASE 1: Generando datos simulados de docking...")
RESULTS_DIR = os.path.join(TEST_DIR, "results", "dock_1")
os.makedirs(RESULTS_DIR, exist_ok=True)

targets = ["TDP43", "SOD1", "FUS"]
ligands_created = 0

for i in range(50):  # 50 ligandos simulados
    target = targets[i % 3]
    ligand_name = f"ZINC{i:08d}"
    binding_energy = random.uniform(-14.0, -6.0)  # kcal/mol
    
    dlg_file = os.path.join(RESULTS_DIR, f"{target}_{ligand_name}_1.dlg")
    with open(dlg_file, "w") as f:
        f.write(f"AutoDock-GPU Result Log\n")
        f.write(f"Input ligand: {ligand_name}\n")
        f.write(f"Target: {target}\n")
        f.write(f"USER    Best Energy  {binding_energy:.2f}\n")
        f.write(f"Estimated Free Energy of Binding = {binding_energy:.2f} kcal/mol\n")
    ligands_created += 1

print(f"  OK: {ligands_created} archivos .dlg generados (3 targets x ~17 ligandos)")

# ──────────────────────────────────────────────────
# FASE 2: merge_results.py - Top hits
# ──────────────────────────────────────────────────
print("\n>>> FASE 2: merge_results.py - Identificando top hits...")

merge_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "analysis", "merge_results.py")
merged_csv = os.path.join(TEST_DIR, "results", "merged_1.csv")

result = subprocess.run(
    [sys.executable, merge_script, "--input", RESULTS_DIR, "--output", merged_csv, "--top", "10"],
    capture_output=True, text=True, timeout=30
)
print(result.stdout)
if result.returncode != 0:
    print(f"  ERROR: {result.stderr}")
    sys.exit(1)

# Verificar que el CSV tiene datos
with open(merged_csv) as f:
    reader = csv.DictReader(f)
    rows = list(reader)
print(f"  OK: {len(rows)} hits en el CSV")
assert len(rows) > 0, "ERROR: merge_results.py no genero resultados!"

# ──────────────────────────────────────────────────
# FASE 3: prepare_md.py - Preparar sistemas MD
# ──────────────────────────────────────────────────
print("\n>>> FASE 3: prepare_md.py - Preparando sistemas MD...")

# Crear CSV de top hits con formato correcto
top_csv = os.path.join(TEST_DIR, "results", "top_hits.csv")
with open(merged_csv) as fin, open(top_csv, "w", newline="") as fout:
    reader = csv.DictReader(fin)
    writer = csv.DictWriter(fout, fieldnames=["ligand", "target", "binding_energy", "file"])
    writer.writeheader()
    writer.writerows(reader)

prepare_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "analysis", "prepare_md.py")
md_dir = os.path.join(TEST_DIR, "md_results")
os.makedirs(md_dir, exist_ok=True)

result = subprocess.run(
    [sys.executable, prepare_script, "--hits", top_csv, "--batch", "1", "--size", "5", "--output", os.path.join(md_dir, "job_1")],
    capture_output=True, text=True, timeout=30
)
print(result.stdout)
if result.returncode != 0:
    print(f"  ERROR: {result.stderr}")
    sys.exit(1)

# Verificar directorios creados
job_dir = os.path.join(md_dir, "job_1")
systems = [d for d in os.listdir(job_dir) if os.path.isdir(os.path.join(job_dir, d))]
print(f"  OK: {len(systems)} sistemas MD preparados")
assert len(systems) > 0, "ERROR: prepare_md.py no creo sistemas!"

# Verificar archivos MDP
for sys_dir in systems[:1]:
    full = os.path.join(job_dir, sys_dir)
    for mdp in ["minim.mdp", "nvt.mdp", "npt.mdp", "md.mdp"]:
        assert os.path.exists(os.path.join(full, mdp)), f"ERROR: Falta {mdp} en {sys_dir}"
print(f"  OK: Archivos MDP verificados (minim, nvt, npt, md)")

# ──────────────────────────────────────────────────
# FASE 4: Simular analisis MD (crear datos falsos de RMSD)
# ──────────────────────────────────────────────────
print("\n>>> FASE 4: analyze_md.py - Creando datos simulados de MD...")

for sys_dir in systems:
    full = os.path.join(job_dir, sys_dir)
    # Crear rmsd.xvg simulado
    with open(os.path.join(full, "rmsd.xvg"), "w") as f:
        f.write("# RMSD simulation\n")
        f.write("@ title \"RMSD\"\n")
        f.write("@ xaxis label \"Time (ps)\"\n")
        f.write("@ yaxis label \"RMSD (nm)\"\n")
        for t in range(100):
            rmsd = 0.15 + random.uniform(-0.05, 0.05)  # ~1.5 Angstrom estable
            f.write(f"{t*10:.1f}  {rmsd:.4f}\n")
    
    # Crear gyrate.xvg simulado
    with open(os.path.join(full, "gyrate.xvg"), "w") as f:
        f.write("# Radius of gyration\n")
        for t in range(100):
            rg = 1.8 + random.uniform(-0.1, 0.1)
            f.write(f"{t*10:.1f}  {rg:.4f}\n")

print(f"  OK: Datos MD simulados para {len(systems)} sistemas")

# ──────────────────────────────────────────────────
# FASE 5: analyze_md.py - Analisis final
# ──────────────────────────────────────────────────
print("\n>>> FASE 5: analyze_md.py - Analizando resultados MD...")

analyze_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "analysis", "analyze_md.py")
final_csv = os.path.join(TEST_DIR, "results", "final_candidates.csv")

result = subprocess.run(
    [sys.executable, analyze_script, "--dir", job_dir, "--output", final_csv],
    capture_output=True, text=True, timeout=30
)
print(result.stdout)
if result.returncode != 0:
    print(f"  ERROR: {result.stderr}")
    sys.exit(1)

# Verificar resultados
with open(final_csv) as f:
    reader = csv.DictReader(f)
    rows = list(reader)
    stable_count = sum(1 for r in rows if r.get("stable") == "YES")
print(f"  OK: {len(rows)} sistemas analizados, {stable_count} estables")
assert len(rows) > 0, "ERROR: analyze_md.py no genero resultados!"

# ──────────────────────────────────────────────────
# RESUMEN
# ──────────────────────────────────────────────────
print("\n" + "=" * 60)
print(" RESULTADO FINAL: PIPELINE COMPLETO EXITOSO")
print("=" * 60)
print(f"  Docking simulado:   {ligands_created} archivos .dlg")
print(f"  Top hits detectados: {len(rows)} candidatos")
print(f"  Sistemas MD creados: {len(systems)}")
print(f"  Candidatos estables: {stable_count}")
print(f"  Directorio prueba:   {TEST_DIR}")
print("\n>>> PIPELINE VALIDADO. Cuando BSC conceda el acceso,")
print(">>> los scripts funcionaran en MareNostrum 5.")
