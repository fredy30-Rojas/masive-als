#!/usr/bin/env python3
"""
PIPELINE CONTINUO 24/7 - Oracle Cloud Free Tier
=================================================
- Ejecuta cribado virtual sin parar
- Guarda resultados incrementales
- Disenado para ARM64 (Oracle Ampere A1)
- 4 OCPU, 24 GB RAM, 200 GB disco
"""
import os, sys, csv, subprocess, time, json
from pathlib import Path
from datetime import datetime

PROJECT_ROOT = "/opt/masive-als"
COMPOUNDS_DIR = os.path.join(PROJECT_ROOT, "compounds")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "oracle_run")
PROTEINS_DIR = os.path.join(PROJECT_ROOT, "proteins")

TARGETS = ["TDP43", "SOD1", "FUS"]
BINDING_SITES = {
    "TDP43": {"center": (12.0, 8.0, 15.0), "size": (25, 25, 25)},
    "SOD1":  {"center": (2.0, 15.0, -4.0), "size": (25, 25, 25)},
    "FUS":   {"center": (5.0, 3.0, 8.0),  "size": (25, 25, 25)},
}

os.makedirs(RESULTS_DIR, exist_ok=True)
LOG_FILE = os.path.join(RESULTS_DIR, "pipeline.log")

def log(msg):
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {msg}"
    print(line)
    with open(LOG_FILE, "a") as f:
        f.write(line + "\n")

def generate_3d(smiles, output_pdbqt):
    """Genera estructura 3D - usa obabel CLI en Linux."""
    try:
        result = subprocess.run(
            ["obabel", f"-:{smiles}", "-O", output_pdbqt, "--gen3d", "-h"],
            capture_output=True, text=True, timeout=30
        )
        return result.returncode == 0 and os.path.exists(output_pdbqt)
    except:
        return False

def dock_single(lig_path, target, prot_path, lig_name):
    """Docking individual con Vina."""
    site = BINDING_SITES[target]
    out_dir = os.path.join(RESULTS_DIR, f"{target}_docked")
    os.makedirs(out_dir, exist_ok=True)
    out_pdbqt = os.path.join(out_dir, f"{lig_name}.pdbqt")
    
    if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 200:
        return None
    
    cmd = [
        "vina",
        "--receptor", prot_path,
        "--ligand", lig_path,
        "--out", out_pdbqt,
        "--center_x", str(site["center"][0]),
        "--center_y", str(site["center"][1]),
        "--center_z", str(site["center"][2]),
        "--size_x", str(site["size"][0]),
        "--size_y", str(site["size"][1]),
        "--size_z", str(site["size"][2]),
        "--exhaustiveness", "4",
        "--num_modes", "5",
    ]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        # Parse energy
        energy = None
        in_table = False
        for line in result.stdout.split("\n"):
            s = line.strip()
            if "mode" in s and "affinity" in s:
                in_table = True; continue
            if not in_table or not s or s.startswith("---"): continue
            parts = s.split()
            if len(parts) >= 2 and parts[0].isdigit():
                try:
                    if int(parts[0]) >= 2:
                        energy = float(parts[1]); break
                    if energy is None:
                        energy = float(parts[1])
                except ValueError:
                    pass
        return energy
    except subprocess.TimeoutExpired:
        return "TIMEOUT"
    except Exception:
        return "ERROR"

def main():
    log("=" * 50)
    log("PIPELINE MASIVE-ALS - Oracle Cloud 24/7")
    log(f"Inicio: {datetime.now()}")
    log("=" * 50)
    
    # Cargar compuestos
    csv_path = os.path.join(COMPOUNDS_DIR, "full_library.csv")
    if not os.path.exists(csv_path):
        log("ERROR: No hay libreria de compuestos")
        sys.exit(1)
    
    with open(csv_path) as f:
        compounds = list(csv.DictReader(f))
    
    log(f"Compuestos: {len(compounds)}")
    log(f"Proteinas: {len(TARGETS)}")
    
    # Preparar proteinas
    prot_paths = {}
    for t in TARGETS:
        pdbqt = os.path.join(RESULTS_DIR, f"{t}.pdbqt")
        if not os.path.exists(pdbqt):
            src = list(Path(os.path.join(PROTEINS_DIR, t)).glob("*.pdb"))
            if src:
                subprocess.run(["obabel", "-ipdb", str(src[0]), "-opdbqt", "-O", pdbqt, "-xr"],
                             capture_output=True, timeout=30)
        if os.path.exists(pdbqt):
            prot_paths[t] = pdbqt
    
    # Generar 3D
    ligands_3d = {}
    lig_dir = os.path.join(RESULTS_DIR, "ligands_3d")
    os.makedirs(lig_dir, exist_ok=True)
    
    for comp in compounds:
        name = comp.get("name", "")
        smiles = comp.get("smiles", "")
        if not smiles or smiles == "N/A":
            continue
        lp = os.path.join(lig_dir, f"{name}.pdbqt")
        if os.path.exists(lp):
            ligands_3d[name] = lp
        elif generate_3d(smiles, lp):
            ligands_3d[name] = lp
    
    log(f"Ligandos 3D: {len(ligands_3d)}")
    
    # Docking
    results_csv = os.path.join(RESULTS_DIR, "all_results.csv")
    results = []
    processed = set()
    
    if os.path.exists(results_csv):
        with open(results_csv) as f:
            for r in csv.DictReader(f):
                results.append(r)
                processed.add((r["ligand"], r["target"]))
    
    total = len(ligands_3d) * len(TARGETS)
    done = len(processed)
    t_start = time.time()
    
    for lig_name, lig_path in ligands_3d.items():
        for target in TARGETS:
            if (lig_name, target) in processed:
                continue
            
            prot = prot_paths.get(target)
            if not prot:
                continue
            
            energy = dock_single(lig_path, target, prot, lig_name)
            done += 1
            
            if energy is not None and energy != "TIMEOUT" and energy != "ERROR":
                results.append({"ligand": lig_name, "target": target, "energy": energy})
            
            if done % 20 == 0 or done == total:
                elapsed = time.time() - t_start
                rate = done / (elapsed / 60) if elapsed > 0 else 0
                eta = (total - done) / rate if rate > 0 else 0
                
                results.sort(key=lambda x: float(x["energy"]))
                with open(results_csv, "w") as f:
                    w = csv.DictWriter(f, fieldnames=["ligand", "target", "energy"])
                    w.writeheader()
                    w.writerows(results)
                
                log(f"[{done}/{total}] {done*100//total}% | {rate:.1f} cpm | ETA: {eta:.0f}min")
    
    log(f"COMPLETADO: {len(results)} resultados en {done} dockings")
    log(f"Tiempo: {(time.time()-t_start)/60:.1f} min")

if __name__ == "__main__":
    while True:
        try:
            main()
        except Exception as e:
            log(f"ERROR: {e}")
            log("Reiniciando en 60s...")
            time.sleep(60)
