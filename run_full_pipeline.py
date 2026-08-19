#!/usr/bin/env python3
"""
PIPELINE COMPLETO CONTINUO - MASIVE-ALS
=========================================
Ejecuta el cribado virtual completo:
1. Genera estructuras 3D para TODOS los compuestos
2. Ejecuta docking contra TDP-43, SOD1, FUS
3. Ranking de hits por energia
4. Guarda resultados incrementales cada 10 compuestos
"""
import os, sys, csv, subprocess, time, json
from pathlib import Path
from datetime import datetime
from openbabel import openbabel as ob

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
COMPOUNDS_DIR = os.path.join(PROJECT_ROOT, "compounds")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "full_run")
TOOLS_DIR = os.path.join(PROJECT_ROOT, "tools")
VINA_EXE = os.path.join(TOOLS_DIR, "vina.exe")

# Encontrar obabel
import shutil as _shutil
def _find_obabel():
    for c in ["obabel", 
              os.path.join(os.path.dirname(sys.executable), "Scripts", "obabel.exe"),
              os.path.join(os.path.dirname(sys.executable), "Lib", "site-packages", "openbabel", "bin", "obabel.exe")]:
        if os.path.exists(c) or _shutil.which(c):
            return c if os.path.exists(c) else _shutil.which(c)
    return "obabel"
OBABEL = _find_obabel()

TARGETS = ["TDP43", "SOD1", "FUS"]
PROTEINS_DIR = os.path.join(PROJECT_ROOT, "proteins")

BINDING_SITES = {
    "TDP43": {"center": (12.0, 8.0, 15.0), "size": (25, 25, 25)},
    "SOD1":  {"center": (2.0, 15.0, -4.0), "size": (25, 25, 25)},
    "FUS":   {"center": (5.0, 3.0, 8.0),  "size": (25, 25, 25)},
}

os.makedirs(RESULTS_DIR, exist_ok=True)


def generate_3d_structure(smiles, output_pdbqt):
    """Genera estructura 3D con OpenBabel Python API."""
    try:
        conv = ob.OBConversion()
        conv.SetInFormat("smi")
        mol = ob.OBMol()
        if not conv.ReadString(mol, smiles):
            return False
        
        mol.AddHydrogens()
        builder = ob.OBBuilder()
        builder.Build(mol)
        
        ff = ob.OBForceField.FindForceField("UFF")
        if ff:
            ff.Setup(mol)
            ff.SteepestDescent(200)
            ff.ConjugateGradients(100)
            ff.GetCoordinates(mol)
        
        conv.SetOutFormat("pdbqt")
        conv.WriteFile(mol, output_pdbqt)
        return True
    except:
        return False


def dock_compound(ligand_path, target_name, protein_pdbqt, lig_name):
    """Ejecuta Vina para un compuesto contra una proteina."""
    site = BINDING_SITES[target_name]
    out_dir = os.path.join(RESULTS_DIR, f"{target_name}_docked")
    os.makedirs(out_dir, exist_ok=True)
    
    out_pdbqt = os.path.join(out_dir, f"{lig_name}.pdbqt")
    
    if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 200:
        return None  # Ya procesado
    
    cmd = [
        VINA_EXE,
        "--receptor", protein_pdbqt,
        "--ligand", ligand_path,
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
        
        # Extraer mejor energia del modo 2 (modo 1 es pose de entrada, afinidad 0)
        best_energy = None
        in_table = False
        for line in result.stdout.split("\n"):
            stripped = line.strip()
            if "mode" in stripped and "affinity" in stripped:
                in_table = True
                continue
            if not in_table or not stripped or stripped.startswith("---"):
                continue
            parts = stripped.split()
            if len(parts) >= 2 and parts[0].isdigit():
                try:
                    mode_num = int(parts[0])
                    energy = float(parts[1])
                    if mode_num >= 2:
                        best_energy = energy
                        break
                    if mode_num == 1 and best_energy is None:
                        best_energy = energy
                except (ValueError, IndexError):
                    pass
        
        return best_energy
    except subprocess.TimeoutExpired:
        return "TIMEOUT"
    except Exception:
        return "ERROR"


def main():
    print("=" * 60)
    print(" PIPELINE COMPLETO MASIVE-ALS")
    print(f" Inicio: {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 60)
    
    # Cargar compuestos
    csv_path = os.path.join(COMPOUNDS_DIR, "full_library.csv")
    if not os.path.exists(csv_path):
        csv_path = os.path.join(COMPOUNDS_DIR, "fda_subset.csv")
    
    with open(csv_path) as f:
        reader = csv.DictReader(f)
        compounds = list(reader)
    
    print(f"\n  Compuestos: {len(compounds)}")
    print(f"  Proteinas:  {len(TARGETS)} ({', '.join(TARGETS)})")
    print(f"  Total dockings: {len(compounds) * len(TARGETS)}")
    
    # Preparar proteinas
    print("\n  Preparando proteinas...")
    protein_paths = {}
    for target in TARGETS:
        pdbqt = os.path.join(RESULTS_DIR, f"{target}.pdbqt")
        if not os.path.exists(pdbqt):
            # Convertir desde proteins dir
            src_pdb = list(Path(os.path.join(PROTEINS_DIR, target)).glob("*.pdb"))[0]
            subprocess.run(
                [OBABEL, "-ipdb", str(src_pdb), "-opdbqt", "-O", pdbqt, "-xr"],
                capture_output=True, timeout=30
            )
        if os.path.exists(pdbqt):
            protein_paths[target] = pdbqt
            print(f"    {target}: {os.path.getsize(pdbqt)/1024:.0f} KB")
    
    # Generar 3D para todos los compuestos
    print("\n  Generando estructuras 3D...")
    ligands_3d = {}
    ligands_dir = os.path.join(RESULTS_DIR, "ligands_3d")
    os.makedirs(ligands_dir, exist_ok=True)
    
    for i, comp in enumerate(compounds):
        name = comp.get("name", f"Compuesto{i}")
        smiles = comp.get("smiles", "")
        if not smiles or smiles == "N/A":
            continue
        
        lig_path = os.path.join(ligands_dir, f"{name}.pdbqt")
        if os.path.exists(lig_path) and os.path.getsize(lig_path) > 200:
            ligands_3d[name] = lig_path
            continue
        
        if generate_3d_structure(smiles, lig_path):
            ligands_3d[name] = lig_path
    
    print(f"    {len(ligands_3d)} compuestos 3D generados")
    
    # Ejecutar docking
    print(f"\n  Iniciando docking...")
    print(f"  {'='*50}")
    
    results = []
    total = len(ligands_3d) * len(TARGETS)
    done = 0
    t_start = time.time()
    
    # Cargar resultados previos si existen
    results_csv = os.path.join(RESULTS_DIR, "all_results.csv")
    if os.path.exists(results_csv):
        with open(results_csv) as f:
            results = list(csv.DictReader(f))
        print(f"    Cargados {len(results)} resultados previos")
    
    processed = set((r["ligand"], r["target"]) for r in results)
    
    for lig_name, lig_path in ligands_3d.items():
        for target in TARGETS:
            if (lig_name, target) in processed:
                done += 1
                continue
            
            prot_path = protein_paths.get(target)
            if not prot_path:
                continue
            
            energy = dock_compound(lig_path, target, prot_path, lig_name)
            done += 1
            
            if energy is not None and energy != "TIMEOUT" and energy != "ERROR":
                results.append({
                    "ligand": lig_name,
                    "target": target,
                    "energy": energy,
                    "category": compounds[list(ligands_3d.keys()).index(lig_name)].get("category", "") if lig_name in ligands_3d else ""
                })
            
            # Guardar cada 10 dockings
            if done % 10 == 0 or done == total:
                elapsed = time.time() - t_start
                rate = done / (elapsed / 60) if elapsed > 0 else 0
                eta_min = (total - done) / rate if rate > 0 else 0
                
                # Guardar CSV
                results.sort(key=lambda x: float(x["energy"]))
                with open(results_csv, "w", newline="") as f:
                    writer = csv.DictWriter(f, fieldnames=["ligand", "target", "energy", "category"])
                    writer.writeheader()
                    writer.writerows(results)
                
                # Mostrar progreso
                top3 = results[:3] if results else []
                top_str = " | ".join(f"{r['ligand'][:12]}:{float(r['energy']):.1f}" for r in top3)
                print(f"    [{done}/{total}] {done*100//total}% | {rate:.1f} comp/min | ETA: {eta_min:.0f}min | Top: {top_str}")
    
    # Informe final
    elapsed = time.time() - t_start
    print(f"\n{'='*60}")
    print(f" PIPELINE COMPLETADO")
    print(f"{'='*60}")
    print(f"  Tiempo: {elapsed/60:.1f} min")
    print(f"  Dockings: {done}")
    print(f"  Resultados: {len(results)}")
    
    if results:
        results.sort(key=lambda x: float(x["energy"]))
        print(f"\n  TOP 10 CANDIDATOS:")
        print(f"  {'#':<4} {'Compuesto':<20} {'Target':<8} {'Energia':<10} {'Categoria'}")
        print(f"  {'-'*60}")
        for i, r in enumerate(results[:10]):
            print(f"  {i+1:<4} {r['ligand'][:20]:<20} {r['target']:<8} {float(r['energy']):>8.2f}   {r.get('category', '')[:15]}")
    
    print(f"\n  Resultados: {results_csv}")


if __name__ == "__main__":
    main()
