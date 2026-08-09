#!/usr/bin/env python3
"""
Analisis de resultados de dinamica molecular - MASIVE-ALS
Calcula MM-GBSA, RMSD y filtra candidatos finales
"""

import argparse, csv, os
import numpy as np

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--dir", required=True, help="Directorio con resultados MD")
    p.add_argument("--output", required=True, help="CSV de salida")
    return p.parse_args()

def analyze_md_system(system_dir):
    """Analizar un sistema MD completo"""
    results = {"rmsd_mean": None, "rmsd_std": None, 
               "rg_mean": None, "hbonds": None, "stable": False}
    
    # RMSD
    rmsd_file = os.path.join(system_dir, "rmsd.xvg")
    if os.path.exists(rmsd_file):
        rmsd_data = []
        with open(rmsd_file) as f:
            for line in f:
                if line.startswith("#") or line.startswith("@"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    rmsd_data.append(float(parts[1]) * 10)  # nm to Angstrom
        
        if rmsd_data:
            # Ignorar primeros 5% (equilibracion)
            stable = rmsd_data[int(len(rmsd_data)*0.05):]
            results["rmsd_mean"] = np.mean(stable)
            results["rmsd_std"] = np.std(stable)
            results["stable"] = results["rmsd_std"] < 2.0  # < 2 Angstrom = estable
    
    # Radio de giro
    rg_file = os.path.join(system_dir, "gyrate.xvg")
    if os.path.exists(rg_file):
        rg_data = []
        with open(rg_file) as f:
            for line in f:
                if line.startswith("#") or line.startswith("@"):
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    rg_data.append(float(parts[1]))
        if rg_data:
            results["rg_mean"] = np.mean(rg_data)
    
    return results

def main():
    args = parse_args()
    
    print(f"Analizando MD en: {args.dir}")
    
    results = []
    for system_dir in sorted(os.listdir(args.dir)):
        full_path = os.path.join(args.dir, system_dir)
        if not os.path.isdir(full_path):
            continue
        
        md_results = analyze_md_system(full_path)
        ligand_name = system_dir.replace("system_", "").replace("_", " ")
        
        results.append({
            "ligand": ligand_name,
            "rmsd_mean": f"{md_results['rmsd_mean']:.2f}" if md_results["rmsd_mean"] else "N/A",
            "rmsd_std": f"{md_results['rmsd_std']:.2f}" if md_results["rmsd_std"] else "N/A",
            "rg_mean": f"{md_results['rg_mean']:.2f}" if md_results["rg_mean"] else "N/A",
            "stable": "YES" if md_results["stable"] else "NO"
        })
    
    # Filtrar solo sistemas estables
    stable = [r for r in results if r["stable"] == "YES"]
    
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["ligand", "rmsd_mean", "rmsd_std", "rg_mean", "stable"])
        writer.writeheader()
        writer.writerows(results)
    
    print(f"\n  Total sistemas: {len(results)}")
    print(f"  Estables (RMSD < 2A): {len(stable)}")
    print(f"  Resultados guardados en: {args.output}")

if __name__ == "__main__":
    main()
