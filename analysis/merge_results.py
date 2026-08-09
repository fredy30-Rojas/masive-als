#!/usr/bin/env python3
"""
Analisis de resultados de docking - MASIVE-ALS
Identifica los mejores candidatos para validacion MD
"""

import argparse, csv, os
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, help="Directorio con resultados de docking")
    p.add_argument("--output", required=True, help="CSV con top hits")
    p.add_argument("--top", type=int, default=1000, help="Numero de hits a seleccionar")
    return p.parse_args()

def read_docking_scores(result_dir):
    """Extraer energias de union de los archivos de resultado"""
    scores = []
    
    for fname in sorted(os.listdir(result_dir)):
        if not fname.endswith(".dlg"):
            continue
        
        filepath = os.path.join(result_dir, fname)
        with open(filepath) as f:
            current_ligand = None
            best_score = float("inf")
            
            for line in f:
                if line.startswith("Input ligand:"):
                    current_ligand = line.split(":")[1].strip()
                
                # AutoDock4: energia de union estimada
                if "Estimated Free Energy of Binding" in line:
                    try:
                        score = float(line.split("=")[1].split()[0].strip())
                        if score < best_score:
                            best_score = score
                    except:
                        pass
                
                # AutoDock-GPU: mejor energia del cluster
                if "USER    Best Energy" in line:
                    try:
                        parts = line.split()
                        score = float(parts[-1])
                        if score < best_score:
                            best_score = score
                    except:
                        pass
            
            if current_ligand and best_score < float("inf"):
                scores.append({
                    "ligand": current_ligand,
                    "target": extract_target(fname),
                    "binding_energy": best_score,
                    "file": fname
                })
    
    return sorted(scores, key=lambda x: x["binding_energy"])

def extract_target(filename):
    """Extraer nombre de proteina diana del nombre de archivo"""
    for target in ["TDP43", "SOD1", "FUS"]:
        if target in filename:
            return target
    return "UNKNOWN"

def main():
    args = parse_args()
    
    print(f"Analizando resultados en: {args.input}")
    scores = read_docking_scores(args.input)
    print(f"  Total de docks analizados: {len(scores)}")
    
    # Seleccionar top N globales
    top_hits = scores[:args.top]
    
    # Guardar
    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["ligand", "target", "binding_energy", "file"])
        writer.writeheader()
        writer.writerows(top_hits)
    
    # Estadisticas por diana
    by_target = defaultdict(lambda: {"count": 0, "best": float("inf")})
    for h in top_hits:
        by_target[h["target"]]["count"] += 1
        by_target[h["target"]]["best"] = min(by_target[h["target"]]["best"], h["binding_energy"])
    
    print(f"\nTop {args.top} hits guardados en: {args.output}")
    print("\nDistribucion por diana:")
    for target, stats in sorted(by_target.items()):
        print(f"  {target}: {stats['count']} hits, mejor energia = {stats['best']:.2f} kcal/mol")

if __name__ == "__main__":
    main()
