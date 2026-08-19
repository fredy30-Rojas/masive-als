#!/usr/bin/env python3
"""MASIVE-ALS: Pipeline 24/7 - ChEMBL 2.4M compuestos, 4 cores paralelo"""
import csv, os, subprocess, time, gzip, sys, shutil
from datetime import datetime
from pathlib import Path

BASE = "/opt/masive-als"
CHEMBL = f"{BASE}/compounds/chembl_34_chemreps.txt.gz"
RES = f"{BASE}/results/massive_run"
PROT = f"{BASE}/proteins"
LOG = f"{RES}/pipeline.log"

TARGETS = ["TDP43", "SOD1", "FUS"]
SITES = {
    "TDP43": (12,8,15,25,25,25),
    "SOD1": (2,15,-4,25,25,25),
    "FUS": (5,3,8,25,25,25),
}

os.makedirs(RES, exist_ok=True)

def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG, "a") as f:
        f.write(line + "\n")

def run_one_dock(smiles, chembl_id, target, prot_pdbqt):
    """Ejecuta un docking completo: SMILES -> 3D -> Vina."""
    lig_dir = f"{RES}/ligands"
    out_dir = f"{RES}/{target}_docked"
    os.makedirs(lig_dir, exist_ok=True)
    os.makedirs(out_dir, exist_ok=True)
    
    lig_pdbqt = f"{lig_dir}/{chembl_id}.pdbqt"
    out_pdbqt = f"{out_dir}/{chembl_id}.pdbqt"
    
    # Saltar si ya existe
    if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 200:
        return None
    
    # Generar 3D
    if not os.path.exists(lig_pdbqt) or os.path.getsize(lig_pdbqt) < 100:
        try:
            r = subprocess.run(["obabel", "-ismi", "-opdbqt", "-O", lig_pdbqt, "--gen3d"],
                              input=smiles, capture_output=True, text=True, timeout=30)
            if not os.path.exists(lig_pdbqt) or os.path.getsize(lig_pdbqt) < 100:
                return None
        except:
            return None
    
    # Docking
    cx,cy,cz,sx,sy,sz = SITES[target]
    try:
        r = subprocess.run([
            "vina", "--receptor", prot_pdbqt, "--ligand", lig_pdbqt,
            "--out", out_pdbqt,
            "--center_x", str(cx), "--center_y", str(cy), "--center_z", str(cz),
            "--size_x", str(sx), "--size_y", str(sy), "--size_z", str(sz),
            "--exhaustiveness", "4", "--num_modes", "5"
        ], capture_output=True, text=True, timeout=120)
        
        # Parse energy
        energy = None
        in_table = False
        for line in r.stdout.split("\n"):
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
                except: pass
        
        return energy
    except:
        return None

def main():
    log("=" * 50)
    log("MASIVE-ALS: Pipeline ChEMBL 2.4M")
    log(f"Inicio: {datetime.now()}")
    
    # Preparar proteinas
    prots = {}
    for t in TARGETS:
        pq = f"{RES}/{t}.pdbqt"
        if not os.path.exists(pq):
            src = list(Path(f"{PROT}/{t}").glob("*.pdb"))[0]
            subprocess.run(["obabel","-ipdb",str(src),"-opdbqt","-O",pq,"-xr"], capture_output=True)
        if os.path.exists(pq):
            prots[t] = pq
    log(f"Proteinas: {len(prots)}")
    
    # Cargar resultados previos
    res_csv = f"{RES}/all_results.csv"
    results = []
    done_set = set()
    if os.path.exists(res_csv):
        with open(res_csv) as f:
            for r in csv.DictReader(f):
                results.append(r)
                done_set.add((r["ligand"], r["target"]))
    log(f"Resultados previos: {len(results)}")
    
    # Contar compuestos en ChEMBL
    log("Contando compuestos en ChEMBL...")
    total_compounds = 0
    with gzip.open(CHEMBL, "rt", encoding="utf-8", errors="replace") as f:
        next(f)  # skip header
        for _ in f:
            total_compounds += 1
    log(f"Compuestos totales: {total_compounds:,}")
    
    # Procesar
    done = 0
    skipped = 0
    t0 = time.time()
    batch = []
    
    with gzip.open(CHEMBL, "rt", encoding="utf-8", errors="replace") as f:
        next(f)  # skip header
        
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue
            chembl_id = parts[0]
            smiles = parts[1]
            
            if not smiles or len(smiles) < 5 or len(smiles) > 500:
                continue
            
            # Procesar contra cada proteina
            for target in TARGETS:
                if (chembl_id, target) in done_set:
                    skipped += 1
                    continue
                
                energy = run_one_dock(smiles, chembl_id, target, prots[target])
                done += 1
                
                if energy is not None:
                    results.append({"ligand": chembl_id, "target": target, "energy": energy})
                
                # Guardar cada 100 dockings
                if done % 100 == 0:
                    results.sort(key=lambda x: float(x["energy"]))
                    with open(res_csv, "w") as fout:
                        w = csv.DictWriter(fout, fieldnames=["ligand","target","energy"])
                        w.writeheader()
                        w.writerows(results)
                    
                    elapsed = time.time() - t0
                    rate = done / (elapsed / 60) if elapsed > 0 else 0
                    eta_total = (total_compounds * 3 - done) / rate if rate > 0 else 0
                    top = results[:3] if results else []
                    if top:
                        parts = []
                        for t2 in top:
                            parts.append(f"{t2['ligand']}={float(t2['energy']):.2f}")
                        log(f"{done:,} docks | {rate:.1f}/min | ETA: {eta_total/60/24:.0f}d | Top: {' | '.join(parts)}")
    
    log(f"COMPLETADO: {len(results)} resultados en {done:,} dockings")

if __name__ == "__main__":
    main()
