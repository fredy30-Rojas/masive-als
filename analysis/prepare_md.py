#!/usr/bin/env python3
"""Preparar archivos de entrada para validacion MD - MASIVE-ALS"""

import argparse, csv, os, subprocess, sys

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--hits", required=True)
    p.add_argument("--batch", type=int, required=True)
    p.add_argument("--size", type=int, default=50)
    p.add_argument("--output", required=True)
    return p.parse_args()

def prepare_system(ligand_id, target, output_dir):
    """Crear archivos de entrada para MD de un complejo"""
    system_dir = os.path.join(output_dir, f"system_{ligand_id}_{target}")
    os.makedirs(system_dir, exist_ok=True)
    
    # Crear parametros MDP para GROMACS
    mdp_templates = {
        "minim.mdp": """integrator = steep
nsteps = 50000
emtol = 100.0
emstep = 0.01
cutoff-scheme = Verlet
nstlist = 10
rlist = 1.2
coulombtype = PME
rcoulomb = 1.2
vdw-type = Cut-off
rvdw = 1.2
pbc = xyz""",
        
        "nvt.mdp": """integrator = md
nsteps = 50000
dt = 0.002
nstxout = 1000
nstvout = 1000
nstenergy = 1000
tcoupl = V-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1 0.1
ref-t = 300 300
pcoupl = no
pbc = xyz
gen-vel = yes
gen-temp = 300""",
        
        "npt.mdp": """integrator = md
nsteps = 50000
dt = 0.002
nstxout = 1000
nstvout = 1000
nstenergy = 1000
tcoupl = V-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1 0.1
ref-t = 300 300
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 2.0
ref-p = 1.0
compressibility = 4.5e-5
pbc = xyz""",
        
        "md.mdp": """integrator = md
nsteps = 500000000
dt = 0.002
nstxout = 10000
nstvout = 10000
nstenergy = 10000
tcoupl = V-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1 0.1
ref-t = 300 300
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 2.0
ref-p = 1.0
compressibility = 4.5e-5
pbc = xyz
constraints = h-bonds
constraint-algorithm = LINCS"""
    }
    
    for fname, content in mdp_templates.items():
        with open(os.path.join(system_dir, fname), "w") as f:
            f.write(content)
    
    return system_dir

def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)
    
    with open(args.hits) as f:
        reader = csv.DictReader(f)
        all_hits = list(reader)
    
    start = (args.batch - 1) * args.size
    batch_hits = all_hits[start:start + args.size]
    
    print(f"Preparando {len(batch_hits)} sistemas MD (batch {args.batch})...")
    
    for i, hit in enumerate(batch_hits):
        ligand = hit.get("ligand", f"LIG{i:04d}")
        target = hit.get("target", "TDP43")
        system_dir = prepare_system(ligand, target, args.output)
        print(f"  [{i+1}/{len(batch_hits)}] {system_dir}")

if __name__ == "__main__":
    main()
