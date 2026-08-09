#!/bin/bash
#==============================================================================
# SLURM: Dinamica molecular de validacion - MASIVE-ALS
# MareNostrum 5 - Particion ACC (NVIDIA Hopper H100)
#==============================================================================
#SBATCH --job-name=masive-als-md
#SBATCH --partition=acc
#SBATCH --qos=acc_res
#SBATCH --nodes=10
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --time=72:00:00
#SBATCH --output=logs/md_%A_%a.out
#SBATCH --error=logs/md_%A_%a.err
#SBATCH --array=1-20
#SBATCH --account=masive-als

# Cada array job valida 50 compuestos x 1 microsegundo de MD

set -e

echo "============================================"
echo " MASIVE-ALS: Validacion MD - Tarea ${SLURM_ARRAY_TASK_ID}"
echo " Fecha: $(date)"
echo "============================================"

module purge
module load gromacs/2024.3-gpu
module load python/3.11

WORKDIR="/gpfs/projects/masive-als"
TOPHITS="${WORKDIR}/results/top_hits.csv"
MDDIR="${WORKDIR}/md_results"

mkdir -p ${MDDIR}/job_${SLURM_ARRAY_TASK_ID}

# Leer top hits para esta tarea (50 compuestos)
BATCH_SIZE=50
python3 ${WORKDIR}/analysis/prepare_md.py \
    --hits ${TOPHITS} \
    --batch ${SLURM_ARRAY_TASK_ID} \
    --size ${BATCH_SIZE} \
    --output ${MDDIR}/job_${SLURM_ARRAY_TASK_ID}

# Para cada compuesto, ejecutar MD de 1 microsegundo
cd ${MDDIR}/job_${SLURM_ARRAY_TASK_ID}

for system in system_*; do
    cd ${system}
    
    # Minimizacion
    gmx grompp -f minim.mdp -c complex.gro -p topol.top -o minim.tpr
    gmx mdrun -deffnm minim -nb gpu -bonded gpu -pme gpu
    
    # Equilibracion NVT
    gmx grompp -f nvt.mdp -c minim.gro -p topol.top -o nvt.tpr
    gmx mdrun -deffnm nvt -nb gpu -bonded gpu -pme gpu
    
    # Equilibracion NPT
    gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr
    gmx mdrun -deffnm npt -nb gpu -bonded gpu -pme gpu
    
    # Produccion 1 microsegundo
    gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr
    gmx mdrun -deffnm md -nb gpu -bonded gpu -pme gpu -nsteps 500000000
    
    # Analisis MM-GBSA
    gmx rms -s md.tpr -f md.xtc -o rmsd.xvg
    gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg
    
    cd ..
done

# Consolidar resultados MD de esta tanda
python3 ${WORKDIR}/analysis/analyze_md.py \
    --dir ${MDDIR}/job_${SLURM_ARRAY_TASK_ID} \
    --output ${MDDIR}/md_results_${SLURM_ARRAY_TASK_ID}.csv

echo "Validacion MD ${SLURM_ARRAY_TASK_ID} completada: $(date)"
