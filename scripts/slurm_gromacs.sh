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
#SBATCH --output=/gpfs/projects/masive-als/logs/md_%A_%a.out
#SBATCH --error=/gpfs/projects/masive-als/logs/md_%A_%a.err
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
JOB_DIR="${MDDIR}/job_${SLURM_ARRAY_TASK_ID}"
cd "${JOB_DIR}" || { echo "ERROR: No se encuentra ${JOB_DIR}"; exit 1; }

for system_dir in system_*; do
    if [ ! -d "${system_dir}" ]; then continue; fi
    cd "${JOB_DIR}/${system_dir}" || continue
    
    echo "  MD: ${system_dir}"
    
    # Verificar que existen los archivos necesarios
    if [ ! -f "complex.gro" ] || [ ! -f "topol.top" ]; then
        echo "  WARNING: Faltan archivos en ${system_dir}, saltando"
        cd "${JOB_DIR}"
        continue
    fi
    
    # Minimizacion
    gmx grompp -f minim.mdp -c complex.gro -p topol.top -o minim.tpr -maxwarn 1
    gmx mdrun -deffnm minim -nb gpu -bonded gpu -pme gpu -ntmpi 1
    
    # Equilibracion NVT
    gmx grompp -f nvt.mdp -c minim.gro -p topol.top -o nvt.tpr -maxwarn 1
    gmx mdrun -deffnm nvt -nb gpu -bonded gpu -pme gpu -ntmpi 1
    
    # Equilibracion NPT
    gmx grompp -f npt.mdp -c nvt.gro -p topol.top -o npt.tpr -maxwarn 1
    gmx mdrun -deffnm npt -nb gpu -bonded gpu -pme gpu -ntmpi 1
    
    # Produccion 1 microsegundo (500M steps x 0.002 ps = 1 us)
    gmx grompp -f md.mdp -c npt.gro -p topol.top -o md.tpr -maxwarn 1
    gmx mdrun -deffnm md -nb gpu -bonded gpu -pme gpu -nsteps 500000000 -ntmpi 1
    
    # Analisis post-MD
    echo "0 0" | gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns 2>/dev/null || true
    echo "0 0" | gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg 2>/dev/null || true
    
    cd "${JOB_DIR}"
done

# Consolidar resultados MD de esta tanda
python3 ${WORKDIR}/analysis/analyze_md.py \
    --dir ${MDDIR}/job_${SLURM_ARRAY_TASK_ID} \
    --output ${MDDIR}/md_results_${SLURM_ARRAY_TASK_ID}.csv

echo "Validacion MD ${SLURM_ARRAY_TASK_ID} completada: $(date)"
