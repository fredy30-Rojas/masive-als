#!/bin/bash
#==============================================================================
# SLURM: Cribado masivo AutoDock-GPU - MASIVE-ALS
# MareNostrum 5 - Particion ACC (NVIDIA Hopper H100)
#==============================================================================
#SBATCH --job-name=masive-als-dock
#SBATCH --partition=acc
#SBATCH --qos=acc_res
#SBATCH --nodes=50
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
#SBATCH --time=48:00:00
#SBATCH --output=logs/dock_%A_%a.out
#SBATCH --error=logs/dock_%A_%a.err
#SBATCH --array=1-50
#SBATCH --account=masive-als

# Total: 50 nodos x 4 GPUs = 200 GPUs simultaneas
# Cada array job procesa 200,000 compuestos (~5 horas)

set -e

echo "============================================"
echo " MASIVE-ALS: Cribado molecular - Tarea ${SLURM_ARRAY_TASK_ID}"
echo " Nodo: $(hostname)"
echo " Fecha: $(date)"
echo "============================================"

# Cargar modulos MareNostrum 5
module purge
module load autodock-gpu/1.6
module load openbabel/3.1.1
module load python/3.11

# Directorios
WORKDIR="/gpfs/projects/masive-als"
COMPOUNDS="${WORKDIR}/compounds"
PROTEINS="${WORKDIR}/proteins"
RESULTS="${WORKDIR}/results/dock_${SLURM_ARRAY_TASK_ID}"

mkdir -p "${RESULTS}"

# Calcular rango de compuestos para esta tarea
CHUNK_SIZE=200000
START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))

echo "Procesando compuestos ${START} a ${END}"

# Para cada proteina diana
for target in TDP43 SOD1 FUS; do
    echo "--- Diana: ${target} ---"
    
    PROTEIN_PDB="${PROTEINS}/${target}/conformations"
    
    # Seleccionar las 500 conformaciones mas representativas
    for conf in $(ls ${PROTEIN_PDB}/*.pdbqt | head -500); do
        conf_name=$(basename ${conf} .pdbqt)
        
        # Ejecutar AutoDock-GPU en las 4 GPUs del nodo
        srun --ntasks=1 --gres=gpu:1 autodock_gpu \
            -lfile ${COMPOUNDS}/ligands_${target}.pdbqt \
            -ffile ${conf} \
            -nrun 10 \
            -ngen 27000 \
            -npdb 10 \
            -resnam ${RESULTS}/${target}_${conf_name}_${SLURM_ARRAY_TASK_ID} &
    done
    wait
done

# Consolidar resultados
python3 ${WORKDIR}/analysis/merge_results.py \
    --input ${RESULTS} \
    --output ${WORKDIR}/results/merged_${SLURM_ARRAY_TASK_ID}.csv

echo "Tarea ${SLURM_ARRAY_TASK_ID} completada: $(date)"
