#!/bin/bash
#==============================================================================
# SLURM: Cribado masivo AutoDock-GPU - MASIVE-ALS
# MareNostrum 5 - Particion ACC (NVIDIA Hopper H100)
#==============================================================================
#SBATCH --job-name=masive-als-dock
#SBATCH --partition=acc
#SBATCH --qos=acc_res
#SBATCH --nodes=50
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:4
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/projects/masive-als/logs/dock_%A_%a.out
#SBATCH --error=/gpfs/projects/masive-als/logs/dock_%A_%a.err
#SBATCH --array=1-50
#SBATCH --account=masive-als

# Total: 50 jobs x 4 GPUs = 200 GPUs simultaneas
# Cada array job procesa 200,000 compuestos

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
mkdir -p "${WORKDIR}/logs"

# Calcular rango de compuestos para esta tarea
CHUNK_SIZE=200000
START=$(( (SLURM_ARRAY_TASK_ID - 1) * CHUNK_SIZE + 1 ))
END=$(( SLURM_ARRAY_TASK_ID * CHUNK_SIZE ))

echo "Procesando compuestos ${START} a ${END}"

# Para cada proteina diana
for target in TDP43 SOD1 FUS; do
    echo "--- Diana: ${target} ---"
    
    PROTEIN_PDB="${PROTEINS}/${target}/conformations"
    LIGAND_FILE="${COMPOUNDS}/${target}_batch_${SLURM_ARRAY_TASK_ID}.pdbqt"
    
    if [ ! -f "${LIGAND_FILE}" ]; then
        echo "ERROR: No se encuentra ${LIGAND_FILE}. Saltando ${target}."
        continue
    fi
    
    # Seleccionar 100 conformaciones mas representativas
    conf_count=0
    for conf in $(ls ${PROTEIN_PDB}/*.pdbqt 2>/dev/null | head -100); do
        conf_name=$(basename ${conf} .pdbqt)
        conf_count=$((conf_count + 1))
        
        # AutoDock-GPU: cada instancia usa 1 GPU
        CUDA_VISIBLE_DEVICES=$(( (conf_count - 1) % 4 )) autodock_gpu \
            -lfile "${LIGAND_FILE}" \
            -ffile "${conf}" \
            -nrun 10 \
            -ngen 27000 \
            -npdb 10 \
            -resnam "${RESULTS}/${target}_${conf_name}_${SLURM_ARRAY_TASK_ID}" &
        
        # Limitar a 4 procesos simultaneos (una por GPU)
        if (( conf_count % 4 == 0 )); then
            wait
        fi
    done
    wait
done

# Consolidar resultados
python3 ${WORKDIR}/analysis/merge_results.py \
    --input "${RESULTS}" \
    --output "${WORKDIR}/results/merged_${SLURM_ARRAY_TASK_ID}.csv"

echo "Tarea ${SLURM_ARRAY_TASK_ID} completada: $(date)"
