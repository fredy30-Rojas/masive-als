#!/bin/bash
#==============================================================================
# Preparacion del entorno en MareNostrum 5 - MASIVE-ALS
# Ejecutar una vez al recibir el acceso
#==============================================================================

set -e

echo "============================================"
echo " MASIVE-ALS: Preparacion del entorno"
echo " Fecha: $(date)"
echo "============================================"

PROJECT_DIR="/gpfs/projects/masive-als"

# 1. Estructura de directorios
echo ">>> Creando estructura de directorios..."
mkdir -p ${PROJECT_DIR}/{compounds,proteins,results,md_results,logs,scripts,analysis,references}
echo "OK"

# 2. Cargar modulos disponibles
echo ">>> Verificando modulos..."
module purge
module load autodock-gpu/1.6 2>/dev/null && echo "  AutoDock-GPU: OK" || echo "  AutoDock-GPU: PENDIENTE"
module load gromacs/2024.3-gpu 2>/dev/null && echo "  GROMACS: OK" || echo "  GROMACS: PENDIENTE"
module load python/3.11 2>/dev/null && echo "  Python 3.11: OK" || echo "  Python 3.11: PENDIENTE"
module load alphafold/2.3.2 2>/dev/null && echo "  AlphaFold: OK" || echo "  AlphaFold: PENDIENTE (requiere instalacion)"
module load openbabel/3.1.1 2>/dev/null && echo "  OpenBabel: OK" || echo "  OpenBabel: PENDIENTE"

# 3. Copiar scripts al proyecto (asumiendo que se ejecuta desde el repo clonado)
echo ">>> Copiando scripts..."
SCRIPT_SRC="$(dirname "$(readlink -f "$0")")"
cp "${SCRIPT_SRC}/scripts/slurm_docking.sh" "${PROJECT_DIR}/scripts/" 2>/dev/null || echo "  slurm_docking.sh: copiar manualmente"
cp "${SCRIPT_SRC}/scripts/slurm_gromacs.sh" "${PROJECT_DIR}/scripts/" 2>/dev/null || echo "  slurm_gromacs.sh: copiar manualmente"
cp -r "${SCRIPT_SRC}/analysis/"* "${PROJECT_DIR}/analysis/" 2>/dev/null || echo "  analysis/: copiar manualmente"
chmod +x ${PROJECT_DIR}/scripts/*.sh 2>/dev/null || true

# 4. Crear entorno virtual Python
echo ">>> Configurando Python..."
python3 -m venv ${PROJECT_DIR}/venv
source ${PROJECT_DIR}/venv/bin/activate
pip install --upgrade pip -q
pip install numpy pandas scipy matplotlib seaborn rdkit biopython -q 2>/dev/null || pip install numpy pandas scipy matplotlib seaborn rdkit-pypi biopython -q
echo "Python OK"

# 5. Verificar espacio
echo ">>> Espacio disponible:"
df -h /gpfs/projects/masive-als

echo ""
echo "============================================"
echo " Entorno preparado. Proyecto listo en:"
echo " ${PROJECT_DIR}"
echo "============================================"
echo ""
echo "Siguientes pasos:"
echo "  1. Ejecutar: bash prep/download_compounds.sh"
echo "  2. Ejecutar: bash prep/prepare_proteins.sh"  
echo "  3. Lanzar cribado: sbatch scripts/slurm_docking.sh"
