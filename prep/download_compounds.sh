#!/bin/bash
#==============================================================================
# Descarga de librerias de compuestos - MASIVE-ALS
# Fuentes: ZINC20, DrugBank, Enamine REAL
#==============================================================================

set -e

WORKDIR="/gpfs/projects/masive-als"
COMPOUNDS="${WORKDIR}/compounds"

echo ">>> Descargando librerias de compuestos..."

# 1. ZINC20 - compuestos drug-like (gratuito, ~5M compuestos)
echo ">>> ZINC20: descargando compuestos drug-like..."
mkdir -p ${COMPOUNDS}/zinc20
# ZINC20 se descarga por tranches desde zinc20.docking.org
# Usamos wget con patrones correctos para cada tranche
for tranche in {A..Z}; do
    for subtranche in {A..Z}; do
        url="https://zinc20.docking.org/tranches/${tranche}${subtranche}/${tranche}${subtranche}.mol2.gz"
        wget -q -P "${COMPOUNDS}/zinc20" \
            "${url}" \
            -A "*.mol2.gz" \
            --limit-rate=100M \
            --timeout=30 \
            --tries=2 \
            --continue 2>/dev/null &
        
        # Limitar a 10 descargas simultaneas
        if (( $(jobs -r | wc -l) >= 10 )); then
            wait -n 2>/dev/null || true
        fi
    done
done
wait
echo "ZINC20: descarga completada ($(ls ${COMPOUNDS}/zinc20/*.mol2.gz 2>/dev/null | wc -l) archivos)"

# 2. DrugBank - farmacos aprobados y experimentales (~15K compuestos)
echo ">>> DrugBank: descargando..."
# DrugBank requiere registro gratuito. Usar XML descargado.
if [ -f "${COMPOUNDS}/drugbank_all_structures.sdf" ]; then
    echo "DrugBank: ya descargado"
else
    wget -q -P ${COMPOUNDS} \
        "https://go.drugbank.com/structures/small_molecule_drugs.sdf" \
        --continue
    mv ${COMPOUNDS}/small_molecule_drugs.sdf ${COMPOUNDS}/drugbank_all_structures.sdf
    echo "DrugBank: OK"
fi

# 3. Enamine REAL - compuestos sintetizables (~2M drug-like)
echo ">>> Enamine REAL: descargando..."
# Enamine ofrece acceso gratuito a su libreria REAL para docking
wget -q -P ${COMPOUNDS}/enamine \
    "https://enamine.net/download/real-database" \
    --continue 2>/dev/null || echo "Enamine: Requiere registro manual en enamine.net"

# 4. Convertir todo a PDBQT para AutoDock
echo ">>> Convirtiendo a formato PDBQT..."
module load openbabel/3.1.1

for f in ${COMPOUNDS}/zinc20/*.mol2.gz; do
    zcat "$f" | obabel -imol2 -opdbqt -O ${COMPOUNDS}/pdbqt/$(basename ${f%.mol2.gz}).pdbqt &
done
wait

echo ">>> Librerias descargadas:"
echo "  ZINC20: $(ls ${COMPOUNDS}/zinc20/*.mol2.gz 2>/dev/null | wc -l) tranches"
echo "  DrugBank: $(test -f ${COMPOUNDS}/drugbank_all_structures.sdf && echo 'OK' || echo 'PENDIENTE')"
echo "  Enamine: $(ls ${COMPOUNDS}/enamine/ 2>/dev/null | wc -l) archivos"
echo "  PDBQT: $(ls ${COMPOUNDS}/pdbqt/*.pdbqt 2>/dev/null | wc -l) ligandos"

echo "Completado: $(date)"
