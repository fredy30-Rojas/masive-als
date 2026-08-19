#!/bin/bash
# ============================================================================
# MASIVE-ALS: Instalacion completa en Oracle Cloud Free Tier (ARM64/Ubuntu)
# Recursos gratis: 4 OCPU Ampere A1, 24 GB RAM, 200 GB disco
# ============================================================================
set -e

echo "============================================"
echo " MASIVE-ALS - Instalacion Oracle Cloud"
echo "============================================"
echo ""

# ──── 1. Actualizar sistema ────
echo "[1/6] Actualizando sistema..."
sudo apt-get update -qq
sudo apt-get upgrade -y -qq

# ──── 2. Instalar dependencias ────
echo "[2/6] Instalando dependencias..."
sudo apt-get install -y -qq \
    python3 python3-pip python3-venv \
    openbabel \
    wget curl git \
    build-essential cmake \
    libboost-all-dev \
    htop tmux

# ──── 3. Instalar AutoDock Vina ────
echo "[3/6] Instalando AutoDock Vina..."
cd /tmp
wget -q "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_aarch64" -O vina || true
if [ ! -f vina ] || [ $(stat -c%s vina) -lt 1000 ]; then
    echo "  Compilando Vina desde fuente..."
    sudo apt-get install -y -qq libeigen3-dev
    git clone --depth 1 https://github.com/ccsb-scripps/AutoDock-Vina.git /tmp/vina-src
    cd /tmp/vina-src/build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make -j$(nproc)
    sudo cp vina /usr/local/bin/
else
    chmod +x vina
    sudo mv vina /usr/local/bin/
fi
echo "  Vina: $(which vina)"

# ──── 4. Instalar Python packages ────
echo "[4/6] Instalando paquetes Python..."
python3 -m venv /opt/masive-als/venv
source /opt/masive-als/venv/bin/activate
pip install --upgrade pip -q
pip install numpy matplotlib biopython pdbfixer openmm -q

echo "  Python packages OK"

# ──── 5. Crear estructura del proyecto ────
echo "[5/6] Creando estructura MASIVE-ALS..."
sudo mkdir -p /opt/masive-als/{proteins,compounds,results,tools,src}
sudo chown -R $USER:$USER /opt/masive-als

# ──── 6. Copiar datos del proyecto ────
echo "[6/6] Descargando datos..."
cd /opt/masive-als

# Proteinas
mkdir -p proteins/{TDP43,SOD1,FUS}
wget -q "https://files.rcsb.org/download/6b1n.pdb" -O proteins/TDP43/PDB-6b1n.pdb || echo "  TDP43: descarga manual"
wget -q "https://files.rcsb.org/download/1hl5.pdb" -O proteins/SOD1/PDB-1hl5.pdb || echo "  SOD1: descarga manual"  
wget -q "https://files.rcsb.org/download/6g99.pdb" -O proteins/FUS/PDB-6g99.pdb || echo "  FUS: descarga manual"

echo ""
echo "============================================"
echo " INSTALACION COMPLETA"
echo "============================================"
echo "  Proyecto: /opt/masive-als"
echo "  Python:   $(python3 --version)"
echo "  Vina:     $(vina --version 2>&1 | head -1)"
echo "  OpenBabel: $(obabel --version 2>&1 | grep -i 'open babel' || echo 'OK')"
echo ""
echo "  Siguiente paso: copiar scripts y lanzar"
echo "  python3 /opt/masive-als/src/run_full_pipeline.py"
