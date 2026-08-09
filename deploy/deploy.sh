#!/bin/bash
# ============================================================
# Deploy script — subir webs al servidor Oracle
# IP: 79.72.57.253
# ============================================================

set -e

SERVER="79.72.57.253"
SSH_KEY="$HOME/.ssh/id_rsa"
LOCAL_DIR="$(dirname "$0")/.."

echo "=== MASIVE-ALS + Conmimirada Deploy ==="
echo ""

# Crear carpetas en el servidor
echo "[1/5] Creando estructura en el servidor..."
ssh -i "$SSH_KEY" "fredy@$SERVER" "
    sudo mkdir -p /var/www/portal
    sudo mkdir -p /var/www/masive-als
    sudo mkdir -p /var/www/accesibilidad
    sudo chown -R fredy:fredy /var/www
"

# Subir portal principal
echo "[2/5] Subiendo portal principal..."
scp -i "$SSH_KEY" "$LOCAL_DIR/web/portal.html" "fredy@$SERVER:/var/www/portal/index.html"

# Subir MASIVE-ALS
echo "[3/5] Subiendo MASIVE-ALS..."
scp -i "$SSH_KEY" "$LOCAL_DIR/web/index.html" "fredy@$SERVER:/var/www/masive-als/index.html"
scp -i "$SSH_KEY" "$LOCAL_DIR/web/presentacion.html" "fredy@$SERVER:/var/www/masive-als/presentacion.html"

# Subir accesibilidad (si hay cambios locales)
echo "[4/5] Verificando página de accesibilidad..."
# La página de accesibilidad ya existe en el servidor, la dejamos como está

# Configurar nginx
echo "[5/5] Configurando nginx..."
scp -i "$SSH_KEY" "$LOCAL_DIR/deploy/nginx-sites.conf" "fredy@$SERVER:/tmp/nginx-sites.conf"
ssh -i "$SSH_KEY" "fredy@$SERVER" "
    sudo cp /tmp/nginx-sites.conf /etc/nginx/sites-available/masive-als
    sudo ln -sf /etc/nginx/sites-available/masive-als /etc/nginx/sites-enabled/masive-als
    sudo nginx -t && sudo systemctl reload nginx
    echo 'Nginx recargado OK'
"

echo ""
echo "=== Deploy completado ==="
echo "Portal:   http://$SERVER/"
echo "MASIVE:   http://$SERVER/masive-als/"
echo "Conmimirada: http://$SERVER/accesibilidad/"
