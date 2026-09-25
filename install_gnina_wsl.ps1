Script de Instalación de WSL2 + Ubuntu + GNINA para Fredy
========================================================

# PASO 1: Instalar WSL2 (como administrador)
# Ejecutar en PowerShell como administrador:
wsl --install -d Ubuntu-22.04
Restart-Computer

# PASO 2: Configurar Ubuntu (después de reiniciar)
# Abrir Ubuntu y crear usuario:
sudo apt update && sudo apt upgrade -y

# PASO 3: Verificar CUDA en WSL
nvidia-smi

# PASO 4: Instalar GNINA (Opción A - binario listo)
cd ~ && mkdir gnina && cd gnina

# Descargar el release más reciente de GNINA
wget https://github.com/gnina/gnina/releases/download/v1.1/gnina-1.1-linux.zip
unzip gnina-1.1-linux.zip
chmod +x gnina

# Verificar instalación
./gnina --version

# PASO 5: Prueba real con datos de Fredy
mkdir -p /mnt/c/Users/Fredy/masive-als/gpu_dock/gnina_test
cd /mnt/c/Users/Fredy/masive-als/gpu_dock/gnina_test

# Tomar 100 poses reales ya generadas por Vina-GPU
ls /mnt/c/Users/Fredy/masive-als/gpu_dock/resultados_libreria/results_TDP43/*_out.pdbqt | head -100 > lista.txt
cat lista.txt | xargs cat > poses_prueba.pdbqt

# Re-puntuar con la CNN
~/gnina/gnina -r /mnt/c/Users/Fredy/masive-als/gpu_dock/TDP43.pdbqt -l poses_prueba.pdbqt --cnn_scoring rescore --out rescored.sdf --num_modes 1

# PASO 6: Guardar resultados
echo "=== GNINA Installation Test Results ===" > ~/gnina/log_preparacion.txt
echo "Date: $(date)" >> ~/gnina/log_preparacion.txt
echo "GNINA Version: $(./gnina --version)" >> ~/gnina/log_preparacion.txt
echo "nvidia-smi output:" >> ~/gnina/log_preparacion.txt
nvidia-smi >> ~/gnina/log_preparacion.txt
echo "Number of poses re-scored: 100" >> ~/gnina/log_preparacion.txt
echo "Top 5 CNN affinity scores:" >> ~/gnina/log_preparacion.txt
grep "minimizedAffinity\|cnn_affinity" rescored.sdf | head -5 >> ~/gnina/log_preparacion.txt

echo "Installation completed successfully!"