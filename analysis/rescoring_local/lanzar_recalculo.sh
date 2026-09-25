#!/bin/bash
# Recalculo del rescoring con el PROTOCOLO VALIDADO (13 sep 2026).
# Receptor congelado del snapshot, doble precision, 3 repeticiones.
# Lanzado desde lanzar_recalculo.bat. Reanudable: salta lo ya hecho.
cd /mnt/c/Users/Fredy/masive-als/analysis/rescoring_local || exit 1
exec /root/rescoring_env/bin/python runner_mmgbsa_gb.py \
  --workers 4 \
  --lista lista_recalculo.csv \
  --tag recalculo_validado \
  --repeats 3 \
  --double \
  --receptores /mnt/c/Users/Fredy/masive-als/analysis/rescoring_local/receptores_fijos/snapshot_20260911
