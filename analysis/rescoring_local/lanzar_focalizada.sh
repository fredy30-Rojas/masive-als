#!/bin/bash
# Rescoring MM-GBSA de la LISTA FOCALIZADA (17 sep 2026).
#
# Que es: las 179 moleculas de la libreria cuya quimica se parece a la de los
# activos conocidos de TDP-43 (bencilisoquinolinas planas y cationicar), que son
# las que el sitio de union de ARN reconoce. Medido hoy: la similitud ordena
# mejor que la energia (AUC 0,92 frente a 0,66).
#
# Protocolo: el VALIDADO del 11 sep, el mismo del recalculo de los 240.
#   - receptor congelado del snapshot (no se vuelve a preparar con PDBFixer)
#   - doble precision en CUDA
#   - 3 repeticiones por compuesto (se reporta media y desviacion)
#
# Poses: las acopladas en el sitio de union de ARN (caja sacada del 4BS2), con
# --poses. No se usan las de la libreria porque esas son de la caja vieja, y son
# la prueba de la comparacion de cajas.
#
# Reanudable: si se corta, al volver a lanzarlo salta lo que ya tiene dG en
# rescoring_focalizada_arn.csv. La lista va ordenada de mejor quimica a peor, asi
# que lo primero que hace es lo mas valioso.
#
# Lanzado desde rescoring_focalizada.bat. Registro: runner_focalizada_arn.log
cd /mnt/c/Users/Fredy/masive-als/analysis/rescoring_local || exit 1
exec /root/rescoring_env/bin/python runner_mmgbsa_gb.py \
  --workers 4 \
  --lista lista_focalizada_rescoring.csv \
  --poses /mnt/c/Users/Fredy/masive-als/gpu_dock/rescoring_caja_arn \
  --tag focalizada_arn \
  --repeats 3 \
  --double \
  --receptores /mnt/c/Users/Fredy/masive-als/analysis/rescoring_local/receptores_fijos/snapshot_20260911
