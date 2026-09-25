#!/bin/bash
# Baja la prioridad (nice 19) del recalculo del rescoring para que Windows,
# OptiKey y el resto de programas de Fredy vayan primero.
# El hijo hereda el nice del padre, asi que basta con el runner; se aplica
# igualmente a los preparadores ya vivos por si acaso.
R=$(pgrep -f runner_mmgbsa_gb.py | head -1)
if [ -n "$R" ]; then
  renice -n 19 -p "$R" >/dev/null 2>&1
  echo "runner $R -> nice 19"
else
  echo "AVISO: no encontre el runner en marcha"
fi
for P in $(pgrep -f mmgbsa_openff_gb.py); do
  renice -n 19 -p "$P" >/dev/null 2>&1
done
echo "--- estado ---"
ps -eo pid,ni,pcpu,args --sort=-pcpu | grep -E 'runner_mmgbsa|mmgbsa_openff' | grep -v grep
