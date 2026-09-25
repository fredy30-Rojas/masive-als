#!/bin/bash
# Para el recalculo del rescoring antes de relanzarlo con otra configuracion.
# El trabajo es reanudable: lo ya escrito en el CSV no se pierde.
pkill -f runner_mmgbsa_gb.py
pkill -f mmgbsa_openff_gb.py
sleep 3
echo "procesos que quedan:"
pgrep -af 'runner_mmgbsa|mmgbsa_openff' | grep -v pgrep
echo "fin de la parada"
