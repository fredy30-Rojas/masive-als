@echo off
rem Recalculo del rescoring MM-GBSA con el protocolo validado (13 sep 2026).
rem Entra en WSL como root, con el entorno del rescoring, y deja el registro
rem en recalculo_stdout.log. Reanudable si se corta.
cd /d "C:\Users\Fredy\masive-als\analysis\rescoring_local"
wsl.exe -d Ubuntu -u root bash /mnt/c/Users/Fredy/masive-als/analysis/rescoring_local/lanzar_recalculo.sh >> "C:\Users\Fredy\masive-als\analysis\rescoring_local\recalculo_stdout.log" 2>&1
