@echo off
rem ============================================================
rem  Rescoring MM-GBSA de la LISTA FOCALIZADA (17 sep 2026)
rem ============================================================
rem  Que es: las 179 moleculas de la libreria cuya quimica se
rem  parece a la de los activos conocidos de TDP-43. La
rem  similitud ordena mejor que la energia (AUC 0,92 vs 0,66).
rem
rem  Protocolo VALIDADO del 11 sep, el mismo del recalculo de
rem  los 240: receptor congelado del snapshot (no se prepara
rem  otra vez con PDBFixer), doble precision en CUDA, y 3
rem  repeticiones por compuesto.
rem
rem  Poses: las acopladas en el sitio de union de ARN (caja
rem  sacada del 4BS2), copiadas a rescoring_caja_arn.
rem
rem  Reanudable: si se corta o hay que pararlo, al volver a
rem  lanzarlo salta lo que ya tenga dG. Va ordenado de mejor
rem  quimica a peor, asi que hace primero lo mas valioso.
rem
rem  Salida: rescoring_focalizada_arn.csv
rem  Registro: runner_focalizada_arn.log
rem  ============================================================
cd /d "C:\Users\Fredy\masive-als\analysis\rescoring_local"
wsl.exe -d Ubuntu -u root bash /mnt/c/Users/Fredy/masive-als/analysis/rescoring_local/lanzar_focalizada.sh >> rescoring_focalizada_stdout.log 2>&1
