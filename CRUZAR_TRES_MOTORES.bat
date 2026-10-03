@echo off
REM Cruza Vina, GNINA y MM-GBSA en resultados\consenso_tres_motores.csv
setlocal
cd /d "%~dp0"
python consenso_tres_motores.py --gnina-top 5000
echo.
echo Escrito en resultados\consenso_tres_motores.csv
pause
endlocal
