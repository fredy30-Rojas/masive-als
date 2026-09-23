@echo off
REM Espera a que TERMINE el acople del fondo duro R-BIND y entonces lanza la
REM validacion limpia de TDP-43 (fondo blando + fondo duro), que es el paso 3
REM que quedo pendiente en INFORME_RESCORING_ESTRATOS_2026-09-22.md.
REM Uso: wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\acople_y_validar_tdp43.bat"
cd /d C:\Users\Fredy\masive-als\analysis
python -X utf8 esperar_acople_y_validar.py > _esperar_acople_out.log 2>&1
