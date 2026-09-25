@echo off
REM Validacion limpia de TDP-43 contra el receptor corregido 4BS2_ph74.
REM Desde el 23 sep 2026 incluye el fondo duro (R-BIND 2.0): tres bloques.
REM Uso: wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\lanzar_validar_tdp43_limpia.bat"
cd /d C:\Users\Fredy\masive-als\analysis
python -X utf8 validar_diana_limpia.py --target TDP43 > _validar_tdp43_limpia_out.log 2>&1
