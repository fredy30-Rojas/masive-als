@echo off
REM Barrido de exhaustividad en CPU: el mismo conjunto de TDP-43 a exhaustividad
REM 8 y 32, con el mismo Vina, el mismo receptor, la misma caja y los mismos
REM ficheros de ligando. El punto medio (exhaustividad 16) ya existe: es la
REM validacion limpia. Reanudable: lo que ya tiene pose no se repite.
REM Uso: wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\barrido_exhaustividad_tdp43.bat"
cd /d C:\Users\Fredy\masive-als\analysis
python -X utf8 barrido_exhaustividad_tdp43.py --todo > _barrido_exhaustividad_out.log 2>&1
