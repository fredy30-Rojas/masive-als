@echo off
REM CONTROL con Vina-GPU: el conjunto completo de TDP-43 (positivos + señuelos +
REM fondo duro de R-BIND) acoplado con la GPU, para ver si el motor cambia el orden.
REM Uso: wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\control_gpu_tdp43.bat"
cd /d C:\Users\Fredy\masive-als\analysis
python -X utf8 control_gpu_tdp43.py --todo > _control_gpu_tdp43_out.log 2>&1
