@echo off
REM Espera al banco de CDK2 y luego pasa el validador para leer el AUC de punta a punta.
REM
REM Por que existe: Fredy pidio el 29 de septiembre que, en cuanto CDK2 terminara, se
REM pasara el validador y se le leyera el AUC completo. La secuencia de la noche es
REM larga (TBK1, luego el paso previo, luego cuatro horas de banco) y no puede quedar
REM pendiente de que alguien este mirando el reloj.
REM
REM Espera a que `resultados_cdk2.csv` exista y lleve 10 min quieto, que es la senal de
REM que el lanzador ha terminado. Si el registro del lanzador se queda quieto una hora,
REM da el banco por parado y mide igual: un banco muerto no escribe el CSV nunca, y
REM quedarse esperando seria perder la noche.
REM
REM Al terminar, habla el AUC crudo, el AUC por atomo pesado y el EF1, leidos del
REM informe que escribe el validador. Si las cifras no se pueden leer, lo dice en voz
REM alta en vez de inventar un numero.
REM
REM Uso:  wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\calibracion_dude_cdk2\encadenar_validacion.bat"
setlocal
cd /d "%~dp0"

REM OJO: el .py ya escribe su propio registro en `encadenar_validacion.log`. Aqui la
REM consola va a OTRO fichero a proposito: si los dos apuntan al mismo, cmd lo tiene
REM abierto y python no puede abrirlo (PermissionError).
if "%~1"=="--ya" (
    python "%~dp0encadenar_validacion.py" --ya >> "%~dp0encadenar_validacion_consola.log" 2>&1
) else (
    python "%~dp0encadenar_validacion.py" >> "%~dp0encadenar_validacion_consola.log" 2>&1
)
exit /b %ERRORLEVEL%
