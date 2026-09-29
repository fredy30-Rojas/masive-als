@echo off
REM Encadena la validacion de CDK2 detras del banco de TBK1 (29 sep 2026).
REM
REM Por que existe: las dos dianas se acoplan en la misma RTX 4080 y Vina-GPU no
REM comparte tarjeta (TBK1 usa casi los 12 GB). El que va segundo tiene que
REM esperar a que el primero suelte, y eso no se puede quedar pendiente de que
REM alguien mire el reloj a las once de la noche.
REM
REM Lanza `encadenar_cdk2.py`, que espera a que el banco de TBK1 termine (o se
REM pare) y entonces acopla la submuestra 1:10 de CDK2: 474 activos contra 4.740
REM senuelos. Todo lo que hace es reanudable.
REM
REM Uso:  wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\calibracion_dude_cdk2\encadenar_cdk2.bat"
setlocal
cd /d "%~dp0"

REM OJO: el .py ya escribe su propio registro en `encadenar_cdk2.log`. Aqui se
REM manda la consola a OTRO fichero a proposito: si los dos apuntan al mismo,
REM cmd lo tiene abierto y python no puede abrirlo (PermissionError).
if "%~1"=="--ya" (
    python "%~dp0encadenar_cdk2.py" --ya %2 >> "%~dp0encadenar_cdk2_consola.log" 2>&1
) else (
    python "%~dp0encadenar_cdk2.py" >> "%~dp0encadenar_cdk2_consola.log" 2>&1
)
exit /b %ERRORLEVEL%
