@echo off
REM Espera al banco de CDK2 y luego hace el control de redocking de ANDR.
REM
REM Por que existe: `encadenar_cdk2.py` ya estaba corriendo cuando se anadio el
REM control de andr al final de la cadena, asi que no le llega el cambio (un proceso
REM no recarga su codigo). Este script hace SOLO esa parte nueva y espera por su
REM cuenta a que la tarjeta quede libre, para no pelearse con el banco de CDK2 ni
REM con el de TBK1.
REM
REM Uso:  wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\calibracion_dude_cdk2\encadenar_andr.bat"
setlocal
cd /d "%~dp0"

REM OJO: el .py ya escribe su propio registro en `encadenar_andr.log`. Aqui la
REM consola va a OTRO fichero a proposito: si los dos apuntan al mismo, cmd lo tiene
REM abierto y python no puede abrirlo (PermissionError).
if "%~1"=="--ya" (
    python "%~dp0encadenar_andr.py" --ya >> "%~dp0encadenar_andr_consola.log" 2>&1
) else (
    python "%~dp0encadenar_andr.py" >> "%~dp0encadenar_andr_consola.log" 2>&1
)
exit /b %ERRORLEVEL%
