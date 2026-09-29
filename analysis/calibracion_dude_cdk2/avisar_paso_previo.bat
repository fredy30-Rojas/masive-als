@echo off
REM Avisa a Fredy por voz cuando el paso previo del encadenador de CDK2 termine.
REM
REM Por que existe: Fredy pidio el 29 de septiembre que le avisaran en cuanto
REM acabara el paso previo (la busqueda agotada con depth 128), para leer el
REM resultado. Eso puede ser dentro de unas horas, a las once de la noche, y no
REM puede quedar pendiente de que alguien este mirando.
REM
REM Mira `encadenar_cdk2.log` y espera a la linea "remedicion hecha", que es lo
REM ULTIMO que hace el paso previo antes de que arranque el banco de cuatro horas.
REM Si el log se queda quieto hora y media sin dar la senal, tambien avisa: si no,
REM Fredy se queda esperando toda la noche un aviso que no va a llegar.
REM
REM Uso:  wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\calibracion_dude_cdk2\avisar_paso_previo.bat"
setlocal
cd /d "%~dp0"

REM OJO: el .py ya escribe su propio registro en `aviso_paso_previo.log`. Aqui la
REM consola va a OTRO fichero a proposito: si los dos apuntan al mismo, cmd lo tiene
REM abierto y python no puede abrirlo (PermissionError).
python "%~dp0avisar_paso_previo.py" >> "%~dp0avisar_paso_previo_consola.log" 2>&1
exit /b %ERRORLEVEL%
