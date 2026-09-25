@echo off
REM Acopla el fondo duro de unidores de ARN (R-BIND 2.0) contra TDP-43.
REM Reanudable: lo que ya tenga pose no se repite.
REM Uso: wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\acoplar_fondo_rbind.bat"
cd /d C:\Users\Fredy\masive-als\analysis
REM %* reenvia los argumentos: con dos procesos, el segundo va con --reversa
REM --log acoplar_rev.log (sin el %* los dos hacian lo mismo y se pisaban el log).
python -X utf8 acoplar_fondo_rbind.py %* > _acoplar_fondo_rbind_out.log 2>&1
