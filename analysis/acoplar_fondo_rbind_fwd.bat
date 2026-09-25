@echo off
REM Tercer proceso del fondo duro: recorre la lista por el PRINCIPIO, para
REM encontrarse con el que va por el final. Los dos son reanudables y se saltan
REM lo que ya tiene pose, asi que no repiten trabajo.
REM Uso: wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\acoplar_fondo_rbind_fwd.bat"
cd /d C:\Users\Fredy\masive-als\analysis
python -X utf8 acoplar_fondo_rbind.py --log acoplar_fwd.log --hilos 6 > _acoplar_fondo_rbind_fwd_out.log 2>&1
