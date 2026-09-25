@echo off
REM Segundo proceso del fondo duro: recorre la lista por el final, para repartir
REM el acople entre dos procesos (la maquina tiene 20 nucleos y cada vina va con
REM --cpu 1). Es reanudable, igual que el otro.
REM Uso: wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\acoplar_fondo_rbind_rev.bat"
cd /d C:\Users\Fredy\masive-als\analysis
python -X utf8 acoplar_fondo_rbind.py --reversa --log acoplar_rev.log --hilos 8 > _acoplar_fondo_rbind_rev_out.log 2>&1
