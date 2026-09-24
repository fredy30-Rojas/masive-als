@echo off
REM La regla de decision con su error medido: AUC por atomo pesado contra el fondo
REM duro, con intervalo bootstrap, aplicada a las cuatro corridas de TDP-43 y a
REM SOD1. Solo lee CSV ya puntuados: no acopla nada y tarda menos de un minuto.
REM Uso: wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\regla_decision.bat"
cd /d C:\Users\Fredy\masive-als\analysis
python -X utf8 regla_decision_bootstrap.py > _regla_decision_out.log 2>&1
