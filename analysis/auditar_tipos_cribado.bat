@echo off
REM Lupa de tipos de atomo sobre los ligandos ya preparados del cribado:
REM busca descartes silenciosos (tipo que Vina no acepta) y elementos criticos
REM nombrados pero tipados de otra cosa. No modifica nada.
REM Uso: wscript.exe ejecutar_bat_oculto.vbs "C:\Users\Fredy\masive-als\analysis\auditar_tipos_cribado.bat"
cd /d C:\Users\Fredy\masive-als\analysis
python -X utf8 auditar_tipos_cribado.py > _auditoria_tipos_out.log 2>&1
