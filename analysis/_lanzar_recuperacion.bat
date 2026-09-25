@echo off
REM Espera a que el re-docking masivo principal termine y luego recupera los timeouts.
cd /d C:\Users\Fredy\masive-als\analysis
python -c "import time, os, csv; \
import subprocess, sys; \
out='_redock_tdp43_masivo/checkpoint.csv'; \
t0=time.time(); \
last=0; \
while True: \
    n=sum(1 for _ in open(out)) if os.path.exists(out) else 0; \
    if n>0 and n==last and time.time()-t0>120: \
        break; \
    last=n; \
    time.sleep(60); \
print('Pool principal terminado en %.0fs (%d lineas). Lanzando recuperacion...' % (time.time()-t0, last))"
python _recuperar_timeouts_tdp43.py
