# -*- coding: utf-8 -*-
"""Mantiene el PC despierto mientras corre el recalculo del rescoring.

Por que existe: la noche del 5 al 6 de septiembre el equipo entro en Modern
Standby (apagado de pantalla sin proteccion) y congelo la GPU a mitad de un
chunk del cribado. Aqui se repite el mismo patron: un calculo largo en la
tarjeta mientras Fredy duerme.

Uso:  python mantener_despierto.py [minutos]     (por defecto 130)
La pantalla SI puede apagarse: solo se pide que el SISTEMA no se suspenda.
"""
import ctypes
import sys
import time

ES_CONTINUOUS = 0x80000000
ES_SYSTEM_REQUIRED = 0x00000001

minutos = float(sys.argv[1]) if len(sys.argv) > 1 else 130.0
k = ctypes.windll.kernel32
k.SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED)
print("manteniendo el sistema despierto %.0f minutos (inicio %s)"
      % (minutos, time.strftime("%H:%M:%S")), flush=True)
time.sleep(minutos * 60)
k.SetThreadExecutionState(ES_CONTINUOUS)
print("fin del bloqueo de suspension", flush=True)
