@echo off
cd /d C:\Users\Fredy\masive-als\analysis
python validar_senuelos.py --activos activos_tdp43.csv --libreria ..\compounds\decoys_library.smi --target TDP43 --receptor ..\gpu_dock\TDP43.pdbqt --centro 16.3,41.1,48.5 --tamano 24 --decoys-por-activo 30 --cpu --salida validacion_TDP43_v4.csv > validacion_tdp43_v4_stdout.log 2> validacion_tdp43_v4_err.log
