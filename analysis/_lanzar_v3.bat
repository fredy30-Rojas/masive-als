@echo off
cd /d C:\Users\Fredy\masive-als\analysis
python validar_senuelos.py --activos activos_tdp43.csv --libreria ..\compounds\decoys_library.smi --target TDP43 --receptor ..\gpu_dock\TDP43.pdbqt --centro 24.4,44.9,53.4 --tamano 18 --decoys-por-activo 30 --cpu --salida validacion_TDP43_v3.csv > validacion_tdp43_v3_stdout.log 2> validacion_tdp43_v3_err.log
