@echo off
cd /d C:\Users\Fredy\masive-als\analysis
echo === VALIDACION TDP43 RRM (receptor 4BS2 modelo 1, cadena A 96-269) === %date% %time%
echo Activos RRM: rTRD01, 5FUrd (sitios de union a RNA validados por literatura)
echo [1/2] Sub-sitio RRM1: centro 26.4,20.2,-2.7 tamano 22
python validar_senuelos.py --activos activos_tdp43_rrm.csv --libreria full_library_solo.smi --target TDP43 --receptor ../gpu_dock/proteins_TDP43/4BS2_chainA.pdb --centro 26.4,20.2,-2.7 --tamano 22 --cpu --decoys-por-activo 30 --salida validacion_TDP43_RRM1.csv
echo [2/2] Sub-sitio RRM2: centro 20.3,33.0,-14.3 tamano 22
python validar_senuelos.py --activos activos_tdp43_rrm.csv --libreria full_library_solo.smi --target TDP43 --receptor ../gpu_dock/proteins_TDP43/4BS2_chainA.pdb --centro 20.3,33.0,-14.3 --tamano 22 --cpu --decoys-por-activo 30 --salida validacion_TDP43_RRM2.csv
echo === FIN VALIDACION RRM === %date% %time%
