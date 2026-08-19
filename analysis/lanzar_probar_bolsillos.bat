@echo off
start "" /min "C:\Users\Fredy\AppData\Local\Python\pythoncore-3.14-64\python.exe" "C:\Users\Fredy\masive-als\analysis\probar_bolsillos_sod1.py" --ligandos "C:\Users\Fredy\masive-als\analysis\_validacion_SOD1\ligands" --receptor "C:\Users\Fredy\masive-als\gpu_dock\SOD1.pdbqt" --exhaustividad 8 --cpu-workers 6 --solo trp32,metal,dimer
