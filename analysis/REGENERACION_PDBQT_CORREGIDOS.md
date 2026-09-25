# Regeneración de PDBQT para los 916 ligandos corregidos

Fecha: 6 de septiembre de 2026 (23:26 local)

## Objetivo
Regenerar los PDBQT de los 916 ligandos con SMILES corregido (de la auditoría de
`ligandos_fallidos_clasificados.csv`, acción `CORREGIR_SMILES`) usando **el mismo
preparador de ligandos del pipeline local**.

## Preparador replicado (idéntico al pipeline)
El pipeline local convierte SMILES → PDBQT 3D con OpenBabel
(`convertir_tandas_v3.py` / `convertir_fda_full.py`):

1. Leer SMILES con `OBConversion.SetInFormat("smi")`
2. `mol.AddHydrogens()`
3. `OBBuilder.Build(mol)` → genera coordenadas 3D
4. `OBForceField UFF`: `SteepestDescent(150)` + `ConjugateGradients(50)` + `GetCoordinates`
5. Escribir PDBQT (`SetOutFormat("pdbqt")`, tipos AutoDock de OpenBabel)

Paralelismo: 10 workers (ProcessPoolExecutor), resumible (salta archivos ya generados).

## Nota sobre el entorno
El paquete pip de openbabel (3.1.0) no incluye los archivos de datos (`UFF.prm`,
`ring-fragments.txt`, `types.txt`, `atomtyp.txt`, `rigid-fragments-*`); se descargaron
de la rama `openbabel-3-1-0` del repositorio oficial al directorio de datos del paquete
para que el preparador funcione igual que en el pipeline.

## Resultados
| Métrica | Valor |
|---|---|
| SMILES de entrada | 916 |
| PDBQT generados (OpenBabel) | 901 |
| Fallos `no_genera_3d` (OBBuilder no puede con sistemas fusionados) | 15 |
| Rescatados con RDKit ETKDG + OpenBabel | 15 |
| **Total PDBQT válidos** | **916 / 916** |
| Tiempo total | ~1 min |

Validación de cada archivo generado:
- 916/916 contienen líneas ATOM/HETATM con coordenadas 3D no nulas.
- 916/916 usan solo tipos de átomo válidos de AutoDock (C, A, N, NA, OA, S, SA, HD, H, P, F, Cl, Br, I).
- Tamaño 886 – 13.099 bytes (ninguno vacío o truncado).

### Sobre los 15 rescatados con RDKit
Son ligandos con sistemas bicíclicos/policíclicos fusionados (p. ej. CHEMBL1186247,
CHEMBL1186304) que el builder 3D de OpenBabel no puede construir. Para ellos se generó
la conformación 3D con RDKit (ETKDGv3, torsiones de anillo pequeño/macrociclo) y se
convirtió a PDBQT con OpenBabel (que conserva las coordenadas y añade los tipos AD).
Para 13 de ellos fue necesario `enforceChirality=False` en el embedding: la
estereoquímica absoluta se representa pero puede no respetar el isómero original en
todos los centros; se recomienda revisar esos 15 con un vistazo visual antes del
rescoring final si alguno aparece entre los top candidatos.

## Salida
```
analysis/_redock_corregidos/pdbqt/<ligand>.pdbqt   (916 archivos)
```
Entrada: `analysis/ligandos_a_reacoplar.smi` (916 SMILES corregidos)

## Siguiente paso
Re-acoplar estos 916 contra TDP-43 v2 (centro 24.23, 16.89, -15.87; tamaño 26;
semillas 42/2026/777) con Vina-GPU local o Kaggle. La caja es la misma del cribado
v4 para que los resultados sean comparables con los 8.001 ya validados.

## Scripts
- `analysis/_regenerar_pdbqt_corregidos.py` — conversor principal (OpenBabel, 10 workers)
- `analysis/_rescatar_fallos_3d_rdkit.py` — rescate de los 15 con RDKit ETKDG
- Log: `analysis/_redock_corregidos/log.txt`