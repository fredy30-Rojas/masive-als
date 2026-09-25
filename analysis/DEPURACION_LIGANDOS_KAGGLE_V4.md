# Depuración de ligandos que fallan en conversión PDBQT — Kaggle v4 (TDP-43)

Fecha: 6 de septiembre de 2026
Fuente de verdad: `kaggle_outputs/masive-als-cribado-master-tdp43-v4/resultados_master_tdp43v4.csv`
Librería: `analysis/_lib_consolidada.csv` (56.772) + `analysis/full_library_solo.smi` (66.478)

---

## 1. Corrección de la estimación previa

En un informe anterior se mencionó "~5.700 ligandos que siempre fallan". **Esa cifra era incorrecta**:
el contador `skipped=5692` del log PROGRESS de Kaggle v4 **no son fallos**, son ligandos que ya
figuraban en `ligandos_completados.txt` (checkpoint restaurado de la sesión anterior) y se saltaron
por reanudación. No son conversiones fallidas.

La cifra real de fallos de Kaggle v4 es **1.080 ligandos únicos**, extraída de la columna `status`
del CSV (todo ligando con `error:...` en sus 3 semillas). De 9.081 filas parseadas: **8.001 ok** y
**1.080 fallos** (tasa de fallo ~11,9 % de lo procesado en esa sesión).

## 2. Categorías de fallo (Kaggle v4, 1.080 ligandos)

| Categoría | Cantidad | Interpretación |
|---|---|---|
| `etiqueta_desconocida` | 992 | PDBQT con tags que Vina no reconoce (casi todas sales: `_SULFATE`, `_MALEATE`, `_SODIUM`, `_HYDROCHLORIDE`, `_BROMIDE`...) |
| `afinidad_no_encontrada` | 32 | Vina corrió pero no devolvió afinidad parseable |
| `multi_modelo_pose` | 30 | Archivos de **pose de acoplamientos previos** (`*_out.pdbqt`, `*_SOD1_out.pdbqt`, `*_TDP43_v2_out.pdbqt`) colados como ligandos |
| `metal_tipo_AD:*` | 24 | Átomos no tipables por AutoDock (B, As, Al, Li, Si, Rb, Sr, Ag) |
| otro | 2 | Misceláneo |

## 3. Clasificación final con acción (1.080)

Archivo: `analysis/ligandos_fallidos_clasificados.csv`

| Acción | Cantidad | Criterio |
|---|---|---|
| `CORREGIR_SMILES` | **916** | Sal/contraion: se resuelve el SMILES del fármaco base (p. ej. `ABACAVIR_SULFATE` → SMILES de `ABACAVIR`) o se quita el contraión con RDKit SaltRemover. **Hay que regenerar el PDBQT desde ese SMILES limpio y re-acoplar.** |
| `EXCLUIR_METAL` | **51** | Organometálicos / sales metálicas: Na, K, Ca, Mg, Al, Zn, Li, Se, As, Te, Sr, Ag, Rb, Ba... AutoDock no puede tipar esos átomos; **no son acoplables con Vina** (ni siquiera regenerando). Se excluyen del cribado. |
| `EXCLUIR_POSE` | **30** | Archivos multi-MODEL de poses de acoplamientos previos; **no son ligandos**; se eliminan de la librería. |
| `REVISAR` | **83** | Sin SMILES resoluble (86 sin SMILES; de ellos 83 caen aquí) o afinidad no encontrada; requiere revisión manual del origen del archivo. |

Nota: un mismo ligando puede aparecer en varias categorías de `categoria`; la acción es única
(primera regla que aplica: metal → pose → sal → revisar).

## 4. Barrido RDKit de TODA la librería (adelanto de fallos futuros)

Archivo: `analysis/libreria_problematicos_rdkit.csv`

Sobre **119.090 entradas SMILES** consolidadas (56.772 de `_lib_consolidada` + 66.478 de
`full_library_solo.smi`), RDKit detectó **1.340 SMILES que contienen metales/metaloides**:

| Metal | Entradas | Metal | Entradas |
|---|---|---|---|
| Na | 868 | Li | 20 |
| K | 124 | Te | 16 |
| Se | 83 | As | 13 |
| Ca | 81 | Sr | 9 |
| Mg | 50 | Ag | 8 |
| Al | 22 | Rb | 6 |
| Zn | 22 | Ra | 6 |
| | | Ba | 3 |

Interpretación científica:
- **Na/K/Ca/Mg/Li/Zn (1.185)** son casi todas **sales** (formas farmacéuticas: `_SODIUM`,
  `_POTASSIUM`, `_CALCIUM`, `_MAGNESIUM`, `_ZINC`): **CORREGIBLES** quitando el contraión y
  regenerando el PDBQT. Son las que producen `etiqueta_desconocida` en Vina.
- **Se/As/Te/Al/Sr/Ag/Rb/Ba (160)** son organometálicos/metaloides **reales en el esqueleto**:
  **NO tipables por AutoDock** → excluir del cribado con Vina (no tienen tipo AD; Vina-GPU
  rechaza el archivo). Si un candidato de este grupo fuera farmacológicamente relevante, el
  acoplamiento requeriría otro motor (p. ej. GNINA con tipos de átomo extendidos).

Esto implica que, de la librería completa, **~1.185 ligandos son corregibles por SMILES limpio**
y **~160 son excluibles de forma permanente** por metal en esqueleto. El resto de fallos
(poses y sin SMILES) se resuelve por higiene de la librería.

## 5. Acción recomendada (orden)

1. **Regenerar PDBQT limpio** para los 916 `CORREGIR_SMILES` desde `smiles_corregido`
   (obtener 3D con el mismo preparador de ligandos usado localmente, p. ej. RDKit ETKDG +
   conversión a PDBQT con el script de la tanda), y **re-acoplarlos** contra TDP-43 v2.
2. **Eliminar los 30 `EXCLUIR_POSE`** de los archivos de ligandos (higiene: no subir poses
   `*_out.pdbqt` como ligandos).
3. **Excluir los 51 `EXCLUIR_METAL`** del cribado Vina y documentarlos (lista para posible
   cribado con GNINA si interesa).
4. **Revisar los 83 `REVISAR`**: verificar de qué tanda/archivo salieron los PDBQT sin SMILES.
5. Actualizar el `ligandos_completados.txt` del checkpoint con los excluidos para que el
   relanzador no vuelva a intentarlos.

## 6. Archivos generados

| Archivo | Contenido |
|---|---|
| `analysis/ligandos_fallidos_kaggle_v4.csv` | Los 1.080 fallos con categoría, SMILES y motivo |
| `analysis/ligandos_fallidos_clasificados.csv` | Mismos 1.080 con **acción** (`CORREGIR_SMILES` / `EXCLUIR_METAL` / `EXCLUIR_POSE` / `REVISAR`), SMILES corregido y razón |
| `analysis/libreria_problematicos_rdkit.csv` | 1.340 SMILES con metales de toda la librería (adelanto) |

Scripts: `analysis/_auditar_fallos_kaggle_v4.py`, `analysis/_clasificar_fallos_final.py`