# Los primeros puestos de TODAS las validaciones eran pseudo-átomos «glue»

**Fecha:** 21 de septiembre de 2026
**Ámbito:** conjuntos de validación de TDP-43, FUS y SOD1
**Herramientas:** `reparar_ligandos_glue.py`, `reparar_validaciones_glue.py`,
`reparar_sod1_v3_glue.py`, `confirmar_glue_controles.py`
**Evidencia:** `reparar_validaciones_glue.log`, `reparar_sod1_v3_glue.log`,
`confirmar_glue_controles.log`
**Corrige:** el diagnóstico de `REVISION_VALIDACION_BOLSILLOS_2026-08-20.md`
(«artefacto de tamaño de Vina») y el aviso de `INFORME_VALIDACION_SOD1_V4_2026-09-20.md`

---

## 1. El defecto

Meeko, con los ajustes por defecto, **no cierra los anillos de 7 eslabones en
adelante**. Los abre en dos ramas y pega los extremos con pseudo-átomos de
pegamento (tipos `CG0` y `G0`, elementos `Du`). Cada anillo abierto añade dos
átomos fantasma al fichero.

Vina reporta `score = intermolecular + intramolecular − unbound`. En el ligando
aislado (columna `REMARK UNBOUND` de la pose) el anillo abierto tiene una
energía interna de **+7 a +16 kcal/mol** en lugar de ≈ 0. Ese exceso se le
resta al score, así que el compuesto recibe un regalo de varios kcal/mol que no
existe. Además el ligando tiene 2 átomos de más, y eso rompe cualquier
normalización por tamaño.

El proyecto **ya conocía este fallo para la librería** (`gpu_dock/reparar_libreria_glue.py`,
3.803 afectados, re-acoplados solo contra TDP43_v2). Los ficheros de control de
las validaciones se prepararon con el Meeko viejo y nunca se repararon.

## 2. Alcance

| Conjunto | Afectados | Papel en el ranking | AUC  antes → después |
|---|---|---|---|
| TDP43 v6 (4BS2) | 16 | los 14 primeros señuelos + **cefarantina** (la mejor activa) | 0,517 → **0,534** crudo · 0,534 → 0,554 por átomo · 0,556 → 0,539 por raíz |
| FUS v2 | 3 | puestos 1-3 (3 señuelos) | 0,609 → **0,626** |
| SOD1 trp32 (v1) | 2 | puestos 1-2 (2 señuelos) | 0,797 → **0,803** |
| SOD1 v3 | 11 | todo el fondo (5 duros + 6 emparejados) | 0,671 → **0,686** crudo · 0,486 → 0,532 sin tamaño (familia «independientes») |

Debajo de los reparados, el patrón es el mismo en los cuatro: **cada uno pierde
entre 3,2 y 6,1 kcal/mol** y su `UNBOUND` cae de +7,4…+15,8 a ≈ −0,4.

### TDP43 v6 (16)

| Ligando | Rol | Antes | Después | Cambio | UNBOUND |
|---|---|---|---|---|---|
| ACT_cepharanthine | **activo** | −10,72 | −7,47 | +3,25 | 8,16 |
| DEC_CHEMBL607833 | señuelo | −13,63 | −9,39 | +4,23 | 15,76 |
| DEC_CHEMBL452 | señuelo | −12,75 | −7,02 | +5,73 | 9,59 |
| DEC_CHEMBL13280 | señuelo | −12,71 | −6,85 | +5,86 | 9,63 |
| DEC_CHEMBL1173055 | señuelo | −12,47 | −8,59 | +3,89 | 8,70 |
| DEC_CHEMBL580 | señuelo | −12,29 | −6,22 | +6,07 | 8,99 |
| DEC_CHEMBL1213252 | señuelo | −12,11 | −6,88 | +5,23 | 10,13 |
| DEC_CHEMBL277062 | señuelo | −12,02 | −7,25 | +4,78 | 9,34 |
| DEC_CHEMBL1435671 | señuelo | −11,53 | −7,85 | +3,68 | 8,59 |
| DEC_CHEMBL22097 | señuelo | −11,19 | −6,66 | +4,53 | 9,49 |
| DEC_CHEMBL1697737 | señuelo | −11,15 | −6,20 | +4,95 | 9,18 |
| DEC_CHEMBL2356527 | señuelo | −10,94 | −7,54 | +3,40 | 8,91 |
| DEC_CHEMBL1310614 | señuelo | −10,17 | −6,24 | +3,94 | 8,23 |
| DEC_CHEMBL538188 | señuelo | −9,46 | −6,28 | +3,18 | 7,38 |
| DEC_CHEMBL10183 | señuelo | −6,63 | −8,32 | −1,69 | −1,45 |
| DEC_CHEMBL8960 | señuelo | −13,49 | −8,03 | +5,46 | 9,47 |

Los 14 señuelos que ocupaban los puestos 1 a 14 eran **moléculas normales de 22 a
54 átomos pesados** (el mayor, CHEMBL607833, tiene 54). No eran «compuestos
grandes»: el puesto 15, el primer señuelo sin glue, tenía −9,43 kcal/mol, cinco
grandes por debajo del tope.

Antes de reparar, los 14 señuelos glue y cefaranthine ocupaban los puestos 1 a 15
de 131. Después de reparar, cefaranthine (−7,47) pasa a ser lo mejor de la lista
global y el primer señuelo es CHEMBL607833 con −9,39.

## 3. Qué cambia y qué no

**No cambia:** las conclusiones de fondo.

- TDP-43 v6 sigue en ≈ azar (AUC 0,534 crudo). La promesa del receptor 4BS2 se
  mantiene como estaba: los activos acoplan bien, pero el ranking crudo no
  discrimina; hace falta rescoring.
- SOD1 trp32 (v1) sigue con AUC 0,803 **crudo**, pero era y sigue siendo
  mayormente el sesgo de tamaño a favor de activos grandes: la normalización por
  átomos lo baja a ~0,5.
- SOD1 v3 sigue sin pasar: 0,686 crudo / 0,532 sin tamaño en la familia
  «independientes»; las familias «serie» y «co-cristalizados» no son
  independientes y no cuentan como validación.

**Sí cambia:**

1. **La causa del ranking de TDP-43 no era el tamaño.** Era un fichero de
   ligando mal preparado. El sesgo de tamaño existe en el fondo
   (≈ −0,09 kcal/mol por carbono pesado en el bolsillo 4BS2), pero no explicaba
   la cabeza.
2. **El score de Vina tenía señal en TDP-43, tapada por el glue.** Tras reparar,
   la activa cefarantina (−7,47) queda por delante del mejor señuelo (−9,39).
   Antes, 14 señuelos inflados artificialmente la superaban.
3. **El defecto estaba en todas las validaciones, no solo en TDP-43.** Cualquier
   tabla de esta carpeta que cite puestos 1-3 debe revisarse.

## 4. Cómo se reparó

1. Detectar los ficheros con tipos `CG0`/`G0`.
2. Re-preparar el ligando desde su SMILES con `rigid_macrocycles=True`
   (el anillo se cierra) y **comprobar que el número de átomos pesados cuadra con
   el SMILES**. Cada reparado pierde exactamente sus 2 átomos duplicados.
3. Re-acoplar **solo esos ligandos** con el mismo receptor, caja y exhaustividad
   del conjunto original, para que el antes y el después sean comparables.
4. Nada se borra: poses viejas en `<out>/_antes_glue/`, ligandos viejos en
   `<ligands>/_roto_glue/`.

La prueba de que la reparación no toca la química real está en los propios
términos: en 4MQ y ZO0 la energía **intermolecular es idéntica** antes y después
(−6,17 vs −6,19 kcal/mol). Lo que desaparece es el término interno del anillo
abierto (+8,94 → −0,36 y +8,83 → −0,57). Detalle en
`confirmar_glue_controles.log`.

## 5. Barrido de las validaciones históricas

Antes de dar por cerrado el asunto, se revisó el resto de tablas de la carpeta:

| Tabla | Filas | Nombres con glue | Veredicto |
|---|---|---|---|
| `validacion_TDP43.csv` | 78 | 0 | limpia (top 5 sin glue) |
| `validacion_TDP43_RRM1.csv` | 42 | 0 | limpia |
| `validacion_TDP43_RRM2.csv` | 42 | 0 | limpia |
| `validacion_TDP43_v2.csv` | 198 | 0 | limpia |
| `validacion_TDP43_v3.csv` | 198 | 0 | limpia |
| `validacion_TDP43_v4.csv` | 198 | 0 | limpia |
| `validacion_TDP43_v5_4BS2.csv` | 58 | **6** | afectada (subconjunto de la v6) |
| `validacion_FUS.csv` (v1) | 35 | **1** | afectada (Sertralina, la nº 1) |
| `validacion_SOD1.csv` (v1) | 219 | **2** | afectada (los mismos dos del trp32) |

Las siete primeras usan conjuntos de ligandos que **ya no están en disco** (los
ficheros de `_validacion_TDP43/ligands` son otros: 133 nombres que no aparecen en
sus CSV). No se pueden reverificar; lo que sí se puede afirmar es que su top 5 no
contiene ninguno de los 30 nombres afectados conocidos, así que **no se deben
citar como contaminadas**. Si alguna se vuelve a usar, hay que re-preparar sus
ligandos desde SMILES con Meeko corregido.

### Conclusión del barrido
Las tablas contaminadas y ya corregidas son las **cinco vivas**: TDP43 v6, FUS v2,
SOD1 trp32 (v1), SOD1 v3 y, por herencia de ligandos, TDP43 v5.

## 6. Pendiente

- [ ] Arreglar la raíz en la preparación de ligandos para que no vuelva a pasar:
      llamar a Meeko con `rigid_macrocycles=True` en todos los scripts que
      preparen controles (hoy solo `reparar_ligandos_glue.py` lo hace).
- [ ] El rescoring MM-GBSA del top de TDP-43 sigue siendo el paso que decide.
- [ ] Decidir si v5 y las tablas viejas de nombre distinto se archivan o se
      rehacen; mientras tanto, la referencia viva de TDP-43 en el bolsillo 4BS2
      es la v6 reparada.
