# Validación de SOD1 rehecha: conjunto limpio de controles y dos fondos

**20 de septiembre de 2026.** Script: `validar_sod1_v3.py`. Datos:
`validacion_SOD1_v3/analisis_sod1_v3.csv`, `detalle_activos_sod1_v3.csv`,
`listas_fondos.csv`, `fondo_duro_motivos.csv`; registro
`validacion_SOD1_v3.log`. Controles de partida: `controles_sod1_v3.csv` (con sus
notas de corrección en `controles_sod1_v3_NOTAS.md`).

---

## 0. Resumen

1. **Se rehízo la validación con el conjunto limpio y dos fondos**: señuelos
   emparejados en propiedades (140 nuevos + los 199 antiguos, ya acoplados) y un
   **fondo duro** nuevo de 143 compuestos quelantes de metales o redox-activos
   (el tipo de química que de verdad se pegaría a una metaloenzima).
2. **El sitio y el protocolo están validados**: el control de redocking pasa en
   isoproterenol (0,48–1,28 Å) y adrenalina (0,69–0,72 Å) — ver
   `redocking_trp32/INFORME_CONTROLES_CORREGIDOS_2026-09-20.md`.
3. **El ranking NO está validado, y ahora se sabe por qué.** Los cuatro ligandos
   co-cristalizados —los que definen el bolsillo— puntúan **peor que el azar**
   frente a los señuelos (AUC 0,32–0,43). La causa medida es el **sesgo de
   tamaño**: en esta caja cada átomo pesado adicional vale **0,096 kcal/mol**.
   Descontado el tamaño, esos mismos ligandos suben a 0,62–0,75.
4. **El 0,815 se explica entero por eso**: 16 de sus 20 activos eran la misma
   serie pirazolona, moléculas grandes y grasas que heredan la ventaja de tamaño.
   Un solo representante de esa serie da AUC 0,93–0,98 en cualquier fondo.
5. **Los dos únicos activos independientes están en el azar**: LCS-1 y PRG-A01
   dan AUC 0,60–0,67 sin normalizar y **0,46–0,57 normalizado por tamaño**; LCS-1
   queda en el puesto 313 de 482.
6. **Ningún activo entra en el 1% superior** (EF1% = 0 en todas las familias).

---

## 1. Activos: el conjunto limpio, separado por tipo de evidencia

De `controles_sod1_v3.csv`, sin contar series congénicas más de una vez:

| Familia | Compuestos | Evidencia |
|---|---|---|
| independientes | LCS-1 (piridazinona), PRG-A01 (cumarina) | actividad celular; las dos únicas quimias independientes |
| serie colapsada | CHEMBL2165613 (representante pirazolona) | serie ChEMBL2165601-2165614 + 1643541/56/57, colapsada a un compuesto |
| co-cristalizados | isoproterenol, adrenalina, dopamina, 5-fluorouridina | **co-cristalización** (PDB 4A7T, 4A7U, 4A7V, 4A7S): evidencia estructural directa |

La separación no es cosmética: cada familia responde a una pregunta distinta.
Los co-cristalizados validan el **sitio** (¿el bolsillo es el que dice la
literatura y el protocolo lo reproduce?). Los independientes preguntan si el
método **ordena quimias distintas**.

## 2. Dos fondos

| Fondo | Cómo se construyó | N |
|---|---|---|
| emparejado, nuevo | `generar_decoys` (mismo peso ±30, logP ±1, rotables ±2, Tanimoto < 0,35) para los siete activos limpios | 140 |
| emparejado, antiguo | los 199 señuelos de la validación anterior, ya acoplados en la misma caja (mismos parámetros) | 199 |
| **duro** | compuestos de la misma librería con grupos que **quelan metales o son redox-activos** (catecol, pirogalol, quinona, o-aminofenol, hidroxamato, tiol, sulfonato, fosfonato, hidrazida, tiosemicarbazona), ordenados por parecido de peso al de los activos | 143 |

Un ligando no llegó a acoplarse (`DECM_CHEMBL443052`), así que el fondo
emparejado nuevo quedó en 140 de 141.

Todo se acopló en la caja del proyecto (Trp32, 46,5 80,0 73,3; 22 Å),
exhaustividad 8, Vina 1.2.3 en CPU, y los ligandos preparados con el mismo
pipeline que los señuelos (ETKDGv3 + MMFF + meeko).

## 3. Resultado: AUC, y AUC sin el tamaño

**El sesgo de tamaño, medido**: ajustando la afinidad frente al número de átomos
pesados en el fondo entero (N = 482),

```
afinidad = -3,42 - 0,096 * (átomos pesados)      r² sobre el fondo
```

es decir, **cada átomo pesado de más vale 0,096 kcal/mol** en esta caja. Con eso
se define el «AUC sin tamaño»: el mismo cálculo sobre el residuo del ajuste (la
parte del score que no explica el tamaño). El ajuste se hace **solo con el
fondo**, nunca con los activos.

| Familia | Fondo | AUC | AUC sin tamaño |
|---|---|---|---|
| independientes (2) | emparejado nuevo | 0,671 | **0,486** |
| independientes (2) | emparejado antiguo | 0,628 | **0,566** |
| independientes (2) | duro | 0,598 | **0,462** |
| independientes (2) | los tres juntos | 0,632 | **0,509** |
| serie (1) | los tres juntos | 0,932 | 0,950 |
| **co-cristalizados (4)** | emparejado nuevo | **0,425** | **0,661** |
| **co-cristalizados (4)** | emparejado antiguo | **0,328** | **0,746** |
| **co-cristalizados (4)** | duro | **0,316** | **0,624** |
| **co-cristalizados (4)** | los tres juntos | **0,353** | **0,682** |
| todos (7) | los tres juntos | 0,515 | **0,671** |

EF1% = 0 y EF5% = 0 (o 5,71 en un caso) en todas las combinaciones: con estos
activos no hay enriquecimiento en la cabeza de la lista.

### Dónde queda cada control (frente a los 482 señuelos juntos)

| Control | Afinidad | Puesto | Mejor que |
|---|---|---|---|
| CHEMBL2165613 (serie) | −6,37 | 34 | 93,0 % del fondo |
| PRG-A01 (cumarina) | −6,19 | 44 | 90,9 % |
| 5-fluorouridina | −5,70 | 113 | 76,6 % |
| LCS-1 (piridazinona) | −4,96 | **313** | 35,2 % |
| adrenalina | −4,78 | 361 | 25,3 % |
| dopamina | −4,71 | 380 | 21,3 % |
| isoproterenol | −4,64 | **398** | 17,6 % |

**Los tres últimos son ligandos con estructura cristalográfica resuelta en este
mismo bolsillo**, y el ranking los coloca por debajo del 75–82 % de los señuelos.

---

## 4. Qué significa

1. **El sitio está bien elegido y el protocolo coloca bien a las catecolaminas**
   (redocking 0,5–1,3 Å). Eso no está en duda.
2. **La puntuación de Vina en esta caja no mide afinidad.** Antes de descontar el
   tamaño, los ligandos verificados puntúan **peor que moléculas al azar
   emparejadas por propiedades**; después de descontarlo, quedan en 0,62–0,75, que
   es «débil», no «validado». Y los dos independientes se quedan en el azar en
   ambos casos.
3. **El 0,815 tiene ahora una explicación cuantitativa**: es una serie congénica
   de moléculas grandes, y en esta caja el tamaño se paga a 0,096 kcal/mol por
   átomo pesado. El número no mide reconocimiento molecular: mide tamaño.
4. **Ningún activo en el 1 % superior** significa que el corte «top 1 %» del
   cribado no contiene ni a los ligandos cuya pose está resuelta. Cualquier lista
   ordenada con esta caja y este protocolo sigue siendo **provisional**, y ahora
   por dos motivos independientes y medidos: controles positivos sesgados y
   puntuación sesgada por tamaño.
5. **Lo que sí sirve**: la caja como **región de triaje** (los ligandos reales
   están en el bolsillo, y el protocolo los coloca), y el fondo duro como
   herramienta reutilizable para cualquier validación de una metaloenzima.

---

## 5. Qué sigue

1. **Re-puntuar con una métrica sin tamaño** antes de volver a mirar cualquier
   lista: eficiencia de ligando (score / átomos pesados) o el residuo del ajuste
   que se ha usado aquí. Ya está implementado y es gratis sobre los datos que hay.
2. **Volver a etiquetar los candidatos de SOD1** con esa métrica y comprobar si
   alguno de los controles verificados entra en la cabeza.
3. **Buscar positivos independientes** (quimias fuera de la serie pirazolona) con
   ensayo de unión directo, no celular: con dos activos independientes el AUC no
   se puede sostener, dé lo que dé.
4. **No presentar ninguna lista de SOD1 como validada** hasta que (1) y (3)
   estén hechos.
