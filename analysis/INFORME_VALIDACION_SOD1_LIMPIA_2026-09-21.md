# SOD1 contra el receptor limpio: el embudo reconoce química, pero no ordena por energía

> **Matizado el 24 de septiembre de 2026:** con la regla nueva (intervalo bootstrap sobre
> el AUC por átomo pesado) **SOD1 queda sin evidencia**, y el hallazgo que se sostiene es
> otro: ordena mejor que el azar contra el fondo duro de quelantes y redox, y **no** contra
> los señuelos de su mismo tamaño. Ver la corrección al principio de este informe.

**21 de septiembre de 2026.**
**Script:** `analysis/validar_sod1_limpia.py` · **Datos:** `validar_sod1_limpia.csv`,
`validar_sod1_limpia_resumen.csv`, `validar_sod1_limpia.log`

Es la validación que el criterio exigía antes de reportar nada: receptor corregido,
positivos de **unión medida** de verdad, series colapsadas y el criterio fijado por
adelantado. Y da un resultado partido, que es más útil que un número redondo.

---

## Corrección del 24 de septiembre de 2026 (regla nueva)

Este informe daba **PASA** apoyándose en la medida emparejada por tamaño (5 de 8
quimiotipos, §4). Esa medida **ya no vale como veredicto**: es la misma familia de
criterios que en TDP-43 dio 1, 2 o 3 quimias con el mismo conjunto, solo por cambiar el
motor o el esfuerzo (`INFORME_BARRIDO_EXHAUSTIVIDAD_TDP43_2026-09-24.md`). El sustituto es
el **AUC por átomo pesado con intervalo del 95 % por bootstrap**, exigiendo que el límite
inferior pase de 0,5 al remuestrear solo el fondo y al remuestrear fondo y positivos, y
que el veredicto no dependa de la corrida (`REGLA_DECISION_2026-09-24.md`).

**Con la regla nueva, SOD1 queda SIN EVIDENCIA, y no en un bloque: en todos.** Además, el
fondo de 482 no es un fondo, son tres pegados: al partirlos por el prefijo del fichero
(`separar_fondos_sod1.py`; `DEC_` y `DECM_` emparejados, `DECH_` duros) se ve qué pregunta
tenía respuesta y cuál no.

| Bloque | AUC/átomo | IC 95 % solo fondo | IC 95 % fondo + positivos | Residual | Veredicto |
|---|---|---|---|---|---|
| fondo entero (482, como en §2) | 0,563 | [0,540; 0,586] | [0,369; 0,749] | 0,478 | SIN EVIDENCIA |
| señuelos emparejados (339) | **0,530** | [0,501; 0,558] | [0,342; 0,720] | 0,454 | SIN EVIDENCIA |
| **fondo duro (143)** | **0,643** | [0,608; 0,676] | [0,435; 0,837] | **0,567** | SIN EVIDENCIA |

### Qué se puede reportar de SOD1, y qué no

**Sí se puede decir, con su error:**

- **La lista ordenada por energía no sirve.** AUC crudo 0,439, y su intervalo
  [0,419; 0,461] queda **entero por debajo de 0,5**. Antes se decía el número; ahora se
  puede decir con el error.
- **Contra el fondo duro (quelantes de metales y redox) el embudo ordena mejor que el
  azar**: 0,643, con el intervalo solo-fondo [0,608; 0,676] entero por encima de 0,5. Es
  el tercer mejor bloque de los 16 medidos, detrás de las cuatro corridas de TDP-43 sobre
  su fondo duro. Y es el **único de los 16** donde el AUC **residual** también supera el
  umbral con ese intervalo (0,567; [0,510; 0,616]).
- **Contra los señuelos emparejados por tamaño, no**: 0,530, al borde (límite inferior
  0,501, el mínimo de los 16 bloques). Ahí el embudo no distingue — y era justo el bloque
  que decidía el veredicto viejo.
- **El confusor de tamaño sigue medido**, y no depende de ningún criterio: −0,10 kcal/mol
  por átomo pesado contra los emparejados y −0,14 contra los duros, correlación −0,72, y
  la cabeza de la lista sigue siendo grande (§2).

**Ya NO se puede decir:**

- *"El embudo reconoce química cuando se le pregunta en igualdad de tamaño, y es
  reportable"* (§5). Lo que hay es **el número de un bloque, sin intervalo que lo cierre**:
  con 11 positivos, para que el límite inferior tocara el umbral habría que llegar a
  **~22 positivos** contra el fondo duro, **~99** contra el fondo entero y **~444** contra
  los emparejados. El cuello de botella son los positivos, no la exhaustividad.
- El **veredicto PASA** de §4. Queda como lo que fue: un criterio permisivo que se cumplía
  con 5 de 8 quimias y que no aguanta el cambio de condiciones.

**Lo que no cambia:** la lista de candidatos de SOD1 sigue sin publicarse, y el defecto
localizado —el término de tamaño de la función de puntuación— sigue siendo el mismo y
sigue siendo lo que hay que atacar. Lo que cambia es que SOD1 **no tiene ninguna
afirmación demostrada**: tiene una prometedora (el fondo duro, 0,643) que pide once
positivos más para poder decirse.

---

## 1. Lo que se montó

| Pieza | Qué es |
|---|---|
| Receptor | `gpu_dock/SOD1_limpio.pdbqt` — un dímero biológico (A–H), sin las 1.763 aguas del cristal ni las 18 copias de 1HL5 |
| Caja | Trp32, centro 46,5/80,0/73,3; 22 Å; exhaustividad 8 (la misma de toda la serie) |
| Positivos | **11 de unión medida** y **8 quimiotipos** distintos (`verdad_de_referencia.csv`) |
| Fondo | **482 señuelos**: 199 emparejados antiguos + 140 emparejados nuevos + 143 duros (quelantes de metales y redox-activos) |
| Preparación | receta canónica (`preparar_ligando.py`) |

**Comprobación previa: los 11 positivos pasan el control de integridad** — átomos
pesados del fichero, de la pose y del SMILES coinciden en los once (15/15, 13/13,
11/11, 18/18, 17/17, 18/18, 22/22, 10/10, 31/31, 33/33, 23/23). Nada que reparar.

Los ocho quimiotipos: catecolamina (3), quinazolina (2), nucleósido, anilina,
quinazolina-CF₃, benzisoxazol-piperidina, fenantridinona y aminoalcohol naftalénico.

## 2. Las métricas, y la trampa del tamaño

| Medida | Valor |
|---|---|
| AUC crudo | **0,439** |
| AUC por átomo pesado | 0,564 |
| AUC residual (descontada la recta del tamaño) | **0,477** |
| EF1 % / EF5 % | 0,00 / 3,64 |

Un AUC crudo **por debajo del azar** con once positivos no es un accidente: es un
síntoma. Y la causa está medida:

- **Dentro del propio fondo, la afinidad correlaciona −0,720 con el número de átomos
  pesados.** Un carbono pesado vale −0,106 kcal/mol. O sea: en esta caja, Vina ordena
  tamaño, casi literalmente.
- **La cabeza de la lista son moléculas grandes.** Los 25 mejores tienen **27,1** átomos
  pesados de media, frente a **19,7** del fondo y 19,2 de los positivos.

Y los positivos son **bimodales**: los tres grandes (33, 31 y 23 átomos) quedan en el
2 %, 3 % y 2 % de cabeza, y los ocho restantes se van al fondo de la lista —
isoproterenol 412, adrenalina 421, dopamina 422, anilina 459. Ninguno de los pequeños
es malo por su química: son pequeños.

## 3. El criterio, medido de tres maneras

El criterio —*positivos de al menos dos quimiotipos distintos por delante de los
señuelos emparejados*— no fija un corte, así que se midió de las tres formas posibles.
**No dan lo mismo, y ahí está la información:**

| Forma de medirlo | Quimiotipos que pasan |
|---|---|
| A) corte fijo de cabeza, 5 % | 2 de 8 (y **0** al 1 %, 1 al 2 %, 3 al 10 %) |
| B) AUC residual (descuenta el tamaño con una recta) | AUC 0,477 → **ninguno** |
| C) **emparejado por tamaño** (cada positivo contra señuelos de su mismo tamaño) | **5 de 8** |

La A no sirve: el veredicto salta de 0 a 3 según dónde se ponga el corte, y un
criterio así no mide nada. La B tampoco: una recta no captura el efecto, que es
bimodal. **La C es la que responde a la pregunta**: si el ranking ordenara química, un
positivo tiene que ganar a los señuelos **del mismo número de átomos**; si solo gana a
los pequeños, lo que ordena es el tamaño.

## 4. El resultado

| Quimiotipo | Positivo (átomos) | AUC emparejada por tamaño | ¿Gana a los de su tamaño? |
|---|---|---|---|
| aminoalcohol naftalénico | 946 (23) | **0,993** | sí |
| benzisoxazol-piperidina | K4I (31) | **0,806** | sí |
| catecolamina | dopamina (11) | **0,686** | sí |
| anilina | ZZT (10) | **0,606** | sí |
| fenantridinona | 6B3 (33) | **0,548** | sí |
| quinazolina-CF₃ | ZO0 (22) | 0,486 | no |
| quinazolina | 12I (17) | 0,394 | no |
| quinazolina | 4MQ (18) | 0,206 | no |
| nucleósido | 5-fluorouridina (18) | 0,243 | no |
| catecolamina | adrenalina (13) | 0,437 | no |
| catecolamina | isoproterenol (15) | 0,244 | no |

**VEREDICTO: PASA — 5 quimiotipos distintos de 8** baten a los señuelos de su mismo
tamaño *(superado el 24 de septiembre de 2026: este criterio ya no vale como veredicto y,
con la regla nueva, SOD1 queda SIN EVIDENCIA — ver la corrección al principio del
informe)* (aminoalcohol naftalénico, benzisoxazol-piperidina, catecolamina, anilina y
fenantridinona). El criterio pedía dos; se cumplen cinco, y con once positivos de ocho
quimias la exigencia no se puede acusar de tramposa.

### Y lo que NO dice, en el mismo párrafo

1. **La lista ordenada por energía, tal cual, no sirve para elegir candidatos.** Su
   cabeza es grande, no buena (AUC crudo 0,439, por debajo del azar). Si se usa para
   triaje, hay que leerla **dentro de estratos de tamaño**, nunca por el puesto
   absoluto y nunca con un corte fijo.
2. El quimiotipo de **catecolamina** pasa **por un solo miembro de tres**: dopamina sí
   (0,686), adrenalina (0,437) e isoproterenol (0,244) no. Si se exigiera la mediana en
   vez del mejor miembro, pasarían 4 de 8 — se sigue cumpliendo el criterio, pero el
   dato hay que darlo entero.
3. **Las tres quinazolinas fallan**, y no por casualidad: son exactamente las que
   tampoco recuperan su pose cristalina (ZO0 recupera en el modo 9, no en el que gana
   por energía; 12I y 4MQ no la recuperan en ninguno). Falla la misma cosa en los dos
   frentes: la función de puntuación.
4. **El 0,80 histórico era la serie pirazolona**, no SOD1. Con los 18 análogos dentro,
   el AUC era 0,797–0,803; con once positivos independientes y de ocho quimias es
   0,439 crudo. La diferencia entre esos dos números es exactamente lo que la auditoría
   predijo.

## 5. Qué significa para el proyecto

- **La diana no está muerta, y el sitio tampoco.** *(Corregido el 24 sep 2026.)* Lo que
  se sostiene hoy es más estrecho y está medido con intervalo: el embudo ordena mejor que
  el azar contra el **fondo duro** (0,643; IC [0,608; 0,676]) y **no** contra los
  señuelos de su mismo tamaño (0,530; IC [0,501; 0,558]). Decir que "reconoce química"
  era apoyarse en el criterio emparejado, que ya no decide.
- **El defecto de ordenación es el que está localizado** —el término de tamaño de la
  función de puntuación—, y es lo que hay que atacar ahora. *(Corrección del 24 sep 2026:
  el reconocimiento, además, **no está demostrado** — ver la corrección.)* Hay dos vías
  concretas:
  1. **Rescoring MM-GBSA dentro de estratos de tamaño**, que es donde el tamaño ya no
     confunde; o
  2. una normalización no lineal (la recta no basta: el residual baja a 0,477).
- **Y lo que NO se hace:** presentar el top de la lista como candidatos. Con el AUC
  crudo por debajo del azar, esa lista elige moléculas grandes.

## 6. Siguiente paso, concreto

1. **Repetir esta misma validación en TDP-43** con los 7 positivos de unión medida y un
   fondo difícil (R-BIND 2.0), aplicando la medida emparejada por tamaño como criterio
   principal. Es el mismo script con otro receptor y otra tabla.
2. **Atacar el término de tamaño** antes de volver a ordenar cualquier lista: probar el
   rescoring dentro de estrato y ver si el AUC emparejado se mantiene.
   **Medido el 22 sep 2026** (`INFORME_RESCORING_ESTRATOS_2026-09-22.md`): el MM-GBSA pasa
   la pregunta emparejada mejor que las funciones clásicas (AUC 0,566 crudo y 0,633
   estratificado sobre los mismos 82 ligandos, 6 de 8 quimias, y recupera 2 de las 3
   quinazolinas que fallaban), pero **no quita el sesgo de tamaño**: su pendiente por
   átomo pesado (−0,482 kcal/mol) es la más fuerte de las cuatro.
3. La lista de candidatos de SOD1 **no se publica** hasta que la ordenación no sea la
   del tamaño.
