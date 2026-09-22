# SOD1 contra el receptor limpio: el embudo reconoce química, pero no ordena por energía

**21 de septiembre de 2026.**
**Script:** `analysis/validar_sod1_limpia.py` · **Datos:** `validar_sod1_limpia.csv`,
`validar_sod1_limpia_resumen.csv`, `validar_sod1_limpia.log`

Es la validación que el criterio exigía antes de reportar nada: receptor corregido,
positivos de **unión medida** de verdad, series colapsadas y el criterio fijado por
adelantado. Y da un resultado partido, que es más útil que un número redondo.

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
tamaño (aminoalcohol naftalénico, benzisoxazol-piperidina, catecolamina, anilina y
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

- **La diana no está muerta, y el sitio tampoco.** El embudo reconoce química cuando se
  le pregunta en igualdad de tamaño. Eso es un resultado positivo y es reportable.
- **El defecto que queda es de ordenación, no de reconocimiento**, y está localizado en
  una sola cosa: el término de tamaño de la función de puntuación. Es lo que hay que
  atacar ahora, y hay dos vías concretas:
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
