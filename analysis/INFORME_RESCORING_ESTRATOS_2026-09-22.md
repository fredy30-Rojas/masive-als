# MM-GBSA por estratos de tamaño: ordena mejor, pero no quita el sesgo de tamaño

**22 de septiembre de 2026.**
**Scripts:** `analysis/preparar_mmgbsa_estratos.py`, `analysis/mmgbsa_runner_estratos.py`
(en Oracle), `analysis/rescoring_estrategia_estratos.py` ·
**Datos:** `analysis/rescoring_estrategia_estratos.csv` (lista completa, N=513),
`rescoring_estrategia_estratos_muestra.csv` (los mismos 82 ligandos),
`mmgbsa_estratos/resultado_oracle.csv`

Responde a la pregunta que quedó abierta en `INFORME_VALIDACION_SOD1_LIMPIA_2026-09-21.md`:
si la ordenación por energía no sabe ordenar porque su energía es casi una función del
tamaño, **¿el rescoring físico MM-GBSA ordena mejor, y le queda el mismo vicio?**

---

## 1. Lo que se montó

| Pieza | Qué es |
|---|---|
| Receptor | `SOD1_limpio.pdb` (dímero A–H, sin aguas ni copias de 1HL5), el mismo de la validación |
| Poses | las del acoplado original (`validacion_SOD1_v5/out`), **no** se volvieron a acoplar |
| Positivos | los **11 de unión medida** de `verdad_de_referencia.csv` |
| Fondo | **79 señuelos**: 66 emparejados por tamaño (±2 átomos pesados con un positivo) + 13 repartidos por cuantiles de tamaño |
| Método | `rescoring_mmgbsa_robusto.py`: GAFF 2.11 + AM1-BCC (antechamber), openmm, OBC2, una pose por ligando |
| Máquina | Oracle, 3 procesos en paralelo, **108,5 min** para los 90 (antes eran ~83 s por ligando en serie) |

**Cobertura: 11 positivos + 71 señuelos con dG (60 emparejados y 11 repartidos), de 90
intentados.** Los 8 fallos son
**todos señuelos**, ninguno positivo, así que el veredicto no está sesgado; pero son 8
defectos que conviene arreglar (ver §5).

Cada positivo tiene entre **10 y 21 señuelos de su mismo tamaño** dentro de la muestra:
el criterio emparejado se puede medir con apoyo, no con anécdotas.

## 2. La comparación, sobre los MISMOS ligandos

Comparar un AUC medido sobre 82 ligandos emparejados con otro medido sobre 513 no vale:
la muestra emparejada es más fácil para cualquier función. Así que las cuatro funciones
se pasaron por los mismos 82 ligandos:

| Función | AUC crudo | AUC residual | Corr. con tamaño | EF5 % | Quimiotipos |
|---|---|---|---|---|---|
| Acoplado (Vina, control) | 0,453 | 0,481 | −0,825 | 3,64 | 4 de 8 |
| Rescoring Vina (`--score_only`) | 0,458 | 0,507 | −0,821 | 3,64 | 3 de 8 |
| Rescoring Vinardo | 0,465 | 0,510 | −0,633 | 1,82 | 6 de 8 |
| **Rescoring MM-GBSA** | **0,566** | **0,624** | −0,510 | 1,82 | 6 de 8 |

**El MM-GBSA es la mejor de las cuatro, y no por poco:** gana 0,10 de AUC crudo y 0,12
de residual a la que más se le acerca, y es la única que pasa de 0,6 en residual.

Y por estratos de tamaño, que es lo que se quería ver:

| Estrato | Positivos | Señuelos | MM-GBSA | Vinardo | Acoplado |
|---|---|---|---|---|---|
| 0–12 átomos | 2 | 10 | **0,950** | 0,850 | 0,650 |
| 13–17 átomos | 3 | 17 | 0,471 | 0,510 | 0,314 |
| 18–22 átomos | 3 | 17 | 0,373 | 0,196 | 0,216 |
| 23–27 átomos | 1 | 15 | 0,933 | 1,000 | 1,000 |
| 28–99 átomos | 2 | 12 | 0,625 | 0,542 | 0,667 |
| **Estratificado (ponderado)** | | | **0,633** | 0,584 | 0,533 |

## 3. El sesgo de tamaño NO se va: se hace más grande

Es el dato incómodo del día y hay que darlo entero. Dentro del fondo:

- Vina (acoplado): **−0,093 kcal/mol por átomo pesado**;
- Vinardo: −0,075;
- **MM-GBSA: −0,482 kcal/mol por átomo pesado.**

El MM-GBSA tiene **el término de tamaño más fuerte de las cuatro funciones**, casi cinco
veces el de Vina en kcal/mol por carbono. Lo que ocurre es que, aun arrastrándolo, dentro
de cada pareja de igual tamaño ordena mejor. O sea: **el MM-GBSA no arregla el defecto de
ordenación, lo compensa**. No es una normalización; es un filtro que acierta más.

## 4. El veredicto

**Sí se mantiene el AUC emparejado, y además mejora.** Con el MM-GBSA pasan **6 de 8
quimiotipos** (el criterio pedía dos):

| Quimiotipo | Positivo | dG (kcal/mol) | Átomos | ¿Gana a los de su tamaño? |
|---|---|---|---|---|
| aminoalcohol naftalénico | 946 | −23,24 | 23 | sí |
| anilina | ZZT | −16,87 | 10 | sí |
| benzisoxazol-piperidina | K4I | −26,44 | 31 | sí |
| catecolamina | isoproterenol / dopamina | −15,65 / −12,00 | 15 / 11 | sí |
| quinazolina | 12I | −15,86 | 17 | **sí** |
| quinazolina-CF₃ | ZO0 | −23,86 | 22 | **sí** |
| fenantridinona | 6B3 | −19,02 | 33 | **no** |
| nucleósido | 5-fluorouridina | −3,84 | 18 | **no** |

**Lo más informativo:** el MM-GBSA **recupera dos de las tres quinazolinas** — precisamente
las que fallaban el control de pose cristalina y que Vina no ordenaba. Eso refuerza el
diagnóstico anterior: el problema estaba en la **función de puntuación**, no en la química
de esos positivos. En contrapartida **pierde la fenantridinona**, que Vina sí ordenaba: no
es un barrido limpio, hay intercambio.

## 5. Lo que NO dice, en el mismo párrafo

1. **No se publica nada con esto.** Once positivos dan para decidir si una idea se sigue
   o se tira, no para reportar un número. La muestra se diseñó para la pregunta, no para
   dar una cifra.
2. **Los estratos son finos:** de 1 a 3 positivos en cada uno. El 0,633 estratificado es
   más creíble que el 0,566 crudo, pero los estratos grandes son los que pesan y el orden
   de los chicos es anécdota (el "1,000" del estrato 23–27 es **un** positivo).
3. **El MM-GBSA sigue teniendo el vicio de tamaño más fuerte de todas** (§3).
4. **Es un MM-GBSA de un punto:** una pose por ligando, sin minimización ni muestreo. No
   es el protocolo de referencia (que promedia), así que los dG no se comparan con la
   literatura; sirven para ordenar dentro de esta tabla y nada más.
5. **Los 8 fallos de preparación** (todos señuelos): 4 por **desacuerdo entre el número
   de átomos del SMILES y el de la pose**, 2 porque **obabel no leyó la pose**, y 2 por
   **"No template found for residue 35/47 (THR)"** — un defecto del receptor, no del
   ligando (mismo tipo de problema que el del calcio que se parcheó en agosto, y aquí en
   el residuo terminal). Los dos últimos son arreglables sin tocar el método.

## 6. Qué significa para el proyecto

- **El MM-GBSA se queda como segunda vuelta del embudo**, pero **dentro de cada estrato
  de tamaño**, no sobre la lista entera. La forma de usarlo: el ranking de Vina elige
  dentro de cada estrato, y el MM-GBSA desempata dentro del estrato. Al revés —MM-GBSA
  sobre la lista completa— vuelve a ordenar por tamaño, con el doble de pendiente.
- **Ninguna de las cuatro funciones ordena sola.** Las cuatro están por debajo o al
  borde del azar en crudo (0,453–0,566) y las cuatro pasan en la pregunta emparejada.
  El problema del proyecto no es "qué función se usa", es **cómo se pregunta**.
- **Siguiente paso concreto, en orden:**
  1. **Repetir la validación en TDP-43** con los 7 positivos de unión medida y el fondo
     duro de R-BIND 2.0, con la medida emparejada como criterio principal
     (`validar_diana_limpia.py` ya está preparado para eso). Es lo que falta para saber
     si esto es de SOD1 o del embudo entero.
  2. **Arreglar los 8 fallos** de preparación (el residuo THR del receptor y el desacuerdo
     SMILES/pose) antes de ampliar la muestra.
  3. Solo si las dos anteriores salen bien: pensar en la lista de candidatos. Hasta hoy,
     **sigue sin publicarse**.
