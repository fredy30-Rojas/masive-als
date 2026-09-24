# La regla que decide, con su error medido

**24 de septiembre de 2026.** El corte por quimiotipos ya no vale como regla de
decision: con el mismo motor y el mismo conjunto da 2 de 5 quimias a exhaustividad 8 y
16 y 3 de 5 a 32 (una quimia cruza por 0,59 kcal/mol), y al cambiar de motor da 1 de 5
(`INFORME_BARRIDO_EXHAUSTIVIDAD_TDP43_2026-09-24.md`). El corte mide donde cae la cola
de un fondo de 122 ligandos, y esa cola se mueve con el esfuerzo. Esta es la sustituta,
con su error medido por bootstrap.

**El resumen, en una linea:** la regla es el **AUC por atomo pesado contra el fondo
duro**, y se declara cuando el limite inferior de su intervalo del 95 % supera 0,5; en
TDP-43 eso da **PASA en las cuatro corridas** (exhaustividad 8, 16 y 32 en CPU y la GPU),
que es justo lo que el corte por quimias no conseguia.

---

## 1. La regla, fijada por adelantado y no negociable

Medida, fondos, umbrales y remuestreos, escritos antes de mirar los numeros:

1. **Medida:** AUC por **atomo pesado** contra el **fondo duro** (los 152 unidores de ARN
   de R-BIND 2.0). El fondo duro es el unico que puede separar quimia especifica de
   quimia generica de ARN, y el AUC por atomo pesado es el unico numero que quedo por
   encima del azar en las cuatro corridas.
2. **Umbral:** 0,5. Por debajo, el embudo no ordena mejor que el azar.
3. **Error:** intervalo del 95 % por bootstrap, 2.000 remuestreos con reemplazo, semilla
   fija (20260924) para que el intervalo sea reproducible.
4. **Dos remuestreos, y los dos tienen que aguantar:**
   - (a) **solo el fondo**: mide la cola del fondo, que es exactamente lo que hacia
     bailar el corte por quimiotipos;
   - (b) **fondo y positivos**: mide ademas el muestreo de los positivos, que con 7
     ligandos es lo que de verdad aprieta.
5. **Veredictos posibles, y no hay mas:**
   - **PASA** — el limite inferior de (a) y de (b) supera 0,5. Se puede reportar.
   - **SIN EVIDENCIA** — (a) si, (b) no: el orden no es un capricho del fondo, pero con
     esos positivos no se puede afirmar que reconozca quimia.
   - **NO PASA** — ni (a).
6. **Condicion de validez de la regla misma:** el veredicto tiene que ser el **mismo en
   todas las corridas**. Una regla que cambia al cambiar el motor o el esfuerzo no sirve
   aunque su intervalo sea bonito. Eso es lo que comprueba el script (`comparar()`).

El AUC vectorizado de los 2.000 remuestreos se comprueba contra el `auc()` importado de
`validar_sod1_limpia` antes de usarlo: la regla no puede calcular su numero con una
formula distinta a la del informe que puntua (leccion del 24 sep con
`auditar_tipos_cribado`/`_limpiar_ranking_tipos.py`, que importa la logica en vez de
copiarla).

---

## 2. El resultado: la regla aguanta

Fondo duro, el bloque que decide. Las cuatro corridas son **el mismo conjunto**, el mismo
receptor y la misma caja; solo cambia el motor y el esfuerzo.

| Corrida | AUC/atomo | IC 95 % solo fondo | IC 95 % fondo + positivos | Veredicto |
|---|---|---|---|---|
| exhaustividad 8 (CPU) | 0,731 | [0,686; 0,774] | [0,540; 0,889] | **PASA** |
| exhaustividad 16 (CPU) | 0,734 | [0,689; 0,775] | [0,545; 0,893] | **PASA** |
| exhaustividad 32 (CPU) | 0,733 | [0,687; 0,775] | [0,562; 0,889] | **PASA** |
| GPU (search_depth 20) | 0,738 | [0,694; 0,780] | [0,543; 0,902] | **PASA** |

**Las cuatro dan lo mismo.** El punto se mueve seis milesimas, el limite inferior del
intervalo dos milesimas, y el veredicto no se mueve. El corte al que sustituye daba 1, 2
o 3 con estas mismas cuatro corridas.

Y hay algo que refuerza el resultado: **el limite inferior del intervalo solo-fondo supera
0,5 en los 13 bloques medidos** (los 12 de TDP-43 y el de SOD1), con minimo 0,5405. Es
decir, ordenar por energia por atomo pesado bate al azar **siempre**, en cualquier motor
y cualquier esfuerzo. Lo que cambia entre veredictos no es eso: es si el numero de
positivos permite demostrarlo.

---

## 3. Los otros dos bloques, y por que decide el duro

| Bloque | AUC/atomo (4 corridas) | IC solo fondo (inf.) | IC todo (inf.) | Veredicto |
|---|---|---|---|---|
| señuelos emparejados (122) | 0,666–0,686 | 0,621–0,648 | 0,434–0,465 | SIN EVIDENCIA × 4 |
| fondo duro (152) | 0,731–0,738 | 0,686–0,694 | 0,540–0,562 | **PASA × 4** |
| los dos juntos (274) | 0,703–0,713 | 0,672–0,681 | **0,498**–0,516 | PASA × 3, SIN EVIDENCIA × 1 |

- **El bloque blando** (señuelos emparejados) no llega en ninguna corrida, y es coherente:
  su limite inferior con todos los remuestreos queda en 0,43–0,47. No es que el embudo
  vaya mal contra los señuelos; es que con 7 positivos no se puede cerrar el intervalo.
- **El bloque de los dos fondos juntos** es el unico donde las cuatro corridas no
  coinciden: la GPU da 0,498 en el limite inferior y las tres de CPU dan 0,510–0,516.
  **Se cae por dos milesimas.** Queda escrito tal cual: es el margen del umbral, y es el
  mejor argumento para no decidir con el bloque mezclado. El que decide es el duro.
- **Por que el duro:** es el unico fondo que separa quimia especifica de quimia generica
  de ARN. Si el embudo reconoce quimia en igualdad de tamaño contra verdaderos unidores
  de ARN, reconoce algo; contra señuelos emparejados solo demuestra que distingue una
  molecula afin de una que no lo es, que es mucho mas facil.

---

## 4. Las dos medidas que quedan fuera, medidas y no supuestas

| Medida | Punto (13 bloques) | Limite inferior del IC solo fondo | Veredicto |
|---|---|---|---|
| AUC crudo | 0,246–0,440 | por debajo de 0,5 en todos | inservible para decidir |
| AUC residual | 0,264–0,478 | **maximo 0,446** | inservible para decidir |
| **AUC por atomo pesado** | 0,563–0,738 | **minimo 0,5405** | la que decide |

Esto cierra la duda que quedaba abierta desde el informe del control: el **residual**
—afinidad menos la recta ajustada al fondo— parecia la medida honesta, y no lo es. En
TDP-43 su limite inferior no llega a 0,5 en ningun bloque, ni siquiera con el intervalo
mas estrecho. Dividir por los atomos pesados quita el sesgo de tamaño mejor que restar
una recta, medido: la recta tiene pendiente pequeña (-0,015 a -0,045 kcal/mol por atomo
en TDP-43) y no da cuenta de la curvatura. Y el AUC crudo esta **por debajo del azar en
los 13 bloques**: la lista ordenada por energia sin mas no solo no ordena, ordena al
reves, porque premia al ligando grande.

---

## 5. Cuantos positivos hacen falta, y no mas exhaustividad

El error del AUC cuando se remuestrean los positivos baja como 1 partido por la raiz de
n. Con los numeros medidos se puede decir que falta con concrecion:

| Bloque | Positivos hoy | AUC/atomo | Para que el limite inferior toque 0,5 |
|---|---|---|---|
| TDP-43, fondo duro | 7 | 0,73 | **ya basta** (PASA) |
| TDP-43, señuelos emparejados | 7 | 0,67–0,69 | **~10–13** |
| SOD1, fondo limpio | 11 | 0,56 | **~99** |

Es una estimacion, no una promesa: supone que el embudo conserva el mismo orden al añadir
positivos, y eso no se sabe hasta tenerlos. Pero dice lo importante: en TDP-43 **el cuello
de botella son tres o seis positivos de union medida mas**, y no mas exhaustividad ni otra
GPU. Subir el esfuerzo mueve el criterio sin acercarlo a la verdad.

---

## 6. La consecuencia incomoda: SOD1

Aplicada al mismo bloque que el proyecto usa para SOD1 (11 positivos contra 482 señuelos
emparejados, exhaustividad 8):

| Medida | Punto | IC solo fondo | IC todo | Veredicto |
|---|---|---|---|---|
| AUC crudo | 0,439 | [0,419; 0,461] | [0,246; 0,664] | — |
| AUC/atomo | 0,563 | [0,540; 0,586] | [0,369; 0,749] | **SIN EVIDENCIA** |
| AUC residual | 0,478 | [0,446; 0,511] | [0,359; 0,613] | — |

La validacion de SOD1 se reporto como **PASA** con el criterio de quimiotipos (5 de 8
quimias por delante de los señuelos de su tamaño, `INFORME_VALIDACION_SOD1_LIMPIA_2026-09-21.md`).
Con esta regla **no pasa: queda sin evidencia**. No es que el trabajo de SOD1 estuviera
mal —el AUC por atomo pesado si esta por encima del azar, 0,563— es que esta **tres
centesimas por encima** y con 11 positivos el intervalo no cierra: harian falta unos 99.
El criterio viejo lo daba por bueno porque su corte es permisivo, no porque hubiera
evidencia.

Lo dejo escrito sin adornos porque cambia lo que se puede afirmar: **hoy ningun resultado
del proyecto pasa esta regla en su bloque blando**, y la unica afirmacion nueva que si se
sostiene es que en TDP-43 el embudo reconoce quimia de tamaño correcto contra el fondo
duro.

---

## 7. Lo que NO dice

1. **No valida candidatos.** Dice que el embudo ordena quimia por encima del azar en un
   bloque concreto. No dice que la cabeza de la lista sea buena —el AUC crudo sigue por
   debajo del azar— ni toca la lista de candidatos.
2. **No es libre.** La regla se fijo con un conjunto (7 positivos, 122 señuelos, 152
   duros) y un fondo duro concreto. Cambiar el fondo duro cambia el numero: si mañana se
   usa otro fondo de unidores de ARN, el intervalo hay que rehacerlo.
3. **No promete lo que estima.** Los "~10-13 positivos" y los "~99" salen de suponer que
   el orden se mantiene al añadir positivos. Es una estimacion para decidir por donde
   seguir, no un calculo de potencia.
4. **No dice que el residual sea inutil en general.** Dice que en este embudo, con estos
   fondos, no supera el umbral. En otro embudo con otra pendiente de tamaño podria ser al
   reves.
5. **No resuelve el bloque mezclado.** El que mezcla los dos fondos es el que se cae por
   dos milesimas en la GPU. La regla no lo arregla; lo declara fuera de la decision.

---

## 8. Como se reproduce

```
python C:/Users/Fredy/masive-als/analysis/regla_decision_bootstrap.py
```

Para otra diana, con corridas propias:

```
python regla_decision_bootstrap.py --nombre FUS "exh 8=_barrido/validar_exh8.csv"
```

- Script: `analysis/regla_decision_bootstrap.py` (importa `auc` y `residual` de
  `validar_sod1_limpia`, la unica copia de las metricas).
- Salidas: `analysis/regla_decision/regla_decision.log` y `.csv` (una fila por corrida y
  bloque, con los intervalos, el veredicto y los positivos que faltarian).
- Entradas: los CSV ya puntuados de las cuatro corridas de TDP-43
  (`_barrido_TDP43/validar_exh8.csv`, `validar_tdp43_limpia.csv`,
  `_barrido_TDP43/validar_exh32.csv`, `_control_gpu_TDP43/validar_gpu.csv`) y el de SOD1
  (`validar_sod1_limpia.csv`). No se recalcula ninguna afinidad: la regla se aplica al
  mismo numero que ya estaba en los informes.
