# Vinardo contra los fragmentos de RRM2 (TDP-43) y SOD1 — 29 de septiembre de 2026

## La pregunta

Por la manana se concluyo que los tres fragmentos de Nshogoza (11, 13 y 15 atomos
pesados) caen al final del ranking de TDP-43, y que la causa probable era la funcion de
puntuacion: **Vina premia la superficie enterrada**, y un fragmento pequeno no puede
competir con un farmaco de 25 a 34 atomos.

La prueba que faltaba era repetir la puntuacion con **Vinardo**, que no premia la
superficie enterrada. Si el fallo fuera de la funcion, los fragmentos recuperarian
puesto.

## Respuesta corta: no recuperan nada

## Como se ha hecho

**No se ha re-acoplado nada.** `repunuar_vinardo.py` toma las poses ya existentes y
las puntua dos veces con `vina.exe` v1.2.3, `--score_only` con `--scoring vina` y con
`--scoring vinardo`. Son las mismas coordenadas: lo unico que cambia es el termino.

| diana | positivos | fondo blando | fondo duro | tareas | fallos | tiempo |
|---|---|---|---|---|---|---|
| TDP-43 | 7 | 122 | 152 | 562 | 0 | 107 s |
| SOD1 | 11 | 482 | 152 | 1.290 | 0 | 151 s |

`metricas_vinardo.py` calcula las metricas **importando `evaluar()`** de
`validar_sod1_limpia.py`, que es la unica copia de AUC, EF, residual y del criterio.
Copiar las cuentas seria volver a tener dos recetas, que es como este proyecto ha
perdido la mitad de sus conclusiones.

Nota sobre lo que Vinardo **no** se ha hecho: no se ha cambiado la funcion del embudo.
Vinardo esta documentada como peor que Vina en pruebas de redocking, asi que cambiar
el pipeline entero por ella solo porque mejora tres cifras no es una decision que se
tome de un solo CSV. Aqui se usa como **experimento diagnostico**.

## Resultado 1: los puestos de los fragmentos (TDP-43)

Contra el **fondo duro** (los 152 unidores de ARN de R-BIND, lista de 159):

| ligando | atomos pesados | puesto con Vina | puesto con Vinardo | cambio |
|---|---|---|---|---|
| fragmento_1 | 13 | 155 de 159 | **155 de 159** | 0 |
| fragmento_2 | 11 | 159 de 159 | **158 de 159** | +1 |
| fragmento_3 | 15 | 157 de 159 | **159 de 159** | −2 |

Contra el fondo blando tampoco: 127→124, 129→128, 128→129.

**No recuperan nada.** El fragmento_3 va dos puestos PEOR con Vinardo.

## Resultado 2: el puesto dentro de su mismo tamano

El puesto global mezcla dos cosas: si la quimia se reconoce y si la molecula es grande.
Para quitar lo segundo, cada positivo se compara solo con los ligandos que tienen
**menos de 4 atomos pesados de diferencia** (la tolerancia que usa `pareado_por_tamano`
en el validador). Fondo duro de TDP-43:

| ligando | pesados | puesto global (Vinardo) | **puesto en su tamano** | con Vina |
|---|---|---|---|---|
| fragmento_1 | 13 | 155 de 159 | **11 de 14** | 10 de 14 |
| fragmento_2 | 11 | 158 de 159 | **10 de 11** | 11 de 11 |
| fragmento_3 | 15 | 159 de 159 | **14 de 14** | 13 de 14 |
| rTRD01 | 25 | 68 de 159 | 24 de 80 | 48 de 80 |
| nTRD22 | 29 | 25 de 159 | 5 de 68 | 1 de 68 |
| PE859 | 34 | 28 de 159 | 10 de 40 | 17 de 40 |
| berberrubine | 24 | 133 de 159 | 55 de 69 | 44 de 69 |

Esta es la tabla que contesta a Fredy. **El fragmento pierde tambien entre los suyos.**
No hay ningun subconjunto del fondo del que salga adelante: contra los ligandos de su
mismo tamano queda 10 de 11, 11 de 14 y 14 de 14. Y Vinardo les hace un poco mas de
daño, no menos.

## Resultado 3: a quien si ayuda Vinardo

Los medianos y grandes, mucho. En su propio grupo de tamano, en TDP-43:

- **rTRD01**: de 48 de 80 a **24 de 80**
- **PE859**: de 17 de 40 a **10 de 40**
- **nTRD22**: sigue siendo de cabeza (1 de 68 -> 5 de 68)
- **berberrubina**: de 44 de 69 a 55 de 69 (empeora)

Eso confirma que Vinardo quita parte del sesgo de tamano **cuando la molecula tiene
superficie que tapar**. El fragmento no la tiene: quitarsela no le devuelve nada,
porque nunca se le contabilizo como premio.

## Resultado 4: las metricas de las dos dianas

La columna Vina usa la afinidad que Vina **escribio al acoplar** (la misma que lee el
validador), para que las cifras sean comparables con las publicadas. Comprobacion de
que la receta no se ha movido: las tres celdas de Vina de TDP-43 salen **0,301 / 0,686**,
**0,357 / 0,734** y **0,332 / 0,713**, y las de SOD1 **0,440 / 0,564** y
**0,313 / 0,838**, que son exactamente las publicadas el 29 de septiembre por la
manana.

**TDP-43**

| fondo | funcion | AUC crudo | AUC/atomo | residual | EF5% | quimias 5% | veredicto |
|---|---|---|---|---|---|---|---|
| blando | Vina | 0,301 | 0,686 | 0,306 | 0,00 | 0 | NO PASA |
| blando | Vinardo | 0,384 | 0,684 | 0,502 | 0,00 | 0 | NO PASA |
| duro | Vina | 0,357 | 0,734 | 0,391 | 0,00 | 0 | NO PASA |
| duro | Vinardo | 0,344 | 0,601 | 0,475 | 0,00 | 0 | NO PASA |
| los dos | Vina | 0,332 | 0,713 | 0,356 | 0,00 | 0 | NO PASA |
| los dos | Vinardo | 0,362 | 0,638 | 0,493 | 0,00 | 0 | NO PASA |

**SOD1**

| fondo | funcion | AUC crudo | AUC/atomo | residual | EF5% | veredicto |
|---|---|---|---|---|---|---|
| blando | Vina | 0,440 | 0,564 | 0,478 | 3,64 | NO PASA |
| blando | Vinardo | 0,432 | 0,518 | 0,447 | 1,82 | NO PASA |
| duro | Vina | 0,313 | 0,838 | 0,461 | 0,00 | NO PASA |
| duro | Vinardo | 0,274 | 0,710 | 0,575 | 0,00 | NO PASA |
| los dos | Vina | 0,409 | 0,629 | 0,486 | 3,64 | NO PASA |
| los dos | Vinardo | 0,394 | 0,564 | 0,469 | 0,00 | NO PASA |

**Ninguna de las doce celdas pasa el criterio.** En SOD1 Vinardo va ademas peor en EF5:
3,64 pasa a 1,82 contra el blando y a 0 contra los dos fondos. En TDP-43 el AUC crudo
sube un poco contra el blando, pero **el AUC por atomo y el residual bajan**, que son
las dos medidas que quitan el tamano.

## Un hallazgo del camino: `--score_only` NO da la afinidad del docking

Al comparar la columna Vina del re-puntado con la que Vina escribio al acoplar, la
diferencia es de **1,96 kcal/mol de mediana** en TDP-43 (maximo 8,79) y **1,23** en SOD1
(maximo 8,09). Se descarto la causa facil antes de aceptar el numero:

- **No es la caja.** `diagnostico_score_only.py` puntua las mismas poses con la caja
  del docking y con la caja centrada en la pose: sale **exactamente lo mismo**
  (diferencia 0,000 en las 8 poses probadas).
- **No es un fallo del script.** El desplazamiento es parejo para los positivos y para
  los dos fondos, asi que no altera el ranking; altera las cifras absolutas.

La explicacion es que al acoplar Vina construye las mapas en una rejilla gruesa y
relaja la pose antes de puntuarla, mientras que `--score_only` evalua la funcion sobre
la pose guardada. Por eso el script **usa la afinidad del docking para la columna
Vina** de las metricas (y asi reproduce las cifras publicadas al milímetro) y se queda
con la de `--score_only` solo para la comparacion entre funciones, que es la pregunta
de los fragmentos y para la que las dos columnas tienen que salir del mismo sitio.

## Conclusion

**El problema no es la funcion de puntuacion, es el rango.** Queda medido, no supuesto:
cambiar a una funcion que no premia la superficie enterrada deja a los tres fragmentos
en el ultimo puesto y empeora uno de ellos, en las dos dianas y en los tres fondos.
No hay una funcion de puntuacion que haga ganar a un fragmento de 11 atomos contra 60
unidores de ARN de 29 a 45 atomos, porque la pregunta que hace el cribado (¿esto separa
TDP-43 de un unidor generico de ARN?) la ganan las moleculas con superficie.

Una observacion lateral que si es util: en SOD1 los ligandos **pequenos** (adrenalina de
13 pesados, dopamina de 11) quedan bien dentro de su grupo de tamano (1 de 16, 3 de 12).
O sea que un farmaco pequeno no tiene por que perder siempre: pierde cuando enfrente
hay algo del mismo tamano que de verdad une otra cosa. En TDP-43, al lado de los
fragmentos medidos hay 14 ligandos del fondo de 11 a 17 atomos que tambien son
unidores de ARN, y ahi el fragmento no tiene nada que ganar.

Esto tiene una consecuencia de estrategia que conviene decir sin rodeos: **el cribado
no puede usarse como filtro de fragmentos con este fondo**. Si el proyecto quiere
trabajar con fragmentos, tiene que cambiar el criterio antes que la funcion: por
ejemplo puntuar contra el **tercer fondo** (los compuestos medidos que NO unen, en
`_medidos_no_unen/`, hoy vacio porque faltan los ensayos) en vez de contra un fondo de
unidores de ARN, porque un fragmento pequeno solo tiene sentido si gana a lo que no une
nada, no a lo que une otra cosa.

## Ficheros

- `analysis/repunuar_vinardo.py` — el re-puntador, parametrizado por diana.
- `analysis/metricas_vinardo.py` — las metricas con las dos funciones, importando `evaluar()`.
- `analysis/comprobar_score_only.py` — mide la diferencia entre docking y `--score_only`.
- `analysis/diagnostico_score_only.py` — descarta la caja como causa de esa diferencia.
- `analysis/vinardo_tdp43.csv`, `analysis/vinardo_sod1.csv` — ligandos con `vina`, `vinardo` y su puesto contra cada fondo.
- `analysis/metricas_vinardo_tdp43.csv`, `analysis/metricas_vinardo_sod1.csv` — las doce celdas.
- `analysis/metricas_vinardo_tdp43.log`, `analysis/metricas_vinardo_sod1.log` — la salida completa.
