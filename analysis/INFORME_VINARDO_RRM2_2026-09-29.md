# Vinardo contra los fragmentos de RRM2 (TDP-43) — 29 de septiembre de 2026

## La pregunta

El 29 de septiembre por la manana se concluyo que los tres fragmentos de Nshogoza
(11, 13 y 15 atomos pesados) caen al final del ranking de TDP-43, y que la causa
probable era la funcion de puntuacion: **Vina premia la superficie enterrada**, y un
fragmento pequeno no puede competir con un farmaco de 25 a 34 atomos.

La prueba que faltaba era repetir la puntuacion con **Vinardo**, que no premia la
superficie enterrada. Si el fallo fuera de la funcion, los fragmentos recuperarian
puesto. Este informe es el resultado de esa prueba.

## Como se ha hecho

**No se ha re-acoplado nada.** `_repunuar_vinardo.py` tomo las 576 poses ya
existentes de TDP-43 (5 positivos + 131 del fondo blando + 152 del fondo duro) y las
punto dos veces con `vina.exe` v1.2.3, `--score_only --scoring vina` y
`--scoring vinardo`. Son las mismas coordenadas: lo unico que cambia es el termino.
576 de 576, 0 fallos, 128 s.

`metricas_vinardo.py` lee ese CSV y calcula las metricas **importando `evaluar()`**
de `validar_sod1_limpia.py`, que es la unica copia de AUC, EF, residual y del
criterio. Copiar las cuentas seria volver a tener dos recetas, que es como este
proyecto ha perdido la mitad de sus conclusiones.

Nota sobre lo que Vinardo **no** se ha hecho: no se ha cambiado la funcion del
embudo. Vinardo esta documentada como peor que Vina en pruebas de redocking, asi que
cambiar el pipeline entero por ella solo porque mejora tres cifras no es una decision
que se tome de un solo CSV. Aqui se usa como **experimento diagnostico**, que es lo
que pedia la pregunta.

## Resultado 1: los puestos de los fragmentos

Contra el **fondo duro** (los 152 unidores de ARN de R-BIND, la lista de 159):

| ligando | atomos pesados | puesto con Vina | puesto con Vinardo | cambio |
|---|---|---|---|---|
| fragmento_1 | 13 | 155 de 159 | **155 de 159** | 0 |
| fragmento_2 | 11 | 159 de 159 | **158 de 159** | +1 |
| fragmento_3 | 15 | 157 de 159 | **159 de 159** | −2 |

Contra el fondo blando (122 senuelos, lista de 136) tampoco: 134→131, 136→135, 135→136.

**No recuperan nada.** El fragmento_3 va dos puestos PEOR con Vinardo.

## Resultado 2: el puesto dentro de su mismo tamano

El puesto global mezcla dos cosas: si la quimia se reconoce y si la molecula es
grande. Para quitar lo segundo, cada positivo se compara solo con los ligandos que
tienen **menos de 4 atomos pesados de diferencia** (la misma tolerancia que usa
`pareado_por_tamano` en el validador). Fondo duro:

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

Los medianos y grandes, mucho. En su propio grupo de tamano:

- **rTRD01**: de 48 de 80 a **24 de 80**
- **nTRD22**: de 1 de 68 a 5 de 68 (sigue siendo el primero de los grandes)
- **PE859**: de 17 de 40 a **10 de 40**
- **berberrubina**: de 44 de 69 a 55 de 69 (empeora)

Eso confirma que Vinardo quita parte del sesgo de tamano **cuando la molecula tiene
superficie que tapar. El fragmento no tiene esa superficie: quitarsela
no le devuelve nada, porque nunca se la會計 como premio.

## Resultado 4: las metricas no cambian de veredicto

Las seis celdas, con las mismas poses:

| fondo | funcion | AUC crudo | AUC/atomo | residual | EF5% | quimias al 5% | veredicto |
|---|---|---|---|---|---|---|---|
| blando | Vina | 0,350 | 0,746 | 0,417 | 0,00 | 0 | NO PASA |
| blando | Vinardo | 0,388 | 0,693 | 0,505 | 0,00 | 0 | NO PASA |
| duro | Vina | 0,313 | 0,724 | 0,439 | 0,00 | 0 | NO PASA |
| duro | Vinardo | 0,344 | 0,601 | 0,475 | 0,00 | 0 | NO PASA |
| los dos | Vina | 0,330 | 0,734 | 0,428 | 0,00 | 0 | NO PASA |
| los dos | Vinardo | 0,364 | 0,643 | 0,496 | 0,00 | 0 | NO PASA |

El AUC crudo sube un poco (0,313 → 0,344 contra el fondo duro), pero **el AUC por
atomo y el residual bajan**, que son las dos medidas que justamente quitan el tamano.
Ninguna de las seis pasa el criterio. TDP-43 sigue sin reconocerse por encima del
azar en ninguna funcion.

## Conclusion

**El problema no es la funcion de puntuacion, es el rango.** Queda medido, no supuesto:
cambiar a una funcion que no premia la superficie enterrada deja a los tres fragmentos
en el ultimo puesto y empeora uno de ellos. No hay una funcion de puntuacion que haga
ganar a un fragmento de 11 atomos contra 60 unidores de ARN de 29 a 45 atomos, porque
la pregunta que hace el cribado (¿esto separa TDP-43 de un unidor generico de ARN?)
la ganan las moleculas con superficie.

Esto tiene una consecuencia de estrategia que conviene decir sin rodeos: **el cribado
no puede usarse como filtro de fragmentos** con este fondo. Si el proyecto quiere
trabajar con fragmentos, tiene que cambiar el criterio antes que la funcion: por
ejemplo puntuar contra el **tercer fondo** (los compuestos medidos que NO unen, en
`_medidos_no_unen/`, hoy vacio porque faltan los ensayos) en vez de contra un fondo
de unidores de ARN, porque un fragmento pequeño solo tiene sentido si gana a lo que
no une nada, no a lo que une otra cosa.

## Ficheros

- `analysis/_repunuar_vinardo.py` — el re-puntador (Vina y Vinardo sobre las mismas poses).
- `analysis/metricas_vinardo.py` — las metricas con las dos funciones, importando `evaluar()`.
- `analysis/vinardo_tdp43.csv` — 288 ligandos con `vina`, `vinardo` y su puesto contra cada fondo.
- `analysis/metricas_vinardo_tdp43.csv` — las seis celdas de la tabla de arriba.
- `analysis/metricas_vinardo_tdp43.log` — la salida completa, con los seis bloques.
