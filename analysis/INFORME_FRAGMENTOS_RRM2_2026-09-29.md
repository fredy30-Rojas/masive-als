# Por qué los tres fragmentos de RRM2 caen al fondo en TDP-43

**29 de septiembre de 2026.** En la puntuacion de TDP-43 contra el fondo duro, los tres
fragmentos de RRM2 salen en los puestos **152, 155 y 157 de 159**. Son la unica
evidencia estructural directa del sitio de union, asi que la pregunta era si el embudo
los hunde por su tamano o porque la caja no cubre donde se unen. Se han medido las dos
cosas, y **las dos quedan descartadas**.

## 1. La caja SI cubre el sitio de RRM2

Centro de la caja `(24.23, 16.89, -15.87)`, lado 26 A, o sea cubre +/-13 A. Los cinco
residuos que la verdad de referencia cita como sitio de los fragmentos estan **todos
dentro**:

| residuo | distancia al centro de la caja | dentro |
|---|---|---|
| Gly245 | 11,5 A | SI |
| Glu246 | 11,4 A | SI |
| His256 | 7,3 A | SI |
| Ile257 | 9,6 A | SI |
| Ser258 | 12,7 A | SI |

El centro de RRM2 esta a 9,5 A del centro de la caja. Es decir, la caja contiene el
sitio, pero **no esta centrada en el**: esta centrada en la interfaz RRM1-RRM2, que es
donde se unen los ligandos grandes (PE859, nTRD22, berberrubine).

## 2. Las poses de los fragmentos SI tocan RRM2

Midiendo contactos a menos de 4 A en el primer modelo de cada pose:

| ligando | RRM2 | RRM1 | cadena B (Phe74) |
|---|---|---|---|
| fragmento_1 | si (His256, Ile253) | si | si |
| fragmento_2 | si (His256, Ile253) | si | si |
| fragmento_3 | casi nada (6 atomos) | si | si |

El motor **encuentra el sitio**. Los fragmentos no se pierden: se colocan en el
centro de la caja, tocan RRM2 y tocan la cadena B, que es el fragmento de fenilalanina
de la otra proteina.

## 3. Y sin embargo, no es que esten penalizados

Comparados con los 10 ligandos del fondo duro que tienen **su mismo tamano** (10-16
atomos pesados):

| rango de tamano | n del fondo | afinidad media |
|---|---|---|
| **10-16** | **10** | **-5,84** |
| 17-22 | 17 | -7,14 |
| 23-28 | 54 | -7,09 |
| 29-34 | 31 | -7,55 |
| 35-45 | 29 | -7,41 |

Los fragmentos dan -5,51, -4,88 y -5,25: **0,63 kcal/mol por debajo de la media de sus
pares de tamano**, y por dentro del rango que ocupa ese fondo (de -4,42 a -5,84). Es
decir, **se puntuan exactamente como cualquier otra molecula pequena**. No hay
penalizacion especifica: sencillamente, Vina da puntuaciones debiles a los ligandos
pequenos porque puntua, en lo de aqui, sobre todo superficie enterrada.

Y el fondo de tamano grande (29-45 atomos, 60 de los 152) promedia -7,5. Ahi esta la
asimetria: los fragmentos compiten contra uniones de ARN que entierran mucha mas
superficie.

## El diagnostico

No es la caja (la cubre) ni una penalizacion especifica (se puntuan como sus pares de
tamano). Es que **el rango de puntuacion de Vina esta dominado por el tamano**, y los
fragmentos de RRM2 estan en el extremo mas pequeno de la distribucion de positivos:
11, 13 y 15 atomos pesados, frente a 25-34 de los demas.

Eso tiene una consecuencia que no es solo de TDP-43: **un fragmento pequeno con union
medida no puede ganar en este embudo, por bien colocado que este**, porque la funcion
de puntuacion le da poca cuenta. Y ahi se cae la fragmentacion como estrategia: el
proyecto tiene tres positivos que son fragmentos pequenos con union medida, y el
embudo los declara, uno por uno, como no activos.

## Que se puede hacer (y que no)

- **Se puede** repetir la puntuacion con una funcion que no premie la superficie
  enterrada (Vinardo, o el perfil de interaccion 3D que ya existe). Es CPU y esta
  medido lo suficiente para saber lo que se espera: el perfil 3D ya dio 0,566 en SOD1
  con el receptor limpio.
- **No se puede** arreglar con la caja: ya cubre el sitio, y recentrarla en RRM2 sacaria del centro los ligandos grandes, que son la mayoria de los positivos.
- **No se puede** arreglar con mas exhaustividad: no es un problema de busqueda, es de
  funcion de puntuacion.
