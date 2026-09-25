# Lista corta de la libreria TDP-43: las tres formas de ordenarla, y por que ninguna sirve sola

25 de septiembre de 2026

## Que se ha hecho

Los resultados de la libreria contra TDP43_v2 ya estaban acoplados desde el 18 de
septiembre: 108.841 poses, y no faltaba ninguna por acoplar (las 4.471 que
aparentemente faltaban son exactamente los ligandos que el motor rechaza por
tipos de atomo invalidos, que estan en `ligandos_excluidos.txt`). Lo que no estaba
hecho era leerlos con la medida que la regla de decision del 24 de septiembre
declara valida. Eso es lo que se ha hecho hoy, sin tocar la tarjeta grafica: solo
leer ficheros que ya existian.

Se han sacado tres ordenaciones de los mismos 108.833 compuestos (los 8 restantes
no tienen fichero en `libreria_ligands`), para poder compararlas:

1. Energia dividida por atomos pesados — `resultados_libreria/ranking_TDP43_v2_por_atomo.csv`
2. Energia menos la recta del fondo duro — `resultados_libreria/residual_TDP43_v2.csv`
3. Emparejando por numero de atomos — `resultados_libreria/emparejado_TDP43_v2_15_50atomos.csv`

Mediana de tamaño de la libreria: 27 atomos pesados, de 1 a 112. Los positivos
medidos de TDP-43 tienen 21,6 de media (nTRD22 29, berberrubine 24, rTRD01 25,
PE859 34, y tres fragmentos de 11, 13 y 15).

## 1. Por atomo pesado: la cabeza son fragmentos

El primer puesto es CHEMBL562972, -8,40 kcal/mol con 15 atomos (-0,560 por atomo).
Barriendo la cabeza de la lista, todos los puestos hasta el 40 tienen entre 15 y 17
atomos. No es quimica: es aritmetica. La afinidad crece menos que proporcionalmente
con el tamaño, asi que al dividir siempre gana el compuesto mas pequeño que se deje
entrar en la ventana. Mover la ventana mueve el resultado, y eso es la definicion de
artefacto.

Esta es, sin embargo, la medida que manda la regla de decision del 24 de septiembre,
y en el control da AUC 0,742 contra el fondo duro. Las dos cosas son verdad a la vez,
y la explicacion esta en el punto 3.

## 2. Por residual: sale casi lo mismo que la energia cruda

La recta ajustada sobre los 152 unidores de ARN de R-BIND 2.0, acoplados contra este
mismo receptor y esta misma caja, es:

    afinidad = -6,686 - 0,0135 * (atomos pesados)

Un carbono pesado de mas vale una centesima de kcal/mol. Con una pendiente tan plana,
restar la recta y no restar nada es casi lo mismo, de modo que esta lista es casi la
lista por energia cruda. Y la energia cruda, el propio validador la hunde: AUC 0,352
contra el fondo duro, por debajo del azar, y lo dice con todas las letras — "la lista
ordenada por energia SIN mas no sirve para elegir candidatos: su cabeza es grande, no
buena".

Cabeza: CHEMBL4581197 (-10,90 con 36 atomos), CHEMBL21197 (-10,80 con 35),
CHEMBL4553709 (-10,80 con 39). Es la misma cabeza que ya salia en `run_libreria.log`
por energia cruda, lo cual confirma el diagnostico.

## 3. Emparejando por tamaño: el criterio que menos trampa tiene, y el que da la mala noticia

Aqui cada compuesto se compara solo contra los de su mismo numero de atomos pesados,
que es la parte C de la validacion, la que decidio el PASA. Cabeza:

    puesto 1   CHEMBL41430    -9,90   20 atomos   z -3,96
    puesto 2   CHEMBL515619  -10,00   21 atomos   z -3,93
    puesto 3   CHEMBL25201   -10,10   26 atomos   z -3,88
    puesto 4   CHEMBL4581197 -10,90   36 atomos   z -3,87
    puesto 5   CHEMBL27251    -9,80   20 atomos   z -3,84

En la libreria, el 1% mejor tiene z -2,33 y el 0,1% mejor z -3,15. La cabeza esta
por debajo de -3,8: son compuestos que destacan muchisimo entre sus iguales.

La mala noticia es donde caen los positivos conocidos en esa misma lista:

    nTRD22        29 atomos   -7,90   z -1,04   le ganan 16.223 compuestos
    berberrubine  24 atomos   -7,80   z -1,09   le ganan 14.969 compuestos
    rTRD01        25 atomos   -7,40   z -0,59   le ganan 29.568 compuestos
    PE859         34 atomos   -7,20   z -0,05   le ganan 52.345 compuestos

Es decir: **la berberrubina y el nTRD22, que son unidores medidos de TDP-43, caen en
el 14% de cabeza de la libreria.** Cualquier corte que deje menos de quince mil
compuestos deja fuera al mejor unidor conocido del proyecto. Un enriquecimiento de
siete veces sobre el azar suena a algo, pero sobre 108.833 compuestos significa que
la lista corta no es corta: para no perder al bueno hay que arrastrar quince mil.

## La conclusion, dicha sin adornos

El cribado de la libreria esta bien hecho y sirve para lo que sirve: saber que esos
108.841 compuestos se han acoplado y con que energia. Lo que no sirve es como
selector de candidatos. Las tres ordenaciones dicen lo mismo por caminos distintos:
la energia de acoplamiento contra este receptor ordena sobre todo el tamaño, y una
vez que se quita el tamaño (punto 3), los unidores conocidos solo suben al 14%.

Esto no es un fallo de hoy ni un fallo del montaje. Es lo que ya decian los propios
informes del proyecto: el del 19 de septiembre avisaba de que el AUC 0,92 de la
similitud era una tautologia, el del 22 de septiembre de que el MM-GBSA tiene el
sesgo de tamaño mas fuerte de las cuatro funciones probadas, la auditoria del 11 de
septiembre de que el rescoring no era validable. Hoy se anade el penultimo eslabon
que quedaba: tampoco el orden por acoplamiento aguanta un examen por tamaño.

## Lo que si se puede hacer con esto, sin mentirse

- Usar el acoplamiento como **prefiltro grueso y nada mas**: quedarse con el 15% de
  cabeza por z emparejado (unas 15.000 moleculas), que es el unico corte que no
  tira a los positivos conocidos. A partir de ahi, filtrar por otra cosa que no sea
  la energia de acoplamiento.
- **Vinardo** como segunda vuelta, que es lo que recomienda el informe del 22 de
  septiembre frente al MM-GBSA, porque conserva mucho menos sesgo de tamaño.
- La **puntuacion de interacciones** que ya existe (`analysis/puntuacion_interaccion.py`,
  `combinar_puntuaciones.py`), que mira contactos y no solo energia.
- El MM-GBSA, si se usa, solo dentro de una ventana estrecha de tamaño y sabiendo
  que su relacion señal/ruido medida era de 1,5 y que dos corridas identicas se
  llevaban 2 kcal/mol.
- Y por encima de todo: **la prueba experimental**, que es lo unico que ningun
  calculo de esta casa puede sustituir.

## Aviso sobre el fichero de anotaciones

Los nombres de los compuestos con mas peso aqui son identificadores de ChEMBL
(CHEMBL4581197, CHEMBL21197, ...). Uno de los que aparece en la lista por residual,
en el puesto 29, tiene nombre comun: CAPMATINIB, un inhibidor de c-Met aprobado.
Aparecer en esta lista no es ninguna evidencia de que se una a TDP-43; el capmatinib
sale por ser grande y puntuar bien, que es justo el sesgo que se ha medido.

## Ficheros que deja este trabajo

- `gpu_dock/rankear_libreria_por_atomo.py` — la lista por atomo pesado
- `gpu_dock/residual_libreria.py` — la lista por residual contra R-BIND
- `gpu_dock/rankear_por_tamaño.py` — la lista emparejando por tamaño
- `gpu_dock/resultados_libreria/_pesados_libreria.json` — atomos pesados de los 113.304
  ligandos, contados una vez para los cuatro destinos
- `gpu_dock/resultados_libreria/ranking_TDP43_v2_por_atomo.csv` (108.833 filas)
- `gpu_dock/resultados_libreria/residual_TDP43_v2.csv` (108.833 filas)
- `gpu_dock/resultados_libreria/emparejado_TDP43_v2_15_50atomos.csv` (102.378 filas)
