# El proyecto no tiene poder estadistico, y eso cambia la pregunta — 29 de septiembre de 2026

## Que se ha medido

Todos los AUC de las dos validaciones, con su **intervalo de confianza al 95 %**
(`bootstrap_auc.py`, 5.000 re-muestreos de los positivos con reemplazo, semilla fija,
fondo fijo porque su error es despreciable al lado del de siete ligandos).

Comprobacion previa: los puntos salen **exactamente** iguales a los que lleva dias
publicando `validar_diana_limpia.py` (0,301 / 0,357 / 0,332 y 0,440 / 0,313 / 0,409).
Esa coincidencia es la garantia de que el intervalo se ha calculado sobre lo mismo que
se venia reportando. (Tuve que invertir el signo del AUC una vez: la version clasica
daba 0,699 donde el validador da 0,301, su complemento exacto. Queda documentado en el
script porque es el error mas facil de cometer al reimplementar la metrica.)

| diana | fondo | n+ | AUC | IC 95 % | P(AUC > 0,5) |
|---|---|---|---|---|---|
| TDP-43 | blando | 7 | 0,301 | **[0,089 , 0,541]** | 0,05 |
| TDP-43 | duro | 7 | 0,357 | **[0,125 , 0,604]** | 0,14 |
| TDP-43 | los dos | 7 | 0,332 | **[0,104 , 0,582]** | 0,10 |
| SOD1 | blando | 11 | 0,439 | **[0,245 , 0,657]** | 0,29 |
| SOD1 | duro | 11 | 0,313 | **[0,116 , 0,526]** | 0,05 |
| SOD1 | los dos | 11 | 0,409 | **[0,223 , 0,622]** | 0,20 |

**Las seis celdas cruzan el 0,5.** El intervalo tiene entre 0,40 y 0,48 de ancho.

## Y el otro criterio, el que DECIDE, tampoco (29 sep 2026, tarde)

El veredicto del proyecto no lo da el AUC crudo: lo da el **criterio C**, el emparejado
por tamaño, que decia "PASA" con dos de cinco químicas en TDP-43 y cinco de ocho en
SOD1. Ese criterio **nunca se le habia calculado la probabilidad de que pase con
moleculas que no unen nada**. Con pocas químicas, quedarse en la mitad buena de su grupo
de tamaño es facil por casualidad.

Se ha añadido `permutacion_ganadoras()` al validador: bajo la hipótesis nula, cada
positivo se sustituye por una molécula del **fondo** de su mismo tamaño, y se cuenta
cuántas químicas ganarían con ese mismo criterio. 2.000 permutaciones:

| diana | fondo | químicas que ganan | lo que daría el azar (mediana) | p |
|---|---|---|---|---|
| TDP-43 | blando | 2 de 5 | 2 | **0,806** |
| TDP-43 | duro | 2 de 5 | 2 | **0,715** |
| TDP-43 | los dos | 2 de 5 | 2 | **0,724** |
| SOD1 | blando | 5 de 8 | 3 | 0,180 |
| SOD1 | duro | 5 de 8 | 3 | 0,192 |
| SOD1 | los dos | 5 de 8 | 3 | 0,183 |

En TDP-43, ganar dos de cinco químicas es **exactamente lo que da el azar**: pasa en el
81 % de las permutaciones. En SOD1 es mejor que el azar (la mediana nula son 3 de 8 y se
llega a 5), pero no llega a significación.

**Ninguna de las seis celdas de las dos dianas es distinguible de la casualidad**, ni
por el AUC ni por el criterio que decide. El log del validador ahora lo dice junto al
"PASA" en todas las corridas, sin cambiar la vara que el proyecto fijó por adelantado
(moverla después de ver el resultado es la forma más rápida de convertir un proyecto en
una historia).

## Un tercer fallo encontrado de paso

`validar_diana_limpia.py --target SOD1` **perdía el fondo duro**: para SOD1, `DIANAS`
deja `poses2` a None porque el fondo lo creó `validar_sod1_limpia.py` en su propio
`POSES_FONDO2`, y el wrapper lo sobrescribía con None. La corrida se quedaba con **un
solo bloque** y no fallaba: devolvía menos y parecía que todo iba bien. Corregido, y
SOD1 vuelve a sacar sus tres bloques.

## Lo que esto significa, y no es lo que parece

El proyecto lleva semanas diciendo "el embudo NO PASA". Con estos intervalos, esa frase
queda vacia de contenido: **ninguna de las seis celdas es distinguible del azar**, ni
por encima ni por debajo. Un AUC de 0,301 con un intervalo que llega hasta 0,541 no es
"el embudo va mal", es **"no hay datos suficientes para saber si va bien o mal"**.

Eso invalida, y hay que decirlo claro, tres cosas que se han dado por buenas:

1. **La comparacion entre funciones de puntuacion.** Doce celdas con Vinardo y con Vina
   que se han comparado punto a punto: con un error de esa anchura, la diferencia entre
   0,357 y 0,344 no es una diferencia, es ruido. **Ninguna conclusion sobre la funcion
   de puntuacion de este proyecto es defendible**, ni a favor de Vina ni en contra.
   Y eso incluye la del informe de Vinardo de esta misma tarde: la conclusion de que los
   fragmentos no recuperan puesto **sigue en pie**, porque es una observacion de puesto
   (155 de 159), no una estimacion de AUC; pero la parte de "Vinardo no cambia nada"
   no se puede defender con estos numeros.

2. **La conclusion sobre los fragmentos se apoya en tres moleculas.** Que los tres
   fragmentos caigan al final es un hecho medido. La generalizacion "un fragmento
   pequeno con union medida no puede ganar en este embudo" sale de **n = 3**. Con tres
   moleculas, tres quiza, o una sola. Es la conclusion mas interesante del proyecto y
   es la menos solida. Habria que decirlo asi en cualquier escrito.

3. **El techo del embudo no se puede fijar.** Decir "el docking clasico no llega de
   0,62 en TBK1" es un numero; con 300 activos si tiene intervalo (TBK1 tiene 2.313, ahi
   si hay poder). Pero en TDP-43 y SOD1, con 7 y 11, no.

Y el que mas duele, porque es el que decidia:

4. **El "PASA" no era evidencia.** El criterio C decia que el embudo reconocia quimia
   en ambas dianas. La permutacion dice que en TDP-43 eso pasa por casualidad el 81 %
   de las veces. La conclusion ya no es "el embudo no pasa" sino "**no hay forma de
   saber si pasa** con estos datos".

## Cuantos positivos harian falta

El margen escala como 1/sqrt(n). Multiplicando por (margen actual / margen deseado)^2:

| diana | para un IC de +-0,10 | para un IC de +-0,05 | hay |
|---|---|---|---|
| TDP-43 | ~40 | ~160 | **7** |
| SOD1 | ~45 | ~185 | **11** |

Para poder afirmar con un 95 % de confianza que el embudo separa una diana del azar
hace falta del orden de **150 a 190 activos medidos**. TDP-43 tiene siete.

## El diagnostico, dicho sin rodeos

**El cuello de botella de este proyecto no es la funcion de puntuacion. Es la
cantidad de verdad de referencia.**

Se han probado, en orden: docking clasico, mas exhaustividad, mas exhaustividad otra
vez, MM-GBSA, el CNN de GNINA, y ahora Vinardo. **Tres re-puntuaciones distintas, tres
numeros que no se distinguen entre si.** Anadir una cuarta no va a cambiar nada, y ya
esta dicho y medido que anadir una quinta tampoco.

Y hay una razon peor, que es de datos y no de metodo: **TDP-43 no puede tener conjunto
de validacion en este proyecto.** No hay estructura con farmaco pequeno (4BS2 es una
interfaz proteina-proteina) y no hay serie de potencia medida de tamano utilizable. No
es que falte computo: falta el sustrato experimental. Uncribado sin serie medida no es
un problema de potencia estadistica, es un problema de existencia.

## Que haria yo, como la cientifica del proyecto

En este orden, y por una razon distinta cada uno:

**1. Dejar de producir AUC sin intervalo de confianza.** A partir de ahora, ningun
numero de este proyecto se escribe sin su IC. `bootstrap_auc.py` se ejecuta dentro de
los validadores, no aparte. Es lo mas barato de todo lo de esta lista y es lo que mas
protege al proyecto de escribir una conclusion falsa por tercera vez.

**HECHO (29 sep, tarde).** `intervalo_auc()` y `permutacion_ganadoras()` ya estan
dentro de `validar_sod1_limpia.py`, asi que las dos medidas salen en el log de TODA
corrona de las dos dianas y en el diccionario de resultados. Comprobado que el punto
del AUC no se mueve: 0,301 / 0,357 / 0,332 en TDP-43 y 0,439 / 0,313 / 0,409 en SOD1,
identicos a los de antes.

**2. No volver a tocar la funcion de puntuacion.** Es la decision que mas horas de GPU
ahorra. El problema ya no es donde poner los numeros: es que no hay con que compararlos.

**3. Construir un conjunto de validacion donde exista la verdad de referencia.** Esta es
la via que si convierte el proyecto de infalsable en falsable. Concretamente: buscar en
ChEMBL una diana con **cincuenta o mas activos medidos** (IC50 o Kd) **y** al menos un
ligando co-cristalizado con estructura. Para esa diana, el control de redocking ya esta
resuelto —`andr` sirve de modelo — y los senuelos se generan con la receta que ya esta
validada con los bancos limpios de CDK2 y TBK1. Con cincuenta positivos el IC baja a
+-0,10 y las preguntas sobre la funcion **empiezan a tener respuesta**.

**4. El tercer fondo pasa de "pendiente" a prioridad.** Sin compuestos medidos que NO
unen, ningun AUC es interpretable: se esta midiendo enrichment contra una lista que solo
contiene binders, y eso es un numero sin contra. Es lo unico que necesita el laboratorio
y lo unico que no se puede resolver con CPU de sobra. **Es la peticion que hay que
llevar al laboratorio.**

**5. Los fragmentos no se resuelven con mas cribado.** Si el proyecto quiere
fragmentos, el experimento es un experimento de fragmento: soaking o SPR sobre los tres
cristalizados, que ya sabemos que unen, y ver si amplian el sitio. Comprar mas
quinolinas y meterlas en el embudo no va a contestar la pregunta.

## Lo unico que hay que decidir

Si el laboratorio puede gastar en ~cincuenta activos medidos de una diana con serie
conocida. Si puede, el proyecto tiene un futuro medible y las tres funciones de
puntuacion que ya se han probado empiezan a poder compararse. Si no puede, lo honesto es
escribir el resultado negativo **con el analisis de poder** —que es un trabajo real y
publicable, y mucho mas defendible que un AUC de 0,30 sin intervalo— y dejar el embudo
como herramienta de priorizacion, que es para lo que los datos actuales muestran que
sirve.

## Ficheros

- `analysis/bootstrap_auc.py` — el calculo, con la comprobacion contra el validador.
- `analysis/bootstrap_auc.log` — la salida con los seis intervalos y la cuenta de cuantos positivos harian falta.
