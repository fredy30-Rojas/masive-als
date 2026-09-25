# Cómo avanzar en MASIVE-ALS — 20 de septiembre de 2026

## El diagnóstico, en una frase

**El cuello de botella no es el cómputo ni el tamaño de la librería: es que el embudo
no ordena, y hasta hoy no estaba validado con quimias independientes.**

Todo lo demás sale de ahí. Con las auditorías de hoy sabemos que:

- **SOD1** (la diana «validada», AUC 0,815): 16 de sus 20 activos son el mismo núcleo
  pirazolona y 18 tienen otro activo a similitud ≥ 0,5. Sin la serie, el AUC queda en
  **0,601 con dos compuestos**, y **LCS-1 —el ligando de referencia real— cae en el
  puesto 152 de 219**.
- **TDP-43**: su métrica de química da **AUC 1,000 dentro de su propia familia** y
  0,538 fuera de ella. Es un detector de familia, no un predictor de unión.
- **FUS**: sin ligandos propios en la literatura. No es evaluable.

Y —esto es lo importante— **el propio paper ya lo dice a medias**: trata la puntuación
de Vina como «filtro de triaje, no como ordenador», admite el sesgo de tamaño y marca
TDP-43 y FUS como provisionales. Lo que falta es la pieza de hoy: que **el número de
SOD1 tampoco es general**, y que los conjuntos de positivos son series.

## Por qué esto importa para los objetivos reales del proyecto

No es un purismo académico. Toca las tres cosas que el proyecto necesita para avanzar:

1. **Las horas de supercomputadora.** Pedir 200.000 horas de GPU para cribar 10
   millones de compuestos con un embudo que no ordena es exactamente lo que un
   revisor rechaza. Escalar mil veces un embudo roto da mil veces más resultados sin
   orden. Un revisor del RES o de EuroHPC mira antes que nada **cómo están validados
   los controles positivos**.
2. **Los centros experimentales.** IRB, VHIR y Bellvitge solo van a probar compuestos
   con una razón creíble. Una lista ordenada por un número que no significa nada no
   les sirve; una lista con un criterio validado, sí.
3. **La credibilidad, que es el activo principal de un proyecto sin financiación y
   llevado por un paciente.** Es mejor llegar con nuestra propia auditoría que ser
   corregidos por ella. Y a los químicos medicionales se les mira primero lo que hoy
   sabemos que no estaba mirado.

## Los seis pasos, por orden de palanca

### 1. Corregir lo entregado y el paper *(sin GPU, esta semana)*

- El paquete de `entregable_investigadoras/` y el paper afirman que SOD1 está
  «VALIDADO (en uso)» con 0,815, y concluyen que la lista de SOD1 es la de mayor
  confianza. Eso hay que matizarlo con la serie congénica: **el 0,839 lo da la serie
  pirazolona; con los dos positivos independientes el número es 0,601**, y uno de
  ellos (LCS-1) no lo encuentra el acoplamiento.
- Añadir la auditoría como anexo. Un proyecto que encuentra y publica su propio sesgo
  gana credibilidad; uno al que se lo encuentran, la pierde.
- Revisar también la afirmación de que el bolsillo TDP-43 v2 está «validado»: lo que
  está validado es que **reproduce las energías publicadas** (PE859 −7,60 frente a
  −8,49), no que el ranking ordene. Ya está dicho a medias en el paper; hay que
  decirlo entero.

### 2. Arreglar la verdad de referencia *(sin GPU)*

- **Colapsar las series a un representante.** Contar 18 veces la misma química no son
  18 positivos. SOD1 pasa de 20 «activos» a ~3 independientes.
- **Buscar positivos de otras quimias, con cita y tipo de ensayo.** Para TDP-43 ya
  están localizados: rTRD01 (Kd 89 µM), nTRD22 (Kd 145 µM), los tres fragmentos de
  RRM2 y AIM4 y bis-ANS. Para SOD1, su propia literatura de ligandos fuera de la serie.
- Criterio de entrada: **sin cita y sin tipo de ensayo, el compuesto no entra**.

### 3. Rehacer las validaciones con el criterio fijado por adelantado *(CPU, horas)*

El criterio, y es uno solo: **el ranking debe colocar positivos de al menos dos
quimiotipos distintos por delante de los señuelos emparejados.** Con un solo
quimiotipo de positivos, el resultado no se reporta. Y con **dos fondos**: señuelos
emparejados en propiedades y, para TDP-43, un fondo difícil de unidores de ARN reales
(R-BIND 2.0).

Un dato que ya existe y que anima: la ronda del 20 de agosto con activos
**químicamente diversos** (rTRD01, nTRD22, bis-ANS, Congo Red…) dio AUC 0,612, mejor
que la de los nueve alcaloides (0,517). La diversidad de los positivos mejora la
señal. Con los positivos buenos y las series colapsadas, esa es la medida que vale.

### 4. Construir la métrica que no está rota: el perfil de interacción 3D *(INTENTADO el 20 sep — ver resultado)*

> **Resultado del primer intento** (`RESULTADO_PERFIL_INTERACCION_2026-09-20.md`):
> la huella **no supera a Vina** (0,65 frente a 0,815) y, con los positivos
> independientes como referencia, da **0,000**. La causa está medida y es más
> gorda que la métrica: **LCS-1 y PRG-A01 comparten 0,038 de los residuos que
> contactan**, menos que dos señuelos al azar (0,297), mientras que dentro de la
> serie pirazolona se solapan 0,397. **En esta caja no hay un modo de unión común
> entre los ligandos independientes**, así que ninguna métrica basada en las
> poses puede funcionar todavía — ni la huella ni un modelo entrenado. Antes de
> volver a tocar métricas hay que averiguar **dónde se coloca de verdad cada
> ligando**: el paso concreto está en la sección 6 de ese documento.
>
> **HECHO, y corrige lo anterior (21 sep 2026) — ver
> `INFORME_DONDE_CAEN_POSITIVOS_SOD1_2026-09-21.md`.** Ese paso se ejecutó y el
> diagnóstico de arriba **era un artefacto del receptor**: se había medido con
> `gpu_dock/SOD1.pdbqt`, que lleva las aguas del cristal y las 18 copias de 1HL5.
> Con el receptor corregido, **los siete positivos caen en Trp32** (los cuatro
> co-cristalizados incluidos), **Trp32 es el único residuo que tocan los siete**, y
> **LCS-1 contra PRG-A01 sube de 0,038 a 0,444**. La caja no hay que rehacerla. La
> huella, re-testeada con el receptor limpio, pasa de 0,000 a **0,566** (Vina 0,530):
> mejor, pero todavía por debajo, y con el vicio de tamaño intacto (+0,42).
> Lo que toca ahora es **rehacer la validación contra el receptor limpio con la
> verdad de referencia corregida** (`verdad_de_referencia.csv`), no volver a tocar
> la métrica sobre las poses viejas.

La energía está sesgada por el tamaño. La similitud 2D es una tautología de familia.
Pero el proyecto tiene algo que todavía no ha explotado: **las poses ya calculadas de
151.382 pares**. De una pose se puede leer **qué residuos se tocan y cómo** —puente de
hidrógeno donador/aceptor, apilamiento aromático, puente salino, contacto hidrófobo—,
y eso **no depende del tamaño ni del parecido 2D**.

- Construir el perfil de los positivos independientes a partir de sus poses.
- Puntuar la librería por **acuerdo de perfil**, no por energía.
- Evaluarlo con el criterio del paso 3: ¿recupera quimias distintas? ¿recupera a
  LCS-1 en SOD1?

Además queda una medida nueva y útil para cualquier diana: **la superposición de
residuos contactados entre los positivos**. Es un test de coherencia del conjunto
de positivos que se puede aplicar antes de gastar horas en métricas — si dos
ligandos que deberían unirse al mismo sitio comparten 0,038, cualquier consenso
de ellos es ruido.

### 5. Solo entonces, escalar

Con el embudo validado, la solicitud a MareNostrum/Frontier deja de ser «200.000 horas
para ordenar 10 millones de compuestos» y pasa a ser «200.000 horas para aplicar un
embudo **ya validado** contra quimias independientes a 10 millones». Eso sí se concede.

### 6. Y la decisión de la diana: ir a la agregación

Si tras el paso 3 el sitio de unión de ARN sigue sin ordenar, se deja de usar para
ordenar y el esfuerzo va al **núcleo de agregación del C-terminal** (segmentos
246-EDLIIKGISV-255 y 311-MNFGAFSINP-320), con la plantilla metodológica que funcionó
en tau (Seidler 2018): la cremallera estérica sí es un bolsillo.

## Lo que NO se hace

- **No se relanza el cribado** ni se compra cómputo hasta que el embudo pase el paso 3.
- **No se persiguen los números grandes** (−71,76, −87,35): son enterramiento.
- **No se presenta ninguna lista como candidatos** sin el criterio cumplido.
- **No se pide más supercomputadora para escalar un embudo no validado.**

## Lo que hay que decidir, y es de Fredy

1. ¿Se corrige y se hace público el sesgo ya, antes de que lo mire nadie más?
2. ¿Se avisa a Ana Martínez y Carmen Gil (CIB-CSIC), que dijeron que revisarían en
   septiembre, con la auditoría hecha?
3. ¿Seguimos con las tres dianas a la vez, o se concentra el esfuerzo en SOD1 —que es
   donde hay un ligando de referencia real, LCS-1, con el que medir de verdad—?
