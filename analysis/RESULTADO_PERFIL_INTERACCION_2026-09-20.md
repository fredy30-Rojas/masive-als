# La huella de interacción 3D: primer intento, y el diagnóstico que sacó

> **CORREGIDO EL 21 DE SEPTIEMBRE DE 2026. Leer antes de usar este documento.**
>
> El diagnóstico de las secciones 4 y 5 («los ligandos independientes no ocupan el
> mismo sitio, no hay modo común») **era un artefacto del receptor, y queda
> retirado**. Esta medida se hizo con `gpu_dock/SOD1.pdbqt`, que no es una proteína
> limpia: lleva las 1.763 aguas del cristal y las 18 copias de 1HL5. Acoplando los
> siete positivos contra el receptor corregido (`gpu_dock/SOD1_limpio.pdbqt`),
> **los siete caen en Trp32** y **LCS-1 contra PRG-A01 pasa de 0,038 a 0,444** de
> solapamiento de residuos.
>
> Lo que **sigue en pie**: la huella no bate a Vina, y re-testeada con el receptor
> limpio da AUC 0,566 (antes 0,000) frente a 0,530 de Vina — mejor, pero insuficiente,
> y con el mismo vicio de tamaño (+0,42 con el número de contactos).
>
> Detalle: `INFORME_DONDE_CAEN_POSITIVOS_SOD1_2026-09-21.md`.

**20 de septiembre de 2026.** Script: `analysis/perfil_interaccion.py`.
Datos: las 219 poses ya calculadas de la validación de SOD1 (20 activos + 199
señuelos, `_validacion_SOD1/out`). No se acopló nada nuevo.

---

## 1. La idea, y por qué parecía la buena

Las dos medidas que el proyecto tenía están rotas, y por motivos distintos: la
**energía** está dominada por el tamaño (−0,59 a −1,26 kcal/mol por átomo pesado)
y la **similitud 2D** es un detector de familia (AUC 1,000 dentro de su familia,
0,538 fuera). Ninguna de las dos mira lo que de verdad describe una unión: **dónde
y cómo se apoya la molécula en la proteína.**

Eso sí está en las poses, que ya están calculadas. Así que se construyó la huella
de interacción (IFP): para cada pose, la lista de **átomos del receptor que el
ligando toca**, con su residuo, su nombre de átomo y los tipos de los dos
extremos; y se puntúa cada molécula por su acuerdo con el consenso de los activos.

## 2. Resultado: no supera a Vina

| Medida | AUC (20 activos vs 199 señuelos) |
|---|---|
| Vina (referencia) | **0,815** † |
| Huella, *jaccard* | 0,650 |
| Huella, *recall* | 0,609 |
| Huella, *coseno* | 0,465 |

† Y ese 0,815 ya sabemos que lo sostiene una serie congénérica
(`AUDITORIA_POSITIVOS_2026-09-20.md`).

**El primer intento no mejora la energía.** Se reporta porque es lo que hay.

## 3. Y al probarlo con la referencia limpia, se cayó: AUC 0,000

Se repitió definiendo el consenso **sólo con los dos activos de quimias
independientes** (LCS-1 y PRG-A01). Resultado: **AUC 0,000**, con los dos en los
puestos 218 y 219 de 219. No es un error de código; es un hallazgo, y el número
apunta a algo más gordo que la métrica.

## 4. La causa, medida: los ligandos independientes no ocupan el mismo sitio

Se midió la superposición de **residuos contactados** entre pares de ligandos
(conjunto de residuos del receptor a los que cada pose se acerca):

| Comparación | Jaccard de residuos |
|---|---|
| Dentro de la serie pirazolona (18 compuestos) | **0,397** |
| Dentro de los señuelos, al azar | 0,297 |
| Serie contra señuelos | 0,243 |
| **LCS-1 contra PRG-A01 (los dos independientes)** | **0,038** |

Los dieciocho pirazolonas van todos al mismo sitio y se solapan entre ellos casi
el doble que dos señuelos al azar. **Los dos activos independientes no comparten
prácticamente ningún residuo entre ellos** — menos incluso que dos señuelos
al azar.

Es decir: **en esta caja no hay un modo de unión común entre los ligandos
independientes.** O no se unen ahí, o el acoplamiento no los coloca ahí. Y la
segunda es coherente con lo que ya sabíamos: LCS-1, el ligando de SOD1 con más
respaldo independiente, cae en el puesto 152 de 219.

## 5. Qué significa esto

1. **El 0,815 no medía lo que parecía.** Mide que la caja encuentra la serie
   pirazolona, que ocupa un sitio concreto y consistente. No mide que la caja
   encuentre ligandos de SOD1 en general, porque los ligandos independientes **no
   reproducen ese modo**.
2. **Ninguna métrica basada en las poses puede funcionar aquí todavía**, ni la
   energía, ni la huella, ni un modelo entrenado: si los positivos independientes
   no comparten colocación, no hay señal común que aprender. El problema es
   anterior a la métrica.
3. **Lo que hay que averiguar antes de tocar ninguna métrica:** dónde se colocan
   de verdad los ligandos independientes. Es un trabajo de inspección de poses y
   de comparación con los sitios documentados, y no necesita GPU.

## 6. El siguiente paso concreto, y es barato

1. **Dockear los siete positivos independientes** de `activos_sod1_v2.csv`
   (LCS-1, PRG-A01, el representante de la pirazolona y los cuatro ligandos
   co-cristalizados de Trp32: 5-fluorouridina, isoproterenol, dopamina y
   epinefrina). Los cuatro cristalográficos son el mejor patrón de medida que se
   puede tener: su sitio está resuelto por cristalografía, no inferido.
2. **Mirar dónde cae cada uno** y compararlo con el sitio documentado. Si los
   cristalográficos caen en Trp32 y LCS-1 no, el problema es de colocación de
   LCS-1; si no cae ninguno, la caja Trp32 no es el sitio de estos ligandos y hay
   que replantearlo — con el mismo criterio que llevó a corregir la caja de
   TDP-43 cuando se descubrió que estaba a 17 Å del sitio real.
3. **Sólo después** tiene sentido volver a intentar una métrica de colocación.
4. Y en TDP-43, lo mismo: la huella se probará allí cuando los positivos estén
   preparados y colocados, porque allí la caja sí reproduce las energías
   publicadas.

## 7. Lo que este intento deja dicho, y sirve

- La huella de interacción **no es un atajo** que arregle esto sin tocar la
  verdad de referencia. Se descarta como atajo, con números.
- El motivo por el que falla **no es la métrica: es que no hay un modo común**.
  Saber eso ahorra semanas de trabajo sobre una base que no lo permitía.
- Y queda una medida nueva y útil para el proyecto: **la superposición de
  residuos contactados** entre positivos. Si dos ligandos que deberían unirse al
  mismo sitio comparten 0,038, cualquier consenso de ellos es ruido. Es un test
  de coherencia del conjunto de positivos que se puede aplicar a cualquier diana
  antes de gastar horas en métricas.
