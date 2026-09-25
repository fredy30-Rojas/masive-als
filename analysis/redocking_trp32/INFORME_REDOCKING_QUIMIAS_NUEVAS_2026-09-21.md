# Las siete quimias nuevas NO pasan el control de redocking, y el único positivo del proyecto no se sostiene

**21 de septiembre de 2026.** Scripts: `redock_trp32.py` (ampliado),
`analizar_modos_quimias.py` (nuevo), `anomalia_quinazolinas.py` (nuevo),
`pose_nativa_en_cribado.py` (nuevo). Registros: `redocking_quimias_nuevas.log`,
`modos_quimias_nuevas.log`, `anomalia_quinazolinas.log`,
`pose_nativa_en_cribado.log`. Datos: `resultado_redocking_quimias_nuevas.json`,
`modos_quimias_nuevas.csv`, `anomalia_quinazolinas.csv`,
`pose_nativa_en_cribado.csv`.
Antecedente directo: `../INFORME_VALIDACION_SOD1_V4_2026-09-20.md`.

---

## 0. Resumen, en tres frases

1. **Ninguna de las siete quimias nuevas reproduce su pose cristalina.** La pose
   mejor puntuada queda entre 4,3 y 10,5 Å de la del cristal; en 27 poses
   generadas por estructura, solo una (la ZO0, a 1,43 Å) aparece cerca del
   cristal, y aun así la función prefiere otra a 5,8 Å.
2. **No es un problema del buscador, es de la función de puntuación.** La pose
   cristalina es un mínimo local de verdad (al minimizarla se queda quieta, deriva
   0,22–0,4 Å) y **aun así** la función la castiga entre **1,0 y 2,8 kcal/mol**
   frente a una colocación distinta. Con el ligando **rígido** (solo traslación y
   giro, cero torsiones) tampoco se recupera: 3,9–11,9 Å.
3. **El único resultado positivo del proyecto (el ranking de quinazolinas de la
   v4) no sobrevive.** Los −10,0 y −9,8 kcal/mol que ponían a 4MQ y ZO0 en el 1 %
   superior no son la pose del cristal. En el receptor de sus propias estructuras,
   4MQ da −5,6 y 12I −5,9: **la diferencia de 5 kcal/mol por un CH₂ solo existe
   dentro del receptor del cribado**, y las dos poses buenas se van a meter en una
   grieta entre dos copias de la proteína del modelo (cadenas A y J de 1HL5).

---

## 1. Qué se pidió y por qué

El informe de la v4 dejaba como primer paso obligado el que faltaba desde el
principio: **«Redocking de las 7 quimias nuevas en su propia estructura
(RMSD ≤ 2 Å) antes de usar ninguna como patrón»**. Tenía razón de ser: el barrido
del PDB había encontrado siete quimias nuevas en Trp32 (4MQ, 12I, ZO0, ZZT, K4I,
6B3, 946) y eran justo las que sostenían el único «positivo» del proyecto — dos
quimias independientes por delante del fondo, en el 1 % superior. Ninguna se había
comparado todavía con su pose depositada.

## 2. Método

Se amplió `redock_trp32.py` (sin tocar la ronda 1):
`SISTEMAS` ahora incluye las siete entradas nuevas y, cuando la estructura no es
el mutante I113T de Wright (8GSQ, 6A9O y 5YTO son otras construcciones y otra
numeración), la cadena se elige sola con `elegir_copia_trp32`: la copia del
ligando que está más cerca del **anillo indol de Trp32** de su misma cadena. Las
distancias que devuelve (3,41 / 3,29 / 3,18 / 3,33 / 3,33 / 3,39 / 3,60 Å)
coinciden con las del barrido original, o sea que la copia elegida es la buena.

Cada control usa **el receptor de su propia estructura** (con el ligando nativo
fuera y todos los átomos donde el depósito dice), caja de 24 Å centrada en el
centroide del ligando cristalino, exhaustividad 8 y **tres semillas fijas**
(42, 2026, 777), con el medidor de RMSD corregido el 20 de septiembre. El
acoplamiento repite la fórmula del cribado, así que lo que se mide es el protocolo
del proyecto, no otro.

### Tres defectos del camino, encontrados y corregidos

1. **Meeko no preparaba 3 de las 7 estructuras** (2WZ6, 8GSQ, 6A9O, 5YTO). Dos
   causas distintas: ambigüedad de conformaciones alternas —en 2WZ6 hay residuos
   que existen con altloc B y no con A— y residuos que no emparejan con ninguna
   plantilla (restos terminales «residuo 0» de 10 cadenas, una serina incompleta
   en 5YTO C:142SER). Se arregla en `extraer_receptor`, que ahora decide átomo por
   átomo: si el átomo tiene un registro sin altloc, ese manda; si solo existe con
   altloc, el de **mayor ocupación**. Es más fiel que `--default_altloc A`, que
   reconstruye cadenas laterales del bolsillo en vez de conservar las coordenadas
   del depósito.
2. **2WZ6 tenía además un residuo que rompía la cadena peptídica** («Expected 4
   paddings for (A:132, A:133)»). El par Glu132–Glu133 de las dos cadenas está a
   media ocupación (0,50 y 0,30) y meeko detecta dos enlaces entre ellos. Está a
   **40 Å de la caja** y no se ve desde el bolsillo, así que se descarta
   (`RESIDUOS_DESCARTADOS`) y queda documentado.
3. **El ligando de 6A9O está a medio resolver**: el depósito solo modela **24 de
   los 33 átomos pesados** del Lig9 (falta toda la cola). No se puede comparar
   átomo a átomo con el ligando ideal, así que su RMSD se mide **solo sobre los
   átomos resueltos**, emparejando el cristal como subconjunto de la pose
   (`rmsd_sobre_resueltos`), y así se reporta.

**Comprobación de que el arreglo no rompió lo que funcionaba:** se repitió el
control de la ronda 1 con el script nuevo. 4A7T (isoproterenol) sigue **PASANDO**:
0,48 / 2,50 / 0,48 Å por semilla, mediana 1,28 Å, 2 de 3 semillas por debajo de
2 Å. Igual que lo publicado.

## 3. Resultado: el control de redocking de las siete quimias nuevas

| Entrada | Ligando | Química | RMSD de la pose **mejor puntuada** | Mejor de las 27 poses | Poses ≤ 2 Å | Veredicto |
|---|---|---|---|---|---|---|
| 4A7Q | 4MQ | diazepanoquinazolina | 10,53 Å | 3,15 Å | 0 | **FALLA** |
| 4A7G | 12I | piperazinaquinazolina | 5,23 Å | 2,70 Å | 0 | **FALLA** |
| 2WZ6 | ZO0 | CF₃-diazepanoquinazolina | 5,77 Å | 1,43 Å | 1 | **FALLA** (recupera en el modo 9, no en el 1) |
| 2WZ0 | ZZT | anilina | 4,29 Å | 3,60 Å | 0 | **FALLA** |
| 8GSQ | K4I | paliperidona | 6,84 Å | 3,25 Å | 0 | **FALLA** |
| 6A9O | 6B3 | Lig9 (fenantridinona) | 5,03 Å\* | 4,75 Å\* | 0 | **FALLA** |
| 5YTO | 946 | aminoalcohol naftalénico | 4,86 Å | 2,43 Å | 0 | **FALLA** |

\* solo sobre los 24 átomos que el cristal tiene resueltos.

Y los cuatro controles de la ronda 1, medidos con el mismo script para tener la
comparación en la misma tabla:

| Entrada | Ligando | Mejor de las 27 poses | Poses ≤ 2 Å | Veredicto |
|---|---|---|---|---|
| 4A7T | isoproterenol | **0,48 Å** | 6 | **RECUPERA** |
| 4A7U | adrenalina | **0,69 Å** | 4 | **RECUPERA** |
| 4A7V | dopamina | **0,74 Å** | 4 | **RECUPERA** |
| 4A7S | 5-fluorouridina | 9,56 Å | 0 | no recupera |

O sea: el protocolo **sí** reproduce las catecolaminas y **no** reproduce ninguna
de las siete quimias nuevas. El patrón es limpio y no es ruido.

## 4. ¿Falla el buscador o falla la puntuación? Falla la puntuación

Esta es la pregunta que separa «no lo encuentra» de «lo encuentra y lo tira», y
se responde con tres medidas independientes.

**(a) Puntuar la pose del cristal sin moverla (`--score_only`).**

| Entrada | E(pose cristal, quieta) | E(pose cristal minimizada) | Deriva al minimizar | E(mejor pose acoplada) | ΔE |
|---|---|---|---|---|---|
| 4A7Q | −4,15 | −4,26 | 0,22 Å | −5,59 | +1,44 |
| 4A7G | −3,85 | −4,24 | 0,30 Å | −5,85 | +2,01 |
| 2WZ6 | −4,68 | −5,18 | 0,37 Å | −5,93 | +1,25 |
| 2WZ0 | −1,77 | −2,62 | 1,19 Å | −4,47 | +2,70 |
| 8GSQ | −6,43 | −7,12 | 0,28 Å | −9,93 | +3,50 |
| 5YTO | −6,04 | −6,86 | 0,33 Å | −8,95 | +2,91 |

La pose cristalina **se queda donde está** al minimizarla (0,2–0,4 Å en cinco de
seis casos): es un mínimo local de verdad. Y aun así la función la coloca entre
**1,25 y 3,50 kcal/mol por debajo** de una colocación que no es la del cristal.

**(b) Acoplar el ligando RÍGIDO, en su conformación cristalina.** Sin torsiones
que muestrear, solo traslación y giro: si el problema fuera de muestreo, aquí la
pose tendría que volver sola a su sitio.

| Entrada | E(cristal quieta) | E(mejor rígida) | ΔE | Deriva al minimizar | Mejor RMSD | Veredicto |
|---|---|---|---|---|---|---|
| 4A7Q | −4,40 | −5,59 | +1,20 | 0,11 Å | 10,85 Å | FALLA |
| 4A7G | −4,07 | −5,64 | +1,57 | 0,24 Å | 10,85 Å | FALLA |
| 2WZ6 | −5,23 | −6,28 | +1,05 | 0,34 Å | 5,46 Å | FALLA |
| 2WZ0 | −1,92 | −4,76 | +2,83 | 1,41 Å | 4,34 Å | FALLA |
| 8GSQ | −8,12 | −9,77 | +1,65 | 0,48 Å | 11,86 Å | FALLA |
| 5YTO | −8,34 | −9,88 | +1,54 | 0,23 Å | 3,87 Å | FALLA |

Con el ligando congelado en la forma exacta del cristal, la función **sigue
prefiriendo otra colocación** y la manda a 3,9–11,9 Å. No hay nada que muestrear:
el que decide mal es el que puntúa.

**(c) Las 27 poses.** Solo la ZO0 tiene una pose a menos de 2 Å (1,43 Å, modo 9),
y la función la descarta en favor de una a 5,77 Å.

**Conclusión:** para estas siete quimias el protocolo no está validado, y la causa
está en el término de puntuación, no en la caja, ni en la exhaustividad, ni en la
preparación del ligando. Es la misma clase de fallo que ya se midió con la
5-fluorouridina, pero aquí es sistemático: **siete de siete**.

## 5. La anomalía de las quinazolinas: de dónde venían los −10 kcal/mol

La v4 se quedó en «cinco kcal/mol por un CH₂ no es físico». Ya se sabe de dónde
salen.

### (a) No es ruido: se reproduce

La validación v4 acopló **sin semilla fija**, o sea que su −10,02 era de una sola
corrida estocástica. Repetido con tres semillas fijas:

| Ligando | v4 (sin semilla) | semilla 42 | semilla 2026 | semilla 777 | Dispersión |
|---|---|---|---|---|---|
| 4MQ (diazepano) | −10,02 | −10,07 | −9,99 | −9,70 | 0,38 |
| ZO0 (diazepano + CF₃) | −9,77 | −9,83 | −9,76 | −9,83 | 0,07 |
| 12I (piperazina) | −4,96 | −4,89 | −4,91 | −4,91 | 0,02 |

Se reproduce. El número era real **en ese receptor**.

### (b) No son las aguas cristalográficas

Aquí apareció un defecto que nadie había mirado: **`SOD1.pdbqt` no es una
proteína limpia**. Lleva dentro las **1.763 aguas cristalográficas** de 1HL5
(tipadas como OA, o sea que Vina las trata como parte rígida del receptor y **no
las puede desplazar**), los iones de zinc, y **las 18 copias de la proteína** del
fichero depositado (cadenas A–Q, S; 21.585 átomos). Cinco de esas aguas caen
dentro del bolsillo de Trp32, a menos de 10 Å del centro de la caja. Un ligando
que se aprieta entre la proteína y aguas congeladas cobra términos que en la
realidad no existirían.

Se repitió todo con el mismo receptor **sin las aguas** (19.822 átomos):

| Ligando | Con aguas | Sin aguas | Cambio |
|---|---|---|---|
| 4MQ | −10,07 | **−10,45** | −0,38 |
| ZO0 | −9,83 | **−10,89** | −1,06 |
| 12I | −4,89 | **−5,52** | −0,63 |

Sin aguas hasta puntúan algo mejor. **Las aguas no son la causa**, pero el defecto
del receptor se queda anotado: hay que decidir si se rehace el receptor sin
disolvente y con una sola copia del dímero, y si eso obliga a repetir el cribado
SOD1 entero.

### (c) La causa: el receptor del cribado no es el bolsillo donde se une el ligando

Se lleva la **pose cristalina** de cada quinazolina al receptor del cribado
(superponiendo el bolsillo de Trp32 de su cadena sobre el de la cadena A de 1HL5
con 140 átomos; desvío de la superposición 0,75–0,86 Å) y se puntúa:

| Ligando | E(pose nativa) en el receptor del cribado | E(pose nativa minimizada) | E(mejor pose acoplada) |
|---|---|---|---|
| 4MQ | +0,75 | −0,09 (deriva 9,95 Å) | −10,07 |
| ZO0 | +7,03 | −3,88 (deriva 4,93 Å) | −9,83 |
| 12I | +7,66 | −1,47 (deriva 6,31 Å) | −4,89 |

La pose nativa **no cabe**: llega con choques (0,75 a 7,66 kcal/mol positivos) y al
minimizarla se escapa 5–10 Å. Hay que ser honesto con la limitación: el desvío de
0,75–0,86 Å del bolsillo entre dos formas cristalinas distintas hace que una pose
trasplantada empiece en contacto malo, así que este test **no prueba** que la pose
nativa valga −5 en ese receptor; solo prueba que **no es la pose que Vina elige**.

Lo que sí es concluyente, y no depende de ninguna superposición:

- **En el receptor de su propia estructura**, 4MQ da −5,59 y 12I da −5,85. Casi
  idénticos. El CH₂ no vale 5 kcal/mol.
- **Las dos poses «de −10» van a parar a una grieta entre dos copias de la
  proteína del modelo**: los contactos a menos de 4 Å se reparten entre la cadena
  **A (23 pares) y la cadena J (18 pares)**. La interfaz A–J de 1HL5 es un
  contacto débil de empaquetamiento cristalino (48 pares de átomos a 4 Å, frente a
  los cientos de un dímero real), no el sitio de unión de nadie. La caja de 22 Å
  centrada en Trp32 (que está bien puesta: Lys30 a 1,1 Å, Glu21 a 3,3, Val29 a
  3,8, Trp32 a 4,1) se extiende hasta esa grieta, y allí el anillo de diazepano se
  entierra a gusto.

En resumen: **el −10 no describe una unión real; describe dónde cabe mejor una
molécula en un modelo con 18 copias de la proteína y una grieta abierta entre
ellas.**

## 6. Qué cambia en el proyecto

1. **El criterio del 20 de septiembre hay que ampliarlo con un paso previo.** El
   criterio actual («un ranking solo cuenta si coloca positivos de al menos dos
   quimias distintas por delante de los señuelos, y sin fondo difícil no se
   reporta») da por supuesto que el protocolo reproduce la pose de los positivos
   que usa. Ya sabemos que eso no es cierto para estas quimias. La regla debería
   ser: **ninguna quimia puede entrar en la lista de positivos de un ranking sin
   haber pasado antes su propio control de redocking** (RMSD ≤ 2 Å de la pose
   mejor puntuada, no de «alguna» pose).
2. **El resultado de la v4 no se puede usar.** El conjunto de controles quedó
   limpio y con 10 quimias independientes, y eso sigue siendo un avance real del
   trabajo — pero la señal que produjo (AUC 0,69 de las cristalizadas, AUC 0,72
   con el fondo emparejado nuevo) depende de dos compuestos cuyas poses son
   artefactos del receptor. Los números hay que corregirlos a la baja y explicarlo.
3. **SOD1 tiene ahora tres problemas abiertos propios:** el receptor lleva aguas
   y 18 copias; la 5-fluorouridina no se recupera (10 de 9 poses); y ninguna de
   las siete quimias nuevas se recupera. Con eso, **el cribado de SOD1 con esa
   caja no es interpretable** hasta rehacer el receptor y volver a medir.
4. **Lo que sí queda en pie:** las catecolaminas se recuperan (0,48 / 0,69 /
   0,74 Å) con ligando flexible y con la caja de 24 Å. El protocolo **está validado
   para ese quimiotipo**, y solo para ese. La dopamina, además, sigue siendo un
   fallo de puntuación puro (la pose buena está en el modo 2).

## 7. Qué hacer a continuación, en orden

1. **Rehacer el receptor SOD1**: una sola unidad biológica (el dímero, no las 18
   copias), sin aguas, con los iones de Zn/Cu en su sitio o explícitamente fuera, y
   documentar la decisión. Volver a medir el control de redocking de las
   catecolaminas contra ese receptor antes de tocar nada más.
2. **Repetir el control de redocking de las siete quimias nuevas contra el
   receptor rehecho**, ya sin grieta artificial, para saber si el fallo de
   puntuación era del bolsillo o del modelo.
3. **Solo después**, y si el control pasa, volver a mirar el ranking. Si no pasa,
   la conclusión honesta es que el cribado de SOD1 con Vina en esa caja no es
   interpretable y hay que cambiar de metodología (rescoring MM-GBSA focalizado,
   como ya se hizo con TDP-43, o un protocolo con receptor flexible).
4. **TDP-43 y FUS** siguen con sus propios pendientes y no se tocan en este
   informe.

---

## Anexo: archivos de esta ronda

| Archivo | Qué es |
|---|---|
| `redock_trp32.py` | ampliado: 7 entradas nuevas, selección de copia por Trp32, resolución de altlocs, `-a` en meeko, `--rehacer-receptor`, `--salida` |
| `analizar_modos_quimias.py` | mide las 27 poses de cada estructura y separa muestreo de puntuación |
| `anomalia_quinazolinas.py` | reproducibilidad con semillas fijas, receptor sin aguas, perfiles de contacto |
| `pose_nativa_en_cribado.py` | trasplanta la pose cristalina al receptor del cribado (Kabsch sobre el bolsillo) y la puntúa |
| `resultado_redocking_quimias_nuevas.json` | resultado por semilla de las 7 entradas |
| `modos_quimias_nuevas.csv` | las 27 poses × 11 entradas con su RMSD |
| `anomalia_quinazolinas.csv` | afinidades y contactos por ligando y variante |
| `pose_nativa_en_cribado.csv` | pose nativa, minimizada y acoplada, por ligando |
| `redocking_quimias_nuevas.log` | traza completa de la ronda 2 |
| `out/anomalia_*_s*.pdbqt` | poses de la prueba de anomalía |
| `receptores/SOD1_sin_agua.pdbqt` | el receptor del cribado sin las 1.763 aguas |
