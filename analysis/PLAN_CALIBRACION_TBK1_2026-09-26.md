# Plan de calibracion del embudo: la diana elegida es TBK1

26 de septiembre de 2026

## De donde viene esto

La auditoria de ayer (`AUDITORIA_ACTIVOS_TDP43_2026-09-26.md`) dejo cerrado que
contra el bolsillo de RRM2 de TDP-43 **no hay con que calibrar**: no existe un solo
ligando potente, medido y publicado contra ese sitio. Con cuatro ligandos de
micromolar alto, cualquier AUC que salga es ruido.

Lo que queda por decidir es que se hace con esa conclusion. Este documento es esa
decision, y esta tomada con datos, no con criterio.

## 1. Que dianas de ELA tienen quimica de verdad publicada

Se han contado, gen por gen, las medidas de afinidad (IC50, Ki, Kd, EC50, con valor
pChEMBL) que cuelgan de la diana **humana y de proteina unica** en ChEMBL. El script
es `buscar_diana_calibracion.py` y este es su resultado.

| gen | ChEMBL | actividades | con afinidad | compuestos | nombre |
|---|---|---|---|---|---|
| HDAC6 | CHEMBL1865 | 29.560 | 3.586 | 2.745 | Protein deacetylase HDAC6 |
| TBK1 | CHEMBL5408 | 13.123 | 1.342 | 1.175 | Serine/threonine-protein kinase TBK1 |
| VCP | CHEMBL1075145 | 1.486 | 721 | 595 | Transitional endoplasmic reticulum ATPase |
| NEK1 | CHEMBL5855 | 2.003 | 305 | 281 | Serine/threonine-protein kinase Nek1 |
| SOD1 | CHEMBL2354 | 175 | 28 | 28 | Superoxide dismutase [Cu-Zn] |
| STMN2, SQSTM1, PFN1, FUS | — | — | — | 2 cada una | — |
| UBQLN2, KIF5A | — | — | — | 1 cada una | — |
| **TARDBP** | CHEMBL2362981 | 40.125 | **0** | **0** | TAR DNA-binding protein 43 |
| ATXN2, ANXA11 | — | — | **0** | **0** | — |
| OPTN, C9orf72 | — | — | sin diana humana de proteina unica | | |

La fila que importa es la de TARDBP: **de 40.125 actividades, cero son una medida de
afinidad**. No es que sean pocas: es que no hay ninguna. Toda la quimica que cuelga de
TDP-43 en ChEMBL son porcentajes de un ensayo de crecimiento en levadura, que no mide
union. Eso confirma, por una via independiente, lo que ya decia la auditoria.

Y las demas dianas de ELA tampoco salvan la situacion: SOD1 tiene 28 compuestos con
afinidad, pero el propio `verdad_de_referencia.py` documenta que 16 de sus 20 activos
son el mismo nucleo pirazolona, o sea una sola serie, y fue precisamente la que
inflaba el AUC de 0,601 a 0,815. SOD1 no calibra nada. OPTN, C9orf72, ATXN2 y ANXA11
no tienen ni diana de proteina unica. FUS, STMN2, SQSTM1, PFN1, UBQLN2 y KIF5A tienen
dos compuestos o uno.

Quedan cuatro candidatas de verdad: HDAC6, TBK1, VCP y NEK1.

## 2. Por que TBK1 y no las otras tres

**TBK1 es una diana de ELA, no una diana de conveniencia.** Las mutaciones de perdida
de funcion de TBK1 causan ELA (Cirulli et al., Science 2015, el exoma de ELA que
identifico TBK1). Es la misma proteina que el proyecto ya estudia por su papel en la
autofagia y en la via de NF-kB.

**Es el caso en el que el acoplamiento deberia funcionar mejor, y eso es justo lo que
se quiere.** El bolsillo de ATP de una quinasa es lo mas parecido a un caso facil que
hay para Vina: un hueco profundo, cerrado, sin metales, con quimica de tipo I bien
conocida. Si el embudo del proyecto **falla incluso ahi**, la conclusion es dura y
limpia: el limite es la herramienta, no la diana. Si **acierta ahi y falla en
TDP-43**, entonces el problema es la forma del bolsillo de RRM2, que es una superficie
poco profunda, y eso es un resultado publicable y util. Cualquiera de las dos
respuestas sirve. Esa es la razon de elegir el caso facil: para que el experimento
tenga las dos salidas cargadas.

**Tiene la estructura con ligando que hace falta, y a buena resolucion.** Se han
buscado en el PDB todas las estructuras humanas de TBK1: hay 25, y **20 con un ligando
de verdad** cocristalizado. La mejor es la **4EUU, a 1,80 angstrom**, que es el
dominio quinasa de TBK1 humano (residuos 2 a 308) con **Ser172 fosforilado** y con
**BX-795** en el bolsillo de ATP (`TITLE: STRUCTURE OF BX-795 COMPLEXED WITH HUMAN TBK1
KINASE DOMAIN`). Eso significa que la caja se puede definir **desde el cristal**, y
esto merece parrafo aparte, mas abajo.

**Por que no las otras tres:**

- **HDAC6** tiene mas quimica que ninguna (2.745 compuestos) y estructuras a 1,05
  angstrom, pero **tiene zinc en el fondo del bolsillo**. Vina y Vinardo no tienen
  termino de metal en su funcion de puntuacion, y casi todos los inhibidores de HDAC
  quelan ese zinc. Un fallo ahi no distingue entre "el embudo no ordena" y "el motor no
  sabe de metales". Mezclaria dos causas y no serviria para calibrar.
- **VCP** tiene 595 compuestos y 54 estructuras con ligando, y es buena candidata,
  pero sus inhibidores se unen sobre todo al bolsillo de nucleotido del dominio D2, que
  es un sitio alosterico y mas dificil de definir. **Se deja como replicacion
  independiente**, no como primera prueba.
- **NEK1** tiene 281 compuestos, pero **solo 2 estructuras en todo el PDB y una sola
  con ligando** (4B9D, 1,90 angstrom). Es poco sitio para definir una caja con criterio.
  Queda descartada por falta de estructura, no por falta de quimica.

## 3. La caja, sacada del cristal y no a ojo

Esto es lo mas importante del plan, porque es exactamente donde el proyecto se equivoco
con TDP-43. Alli habia tres cajas (A, B y C) y se eligio por AUC: la caja C era la que
daba 0,742 y por eso se adopto. **Eso es circular**: se elige la caja con la regla que
luego se presume validada por esa misma caja. Aqui no se puede hacer eso, y no se va a
hacer.

La caja de TBK1 queda fijada por el ligando cocristalizado, con un criterio que no
mira ningun resultado de acoplamiento:

- centro del ligando BX-795 en la cadena A: **-0.45, -9.90, 9.00**
- centro de los 15 residuos del bolsillo a menos de 8 angstrom: **-1.44, -10.01, 8.66**
- residuos que forman el bolsillo, por distancia al ligando: Leu15, Gly16, Gln17,
  Val23, Phe24 (asa P), Phe88, Cys89, Gly92, Ser93 (bisagra), Ala36, Gly139, Asn140,
  Met142, Thr156, Ile14
- la extension de esos residuos es de 11,4 por 12,0 por 13,8 angstrom

Dos comprobaciones de que el sitio es el correcto: la **Cys89 es la cisteina de la
bisagra** de TBK1, y la **Met142 es la puerta** (gatekeeper), que en TBK1 es metionina
en vez de la treonina habitual, y es la razon de que TBK1 tenga un bolsillo mas ancho
de lo normal. Los dos residuos estan en la lista, o sea que el ligando esta donde tiene
que estar.

El centro recomendado es **-1.44, -10.01, 8.66**, que es el del bolsillo y no el del
ligando entero, porque BX-795 asoma una cola al disolvente y su centro se corre hacia
fuera. Tamaño de caja recomendado: **24 angstrom** de lado, que deja cuatro o cinco
angstrom de margen alrededor del ligando. El resto del protocolo se copia tal cual del
que ya esta validado: SEARCH_DEPTH 20, NUM_MODES 3, THREAD 8000.

## 4. El banco de prueba

La quimica medida contra TBK1 en ChEMBL (CHEMBL5408), quedandose con el mejor valor por
compuesto y solo con afinidades de verdad (IC50, Ki, Kd, EC50), son **561 compuestos**.
Repartidos por potencia:

| franja | compuestos |
|---|---|
| 10 nanomolar o mejor (pChEMBL >= 8) | 82 |
| de 10 a 100 nanomolar | 117 |
| de 0,1 a 1 micromolar | 232 |
| de 1 a 10 micromolar | 110 |
| peor que 10 micromolar | 20 |

La comparacion con TDP-43 no admite discusion: alli hay 0 compuestos con afinidad
medida, aqui hay **199 por debajo de 100 nanomolar**. Ese es todo el cambio.

El lado de los negativos hay que decirlo con honestidad: por el criterio estricto solo
20 compuestos se han ensayado contra TBK1 y han salido peores que 10 micromolar. No
son muchos, y por eso **el fondo del banco no van a ser negativos medidos, sino
decoys emparejados en propiedades**, exactamente como se hace en DUD-E. Lo que cambia
respecto a TDP-43 no son los negativos: son los positivos, que ahora existen.

Una advertencia que hay que respetar al leer los resultados: la quimica de TBK1 es
mucho mas congenerica que la de una libreria diversa. Es decir, es muy posible que
muchos de los 199 activos sean el mismo esqueleto con grupos cambiados. Si se mide
enriquecimiento sin tener eso en cuenta, sale inflado. **El remedio ya esta escrito en
el proyecto**: la regla de decision de `analysis/regla_decision/` no da un AUC global,
da un AUC por quimiotipo y exige que pasen varios. Se usa esa misma, sin inventar otra.

## 5. Que se mide, y que significa cada resultado

Se corre el mismo embudo (Vina-GPU, mismo protocolo) sobre el banco de TBK1 y se
calculan **las mismas metricas que el proyecto ya tiene escritas**, para que sean
comparables una a una con TDP-43: AUC por atomo pesado contra el fondo duro, AUC del
residuo, y la ordenacion emparejada por tamaño con su z. Se añaden **factores de
enriquecimiento (EF1%, EF5%) y BEDROC**, que es lo que mide "puedo sacar una lista
corta", que es la pregunta de verdad.

La lectura es esta, y conviene dejarla escrita antes de ver el resultado para no
torcerla despues:

- **Si TBK1 sale bien** (activos arriba, EF1% claramente mayor que 1, y sobre todo
  varios quimiotipos pasando la regla de decision), entonces el embudo funciona y el
  problema de TDP-43 es que su bolsillo no tiene quimica de referencia. La conclusion
  es que **el proyecto puede seguir con TDP-43 solo como catalogo**, y que la linea de
  trabajo pasa a ser conseguir un positivo real, no afinar nada mas.
- **Si TBK1 sale mal**, entonces el embudo no ordena ni en el caso facil, y lo que hay
  que rehacer no es la diana sino la herramienta: pasaria a tener sentido probar
  Vinardo, el rescoring por interacciones, o los dos juntos, **contra este banco**, que
  ahora si tiene respuesta conocida. Se afinaria contra TBK1 y se volveria a TDP-43
  cuando hubiera referencia.
- **Si sale a medias** (bien en AUC global y mal en EF1%), el diagnostico es que ordena
  en grueso pero no discrimina arriba, y entonces el uso honesto del acoplamiento es
  como prefiltro amplio, tal como decia el informe de la lista corta, nunca como
  selector de candidatos.

Ninguna de las tres salidas deja al proyecto donde estaba. Eso es lo que se buscaba.

## 6. Lo que este plan no hace

- **No acopla los dos millones de ZINC.** Sigue en pie lo de la auditoria: multiplicar
  por veinte un selector que no selecciona no arregla nada. Aquello se retoma cuando el
  embudo tenga respuesta conocida, y entonces tendra sentido.
- **No tira el cribado de 108.841.** Se queda como catalogo, con las tres ordenaciones
  guardadas (`ranking_TDP43_v2_por_atomo.csv`, `residual_TDP43_v2.csv`,
  `emparejado_TDP43_v2_15_50atomos.csv`).
- **No afina la metrica.** Ese era el vicio de fondo: se han escrito reglas de
  decision, barridos de exhaustividad, controles de motor y un modelo ExtraTrees con
  R2 de 0,26, todo sobre una vara que no medía. Aqui se para hasta tener el banco.
- **No empieza el acoplamiento hoy.** El plan queda escrito; la maquina se apaga.

## 7. Ficheros

- `analysis/PLAN_CALIBRACION_TBK1_2026-09-26.md` — este documento
- `analysis/buscar_diana_calibracion.py` — el recuento por gen que decide la diana
- `analysis/AUDITORIA_ACTIVOS_TDP43_2026-09-26.md` — por que TDP-43 no vale
- `analysis/fragmentos_nshogoza.csv` — los tres fragmentos, verificados contra la figura

Los scripts de las tres ordenaciones de la libreria quedan en `gpu_dock/`
(`rankear_libreria_por_atomo.py`, `residual_libreria.py`, `rankear_por_tamaño.py`),
que el proyecto mantiene fuera de git a proposito por ser datos masivos. Los resultados
que producen son los tres CSV de `gpu_dock/resultados_libreria/`.

Datos de partida verificados hoy: PDB 4EUU (dominio quinasa de TBK1 humano, residuos
2-308, Ser172 fosforilado, BX-795 cocristalizado, 1,80 angstrom); las 25 estructuras
humanas de TBK1 y las 20 con ligando; ChEMBL CHEMBL5408 con 561 compuestos con
afinidad medida y el reparto por potencia de la tabla de arriba.
