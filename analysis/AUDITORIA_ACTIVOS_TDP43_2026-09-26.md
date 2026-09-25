# Auditoria de los activos de TDP-43 y de lo que hay publicado

26 de septiembre de 2026

## Por que se hace esta auditoria

Ayer se comprobo que el cribado de la libreria (108.841 compuestos de TDP43_v2) no
sirve para elegir candidatos: los unidores conocidos de TDP-43 caen todos en el 14%
de cabeza, y cualquier lista corta que deje menos de quince mil compuestos deja fuera
al mejor de ellos. La explicacion que se dio fue que un embudo que acierta con AUC
0,742 contra el fondo duro no puede rankear una libreria. Pero eso deja una pregunta
mas de fondo sin contestar: **¿es que la herramienta es floja, o es que la vara de
medir esta mal?**

Esta auditoria contesta esa pregunta. Se ha ido a las fuentes primarias:
ChEMBL, PubMed, Europe PMC y la figura original del articulo de los fragmentos.

## 1. Lo que hay publicado sobre TDP-43 en ChEMBL

La diana humana esta en ChEMBL como **CHEMBL2362981**. Cuelgan de ella 40.125
actividades, que suenan a muchisimo para una diana tan dificil. No lo son: el
desglose por ensayo es este.

| ensayo | que mide | actividades |
|---|---|---|
| CHEMBL2354287 | qHTS de inhibidores de TDP-43 (PubChem 652104) | 40.105 |
| CHEMBL4627527 | inhibicion de la union de TDP-43 al ARN | 5 |
| CHEMBL4627528 | afinidad de union a TDP-43 | 1 |
| CHEMBL4627529 | union al ADN en cuadruple helice del repetido G4C2 | 6 |
| CHEMBL5392564 | morfologia de anisosomas en celulas U2OS | 3 |
| CHEMBL5392565 | empalme del ARNm de STMN2 en SH-SY5Y | 5 |
| CHEMBL5156952 | fosforilacion de TDP-43 en linfoblastos | 1 |
| CHEMBL5156953 | nivel de TDP-43 en linfoblastos | 1 |
| CHEMBL5652589 | afinidad por Kinobead pull-down | 4 |

Hay que mirar cada uno antes de creerselo, y ninguno aguanta como fondo o como
positivo de un cribado contra el bolsillo de RRM:

**El cribado masivo no vale.** El ensayo CHEMBL2354287 es el PubChem AID 652104, de
Susan Lindquist en el Whitehead Institute. Es un **ensayo de crecimiento en levadura**
que expresa TDP-43: se mide la luminiscencia de un reactivo de viabilidad, no una
union. Un compuesto inactivo ahi es un compuesto que no rescato el crecimiento de la
levadura, y eso puede ser por mil motivos que no tienen nada que ver con el bolsillo
de RRM. Sus 40.105 registros parecen un fondo experimental enorme y no lo son.

**El ensayo de 900 nanomolar tampoco.** CHEMBL4627529 mide la union a la
**cuadruple helice de ADN del repetido G4C2**, o sea al acido nucleico, no a la
proteina. Sus Kd de 900 nM a 16 uM (elliptina, datelliptium y compania) son afinidades
por el ADN repetido. Meterlas como positivos de TDP-43 seria confundir la diana con
su ligando.

**Lo que queda de union real a la proteina es esto:**

| compuesto | afinidad | tecnica | sitio |
|---|---|---|---|
| rTRD01 | Kd 89,4 uM; IC50 150 uM | MST y RMN | RRM1 |
| nTRD22 | Kd 145 +- 3 uM (1-260) y 96 +- 36 uM (1-102) | MST | dominio N-terminal |
| CHEMBL4638490 | Kd 89 uM | — | (es rTRD01) |
| CHEMBL5653589 | Kd 176 nM | Kinobead pull-down | no especificado |
| CHEMBL3752910 | Kd 813 nM | Kinobead pull-down | no especificado |

Las dos entradas de Kinobead son de un ensayo de arrastre de quinasas: que una
molecula aparezca ahi no significa que se una a un sitio concreto de TDP-43, y menos
al bolsillo de RRM.

**Y hay una entrada que no cuadra.** CHEMBL4635203 aparece con IC50 = 100 nM en el
ensayo de inhibicion de la union al ARN, del mismo articulo de revision que da los
otros valores. En la misma tabla convive con un ">100 nM" del compuesto
CHEMBL4647070, con cafeina a 10 uM y con acido tioctico a 19,95 uM. Que la unica
entrada nanomolar salga junto a un ">100" del mismo orden, y rodeada de compuestos
que son interferentes clasicos de ensayo (la cafeina y el acido lipoico unen acidos
nucleicos de forma inespecifica), apunta a un error de unidades en la extraccion de
la revision: lo mas probable es que la tabla original estuviera en micromolar. No se
ha podido comprobar contra la tabla original porque la revision esta de pago.
**Queda marcado como no verificado y no se usa para nada.** Si alguien quiere
perseguirlo: seria el unico inhibidor nanomolar publicado de la union de TDP-43 al
ARN, y merece una carta a los autores o la busqueda del articulo primario.

## 2. La comprobacion que si se ha podido cerrar: los tres fragmentos

De los siete activos que el proyecto usa, los tres **fragmentos son los unicos que
estan medidos contra el bolsillo que de verdad se acopla** (RRM2, residuos G245,
E246, H256, I257 y S258, por desplazamiento quimico de RMN). El proyecto tenia dos
reparos escritos sobre ellos: que sus SMILES no estaban publicados y se habian
sacado por OCR quimico de un dibujo, y que nadie los habia revisado con ojos de
quimico.

**Eso queda resuelto hoy.** Se ha bajado la Figura 1 original del articulo
(Nshogoza et al., IJMS 2019, 20:3230, DOI 10.3390/ijms20133230) desde el servidor de
figuras de MDPI, que es abierto, y se ha mirado el panel (c) a resolucion completa.
Las tres estructuras coinciden con las que el proyecto tiene:

- **hit 1 / fragmento_1**: piperidin-3-amina unida por el nitrogeno a una pirimidina.
  SMILES `NC1CCCN(c2ncccn2)C1`, C9H14N4. Coincide.
- **hit 2 / fragmento_2**: un biciclo de cinco y cinco con **un imidazol fusionado a
  un 1,3,4-tiadiazol**, y un grupo 1-aminoetilo en el carbono del imidazol. SMILES
  `CC(N)c1cn2ncsc2n1`, C6H8N4S. Coincide, y la formula cuadra: contando el dibujo
  salen 6 carbonos, 4 nitrogenos y 1 azufre, que es exactamente C6H8N4S. Este era el
  unico de los tres que ofrecia duda, porque el OCR podia haber confundido el orden de
  los nitrogenos del anillo; con la figura a resolucion completa delante, el azufre
  esta arriba a la izquierda, un nitrogeno arriba a la derecha, y los otros dos
  abajo, uno en el centro y otro a la izquierda, que es justo lo que dibuja el SMILES.
- **hit 3 / fragmento_3**: 2-(4-aminopiperidin-1-il)nicotinonitrilo, con el nitrilo
  en la posicion contigua a la del enlace, no en la opuesta. SMILES
  `N#Cc1cccnc1N1CCC(N)CC1`, C11H14N4. Coincide.

Las tres imagenes de cotejo que se guardaron en `_nshogoza_fig1c/` eran correctas.
Las copias de la figura original y los recortes ampliados de cada estructura quedan
en esa misma carpeta, con los nombres `original_g001.jpg`, `original_general.png` y
`publicado_hit1.png`, `publicado_hit2.png`, `publicado_hit3.png`, para que cualquiera
pueda rehacer el cotejo sin volver a bajar nada.

**Lo que la lectura del articulo anade, y es importante:** los tres fragmentos **no
tienen ninguna afinidad medida**. El articulo no da Kd ni IC50 de ninguno de ellos.
Lo que da son desplazamientos quimicos, y los titula con el fragmento en exceso de
**ocho a uno** sobre la proteina para verlos. Eso, en el propio lenguaje del articulo,
es "three hits weakly binding"; la resolucion de las tecnicas que usan llega hasta
milimolar, y el fragmento esta a 0,4 mM en el cribado de 890 fragmentos. Los
fragmentos son señales de RMN, no ligandos con afinidad. Que queden en los puestos
125, 127 y 129 de 129 al acoplarlos no es un fallo del embudo: es que no hay nada que
ordenar.

## 3. La conclusion

Juntando las dos cosas, la respuesta a la pregunta de la que se partia es que **no es
la herramienta, es la vara**. El conjunto de referencia con el que el proyecto valida
su embudo contra el bolsillo de RRM de TDP-43 esta formado por:

- tres fragmentos sin ninguna afinidad medida, de 9, 6 y 11 atomos pesados, medidos
  solo por desplazamiento de RMN y a exceso ocho a uno;
- rTRD01, con Kd 89 uM, y medido contra **RRM1**, que no es el bolsillo de los
  fragmentos;
- nTRD22, que segun la fuente primaria se une al **dominio N-terminal** y es
  alosterico, y con el constructo de solo RRM los autores **no ven ningun
  desplazamiento**;
- PE859 y berberrubina, descritos en la interfaz RRM1-RRM2 y medidos en celula y en
  nematodo, sin Kd ni IC50, y que el propio proyecto caracteriza como cationes planos
  intercaladores.

No hay un solo ligando con afinidad submicromolar conocido, medido y publicado contra
el bolsillo de RRM2 de TDP-43. Y no es un descuido de la busqueda: la diana humana
entera, en ChEMBL, tiene menos de veinte medidas de union directa, y las mejores son
de micromolar alto. El NIH si tiene un proyecto financiado (5U01AG068823) buscando
justamente eso, moleculas que se unan a RRM2 e impidan la union al ARN, lo cual dice
que la cosa esta por hacer, no que este hecha.

Esto explica de una vez todos los resultados raros que el proyecto lleva meses
encontrando: que ninguna metrica cruda ordene, que el corte fijo de cabeza no
recoja ninguna quimia, que los fragmentos queden en el fondo de la lista, que los dos
unicos que pasan el criterio no se unan a ese bolsillo. No era el receptor, ni la
caja, ni la protonacion, ni el motor, ni la metrica. Era que **no hay con que
calibrar**. Con cuatro ligandos de micromolar alto, uno de ellos alosterico de otro
dominio, cualquier AUC que salga es ruido con dos decimales.

## 4. Lo que esto implica para el trabajo que viene

- **Acoplar los dos millones en la caja corregida no tiene sentido.** No porque este
  mal hecho, sino porque multiplica por veinte un selector que no selecciona, y lo
  medido es que no selecciona por falta de referencia, no por falta de compuestos.
- **El cribado de la libreria conserva un valor, que es el de un catalogo.** Saber
  que esos 108.841 compuestos se acoplaron y con que energia esta bien y es util como
  archivo, para volver a el el dia que exista una referencia de verdad. No es una
  lista de candidatos y no debe presentarse como tal.
- **Lo unico que desbloquea el proyecto es conseguir un conjunto de referencia
  honesto para ese bolsillo.** Hay tres caminos, y los tres son de meses, no de dias:
  buscar el articulo primario del dato de 100 nM y, si es real, comprarlo o sintetizar
  lo necesario para tener un positivo potente; esperar a que el proyecto del NIH
  publique; o cambiar de diana por una que si tenga quimica conocida, y dejar TDP-43
  para cuando la tenga.
- **La diana de calibracion ya esta elegida: TBK1.** Se ha pasado revista a los genes
  de ELA y se ha visto cual tiene quimica de verdad publicada: TBK1 tiene 199
  compuestos con afinidad medida por debajo de 100 nanomolar y una estructura a 1,80
  angstrom con inhibidor cocristalizado (4EUU, BX-795). El plan completo, con la caja
  sacada del cristal y no a ojo, esta en `PLAN_CALIBRACION_TBK1_2026-09-26.md`.
- **Y lo que no hay que hacer:** seguir afinando la metrica. Se han escrito reglas de
  decision, barridos de exhaustividad, controles de motor, auditorias de rescoring y
  un modelo ExtraTrees con R2 de 0,26. Todo eso es trabajo bien hecho sobre una vara
  que no mide. Afinar la regla no arregla el problema, porque el problema no esta en
  la regla.

## 5. Nota sobre el coste medido, para que no se cite mal

En la carpeta del cribado de la caja corregida quedo un log con la medida de ritmo:
28,47 segundos por ligando, que extrapolado da 896 horas (37 dias) por diana. Esa
medida se tomo con la tarjeta grafica ocupada por un modelo de Ollama que retenia
5,58 GB de memoria y frenaba el acoplamiento veinte veces. Con la tarjeta libre la
medida que se tomo despues fue de **1,40 segundos por ligando**, o sea unas 44 horas
por diana. La cifra de 896 horas no debe usarse.

## 6. Ficheros de esta auditoria

- `analysis/AUDITORIA_ACTIVOS_TDP43_2026-09-26.md` — este informe
- `analysis/_nshogoza_fig1c/original_g001.jpg` — la Figura 1 original del articulo
- `analysis/_nshogoza_fig1c/publicado_hit1.png`, `publicado_hit2.png`,
  `publicado_hit3.png` — cada estructura ampliada, para el cotejo
- `analysis/fragmentos_nshogoza.csv` — actualizado con la verificacion
- `analysis/INFORME_LISTA_CORTA_LIBRERIA_2026-09-25.md` — el informe de ayer
