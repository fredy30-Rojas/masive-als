# VALIDACION ESTANDAR: DUD-E CDK2 — estado del 29 de septiembre de 2026

## 1. La pregunta del dia, que no es la del 28

El 28 de septiembre el banco de TBK1 dejo un techo: AUC crudo 0,602, AUC por atomo
pesado 0,617, 0,637 con la media de las nueve poses, y el MM-GBSA en 0,498. Con
esos numeros solos **no se sabe que significa el techo**: puede ser que TBK1 sea
una diana dificil, o puede ser que el embudo (Vina-GPU con `search_depth 20`) no
de mas de si.

Para distinguirlo hace falta una diana donde la respuesta no la pongamos nosotros.
CDK2 de DUD-E es eso: 474 activos y 27.850 señuelos emparejados por fisicoquimica y
por forma, publicados con la estructura 1h00, y con **un numero de referencia
encima de la mesa: AUC 0,791 con el DOCK 3.6 de los propios autores** (Mysinger y
col., *Directory of Useful Decoys, Enhanced*, J. Med. Chem. 2012).

La lectura, decidida hoy antes de ver nada:

* CDK2 en ~0,79 → el embudo esta bien y el techo de 0,62 es de TBK1;
* CDK2 en ~0,60 → el techo es del MOTOR, no de la diana, y lo que hay que cambiar
  es el motor;
* CDK2 por debajo de 0,55 → el embudo no separa ni donde ya esta hecho y otros
  separan, y todo lo que se apoya encima (rescoring, MM-GBSA) se apoya en nada.

## 2. Lo que quedo montado

* **Receptor** (`preparar_receptor_cdk2.py`): el `receptor.pdb` de DUD-E es de
  formato viejo (54 columnas, sin columna de elemento) y trae el residuo 83
  rotulado **LEV** cuando es una leucina entera (N CA C O CB CG CD1 CD2). Se
  normaliza el elemento, se quitan los 479 hidrogenos viejos y se renombra LEV a
  LEU. Quedan 2.244 atomos pesados y el PDBQT de 2.723 atomos.
* **Caja**: 24 A centrada en el ligando del cristal, `[2,60, 26,32, 8,59]`. El
  centro del bolsillo (residuos a menos de 8 A) cae a 0,93 A de ahi: esta bien
  puesta.
* **Banco** (`preparar_banco_cdk2.py`): **28.301 PDBQT de 28.324** en 18 minutos
  (474 activos y 27.846 señuelos; 4 fallos, los cuatro señuelos, por geometria 3D
  imposible). Sin un solo descartado por tipo de atomo: todo lo entiende AutoDock.
* **Lanzador** (`lanzar_banco_cdk2.py`): el de TBK1 copiado, con `--decoys N` para
  acotar los señuelos con semilla fija. Reanudable, lotes de 500 con volcado cada
  15 s.
* **Validador** (`validar_banco_cdk2.py`): las mismas medidas que el de TBK1 (AUC
  crudo y por atomo pesado, EF1 %, EF5 %, BEDROC con su base al azar medida,
  bootstrap de 2.000 remuestreos, reparto por esqueletos de Murcko) mas la
  comparacion contra el 0,791 publicado y contra el techo de TBK1 recalculado con
  el mismo codigo. Comprobado que los 28.301 identificadores del manifiesto
  empalman con su SMILES del `.ism`: **cero sin empalmar**.

## 3. El hallazgo del dia: el control estaba mal medido

El control de redocking (acoplar el ligando del cristal y ver a que distancia cae)
daba 6,6 A en las seis variantes del barrido. Se dio por fallado. **La medida
estaba rota**, y por esto:

El script suponia —lo decia en su propio docstring— que "meeko respeta el orden de
atomos de la molecula y Vina-GPU mantiene ese orden en las poses, asi que el
emparejamiento por indice es una resta directa". **No es verdad.** Meeko escribe
los atomos pesados en otro orden: el PDBQT del ligando del cristal, leido en orden
de fichero, empieza `C,C,C,C,C,C,O,C,C,C,N,...` y el mol2 empieza
`C,C,C,C,C,C,C,C,C,...`. Restando por indice, cada coordenada caia en el atomo
equivocado: lo que se medía no era la pose, era una molecula con la geometria
desordenada.

El arreglo no es suponer mejor: es **deducir** la correspondencia. El PDBQT que se
acopla es el ligando del cristal, asi que se busca la asignacion optima entre sus
atomos y los del mol2 exigiendo que cada uno caiga en uno del mismo elemento. Da
**0,0000 A** y devuelve la permutacion exacta
(`5 6 7 8 13 14 27 12 20 11 25 0 1 26 23 17 9 21 18 10 22 24 19 15 16 3 4 2 28 29`).
Con ella, la resta por indice vuelve a significar lo que dice. Si esa asignacion no
diera ~0, el control se para: significaria que lo que se acoplo no era la pose del
cristal.

### La segunda trampa, que casi cuela: `GetBestRMS` superpone antes de medir

Al arreglar el orden aparecieron dos numeros que no cuadraban entre si: 6,87 A por
indice y 1,55 A por `GetBestRMS`. La explicacion es que **`GetBestRMS` de RDKit
alinea la pose sobre el cristal y despues mide**. Eso responde "se conserva la
forma?", no responde "esta en el mismo sitio del bolsillo?", que es lo que se
pregunta un redocking. En CDK2 la forma se conserva (1,55-2,30 A) y el sitio no
(3,81-6,87 A).

Asi que la medida que decide, y la que usa el liston de 2 A, es la del **bolsillo**,
sin alinear. El RMSD alineado se deja al lado, etiquetado como lo que es.

### Y el veredicto, medido bien, no cambia

Barrido remedido desde disco, sin volver a usar la GPU (`diagnostico_orden_poses.py
--remedir-barrido`, `barrido_redocking_cdk2_corregido.txt`):

| variante | entrada | bolsillo | **piso** | alineado | veredicto |
|---|---|---|---|---|---|
| caja 20 A, depth 20 | 0,0000 A | 3,81 A | 3,67 A | 2,30 A | no pasa |
| caja 20 A, depth 32 | 0,0000 A | 6,55 A | 5,39 A | 1,83 A | no pasa |
| caja 22 A, depth 20 | 0,0000 A | 6,87 A | 5,69 A | 1,55 A | no pasa |
| caja 22 A, depth 32 | 0,0000 A | 3,99 A | 3,70 A | 2,30 A | no pasa |
| caja 24 A, depth 20 | 0,0000 A | 6,41 A | 5,52 A | 1,88 A | no pasa |
| caja 24 A, depth 32 | 0,0000 A | 6,85 A | 4,60 A | 1,88 A | no pasa |

`piso` es el RMSD mas pequeño que puede dar **cualquier** emparejamiento de atomos
que respete el elemento (asignacion optima libre). No es la medida del proyecto: es
un suelo, y sirve para decir lo mas fuerte que se puede decir. Aqui dice que **ni
con el mejor emparejamiento imaginable la pose baja de 3,67 A**. El veredicto no
depende del orden de los atomos: el motor no encuentra esa pose.

## 4. El control de TBK1 no tenia este problema

Importa saberlo, porque si el liston del proyecto estuviera roto por todas partes,
nada de lo de TBK1 valdria. `control_redocking_bx795.py` casa **por nombre de
atomo** (`N07`, `C26`...) y lo explica en su codigo: los nombres del cristal si se
conservan, y el emparejamiento por nombre es exacto. Comprobado hoy: la entrada de
ese control da 5,72 A por indice y **0,0000 A** por correspondencia medida, o sea
que la suposicion del indice tambien fallaba alli, pero el script no la usaba. Su
1,23 A esta bien medido.

La diferencia de fondo: el ligando de CDK2 sale de un **mol2**, que no trae nombres
de atomo y deja que meeko los sustituya por el elemento (todos los carbonos se
llaman `C`). El de TBK1 sale de un **PDB**, que los trae. Por eso uno se puede
casar por nombre y el otro no.

## 5. Lo que queda por hacer y que se leera al final

1. Acoplar la submuestra 1:10 (474 activos + 4.740 señuelos, unas 4 horas) y,
   despues, los 27.846 señuelos enteros, sobre lo ya hecho.
2. Pasar `validar_banco_cdk2.py` y leer si el AUC por atomo pesado se parece al
   0,79 publicado, al 0,617 de TBK1, o a nada.
3. **Escribir el control de CDK2 en el informe junto al AUC.** El control no pasa y
   eso no invalida la medida del AUC (son dos preguntas distintas: "encuentra esta
   pose?" y "ordena activos sobre señuelos?"), pero no puede aparecer una sin la
   otra. El protocolo del proyecto lo tiene como requisito y aqui hay que decir
   por que se sigue adelante de todas formas.

## 6. Estado de la maquina

* **TBK1**: 30.334 poses de 35.207 a las 18:03, unos 2,70 s por ligando, ~3,7 h por
  delante. Vivo (PID 28900) y con la GPU al 100 %.
* **Encadenador** (`encadenar_cdk2.py` + `encadenar_cdk2.bat`, lanzado a las
  18:23 oculto y esperando): cuando el lanzador de TBK1 reescriba
  `resultados_tbk1.csv` —o si su registro lleva una hora parado, que es como se
  detecta un banco muerto— lanza `lanzar_banco_cdk2.py --decoys 4740` y deja la
  submuestra acoplada sin que nadie tenga que estar mirando la hora.
* **GPU**: RTX 4080, 8.128 MiB en uso.

## 7. Por que el motor no devuelve esa pose (investigado el mismo dia)

El control falla y hay que saber por que antes de tocar nada. Son tres las
sospechosas, y solo una es cierta: **busqueda** (la pose esta pero no la
encuentra), **puntuacion** (la encuentra y la descarta) o **preparacion** (la pose
no cabe en el receptor que hemos preparado). `diagnostico_pose_cristal.py` las
separa sin usar la GPU, y el resultado es tajante.

**La preparacion queda descartada.** La pose del cristal cabe: 17 residuos en
contacto a menos de 4 A, distancia minima de atomos pesados 2,12 A (un enlace de
hidrogeno, no un choque), y el anclaje canonico del sitio de ATP esta ahi entero,
LEU83 a 2,60 A, GLU81 a 2,85 A, HIE84 y LEU134 tambien. Seis contactos polares. Si
el motor huyera de esa zona por un choque, se veria aqui, y no se ve.

**Y el motor no pierde el bolsillo.** El nucleo de diaminopirimidina se queda
anclado al enganche en todas las variantes (LEU83 y LEU134 se mantienen siempre, y
el centroide del ligando se queda a 1,2-3,3 A), pero el brazo saturado con el
amonio se va de **5,9 a 9,4 A**. No es un ligando perdido: es un ligando anclado
por un lado y con el otro brazo puesto en otro sitio.

| fragmento | se mueve (peor caso de las 6 variantes) |
|---|---|
| nucleo pirimidina | 1,45 a 4,56 A |
| fenilo difluorado | 2,34 a 8,16 A |
| anillo fusionado | 2,36 a 6,58 A |
| **resto saturado + amonio** | **5,87 a 9,39 A** |

**Y el sitio donde lo pone no es peor en contactos.** El cristal hace 6 contactos
polares; las poses del motor hacen entre 4 y 9, con el mas corto a 2,79-3,05 A
frente a los 2,60 A del cristal. Es decir: el motor encuentra acomodos **iguales o
con mas contactos**, y los prefiere. La funcion de puntuacion no esta echando de
menos una interaccion que no ve; esta pesando distinto las que si ve.

**Tampoco es un conformero forzado.** Energia MMFF94 del mismo esqueleto, con los
atomos pesados fijos y los hidrogenos regenerados y relajados: el conformero del
cristal da **49,19 kcal/mol** y las poses del motor **58,67 a 70,25**. El cristal
esta entre 9 y 21 kcal/mol **menos** tenso que lo que devuelve el motor. No se
descarta la pose del cristal por estar forzada.

Con todo eso, la lectura es **(b) PUNTUACION**, con una salvedad honesta que hay
que dejar escrita: la busqueda no se ha agotado del todo (se probo
`search_depth` 20 y 32). Queda un matiz mas que no se puede comprobar con lo que
hay en disco y que conviene no perder de vista: el brazo del amonio es un cation, y
el `receptor.pdb` de DUD-E **no trae ni una sola molecula de agua** (cero HETATM),
asi que si en el cristal ese cation esta sujeto por una red de aguas, ese anclaje
se ha borrado antes de empezar. Es una hipotesis, no una medida, y va anotada como
tal.

### Lo que el propio cristal no decidió (segunda pasada del mismo día)

Se bajó el PDB original entero (`1H00_cristal_completo.pdb`) para comprobar una
hipótesis: que el amonio estuviera sujeto por una red de aguas que DUD-E quitó. La
hipótesis es **falsa**, y al buscarla apareció algo bastante mas gordo.

**El cristal trae el ligando DOS veces.** `FAP` (residuo 1300, altloc A) y `FCP`
(residuo 1400, altloc B), las dos con **ocupación 0,50** y B medio ~47, a 0,40 Å de
centroide. Comprobado: FAP es **exactamente** el mol2 de DUD-E (0,000 Å) y FCP está
a 1,41 Å. No son dos moléculas: son **dos conformaciones alternativas del mismo
ligando**, que el cristal no pudo decidir.

Y la diferencia entre ellas se concentra en un sitio muy concreto:

| parte | FAP contra FCP |
|---|---|
| núcleo de pirimidina | **0,00 Å** (idéntico) |
| anillo fusionado | **0,00 Å** (idéntico) |
| fenilo difluorado | 1,72 Å |
| brazo saturado + amonio | 1,90 Å |

**El amonio no está sujeto por nadie.** Medido sobre el PDB original: no tiene
ningún átomo de proteína a menos de 4 Å (el más cercano, a 4,17 Å) y **ninguna agua
a menos de 4,5 Å** (la más cercana, a 5,90 Å). Es un catión desolvatado, sin ningún
enganche, en la zona que el propio cristal dejó sin determinar. La hipótesis del
agua queda refutada, pero por el camino se entendió lo importante: **es justo el
grupo que el motor mueve de 6 a 9 Å**.

### Entonces el listón de 2 Å sobre los 30 átomos no vale para este ligando

Pedirle al motor que clave con 2 Å un grupo que (a) el cristal modela en dos
posiciones al 50 %, (b) no toca a nadie y (c) está a 4,17 Å de lo más cercano, es
pedirle que acierte algo que la estructura no sabe. El control hay que leerlo por
partes —núcleo anclado contra brazo suelto—, que es lo que hacen
`cristal_1h00_dos_conformaciones.py` y su informe `cristal_1h00.txt`:

| variante | pose | **núcleo (12 átomos)** | brazo (18 átomos) |
|---|---|---|---|
| caja 20 Å, depth 20 | 6 | **1,95 Å — PASA** | 4,65 Å |
| caja 20 Å, depth 32 | 6 | 4,97 Å | 9,47 Å |
| caja 22 Å, depth 20 | 7 | 5,06 Å | 9,41 Å |
| caja 22 Å, depth 32 | 9 | **2,00 Å — PASA** | 4,88 Å |
| caja 24 Å, depth 20 | 2 | 5,25 Å | 8,24 Å |
| caja 24 Å, depth 32 | 9 | 4,91 Å | 7,88 Å |

Sobre el núcleo, que es lo único que el cristal tiene clavado y lo que se ancla al
enganche de LEU83 y GLU81, el motor **sí lo reproduce**: 1,95 y 2,00 Å en dos de las
seis variantes. Sobre el brazo, ninguna. Así que el control no se puede leer como un
sí o un no sobre 30 átomos: leído por partes, el núcleo pasa.

Esto no borra lo de antes, lo ordena: la preparación está descartada, el motor no
pierde el bolsillo, y lo que no acierta es la colocación de un catión que el cristal
tampoco tiene resuelto. Queda por ver si con la búsqueda agotada aparece el núcleo en
mas variantes.

### La prueba que lo cierra, ya encadenada

Agotar la busqueda en la misma caja del banco: `search_depth 128` (y 9 poses) con
caja de 24 A y de 20 A. Si la pose del cristal aparece, era busqueda. Si no
aparece ni asi, es puntuacion y no hay mas que hablar. Son unos minutos de tarjeta
y ya esta metido en el encadenador, **antes** del banco: despues serian cuatro
horas de espera para tres minutos de prueba.

### Que implica, si se confirma que es puntuacion

Que el control no pase **no es un defecto del receptor ni de la caja**: es que este
motor, con esta diana y con este ligando cationico, prefiere un acomodo distinto
del mismo bolsillo. Para el objetivo del dia (medir el embudo contra un banco
publico) eso es informacion, no un bloqueo: el AUC se mide igual y se lee con esto
al lado. Para arreglarlo haria falta un rescoring, y aqui hay que ser claro: el
proyecto ya sabe que su MM-GBSA no ordena (AUROC 0,498 en TBK1), asi que no es el
rescoring que se va a probar por probar. Lo que esta prueba pide, si se confirma,
es un motor con desolvacion de verdad (el CNN de GNINA, que ya esta en el proyecto,
o un campo de fuerzas explicito), no otra capa encima de Vina.

## 8. Archivos del dia

| fichero | que es |
|---|---|
| `preparar_receptor_cdk2.py` | receptor y caja de CDK2, con el arreglo del PDB viejo y del LEV |
| `preparar_banco_cdk2.py` | los 28.301 PDBQT del banco |
| `lanzar_banco_cdk2.py` | acopla el banco, reanudable, con `--decoys` |
| `validar_banco_cdk2.py` | las medidas del banco y la comparacion con DUD-E y TBK1 |
| `control_redocking_cdk2.py` | el control, con la correspondencia deducida y el piso |
| `diagnostico_orden_poses.py` | la prueba del orden de atomos y la remedicion del barrido |
| `encadenar_cdk2.py` / `.bat` | el encadenado detras de TBK1 |
| `barrido_redocking_cdk2_corregido.txt` | el barrido remedido |
| `diagnostico_orden_poses.txt` | el detalle de la prueba del orden |
| `diagnostico_pose_cristal.py` | por que el motor no devuelve la pose: fragmentos, contactos, tension |
| `cristal_1h00_dos_conformaciones.py` | las dos conformaciones del cristal y el nucleo separado del brazo |
| `1H00_cristal_completo.pdb` | el PDB original entero, con las 209 aguas y las dos copias del ligando |
