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

## 7. Archivos del dia

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
