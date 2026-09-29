# La diana donde el embudo se pueda validar, y por qué hCA2 — 29 de septiembre de 2026

## Qué se buscaba

Una diana con (1) **activos medidos** suficientes para que el intervalo del AUC sirva
para algo, y (2) un **farmaco co-cristalizado** para el control de redocking. Y una
condición que el proyecto ya se había ganado a pulso: que la serie sea **ganable** por
un cribado, que ordena por tamaño.

## Las que pasan los dos primeros requisitos

`elegir_diana_validable.py` lo cuenta por API de ChEMBL, sin fiarse de la memoria: los
nombres de la lista corta son hipótesis de partida y todo lo que se cuenta sale de
ChEMBL (actividades con `pchembl_value`, o sea potencia estandarizada) y de PDBe
(entidades de tipo `bound`, que es el farmaco en el bolsillo).

| diana | ChEMBL | activos medidos | entradas PDB |
|---|---|---|---|
| Beta-secretasa 1 | CHEMBL4822 | 15.327 | 431 |
| Anhidrasa carbonica 2 | CHEMBL205 | 13.051 | 1.135 |
| Anhidrasa carbonica 1 | CHEMBL261 | 10.922 | 54 |
| Acetilcolinesterasa | CHEMBL220 | 8.144 | 78 |
| Beta-secretasa 2 | CHEMBL2525 | 3.181 | 19 |
| Anhidrasa carbonica 14 | CHEMBL3510 | 1.031 | 2 |
| Anhidrasa carbonica 13 | CHEMBL3912 | 723 | 2 |
| Quimotripsina C | CHEMBL2386 | 141 | 7 |

Ocho de doce. Descartadas por la vía rápida: aldosa reductora (141 y sin cristal),
AChE de Drosophila y otras sin serie o sin estructura.

## El segundo filtro, que es el que decide

Tener muchos activos medidos no sirve si la serie es grande y flexible: el cribado
ordena por tamaño y la señal se la come el tamaño. Es exactamente lo que se midió hoy
con los fragmentos de TDP-43. Así que para cada candidata se midió el perfil de **su**
serie (mediana de una muestra de 250 activos reales):

| diana | muestra | atomos pesados | enlaces rotatorios | cLogP |
|---|---|---|---|---|
| Beta-secretasa 1 | 223 | **51** | **20** | 2,07 |
| Beta-secretasa 2 | 235 | 40 | 11 | 3,59 |
| Quimotripsina C | 135 | 34 | 9 | 2,96 |
| Anhidrasa carbonica 1 | 222 | 24 | 6 | −0,26 |
| **Anhidrasa carbonica 2** | 183 | **22** | **5** | 1,17 |
| Anhidrasa carbonica 13 | 176 | 23 | 5 | 2,04 |
| Anhidrasa carbonica 14 | 186 | 19 | 3 | 0,70 |
| Acetilcolinesterasa | 219 | 25 | 4 | 2,94 |

**Beta-secretasa 1 es la que más activos tiene y la peor candidate**: mediana de 51
atomos pesados y 20 enlaces rotatorios. Son peptidomiméticos grandes; un cribado sobre
ellos mide tamaño, no quimia, y se acaba donde se empezó hoy.

## La elegida: anhidrasa carbonica 2 (CHEMBL205)

- **13.051 activos medidos** con pchembl. Muy por encima de los ~40 que hacen falta
  para un intervalo de ±0,10, así que el margen deja de ser el problema.
- **Mediana de 22 atomos pesados y 5 enlaces rotatorios**: una serie pequeña y semirrigida
  (sulfonamidas), en el mismo rango de tamaño que los positivos de SOD1, que sí
  funcionaban dentro de su grupo de tamaño.
- **1.135 estructuras con farmaco**: el control de redocking tiene de sobra para elegir
  uno con la pericia que se usó con `andr`.

Respaldo: anhidrasa carbonica 14 (19 pesados, 3 rotatorios), con serie cuatro veces
peor pero más limpia; y acetilcolinesterasa, si hiciera falta una segunda diana.

## Tres avisos técnicos, todos aprendidos en este proyecto

1. **La serie de hCA2 es congenerica, y eso es un peligro conocido.** Miles de
   sulfonamidas distintas: es la forma exacta de construir un conjunto de validación donde casi
   todo el mundo comparte el mismo anillo, que es lo que se refutó en CDK2 (474 activos,
   474 esqueletos distintos, y aun así había que mirarlo). **Antes de acoplar nada hay
   que auditar la serie con `auditar_activos_cdk2.py`**: esqueletos de Murcko, MCS y
   solape con los señuelos. sin esa cifra, el AUC puede volver a ser un artefacto de la
   construcción del banco.
2. **El zinc se queda en el receptor.** hCA2 es una enzima de Zn y el sulfonamida se une
   al metal. Si el metal se cae al limpiar el receptor, todo lo medido carece de
   sentido. Es el mismo caso que el zinc de SOD1, que ya se resolvio.
3. **La sulfonamida de un anillo de mas o un tautomero cambia la union.** Con 13.000
   activos, la parte de la verdad de referencia que se prepare tiene que usar la misma
   receta de protonación y tautómeros para todos, o la comparación deja de ser justa.
   `preparar_ligando.py` es la que hay que usar, sin excepciones.

## El orden que yo propongo a partir de aquí

1. Auditar la serie de hCA2 con la receta de CDK2 (Murcko, MCS, solape). Si la serie
   resulta congenerica, el numero de activos independientes baja y hay que decirlo antes
   de acoplar.
2. Elegir el co-cristal de control con el mismo criterio pericial que se usó para `andr`
   (rigidez, sin polares sueltos, sin aguas, B bajo) y hacer el redocking. Ese control
   tiene que pasar **antes** de generar un solo señuelo, igual que se hizo con CDK2.
3. Solo entonces, banco de señuelos y cribado, con el IC saliendo solo del validador.

Esto es lo que hace al proyecto falsable por primera vez, y es tambien lo que
explicaría por qué las validaciones de TDP-43 y SOD1 no podían decir nada: no era la
técnica, era que se estaba midiendo con siete molecules.

## Ficheros

- `analysis/elegir_diana_validable.py` — la búsqueda por API y los dos filtros.
- `analysis/elegir_diana_validable.log` — la salida con los conteos y los perfiles.
- `analysis/dianas_validables.csv` — las doce candidatas con sus numeros.
