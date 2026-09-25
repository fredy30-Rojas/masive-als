# XL20 y XL23, la pareja del CR: comparados y acoplados (25 sep 2026)

**De dónde viene esto.** El barrido del PDB de ayer dejó escritas dos cosas: el PDB
está agotado para buscar positivos del sitio que usamos, y de la literatura salió
**XL20**, el primer unido medido del **CR del C-terminal (320-340)**, con la
dependencia de sitio medida (borrar el CR quita la unión; mutar el **Trp334** la
baja). De su misma tabla, **XL23 comparte el quimiotipo** y además inhibe la
agregación in vitro. Este informe hace las dos cosas que se pidieron: **comparar la
pareja** con números reproducibles y **construir el receptor del CR** para acoplarlos.

**El resultado, en una línea:** XL20 y XL23 comparten 19 de sus 26 y 34 átomos
pesados (la adenina, el aminociclohexanol y el enlace amida), y acoplados en la
hélice del CR **los dos caen en la misma cara del Trp334**; el que puntúa mejor en
crudo no es el que tiene la unión medida, sino XL23, en **las 18 combinaciones de
modelo y semilla** (42/2026/777) — y sin embargo **por átomo pesado gana XL20 en 17 de
esas mismas 18**. Las dos formas de leer el mismo acoplamiento dicen cosas contrarias,
y eso es lo que hay que contar: la diferencia es el tamaño del ligando, no el ajuste.
En este sitio el acoplamiento dice **dónde** se sientan, no **quién** se une.

---

## 1. La pareja, medida

`analysis/comparar_xl20_xl23.py` (nuevo) lee los dos SMILES de donde ya viven en el
repositorio (`controles_tdp43_xl20.csv` y `suplementario_tdp43_xl21_27.csv`) y
escribe la comparación en `analysis/_xl20_xl23/pareja.csv` y `pareja.txt`, con los
dos dibujos y el núcleo común resaltado en `XL20_XL23_comparacion.png`.

| | **XL20** (unión medida) | **XL23** (sin unión medida) |
|---|---|---|
| SMILES | `CN(Cc1ccccc1)[C@H]1CCC[C@H](n2cnc3c(N)ncnc32)[C@H]1O` | `Nc1ncnc2c1ncn2[C@H]1CCC[C@H](NC(=O)Cn2c(C(F)(F)F)nc3ccccc32)[C@@H]1O` |
| Fórmula | C19H24N6O | C21H21F3N8O2 |
| Masa | 352,4 Da | 474,4 Da |
| Átomos pesados | 26 | 34 |
| logP | 1,99 | 2,05 |
| TPSA | 93,1 Å² | 136,8 Å² |
| Donantes / aceptores | 2 / 6 | 3 / 7 |
| Enlaces rotables | 4 | 4 |
| Anillos (aromáticos) | 4 (3) | 5 (4) |
| Centroesteros definidos | 3 | 3 |

**Parecido 2D.** Tanimoto **ECFP4 = 0,444** (0,500 contando la estereoquímica), claves
**MACCS = 0,700**. Es el parecido de dos compuestos de la misma cabeza con colas
distintas, no el de un análogo cercano.

**Qué comparten exactamente.** El **núcleo común (MCS con anillos completos) son 19
átomos y 21 enlaces**, y contiene la **adenina** y el **aminociclohexanol**, con los
tres centroesteros alineados. Ese núcleo es casi todo XL20: 19 de sus 26 átomos
pesados. Lo que cada uno pone de su parte, cortando el núcleo:

- **XL20** añade `[1*]Cc1ccccc1` — el **bencilo sobre el nitrógeno** (7 átomos pesados).
- **XL23** añade `*Cn1c(C(F)(F)F)nc2ccccc21` (14 átomos pesados) más el **oxígeno del
  carbonilo** (1 átomo): un **benzimidazol con trifluorometilo colgando de una amida**,
  que en XL20 no existe.

**El matiz que conviene no saltarse:** el **esqueleto de Murcko NO coincide**. Es
decir, son la misma *cabeza* (misma química de unión, misma parte que se supone que se
apoya en el sitio) con colas que van en direcciones químicas distintas —bencilo
hidrófobo frente a benzimidazol-CF₃ polar, TPSA que sube de 93 a 137 Å²—. Lo que el
artículo llama «mismo quimiotipo» es esto, y así queda medido: 19 átomos comunes, dos
esqueletos distintos, 4 enlaces rotables en los dos.

## 2. El receptor del CR, elegido con números

`analysis/construir_receptor_cr.py` (nuevo) recorre las entradas del PDB que traen esa
región, pide a PDBe la secuencia de cada entidad y **localiza el fragmento dentro de
la TDP-43 humana (UniProt Q13148) por un ancla interna** (`MNFGAFSINP`), que es lo que
fija la numeración real de cada construcción. Comprobado antes de nada: en Q13148 el
**residuo 334 es triptófano** y el CR del artículo es
`PAMMAAAQAALQSSWGMMGML` (320-340). El barrido entero queda en
`analysis/_cr_receptor/candidatos_cr.csv`.

| PDB | Residuos | Modelos | Anillo del Trp334 | CR completo |
|---|---|---|---|---|
| **2N2C** | **307-349 (43)** | **6** | **sí** | **sí** |
| 2N3X / 2N4G / 2N4H | 311-360 (50) | 10 | sí | sí |
| 6N37 / 6N3A / 6N3B | 311-360 (50) | 1 | sí | sí |
| 7KWZ | 267-414 (148) | 1 | sí | sí |
| 7Q3U | 279-360 (82) | 1 | sí | sí |
| 9FOF / 9FOR | 282-345 / 284-345 | 1 | sí | sí (filamentos, con anexina A11) |
| 7PY2 / 8CG3 / 8CGG / 8CGH / 8QX9 / 8QXA / 8QXB | 1-414 (414) | 1 | sí | sí (TDP-43 entera) |
| 6N3C | — | — | — | se cae: su secuencia no lleva el ancla |

**Elegida: 2N2C.** Cumple los dos requisitos (cubre 320-340 entero y trae el anillo
indol resuelto) y es **la construcción más corta que lo hace**: 43 residuos frente a
50, 82, 148 o 414. Como el ligando solo puede tocar lo que hay, cuanta menos proteína
haya fuera del CR, menos se le puede atribuir a otra cosa lo que aparezca en las
poses. 2N2C es además la estructura **del propio CR** («TDP-43 prion-like hydrophobic
helix in DPC»), no una fibra ni la proteína entera. De las 19 estructuras probadas,
**18 cumplen** los dos requisitos y solo se cae 6N3C, cuya secuencia no lleva el ancla
(el recuento lo hace el script y queda en `caja.json`). Y el CR es de verdad una
hélice en ella: la distancia CA(i)–CA(i+4) en 320-340 promedia **6,59 Å** en los seis
modelos (de 6,49 a 6,64), que es lo que da una α-hélice.

**Los seis modelos, uno por uno.** El CR está intrínsecamente desordenado: acoplar en
un solo modelo sería elegir a dedo la conformación que a uno le guste. Se separan los
seis modelos de la RMN, se prepara cada uno con **la receta de acoplamiento del
proyecto** (meeko, `mk_prepare_receptor -p -j -a` — la misma con la que los controles
de redocking del Trp32 dan RMSD ≤ 2 Å) y **cada modelo lleva su propia caja**: el
centroide del anillo indol del Trp334 se separa hasta **4,15 Å** entre modelos
(`caja.json`), así que una caja compartida estaría descentrada en varios.

**No se pasa por PDBFixer, a propósito.** La receta validada del proyecto para acoplar
no lo usa (`redocking_trp32/preparar_receptor_limpio.py` documenta que reconstruía
cadenas laterales del bolsillo y las movía hasta 2,3 Å), y además sus hidrógenos no son
deterministas (medido en `RESCORING_WSL_CUDA_2026-09-11.md`).

**Qué ve la caja.** Caja de 22 Å (el lado del resto del proyecto) sobre el anillo del
Trp334: cubre **13-14 residuos**, entre el 326 y el 339 en casi todos los modelos — o
sea el CR y muy poco más.

## 3. El acoplamiento de la pareja

`analysis/acoplar_xl20_xl23_cr.py` (nuevo) prepara XL20 y XL23 con **la receta
canónica** (`preparar_ligando.escribir`, con la comprobación de pseudo-átomos de
pegamento), los acopla en los seis modelos con Vina (exhaustividad 8, 9 modos y
**tres semillas: 42, 2026 y 777**, las que usa el proyecto) y mide, además de la
afinidad, **los contactos a 4 Å**. Las dos cosas están puestas a propósito: una sola
semilla no permite decir ni una cosa ni la otra cuando las diferencias son de décimas
y Vina es estocástico, y con una caja que desborda una hélice de 43 residuos una pose
flotando en el hueco también puntúa, así que se registra cuántos átomos del ligando
tocan de verdad la proteína.

**Las dos caen en la cara del Trp334.** De 161 poses por compuesto (9 modos × 6 modelos
× 3 semillas; en dos combinaciones Vina devolvió 8 modos y no 9), **137 de XL20 y 140 de
XL23** tienen algún átomo a menos de 4 Å de un átomo del anillo indol del Trp334. Y no
es un sitio distinto: la mejor pose de cada uno toca el mismo tramo
`Q331-W334-G335-M337-G338` (XL20 añade P349; XL23 añade M339, A341 y P349). La misma
grieta.

| Modelo | XL20 (mediana) | XL23 (mediana) | Diferencia | Gana |
|---|---|---|---|---|
| 1 | −4,87 | −5,57 | −0,70 | XL23 |
| 2 | −5,08 | −5,65 | −0,61 | XL23 |
| 3 | −4,79 | −5,59 | −0,75 | XL23 |
| 4 | −5,35 | −6,58 | −1,20 | XL23 |
| 5 | −5,20 | −6,71 | −1,51 | XL23 |
| 6 | −4,56 | −5,30 | −0,73 | XL23 |
| **mejor de todos** | **−5,70** | **−6,92** | | **XL23** |

Mediana de la diferencia (XL23 − XL20) en las 18 combinaciones: **−0,80 kcal/mol**, y
va de −2,03 a −0,11. **XL23 gana en crudo las 18 de 18**, o sea que el orden no es
ruido de semilla: es el mismo en todos los modelos y con todas las semillas.

**Y aquí está lo que hay que leer despacio.** XL23 gana en crudo a favor de la
dirección contraria a la evidencia: XL20 es el único de la pareja con unión medida. La
explicación se puede comprobar, y se comprueba: **por átomo pesado gana XL20**, en 17
de esas mismas 18 combinaciones (mediana −0,190 frente a −0,166 kcal/mol por átomo;
buscando la eficiencia de la mejor pose de todos los modelos, −0,219 frente a −0,203).
La diferencia de 1,6 kcal/mol es del orden de la que Vina regala a un ligando con 8
átomos pesados más: en un bolsillo poco profundo como la cara de una hélice, el que más
abulta más contacta. **Es decir: dos formas igual de razonables de leer el mismo
acoplamiento dan ganadores distintos, y la que uno elija decide la historia. Ese no es
un resultado sobre la pareja: es un resultado sobre la función de puntuación**, y por
eso aquí el acoplamiento no se usa para ordenar.

## 4. Lo que se puede afirmar y lo que no

**Se puede afirmar, y está medido hoy:**

1. **Los dos compuestos son la misma cabeza con dos colas distintas**: 19 átomos
   comunes (adenina, aminociclohexanol y la amida), esqueletos de Murcko distintos,
   Tanimoto ECFP4 0,444. Ya no es una impresión de mirar la tabla del artículo.
2. **El CR tal como está depositado (2N2C) presenta una grieta en la cara del Trp334 y
   los dos compuestos caben en ella.** En los seis modelos de la RMN, con las tres
   semillas, y con los mismos residuos de contacto (Q331, W334, G335, M337, G338). El
   sitio que el artículo señala por mutagénesis es, en la estructura, un sitio que un
   acoplamiento encuentra solo.
3. **El receptor del CR existe ya en el repositorio** y es reproducible
   (`_cr_receptor/`): seis receptores PDBQT, su caja por modelo y las 322 poses con sus
   contactos.

**Lo que NO se puede afirmar, y por qué:**

1. **Que XL23 no se una.** No está medido: el artículo no publica unión para él, solo
   actividad funcional. Un acoplamiento no dice que algo no se una, y menos con 1,6
   kcal/mol de diferencia que son, por átomo, cero.
2. **Que XL20 se una mejor que XL23.** El puntaje en crudo dice lo contrario (XL23,
   18 de 18) y el puntaje por átomo pesado dice que sí (XL20, 17 de 18): cuando dos
   lecturas del mismo cálculo se contradicen, la conclusión honesta es que el cálculo
   no ordena esta pareja, no que una de las dos sea la buena.
3. **Nada sobre afinidad.** El propio artículo llama a su unión «micromolar aparente»
   y desvía del modelo 1:1; un número de Vina no es una Kd.
4. **Que esto sea una validación.** No lo es: **no hay control de redocking posible en
   el CR**, porque no existe ni un solo ligando co-cristalizado con TDP-43 en todo el
   PDB (barrido de ayer) y por tanto no hay pose experimental que reproducir. Sin ese
   control, el acoplamiento del CR es una **primera mirada**, no un protocolo
   validado — y así queda dicho, aunque la receta sea la misma que sí pasa el control
   en el bolsillo de RRM1-RRM2.

**Y una limitación de la estructura, dicha sin adornos:** 2N2C resuelve el CR como
hélice **dentro de una micela de DPC**, y la micela no está. El péptido desnudo es una
foto de una conformación de una región que en solución es desordenada; sirve para
decir dónde cabe un ligando y con qué contactos, no para cuantificar. La caja, además,
desborda el CR: contactos con 341, 347, 348 y 349 son de la hélice de fuera del CR, y
están anotados como tales en el CSV.

## 5. Y los otros seis de la tabla: el sitio no los distingue

Con la maquinaria de la pareja ya montada, la pregunta que sigue se contesta casi
sola: **los que sí tienen actividad funcional en el artículo (XL21 y XL23, que inhiben
la agregación del LCD in vitro a 100 µM) se colocan en el sitio del Trp334 de otra
manera que los que no tienen nada reportado (XL22, XL25, XL26)?**

`analysis/acoplar_familia_cr.py` (nuevo) prepara los siete compañeros de XL20 con la
misma receta canónica, los acopla con **el mismo receptor, la misma caja y las mismas
tres semillas**, y **reutiliza las poses de XL20 y XL23** que ya estaban calculadas: los
ocho quedan hechos exactamente con lo mismo. No reimplementa nada — lee poses,
contactos, distancias y la llamada a Vina del script de la pareja.

| Compuesto | Evidencia funcional (del artículo) | At. pesados | Mejor | Mediana | kcal/mol por átomo | Poses sobre Trp334 |
|---|---|---|---|---|---|---|
| XL20 | **unión medida** (SPR + CETSA) | 26 | −5,70 | −4,99 | −0,192 | 137 de 161 (85 %) |
| XL21 | inhibe la agregación (100 µM) | 23 | −5,38 | −4,93 | −0,215 | 71 de 162 (44 %) |
| XL23 | inhibe la agregación (100 µM) | 34 | **−6,92** | −5,64 | −0,166 | 140 de 161 (87 %) |
| XL22 | sin nada reportado | 22 | −5,31 | −4,84 | −0,220 | 80 de 158 (51 %) |
| XL25 | sin nada reportado | 29 | −6,55 | **−5,78** | −0,199 | 106 de 161 (66 %) |
| XL26 | sin nada reportado | 22 | −6,01 | −5,30 | **−0,241** | 104 de 162 (64 %) |
| XL24 | **empeora** la muerte neuronal | 20 | −5,00 | −4,47 | −0,224 | 58 de 162 (36 %) |
| XL27 | neuroprotege solo a 100 µM | 23 | −5,49 | −4,98 | −0,217 | 69 de 162 (43 %) |

Las medianas de cada combinación de modelo y semilla: inhibidores (XL21, XL23) de
−6,92 a −4,32, mediana −5,30; sin nada reportado (XL22, XL25, XL26) de −6,55 a −4,07,
mediana −5,05. **Se solapan por completo**: el peor de los inhibidores (−4,32) es peor
que el mejor de los que no tienen nada reportado (−6,55). El sitio no los separa.

**Y quitando el efecto del tamaño** (que es lo que Vina premia, y lo que hizo ganar a
XL23 en la pareja) pasa lo siguiente: la afinidad mediana sube −0,077 kcal/mol por
átomo pesado, y el residuo de cada uno —lo que le sobra o le falta respecto a lo que su
tamaño predice— queda así:

| Mejor de lo que le toca | | Peor de lo que le toca | |
|---|---|---|---|
| XL26 (nada reportado) | −0,405 | XL23 (inhibe) | +0,172 |
| XL25 (nada reportado) | −0,342 | **XL20 (unión medida)** | **+0,213** |
| XL27 (neuroprotege a 100 µM) | −0,011 | XL24 (empeora) | +0,273 |

Es decir: **el orden por residuo va al revés de la evidencia** —los dos mejor colocados
no tienen nada reportado y el único con unión medida es el segundo peor—, y el recorrido
entero del residuo (de −0,40 a +0,27, o sea 0,68 kcal/mol) es **del orden del ruido del
propio método**: el reparto de un mismo compuesto entre modelos y semillas tiene 0,39
kcal/mol de desviación típica. Así que ni el orden directo ni el inverso son una señal.
Esto **no** quiere decir que el CR prefiera a los inactivos; quiere decir que, con esta
índole, el sitio no está codificando la química que los separa.

Hay una lectura que sí es informativa y va en la misma dirección: la **fracción de
poses que se apoyan en el Trp334**. XL20 y XL23 la tienen altísima (85 % y 87 %) y
XL25/XL26 media (66 % y 64 %), pero XL21 —que también inhibe la agregación— se queda en
el 44 %, igual que XL27 (43 %) y XL24 (36 %). O sea que incluso la medida puramente
geométrica de «¿prefiere este sitio?» deja a los inhibidores en extremos opuestos. Y
conviene recordar que esos dos inhibidores **no se parecen entre sí** (uno es
uracilo-carboxamida y el otro adenina-aminociclohexanol): es lo esperable si la
inhibición de la agregación a 100 µM no es un dato de unión específica —en este campo
es conocido que buena parte de esos inhibidores lo son por mecanismos inespecíficos o
coloidales, y el artículo no publica unión para ninguno de los dos.

**La conclusión útil es negativa y ahorra trabajo:** el CR tal como está modelado (la
hélice del 2N2C) **no sirve para ordenar los compuestos de esta tabla**, ni por afinidad
ni por ocupación del sitio. Si algún día se quiere priorizar análogos con este receptor,
antes habrá que tener el dato que falta —unión medida de más de uno de ellos— y no al
revés.

## 6. Qué sigue

1. **La pareja solo se convierte en relación estructura-actividad cuando haya unión
   medida de XL23** (SPR o CETSA, como XL20). Mientras eso no exista, la pareja es una
   cabeza común con dos colas y nada más; y el dato que falta es un experimento, no
   más acoplamiento. Con los ocho de la tabla (§5) la conclusión es aún más clara: este
   receptor no ordena ninguno de los dos criterios, así que por esta vía no se puede
   elegir a quién mandar a medir.
2. **No se relanza nada por esta vía.** El CR no entra en la validación del bolsillo
   de RRM (XL20 sigue en `verdad_de_referencia.csv` como unido de otro sitio, marcado
   como no apto) y no hay GPU autorizada para una diana nueva sin criterio previo
   (`PLAN_TDP43_2026-09-20.md`).
3. **Si se quiere seguir la línea del CR**, lo que hay que decidir antes es el
   criterio: qué contaría como éxito (¿reproducir la dependencia del Trp334 con una
   serie de análogos sintéticos?) y con qué fondo de señuelos se compararía. Eso es una
   decisión de proyecto, no un script.

## 7. Cómo se reproduce

```
python C:/Users/Fredy/masive-als/analysis/comparar_xl20_xl23.py     # la pareja
python C:/Users/Fredy/masive-als/analysis/construir_receptor_cr.py  # elige el CR y arma el receptor
python C:/Users/Fredy/masive-als/analysis/acoplar_xl20_xl23_cr.py   # prepara y acopla la pareja, y lee contactos
python C:/Users/Fredy/masive-als/analysis/acoplar_familia_cr.py     # los ocho de la tabla, por grupos
```

- Pareja: `analysis/_xl20_xl23/pareja.csv`, `pareja.txt` y `XL20_XL23_comparacion.png`
  (los dos dibujos con el núcleo común resaltado).
- Receptor: `analysis/_cr_receptor/candidatos_cr.csv` (el barrido de estructuras del
  CR, con la razón de cada descarte), `pdb/2n2c.pdb`, `modelos/cr_modelo{1..6}.pdb` y
  `.pdbqt`, y `caja.json` (centro de caja, tamaño y residuos que ve, por modelo).
- Acoplamiento: `analysis/_cr_receptor/out/<ligando>_modelo<n>_s<semilla>.pdbqt` (las
  poses), `acoplamiento_cr.csv` (una fila por pose, con contactos, residuos, si toca el
  Trp334 y la distancia al anillo) y `acoplamiento_cr.txt` (el resumen legible: la
  mejor pose de cada combinación de modelo y semilla, la tabla por modelo y semilla, y
  el parteo por átomo pesado). Para repetirlo con otras semillas o más
exhaustividad está el `--semillas 42,2026,777` y el `--exhaustividad` del script.
- Los ocho de la tabla: `analysis/_cr_receptor/familia_cr.csv` (una fila por pose, con
  su grupo de evidencia) y `familia_cr.txt` (el resumen legible: la tabla por compuesto,
  la comparación por grupos, el residuo tras quitar el tamaño y el ruido del método).
  Comparte la carpeta de poses con la pareja, así que no repite trabajo.
- El reparto de trabajo con las herramientas que ya existían: la preparación de
  ligandos es `preparar_ligando.py` (la única copia de la receta) y la de receptores y
  el propio Vina son los de `redocking_trp32/redock_trp32.py`; aquí no se reimplementa
  ninguna de las dos.
