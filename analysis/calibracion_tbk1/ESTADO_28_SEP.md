# TBK1: qué estaba parado y qué se arregló — 28 sep 2026

## 1. El docking se había parado en silencio

Arrancó solo la madrugada del 28 a las 00:52 (lo lanzó el vigilante cuando vio
240 s de quietud en la preparación) y el log dejó de escribir a las 04:00. A las
10:00 la RTX 4080 estaba al 0% y no había ningún proceso de Vina. De los 35.349
ligandos solo había 300 poses, las del lote de prueba del 27.

Dos fallos de fondo, ninguno era la GPU:

### a) El lote de 5.000 se perdía entero
`lanzar_banco_tbk1.py` copiaba las poses del lote a la carpeta buena **solo
cuando el lote terminaba**. El lote de 5.000 murió a las tres horas y se llevó
por delante las tres horas de trabajo.

**Cambio:** las poses se vuelcan en vivo cada 15 s, y el lote baja de 5.000 a
500. Una muerte ahora cuesta como mucho el ligando que esté en curso, unos
22 minutos de trabajo en vez de tres horas.

### b) 142 ligandos que Vina no puede leer
El log estaba lleno de `ATOM syntax incorrect: "B" is not a valid AutoDock type`
y también de `"Si"`. Son ligandos con boro, silicio o sodio, tipos que AutoDock 4
no tiene: se perdían en silencio dentro del lote.

**Cambio:** `lanzar_banco_tbk1.py` los detecta antes de acoplar y los salta. El
filtro tarda 2 s para los 35.349 y es exacto:

    descartados 142: tipo de atomo que AutoDock no entiende (B=104, Na=2, Si=38)

Los nombres quedan en `descartados_tbk1.txt` para poder contarlos después.

## 2. El MM-GBSA: cinco fallos, todos reales, todos medidos en Oracle

El paso 4 (tleap) llevaba noches atascado y el diagnóstico que había ("le faltan
4 términos de torsión a GAFF") era solo el último. En total:

1. **Faltaba `source leaprc.water.tip3p`.** Sin él, `TIP3PBOX` no existe como
   argumento y tleap lo leía como texto: *"Argument #2 is of type String must be
   of type: [unit]"*.
2. **`solvateBox` con 4 argumentos.** Esta versión de AmberTools pide 3:
   `solvateBox <solute> <solvent> <buffer> [iso] [closeness]`. Con 4 decía
   *"usage: solvateBox..."*.
3. **`addions` no existe**, es `addIons` con mayúscula. Y no hacía falta: el
   complejo ya sale neutro. Se quitó.
4. **`saveoff3` y `saveamberpdb` no existen** en AmberTools moderno. El
   `desc` los devuelve como `STRING (with no reference)`, o sea que tleap los
   trata como texto. El MM-GBSA solo necesita `.prmtop` y `.rst7`, así que se
   quitaron; el sustituto de `saveoff3` es `saveoff`, con 2 argumentos.
5. **Los 4 términos de torsión de GAFF sí eran un problema**, pero no de donde
   se pensaba. `antechamber` **no** genera el `frcmod` para esas combinaciones
   (`ca-cp-nb-ca`, `cp-cp-nb-ca`): no escribe ninguno. Lo que sí las rellena es
   **`parmchk2 -s gaff`**, llamado aparte, que las cubre como
   *"same as X -ca-nb-X"*. Con su `frcmod` cargado: **0 torsiones faltantes** y el
   complejo pasa a 17 MB de prmtop (antes salía de 0 bytes).

### Y el paso 5, que no se había llegado ni a tocar

`sander` a mano tampoco funciona en AmberTools moderno: `-O` ya no existe
(*"mdfil: Error unknown flag"*) y la carpeta `MDOUT` con las plantillas no viene
en la instalación de conda-forge. Se cambió a **`MMPBSA.py`**, que es la
herramienta para esto y genera la entrada de sander sola. Detalles que costaron
tiempo y quedan en el código:

- `igb`, `saltcon` van en `&gb`, no en `&general` (*"Unknown variable igb"*).
- `prmtop` **no** se puede pasar en la entrada ni como argumento (*"Unknown
  variable prmtop in &general"*). MMPBSA.py espera los ficheros con el nombre
  literal `complex_prmtop` y `complex_inpcrd`, sin ningún sufijo.
- También espera un `mdcrd` con las coordenadas, que no es la inpcrd.
- `surften`/`surftat` **no existen** en esta versión.

### El error de fondo del paso 5: no era de sintaxis

Con todo lo anterior, `mmpbsa_py_energy` fallaba con:

    Error: gb>0 is incompatible with periodic boundary conditions.
    To use this method set IFBOX in the PRMTOP file to 0

Es un problema físico, no de código: `solvateBox` deja el prmtop con `IFBOX=1` y
el Generalized Born es, por definición, un modelo de solvente **no** periódico.
Quitar la caja no es un truco, es lo que el método pide.

**Cambio:** `quitar_caja()` en el script, con **parmed** (que viene con
AmberTools) en vez de editar el prmtop a mano. Editar el flag a mano deja los
contadores de la sección BOX desfasados y sander lee memoria de más.

## Cómo queda

- El docking corre en la RTX 4080 del PC, lotes de 500, volcado en vivo.
  Ritmo medido 2,25 s por ligando; con 34.867 pendientes son unas **22 horas**.
- El MM-GBSA está en Oracle, en `~/prueba_mmgbsa/`, sobre una pose real de Vina
  (ACT_CHEMBL1078178). Al cerrarlo la sesión, el cálculo de GB con 90.570 átomos
  dejó la VM sin responder por RAM: es pesado, pero arranca.

---

# Cierre del 28 de septiembre: el MM-GBSA ya da dG

## El fallo de fondo: el emparejamiento por orden

Durante días se tipó el ligando desde el SMILES con antechamber y luego se
metieron las coordenadas de la pose **emparejando uno con uno por orden**.
Ese emparejamiento era falso: **antechamber reordena los átomos**.

Medido: en el mol2 tipado el C4 y el C9 están a cuatro enlaces de distancia;
al meter las coordenadas de la pose salen a 1,489 Å, que es un enlace y medio.
El C de la pose se quedaba con el tipo y la carga del C equivocado. Todos los
síntomas posteriores (VDWAALS nan, EEL inf, H a 0,5 Å del vecino, tipos de H
que no existían, ángulos imposibles) eran de ese único error, no fallos
aparte.

**La solución no fue arreglar el emparejamiento sino invertir el orden**, para
que sea AmberTools quien reordene sobre unas coordenadas que ya son las suyas:

1. `obabel -ipdbqt` lee el pdbqt de la pose, que **sí** trae conectividad
   (Vina escribe las ramas ROOT/BRANCH). 33 átomos, 36 enlaces.
2. Se le ponen los H que falten, contando por **valencia** (no a ojo).
   Pide 20, y salen 53 átomos y 56 enlaces: exactamente los mismos recuentos
   que daba antechamber desde el SMILES. Que cuadre la cuenta por dos caminos
   distintos es lo que hizo confiar en el método.
3. Antechamber + parmchk2 lo tipan y escriben un mol2 que **ya lleva la pose
   dentro**. De regalo, GAFF se asigna por contexto real y no por elemento.

El script ya **no necesita `--smiles`**: la topología sale de los enlaces que
el propio Vina escribió.

## Colocación de los hidrógenos

AmberTools 21 (conda-forge) no sabe ponerlos, y está todo medido en
`anadir_hidrogenos`: `addHydrogens` y `hadd` no existen, `addH` los crea sin
tipo y deja el prmtop de 0 bytes, `pdb4amber` y `reduce` no hacen nada con un
residuo UNK, y `obabel -h` sobre un mol2 GAFF añade 1 H de 32.

Se ponen con distancia de enlace estándar, en la dirección opuesta a la suma
de los vecinos, en cono de 2 o 3 si hacen falta varios. El criterio de elección
del sitio es **uno solo**: la distancia mínima a todo lo ya colocado (el
ligando y el receptor). Antes había dos criterios que se peleaban y dejaban
los dos H de dos metilenos contiguos a 0,83 Å.

Los H también esquivan al **receptor**: sin eso un H se metía dentro de la
proteína (el H5 a 0,269 Å de un LEU15) y el vdW salía 1,98e15.

## El agua: fuera, y por qué

Se probó `solvateBox` (90.570 átomos, el MM-GBSA se colgaba; no era RAM, era
tamaño) y `solvateOct` (rápido, pero...). El problema de `solvateOct` es que
pone agua alrededor de lo que le digas **sin preguntar si hay proteína
alrededor**: metía 5.685 moléculas de agua dentro del receptor, con pares a
0,28 Å (OG 409 contra 14532 a 0,684 Å).

Y el agua no hacía falta: con `igb=5` el solvente es **implícito**. El agua
explícita no aporta al dG y, al estar en el sistema, GB la trata como soluto.
Un MM-GBSA con GB no lleva agua.

**Sistema ahora:** receptor 9.688 + ligando 53 = complejo 9.741 átomos.

## Dos comprobaciones que convierten el silencio en error

- `comprobar_solapes_sistema()` sobre el `complex.rst7` antes de gastar un
  minuto en el MM-GBSA, con rejilla para que tarde menos de un segundo. Cuando
  dos átomos caen en el mismo sitio, el motor **no protesta**: devuelve
  1,98e15 y sigue. Ahora sale el par y la distancia.
- `leer_energia_mmpbsa()` se reescribió con el formato real de AmberTools 14.
  No existe ninguna línea "DG (Energy terms)"; las filas se llaman `VDWAALS`,
  `EEL`, `EGB`, `ESURF`, `G gas`, `G solv`, `TOTAL`. El lector viejo buscaba
  un encabezado que esta versión no escribe, devolvía `None` sin decir nada y
  el paso 3 reportaba "FALLO" con las energías ya calculadas.

## Resultado

```
3. MM-GBSA complex   G =   -12090.1 kcal/mol
3. MM-GBSA receptor  G =   -17023.0 kcal/mol
3. MM-GBSA ligand    G =     4947.1 kcal/mol
dG_union MM-GBSA = -14.3 kcal/mol
```

Ligando de prueba ACT_CHEMBL1078178, afinidad de Vina -10,2. Tarda 12 s por
sistema. **El MM-GBSA del proyecto funciona.**

## Docking

Va por 13.837 de 35.207 poses, 3,19 s por ligando, unas 19 horas left.

---

# Cribado de 20 y de 300: el MM-GBSA da numero pero no se sabe si ordena

## Lo que se midio con 20 ligandos (10 activos, 10 inactivos)

19 de 20 salieron con dG (el vigesimo fallaba, ver abajo). 36 s por ligando
con los 4 nucleos de Oracle, en vez de los 40 s en serie: doce veces mas
rapido por paralelizar, que es abrir un proceso por nucleo y ya.

| clase | dG min | mediana | dG max |
|---|---|---|---|
| activos (9) | -45,6 | -10,5 | +104,9 |
| inactivos (10) | -52,2 | -12,9 | +255,5 |

- **AUROC 0,567.** Cero coma cincuenta es azar. Con nueve y diez ligandos el
  margen de error es de unas dos décimas, así que con veinte ligandos NO SE
  PUEDE CONCLUIR. Solo queda dicho que el número sale.
- Spearman contra pChEMBL: -0,095.
- El rango de 300 kcal/mol es la pista: eso no es señal, es ruido.

## El fallo que quedaba: dos H a 0,532 A

CHEMBL5758899, en la colocación de H. La poda del bucle de candidatos
comparaba **ángstroms al cuadrado contra ángstroms** (`puntuacion <=
mejor_puntuacion`, con `mejor_puntuacion` ya en Å). Eso no descartaba el
candidato que debía, y se quedaba con uno cuyo mínimo real era 0,25 Å.
Corregido a `mejor_puntuacion ** 2`. Con eso: **20 de 20**.

## La minimizacion: se probo y NO sirve (queda documentado, apagada)

Relaxar la geometria antes de calcular es la via correcta en principio: se
calculaba sobre la pose de Vina sin tocar un atomo, con los H puestos por una
regla de geometria ideal, y de ahi el ruido. Se anadio `minimizar()` con
sander, 200 pasos. Y al medir:

- dG_union pasa de **-14,3 a -30,1 kcal/mol**. En 9.741 atomos eso no es
  fisico.
- **vdW del complejo = -5.420.** Una proteina minimizada da +30 a +80. Un vdW
  muy negativo significa que los atomos se han amontonado: la minimizacion
  COLAPSA en vez de relajar.
- sander avisa "Maximum number of minimization cycles reached": no converge.
- 36 s -> **6 min y medio** por ligando.
- El rst7 que escribe sander en Amber 26 es BINARIO (Fortran nuevo, con
  "spatial" y "atom" dentro). MMPBSA.py lo lee, pero cualquier script que lo
  abra a pelo revienta con ValueError. Para volver a texto hace falta parmed.

Queda `--minimizar N` disponible y **desactivado por defecto (0)**. Lo
correcto seria minimizar SOLO el ligando (53 atomos, converge en segundos) e
inyectar esas coordenadas en el complejo, o minimizar con el receptor
restringido con `nmropt`. Pendiente de medir.

## Un fallo de sander que costo tiempo: la linea de titulo

Un mdin que empieza directamente por `&cntrl` hace que sander responda
"Could not find cntrl namelist" y salga con codigo 1 sin calcular nada. Se
come la primera linea como TITULO y luego busca el namelist a partir de la
segunda. Se probaron ocho variantes (`&end` en vez de `/`, el bloque entero en
la linea 1, dos y seis espacios de sangria, `mdin` en vez de `-i`, con
AMBERHOME, con `amber.sh` cargado) y las ocho fallan igual. **Con una linea de
texto delante, rc=0.** El mdout lo delata: "Here is the input file:" salia
VACIO. El MM-GBSA usa MMPBSA.py, que genera su propio mdin, asi que este
fallo solo aparece en la minimizacion.

## Muestra grande en marcha

`preparar_muestra.py` toma 50 ligandos de CADA una de las cuatro franjas de
activos (no al azar: al azar casi todos caen en la mejor y no hay rango con el
que trabajar) y 100 de los 141 inactivos. Son 300, con las 300 poses ya
empaquetadas. Con 300 el margen del AUROC baja a unas cinco centesimas, que
si es concluyente.

---

# El diagnóstico: el problema no era el método, era la lista de referencia

Este es el hallazgo del 29 de septiembre y cambia la lectura de todo lo anterior.

## El MM-GBSA no aporta nada sobre Vina (medido sobre 277 ligandos con dG)

| criterio de verdad | Vina | MM-GBSA |
|---|---|---|
| lista del banco (mejor valor) | 0,602 ± 0,034 | 0,524 ± 0,036 |
| verdad v2 (mediana, 2+ medidas) | 0,617 ± 0,033 | 0,692 ± 0,097 |

El mejor peso que se le puede dar al MM-GBSA sobre la afinidad de Vina es
**CERO**. Cero coma cincuenta y dos contra cero coma cincuenta es azar, con
un margen de 0,036: no es "no demostrado", es medido como nulo.
**Veredicto: el rescoring con MM-GBSA no se usa más.** El código se queda
documentado, con los seis fallos que hicieron falta para que funcionara, para
no repetir el camino.

## Pero el docking no es malo: la lista lo estaba engañando

El control de redocking del BX-795 da **RMSD 1,23 Å** (listón del proyecto: 2 Å).
Receptor y caja están bien. Y sin embargo:

| verdad de referencia | Vina |
|---|---|
| el mejor pChEMBL de cada compuesto (la del banco) | **0,535** |
| mediana de las medidas, 2+ medidas | **0,617 ± 0,033** |
| mediana, 3+ medidas | **0,730 ± 0,095** |

Y desglosado por cuántas medidas tiene cada compuesto:

| medidas por compuesto | activos | inactivos | AUROC de Vina |
|---|---|---|---|
| 1 sola | 1156 | 703 | **0,515** (azar) |
| 2 o más | 435 | 17 | **0,590** |

**Cuantos más datos hay de un compuesto, mejor predice.** Con un solo dato
publicado, el docking no dice absolutamente nada. Eso era lo que ahogaba la
señal.

## Por qué la lista del banco estaba contaminada (medido)

1. **"El mejor valor de cada compuesto".** Con la mediana, **620 de los 2.172
   activos (29%) dejan de serlo**. Elegir el dato más favorable de todos es
   elegir siempre el que favorece la hipótesis.
2. **El 80% de los compuestos tienen UNA sola medición**: 1863 de 2315.
3. **1.221 medidas son límites** (`<` o `>`), no medidas: 831 `>` y 388 `<`.
   Contarlas como exactas infla los inactivos.
4. **518 actividades las marca ChEMBL como dudosas** ("los valores parecen un
   orden de magnitud distintos de los publicados, las unidades pueden estar
   mal"). Tirarlas está bien hecho; el resto de las perdidas sí se recuperan.
5. **6.056 actividades cinéticas (3.028 kon + 3.028 k_off) están VACÍAS**:
   tienen `value` y `standard_value` a `None`. ChEMBL guardó el registro del
   ensayo y no el número. Es el 46% de las actividades de la diana y no
   aporta ni una cifra.

De 13.123 actividades crudas solo **2.839 tienen pChEMBL utilizable**, y tras
filtrar quedan 393 activos y 77 inactivos con dos o más medidas de unión
(`verdad_tbk1_v2.csv`).

## Lo que queda

Con 77 inactivos el margen del AUROC es de ±0,04: se puede medir, pero no hay
margin para más. **TBK1 no tiene datos suficientes para validar un método de
cribado.** Para validar el pipeline de punta a punta hace falta un conjunto
público (DUD-E, CASF-2016 o equivalente). Para descubrir compuestos con un
AUROC de 0,62 no hay método: el consenso multi-diana (network pharmacology) es
lo que aguanta esa señal.

## 29 sep: veredicto de los tres rescorings (medido, 470 ligandos)

La verdad de referencia v2 (393 activos y 77 inactivos con dos o mas medidas,
IC50/Ki/EC50 en nM, sin limites) se rescoreo tres veces con las mismas poses:

- Vina, mejor pose: AUROC 0,617. Con la media de las 9 poses: 0,637.
- MM-GBSA (igb=5, sin agua, sin minimizar): AUROC 0,498. No discrimina nada.
- GNINA 1.3.3 CNN (score_only, pose 1, en Kaggle): CNNscore 0,536,
  CNNaffinity 0,547, CNN_VS 0,545. La fusion por rangos con Vina no suma
  (0,613-0,617): el CNN y Vina se equivocan en los mismos compuestos.

Conclusion: el techo del docking clasico en TBK1 con este protocolo esta en
0,62-0,64. El MM-GBSA tal cual esta (pose congelada, H por geometria ideal,
sin minimizacion del complejo) no aporta senal y no vale como rescore aqui.
La senal de Vina es real (la normalizacion por tamano no la reproduce) y
mejora poco con las 9 poses.

Fallos que dejo resueltos en el camino del kernel GNINA (binario x86 en
Kaggle, no en Oracle que es ARM): faltaban las librerias CUDA 12 de pip
(cudnn, cudart, cusparse), las lineas MODEL/ENDMDL hay que quitarlas al
partir un pdbqt de 9 modelos (smina da Parse error si las copia), y el
CLI de kaggle en Windows necesita PYTHONUTF8=1 para bajar el log.

La criba MM-GBSA de los 470 dejo 435 ok y 35 fallos, 7 dG absurdos
(|dG|>500, choques) filtrados por umbral.
