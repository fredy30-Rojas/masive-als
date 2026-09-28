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
