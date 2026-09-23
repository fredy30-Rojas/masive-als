# Control con Vina-GPU: el veredicto de quimias depende del motor

**23 de septiembre de 2026.** La validación limpia de TDP-43 va con **Vina de CPU**
(`tools/vina.exe`) porque el fondo duro tiene que salir del mismo motor que los
positivos y los señuelos. Este control hace lo otro: acopla **todo el conjunto**
(positivos, señuelos emparejados y fondo duro de R-BIND 2.0) con **Vina-GPU 2.1** de
una vez y pasa las mismas métricas. La pregunta era una: **¿el motor cambia el orden?**

**Sí lo cambia, y cambia el veredicto.** Con los señuelos emparejados la CPU dice
PASA (2 de 5 quimiotipos) y la GPU dice NO PASA (1 de 5). El cuadro grande sí se
sostiene: en los dos motores el AUC crudo queda **por debajo del azar**, el AUC por
átomo pesado anda entre 0,67 y 0,75, el enriquecimiento de cabeza es **cero**, y el
fondo duro de R-BIND **no ordena mejor** que los señuelos.

---

## 1. Lo que se comparó, y lo que NO es comparable

Los ligandos **son los mismos ficheros `.pdbqt`**: se copiaron los ya preparados con
la receta canónica (`preparar_ligando.py`), sin re-preparar nada. Mismo receptor
(`_tdp43_bolsillo_v2/4BS2_ph74.pdbqt`), misma caja de 26 Å, mismos positivos de
`verdad_de_referencia.csv`, el mismo criterio (quimiotipos que baten a los señuelos
**de su mismo tamaño**). Lo único que cambia es el motor.

**Y hay que decir la diferencia de perilla, porque no son la misma cosa:**
`search_depth` de Vina-GPU **no es** `exhaustiveness` de Vina. Aquí se usó
`search_depth = 20` (lo recomendado para 2.1) frente a `exhaustiveness = 16` en CPU.
Por eso este control juzga **el orden y los veredictos**, no los valores absolutos de
afinidad. Un detalle que ya se ve en la tabla: la pendiente de tamaño baja de
−0,039 a −0,022 kcal/mol por átomo pesado, y eso es parte de `search_depth`, no una
virtud del motor.

| Pieza | CPU | GPU |
|---|---|---|
| Motor | `tools/vina.exe` (AutoDock Vina 1.2.3) | `gpu_dock/Vina-GPU-2.1-win.exe` |
| Perilla | exhaustividad 16 | search_depth 20 |
| Máquina | CPU, hasta 6 procesos `--cpu 1` | RTX 4080 (8.128 MiB, 100 %) |
| Tiempo | 152 ligandos del fondo duro: **2 h 28 min** (19:23 → 21:51) | 281 ligandos: **10,7 min** |

---

## 2. La tabla, bloque a bloque

Mismos 7 positivos de unión medida, 122 señuelos emparejados, 152 unidores de ARN
de R-BIND 2.0.

| Bloque | Motor | AUC crudo | AUC/átomo | Residual | EF5 % | Quimias | Veredicto |
|---|---|---|---|---|---|---|---|
| Señuelos 122 | CPU | 0,301 | 0,686 | 0,306 | 0,00 | **2/5** | **PASA** |
| Señuelos 122 | **GPU** | 0,246 | 0,667 | 0,264 | 0,00 | **1/5** | **NO PASA** |
| Fondo duro 152 | CPU | 0,357 | 0,734 | 0,391 | 0,00 | 2/5 | PASA |
| Fondo duro 152 | **GPU** | 0,401 | 0,738 | 0,414 | 0,00 | **3/5** | PASA |
| Los dos 274 | CPU | 0,332 | 0,713 | 0,356 | 0,00 | 2/5 | PASA |
| Los dos 274 | **GPU** | 0,332 | 0,707 | 0,345 | 0,00 | **3/5** | PASA |

**Los 5 quimiotipos de los 7 positivos, y quién gana a los señuelos de su tamaño:**

| Quimiotipo | Positivo | CPU | GPU |
|---|---|---|---|
| bencilisoquinolina | berberrubine | sí | sí |
| isoxazol-piperidina | nTRD22 | sí | sí |
| piperidinil-pirimidina | rTRD01 | no | **sí** |
| piridil-pirazol | PE859 | no | no |
| fragmento | fragmento_1 / 2 / 3 | no | no |

**El que da la vuelta es rTRD01** (piperidinil-pirimidina): con la CPU no cruza el
corte de su tamaño, con la GPU sí. Los tres fragmentos van todos juntos en un solo
quimiotipo («fragmento»), así que el total son 5 y no 7.

Los dos que ganan en las dos corridas son **exactamente los dos que el informe del
23 de septiembre coloca, según la literatura, fuera de ese bolsillo** (nTRD22 es
modulador alostérico del dominio N-terminal, la berberrubina se describe en la
interfaz RRM1–RRM2). El motor no cambia eso: cambia cuántos cruzan el corte.

**Lo que se mueve y lo que no.** Se mueve el veredicto del bloque blando (2/5 → 1/5,
PASA → NO PASA), una quimia en el fondo duro, y el AUC crudo del fondo duro
(0,357 → 0,401). No se mueve la conclusión de fondo: **ninguna de las dos corridas
ordena en crudo** —las seis AUC están por debajo de 0,5—, el enriquecimiento de
cabeza es cero en las seis, y el fondo duro no separa mejor que los señuelos.
Curiosidad que conviene no confundir con un hallazgo: el bloque de los dos fondos
juntos da **0,332 con los dos motores**, con residuales distintos; es casualidad de
dos decimales, no una medida de acuerdo.

---

## 3. Las dos poses que faltan, dichas enteras

1. **`RB_SM_0085` (DB1273), del fondo duro: no se puede preparar.** Sus dos anillos
   son de **selenio** y Meeko no sabe tipar el Se: `escritura PDBQT: atom number 24
   has None type`. Comprobado que el mismo esqueleto con azufre prepara sin problema
   (28 átomos pesados), pero **no se le cambia el Se por S a la callada**. Entra
   152 de 153, y el que falta se declara.
2. **`DEC_CHEMBL4543460`, señuelo: tiene un ácido borónico y Vina no acepta el boro.**
   El motor de la GPU no escribió su pose y parecía un fallo del control; no lo era.
   Probado a mano con la CPU (ver §7), Vina dice: `PDBQT parsing error: Atom type B
   is not a valid AutoDock type`. **Ninguno de los dos motores puede acoplarlo**, y
   de hecho la corrida de CPU tampoco lo tenía: el bloque de señuelos es de 122 por
   esto, en las dos columnas, y no porque la GPU se quedara corta.

   Se revisó el resto del conjunto con la misma lupa (tipos de átomo fuera de la
   lista de AutoDock en los 132 ficheros del conjunto viejo y en los 152 del fondo
   duro): **solo este**. Un descarte silencioso, no un problema de fondo.

Además, el montaje del control lleva **14 ficheros `ACT_`** (los 9 del conjunto viejo
más los 5 que prepara la corrida de CPU), pero el criterio solo usa **los 7 aptos**
de `verdad_de_referencia.csv`: los otros son positivos de otros conjuntos que viajan
en la carpeta y el script no puntúa.

Los 7 aptos, para que quede claro quién es quién: `rTRD01` (piperidinil-pirimidina),
`nTRD22` (isoxazol-piperidina), `fragmento_1`, `fragmento_2`, `fragmento_3`
(fragmento), `PE859` (piridil-pirazol) y `berberrubine` (bencilisoquinolina).

---

## 4. Una trampa que hay que dejar escrita

La primera pasada del control se quedó **sin los 5 positivos nuevos** (`rTRD01`,
`nTRD22` y los tres fragmentos), porque sus ficheros los prepara la corrida de CPU
(`validar_tdp43_limpia/ligands`) y el montaje se hizo antes. El pre-flight de la
validación los detectó y **los re-acopló con CPU dentro del control**, así que esa
pasada mezclaba motores sin decirlo: dio fondo duro 0,404 / 0,747 / 0,411 con 3 de 5.

Ya está corregido —`control_gpu_tdp43.py` copia también de `validar_tdp43_limpia/ligands`
y los acopla con la GPU— y los números de §2 son de la pasada **toda en GPU**
(289 ligandos en el montaje, 288 con pose, 7 de 7 positivos).

---

## 5. Lo que NO dice

1. **No se publica nada con esto.** Es un control interno de método, con 7 positivos.
2. **No dice que la GPU sea más exacta ni menos.** Dice que entre los dos motores
   **cambia el veredicto**, y que por eso el "PASA por emparejado de tamaño" no se
   puede reportar como si fuera una propiedad de la quimia.
3. **No separa las dos causas.** El cambio puede venir del motor (implementación
   distinta) o de la perilla (`search_depth` 20 frente a exhaustividad 16). Para
   separarlas haría falta un barrido de exhaustividad **dentro de la CPU** (por
   ejemplo 8 frente a 32) con el mismo conjunto: si el veredicto también se mueve,
   entonces lo que no aguanta es el criterio, no el motor.
4. **No arregla el fondo duro como fondo duro:** sigue por debajo del azar en crudo
   (0,357 y 0,401) con los dos motores. Que los señuelos emparejados y los verdaderos
   unidores de ARN den casi lo mismo es el resultado incómodo, no un fallo del control.
5. **No dice nada del descarte por boro de otros conjuntos.** Aquí se revisaron los
   dos que alimentan esta validación y solo había uno; el cribado viejo no se ha
   revisado con esta lupa.

---

## 6. Qué significa para el proyecto

- **El criterio de decisión no es estable.** "Al menos dos quimiotipos por delante de
  los de su tamaño" da 1, 2 o 3 según el motor y según el fondo. Es el número que el
  proyecto usa para decidir si una idea sigue o se tira, y se mueve solo con cambiar
  de motor. Eso hay que arreglarlo antes de volver a decidir nada con él.
- **Lo que sí aguanta en las dos corridas, y por tanto es lo que se puede decir:**
  la energía acoplada **no ordena en crudo** (AUC 0,25–0,40, por debajo del azar), el
  corte fijo de cabeza **no sirve** (EF 0 % en los seis bloques), el AUC por átomo
  pesado es el único número por encima del azar (0,67–0,75), y el fondo de verdaderos
  unidores de ARN **no se separa** del de señuelos. Nada de esto depende del motor.
- **Consecuencia práctica:** la lista de candidatos **no se toca** por esto, y la
  validación de TDP-43 **sigue sin reportarse**. El paso que toca es el barrido de
  exhaustividad dentro de la CPU (§5.3), que es barato y decide si el problema es el
  criterio.

---

## 7. Anexo: cómo se reprodujo

```
python C:/Users/Fredy/masive-als/analysis/control_gpu_tdp43.py --todo
python C:/Users/Fredy/masive-als/analysis/validar_diana_limpia.py --target TDP43
```

Lanzadores: `analysis/control_gpu_tdp43.bat` y
`analysis/lanzar_validar_tdp43_limpia.bat` (los dos por
`wscript.exe ejecutar_bat_oculto.vbs`, que es como corren sin ventana).

- Montaje y poses: `_control_gpu_TDP43/ligands` (289), `out` (positivos y señuelos),
  `out_duro` (fondo duro R-BIND).
- Registros: `_control_gpu_TDP43/acoplar_gpu.log`, `validar_gpu.log`,
  `validar_gpu.csv`, `validar_gpu_resumen.csv`; la salida cruda del motor, en
  `_control_gpu_TDP43/vina_gpu_stdout.log`.
- La corrida de CPU: `validar_tdp43_limpia.log`, `validar_tdp43_limpia.csv`,
  `validar_tdp43_limpia_resumen.csv`; el fondo duro que la alimenta se acopló con
  `analysis/acoplar_fondo_rbind.py` (152/153).

La comprobación del boro, a mano:

```
tools/vina.exe --receptor analysis/_tdp43_bolsillo_v2/4BS2_ph74.pdbqt \
  --ligand analysis/_validacion_TDP43/ligands/DEC_CHEMBL4543460.pdbqt --cpu 1 ...
PDBQT parsing error: Atom type B is not a valid AutoDock type (atom types are case-sensitive).
```
