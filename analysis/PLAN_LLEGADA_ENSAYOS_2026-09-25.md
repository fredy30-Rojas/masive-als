# Plan para cuando lleguen los ensayos (25 sep 2026)

**Qué es esto.** La petición de ensayos está escrita (`PETICION_ENSAYOS_2026-09-25.md`) y
lo que falta está contado con números. Esto es lo que hay que hacer el día que llegue un
resultado, **escrito antes de verlo**, que es la única forma de que valga: si las reglas
se eligen después de mirar el dato, lo que se mide es la elección, no el embudo.

---

## 1. El veredicto de hoy queda congelado

El estado con el que se compara todo lo que venga está copiado en
`analysis/regla_decision/_antes_2026-09-25/`: `regla_decision.csv`,
`sod1_fondos_separados.csv` y el log de la corrida, guardado como
`regla_decision_log_congelado.txt` —el original se llama `.log` y el `.gitignore` del
repositorio ignora esa extensión, así que con ese nombre no se habría subido—. La regla
**sobrescribe** sus salidas cada vez que corre, así que sin esa copia el veredicto de hoy
se perdería al primer re-run.

De hoy, lo que hay que batir: **TDP-43 PASA contra el fondo duro** (0,731–0,738; límite
inferior 0,686–0,694 en las cuatro corridas) y **SIN EVIDENCIA** en señuelos emparejados
(0,666–0,686; límite inferior 0,434–0,465). **SOD1 SIN EVIDENCIA en sus tres bloques**,
con el 0,530 contra señuelos emparejados y el 0,643 contra su fondo duro.

## 2. Lo que está fijado por adelantado y no se toca al ver el resultado

Umbral 0,5. Semilla 20260924. 2.000 remuestreos. Los dos fondos, tal como están hoy (122
señuelos emparejados y 152 unidores de ARN en TDP-43; 339 y 143 en SOD1). La condición de
validez de la regla: **el veredicto tiene que ser el mismo en todas las corridas**, y una
corrida solo se puede caer por una razón técnica declarada por escrito *antes* de mirar el
veredicto (por ejemplo, que la GPU no esté disponible).

Y tres reglas de entrada que no admiten excepción:

1. **El positivo entra completo.** No se saca de la lista porque puntúe mal. Un positivo
   que el embudo no recupera es precisamente el dato que mide el embudo.
2. **El positivo entra con su sitio.** Un unido medido de otro bolsillo no entra en el
   bloque de la caja que no le toca (es el caso de XL20, que es del CR y no del RRM, y así
   está marcado en la verdad de referencia).
3. **Un análogo suma `n`, no independencia.** Si el positivo nuevo es de la misma familia
   que otro que ya está, se declara como congénere: la regla crece con positivos
   independientes y con análogos crece menos de lo que parece. Esto no veta nada, se
   anota.

## 3. El portero: `analysis/incorporar_positivo.py`

Nadie mete una fila a mano en `verdad_de_referencia.csv`. El portero comprueba la entrada
y **no escribe** en la verdad: deja la fila preparada en `analysis/positivos_pendientes.csv`
con su veredicto y su motivo, y una persona la revisa y la mueve.

Comprueba, en este orden: que el ensayo sea de **unión medida** (Kd, SPR, ITC, MST, CSP,
RMN, CETSA) y no actividad celular o de agregación; que el **sitio declarado sea el de la
caja** de esa diana; que el **SMILES** lo construya RDKit y cuántos átomos pesados tiene; y
la **independencia** por Tanimoto ECFP4 y esqueleto de Murcko contra los positivos que ya
están. Devuelve `entra`, `entra_congenere` o `no_entra`, y si entra imprime los comandos
exactos que toca correr.

Probado hoy con tres casos: XL23 declarado contra el RRM (`no_entra`, su sitio es el CR),
CHEMBL2165613 con un Kd por SPR (`entra`, química nueva) y PRG-A01 con un EC50 celular
(`no_entra`). Las filas de prueba no se han dejado en el repositorio.

## 4. Los tres caminos posibles

### A. Llega un positivo nuevo del sitio correcto

Es el caso que se pidió. La secuencia, y en este orden:

```
# 1. el portero (deja la fila preparada y dice si vale)
python analysis/incorporar_positivo.py --target TDP43 --ligando ... --smiles ... \
    --sitio "RRM1" --ensayo "SPR Kd ..." --cita "..." --quimia "..."
# 2. revisar la fila y moverla a verdad_de_referencia.csv  (a mano, con su cita)
# 3. guardar el veredicto de hoy
mkdir -p analysis/regla_decision/_antes_AAAA-MM-DD
cp analysis/regla_decision/regla_decision.* analysis/regla_decision/_antes_AAAA-MM-DD/
# 4. re-puntuar: el validador re-prepara y re-acopla SOLO lo que no cuadre
python analysis/validar_diana_limpia.py --target TDP43
python analysis/barrido_exhaustividad_tdp43.py        # exh 8 y 32
python analysis/control_gpu_tdp43.py                  # la corrida de GPU
# 5. el veredicto nuevo
python analysis/regla_decision_bootstrap.py
```

> **Aviso del 25 sep 2026, antes de re-puntuar TDP-43.** Al montar el fondo de
> inactividad se corrió el validador de TDP-43 para comprobar que no se rompía nada, y
> salieron dos cosas que hay que tener delante el día del re-run. Primera: los ligandos
> preparados de `rTRD01` y `nTRD22` **ya no están** en la carpeta del fondo
> (`analysis/_validacion_TDP43/ligands` está a medias en disco), así que la comprobación
> previa los da por sospechosos y **los re-acopla**: con la receta de hoy salen −6,624 y
> −8,228 en vez de los −6,621 y −8,026 publicados. Las poses re-acopladas **no se han
> dejado en el repositorio** (se restauraron las de la corrida que produjo el CSV)
> precisamente para no mover la tabla de §3.4 del paper por la puerta de atrás: antes de
> dar por buenos unos números nuevos, decidir si se re-acoplan esos dos positivos o se
> recuperan sus ficheros. Segunda: los tres fragmentos de Nshogoza siguen **sin SMILES** en
> `verdad_de_referencia.csv` (dice "pendiente: hay que dibujarlo"), y eso **rompía la
> corrida entera** con un `KeyError` al calcular el AUC por átomo; ya no rompe (su número
> de átomos pesados se lee del fichero preparado y coincide con los publicados: 13, 11 y
> 15), pero el SMILES sigue pendiente.

**El orden importa, y por eso está numerado.** El validador es el que prepara y acopla el
positivo nuevo, y lo deja en `validar_tdp43_limpia/ligands`; el barrido de exhaustividad y
la corrida de GPU **copian de ahí** y solo acoplan lo que no tenga pose válida. Si se
corrieran al revés, la corrida de GPU se dejaría fuera el positivo nuevo y el veredicto se
leería sobre conjuntos distintos.

**Lo que se espera.** Contra los señuelos emparejados de TDP-43 hacen falta 3–6 positivos
más; contra el fondo duro de SOD1, del orden de 11. Con uno o dos no se cierra el
intervalo, y eso está bien: la regla es barata de repetir, así que se corre con cada
positivo que llegue y se lee el intervalo, no el titular. Un solo positivo no cierra nada;
lo que sí hace es empezar a estrecharlo.

**Y lo que también puede pasar, y hay que reportarlo igual:** que con los nuevos el
veredicto **cambie de signo** —que el bloque duro de TDP-43 deje de pasar, por ejemplo, si
el positivo nuevo es un compuesto que el embudo no recupera—. Eso no es un fallo del
experimento: es la regla funcionando, y se escribe tal cual.

### B. Llega un negativo (se midió y no une)

No entra en la verdad de referencia: la verdad es de unidos medidos. Pero **no se tira**,
porque es material que hoy no existe en el proyecto: el primer **señuelo con inactividad
medida**. Hoy los fondos son parecidos en propiedades, no inactivos conocidos, y eso está
escrito como limitación.

**Ya está montado** (25 sep 2026), y a propósito antes de tener el dato: así el día que
llegue no hay que decidir nada con el número delante. Los compuestos que la propia
petición propone medir ya están preparados y acoplados con el mismo receptor, caja y
exhaustividad que la validación de su diana (`analysis/_medidos_no_unen/`, con
`fondo_inactivos.py` y su `manifiesto.csv`), la forma de leerlos está declarada de
antemano en `_medidos_no_unen/DECLARACION.md` —bloque aparte, que no se mezcla con los
otros dos; el compuesto entra con su ensayo, su constructo y su rango; el sitio manda;
piso de 30 compuestos; decide el AUC por átomo pesado, nunca el crudo— y la regla ya
reconoce ese bloque (papel `inactivo`).

Hoy el bloque **no sale** en `regla_decision.csv` y los números publicados no se han
movido ni un decimal, porque no hay ningún negativo: el manifiesto está entero en
`pendiente_medicion` y se comprobó que la salida de la regla es idéntica byte a byte a la
congelada. El día que llegue el primero se cambia su `estado` a `no_une` (con su ensayo y
su cita) y el bloque aparece solo. Lo que no se puede hacer es meterlo en `fondo` sin
decirlo: eso cambiaría los bloques que ya están publicados.

### C. Llega el dato del CR (XL23)

El CR es **otra caja**, así que este dato no entra en la regla del RRM ni la mueve. Lo que
hace es abrir la línea que hoy está cerrada, y conviene tener el criterio escrito antes de
tener el número:

- **Criterio del CR, fijado aquí:** si XL23 une, el control pasa a ser *¿reproduce el
  acoplamiento la dependencia de sitio medida?* — el acoplamiento tiene que ordenar XL20 y
  XL23 igual que los ordena el ensayo, y la lectura que se usa es la **por átomo pesado**,
  que es la que pasó la regla del 24 y la que no premia el tamaño (en crudo XL23 ganaba
  18 de 18 sin tener unión medida: esa lectura está descartada para decidir).
- **Dos puntos es poco, y se dice:** con dos compuestos, coincidir puede ser casualidad.
  El control se puede montar, la afirmación que sale de él es débil, y así se escribe.
- Si **no** une, la pareja XL20–XL23 queda como una cabeza común con dos colas, sin más, y
  el CR sigue cerrado. También es un resultado.

## 5. Qué NO se hace, pase lo que pase

No se mueve el umbral, ni la semilla, ni los remuestreos. No se cambia el fondo ni se parte
en trozos nuevos para encontrar un bloque que pase: eso es exactamente lo que hacía el
corte por quimiotipos al que esta regla sustituyó, y está documentado por qué se retiró. No
se retira un positivo porque puntúe mal. No se cuenta un análogo como química nueva. Y no
se reporta como PASA un bloque cuyo límite inferior no llega: `SIN EVIDENCIA` es un
resultado, no un fracaso.

## 6. Qué queda escrito al terminar

Cada tanda de ensayos que cambie la verdad de referencia deja: su carpeta con fecha en
`analysis/regla_decision/`, la fila nueva con su cita en `verdad_de_referencia.csv`, y el
párrafo en el informe del día con los intervalos **antes y después**. La tabla de la
sección 3.4 del paper se actualiza **solo** con una corrida cuyo veredicto sea el mismo en
todas las corridas; si no lo es, se escribe que no lo es.

## 7. Cómo se reproduce

```
python analysis/incorporar_positivo.py --help        # el portero y sus comprobaciones
python analysis/regla_decision_bootstrap.py          # la regla, con lo que haya hoy
```

- El portero está en `analysis/incorporar_positivo.py`; el veredicto congelado, en
  `analysis/regla_decision/_antes_2026-09-25/` (el log, como `.txt`, porque `.log` está en
  el `.gitignore`).
- La regla importa `auc` y `residual` de `validar_sod1_limpia.py`, que sigue siendo la
  única copia de las métricas.
