# Barrido de exhaustividad en CPU: el criterio tambien se mueve solo

**24 de septiembre de 2026.** El control con Vina-GPU
(`INFORME_CONTROL_VINAGPU_TDP43_2026-09-23.md`, §5.3) dejo escrita una duda que no
podia cerrar: el veredicto de quimiotipos cambio al cambiar de motor, pero el cambio
podia venir del motor **o** de la perilla, porque `exhaustiveness` de Vina no es
`search_depth` de Vina-GPU. Este barrido mueve **solo la perilla, dentro del mismo
motor**: la misma CPU, el mismo Vina, el mismo receptor, la misma caja y los mismos
ficheros de ligando. Si el veredicto se mueve aqui tambien, entonces lo que no aguanta
es el criterio.

**Se mueve.** Con el mismo motor y el mismo conjunto, exhaustividad 8 y 16 dan **2 de 5
quimiotipos**, y exhaustividad 32 da **3 de 5**. El quimiotipo que cruza el corte al
subir el esfuerzo es `PE859` (piridil-pirazol), y lo hace por **0,59 kcal/mol**: su
afinidad pasa de -7,479 a -8,068 con el mismo ligando, el mismo receptor y la misma
caja. El "PASA por emparejado de tamano" no es una propiedad de la quimia: es una
funcion del esfuerzo de busqueda.

---

## 1. Que se movio y que no

Mismos 7 positivos de union medida, 122 señuelos emparejados por propiedades y 152
unidores de ARN de R-BIND 2.0. El punto medio (exhaustividad 16) ya existia y no se
repitio: es la validacion limpia `_validacion_TDP43` + `_validacion_TDP43_rbind`.
Aqui solo se añadieron los extremos, 8 y 32.

| Bloque | Metrica | exh 8 | **exh 16** (limpieza) | exh 32 | GPU (sd 20) |
|---|---|---|---|---|---|
| Señuelos 122 | AUC crudo | 0,301 | 0,301 | 0,288 | 0,246 |
| | AUC/atomo | 0,683 | 0,686 | 0,666 | 0,667 |
| | residual | 0,307 | 0,306 | 0,288 | 0,264 |
| | quimias (C) | 2/5 | 2/5 | **3/5** | 1/5 |
| | veredicto | PASA | PASA | PASA | **NO PASA** |
| R-BIND 152 | AUC crudo | 0,371 | 0,357 | 0,390 | 0,401 |
| | residual | 0,407 | 0,391 | 0,423 | 0,414 |
| | quimias (C) | 2/5 | 2/5 | **3/5** | 3/5 |
| Los dos 274 | AUC crudo | 0,340 | 0,332 | 0,345 | 0,332 |
| | residual | 0,361 | 0,356 | 0,362 | 0,345 |
| | quimias (C) | 2/5 | 2/5 | **3/5** | 3/5 |
| | pendiente de tamaño | -0,045 | -0,039 | -0,039 | -0,022 |

`EF5 %` es **0,00 en los tres** bloques y en las tres exhaustividades.

**Quien gana a los señuelos de su tamaño:**

| Quimiotipo | Positivo | exh 8 | exh 16 | exh 32 |
|---|---|---|---|---|
| bencilisoquinolina | berberrubine | si | si | si |
| isoxazol-piperidina | nTRD22 | si | si | si |
| piridil-pirazol | PE859 | no | no | **si** |
| piperidinil-pirimidina | rTRD01 | no | no | no |
| fragmento | fragmento_1/2/3 | no | no | no |

Tres cosas que se leen de la tabla y conviene no mezclar:

1. **8 y 16 dan lo mismo, casi al decimal** (AUC 0,301 frente a 0,301; mismo veredicto).
   Para ordenar y para el cuadro grande, el ajuste barato basta: no hay que gastar
   cuatro veces mas tiempo para obtener el mismo numero.
2. **El que se mueve es el bloque de los señuelos emparejados**, que es justo el que
   decide el veredicto (el fondo duro no decide porque no separa, §3). Pasar de 16 a 32
   sube a PE859 lo justo para cruzar el corte de su tamaño, y eso cambia el veredicto
   de 2 a 3 quimias.
3. **La direccion del movimiento no es la que uno esperaria.** Mas esfuerzo no endurece
   el criterio: lo ablanda. Aqui subir de 16 a 32 **mete** una quimia mas; en el control
   de GPU, contra una perilla distinta, **salio** una. El corte no mide cuanta quimia
   hay: mide donde cae la cola de un fondo de 122 ligandos, y esa cola se mueve.

**Lo que no se movio, en las cuatro corridas (8, 16, 32 y GPU):** el AUC crudo esta por
debajo del azar en los nueve bloques (0,246 a 0,401), el enriquecimiento de cabeza es
**cero**, el AUC por atomo pesado es el **unico** numero por encima del azar (0,666 a
0,738), y el fondo de verdaderos unidores de ARN **no ordena mejor** que los señuelos
emparejados. Nada de esto depende del esfuerzo.

---

## 2. Como corrio, y un incidente que se declara

Mismo montaje para las tres: 289 ligandos `.pdbqt` ya preparados con la receta canonica
(no se re-preparo nada), receptor `_tdp43_bolsillo_v2/4BS2_ph74.pdbqt`, caja de 26 A,
`num_modes 3`, Vina 1.2.3 de CPU, 14 procesos de un hilo.

| Exhaustividad | Ligandos | Tiempo | Poses |
|---|---|---|---|
| 8 | 287 por acoplar | **79,3 min** (0,3 min/ligando) | 288 (1 sin pose) |
| 32 | 289 por acoplar | ~205 min hasta la ultima linea del diario | 288 (1 sin pose) |

**El incidente.** La corrida a exhaustividad 32 **no llego a su ultima linea: el trabajo
se cayo**. El ultimo ligando, `RB_SM_0145` (un macrociclo del fondo de R-BIND), tardo
mas de dos horas y chococon el tope por ligando del guion (`timeout=7200`), que al saltar tumba el conjunto de
procesos. No se perdio nada: su pose quedo escrita en ese mismo momento y se verifico
completa (3 modelos, `REMARK VINA RESULT` presente), los 288 ficheros tienen resultado,
y la puntuacion se corrio despues sobre los 288. Que un
solo ligando de 289 se coma dos horas mientras la mediana va a 0,7 min no es un fallo
del montaje: es el comportamiento de Vina con un macrociclo grande y flexible, y queda
apuntado porque explica por que la corrida no termino sola.

**Las dos poses que faltan, dichas enteras** (son las mismas del control de GPU, no
nuevas):

1. `RB_SM_0085` (DB1273), del fondo duro: **no se puede preparar**. Sus dos anillos son
   de selenio y Meeko no sabe tipar el Se. Entra 152 de 153, y el que falta se declara.
2. `DEC_CHEMBL4543460`, señuelo: tiene un **acido boronico** y Vina no acepta el boro
   (`Atom type B is not a valid AutoDock type`). Ninguno de los dos motores puede
   acoplarlo, y por eso el bloque de señuelos es de 122 y no de 123.

---

## 3. Lo que NO dice

1. **No valida TDP-43.** Sigue sin reportarse. El conjunto es el mismo de siempre: 7
   positivos de union medida y 122 señuelos, y un AUC crudo por debajo del azar no se
   puede presentar como validacion.
2. **No arregla el criterio; lo descalifica.** El barrido estaba pensado para separar la
   causa (motor o perilla). Separa una parte: la perilla, sola, mueve el veredicto. El
   otro control dice que el motor tambien lo mueve. Con las dos piezas encima de la
   mesa, el resultado es que **"al menos dos quimiotipos por delante de los de su
   tamaño" da 1, 2 o 3 segun con que se corra**, y ese numero es el que el proyecto usa
   para decidir si una idea sigue o se tira.
3. **No dice que la GPU sea mala ni que la CPU sea buena.** Dice que el numero que se
   usa como arbitro no aguanta el cambio de condiciones. El juicio sobre los dos motores
   sigue igual que ayer: el cuadro grande coincide, el veredicto no.
4. **No explica cual de las tres quimias es real.** Las dos que ganan en **las cuatro**
   corridas son, segun la literatura, las dos que caen **fuera** de ese bolsillo:
   `nTRD22` es modulador alosterico del dominio N-terminal y la berberrubina se describe
   en la interfaz RRM1-RRM2. La que entra a exigir mas es `PE859`, que en las corridas
   baratas no cruzaba. Con 5 quimiotipos y 7 positivos no se puede decidir quien tiene
   razon; lo que se puede decir es que **la que decide el veredicto es la mas fragil de
   las cinco**.
5. **No es una ley del barrido.** El rango explorado es 8-32 en un solo conjunto. Que
   16→32 sea mas permisivo que 8→16 es un dato, no una tendencia demostrada: haria falta
   otra diana y otro conjunto de señuelos para hablar de pendiente.

---

## 4. Que toca ahora

- **La lista de candidatos no se toca por esto**, y la validacion de TDP-43 sigue sin
  reportarse. El barrido no cambia ningun candidato; cambia lo que se puede afirmar.
- **El corte por quimiotipos no puede seguir siendo la regla de decision** mientras no
  se defina con un error. Si el proyecto quiere una regla util, las dos opciones que se
  sostienen solas en las cuatro corridas son el **AUC por atomo pesado** (unico por
  encima del azar, 0,666-0,738) y el **residual**. Cualquiera de las dos habria que
  acompañarla de un intervalo (bootstrap sobre los señuelos) antes de decidir nada con
  ella, y con 5 quimiotipos el corte por quimias no da ese intervalo.
- **Para seguir con TDP-43 hace falta un conjunto con mas positivos**, no mas
  exhaustividad: subir el esfuerzo mueve el criterio sin acercarlo a la verdad.

---

## 5. Anexo: como se reprodujo

```
python C:/Users/Fredy/masive-als/analysis/barrido_exhaustividad_tdp43.py --todo
python C:/Users/Fredy/masive-als/analysis/barrido_exhaustividad_tdp43.py --acoplar 32
python C:/Users/Fredy/masive-als/analysis/barrido_exhaustividad_tdp43.py --puntuar 32
```

- Montajes y poses: `_barrido_TDP43/exh8` y `_barrido_TDP43/exh32` (`ligands`, `out`,
  `out_duro`). Punto medio: `_validacion_TDP43` + `_validacion_TDP43_rbind`.
- Puntuacion: `_barrido_TDP43/validar_exh8.csv`, `validar_exh32.csv` y sus
  `_resumen.csv`; los dos extremos juntos, en `_barrido_TDP43/barrido_resumen.csv`.
- Registros: `_barrido_TDP43/barrido.log` y `_barrido_exhaustividad_out.log`.
- La comprobacion de que el ligando caido dejo su pose entera:

```
grep -a -c '^MODEL' _barrido_TDP43/exh32/out_duro/RB_SM_0145_out.pdbqt   # 3
```
