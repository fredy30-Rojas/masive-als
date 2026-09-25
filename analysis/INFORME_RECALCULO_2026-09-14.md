# Recálculo del rescoring MM-GBSA con el protocolo validado

**Fecha:** 14 de septiembre de 2026
**Entrada:** `rescoring_local/rescoring_recalculo_validado.csv` (242 filas)
**Comparado con:** `VALIDACION_RESCORING_2026-09-11.md`
**Script del análisis:** `rescoring_local/analizar_recalculo.py`
**Salida cruda:** `rescoring_local/informe_recalculo_2026-09-14.txt`

---

## 1. Qué se hizo

Se repasaron 240 moléculas distintas con el protocolo que se validó el 11 de
septiembre, para comprobar si el diagnóstico de aquel día se sostenía o había
sido un accidente de aquella muestra concreta.

El resultado cabe en una frase: **se sostiene, y con una fidelidad casi
idéntica.** El dG bruto vuelve a ordenar al revés, y en TDP43_v2 no ordena ni
quitándole el tamaño.

---

## 2. ¿Es comparable con la validación del 11 de septiembre?

Sí, y era la comprobación previa obligatoria. Sin esto no se puede comparar
nada.

Los tres receptores de este recálculo tienen el **mismo md5** que los del
snapshot congelado `receptores_fijos/snapshot_20260911/manifest.json`:

| Diana | md5 | Átomos |
|---|---|---|
| TDP43_v2 | `1a17b49efaf7677d2909f5b45d0b3b4f` | 2067 |
| SOD1 | `f4331421f45c096d929c979dd1a06f43` | 1459 |
| FUS | `6c90f013f59aae2f5a0bf8706af198c8` | 601 |

Las 238 filas con resultado son `plataforma=CUDA`, igual que la validación.

**Caveat metodológico honesto:** el recuento de átomos pesados de este informe
sale del SMILES (RDKit, fragmento mayor), no del fichero de pose como en el
script del 11 de septiembre, porque las poses de este recálculo no están en
`rescoring_local`. Los valores salen sistemáticamente unos 2-3 átomos por
debajo (controles SOD1: 17,5 aquí frente a 19,5 entonces; fondo: 28,3 frente a
31,1), lo cual es coherente con que el fichero de pose cuente algo de más. El
desplazamiento es el mismo para controles y fondo, así que la comparación
interna se mantiene; y las diferencias que importan (AUC bruto 0,14 frente a
0,25; residual 0,70 frente a 0,75) están demasiado lejos de 0,5 como para que
dos átomos las den la vuelta.

---

## 3. La reproducción, número a número

| Métrica | Validación 11 sep | Recálculo 14 sep |
|---|---|---|
| SOD1 — controles / fondo | 18 / 60 | 18 / 68 |
| SOD1 — dG controles (mediana) | −15,8 | −17,7 |
| SOD1 — dG fondo (mediana) | −21,0 | −25,5 |
| SOD1 — **AUC bruto** | **0,253** | **0,143** |
| SOD1 — AUC tamaño-comparable | 0,430 | 0,230 |
| SOD1 — **AUC residual** | **0,754** | **0,703** |
| SOD1 — AUC eficiencia | 0,656 | 0,595 |
| SOD1 — Spearman dG vs tamaño | −0,508 | −0,585 |
| TDP43_v2 — controles / fondo | 8 / 60 | 8 / 71 |
| TDP43_v2 — dG controles (mediana) | −12,6 | −16,9 |
| TDP43_v2 — dG fondo (mediana) | −31,1 | −34,5 |
| TDP43_v2 — **AUC bruto** | **0,140** | **0,130** |
| TDP43_v2 — AUC tamaño-comparable | 0,139 | 0,108 |
| TDP43_v2 — **AUC residual** | **0,073** | **0,032** |
| TDP43_v2 — AUC eficiencia | — | 0,030 |
| FUS — controles | no evaluable | no evaluable |

Todas las métricas se mueven en el mismo sentido y con la misma magnitud. Las
dos que deciden (AUC residual) quedan a 0,05 y a 0,04 de las de septiembre.
Esto ya no es un resultado de una muestra: es el comportamiento del protocolo.

---

## 4. El sesgo de tamaño, medido

La pendiente del ajuste dG ~ átomos pesados es:

- **SOD1: −1,26 kcal/mol por átomo pesado.** Un ligando de 20 átomos y otro de
  35 (ambos corrientes) se llevan casi **19 kcal/mol solo por tamaño**.
- **TDP43_v2: −0,93 kcal/mol por átomo pesado.**

El sesgo es tan grande que se come la señal entera. La prueba más clara: **los
controles de SOD1 son MÁS PEQUEÑOS que el fondo** (17,5 frente a 28,3 átomos),
así que el tamaño juega contra ellos, y el AUC bruto sale 0,143 — muy por
debajo de 0,5, es decir, ordenados al revés de como deberían.

Y esto no es ruido. La dispersión interna entre réplicas (`dG_sd`) tiene
mediana **0,18 kcal/mol** y máximo 4,19. El cálculo es internamente estable:
el problema no es que MM-GBSA sea ruidoso, es que mide tamaño y lo llama
afinidad.

**Valores por debajo de −60 kcal/mol no son creíbles** como energías de
unión de moléculas de este tamaño. En este recálculo hay varios: −87,35
(TDP43_v2), −72,01 (SOD1), −54,65 y −54,51. Son artefactos de enterramiento,
no ligandos extraordinarios.

---

## 5. Dónde caen de verdad los controles

El número más fácil de leer de todo el informe. Si el orden sirviera, los
fármacos conocidos deberían quedar arriba.

**SOD1** (18 controles entre 86 moléculas):
- por dG bruto → percentil **18** (puesto 70 de 86) — abajo del todo
- por residual → percentil **72** (puesto 24 de 86) — arriba, pero no arriba del todo
- por eficiencia → percentil **69** (puesto 26 de 86)

**TDP43_v2** (8 controles entre 79 moléculas):
- por dG bruto → percentil **6**
- por residual → percentil **5**
- por eficiencia → percentil **5**

En TDP43_v2 los ligandos conocidos caen **al fondo de la lista por los tres
caminos**. Eso ya no se arregla con una métrica mejor: apunta a que el
problema está antes, en el receptor o en el bolsillo, no en cómo se puntúa.

---

## 6. La prueba de control que faltaba: ¿cuánto lo explica solo el tamaño?

Antes de concluir que el método está roto hay que descartar lo obvio. Si los
controles fuesen simplemente **más pequeños**, ordenar por tamaño daría el
mismo resultado que ordenar por dG, y no habríamos aprendido nada del método.

Se calculó el AUC usando como puntuación **solo el número de átomos pesados**
(más grande = peor, igual que en dG):

| Diana | AUC solo-tamaño | AUC dG bruto | Controles | Fondo |
|---|---|---|---|---|
| SOD1 | **0,019** | 0,143 | 17,5 átomos | 28,3 átomos |
| TDP43_v2 | **0,518** | 0,130 | 27,4 átomos | 26,7 átomos |

Esto parte el diagnóstico en dos, y es la conclusión más importante de todo el
informe:

**En SOD1 el tamaño lo explica casi todo.** Ordenar por número de átomos ya da
0,019, prácticamente la separación perfecta en el sentido equivocado. El dG
bruto (0,143) apenas mejora eso. Es decir: **el AUC bruto de 0,253 que se
reportó el 11 de septiembre no era evidencia de que el método falle, sino de
que los controles de SOD1 son más pequeños que el fondo.** Sólo el residual
(0,703) mide lo que de verdad aporta el método, y ahí sí hay señal.

**En TDP43_v2 el tamaño no explica nada** (0,518 = azar puro; 27,4 frente a
26,7 átomos). Los controles están bien igualados en tamaño con el fondo y
**aun así caen los últimos por los tres caminos**. Aquí no hay artefacto de
tamaño que valga: **el fallo es real y específico de esta diana.**

### Lo mismo pasa ya en el acoplamiento

El AUC del propio Vina, en las mismas 86 moléculas de SOD1, sale **0,000**:
los 18 controles puntúan peor que los 68 fondos, sin una sola excepción. Como
los controles tienen 17,5 átomos y el fondo 28,3, el sesgo de tamaño ya está
instalado antes de que empiece el MM-GBSA. **El recálculo no puede arreglar lo
que el acoplamiento ya entrega sesgado.**

### Las poses de TDP43_v2 están dispersas

Medido sobre las poses reales: en TDP43_v2 los 8 controles tienen un centro
disperso de **8,1 Å**, mientras los 71 fondos se agrupan en **2,2 Å**. Ninguna
pose se va del receptor (todas tocan a menos de 3 Å). En SOD1 los dos grupos se
solapan (2,5 Å entre centros, 4,9 frente a 3,1 Å de dispersión).

Que los ligandos conocidos de TDP43_v2 acaben cada uno por su lado, mientras el
fondo converge a un sitio, apunta a que **el sitio de acoplamiento que se está
usando no es el bueno para ellos**. Es la pista más concreta que hay para esa
diana, y explica por qué ninguna métrica la salva: no es un problema de
puntuación, es un problema de colocación.

---

## 7. Cómo ordenar, y cómo no

**No ordenar nunca por dG bruto.** Confirmado por segunda vez. En SOD1 pone a
los buenos en el percentil 18; en TDP43_v2 en el 6. Y ahora sabemos por qué en
SOD1: el dG bruto es, en la práctica, un contador de átomos disfrazado.

**SOD1:** el residual es lo único defendible (AUC 0,703). Es el único número
que mide lo que el método aporta de verdad, porque el bruto está dominado por
el tamaño. La eficiencia de ligando (AUC 0,595) es casi azar. El residual sirve
para *descartar*, no para elegir. Y cualquier comparación futura en SOD1 debe
hacerse **dentro de una ventana de tamaño**: de los 68 fondos, sólo 21 caen en
el rango de 15-26 átomos de los controles.

**TDP43_v2:** nada sirve, y aquí no vale echarle la culpa al tamaño (sección
6). El residual (AUC 0,032) ordena **al revés** que el dG bruto. Ordenar los
candidatos por residual aquí sería seleccionar sistemáticamente los peores.
**Hasta que no se arregle el sitio de acoplamiento, ningún número de TDP43_v2
debe usarse para decidir nada.**

**FUS:** sin controles no se puede decir nada. Su ranking existe pero **no está
validado** y no se debe usar para decidir.

### Candidatos que salen arriba (con la cautela puesta)

Por residual, SOD1: CHEMBL520254 (−72,01; 34 átomos; residual −37,6),
CHEMBL4570385 (residual −23,8), CHEMBL521548 (−21,5), CHEMBL561610 (−18,2),
CHEMBL500921 (−12,1), GLASDEGIB (−11,1).

**CHEMBL520254 hay que mirarlo antes de creérselo.** Es el primero por las tres
métricas a la vez, y su residual (−37,6) está a más de 13 kcal/mol del
siguiente. Con un AUC de controles de solo 0,70, un outlier así es más
sospechoso de artefacto que de hallazgo. Toca inspeccionar su pose.

---

## 8. Errores y reintentos

240 moléculas distintas, 238 con resultado, **2 fallos reales (0,8 %)**:
`TDP43_v2` cepharanthine y cepharanthine_limpio, las dos por
`n atomos no coincide: SMILES 45 vs pose 47`.

**Corrección (14 sep, tarde): no es un problema de protonación.** El SMILES es
correcto — la cefarantina es C37H38N2O6, 45 átomos pesados, y el SMILES da
exactamente eso. La que está mal es **la pose**: tiene 39 carbonos donde la
molécula tiene 37, o sea **dos carbonos de más**. La pose que hay guardada como
cepharanthine no es cefarantina. No se arregla reprotonando: hay que reacoplar,
o dar por perdidas esas dos filas y revisar la entrada de la librería.

Además, CHEMBL1206581 y CHEMBL4434647 fallaron primero con `salida vacia` y
**el reintento los resolvió** (−38,08 y −24,36). La reanudación del runner
funciona.

Tiempo de GPU sumado de las 238 filas: **77,2 horas** (media 19,5 min por
molécula, mediana 16 min, de 40 s a 4389 s).

---

## 9. Qué haría falta ahora

> **Los puntos 1, 3 y 4 ya están hechos esa misma noche: ver la sección 11.**
> El 1 se resolvió midiendo el enterramiento (los controles entierran *menos*
> que el fondo, y son intercaladores de ácido nucleico), el 3 se resolvió
> consultando ChEMBL (FUS no tiene ligandos, punto), y el 4 se resolvió
> descartando la pose y encontrando de qué enzima es de verdad el compuesto.

1. **TDP43_v2: revisar el SITIO de acoplamiento, no la métrica.** Es lo
   prioritario, y por primera vez hay una pista concreta: los 8 controles
   acaban dispersos 8,1 Å mientras el fondo converge en 2,2 Å. Hay que mirar si
   el bolsillo que se usa es el sitio real de unión y si las poses de los
   controles son razonables. Mientras eso no se arregle, ningún número de esa
   diana vale para decidir.
2. **SOD1: dejar de comparar controles contra un fondo de otro tamaño.** De 68
   fondos, sólo 21 caen en el rango de tamaño de los controles. Cualquier
   validación futura en SOD1 tiene que ser **dentro de una ventana de tamaño**,
   y juzgarse por el residual (0,703), no por el bruto (0,143).
3. **FUS: buscar controles propios.** Es la única diana grande sin evaluar, y
   sin controles su ranking no vale para nada. **Sin controles no se puede
   decir absolutamente nada de FUS**, y conviene no olvidarlo al leer listas.
4. **CHEMBL520254: inspeccionar la pose** antes de meterlo en cualquier lista,
   que es el primero por las tres métricas pero con un residual a 13 kcal/mol
   del siguiente.
5. **Cefarantina: la pose es de otra molécula** (39 carbonos frente a 37). Hay
   que reacoplar esa entrada o dar por perdidas las dos filas, y **revisar la
   librería de ligandos por si hay más entradas mal etiquetadas**.
6. **Los 2.973 de la lista completa siguen sin poder ordenarse.** Este recálculo
   confirma que ordenar por dG bruto no sirve para decidir candidatos.

---

## 10. Auditoría de toda la librería (14 sep, tarde)

Después de ver que la pose de la cefarantina era de otra molécula, se auditó
**toda** la librería: para cada pose se comparó su composición de elementos
pesados contra la del SMILES de su ligando. Script: `analizar` →
`analysis/auditar_libreria.py`; detalle en `analysis/poses_sospechosas.csv`.

| | |
|---|---|
| Poses revisadas | 430.440 (4 dianas) |
| Comparadas contra su SMILES | 60.750 |
| **Coinciden** | **60.712 (99,94 %)** |
| No coinciden | 38 |
| Sin SMILES en ningún CSV | 369.652 |

**La librería está limpia.** El 99,94 % de las poses son de la molécula que
dicen ser. Los 38 desajustes se reparten así:

- **20 poses degeneradas** (casi vacías, con 1-3 átomos).
- **4 poses fragmentadas** (les falta medio ligando).
- **14 de molécula distinta**, que son en realidad **dos patrones**:
  - **12 = tres compuestos con boro** (CHEMBL15632, CHEMBL157117,
    CHEMBL4553125), cada uno en las cuatro dianas. El SMILES tiene boro y la
    pose tiene un carbono en su lugar: el fichero PDBQT trae el átomo
    **nombrado "B" pero tipado como "C"**. La preparación convirtió el boro en
    carbono, así que esos compuestos se acoplaron y puntuaron **como si fueran
    su análogo de carbono**. No son ellos.
  - **2 = cefarantina** (la de la sección 8).

**Cautela honesta:** sólo se pudieron comprobar 60.750 de las 430.440 poses,
porque las otras 369.652 no tienen SMILES en ningún CSV de la carpeta. En la
parte comprobada el error es de 0,06 %, pero **no se puede afirmar lo mismo de
la parte no comprobada**.

---

## 11. Los tres cabos sueltos, atados (14 sep, noche)

Scripts: `rescoring_local/diag_energia_ligando.py` y
`rescoring_local/diag_bolsillo_candidatos.py`.

### 11.1 CHEMBL520254 no es un candidato: es un falso positivo con nombre

Era el primero por las tres métricas. Se ha mirado por dentro y por fuera.
Primero, las hipótesis de artefacto — **todas descartadas, midiendo**:

| Hipótesis | Medición | Veredicto |
|---|---|---|
| La energía interna del ligando infla el dG | correlación e_ligand↔dG = +0,13 (SOD1), −0,01 (TDP43_v2), +0,03 (FUS). El término varía 434 kcal/mol pero **se cancela** | descartada |
| Está enterrado en una cavidad falsa | 74 contactos frente a mediana 66 en su rango de tamaño (puesto 3 de 19) | descartada |
| Está en otro bolsillo | su centro está a **1,4 Å** del centro del fondo | descartada |
| Choca con el receptor | distancia mínima 2,13 Å | descartada |
| La carga exagera el MM-GBSA | es neutra, y las 86 de SOD1 lo son | descartada |

O sea: **la pose es físicamente normal y está exactamente donde están los
controles.** Ningún artefacto detectable. Y aun así es un falso positivo, por
lo que se ve al buscarlo fuera:

**CHEMBL520254 es un inhibidor de la leucotrieno A4 hidrolasa** (IC50 44 nM por
ensayo enzimático y 240 nM en sangre completa, *Bioorg Med Chem Lett* 2008).
Sin nombre comercial, sin fase clínica, y **nunca probado contra SOD1** en
ChEMBL.

Es un compuesto real y potente — contra otra enzima. Su farmacóforo
(carboxilato + amida + arilo hidrofóbico) es un **quelante de metales**, y el
sitio activo de SOD1 es un canal cargado que lleva a cobre y zinc. Un quelante
puntúa bien en cualquier bolsillo con metal sin ser ligando de verdad: es un
falso positivo clásico del acoplamiento en metaloenzimas.

**Veredicto: no perseguirlo.** No porque su pose esté mal, sino porque no hay
ninguna prueba independiente de que se una a SOD1.

### 11.2 TDP43_v2: el bolsillo premia al que llena hueco, y los buenos no lo llenan

La medida más limpia de todo el informe:

| | contactos por átomo | dG mediana |
|---|---|---|
| **TDP43_v2** controles (8) | **3,41** | −16,9 |
| TDP43_v2 fondo (71) | 4,27 | −34,5 |
| **SOD1** controles (18) | **2,63** | −17,7 |
| SOD1 fondo (68) | 2,19 | −25,5 |

En SOD1 los fármacos conocidos entierran **más** que la librería (2,63 frente a
2,19). En TDP43_v2 entierran **menos** (3,41 frente a 4,27): el fondo llena
mejor la cavidad que los ligandos verdaderos. Ordenando las 79 moléculas de
TDP43_v2 por enterramiento por átomo, los 8 controles caen en los puestos
**41, 63, 67, 72, 75, 77, 78 y 79**. La correlación contactos↔dG es −0,45: en
esta diana el dG mide, sobre todo, cuánto hueco rellenas.

**Y ahora se sabe por qué, mirando quiénes son los controles:** PE859,
berberrubina, **berberina**, **sanguinarina** y **nitidina**. Berberina,
sanguinarina y nitidina son alcaloides planos, aromáticos y con carga
permanente, intercaladores de ácido nucleico. Y **TDP-43 es una proteína que se
une a ARN**: su sitio real es la superficie plana de sus dominios RRM1/RRM2,
que es **somera** y da pocos contactos.

Se está acoplando en una cavidad profunda mientras los ligandos verdaderos se
unen a una superficie plana. El acoplamiento selecciona llenadores de cavidad,
no ligandos de TDP-43. Eso explica las tres cosas a la vez: que los controles
entierren menos, que sus poses se dispersen 8,1 Å (una superficie somera no
ancla) y que ninguna métrica los salve.

**El arreglo no es la métrica: es mover la caja al sitio de unión de ARN**
(RRM1/RRM2), o elegir un bolsillo definido del dominio C-terminal si el
objetivo es la agregación.

### 11.3 FUS no se puede validar, y no es por dejadez: no hay ligandos

La diana «RNA-binding protein FUS» es **CHEMBL5724679**. Tiene **diez registros
de actividad en total**:

- CHEMBL5653589 — Kd 36,6 nM
- CHEMBL3752910 — Kd 17,5 nM
- Molibresib — IC50 ≥ 10.000 nM (es decir, nada)

Los dos primeros salen de un **ensayo de pull-down con perlas de quinasa**
(*Kinobead*), un solo artículo y sin confirmación por otra técnica: el perfil
típico de un falso positivo de perlas promiscuas.

**FUS no tiene ligandos pequeños creíbles.** No es que falten controles por
dejadez — es que no existen, y por tanto **ningún ranking de FUS es
interpretable ni hay trabajo que lo arregle**. Si algún día se quiere intentar,
la vía razonable sería usar como controles alcaloides de unión a ARN (los
mismos berberina o sanguinarina), ya que FUS también es una proteína de ARN.

---

## Conclusión

El recálculo hizo lo que tenía que hacer: **reprodujo la validación**. Mismo
receptor, misma plataforma, mismos números, misma conclusión. El hallazgo del
11 de septiembre era real y no un accidente de muestra.

Y las pruebas de control de hoy afinan bastante más el diagnóstico, partiéndolo
en dos:

**En SOD1 el problema es de diseño experimental, no del método.** Ordenar por
número de átomos ya reproduce casi toda la separación (AUC 0,019), porque los
controles son mucho más pequeños que el fondo. El AUC bruto de 0,253 nunca fue
la prueba de que el método fallara. Lo que el método aporta de verdad es el
residual, 0,703, y es una señal modesta pero real: **sirve para descartar, no
para elegir.**

**En TDP43_v2 el problema es real.** Ahí el tamaño no explica nada (0,518) y
los ligandos conocidos caen los últimos por los tres caminos. Con la pista de
las poses dispersas, esto ya no es "la métrica no sirve", es **"estamos
acoplando en el sitio equivocado"**. Es el bloqueo de verdad y es donde hay que
trabajar.

Para FUS no se sabe nada, y conviene decirlo así de claro.**En una frase: no hay que afinar el ranking, hay que arreglar el sitio de acoplamiento de TDP43_v2 y conseguir controles para FUS.**

---

## 12. Los cabos sueltos, cerrados (17 de septiembre de 2026)

### 12.1 Los dos de FUS: eran reintentos, y están resueltos

Se comprobó fila por fila en `rescoring_recalculo_validado.csv`. Las filas 49
y 51 (`CHEMBL1206581` y `CHEMBL4434647`, ambas con `salida vacia`) **tienen su
fila de reintento en el mismo CSV**, la 61 y la 62, con los valores reales
(−38,08 y −24,36). El runner escribe una fila nueva al reintentar, no reescribe
la fallida, y por eso aparecen duplicadas. El recuento del apartado 8 se sostiene
tal como está: **240 moléculas distintas, 238 con resultado, 2 perdidas**.

### 12.2 Las dos filas de cefarantina: cerradas como perdidas, sin reacoplar

Se decide **no reacoplarlas**. No se reacoplan porque no son candidatas del
cribado: `cepharanthine` no aparece en `resultados_libreria_total.csv`, o sea
que no salió del top-100 por afinidad de ninguna diana, sino que entró por otra
vía. Reacoplar una entrada de la que no consta su origen no arregla nada: deja
dos números sueltos que no pertenecen a ninguna lista.

Quedan cerradas así: **perdidas, con la causa escrita, y el compuesto apartado**.
La causa no es la protonación (el SMILES es correcto, C37H38N2O6, 45 átomos
pesados): **la pose guardada es de otra molécula**, con 39 carbonos donde la
cefarantina tiene 37.

### 12.3 Los doce del boro: no valen, y hay que excluirlos de cualquier lista

Tres compuestos, cada uno en las cuatro dianas, doce filas. El fichero PDBQT
trae el átomo **nombrado «B» pero tipado como «C»**, así que se acoplaron y se
puntuaron como su análogo de carbono. Es decir: **los números que hay en la
librería para esos tres no son de esos tres compuestos.**

| Compuesto | TDP43_v2 | TDP43 | SOD1 | FUS |
|---|---|---|---|---|
| CHEMBL15632 | −8,8 | −7,5 | −6,1 | −5,7 |
| CHEMBL157117 | −8,4 | −6,0 | −5,8 | −6,2 |
| CHEMBL4553125 | −7,8 | −7,5 | −7,1 | −6,7 |

Están en `analysis/compuestos_invalidos.csv`, junto con la cefarantina, para
que **cualquier ranking futuro los descarte sin tener que volver a descubrirlo**.
Si algún día interesan de verdad, son solo tres y la reparación es pequeña:
repreparar con el boro bien tipado y reacoplar. Mientras tanto, sus afinidades
no pueden usarse, y la de CHEMBL15632 contra TDP43_v2 (−8,8) está entre las
buenas de la diana: si alguien leyera la lista sin este aviso, lo elegiría.

### 12.4 Lo que queda abierto, y no por dejadez

- **TDP43_v2 es el único cabo que bloquea una decisión.** Mientras la caja no
  se mueva al sitio de unión de ARN, ningún número de esa diana vale (apartado 11.2).
- **FUS no tiene salida**: no hay ligandos pequeños creíbles con los que
  validar, así que ningún ranking de FUS es interpretable (apartado 11.3).
- **La auditoría de la librería está hecha a medias, y hay que decirlo así:**
  se pudieron comparar 60.750 poses de las 430.440 (las que tienen SMILES). El
  99,94 % de limpieza vale para esa parte, **no para el total**. Las 369.652
  restantes siguen sin comprobar.

### 12.5 La corrida de Oracle: retirada, no terminada

Quedaba una corrida de rescoring MM-GBSA en Oracle, parada desde el 12 de
septiembre en la fila 150 de 193. Lo natural habría sido relanzarla y acabarla.
**No se hace, y el motivo es de peso:** los receptores que usa ese runner no son
los congelados de la validación del 11.

| Diana | Átomos en la corrida de Oracle | Átomos del receptor congelado |
|---|---|---|
| SOD1 | 19.783 | 1.459 |
| TDP43_v2 | 1.401 | 2.067 |
| FUS | 391 | 601 |

El de SOD1 es el receptor viejo con las copias concatenadas (numeración 1..153
repetida). Y el script (`rescoring_mmgbsa_robusto.py`) **no trabaja en doble
precisión ni hace las tres repeticiones** del protocolo validado. Son, por
tanto, números del método que se descartó el 13 de septiembre: terminarlos
habría dejado ochenta y pico resultados falsos con apariencia de resultados.

Los ficheros se guardan en `~/mmgbsa/descartado_protocolo_viejo/` con su
`LEEME.md`; no se borran, para que quede la constancia de qué eran.

De paso apareció por qué el latido llevaba tres días dando esa corrida por
viva: **el `pgrep` se contaba a sí mismo** (su propia línea de comandos lleva el
patrón que busca), así que Oracle salía siempre «activo», y el progreso se leía
del último registro de una corrida muerta. Corregido en `latido_corazon.py`: el
patrón va con corchetes, que no coinciden consigo mismos, y **el progreso solo
se publica si hay un proceso vivo de verdad**.

---

## 13. El sitio de acoplamiento de TDP43_v2, probado — y descartado como causa (17 sep 2026)

Este era el único cabo que bloqueaba una decisión. La conclusión de las
secciones 5 y 11.2 era que se está acoplando en el sitio equivocado, y que por
eso ninguno de los caminos salva a los controles. **Se ha probado, y no es
así.** El experimento destapó además un sesgo en el propio diseño de la prueba.

### 13.1 La caja nueva, sacada de la estructura y no de la literatura

El receptor `TDP43_v2.pdbqt` sale de la estructura **4BS2**, que es el TDP-43
con **ARN UG-rico unido**. O sea que el sitio de unión de ARN no hay que
buscarlo en un artículo: está en el fichero del que salió el receptor.

Primero se comprobó que se puede comparar sin alinear nada: **mismo sistema de
coordenadas** (desviación media del CA de 0,00 Å sobre 174 residuos comunes).
Después se midió qué residuos tocan al ARN (a menos de 4 Å): **39 residuos**,
repartidos por **35,6 × 42,0 × 31,6 Å**, porque el ARN cruza los dos dominios.

Cuatro cajas posibles, **todas del mismo tamaño que la actual (26 Å)** para que
la comparación sea limpia (una caja más grande encierra más receptor por
construcción, y eso confundiría el enterramiento con el mérito):

| Candidata | Centro | Átomos del sitio de ARN dentro | Átomos de receptor dentro | Distancia a la caja actual |
|---|---|---|---|---|
| Interfaz de ARN completa | (22,7, 21,8, −10,2) | **452 de 657** | **511** | 7,6 Å |
| Plataformas aromáticas | (22,3, 25,6, −14,6) | 402 | 535 | 9,0 Å |
| RRM1 (RNP2-RNP1) | (30,6, 13,7, −3,9) | 360 | 476 | 14,0 Å |
| **La que usa el cribado** | (24,2, 16,9, −15,9) | 323 | 611 | — |

Elegida la primera: cubre el 69 % del sitio de unión de ARN **encerrando menos
receptor que la caja actual** (511 frente a 611), que es exactamente lo que se
quería (más sitio real, menos relleno de cavidad).

**Primer hallazgo, y ya corrige al informe:** la caja del cribado **no está en
otro bolsillo**. Está a 7,6 Å del centro de la interfaz de ARN y ya contiene el
49 % de sus átomos. El problema no era mirar a otro sitio, era mirar a medias.
Guardado en `rescoring_local/caja_tdp43v2_arn.json`.

### 13.2 El sesgo que apareció al hacer la prueba

Se acopló el mismo panel (los 8 controles + las 71 moléculas de fondo del
recálculo) en la caja nueva con Vina-GPU. Resultado: AUC **0,092** frente al
**0,000** de la caja actual. Los dos números son catastróficos… y los dos
engañan, por la misma razón:

> **Ese «fondo» son los 100 mejores por afinidad de Vina contra esta misma
diana.** Pedir a los ligandos conocidos que queden por delante de ellos no mide
si el método reconoce a los buenos: mide si los buenos puntúan mejor que los que
mejor puntúan. Es pedirles que ganen un concurso a los ganadores del concurso.

Con esa prueba no se puede concluir nada, ni a favor ni en contra. Así que se
repitió con un fondo honesto: **200 moléculas al azar de la librería entera
(113.304), sin haber pasado filtro ninguno**, con la semilla fijada para que el
sorteo se pueda repetir.

| | Caja actual | Caja de ARN |
|---|---|---|
| **AUC controles vs fondo aleatorio** | **0,657** | **0,655** |
| Percentil medio del control dentro del fondo aleatorio | 64 | 64 |

### 13.3 Qué significa esto, sin adornos

1. **El acoplamiento no está roto.** Contra fondo aleatorio, los ligandos
   conocidos de TDP-43 quedan por delante del 64 % de la librería (AUC 0,66).
   Es modesto, pero está por encima del azar y es lo normal en acoplamiento.
2. **Cambiar la caja al sitio de unión de ARN no cambia nada**: 0,657 → 0,655.
   La hipótesis «estamos en el sitio equivocado» **no se sostiene** como
   explicación de que los controles caigan al final.
3. **Lo que explica que caigan al final es el sesgo de la prueba**, no la
   geometría: caían al final *entre los mejores por Vina de su propia diana*,
   que es el grupo más difícil posible. La afirmación del apartado 5 («caen los
   últimos por los tres caminos») hay que leerla así, con el fondo entre
   comillas: caen los últimos entre los preseleccionados, no en la librería.
4. **Lo que este experimento NO toca:** el dominio del tamaño en MM-GBSA
   (apartado 4) sigue tan medido como estaba, porque aquello no era una
   comparación contra un fondo elegido, sino una pendiente de −1,26 kcal/mol por
   átomo pesado. Esa parte del informe se mantiene entera.

### 13.4 Lo que quedaría por hacer, y por qué NO se lanza solo

- **Rescoring MM-GBSA de las poses de la caja nueva:** 79 moléculas × ~20 min
  ≈ 26 h de GPU. No se lanza por su cuenta: la prueba de acoplamiento dice que
  el sitio no es el problema, así que lo razonable es decidir antes *para qué* se
  quiere el ranking, no gastar un día de GPU por inercia.
- **Y la decisión de fondo, que es suya:** el fondo del que salen los candidatos
  son los 100 mejores por Vina. Si esa selección no es neutral —y la sección
  13.2 demuestra que no lo es—, el sitio donde hay que buscar activos no es ahí,
  sino en un cribado contra fondo aleatorio o contra panel de diversidad.

Ficheros de este apartado: `analysis/tdp43v2_sitio_arn.py` (interfaz de ARN),
`analysis/tdp43v2_caja_arn.py` (elección de caja),
`gpu_dock/acoplar_panel_caja_arn.py` (panel de control),
`gpu_dock/acoplar_fondo_aleatorio_caja_arn.py` (fondo aleatorio),
`rescoring_local/panel_caja_arn.csv`, `rescoring_local/fondo_aleatorio_caja_arn.json`.

---

## 14. La solución: no es la energía, es la química (17 sep 2026, tarde)

Se buscó qué separa de verdad a los ligandos conocidos del fondo aleatorio, con
las poses que ya había calculadas (sin gastar GPU). Todo con el mismo panel: 8
controles contra 200 moléculas al azar.

### 14.1 Lo que se probó primero, y por qué no sirve

| Medida | AUC | Lectura |
|---|---|---|
| Afinidad de acoplamiento | 0,655 | lo mejor de lo clásico |
| Afinidad corregida por tamaño | 0,657 | corregir el tamaño no cambia nada |
| Apilamiento con aromáticas | 0,581 | informa algo, pero es una cuenta de átomos |
| Puentes de hidrógeno | 0,526 | azar |
| Contactos totales | 0,432 | **peor que el azar: es un contador de tamaño** |
| Contactos por átomo | 0,379 | peor aún |
| Afinidad + apilamiento (pesos iguales) | 0,629 | **combinar empeora**: mete ruido de tamaño |

Nada de esto pasa del 0,66. Y el detalle que da la pista está en una tabla de
percentiles: ordenando por afinidad, arriba quedan **sanguinarina (88),
coptisina (83), berberrubina (79), berberina (69), nitidina (68)** —todos
bencilisoquinolinas planas y catiónicas— y abajo el único que no es de esa
familia, el **ketoconazol (20)**. No es que el método no sepa puntuar: es que el
sitio de unión de ARN del TDP-43 reconoce una química concreta, y a eso la
energía sólo llega de refilón.

### 14.2 La vía que no se había probado: parecerse a los activos

Con pocos activos conocidos y una química clara, la búsqueda por similitud es lo
primero que se hace en un proyecto de verdad. Aquí no se había probado. Medida
igual, contra el mismo fondo aleatorio, con huellas ECFP4 (Morgan, radio 2):

| Medida | AUC |
|---|---|
| Similitud máxima con cualquier activo conocido | **0,885** |
| Media de las 3 similitudes más altas | **0,902** |
| Similitud máxima con la familia plana | **0,923** |
| Afinidad de acoplamiento (referencia) | 0,656 |

Cada molécula se compara contra los activos **menos ella misma** (si no, un
control se parecería a sí mismo y el número se inflaría solo). Los SMILES de las
200 del fondo se piden a ChEMBL por su identificador y quedan en caché.

**Cautela honesta, y es grande:** seis de los ocho activos son la misma familia
química, así que esta medida responde sobre todo a «¿es una bencilisoquinolina?».
Sirve para **ampliar una serie conocida** y para rescatar fármacos de esa familia;
**no descubre quimiotipos nuevos**. Y el acoplamiento sigue siendo útil para
otra cosa: desempatar dentro de la familia, que es donde el tamaño no varía tanto.

### 14.3 El embudo, y la lista que sale de él

Orden de trabajo que se propone, cada paso con su porqué medido:

1. **Filtrar por química**, no por energía: similitud con los activos conocidos.
   De las 56.772 moléculas de la librería con SMILES, **213 pasan el umbral de
   0,25** (179 esqueletos distintos, una vez quitadas las sales del mismo
   compuesto).
2. **Acoplar sólo ese conjunto** en el sitio de unión de ARN. Son 179 moléculas:
   minutos de GPU, no días. Hecho.
3. **MM-GBSA sólo al final y dentro de una ventana de tamaño**, si se quiere
   apurar el orden dentro de la familia. Ahí su sesgo de tamaño ya no puede
   mandar, porque todos los candidatos son del mismo tamaño.

Resultado en `rescoring_local/lista_focalizada_tdp43v2.csv` (ordenada por
similitud, con la afinidad como desempate y con el nombre y la fase clínica de
cada una cuando existe), y la versión sin deduplicar en
`..._completa.csv`. De los 179 esqueletos, **7 ya son fármacos con fase
clínica**: son candidatos de reposicionamiento de una serie que ya está en
clínica, no invenciones.

Y hay un detalle que vale la pena: de los primeros de la lista, **muchos no
tenían afinidad nunca calculada** porque no salieron en el top-100 por Vina. Es
decir, la vía clásica los había dejado fuera por definición —el top-100 por
energía no los contiene— y son justamente los más parecidos a los activos.

### 14.4 Qué quedaría por comprobar, dicho claro

Todo esto está medido contra **un fondo aleatorio**, que es honesto pero no es la
prueba de fuego: la prueba de verdad es experimental. Lo que se puede afirmar es
que en esta diana, con estos activos conocidos, la química ordena mucho mejor que
la energía (0,92 frente a 0,66), y que el embudo no inventa nada: filtra, acopla
y ordena con criterios que se pueden explicar uno a uno.

Scripts de este apartado: `analysis/puntuacion_interaccion.py`,
`analysis/combinar_puntuaciones.py`, `analysis/similitud_activos.py`,
`analysis/lista_focalizada.py`, `gpu_dock/acoplar_lista_focalizada.py`.

