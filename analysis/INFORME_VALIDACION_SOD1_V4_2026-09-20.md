# Validación de SOD1 v4: por fin hay quimias independientes en Trp32

> **AVISO DEL 21 DE SEPTIEMBRE DE 2026 — ESTE INFORME QUEDÓ SIN EFECTO EN SU
> CONCLUSIÓN.** El control de redocking que él mismo dejaba pendiente (sección 5,
> punto 1) se hizo al día siguiente y **falla en las siete quimias nuevas**: la
> pose mejor puntuada queda a 4,3–10,5 Å de la del cristal en todas ellas, la pose
> cristalina es un mínimo local de la función y aun así se puntúa 1,0–2,8 kcal/mol
> peor que otra colocación, y con el ligando rígido tampoco se recupera. Además,
> la anomalía de las quinazolinas (sección 4) ya está explicada: los −10,0 y −9,8
> no son una pose de unión real, sino una colocación en una grieta entre dos
> copias de la proteína del modelo de 1HL5 (cadenas A y J). En el receptor de sus
> propias estructuras 4MQ da −5,6 y 12I −5,9. **El ranking de esta validación no
> se puede reportar tal cual.** Ver `redocking_trp32/INFORME_REDOCKING_QUIMIAS_NUEVAS_2026-09-21.md`.
> Lo que sigue siendo válido de este informe: el conjunto de controles limpio de
> 10 quimias independientes, el sesgo de tamaño medido, y el hecho de que las
> catecolaminas puntúan mal.

**20 de septiembre de 2026 (noche).** Scripts: `ligandos_cristal_trp32.py` y
`validar_sod1_v4.py`. Datos: `analysis/trp32_cristal/ligandos_trp32.csv`,
`analysis/validacion_SOD1_v4/analisis_sod1_v4.csv`,
`detalle_controles_sod1_v4.csv`. Controles: `controles_sod1_v4.csv`.
Antecedentes: `INFORME_VALIDACION_SOD1_2026-09-20.md` (v3) y
`redocking_trp32/INFORME_CONTROLES_CORREGIDOS_2026-09-20.md`.

---

## 0. Lo que se buscaba y por qué

El criterio que el proyecto se fijó el 20 de septiembre exige **positivos de al
menos dos quimias independientes** para poder reportar un ranking. Hasta hoy
solo había dos compuestos independientes (LCS-1 y PRG-A01) y **ninguno con
estructura**. Con dos no se sostiene nada, dé lo que dé la métrica.

En vez de buscar más artículos al azar, se hizo lo sistemático: **recorrer las
156 estructuras de SOD1 humana depositadas** en el PDB y medir qué moléculas
pequeñas están a menos de 6 Å del anillo de Trp32. Cada fichero trae su propia
Trp32 y su propio ligando en el mismo sistema de coordenadas, así que no hay que
superponer nada.

## 1. El barrido del PDB: siete quimias nuevas en Trp32

`ligandos_cristal_trp32.py` recorre las 156 entradas (146 sirven en formato PDB),
descarta disolvente, iones y aditivos de cristalización, y mide distancias. En
`trp32_cristal/ligandos_trp32.csv` queda cada contacto. Resultado:

| Código | PDB | Química | Dist. min. a Trp32 | Evidencia |
|---|---|---|---|---|
| 5FW | 4A7T | isoprenalina (catecolamina) | 3,13 Å | co-cristalización |
| ALE | 4A7U | adrenalina (catecolamina) | 3,69 Å | co-cristalización |
| LDP | 4A7V | dopamina (catecolamina) | 5,92 Å | co-cristalización |
| 5UD | 4A7S | 5-fluorouridina (nucleósido) | 2,58 Å | co-cristalización |
| 12I | 4A7G | 4-(4-metilpiperazin-1-il)quinazolina | 3,29 Å | co-cristalización |
| 4MQ | 4A7Q | 4-(4-metil-1,4-diazepan-1-il)quinazolina | 3,41 Å | co-cristalización |
| ZO0 | 2WZ6 | 4-(4-metil-1,4-diazepan-1-il)-2-(CF₃)quinazolina | 3,18 Å | co-cristalización |
| ZZT | 2WZ0 | 2-metoxi-5-metilanilina | 3,33 Å | co-cristalización |
| **K4I** | **8GSQ** | **paliperidona** (fármaco aprobado) | 3,33 Å | co-cristalización + MST |
| **6B3** | **6A9O** | **fenantridinona (Lig9)** | 3,39 Å | co-cristalización + inhibición de oxidación de Trp32 |
| 946 | 5YTO | aminoalcohol naftalénico | 3,45 Å | co-cristalización |

Dos de ellos son especialmente valiosos: **paliperidona** (antipsicótico
aprobado; Aouti et al. 2023, *Acta Cryst* D79:531, con afinidad por MST) y
**Lig9** (Manjula et al. 2019, fenantridinona que inhibe la oxidación de Trp32).
Los demás vienen de la misma serie de estructuras de Wright et al. (2013) que dio
las catecolaminas, pero son quimias completamente distintas.

Con esto el conjunto de controles pasa de **2 quimias independientes** a **10
quimias** (Murcko) entre 14 compuestos: catecolaminas, nucleósido, tres
quinazolinas, anilina, fenantridinona, benzisoxazol-piperidina, aminoalcohol
naftalénico, piridazinona, cumarina y pirazolona.

## 2. El resultado del acoplamiento

Todos los controles se acoplaron en la **misma caja Trp32** del cribado (46,5
80,0 73,3; 22 Å), con el **mismo receptor** (`gpu_dock/SOD1.pdbqt`), la misma
exhaustividad (8) y el mismo pipeline de preparación que los señuelos. Los
fondos son los ya acoplados en la v3: 140 señuelos emparejados nuevos, 143 del
fondo duro (quelantes/redox) y 199 emparejados antiguos.

### Puesto de cada control (frente a los 482 señuelos juntos)

| Control | Química | Afinidad | Átomos pesados | Puesto | Mejor que |
|---|---|---|---|---|---|
| diazepanoquinazolina (4MQ) | quinazolina | **−10,02** | 20 | **8** | 98,3 % |
| cf3quinazolina (ZO0) | quinazolina | **−9,77** | 24 | **10** | 97,9 % |
| Lig9 (6B3) | fenantridinona | −6,45 | 33 | **26** | 94,6 % |
| paliperidona (K4I) | benzisoxazol | −6,14 | 31 | 51 | 89,4 % |
| PRG-A01 | cumarina | −6,12 | 31 | 51 | 89,4 % |
| CHEMBL2165613 | pirazolona | −5,98 | 20 | 62 | 87,2 % |
| 5-fluorouridina | nucleósido | −5,35 | 18 | 200 | 58,6 % |
| naftalenoaminoalcohol (946) | aminoalcohol | −5,34 | 23 | 202 | 58,2 % |
| quinazolina (12I) | quinazolina | −4,97 | 17 | 311 | 35,6 % |
| LCS-1 | piridazinona | −4,86 | 16 | 348 | 28,0 % |
| adrenalina | catecolamina | −4,77 | 13 | 363 | 24,8 % |
| isoprenalina | catecolamina | −4,72 | 15 | 377 | 21,9 % |
| dopamina | catecolamina | −4,70 | 11 | 380 | 21,3 % |
| anilina (ZZT) | anilina | −4,58 | 10 | 407 | 15,7 % |

### AUC por familia y fondo

| Familia | Fondo | AUC | AUC sin tamaño | EF1% | EF5% |
|---|---|---|---|---|---|
| catecolaminas (1 quimia, 3 comp.) | los tres | 0,228 | 0,669 | 0 | 0 |
| **cristalizadas, quimias nuevas (8)** | emparejado nuevo | **0,724** | 0,617 | 0 | 5,0 |
| **cristalizadas, quimias nuevas (8)** | emparejado antiguo | **0,688** | 0,678 | **25,0** | 7,5 |
| **cristalizadas, quimias nuevas (8)** | duro | **0,650** | 0,600 | 0 | 5,0 |
| **cristalizadas, quimias nuevas (8)** | los tres juntos | **0,687** | 0,637 | 0 | 5,0 |
| Lig9 sola | los tres | 0,948 | **0,390** | 0 | 0 |
| paliperidona sola | los tres | 0,896 | **0,328** | 0 | 0 |
| sin estructura (LCS-1, PRG-A01) | los tres | 0,588 | 0,425 | 0 | 0 |
| serie pirazolona (1) | los tres | 0,873 | 0,902 | 0 | 0 |
| cristalizadas + sin estructura (10) | los tres | 0,667 | 0,595 | 0 | 4,0 |
| todos los controles (14) | los tres | 0,588 | 0,633 | 0 | 2,86 |

Sesgo de tamaño medido esta vez **con los 482 señuelos** (contando los átomos
pesados del propio PDBQT acoplado, no de un SMILES emparejado por nombre):

```
afinidad = -3,24 - 0,105 * (átomos pesados)
```

## 3. Lo que cambia

1. **El criterio se cumple por primera vez.** En la lista ordenada por afinidad
   entran en el **1 % superior** dos compuestos (4MQ, ZO0) y en el **5 %
   superior** tres, y hay **dos quimias independientes** por delante del fondo:
   la **quinazolina** (4MQ, ZO0) y la **fenantridinona** (Lig9). Hasta hoy el
   mejor resultado del proyecto tenía **cero** quimias independientes en la
   cabeza.
2. **El conjunto de positivos ya no es el problema.** Es la primera vez que se
   puede decir: los positivos son 10 quimias distintas, **todas con estructura
   resuelta o ensayo directo**, y aun así el número sale.
3. **El AUC de las cristalizadas es 0,69 (0,64 sin tamaño): por encima del azar,
   pero «débil», no «validado».** La diferencia con las catecolaminas solas
   (0,23 bruto) es enorme y muestra que la familia importa: las catecolaminas,
   que son las que definen el bolsillo, puntúan mal; las quimias nuevas, mejor.
4. **Quitando el tamaño, los dos mejores individuales se caen**: Lig9 pasa de
   0,948 a 0,390 y la paliperidona de 0,896 a 0,328. Su buena posición bruta es
   **tamaño** (33 y 31 átomos pesados) — exactamente el sesgo que ya estaba
   medido en la v3. Las quinazolinas **sí** sobreviven al descuento (residuo
   −4,68 y −4,01), que es lo único que se puede considerar señal real hoy.
5. **La serie pirazolona sigue ganando a todo** (0,90 sin tamaño) siendo un solo
   compuesto: es la prueba de que la métrica premia el tamaño y la lipofilia, no
   el reconocimiento.

## 4. Dos cosas que hay que verificar antes de dar esto por bueno

1. **Las dos quinazolinas extremas.** 4MQ y ZO0 puntúan −10,0 y −9,8 kcal/mol,
   mientras 12I —misma química, un anillo de piperazina en vez de diazepano—
   puntúa −4,97. **Cinco kcal/mol por un CH₂ no es físico.** Las dos diazepano
   concuerdan entre sí, así que no es ruido: es geometría.

   **Primera hipótesis, medida y descartada.** El enterramiento no lo explica
   (`enterramiento_quinazolinas.py`, sobre las poses ya calculadas):

   | Pose | Afinidad | Pesados | Receptor ≤4,5 Å | ≤5,5 Å | Residuos |
   |---|---|---|---|---|---|
   | 4MQ | −10,02 | 20 | 25 | 44 | 12 |
   | ZO0 | −9,77 | 24 | 22 | 43 | 14 |
   | 12I | −4,96 | 17 | 22 | 42 | 12 |

   Las tres están **rodeadas por prácticamente el mismo número de átomos de
   receptor** (22–25 a 4,5 Å, 42–44 a 5,5 Å), y sin embargo difieren en 5
   kcal/mol. No es cuánto las rodea, es **cómo**: la diferencia tiene que estar
   en el término de interacción concreto (forma del anillo de diazepano frente a
   la piperazina), y eso solo se ve mirando la pose. Siguiente paso obligado:
   **re-acoplar 4MQ y ZO0 en su propia estructura cristalográfica** (4A7Q y
   2WZ6) y medir el RMSD frente a la pose depositada, que es la única referencia
   independiente que existe.
2. **La pose real.** Ninguno de estos ligandos se ha comparado aún con su pose
   cristalina. El control de redocking está hecho para las catecolaminas y no
   para estas quimias: es el siguiente paso obligado, y ahora es posible porque
   cada una tiene su estructura.

## 5. Qué queda

1. Redocking de las 7 quimias nuevas en su propia estructura (RMSD ≤ 2 Å) antes
   de usar ninguna como patrón.
2. Explicar la anomalía 4MQ/12I con medidas. El enterramiento ya está medido y
   descartado (`enterramiento_quinazolinas.py`); queda comparar las poses y los
   términos de interacción frente a la referencia cristalina.
3. Si la señal de las quinazolinas se sostiene, es la primera familia con la que
   el proyecto puede reportar algo: son dos compuestos cristalizados y una
   familia química concreta con potencial sintético.
