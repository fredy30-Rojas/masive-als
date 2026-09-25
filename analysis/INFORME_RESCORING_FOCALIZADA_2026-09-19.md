# Rescoring de la lista focalizada de TDP-43 — cierre y auditoría

**Fecha:** 20 de septiembre de 2026
**Cubre:** la corrida lanzada el 17 sep, terminada el 19 sep, y la auditoría de la
métrica de química que la seleccionó.

---

## 1. La corrida: terminada

| | |
|---|---|
| Lista | `rescoring_local/lista_focalizada_rescoring.csv` (179 esqueletos distintos) |
| Salida | `rescoring_local/rescoring_focalizada_arn.csv` |
| Terminada | **19 sep 2026, 22:07** |
| Resultado | **179 de 179 con dG, cero fallos** |
| Última tanda (63 restantes) | 711,2 min (11,85 h) con 4 procesos |
| Protocolo | receptor congelado `snapshot_20260911` (TDP43_v2 md5 `1a17b49efaf7`), doble precisión en CUDA, 3 repeticiones |
| Poses | acopladas en el sitio de unión de ARN (`gpu_dock/rescoring_caja_arn/results_TDP43_v2`) |
| Registro | `runner_focalizada_arn.log`, `rescoring_focalizada_stdout.log` |
| Aviso | `vigilante_rescoring.log` — 22:11, Telegram entregado; el vigilante terminó solo |

Se había parado dos veces: por suspensión del equipo la noche del 18 (22:52, en la
fila 116) y se reanudó al día siguiente. Es reanudable y no se perdió ninguna fila.

### Estadística de los 179

```
mínimo   −71,76   (CHEMBL2104842)
mediana  −29,56
media    −31,08
máximo   −10,55   (CHEMBL13045, cheleritrina)
dG_sd mediana 0,29   |   máximo 5,98
```

El cálculo es **estable**: la dispersión entre repeticiones tiene mediana 0,29
kcal/mol. Lo que sigue no es ruido.

Dispersión alta, a revisar por pose antes de creerles nada:
`CHEMBL1395394` (sd 5,53), `CHEMBL508159` (3,50), `CHEMBL2104842` (1,73),
`CHEMBL3921127` (1,65).

### Los diez primeros por energía

```
CHEMBL2104842  −71,76    CHEMBL2314732  −54,97    CHEMBL113229   −52,54
CHEMBL522746   −66,56    CHEMBL449782   −53,42    CHEMBL2334893  −51,06
CHEMBL3921127  −60,42    CHEMBL2334880  −53,02    CHEMBL1395394  −50,36
CHEMBL1739674  −58,01                               CHEMBL4589737  −49,22
CHEMBL2334892  −56,59                               CHEMBL508159   −48,68
CHEMBL169186   −55,17                               CHEMBL3903297  −48,37
```

**Estos números no se creen.** El sesgo de tamaño sigue intacto: pendiente
**−0,592 kcal/mol por átomo pesado**, correlación −0,35, con tamaños que van de
11 a 52 átomos pesados (mediana 25). Los de arriba son sencillamente los más
grandes. El −71,76 de CHEMBL2104842 es del mismo tipo que el enterramiento que ya
se vio en la corrida de los 240 (CHEMBL1094870, −87,35).

---

## 2. La prueba que faltaba: ¿dónde caen los controles?

Los ocho controles conocidos, ordenados por dG entre los 179:

| Control | Puesto | dG |
|---|---|---|
| Nitidina (cloruro / base) | 159 / 164 | −18,24 / −17,19 |
| Sanguinarina (base / cloruro) | 169 / 173 | −14,59 / −13,49 |
| Coptisina | 175 | −12,65 |
| Berberina (cloruro / sulfato) | 176 / 178 | −12,42 / −11,35 |
| Fagaridina | 177 | −11,46 |
| **Cheleritrina** | **179 (última)** | **−10,55** |

**Todos los activos conocidos de la familia caen en el fondo de la tabla.** Con la
caja vieja ya pasaba (percentil 5 por los tres caminos). Con la caja del ARN pasa
igual. Es la **segunda confirmación**: cambiar el sitio de acoplamiento no cambió
el orden.

---

## 3. Auditoría de la métrica que armó la lista

La lista de 179 se seleccionó por **similitud ECFP4 con los activos** (umbral 0,25),
que era la única medida que había dado AUC alto (0,885–0,923). La pregunta obligada
era si esa medida mide **que la molécula se une a TDP-43** o sólo **que la molécula
es de la misma familia química**. Lo segundo convertiría el AUC en una tautología.

Script: `analysis/auditar_similitud_tdp43.py` (sin GPU).
Salida: `rescoring_local/auditoria_similitud_tdp43.json`.
Fondo: el mismo aleatorio de 200 (semilla 20260917, 195 con SMILES).

### 3.1 Un solo compuesto da tanto como los ocho juntos

| Referencia | AUC |
|---|---|
| Los 8 activos | 0,885 |
| Sólo la familia plana (6) | 0,923 |
| **Sólo nitidina** | **0,952** |
| **Sólo sanguinarina** | **0,941** |
| **Sólo berberina** | **0,939** |
| Sólo PE859 | 0,766 |
| Sólo ketoconazol | 0,305 |

Una sola molécula de referencia bate al conjunto entero. No hay información de
"unión": hay **una etiqueta de familia**.

### 3.2 Dentro de la familia el AUC es 1,000 por construcción

| | familia sola | ketoconazol + PE859 |
|---|---|---|
| sim_max con los 8 | **1,000** | 0,538 |
| sim_max con la familia | **1,000** | 0,692 |

Que la métrica encuentre a los seis planos no es un resultado: es lo que la métrica
hace por definición. Para los dos activos que **no** son de la familia, está en el
azar.

### 3.3 Una regla de dos líneas empata con la huella de 2048 bits

| Medida | AUC |
|---|---|
| Huella ECFP4 con la familia plana | 0,923 |
| **Nº de anillos aromáticos + presencia de catión** | **0,921** |
| Sólo catión permanente | 0,837 |
| Sólo nº de anillos aromáticos | 0,771 |
| Fracción de átomos sp2 | 0,631 |
| Átomos pesados | 0,479 |
| Peso molecular | 0,436 |

Contar anillos y mirar si hay carga da lo mismo que la huella. La huella **no aporta
nada** que no estuviera ya en la etiqueta "familia".

### 3.4 Y no generaliza fuera de su familia

```
Referencia sólo familia plana      -> ¿encuentra ketoconazol y PE859?   AUC 0,692
Referencia sólo ketoconazol+PE859  -> ¿encuentra los 6 planos?          AUC 0,576
```

En una dirección roza el azar y en la otra es el azar. **La métrica es un detector
de familia, no un predictor de unión.**

---

## 4. La causa raíz: los ocho controles no son independientes

Procedencia documentada de los ocho (de `REVISION_VALIDACION_BOLSILLOS_2026-08-20.md`,
sección v6):

| Control | Respaldo real |
|---|---|
| PE859 | validado en HEK + *C. elegans* (paper 2026), unión directa |
| Berberrubina | validado en HEK + *C. elegans* (paper 2026), unión directa |
| Berberina, sanguinarina, epiberberina, coptisina, nitidina | **"cluster de berberrubina con actividad reportada"** — es decir, asumidos activos POR PARECERSE a la berberrubina |
| Ketoconazol | "reportado reducir agregación" — efecto funcional, no unión medida |

O sea:

- **Dos** ligandos con unión medida.
- **Seis** que entran en la lista por similitud con uno de ellos.
- **Uno** con efecto funcional.

**Validar una métrica de similitud contra un conjunto definido por similitud no es
validación.** Es la razón por la que el AUC sale 0,92: se está midiendo la métrica
contra sí misma. Y con dos ligandos reales de quimiotipo distinto no se puede
validar nada — dos puntos no hacen una curva.

Un dato independiente que confirma lo mismo: **de los 179, 89 son parecidos a la
berberrubina (sim ≥ 0,25) y 1 solo es parecido al PE859.** La lista es, en la
práctica, una serie de análogos de berberrubina.

Y el remate: **la correlación entre similitud con la berberrubina y el dG es +0,243**,
es decir **positiva** — cuanto más de la familia es la molécula, *peor* puntúa por
energía. Los dos criterios apuntan en direcciones opuestas. No se pueden combinar:
hay que elegir uno, y ninguno de los dos está validado.

---

## 5. Qué se sostiene y qué no

**Se sostiene**
- El cálculo: estable (sd mediana 0,29) y reproducible.
- El sesgo de tamaño del MM-GBSA en este receptor: **−0,592 kcal/mol por átomo
  pesado** (antes −1,26 en el conjunto de las tres dianas). Medición absoluta.
- Que el sitio de unión de ARN de TDP-43 **no se puede ordenar** con acoplamiento ni
  con MM-GBSA. Van dos cajas y dos veces el mismo resultado.
- Que los activos conocidos de la familia caen al fondo por energía (ya iba dos veces).

**No se sostiene**
- Que la similitud química sea una medida de unión en esta diana. Es un detector de
  familia y su AUC 0,92 es una tautología (secciones 3.1–3.4).
- Por tanto, que la lista de 179 sea una lista de **candidatos**. Es una serie de
  **análogos** de la familia de la berberrubina, y sirve para lo que sirve eso:
  ampliar la serie y buscar reposicionamiento dentro de ella.
- El orden por energía de esa lista, que además está **anti-correlacionado** con la
  química.

**Lo que queda por decidir, y no es un cálculo:** ni el acoplamiento ni el MM-GBSA
ordenan este sitio. Eso no es un fallo afinable con otra métrica; es la física de una
superficie plana y somera. La decisión es **cambiar de sitio o cambiar de estrategia**,
no afinar el ranking.

---

## 6. Propuesta, en orden

1. **Cerrar y no gastar más GPU en este sitio.** Ya está cerrado: no hay nada
   corriendo.
2. **Rehacer el conjunto de controles con procedencia y quimiotipos independientes.**
   Cada control con su cita y su tipo de ensayo (unión directa / funcional /
   asumido), y buscar activos de otras familias. Mientras n=2 reales, ningún ranking
   es evaluable.
3. **Fijar el criterio de aceptación por adelantado:** un ranking vale si coloca al
   PE859 y a la berberrubina por delante de señuelos emparejados en propiedades. Hoy
   **ningún sitio lo cumple** (v4 0,465; v6 0,517; MM-GBSA en las dos cajas, no).
4. **Elegir el sitio por ese criterio, no por comodidad.** Si el sitio de ARN no
   pasa la prueba con ligandos de dos quimiotipos distintos, se deja de usar para
   ordenar.
5. **Si se sigue con TDP-43, ir al núcleo de agregación** (dominio C-terminal,
   segmento 311–360; estructuras 5D2N/6B1Z): tiene estructura definida —láminas y
   superficie de crecimiento de fibra— y literatura propia de inhibidores de
   agregación. Ahí la geometría significa algo.
6. **Etiquetar el entregable como lo que es:** la lista de 179 es una serie de
   análogos de bencilisoquinolinas ordenada por similitud, con la energía sólo como
   desempate y nunca como criterio de elección. Y, como pedía la revisión del 20 ago,
   los candidatos de TDP-43 siguen **no validados** en el paper.

---

## 7. El arreglo, con nombre y cita

Está en **`analysis/PLAN_TDP43_2026-09-20.md`**: hay **más ligandos de TDP-43 con
unión medida de los que usamos, y de quimiotipos distintos** (rTRD01, Kd 89 µM por
MST+NMR; nTRD22, Kd 145 µM; tres fragmentos de RMN que se unen a un bolsillo
concreto de RRM2; AIM4; bis-ANS). Con siete positivos de cinco quimiotipos se puede
validar de verdad —y se puede ver caer la métrica de similitud, que es el resultado
honesto que hoy no somos capaces de obtener. Con eso, un fondo difícil de unidores de
ARN de R-BIND 2.0, y un sitio elegido por criterio puesto por adelantado.

## Archivos de este cierre

- `analysis/auditar_similitud_tdp43.py` — la auditoría de la métrica (nuevo)
- `analysis/rescoring_local/auditoria_similitud_tdp43.json` — su salida (nuevo)
- `analysis/rescoring_local/rescoring_focalizada_arn.csv` — los 179 (cerrado)
- `analysis/rescoring_local/runner_focalizada_arn.log` — el registro de la corrida
- `analysis/PLAN_TDP43_2026-09-20.md` — el plan de arreglo (nuevo)
- `analysis/rescoring_local/controles_tdp43_procedencia.csv` — procedencia de los 8 (nuevo)
- `analysis/INFORME_RECALCULO_2026-09-14.md` — el informe anterior, que esto continúa
