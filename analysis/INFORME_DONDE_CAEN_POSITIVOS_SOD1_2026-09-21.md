# Dónde caen de verdad los positivos de SOD1: se retira el «no comparten modo de unión»

**21 de septiembre de 2026.**
**Herramientas:** `donde_caen_positivos_sod1.py`, `verdad_de_referencia.py`,
`preparar_ligando.py`, `perfil_interaccion.py`
**Datos:** `donde_caen_positivos_sod1.csv`, `verdad_de_referencia.csv`

Este informe **retira una conclusión del 20 de septiembre** y la sustituye por la
medida correcta. También cierra la causa raíz del fallo «glue» y deja armada la
verdad de referencia.

---

## 1. Lo que estaba bloqueando todo

El 20 de septiembre se midió que **LCS-1 y PRG-A01 —los dos únicos positivos
independientes de SOD1— comparten 0,038 de los residuos que contactan**, menos que
dos señuelos al azar (0,297). De ahí salió la conclusión de
`RESULTADO_PERFIL_INTERACCION_2026-09-20.md`:

> «En esta caja no hay un modo de unión común entre los ligandos independientes…
> ninguna métrica basada en las poses puede funcionar aquí todavía.»

Esa medida se hizo con el **receptor histórico `gpu_dock/SOD1.pdbqt`**, que no es una
proteína limpia: lleva las 1.763 aguas del cristal, el zinc y **las 18 copias** de
1HL5 (21.585 átomos). Y con las poses acopladas contra ese receptor.

Antes de tocar ninguna métrica había que responder una pregunta de hecho: ¿dónde se
coloca cada positivo?

## 2. El patrón de medida: la cristalografía

Cuatro de los siete positivos de SOD1 están **co-cristalizados en Trp32**, así que su
sitio está resuelto y no inferido. Sus contactos, leídos de su propia estructura:

| Ligando | Estructura | Residuos del bolsillo que toca |
|---|---|---|
| 5-fluorouridina (5UD) | 4A7S | 21, 30, **32**, 97, 98 |
| isoprenalina (5FW) | 4A7T | 21, 22, 23, 28, 29, 30, **32**, 100 |
| epinefrina (ALE) | 4A7U | 21, 22, 23, 28, 29, 30, **32**, 100 |
| dopamina (LDP) | 4A7V | 21, 22, 23, 28, 29, 30, 100 |

Consenso (residuos que toca al menos la mitad): **21, 22, 23, 28, 29, 30, 32, 100**.

## 3. El resultado

Acoplados los siete con la misma caja Trp32 (centro 46,5/80,0/73,3; 22 Å;
exhaustividad 8) contra los dos receptores:

| Ligando | Receptor limpio | Receptor histórico |
|---|---|---|
| **LCS-1** | **−5,00  EN Trp32** (7 contactos) | −4,93  EN Trp32 (5) |
| **PRG-A01** | **−7,46  EN Trp32** (6) | −5,56  EN Trp32 (5) |
| CHEMBL2165613 (pirazolona) | −6,53  EN Trp32 (9) | −6,71  EN Trp32 (7) |
| 5-fluorouridina (co-crist.) | −4,70  EN Trp32 (6) | −5,02  **no toca Trp32** |
| isoproterenol (co-crist.) | −4,44  EN Trp32 (6) | −4,92  **no toca Trp32** |
| dopamina (co-crist.) | −4,57  EN Trp32 (4) | −4,73  EN Trp32 (7) |
| epinefrina (co-crist.) | −4,51  EN Trp32 (4) | −4,57  EN Trp32 (7) |

**Con el receptor corregido, los siete caen en Trp32, y los cuatro co-cristalizados
lo reproducen (4 de 4).** Con el receptor histórico solo lo hacen dos de los cuatro:
el receptor contaminado es el que no reproduce el sitio documentado.

| | Receptor limpio | Receptor histórico |
|---|---|---|
| Positivos que tocan Trp32 | **7 de 7** | 5 de 7 |
| Co-cristalizados que lo reproducen | **4 de 4** | 2 de 4 |
| Residuo común a los siete | **32 (Trp32)**, y solo ese | — |
| **LCS-1 vs PRG-A01 (Jaccard)** | **0,444** | — |

Antes: 0,038. Ahora: **0,444**, por encima del azar entre señuelos (0,297) y en el
mismo orden que la serie pirazolona (0,397).

### Consecuencia

**La conclusión del 20 de septiembre era un artefacto del receptor, y se retira.**
Sí hay un modo de unión común entre los ligandos independientes: **Trp32**, el único
residuo que tocan los siete. Lo que no había era un receptor capaz de encontrarlo.

## 4. La huella de interacción, re-testeada con el receptor limpio

Con las 513 poses de `validacion_SOD1_v5/out` (todas acopladas contra el receptor
limpio) y con `perfil_interaccion.py`:

| Referencia del consenso | AUC por huella (recall) | Antes (receptor histórico) |
|---|---|---|
| los 31 activos de la validación | 0,590 (Vina: 0,530) | — |
| **solo LCS-1 y PRG-A01** | **0,566** | **0,000** |

El 0,000 catastrófico desaparece. Pero la huella **sigue sin batir a Vina**, y su
defecto estructural no ha cambiado: correlaciona **+0,42** con el número de contactos
(el mismo vicio de tamaño que la energía, que correlaciona −0,21). LCS-1 queda en el
puesto 378 de 513 y PRG-A01 en el 63.

Es decir: la métrica no está rota por lo que creíamos, pero **tampoco sirve todavía**,
y ahora se sabe por qué. Reintentarla exige la validación rehecha contra el receptor
limpio y con la verdad de referencia arreglada, no con las poses actuales.

## 5. Causa raíz del fallo «glue», cerrada

El fallo no volverá a producirse por copiar la receta equivocada: **la preparación de
ligandos vive ahora en un solo sitio**, `masive-als/preparar_ligando.py`, que

1. cierra los anillos grandes (`rigid_macrocycles=True`),
2. incrusta con seis semillas en vez de una (con una sola, algunos macrociclos no
   incrustaban y el ligando se caía del conjunto en silencio),
3. **comprueba que no queden pseudo-átomos y que el número de átomos cuadre con el
   SMILES**, y devuelve el motivo en vez de escribir un fichero malo,
4. desalar cuando se le pide (las tres catecolaminas venían con su clorhidrato o
   tartrato; meeko rechaza moléculas de dos fragmentos y se caían del conjunto).

Los cuatro scripts que preparan ligandos delegan ya en ella (`analysis/validar_senuelos.py`,
`gpu_dock/_preparar_libreria_gpu.py`, `gpu_dock/acoplar_controles_tdp43.py`,
`analysis/reparar_ligandos_glue.py`). Y `sanear_tipos` de `acoplar_controles_tdp43.py`,
que traducía `CG0` a carbono y **metía dos carbonos fantasma en el ligando acoplado**,
ahora **aborta** en vez de traducir.

Prueba de regresión incluida: `python preparar_ligando.py` re-prepara los 47 ficheros
reales de los `_roto_glue/` y comprueba que ninguno queda con pseudo-átomos. **47 casos,
0 fallos.**

## 6. Verdad de referencia

`verdad_de_referencia.py` consolida los positivos de las tres dianas con **cita y tipo
de ensayo**, colapsando las series. La regla, escrita por adelantado: sin cita y sin
tipo de ensayo, el compuesto no entra.

| Diana | Entradas | Aptos (unión medida) | Quimias distintas |
|---|---|---|---|
| SOD1 | 14 | **11** | 8 |
| TDP-43 | 17 | **7** | 5 |
| FUS | 2 | **0** | 0 |

- **SOD1**: los 11 aptos son co-cristalizados en Trp32. LCS-1, PRG-A01 y la
  pirazolona quedan marcados **no aptos** — su respaldo es actividad celular, no unión
  medida, y la pirazolona es el representante único de 18 análogos colapsados (no es
  evidencia de generalidad).
- **TDP-43**: 7 aptos de 5 quimiotipos (rTRD01, nTRD22, tres fragmentos de RRM2,
  PE859, berberrubina). Los seis que entraron por parecerse a la berberrubina quedan
  fuera con ese motivo escrito. Faltan los SMILES de los tres fragmentos: hay que
  dibujarlos desde la fuente (Nshogoza 2019).
- **FUS**: no tiene ligandos propios en la literatura. No es evaluable, y así queda
  escrito.

### Control de redocking de los once co-cristalizados

Medido con el mismo protocolo (caja 24 Å, exhaustividad 8, tres semillas, 27 modos),
uniforme para los once:

| Entrada | Ligando | Mejor de los 27 modos | El que gana por energía | Modos ≤ 2 Å | Veredicto |
|---|---|---|---|---|---|
| 4A7T | isoproterenol | **0,48 Å** | 0,48 Å | 6/27 | RECUPERA |
| 4A7U | adrenalina | **0,69 Å** | 0,72 Å | 4/27 | RECUPERA |
| 4A7V | dopamina | **0,74 Å** | 3,15 Å | 4/27 | RECUPERA |
| 2WZ6 | ZO0 | **1,43 Å** | 5,77 Å | 1/27 | RECUPERA |
| 5YTO | 946 | 2,43 Å | 4,89 Å | 0/27 | FALLA |
| 4A7G | 12I | 2,70 Å | 5,23 Å | 0/27 | FALLA |
| 4A7Q | 4MQ | 3,15 Å | 10,55 Å | 0/27 | FALLA |
| 8GSQ | K4I | 3,25 Å | 6,88 Å | 0/27 | FALLA |
| 2WZ0 | ZZT | 3,60 Å | 4,36 Å | 0/27 | FALLA |
| 6A9O | 6B3 | 4,75 Å | 5,03 Å | 0/27 | FALLA |
| 4A7S | 5-fluorouridina | 9,56 Å | 11,07 Å | 0/27 | FALLA |

Dos hechos que importan:

1. **4 de 11 recuperan la pose; 7 no.** El sitio está bien, el acoplador no siempre lo
   encuentra.
2. **En tres de los cuatro que recuperan (dopamina, ZO0, y casi adrenalina), el modo que
   gana por energía NO es la pose cristalina** (0,74 frente a 3,15 Å; 1,43 frente a
   5,77 Å). Es decir: el problema vuelve a ser de **puntuación, no de muestreo** — el
   muestreo genera la pose correcta y la función prefiere otra.
3. La pose depositada de la **5-fluorouridina no es un mínimo del potencial** (contacto
   F···N de 2,07 Å en una estructura de 1,06 Å de resolución). No sirve como patrón de
   medida y está marcado así.

## 7. Qué sigue, en orden

1. **Rehacer la validación de SOD1 con el receptor limpio** y con la verdad de
   referencia corregida: 11 positivos de unión medida, 8 quimias, más los señuelos
   emparejados y el fondo duro. Es lo que el criterio exige antes de reportar nada.
2. **Dibujar los tres fragmentos de RRM2** y completar TDP-43 con un fondo difícil
   (R-BIND 2.0), que es el arreglo 2 de `PLAN_TDP43_2026-09-20.md`.
3. **Solo después**, reintentar una métrica de colocación — ahora con receptores limpios
   y positivos independientes de verdad.
4. El bolsillo de Trp32 ya no está en duda; **la caja no hay que rehacerla**.
