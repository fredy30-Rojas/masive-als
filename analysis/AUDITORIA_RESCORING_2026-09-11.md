# AUDITORÍA CIENTÍFICA DEL RESCORING MM-GBSA — MASIVE-ALS

**Fecha:** 11 de septiembre de 2026
**Objeto:** determinar si los datos de rescoring MM-GBSA (corrida local sobre la lista corta de 2.973 candidatos) son válidos y utilizables.
**Datos auditados:**
- `analysis/rescoring_local/rescoring_lista_corta.csv` (2.973 filas; corrida local 10/09 23:29 → 11/09 07:29)
- `rescoring_datos/*.csv` (corridas previas, 99 ligandos con dG)
- `analysis/rescoring_mmgbsa_robusto.py` + `runner_rescoring_v5.py` (pipeline de Oracle, con OBC2)
- `analysis/lista_corta_candidatos.csv` (2.973 candidatos: TDP43 1.237, SOD1 1.167, FUS 569)
- `analysis/controles_calibracion.csv` (20 controles positivos conocidos)

---

## VEREDICTO

**Los dG del rescoring local NO son utilizables para clasificar ni para decidir candidatos.**

La infraestructura funciona (91,9 % de éxito), pero el número que produce **no es una energía libre de unión**: es una diferencia de energía potencial **en vacío**, no reproducible y sin capacidad de discriminación. Hay 7 hallazgos, todos verificados sobre los datos.

---

## HALLAZGO 1 — No es MM-GBSA: falta el disolvente implícito (causa raíz)

El pipeline local `analysis/rescoring_local/mmgbsa_openff.py` construye los tres sistemas (complejo, receptor, ligando) con SMIRNOFF y los evalúa con **`NonbondedForce.NoCutoff`**, sin ningún modelo de disolvente: no hay `implicit/obc2.xml`, ni `CustomGBForce`, ni término de Born generalizado. **La GBSA no existe en el código.**

En contraste, el pipeline de Oracle sí lo tiene:

```
rescoring_mmgbsa_robusto.py:301    forcefield = ForceField("amber14-all.xml", "implicit/obc2.xml")
```

Consecuencia: sin apantallamiento dieléctrico, la interacción Coulombiana ligando–bolsillo en vacío domina y da valores enormes. Los dos pipelines **no miden la misma magnitud**, por lo que sus números no son comparables (ver hallazgo 6).

---

## HALLAZGO 2 — Magnitudes fuera de rango físico

| Métrica | Valor observado | Rango realista |
|---|---|---|
| dG mediano | **−60,8 kcal/mol** | −5 a −40 (MM-GBSA) |
| dG rango | −13,9 a −137,7 | — |
| Eficiencia de ligando (LE) mediana | **2,16 kcal/mol·átomo** | 0,2 – 0,5 |
| Ligandos con LE > 1,5 | **2.410 (88,2 %)** | ~0 % |
| LE máxima | **4,71** | ~0,6 extremo |

La eficiencia de ligando es ~5–10 veces mayor de lo físicamente posible en cualquier compuesto de esta librería. Un dG de −131,8 kcal/mol para una molécula de 28 átomos (LE 4,7) no corresponde a ninguna unión real.

**Nota:** 1.791 de 2.731 filas (65,6 %) tienen `e_ligand > 0` (mediana +25,3 kcal/mol). No es un error de resta — la aritmética interna `dG = E_complejo − E_receptor − E_ligando` es correcta en las 2.731 filas — pero confirma que se están evaluando energías de vacío, no un ciclo termodinámico solvatado.

---

## HALLAZGO 3 — Cero capacidad de discriminación frente al docking

Correlación de rangos entre afinidad Vina y dG MM-GBSA:

| Target | n | Spearman ρ |
|---|---|---|
| FUS | 493 | **0,018** |
| SOD1 | 1.073 | **0,144** |
| TDP43_v2 | 1.165 | **0,052** |

El rescoring no reproduce ni siquiera el orden del docking. Si el reescalado fuese informativo, ρ debería estar en 0,3–0,6. Valores de 0,02–0,14 significan que **el reescalado es esencialmente ruido respecto al cribado**.

---

## HALLAZGO 4 — El “score” premia el tamaño molecular, no el encaje

Correlación del dG con propiedades del ligando:

| Target | ρ(dG, átomos pesados) | ρ(dG, TPSA) |
|---|---|---|
| SOD1 | **−0,353** | −0,369 |
| TDP43_v2 | **−0,511** | −0,507 |
| FUS | −0,055 | −0,114 |

En TDP-43, media varianza del orden lo explica solo el tamaño de la molécula. Es el artefacto clásico de una función de energía sin desolvación: cuanto más grande el ligando, más contactos sin apantallar y más negativo el número. **Clasificar con esto favorece a las moléculas grandes, no a las que unen bien.**

---

## HALLAZGO 5 — No reproducible: el error numérico ≈ el rango de señal

Mismo compuesto (CHEMBL4584906), misma pose, mismo receptor, repetido:

| Corrida | Plataforma | dG |
|---|---|---|
| Guardado en CSV | OpenCL | −66,74 |
| Repetición 1 | OpenCL | −67,33 |
| Repetición 2 | OpenCL | −76,75 |
| Repetición 3 | OpenCL | −66,48 |
| Repetición 4 | OpenCL | −64,73 |
| Repetición 5 | OpenCL | −66,72 |
| **Repetición en CPU (Reference)** | **Reference** | **−78,03** |

- Ruido entre corridas idénticas (OpenCL): **≈ 2 kcal/mol**
- Diferencia entre plataformas (OpenCL vs Reference): **≈ 11,5 kcal/mol**

El rango total de los 2.731 dG es −13,9 a −137,7 con **desviación típica 17,9** e IQR de 24 kcal/mol. Con un error de 11–12 kcal/mol, la **relación señal/ruido es ≈ 1,5**: el instrumento no distingue dos compuestos cualesquiera de la lista.

Causa técnica identificable: el ligando se **re-embebe** (`ETKDGv3`, `useRandomCoords=True`, `CoordMap` de la pose) en lugar de conservar las coordenadas acopladas, y la minimización del complejo se hace en plataforma GPU con reducciones no deterministas. La geometría final —y con ella E_receptor y E_complejo— depende de la semilla y del hardware.

---

## HALLAZGO 6 — Inconsistencia con el rescoring previo (OBC2)

Solo hay 3 ligandos comparables entre la corrida local y las corridas previas:

| Ligando | Target | Oracle (OBC2) | Local (vacío) | Δ |
|---|---|---|---|---|
| CHEMBL4530588 | TDP43 | −22,56 | −67,99 | **−45,4** |
| CHEMBL564972 | TDP43 | −26,62 | −54,31 | **−27,7** |
| CHEMBL4634995 | TDP43 | −29,37 | −55,48 | **−26,1** |

Diferencias de 26–45 kcal/mol confirman que son magnitudes distintas (una solvatada, otra no) y que ninguno de los dos conjuntos puede mezclarse con el otro.

---

## HALLAZGO 7 — Sin controles positivos: el rescoring no es validable con estos datos

Los 20 controles positivos conocidos (18 SOD1 + berberina y ketoconazol para TDP-43) tienen afinidades Vina de **−4,7 a −6,3** y **todos quedan fuera del corte top-5 %** (−6,2 en SOD1, −8,4 en TDP-43). Por tanto **ningún activo conocido entró en la lista corta ni fue reescalado**.

Sin un solo activo conocido dentro del conjunto reescalado, **no se puede medir enriquecimiento** (AUC, factor de enriquecimiento, top-N) y no hay ninguna evidencia de que el orden del rescoring tenga relación con la afinidad real. Este es el argumento decisivo: no es que el rescoring funcione mal, es que **no hay un único dato que demuestre que funciona bien**.

---

## QUÉ SÍ ES VÁLIDO Y REUTILIZABLE

1. **Toda la infraestructura**: 2.973 candidatos procesados, 91,9 % de éxito, checkpoint reanudable, 242 fallos documentados (embeddings y receptores no estándar de FUS).
2. **El cribado Vina y el embudo de filtros** (top-5 %, PAINS, Lipinski/Veber, CNS MPO): no dependen de este rescoring.
3. **La lista corta de 2.973 candidatos** sigue siendo la entrada correcta para un rescoring bien hecho.
4. **La calibración con controles positivos** (informe `CALIBRACION_CORTE_2026-09-09.md`): sigue siendo el hallazgo más sólido del proyecto.

## QUÉ NO DEBE USARSE

- Los valores `mmgbsa_dG` de `rescoring_lista_corta.csv` para **ordenar, filtrar o elegir** candidatos.
- Mezclar estos dG con los de `rescoring_mmgbsa_v5.csv` / `rescoring_tdp43_mmgbsa.csv` (magnitudes distintas).
- Cualquier cifra derivada: eficiencias de ligando, “mejores candidatos”, figuras del paper basadas en este ranking.

---

## REMEDIOS CONCRETOS (en orden)

1. **Añadir el término GBSA.** Opciones: (a) ejecutar el rescoring en la máquina Linux con AMBER/GAFF + `implicit/obc2.xml` (el pipeline de Oracle ya está probado allí y es el camino corto); (b) en Windows, añadir a mano un `CustomGBForce` con radios OBC2 sobre el sistema OpenFF.
2. **Usar la geometría acoplada, no re-embebida.** Transferir las coordenadas de la pose y minimizar solo el complejo: elimina la fuente principal de ruido.
3. **Fijar una plataforma determinista** (Reference/CPU o `CUDA` con `deterministicForces=True`) y reportar la plataforma en el CSV.
4. **Validar antes de escalar.** Incluir en el conjunto a reescalar los 20 controles positivos (con sus poses) y medir enriquecimiento. Sin eso, no se puede afirmar que el rescoring aporta.
5. **Reportar siempre la incertidumbre.** Publicar el dG con su error (± 2 kcal/mol intra-plataforma) y no ordenar diferencias menores que ese error.

---

## CONCLUSIÓN

El rescoring realizado **no produce datos usables**. La causa es doble y verificable: **no hay modelo de disolvente** (por lo que la cifra no es una energía libre de unión) y **el cálculo no es reproducible** (el error supera la mitad del rango de señal). La consecuencia científica es que el ranking actual **no está validado ni es validable** con los datos disponibles, porque ningún activo conocido entró en el conjunto reescalado.

Nada de esto invalida el cribado ni el embudo de filtros: invalida el paso de rescoring tal como se ejecutó. Con los cinco remedios anteriores, la lista corta de 2.973 candidatos sigue siendo una base sólida para obtener un ranking defendible.

**Auditoría reproducible:** `analysis/_auditoria_rescoring.py`, `_auditoria_rescoring2.py`, `_auditoria_rescoring3.py`.
