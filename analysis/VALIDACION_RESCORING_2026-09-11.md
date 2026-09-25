# VALIDACIÓN DEL RESCORING CORREGIDO — MASIVE-ALS

**Fecha:** 11 de septiembre de 2026
**Objeto:** comprobar si el rescoring MM-GBSA corregido (OBC2 + geometría de
la pose + receptor congelado + GPU determinista) **distingue compuestos
activos conocidos de compuestos cualesquiera**. Es la validación que la
auditoría de esta mañana señaló como imprescindible antes de ordenar
candidatos: ningún activo conocido había entrado en la lista corta, así que
no había forma de saber si el ranking servía.

**Ejecutado por:** Claude, a petición de Fredy.
**Artefactos:** `rescoring_validacion.csv` (200 filas), `runner_validacion.log`,
`receptores_fijos/snapshot_20260911/` (receptores congelados + `manifest.json`
con hashes).

---

## DISEÑO

Se reescalaron **20 controles positivos conocidos** (18 de SOD1 y 2 de
TDP-43, los únicos con pose acoplada) contra **180 compuestos de fondo**
tomados al azar de la propia lista corta, con la misma tubería y el mismo
receptor por diana. Si el score ordena bien, los controles deben quedar por
delante del fondo.

- Plataforma: **CUDA** (GPU), un solo receptor congelado por diana
  (md5 registrado en cada fila), 200/200 compuestos resueltos sin errores.
- Tiempo: **16,8 minutos** para 200 compuestos con 4 procesos (≈25 s cada uno).
- Comprobación previa: **todas** las poses (controles y fondo) caen dentro
  del bolsillo del receptor congelado, a distancias al centro comparables
  (controles 8,7 Å de media; fondo 9,5 Å en SOD1). La comparación es válida.

---

## RESULTADO

### SOD1 (18 controles vs 60 de fondo)

| Métrica | Valor | Lectura |
|---|---|---|
| dG mediano controles | −15,8 kcal/mol | — |
| dG mediano fondo | −21,0 kcal/mol | el fondo puntúa **mejor** |
| AUC bruto | **0,253** (p = 0,0016) | los controles quedan **por detrás**, significativamente |
| AUC con fondo de tamaño comparable (17-28 átomos, n=15) | 0,430 | sigue sin discriminar |
| **AUC sobre residuales (quitando el tamaño)** | **0,754** | aquí sí hay señal |
| Eficiencia de ligando (dG/átomo pesado) | controles −0,825 vs fondo −0,680 | AUC = 0,656 |
| Correlación dG vs tamaño | **−0,508** | el score premia el tamaño |

Por dG bruto, 10 de los 60 compuestos de fondo superan al mejor control. Por
eficiencia de ligando (normalizada por tamaño), solo **1 de 60**.

### TDP43_v2 (8 controles vs 60 de fondo) — RESULTADO NEGATIVO FIRME

Ampliado el 12 de septiembre: se acoplaron 6 controles más (PE859,
berberubina, sanguinarina, epiberberina, coptisina y nitidina) con la misma
tubería Vina-GPU, la misma caja y el mismo receptor, pasando de 2 a 8
controles. **Y con ocho controles el tamaño ya queda apareado** (controles
27,8 átomos pesados de media frente a 29,9 del fondo), así que el sesgo de
tamaño deja de explicar el resultado.

| Métrica | Valor | Lectura |
|---|---|---|
| dG mediano controles | −12,6 kcal/mol | — |
| dG mediano fondo | −31,1 kcal/mol | el fondo puntúa **mucho mejor** |
| AUC bruto | **0,140** (p = 0,0010) | los activos conocidos quedan por detrás |
| AUC con tamaño apareado (n=53) | **0,139** | **no es efecto del tamaño** |
| AUC sobre residuales | **0,073** | sigue sin discriminar |

En TDP-43 el rescoring **no solo no distingue a los activos conocidos: los
coloca sistemáticamente peor que compuestos cualesquiera**, y ahora está
descartado que sea por el tamaño molecular. Es un resultado negativo firme,
no una falta de datos.

### FUS

Sin controles con pose: no evaluable.

---

## CONCLUSIÓN

**El score corregido sigue dominado por el tamaño molecular**, aproximadamente
−1 kcal/mol por átomo pesado. Eso tiene una consecuencia práctica grave: en
SOD1 los activos conocidos son moléculas pequeñas, así que **ordenar los
candidatos por dG bruto coloca a los activos conocidos en la parte baja de la
tabla** (AUC 0,253, por debajo del azar, con p = 0,0016). No es que el
rescoring no aporte información: es que la información queda sepultada bajo el
sesgo de tamaño.

Cuando se quita el tamaño (residuales, o eficiencia de ligando), **aparece
señal moderada en SOD1** (AUC 0,65-0,75; solo 1 de 60 compuestos de fondo
supera al mejor control por eficiencia). Es un indicio, no una validación: con
18 controles el intervalo de confianza del AUC es amplio (±0,10).

En TDP-43 no hay señal: los dos controles quedan por detrás de la mitad del
fondo con cualquier normalización.

### Qué significa para el proyecto

1. **NO ordenar candidatos por dG bruto.** Un ranking por dG bruto de la
   lista corta es, en la práctica, un ranking por tamaño molecular, y
   perjudica justamente a los activos pequeños conocidos.
2. **La pasada completa de los 2.973 no debe lanzarse con el score bruto**:
   serían ~6 h de GPU para producir un orden que ya sabemos sesgado.
3. Si se usa el rescoring, usar la métrica **normalizada por tamaño**
   (residual o eficiencia de ligando) y declararlo explícitamente.
4. **Ampliar el conjunto de validación**, sobre todo en TDP-43 (solo 2
   activos): sin más controles no se puede afirmar ni negar capacidad de
   discriminación.
5. Pendiente técnico menor: la reproducibilidad medida entre corridas
   idénticas es de ~0,1-0,4 kcal/mol, no 0,00 como decía el informe previo
   (ver nota de reproducibilidad abajo).

---

## PENDIENTES QUE ESTO ABRE

- Reunir más activos conocidos de TDP-43 (literatura) y acoplarlos, para
  tener un conjunto de validación con potencia estadística.
- Decidir la métrica final (residual vs eficiencia) y congelarla antes de
  cualquier ranking.
- Añadir señuelos apareados por propiedades (tipo DUD-E) para medir
  enriquecimiento de verdad, no solo separación contra fondo aleatorio.
## REPRODUCIBILIDAD — MEDIDA (corrige al informe del 12:24)

El informe `RESCORING_WSL_CUDA_2026-09-11.md` afirmaba dispersión 0,00 con
receptor congelado. Medido hoy con el protocolo corregido y el receptor del
snapshot, **tres corridas idénticas del mismo compuesto** (berberina contra
TDP43_v2) dan:

| Plataforma | dG de las 3 corridas idénticas | Dispersión | Tiempo |
|---|---|---|---|
| CUDA, precisión mixta (float32) | −8,43 · −8,05 · −6,92 | **1,51 kcal/mol** | ~15 s |
| CUDA, precisión doble | −14,62 · −14,25 · −14,49 | **0,37 kcal/mol** | ~20 s |
| CPU (float64, 1 hilo) | −8,48 · −8,48 · −8,48 | **0,00 (bit a bit)** | ~2-4 min |

Dos consecuencias:

1. **La precisión simple no solo añade ruido: sesga el valor.** El mismo
   compuesto da ≈ −8 kcal/mol en float32 y ≈ −14,5 en float64 sobre CUDA.
2. **Solo la vía CPU (float64, un hilo) es reproducible bit a bit.** La GPU,
   incluso en doble precisión, deja 0,37 kcal/mol de dispersión; y entre
   plataformas (CPU vs CUDA) el mismo cálculo difiere varios kcal/mol,
   porque la minimización acaba en mínimos distintos.

**Alcance sobre esta validación:** la corrida de los 200 se hizo entera en
CUDA float32, así que cada valor lleva ±1,5 kcal/mol de ruido. El veredicto
(el sesgo de tamaño domina; ordenar por dG bruto entierra a los activos
conocidos) es robusto porque los efectos medidos son de decenas de kcal/mol.
Los números finos (AUC 0,754 residual) sí llevan esa incertidumbre encima y
deben rehacerse con el protocolo corregido antes de citarlos.

**Protocolo que se deriva:** una sola plataforma para todo el conjunto, con
precisión doble, y N repeticiones por compuesto con su dispersión declarada
(media ± error). Nunca comparar valores medidos en plataformas distintas.
La opción CPU es la única determinista pero cuesta ~2-4 min por compuesto.

---

**Reproducible con:** `runner_mmgbsa_gb.py` (lanzador), `preparar_validacion.py`
(lista), `analizar_validacion.py` (análisis), `_chequeo_bolsillo.py`
(comprobación de encaje en el bolsillo) y el snapshot de receptores con su
manifiesto de hashes.
