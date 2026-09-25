# REORDENAMIENTO POR TAMAÑO DE LA LISTA CORTA — MASIVE-ALS

**Fecha:** 13 de septiembre de 2026
**Ejecutado por:** Claude, a petición de Fredy
**Script:** `analysis/rescoring_local/reordenar_normalizado.py`
**Salidas:** `rescoring_local/rescoring_ranking_normalizado.csv` (2.731 filas) y
`rescoring_local/ranking_normalizado_resumen.txt`

---

## POR QUÉ SE HIZO

La validación del 11 de septiembre demostró que el dG bruto del MM-GBSA está
dominado por el tamaño molecular: ordenar por energía bruta es, en la práctica,
ordenar por número de átomos, y entierra a los activos pequeños conocidos
(AUC 0,253 en SOD1 y 0,140 en TDP43_v2, los dos **por debajo del azar**).

Este trabajo reordena los **2.731** resultados del rescoring con las métricas
normalizadas que la validación señaló como las únicas con señal:

- **residual** = dG − (a + b · átomos_pesados), con la recta ajustada por diana
  (ajuste con recorte de atípicos a 3 MAD y reajuste). Residual negativo =
  mejor de lo que su tamaño predice.
- **eficiencia de ligando (LE)** = dG / átomos_pesados.

Los átomos pesados se cuentan de la pose acoplada (mismo criterio que
`analizar_validacion.py`), no del SMILES: 2.731 de 2.731 con pose encontrada.

---

## QUÉ SALE

| Diana | n | pendiente dG/tamaño | Spearman dG~tamaño | top-100 bruto vs residual |
|---|---|---|---|---|
| SOD1 | 1.073 | −2,11 kcal/mol por átomo | −0,403 | 88 de 100 |
| TDP43_v2 | 1.165 | −2,02 | −0,571 | **71 de 100** |
| FUS | 493 | −0,34 | −0,069 | 96 de 100 |

En total, **761 de 2.731 candidatos suben más de 50 puestos** al quitar el
tamaño. El reordenamiento cambia poco en FUS, algo en SOD1 y bastante en
TDP43_v2.

### Candidatos que suben en más de una diana (top-50 por residual)

| Ligando | Dianas y puesto |
|---|---|
| CHEMBL521548 | SOD1 p5 · TDP43_v2 p6 |
| CHEMBL4566549 | FUS p4 · TDP43_v2 p41 |
| CHEMBL4573275 | SOD1 p14 · FUS p21 |
| CHEMBL1204619 | SOD1 p25 · FUS p16 |

### Líderes por eficiencia de ligando (dG por átomo pesado)

- **SOD1:** CHEMBL19179 (27 átomos, LE −4,13; puesto bruto 17 → 1 por LE)
- **TDP43_v2:** CHEMBL12225 (30 átomos, LE −4,39)
- **FUS:** CHEMBL28751 (27 átomos, LE −3,29)

Los que más suben en TDP43_v2: CHEMBL519419 (del puesto 104 al 35 por residual),
ESTRONE SULFURIC ACID (86 → 14) y CHEMBL515831 (47 → 10). Es justo lo esperado:
moléculas pequeñas que el score bruto castigaba por su tamaño.

---

## ADVERTENCIA GRAVE — ESTOS NÚMEROS NO SON DEL PROTOCOLO VALIDADO

Al revisar los resultados apareció un problema de fondo que hay que decir sin
adornos:

**Los 2.731 cálculos se hicieron con los receptores viejos de `gpu_dock/`, no
con los receptores congelados de la validación** (creados el 11 sep a las
14:44-14:58, cuando la lista corta ya había terminado a las 07:29).

| Receptor | Archivo viejo (lista corta) | Congelado (validación) |
|---|---|---|
| SOD1 | **2.738 residuos / 21.585 átomos** | 98 residuos / 1.459 átomos |
| TDP43_v2 | 174 residuos / 1.401 átomos | 128 residuos / 2.067 átomos |
| FUS | 41 residuos / 514 átomos | 41 residuos / 601 átomos |

SOD1 tiene unas **18 copias** del monómero (153 residuos): es el receptor sucio
que ya se había detectado. Por eso las energías de la lista corta van de
−80 a −137 kcal/mol, mientras que los controles validados dan de −8 a −16.

Consecuencias:

1. El reordenamiento es **coherente dentro de cada diana** (todas las filas de
   una diana usan el mismo receptor), así que la corrección de tamaño es válida
   como ejercicio interno.
2. Pero **no sabemos si ese score discrimina nada**: la validación de
   discriminación (AUC) se hizo con el receptor congelado, no con este.
3. Por tanto, **este ranking NO debe usarse para decidir candidatos**. Sirve
   para ver a quién favorece quitar el tamaño, y para elegir a quién recalcular.

---

## SIGUIENTE PASO RECOMENDADO

Recalcular con el **protocolo validado** (receptores de
`receptores_fijos/snapshot_20260911/`, doble precisión, 3 repeticiones por
compuesto, plataforma única) un subconjunto: los mejores por residual y por
eficiencia de cada diana.

Coste medido: 84 s por compuesto con 3 repeticiones en CUDA. Con 4 procesos:

- 100 compuestos ≈ 35 minutos
- 300 compuestos ≈ 1 hora 45 minutos

Con eso se tendría un ranking que sí se puede defender, con su dispersión
declarada y sin comparar plataformas distintas.
