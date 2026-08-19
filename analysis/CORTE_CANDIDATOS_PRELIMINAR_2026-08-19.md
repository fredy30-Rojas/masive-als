# Corte de candidatos — PRELIMINAR

**Fecha:** 19/08/2026
**Estado:** ⚠️ PRELIMINAR — el cribado sigue corriendo (400 tandas, ahora en z006). Se re-evaluará a fin de agosto con todos los datos acumulados.

---

## Resumen del lote actual

- **42 candidatos** filtrados y balanceados (`candidatos_filtrados_v4.csv`)
- Distribución por proteína: **SOD1: 10 · FUS: 9 · TDP43: 23**
- Rango de afinidad: **-8.1 a -6.4 kcal/mol**
- **7 compuestos multi-diana** (atacan 2-3 proteínas)
- CNS MPO (penetración cerebral): 3.15 a 5.0
- Ninguno es fármaco ya aprobado

## Criterio de corte propuesto

1. **Afinidad ≤ -7.0 kcal/mol** → 14 candidatos prioritarios
2. Dentro de esos, priorizar:
   - **Multi-diana** (2-3 proteínas) — valor "pan-ELA"
   - **CNS MPO ≥ 4.0** — cruza la barrera hematoencefálica

## Los 14 prioritarios (afinidad ≤ -7.0)

| Ligando | Diana | Afinidad | CNS MPO | Multi-diana |
|---|---|---|---|---|
| CHEMBL3311449 | SOD1 | -8.1 | 3.52 | 2 |
| CHEMBL9010 | SOD1 | -7.4 | 3.55 | 3 |
| CHEMBL9532 | SOD1 | -7.4 | 3.15 | 3 |
| CHEMBL9532 | FUS | -7.4 | 3.15 | 3 |
| CHEMBL3310304 | SOD1 | -7.3 | 4.64 | — |
| CHEMBL7360 | SOD1 | -7.3 | 3.96 | 2 |
| CHEMBL3309775 | TDP43 | -7.2 | 4.64 | — |
| CHEMBL9347 | TDP43 | -7.2 | 4.43 | — |
| CHEMBL1163245 | SOD1 | -7.2 | 5.00 | — |
| CHEMBL8905 | SOD1 | -7.2 | 4.15 | 3 |
| CHEMBL601104 | TDP43 | -7.1 | 3.71 | — |
| CHEMBL3309990 | SOD1 | -7.1 | 5.00 | — |
| CHEMBL6900 | SOD1 | -7.1 | 4.61 | — |
| CHEMBL9440 | SOD1 | -7.1 | 4.13 | — |

## Élite — multi-diana (los más valiosos para ELA)

1. **CHEMBL8905** — ataca las **3 proteínas** + CNS MPO 4.15 (cruza cerebro). El más completo.
2. **CHEMBL9010** — 3 proteínas, SOD1 a -7.4.
3. **CHEMBL9532** — 3 proteínas, SOD1 y FUS a -7.4.
4. **CHEMBL3311449** — mejor afinidad global (-8.1), 2 proteínas.
5. **CHEMBL7360** — 2 proteínas, SOD1 a -7.3.

## Récord global aparte: CHEMBL1162102 (no en la lista CNS)

- **CHEMBL1162102** es el ligando con mejor afinidad global: **SOD1 -8.883**, TDP43 -7.274, FUS -7.212 → ataca las **3 proteínas**.
- **PERO** tiene 2 grupos ácido sulfónico (`S(=O)(=O)O`) en su estructura: es muy polar y con carga negativa a pH fisiológico, por lo que **no cruza la barrera hematoencefálica**. Por eso el filtro CNS lo dejó fuera de los 42 candidatos.
- Se conserva como **compuesto de referencia** (muestra qué prefiere el bolsillo de unión) y posible punto de partida para química médica (modificar los sulfónicos para mejorar penetración cerebral manteniendo el núcleo).

## Nota científica

- El CNS MPO usado es una aproximación de 4 componentes (no los 6 del método original), por lo que los valores **no se comparan directamente** con el punto de corte de 4.0 de la literatura.
- Corte definitivo: fin de agosto, con el total acumulado y validación con controles positivos.
