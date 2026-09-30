# Auditoría de la serie de hCA2 — 30 de septiembre de 2026

Antes de acoplar un solo ligando. Con la receta de `calibracion_dude_cdk2/auditar_activos_cdk2.py`, que es la que demostró que el AUC de CDK2 era interpretable.

## Cifras

| medida | valor |
|---|---|
| diana | Anhidrasa carbonica 2 (hCA2) (CHEMBL205) |
| activos medidos auditados | 1200 |
| esqueletos de Murcko distintos | 281 |
| familias de 5+ activos | 44 |
| esqueleto mas grande real | 45 activos (3.8 %) |
| MCS sobre 80 activos | 4 atomos |
| MCS dentro de las 203 sulfonamidas | 4 atomos |
| con sulfonamida | 203 (16.9 %) |
| Tanimoto mediano entre activos (Morgan con quiralidad, conteos) | 0.185 |
| Tanimoto maximo | 0.935 |
| compuestos con dos entradas de ChEMBL (mismo canonico) | 0 |

## Los diez esqueletos más poblados

| activos | %% | esqueleto |
|---|---|---|
| 314 | 26.2 | `c1ccccc1` |
| 45 | 3.8 | `O=C(CCC[C@H]1CC[C@@H]2C1CC[C@@H]1C3CCCC[C@H]3CC[C@H]12)Nc1ccccc1` |
| 41 | 3.4 | `O=S(=O)(Nc1nncs1)c1ccccc1` |
| 30 | 2.5 | `O=C(C[n+]1ccc(-c2ccccc2)cc1)Nc1ccccc1` |
| 28 | 2.3 | `c1cc2ccsc2s1` |
| 25 | 2.1 | `c1ccc2scnc2c1` |
| 24 | 2.0 | `O=S(=O)(Nc1ccccc1)c1ccccc1` |
| 19 | 1.6 | `O=S(=O)(/N=c1\[nH]ncs1)c1ccccc1` |
| 19 | 1.6 | `O=C(NCCc1ccccc1)c1ccccc1` |
| 18 | 1.5 | `c1nncs1` |

## Veredicto

- LIMPIA por los tres criterios

## Qué significa para el AUC

El número grande de hCA2 (13.051 activos medidos en total) se reduce aquí a **44 independientes**, y ese es el que limita el intervalo del AUC, no los 1200 brutos. Con 44 independientes el semiancho del IC 95 % anda alrededor de 0,1 si el AUC sale cerca de 0,5, y baja de 0,06 si sale alto. Sigue siendo un orden de magnitud mejor que TDP-43 (7 positivos) y SOD1 (11).

## Ficheros

- `analysis/auditar_serie_hca2.py` — la auditoría.
- `analysis/serie_hca2_auditoria.csv` — los 1200 activos con su esqueleto.
- `analysis/auditar_serie_hca2.log` — la salida completa.
