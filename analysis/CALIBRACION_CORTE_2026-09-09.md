# CALIBRACION DEL CORTE CON CONTROLES POSITIVOS — MASIVE-ALS

Generado: 2026-09-09 desde resultados en disco del GPU (RTX 4080).

## 1. Estado del cribado (outputs validos por proteina)

| Proteina | Outputs | Top-5%% (n) | Corte Vina (kcal/mol) | Mejor score |
|---|---|---|---|---|
| TDP43_v2 | 107608 | 5380 | -8.4 | -10.9 |
| SOD1 | 107608 | 5380 | -6.2 | -8.6 |
| FUS | 70863 | 3543 | -6.4 | -9.1 |

Excluidos documentados (ligandos imposibles, sin relanzar): **5696**.

## 2. Controles positivos conocidos vs corte

**TDP43**: 9 controles (2 en libreria cribada). Afinidades Vina: -5.7 a -6.3 kcal/mol.
- Corte top-5%: **-8.4** → pasarian el corte **0 de 2** controles.
  - berberine                -5.7
  - ketoconazole             -6.3

**SOD1**: 20 controles (18 en libreria cribada). Afinidades Vina: -4.7 a -6.0 kcal/mol.
- Corte top-5%: **-6.2** → pasarian el corte **0 de 18** controles.
  - CHEMBL1939222            -4.7
  - CHEMBL1643557            -4.7
  - CHEMBL2165612            -5.0
  - CHEMBL2165611            -5.0
  - CHEMBL2165603            -5.0
  - CHEMBL2165610            -5.1
  - CHEMBL2165607            -5.1
  - CHEMBL1643556            -5.2
  - CHEMBL2165605            -5.3
  - CHEMBL2165604            -5.3
  - CHEMBL2165608            -5.4
  - CHEMBL2165602            -5.4
  - CHEMBL1643541            -5.4
  - CHEMBL2165606            -5.5
  - CHEMBL2165609            -5.6
  - CHEMBL2165613            -5.8
  - CHEMBL2165601            -5.8
  - CHEMBL2165614            -6.0

**FUS**: 2 controles, 0 en libreria; ninguno con afinidad medible.

## 3. Lectura cientifica

- **SOD1** es la unica con calibracion util: 18/20 controles en libreria, afinidades **-4.7 a -6.0** kcal/mol. El corte top-5% en SOD1 es **-6.2**, por debajo (mas estricto) que el mejor control conocido. Ningun control conocido pasa el corte: el docking por si solo NO discrimina a los activos conocidos; el ranking necesita rescoring (MM-GBSA) + filtros de quimica medicinal (CNS, PAINS) antes de elegir candidatos.

- **TDP-43** (v2, bolsillo corregido): solo 2/9 controles en libreria (berberine -5.7, ketoconazole -6.3); los demas son sondas publicadas (PE859, sanguinarine...) ausentes de la libreria de farmacos. Corte top-5%: -8.4.

- **FUS**: 0/2 controles en libreria; calibracion pendiente de incluir ligandos FUS conocidos en futuras tandas. Corte top-5%: -6.4.

## 4. Recomendacion de corte (conservadora)

Con los datos actuales, el corte por percentil NO es comparable entre proteinas (los scores Vina dependen del receptor). Se recomienda:
1. Usar el ranking solo como filtro de enriquecimiento, no como valor absoluto.
2. Aplicar rescoring MM-GBSA sobre el top-5%% (14.303 compuestos).
3. Filtrar por CNS MPO >= 4 (o BBB clasica) — ya calculado en el CSV.
4. Validar los candidatos finales con controles positivos por proteina (ampliar la lista FUS y TDP-43 en libreria).

## 5. Archivos

- `ranking_afinidades.csv` — 286.079 afinidades Vina (las 3 proteinas)
- `ranking_top5_consolidado.csv` — top-5%% con SMILES, propiedades y filtros CNS
- `controles_calibracion.csv` — controles con afinidad real en libreria
- `_excluidos_clasificados.csv` — los 5.696 ligandos imposibles con motivo
