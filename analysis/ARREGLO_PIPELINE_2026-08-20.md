# Arreglo del pipeline de validación y rescoring — 20 agosto 2026

Revisión y correcciones aplicadas por Salomé (Buffy) sobre los tres problemas
pendientes de la auditoría científica del 18/08.

## 1. Re-validación de bolsillos con señuelos (7× más)

Se repitió el protocolo decoy con 6.462 señuelos property-matched (antes 73–33).

| Diana | Caja | Activos | Señuelos | ROC-AUC | EF5% | Veredicto |
|---|---|---|---|---|---|---|
| SOD1 Trp32 | 22 Å (validada) | 20 | 199 | 0.815 | 4.0 | ✅ válida |
| FUS | modelo 1 | 2 | 87 | 0.609 | 0.0 | ⚠️ débil |
| TDP-43 RRM1 (25 Å) | dominio entero | 5 | 193 | 0.465→0.506 | 0.0 | ❌ no discrimina |
| TDP-43 RNA-face (18 Å) | centrada en RNP1/RNP2 | 5 | 193 | 0.506 | 0.0 | ❌ no discrimina |

**Hallazgos clave:**

1. **TDP-43 no es validable hoy**: ni la caja de 25 Å ni la de 18 Å centrada
   en la cara de unión a RNA separan activos de señuelos (AUC ~0.5 = azar).
   La causa no es la caja, sino los controles: los "activos" de literatura
   (rTRD01, nTRD22, bis-ANS, Congo Red, 5-fluorouridina, isoproterenol) son
   heterogéneos y varios no se unen al RRM1 (bis-ANS va al dominio C-terminal;
   5FUrd/isoproterenol son ligandos cristalizados en SOD1, no en TDP-43).

2. **Artefacto de tamaño en Vina**: señuelos grandes puntúan hasta −28 kcal/mol
   frente a −5.2 del mejor activo. Vina infla la afinidad de moléculas grandes
   en cajas pequeñas → el corte fijo sesga hacia compuestos grandes. Confirma
   que el score de Vina es un *filtro de triaje*, no un ranker.

3. **No hay controles positivos públicos**: ChEMBL (targets CHEMBL2362981 y
   CHEMBL5724679) y PubChem no tienen bioactividades curadas para TDP-43 ni FUS.
   Hasta que existan, el AUC de estas dos dianas no se puede calibrar.

## 2. MM-GBSA de TDP-43 — causa raíz y arreglo

**Diagnóstico.** El runner de Oracle fallaba los 15 candidatos TDP-43 con
timeouts de 30 min, `NaN` y "No template found for residue UNK". La causa era
doble:

1. El receptor PDBQT de AutoDock intercala los átomos por torsión activa
   (GLY142: C,O ... HIS143: N,H ... GLY142: CA), no por residuo. PDBFixer no
   agrupa los residuos → química rota → minimización que no converge (NaN).
2. El receptor de **FUS contiene RNA unido** (guanosina G y uridina U del
   complejo 6G99); PDBFixer no puede parametrizar nucleótidos.

**Arreglo** (`_reordenar_receptor.py` + `rescoring_mmgbsa_robusto.py` v5):

- Reordenar los átomos del PDBQT por (cadena, residuo) en orden canónico
  (N-CA-C-O-CB...), conservando **exactamente el frame del docking**.
- Eliminar residuos no aminoácidos (RNA, aguas).
- Recortar el receptor canónico al entorno del ligando (10 Å) por (cadena,
  residuo) para acelerar la minimización.

**Resultado:** CHEMBL3309775/TDP43 pasó de timeout 30 min → **dG calculado en
<2 min**. Los receptores canónicos quedan en Oracle: `receptor_TDP43_rrm1_canon.pdb`
(786 át), `receptor_FUS_canon.pdb` (391 át, sin RNA), `receptor_SOD1_canon.pdb`.

## 3. Auditoría de poses extremas de FUS

Las poses con dG −30/−32 no tienen colisiones estéricas (0–2 contactos
<0.75·vdW). El patrón sospechoso era `e_ligand` muy negativo (−151, −236),
síntoma de minimización inestable del ligando aislado, no de choques. El
receptor FUS canónico (sin RNA) corrige la geometría; el runner v5 recalcula
los 42 candidatos con el mismo protocolo y produce valores comparables.

## 4. Estado del runner v5 (en curso)

Lanzado en Oracle (`runner_rescoring_v5.py`, `nohup`, PID registrado).
Recalcula los 42 candidatos con receptores canónicos, ~3 min/candidato,
reanudable (checkpoint en `rescoring_mmgbsa_v5.csv`). Los resultados
definitivos se consolidarán al terminar (~1.5–2 h).

## 5. Paper

Actualizada la sección 5.3 (Tabla 1 y su interpretación) con la validación
de 20/08: AUC TDP-43 0.506 (RNA-face 18 Å) y FUS 0.609 (87 señuelos), el
artefacto de tamaño de Vina, y la ausencia de bioactividades curadas públicas
para TDP-43/FUS.

## 6. Corrección definitiva de TDP-43 (17:00–17:40)

### 6.1 El error documental: 6B1N no es TDP-43
El PDB citado como referencia del receptor (6B1N) es **Hsc70/HSPA8**, una
proteína de choque térmico sin relación con TDP-43. La estructura real del
receptor de docking es **4IUF** (TDP-43 RRM1 + ADN co-cristalizado; 773/786
átomos coinciden con el pdbqt). Corregido `proteins/TDP43/metadata.json` y
el paper.

### 6.2 La caja estaba 17 Å fuera del sitio real
El ADN nativo de 4IUF (residuo XUA + cadena B) marca el sitio de unión a RNA
(residuos 109–179, RNP1/RNP2). Las cajas viejas (centros 24.4–28.3, 43.7–44.9,
52.5–53.4) quedaban a ~17 Å del sitio. Por eso el AUC era de azar.

**Caja corregida**: centro (16.3, 41.1, 48.5), 24×24×24 Å.

### 6.3 Re-docking del ligando nativo: RMSD 1.75 Å
Re-acoplamiento del fragmento nativo XUA contra la caja corregida →
**RMSD 1.75 Å < 2 Å** (umbral de éxito). El protocolo de docking de TDP-43
es válido; solo la caja estaba mal colocada.

### 6.4 Validación v4: AUC 0.612
Re-validación con 6.462 señuelos (193 usados) y la caja corregida:
**ROC-AUC = 0.612** (vs 0.506 con la caja vieja). bis-ANS (ligando real de
TDP-43) puntúa −6.58, mejor que el 94.8% de los decoys. Queda el artefacto
de tamaño de Vina (decoys grandes hasta −28).

### 6.5 Impacto en el ranking z001
Re-docking de los top-10 TDP-43 con la caja corregida: CHEMBL6356 (era #1,
−8.6) cae a −7.39; CHEMBL503046 sube a −7.71. Los "récords" viejos eran
artefactos de la caja mal puesta.

### 6.6 Re-docking masivo en curso
Los 2.941 hits TDP-43 ≤ −6.3 de las tandas z001–z010 se re-dockean con la
caja corregida (20 workers CPU, checkpoint reanudable en
`_redock_tdp43_masivo/checkpoint.csv`). Estimado 4–6 h. Al terminar se genera
el ranking definitivo `resultados_tdp43_corregido.csv`.

### 6.7 Bug del receptor canónico (HOH)
El script `_reordenar_receptor.py` incluía `HOH` en el set de residuos
permitidos → los receptores finales retenían aguas que rompían PDBFixer
(timeouts MM-GBSA). Corregido: HOH fuera del set. Receptores regenerados:
TDP43 773 átomos, FUS 391, SOD1 19.783. Runner v5 relanzado limpio (un solo
proceso) a las 16:02 UTC.


### 6.8 Runner v5 MM-GBSA completado (Oracle)
El runner v5 terminó los 42 pares (ligand × target) el 20/08 ~18:04 UTC.
Bug residual del embedding (`EmbedMolecule` sin verificar código de retorno →
`NoneType`) parcheado con reintento robusto; 4 filas reprocesadas. De 42 pares:
**36 válidos, 6 degenerados** (dG=0.0 = complejo equivale a partes separadas,
o TIMEOUT). Ranking consolidado en `ranking_mmgbsa_v5.csv`.

**Top MM-GBSA por target (dG kcal/mol):**
- SOD1: CHEMBL3311449 **−33.42**, CHEMBL8905 −25.81, CHEMBL3310304 −23.75
- TDP43: CHEMBL3311449 **−32.23**, CHEMBL8905 −31.91, CHEMBL9010 −30.62
- FUS: CHEMBL7563 **−22.28**, CHEMBL9532 −21.87, CHEMBL1210578 −19.05

**CHEMBL3311449 (pirrolo[2,3-d]pirimidina) domina SOD1 y TDP43** con dG
≈ −33/−32 kcal/mol, coherente con su rol de líder del rescoring MM-GBSA.

Degenerados a revisar: CHEMBL7360 (SOD1), CHEMBL9440 (SOD1),
CHEMBL601104 (TDP43), CHEMBL3309822 (TDP43), CHEMBL7256 (TDP43).
