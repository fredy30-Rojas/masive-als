# REVISIÓN CIENTÍFICA MASIVE-ALS — 18 agosto 2026

Revisión metodológica solicitada por Fredy. Complementa la auditoría de datos del 17-18/08 (bolsillo, embudo, documentación), ya corregida.

## RESUMEN EJECUTIVO

El diseño general es sólido: el repurposing de fármacos aprobados es la vía más realista hacia 3-5 candidatos utilizables, y la corrección del bolsillo fue el paso correcto. Pero el ranking actual tiene 3 fallos metodológicos que pueden invalidar la selección:

1. Sin controles positivos (ligandos conocidos de cada proteína), el corte por percentil no está calibrado.
2. El corte top-5% se aplicó global sobre las 3 proteínas, y los scores de Vina NO son comparables entre receptores distintos. Resultado medible: 138 candidatos SOD1, 17 TDP-43, 16 FUS (81% SOD1).
3. No hay filtro de barrera hematoencefálica (BBB/CNS), crítico en una enfermedad neurodegenerativa.

Ninguno es irreparable. Remedios concretos en cada sección.

---

## 1. CRÍTICO — Ranking no calibrado (sin controles positivos)

- Vina tiene error de afinidad de ~1-2 kcal/mol y RMSD de ~2 Å. El corte top-5% (umbral ≈ -6.7/-6.8 kcal/mol) es arbitrario sin calibrar contra ligandos conocidos.
- **Remediación estándar**: acoplar un set de ligandos CONOCIDOS de cada proteína (TDP-43: ligadores de RRM reportados; SOD1: compuestos anti-agregación reportados; FUS: ligadores de dominio de baja complejidad) y calcular enrichment factor / ROC-AUC del ranking. Si los conocidos no aparecen en el top 5%, el protocolo no discrimina entre binders y no-binders.
- **Redocking de validación**: reproducir una pose cocristalizada conocida (RMSD < 2 Å) en cada caja antes de confiar en cualquier ranking. No consta que las cajas se hayan validado con ligando cocristalizado.

## 2. CRÍTICO — El corte global mezcla scores incomparables entre proteínas

- Los scores de Vina solo son comparables DENTRO de un mismo receptor/bolsillo: dependen del tamaño del receptor, del número de átomos y del tamaño del bolsillo. -8.6 en TDP-43 no es comparable con -8.6 en SOD1.
- El embudo se corrió sin `--proteina` sobre resultados_z001.csv: el umbral global (-6.70) admitió desproporcionadamente filas de SOD1, cuya distribución de energías es más favorable. Resultado en candidatos_filtrados.csv: **138 SOD1 / 17 TDP43 / 16 FUS (81% SOD1)**.
- Umbrales top-5% POR PROTEÍNA (los correctos): FUS -6.40, TDP43 -6.40, SOD1 -7.10.
- **Remediación**: correr el embudo una vez por proteína (`--proteina TDP43`, `--proteina SOD1`, `--proteina FUS`) con su propio percentil, y unir los resultados. Es un cambio de un comando por proteína.

## 3. CRÍTICO — Sin filtro de barrera hematoencefálica (BBB/CNS)

- La ELA es una enfermedad neurodegenerativa: el candidato DEBE cruzar la barrera hematoencefálica para llegar a las motoneuronas. Lipinski/Veber no cubren penetración CNS.
- **Remediación**: añadir criterios CNS al embudo: TPSA ≤ 90 Å², MW ≤ 450 (orientativo), o el score CNS-MPO (Wager 2010, ≥ 4). El filtro actual puede aprobar moléculas que jamás alcanzarán las motoneuronas.

## 4. BIOLÓGICO — Mecanismo de los blancos (declarar en el paper)

- **TDP-43**: el bolsillo está en RRM1 (dominio de unión a RNA), no en el dominio de agregación (LCD C-terminal, ~aa 274-414). Ocupar RRM1 puede modular la unión a RNA — mecanismo de toxicidad discutido en la literatura reciente — pero NO previene directamente la agregación citoplasmática. Hay que declarar el mecanismo postulado.
- **FUS**: el receptor usa el modelo 1 de 20 del NMR — elección arbitraria; usar el ensemble (centroide de cluster) o declarar por qué modelo 1. Además, la caja de 25 Å cubre ~91% de los átomos del receptor (514 átomos): es docking casi ciego, no dirigido. Para un estudio dirigido hay que definir un sitio específico.
- **SOD1**: la especie patogénica es el monómero apo (sin Cu/Zn) misfolded; la enfermedad clásica implica disociación del dímero. Acoplar al tetrámero nativo holometalado con caja en la región de agregación C-terminal es una aproximación pragmática válida para anti-agregación, pero debe declararse explícitamente que no modela el monómero misfolded ni la interfaz del dímero.

## 5. PIPELINE — Huecos operativos

- **El CSV de candidatos no referencia los poses**: candidatos_filtrados.csv exporta ligand/target/affinity/smiles pero NO la ruta al PDBQT del complejo acoplado. Para MM-GBSA se necesitan los complejos (pose + receptor), no los SMILES. Cada candidato debe mapear a su pose de docking.
- **Cobertura del top**: 292/930 filas (31%) del top-5% quedaron sin SMILES y fueron descartadas. Con 56.772 SMILES cargados, hay identificadores (¿ZINC?) sin resolver en las librerías. El top-5% "oficial" excluye un tercio de su propia selección.
- **enriquecer_smiles.py**: la búsqueda por nombre toma el PRIMER hit de ChEMBL; para nombres genéricos el primer hit puede ser una sal o forma distinta de la molécula correcta. Validar manualmente los ~98 nombres añadidos vía búsqueda por nombre.
- **MM-GBSA en GROMACS**: requiere gmx_MMPBSA y decisiones no documentadas: campo de fuerza (ff14SB/GAFF2), modelo de agua (TIP3P), longitud de simulación y número de frames. Además, MM-GBSA correlaciona modestamente con la afinidad experimental: es un ranker adicional, no evidencia de unión.

## 6. REPRODUCIBILIDAD

- Documentar versión exacta de Vina-GPU (2.1), exhaustiveness usado, seed aleatorio, y plataforma (Windows local vs Kaggle P100 vs Oracle ARM): los scores pueden variar entre plataformas.
- El total de pares cribados quedó inconsistente en el RESUMEN (línea de COMPUESTOS decía 27.721 + z001 en curso; el estado actual dice 32.699 únicos tras integración de Buffy). Corregido en esta revisión.

## 7. ESTRATÉGICO — Recomendación honesta

- El entregable realista de este proyecto, sin ensayos experimentales, es una lista corta priorizada de fármacos aprobados (o en fase clínica) con justificación mecanística y perfil de seguridad conocido: eso es lo que un neurólogo puede evaluar y actuar.
- El top actual de candidatos son bioactivos CHEMBL sin perfil clínico. Considerar priorizar por: fármaco aprobado o en fase clínica + permeable a CNS + afinidad. Es decir, añadir dos columnas al CSV de candidatos: `aprobado` (sí/no) y `permeabilidad_CNS` (estimada por TPSA/MW).
