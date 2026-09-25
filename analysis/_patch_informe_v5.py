# -*- coding: utf-8 -*-
"""Añade la seccion 6.8 (runner v5 completado + ranking MM-GBSA) al informe."""
import io

P = r"C:/Users/Fredy/masive-als/analysis/ARREGLO_PIPELINE_2026-08-20.md"
adic = """

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
"""

with open(P, "a", encoding="utf-8") as f:
    f.write(adic)
print("ok")
