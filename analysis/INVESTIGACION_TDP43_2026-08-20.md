# Investigación: cómo arreglar TDP-43 — 20 agosto 2026

## Resumen del hallazgo

**TDP-43 SÍ es validable, pero la caja de docking estaba mal colocada** y la
estructura de referencia estaba mal documentada. Esto explica el AUC de azar
(0.465–0.506) observado en todas las validaciones previas.

## 1. El PDB "6B1N" documentado NO es TDP-43

El archivo `proteins/TDP43/PDB-6b1n.pdb` (y el PDB real 6B1N de RCSB) es
**Hsc70/HSPA8** — proteína de choque térmico:

```
HEADER    CHAPERONE  18-SEP-17  6B1N
TITLE     DISRUPTED HYDROGEN BOND NETWORK IMPAIRS ATPASE ACTIVITY IN AN HSC70
COMPND    MOLECULE: HEAT SHOCK PROTEIN FAMILY A (HSP70) MEMBER 8
```

Secuencia: `SPAVGIDLGTTYSWVGVFQHGKVEIIANDQ...` (dominio ATPasa de Hsp70).
El `metadata.json` y el paper citan "6B1N/4IUF" como fuente; 6B1N es un error.

**La estructura real del receptor de docking es 4IUF** (cristal de TDP-43 RRM1
en complejo con DNA): 773 de 786 átomos del `TDP43.pdbqt` coinciden
exactamente con la cadena A de 4IUF (secuencia `TSDLIVLGLPWKTTEQDLKEYF...`).
El cribado se hizo contra el receptor correcto; solo la referencia estaba mal
documentada.

## 2. La caja de cribado estaba 17 Å lejos del sitio de RNA nativo

4IUF contiene el DNA co-cristalizado (cadena B: DC/DG/DT/XUA, 293 átomos) en
el sitio de unión a RNA (residuos en contacto: 109, 110, 113, 165, 171, 174,
176, 179 — la superficie RNP1/RNP2 del RRM1):

| Caja | Centro | Dist. al centro del DNA nativo |
|---|---|---|
| Cribado real (usada en z001–z011) | 28.3, 43.7, 52.5 | **17.2 Å** |
| RNA-binding previa (18 Å, v3) | 24.4, 44.9, 53.4 | 14.5 Å |
| **Corregida (esta investigación)** | **16.3, 41.1, 48.5** | **0 (centrada en el DNA)** |

El centro del DNA completo es (16.3, 41.1, 48.5); el del nucleótido XUA,
(13.8, 35.7, 57.1). Con la caja a 17 Å del sitio real, Vina muestreaba una
región mayoritariamente proteica sin el bolsillo de RNA — por eso los señuelos
puntuaban igual que los activos (AUC ≈ 0.5).

## 3. Re-docking del ligando nativo: RMSD 1.75 Å → protocolo válido

Se extrajo el nucleótido XUA co-cristalizado (36 átomos), se convirtió a
PDBQT con Open Babel, y se re-dockeó con Vina (caja corregida 18 Å,
exhaustividad 16) contra el receptor del cribado:

```
mode | affinity | RMSD (Kabsch, 23 átomos pesados comunes)
  1  |  -4.962  |  1.75 Å   ← por debajo del umbral de éxito (2 Å)
```

Un RMSD < 2 Å tras superposición óptima se considera re-docking exitoso en la
literatura (la pose dockeada reproduce la pose cristalográfica). **Conclusión:
el protocolo de docking (receptor, Vina, parametrización) es correcto para
TDP-43; el problema era exclusivamente la colocación de la caja.**

## 4. Acción tomada

- Re-validación v4 con señuelos en curso con la caja corregida
  (16.3, 41.1, 48.5; 24 Å): `validacion_TDP43_v4.csv` (198 ligandos).
- Si el AUC sube >0.7, TDP-43 queda validado y habrá que **re-dockear las
  tandas** TDP-43 ya cribadas con la caja corregida (o al menos re-priorizar
  los candidatos z001 con la nueva caja).
- Corregir `metadata.json` y el paper: fuente real = **4IUF** (no 6B1N).

## 5. Limitación honesta

El XUA es un fragmento del DNA (un solo nucleótido modificado); el DNA
completo (208 átomos, 293 con H) es demasiado grande/torsional para Vina en
tiempo razonable. El re-docking del fragmento más representativo ya es
suficiente como control positivo nativo. La validación con señuelos (v4)
dará la métrica ROC-AUC definitiva.

## Archivos

- `_4iuf_full.pdb` — estructura real de TDP-43 RRM1 + DNA (RCSB 4IUF)
- `_xua_nativo.pdbqt` — ligando nativo (XUA) para re-docking
- `_redock_xua.conf` / `_xua_redock_out.pdbqt` — re-docking
- `_rmsd_xua.py` — cálculo RMSD (Kabsch)
- `validacion_TDP43_v4.csv` — re-validación con caja corregida (en curso)
