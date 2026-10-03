# MASIVE-ALS — Informe de progreso

**Fecha:** 4 de octubre de 2026, 00:30
**Máquina:** PC principal (RTX 4080 12 GB) + Oracle Cloud (VM gratis)
**Frentes activos:** cadena ZINC (Vina-GPU), re-puntaje CNN (GNINA), MM-GBSA de FUS (Oracle)

---

## 1. Librería ChEMBL — TERMINADA

| Concepto | Valor |
|---|---|
| Resultados totales | **431.665** |
| TDP43 | 107.608 |
| SOD1 | 107.608 |
| FUS | 107.608 |
| TDP43_v2 | 108.841 |

Fichero: `gpu_dock/resultados_libreria/resultados_libreria_total.csv`

---

## 2. Cribado ZINC 2M — EN CURSO

- **400 tandas** preparadas en `C:/Users/Fredy/zinc2m/tandas_tgz/lote_zNNN_plano.tar.gz` (~1,93 M ligandos en pdbqt).
- **9 tandas terminadas** (`estado_400/tanda_z001..z009`); **z010 en curso**.
- Fichero acumulado: `gpu_dock/resultados_vinagpu_total.csv` = **382.075 resultados** (todas parejas únicas).

| Diana | Resultados | Mejor energía |
|---|---|---|
| TDP43 | 109.305 | — |
| FUS | 109.362 | −9,10 (CHEMBL1209725) |
| TDP43_v2 | 108.841 | −10,90 (CHEMBL4581197) |
| SOD1 | 54.567 | −9,50 (CHEMBL4559945) |

**No acoplables acumulados:** 5.099 parejas únicas (deduplicado). Dos causas, ambas resueltas dentro de la cadena:

1. **Doble bloque**: el ligando traía dos moléculas concatenadas tras el primer `TORSDOF`. Se recorta al primer bloque.
2. **Tipo de átomo inválido**: `B`, `Si`, `Li`, `Te` no son tipos AutoDock. Se remapean a `C`/`S`.

La cadena `cadena_400.py` aplica ambos parches sola antes del reintento, así que ya no se pierden ligandos.

---

## 3. Re-puntaje CNN (GNINA) — EN CURSO

Puntúa con la red neuronal de GNINA las poses que Vina-GPU ya encontró (barato: <1 s por pose).

| Diana | Top 5.000 | Top 10.000 |
|---|---|---|-|
| TDP43 | ✅ 5.001 | ✅ 10.001 |
| TDP43_v2 | ✅ 5.001 | 5.647 (en curso) |
| SOD1 | ✅ 5.001 | pendiente |
| FUS | ✅ 5.001 | pendiente |

Ficheros: `gpu_dock/gnina_top5000/rescore_<diana>.csv`

---

## 4. Docking real GNINA del consenso — COMPLETO

| Diana | Ligandos | Dockings | Fichero |
|---|---|---|---|
| TDP43 | 60 | 120 | `resultados/consenso_gnina_tdp43.csv` |
| SOD1 | 31 | 62 | `resultados/consenso_gnina_sod1.csv` |

Cada ligando se acopló contra **dos conformaciones** del receptor (p. ej. TDP43 y TDP43_v2).
Además: `resultados/rescate8_sod1.csv` (16 filas, ligandos de cristal SOD1).

---

## 5. MM-GBSA de FUS (Oracle) — EN CURSO, 105/119

- Proceso vivo: `runner_mmgbsa_fus.py` en `~/mmgbsa` (VM Oracle).
- **105 de 119** candidatos procesados: **91 con energía**, 13 timeouts (>1800 s), ~14 pendientes.
- Mejores dG hasta ahora: CHEMBL3311273 (−25,62), CHEMBL11378 (−25,05), CHEMBL503752 (−24,99), CHEMBL4584807 (−24,12), CHEMBL508135 (−23,25).
- Parcial traído a `resultados/mmgbsa_fus_oracle_parcial.csv`.
- **Vigilante activo** (`vigilar_mmgbsa_fus.py`): trae el parcial cada 5 min y, al terminar, recalcula el ranking en Oracle y regenera el consenso de tres motores solo.

---

## 6. Cruce de tres motores (MM-GBSA + Vina + GNINA)

Fichero: `resultados/consenso_tres_motores.csv` (202 candidatos)

| Concepto | Valor |
|---|---|
| Candidatos totales | 202 |
| Con los **tres** motores | 35 |
| Buenos en los tres | 34 |

### Top 5 del consenso

| # | Diana | Ligando | MM-GBSA | Vina | CNN affinity |
|---|---|---|---|---|---|
| 1 | TDP43 | **CHEMBL4574530** | −36,01 | −10,8 | 7,769 |
| 2 | TDP43 | CHEMBL4559959 | −38,07 | −9,8 | 7,803 |
| 3 | TDP43 | CHEMBL4584797 | −33,29 | −10,4 | 7,819 |
| 4 | TDP43 | CHEMBL530000 | −34,7 | −9,6 | 7,129 |
| 5 | TDP43 | CHEMBL4581197 | −31,24 | −10,9 | 7,027 |

> El antiguo número 1 por MM-GBSA, CHEMBL1795884 (−48,29), cae al puesto 60 por un Vina flojo.

---

## 7. Infraestructura en marcha

| Proceso | Dónde | Estado |
|---|---|---|
| `cadena_400.py` (Vina-GPU) | PC, RTX 4080 | ✅ corriendo (z010) — 100% GPU |
| `gnina_repuntar.py` (CNN) | WSL2, CPU | ✅ corriendo (TDP43_v2) | 
| `vigilar_mmgbsa_fus.py` | PC | ✅ vigilando Oracle |
| `vigilar_cnn_10k.py` | PC | ✅ vigilará y regenerará el consenso |
| `runner_mmgbsa_fus.py` | Oracle | ✅ corriendo (105/119) |

---

## 8. Qué queda por delante

- **391 tandas ZINC** restantes (z010→z400) a ~40 min/tanda → semanas de cómputo continuo.
- **CNN a 10.000** en las cuatro dianas (TDP43 ya hecho).
- **MM-GBSA de FUS**: ~14 candidatos pendientes.
- Cuando CNN y FUS cierren, el consenso se regenera solo.

---

*Generado automáticamente por Luz (Buffy) — 4 de octubre de 2026.*