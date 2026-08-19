# Proyecto MASIVE-ALS — Informe de progreso (10 de agosto de 2026)
## Para: Equipo investigador — validación y colaboración

---

## Resumen ejecutivo

El pipeline de cribado virtual MASIVE-ALS está operativo 24/7 y ha completado **800 simulaciones de docking molecular** con proteínas reales implicadas en ELA (TDP-43, SOD1, FUS) contra la base de datos ChEMBL de 2.4 millones de compuestos bioactivos.

**Este documento resume el progreso real y solicita su colaboración para los siguientes pasos.**

---

## Estado actual del pipeline

| Métrica | Valor |
|---|---|
| **Pipeline** | AutoDock Vina 1.2.3 + OpenBabel + OpenMM |
| **Proteínas** | TDP-43 (PDB: 6b1n), SOD1 (PDB: 1hl5), FUS (PDB: 6g99) |
| **Compuestos objetivo** | 2,409,270 (ChEMBL 34) |
| **Dockings completados** | 800+ (y corriendo 24/7) |
| **Velocidad** | ~2.4-2.5 docking/minuto |
| **Infraestructura** | Oracle Cloud ARM (gratuita, 24/7) |
| **Estado** | Activo y generando resultados continuamente |

---

## Mejores candidatos preliminares (datos reales)

| Compuesto (ChEMBL) | Proteína | Energía (kcal/mol) |
|---|---|---|
| CHEMBL156224 | SOD1 | -0.10 |
| CHEMBL503549 | SOD1 | -0.09 |
| CHEMBL439520 | TDP-43 | -0.09 |
| CHEMBL154341 | TDP-43 | -0.09 |
| CHEMBL409812 | FUS | -0.09 |

**Nota técnica:** Energías calculadas con exhaustiveness=4 (parámetro mínimo de validación). Para resultados publicables y selección definitiva de candidatos se requiere exhaustiveness≥32, que implica supercomputación (MareNostrum 5 / RES).

---

## Validación técnica completada

- ✅ Docking molecular funcional con 3 proteínas ELA y fármacos reales
- ✅ Preparación de estructuras (PDBFixer + OpenBabel)
- ✅ Dinámica molecular OpenMM configurada y operativa
- ✅ Pipeline autónomo 24/7 en cloud gratuito
- ✅ Datos exportables en CSV estándar (resultados por compuesto, proteína y energía)

---

## Qué solicitamos a los grupos de investigación

1. **Validación experimental**: Ensayos in vitro de los mejores candidatos (IRB Barcelona, VHIR)
2. **Respaldo científico**: Aval institucional para solicitar tiempo de cálculo en MareNostrum 5 / RES
3. **Supercomputación**: Cribado completo a exhaustiveness≥32 sobre 10M compuestos
4. **Colaboración clínica**: Orientación sobre relevancia terapéutica de los candidatos

---

*Pipeline corriendo en Oracle Cloud (IP 79.72.57.253). Todos los resultados estarán disponibles en acceso abierto.*
*Contacto: Fredy Rojas Gutiérrez — fredy_30@hotmail.com — +34 675 31 58 41 — Rubí, Barcelona.*
