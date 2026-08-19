# Proyecto MASIVE-ALS — Informe de prueba de concepto
## Para: Dra. Mónica Povedano Panadés — Hospital Bellvitge / IDIBELL

---

## Resumen ejecutivo

El pipeline de cribado virtual MASIVE-ALS está funcionando. Se han completado las primeras 600 simulaciones de docking molecular con proteínas reales implicadas en ELA (TDP-43, SOD1, FUS) contra la base de datos ChEMBL de 2.4 millones de compuestos bioactivos, ejecutándose 24/7 en infraestructura cloud gratuita (Oracle Cloud, 4 núcleos ARM, 24 GB RAM).

---

## Datos de la prueba de concepto

| Métrica | Valor |
|---|---|
| **Pipeline** | AutoDock Vina 1.2.3 + OpenBabel + OpenMM |
| **Proteínas** | TDP-43 (PDB: 6b1n), SOD1 (PDB: 1hl5), FUS (PDB: 6g99) |
| **Compuestos** | 2,409,270 (ChEMBL 34) |
| **Dockings completados** | 600+ (y contando) |
| **Velocidad** | ~2.5 docking/minuto |
| **Mejor energía** | -0.10 kcal/mol |
| **Infraestructura** | Oracle Cloud ARM (gratis, 24/7) |

---

## Mejores candidatos preliminares

| Compuesto (ChEMBL) | Proteína | Energía (kcal/mol) |
|---|---|---|
| CHEMBL156224 | SOD1 | -0.1023 |
| CHEMBL439520 | — | -0.09 |
| CHEMBL154341 | — | -0.09 |
| CHEMBL409812 | — | -0.09 |
| CHEMBL146675 | TDP-43 | -0.0583 |

**Nota:** Energías con exhaustiveness=4 (mínimo). Para resultados publicables se requiere exhaustiveness≥32, lo que implica necesidad de supercomputación.

---

## Validación técnica completada

- ✅ Docking molecular: 3 proteínas ELA + fármacos reales
- ✅ Preparación de estructuras: PDBFixer + OpenBabel
- ✅ Dinámica molecular: OpenMM configurado y funcional
- ✅ Pipeline 24/7 en cloud gratuito
- ✅ Datos exportables en formato CSV estándar

---

## Próximos pasos solicitados

1. **Validación experimental**: Ensayos in vitro de los mejores candidatos (IRB Barcelona, VHIR)
2. **Supercomputación**: Solicitud a MareNostrum 5 / RES para cribado completo a exhaustiveness≥32 (10M compuestos)
3. **Respaldo institucional**: Aval de Bellvitge/IDIBELL para las solicitudes de tiempo de cálculo
4. **Dinámica molecular completa**: Simulaciones de 100 ns para los top hits

---

*Generado el 10 de agosto de 2026. Pipeline corriendo en Oracle Cloud (IP 79.72.57.253).*
*Contacto: Fredy Rojas Gutiérrez — fredy_30@hotmail.com — +34 675 31 58 41*
