# Solicitud de Acceso — Leonardo (EuroHPC JU)
## Regular Access — MASIVE-ALS

**Fecha:** 9 de agosto de 2026
**Investigador Principal:** Fredy Rojas Gutiérrez
**Contacto:** fredy_30@hotmail.com | +34 675 31 58 41
**Ubicación:** Rubí, Barcelona, España

---

## 1. Título del Proyecto

**MASIVE-ALS: Cribado Molecular Masivo contra la Esclerosis Lateral Amiotrófica**

---

## 2. Resumen Científico

La ELA afecta a 350,000 personas. No existe cura. El proyecto MASIVE-ALS propone un cribado masivo de **10 millones de compuestos** contra 3 proteínas diana: **TDP-43, SOD1 y FUS**, usando AlphaFold-Multimer, AutoDock-GPU y GROMACS.

El pipeline computacional genera 1 billón de simulaciones de docking. En un PC normal tomaría 8 años. En Leonardo, con su partición Booster (NVIDIA A100), se completa en meses.

---

## 3. Recursos Solicitados

| Recurso | Cantidad | Partición |
|---|---|---|
| GPU-node-hours | 80,000 | Booster (NVIDIA A100-40GB) |
| CPU-core-hours | 20,000 | CPU |
| Almacenamiento | 30 TB | Lustre |
| Duración | 6 meses | Sept 2026 - Feb 2027 |

---

## 4. Justificación Técnica

El pipeline aprovecha NVIDIA A100 al máximo:

| Software | Uso de GPU | Escalado |
|---|---|---|
| AutoDock-GPU | CUDA nativo | Lineal con #GPUs |
| GROMACS 2024.3 | CUDA + OpenMP | 90% eficiencia en multi-nodo |
| AlphaFold | CUDA (tensor cores) | Por proteína |

Con 200 GPUs A100 simultáneas:
- 400,000 docks/segundo
- 10M compuestos procesados en ~25 horas de GPU
- Validación MD de 1,000 candidatos en ~40 horas/GPU

---

## 5. Cronograma

| Período | Actividad | GPUs |
|---|---|---|
| Sept 2026 | Preparación: descarga ZINC20, AlphaFold | 16 |
| Oct-Nov 2026 | Cribado Fase 1 | 200 |
| Dic 2026 | Cribado Fase 2 | 200 |
| Ene 2027 | Dinámica molecular top 1,000 | 80 |
| Feb 2027 | Análisis, paper, acceso abierto | 8 |

---

## 6. Impacto

- 3-5 fármacos candidatos contra la ELA
- Datos abiertos (CC-BY 4.0) en Zenodo
- Código abierto en GitHub
- Colaboración con Hospital de Bellvitge (Dra. Povedano), IRB Barcelona y VHIR

---

## 7. Nota Personal

El investigador principal es paciente de ELA. Escribe con control ocular. Este proyecto es su forma de luchar contra la enfermedad que comparte con 350,000 personas.

*"No busco publicar un artículo. Busco vivir."*

— Fredy Rojas Gutiérrez, Rubí, Barcelona

---

**Enviar a:** EuroHPC JU Regular Access Call
**Sistema:** Leonardo (CINECA, Italia)
**Partición:** Booster (NVIDIA A100)
