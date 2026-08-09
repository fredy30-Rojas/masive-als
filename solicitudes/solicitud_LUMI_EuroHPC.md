# Solicitud de Acceso — LUMI (EuroHPC JU)
## Development Access Mode — MASIVE-ALS

**Fecha:** 9 de agosto de 2026
**Investigador Principal:** Fredy Rojas Gutiérrez
**Contacto:** fredy_30@hotmail.com | +34 675 31 58 41
**Ubicación:** Rubí, Barcelona, España

---

## 1. Título del Proyecto

**MASIVE-ALS: Cribado Molecular Masivo contra la Esclerosis Lateral Amiotrófica**

---

## 2. Resumen Científico (máx. 2000 caracteres)

La Esclerosis Lateral Amiotrófica (ELA) afecta a 350,000 personas en el mundo y 4,000 en España. No tiene cura. El proyecto MASIVE-ALS propone un cribado virtual masivo de **10 millones de compuestos** contra las 3 proteínas clave de la enfermedad: **TDP-43** (agregación citoplasmática en el 97% de pacientes), **SOD1** (radicales libres en el 20% de ELA familiar) y **FUS** (transición líquido-sólido patológica).

Usando AlphaFold-Multimer (Nobel de Química 2024) generaremos 100,000 conformaciones proteicas. Con AutoDock-GPU ejecutaremos 1 billón de simulaciones de docking. Los mejores 1,000 candidatos se validarán con dinámica molecular de 1 microsegundo en GROMACS.

Recursos solicitados en LUMI-G (GPU):
- **50,000 GPU-node-hours** en la partición LUMI-G (AMD MI250X)
- **10,000 CPU-core-hours** para preparación
- **20 TB** almacenamiento en Lustre
- Duración: 6 meses

**Nota:** Esta es una solicitud de Development Access, ideal para portar y escalar nuestro pipeline de docking antes de solicitar una asignación mayor en la convocatoria Extreme Scale.

---

## 3. Justificación Técnica para LUMI

El pipeline MASIVE-ALS está optimizado para GPU:

| Componente | GPU requerida | Software |
|---|---|---|
| Docking masivo | AMD MI250X | AutoDock-GPU 1.6 (OpenCL) |
| Dinámica molecular | AMD MI250X | GROMACS 2024.3 (HIP/ROCm) |
| Predicción estructural | AMD MI250X | AlphaFold-Multimer |

AutoDock-GPU escala linealmente con el número de GPUs. Cada GPU MI250X procesa ~2,000 docks/segundo. Con 128 GPUs simultáneas alcanzamos 256,000 docks/segundo, completando el cribado en semanas en lugar de años.

El código está disponible en GitHub: `github.com/fredy30-Rojas/masive-als`

---

## 4. Plan de Trabajo en LUMI

| Mes | Actividad | GPUs |
|---|---|---|
| Sept 2026 | Portar y optimizar AutoDock-GPU para ROCm/MI250X | 8 |
| Oct-Nov 2026 | Cribado Fase 1: 5M compuestos | 128 |
| Dic 2026 | Cribado Fase 2: 5M compuestos | 128 |
| Ene 2027 | Dinámica molecular: top 1,000 hits | 64 |
| Feb 2027 | Análisis y publicación | 4 |

---

## 5. Resultados Esperados

- 3-5 fármacos candidatos contra la ELA listos para validación experimental
- Datos publicados en Zenodo con licencia CC-BY 4.0
- Código en GitHub, ciencia reproducible
- Citación a EuroHPC JU y LUMI en todas las publicaciones

---

## 6. Información Adicional

El investigador principal, Fredy Rojas, es paciente de ELA. Escribe con control ocular Tobii 4C. Este proyecto se construyó sin mover las manos.

"Quien firma no busca publicar un artículo. Busca vivir."

---

**Solicitud enviada a través de:** EuroHPC JU Development Access
**Sistema:** LUMI (CSC, Finlandia)
**Tipo:** Development → preparación para Extreme Scale
