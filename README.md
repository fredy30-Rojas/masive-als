# MASIVE-ALS: Molecular Docking at Exascale for ALS Proteinopathy

**Investigador Principal:** Fredy Rojas Gutiérrez  
**Supercomputador:** MareNostrum 5 (BSC-CNS, Barcelona)  
**Programa:** Red Española de Supercomputación (RES)  
**Horas solicitadas:** 200,000 GPU + 50,000 CPU  
**Duración:** Septiembre 2026 - Febrero 2027  

---

## 🔬 Objetivo

Identificar 3-5 fármacos candidatos contra la ELA mediante cribado virtual masivo de **10 millones de compuestos** contra las 3 proteínas clave de la enfermedad: **TDP-43, SOD1 y FUS**.

---

## 📁 Estructura del proyecto

```
masive-als/
├── README.md                    ← Este archivo
├── setup_env.sh                 ← Preparación del entorno en MareNostrum 5
├── scripts/
│   ├── slurm_docking.sh         ← Cribado masivo (200 GPUs, 50 nodos)
│   └── slurm_gromacs.sh         ← Validación MD (40 GPUs, 10 nodos)
├── prep/
│   ├── download_compounds.sh    ← Descarga ZINC20, DrugBank, Enamine REAL
│   └── prepare_proteins.sh      ← Preparación de conformaciones proteicas
├── analysis/
│   ├── merge_results.py         ← Análisis de resultados de docking
│   ├── prepare_md.py            ← Preparación de sistemas para MD
│   └── analyze_md.py            ← Análisis de dinámica molecular
├── compounds/                   ← Librerías de compuestos (PENDIENTE en Marenostrum)
├── proteins/                    ← Conformaciones proteicas (PENDIENTE en Marenostrum)
└── results/                     ← Resultados de docking y MD
```

---

## 🚀 Flujo de trabajo

### Fase 0 - Preparación (1-2 semanas)
```bash
# En MareNostrum 5, una vez concedido el acceso:
ssh usuario@marenostrum5.bsc.es
cd /gpfs/projects/masive-als
bash setup_env.sh
bash prep/download_compounds.sh
bash prep/prepare_proteins.sh
```

### Fase 1 - Cribado masivo (oct-nov 2026)
```bash
sbatch scripts/slurm_docking.sh
# 50 tareas en paralelo, 200 GPUs simultáneas
# ~5 horas por tarea = 1 semana total
```

### Fase 2 - Segunda ronda (dic 2026)
```bash
# Seleccionar top 1,000 hits del Fase 1
python3 analysis/merge_results.py --input results/ --output results/top_hits.csv --top 1000
sbatch scripts/slurm_docking.sh  # Segunda tanda
```

### Fase 3 - Validación MD (ene 2027)
```bash
python3 analysis/prepare_md.py --hits results/top_hits.csv --batch 1
sbatch scripts/slurm_gromacs.sh
# 1 microsegundo de MD por candidato
python3 analysis/analyze_md.py --dir md_results/ --output results/final_candidates.csv
```

### Fase 4 - Publicación (feb 2027)
- Resultados en Zenodo (CC-BY 4.0)
- Código en GitHub: github.com/fredy30-Rojas/masive-als
- Citación a RES y BSC

---

## 🧬 Proteínas diana

| Proteína | Rol en ELA | Estrategia |
|---|---|---|
| **TDP-43** | 97% pacientes: agregación citoplasmática | Desagregar acúmulos |
| **SOD1** | 20% ELA familiar: radicales libres | Estabilizar dímero nativo |
| **FUS** | Transición líquido-sólido patológica | Bloquear fase aberrante |

## 🛠 Herramientas

| Software | Uso | Licencia |
|---|---|---|
| AlphaFold-Multimer | Predicción de estructura 3D | Open source |
| AutoDock-GPU | Docking molecular masivo | Open source |
| GROMACS | Dinámica molecular | Open source |
| OpenBabel | Conversión de formatos químicos | Open source |

---

## 📊 Estimaciones de rendimiento

| Recurso | PC normal (4 GPUs) | MareNostrum 5 (200 GPUs) |
|---|---|---|
| 10M compuestos | ~8 años | ~6 meses |
| Docks/segundo | ~2,000 | ~400,000 |
| MD 1 µs | ~30 días | ~2 días |

---

*Proyecto MASIVE-ALS — "No busco publicar un artículo. Busco vivir."*
