# Solicitud de Acceso — Frontier (OLCF)
## Director's Discretionary Program — MASIVE-ALS

**Fecha:** 9 de agosto de 2026
**Investigador Principal:** Fredy Rojas Gutiérrez
**Contacto:** fredy_30@hotmail.com | +34 675 31 58 41
**Ubicación:** Rubí, Barcelona, España
**Afiliación sugerida:** Hospital Universitari de Bellvitge (Dra. Mónica Povedano)

---

## 1. Project Title

**MASIVE-ALS: Massive Molecular Screening for Amyotrophic Lateral Sclerosis**

---

## 2. Scientific Summary

ALS affects 350,000 people worldwide with no cure. MASIVE-ALS proposes the largest virtual screening campaign ever conducted against ALS: **10 million compounds** screened against 3 key proteins (**TDP-43**, **SOD1**, **FUS**) using AlphaFold-Multimer, AutoDock-GPU, and GROMACS.

Total computational scale: 1 trillion docking simulations. On a standard workstation: 8 years. On Frontier: **weeks**.

---

## 3. Resources Requested

| Resource | Amount | Partition |
|---|---|---|
| GPU-node-hours | 80,000 | AMD MI250X |
| CPU-core-hours | 15,000 | CPU |
| Storage | 30 TB | Orion filesystem |
| Duration | 6 months | Sep 2026 - Feb 2027 |

---

## 4. Technical Justification

AutoDock-GPU runs natively on AMD GPUs via OpenCL. MI250X provides 2 GCDs per GPU, effectively doubling throughput. GROMACS 2024 supports AMD HIP/ROCm backend.

| Benchmark | Per MI250X GPU | 200 GPUs |
|---|---|---|
| Docks/second | ~2,500 | 500,000 |
| MD performance | ~120 ns/day | ~24 µs/day |

Frontier's architecture (AMD EPYC + MI250X) is ideal for this workload.

---

## 5. Timeline

| Phase | Activity | GPUs |
|---|---|---|
| Sep 2026 | Setup, ZINC20 download, AlphaFold | 16 |
| Oct-Nov 2026 | Phase 1 screening (5M compounds) | 200 |
| Dec 2026 | Phase 2 screening (5M compounds) | 200 |
| Jan 2027 | MD validation (top 1,000 hits) | 80 |
| Feb 2027 | Analysis and open-access publication | 8 |

---

## 6. Deliverables

- 3-5 ALS drug candidates ready for experimental validation
- Open-access data (CC-BY 4.0) on Zenodo
- Open-source code on GitHub
- Citation to OLCF and DOE in all publications

---

## 7. Collaborating Institutions

- Hospital de Bellvitge, Barcelona (Dra. Mónica Povedano — ALS clinical trials)
- IRB Barcelona (in vitro validation)
- VHIR Vall d'Hebron (in vivo models)

---

## 8. Personal Statement

The PI is an ALS patient who writes with eye-tracking technology (Tobii 4C). This project was built without moving his hands. He lives 20 km from his first-choice supercomputer (MareNostrum 5), but is pursuing all available paths.

*"I am not trying to publish a paper. I am trying to live."*

— Fredy Rojas Gutiérrez

---

**Submit to:** OLCF Director's Discretionary Program
**System:** Frontier (Oak Ridge National Laboratory)
**Partition:** AMD MI250X GPU
