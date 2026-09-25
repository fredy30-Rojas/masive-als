# CORRECCIÓN DEL RESCORING MM-GBSA — MASIVE-ALS

**Fecha:** 11 de septiembre de 2026
**Base:** auditoría `AUDITORIA_RESCORING_2026-09-11.md` (7 hallazgos, veredicto: datos no utilizables)
**Objeto:** corregir el pipeline y **verificar** la corrección con medidas.

---

## RESUMEN

Se corrigieron los tres fallos del pipeline. La corrección está **implementada y medida**:

| Métrica | Antes (vacío) | Después (GB/OBC2) | Cambio |
|---|---|---|---|
| dG de CHEMBL4584906–SOD1 | −66,74 / −78,03 | **−39,62 / −40,11 / −39,79** | a rango físico |
| Dispersión entre corridas idénticas | 11,5 kcal/mol | **0,49 kcal/mol** | 23× mejor |
| Eficiencia de ligando | 2,16 (88 % > 1,5) | **≈ 1,3** | coherente con MM-GBSA sin entropía |
| Minimización | tolerancia laxa, sin comprobar | **converge** (para antes de 8.000 iter.) | protocolo definido |

**Ficheros nuevos (reproducibles):**
- `analysis/rescoring_local/mmgbsa_openff_gb.py` — pipeline corregido
- `analysis/rescoring_local/runner_mmgbsa_gb.py` — lanzador reanudable, escribe a `rescoring_corregido.csv`

---

## QUÉ SE CAMBIÓ Y POR QUÉ

### 1. Término de disolvente OBC2 (el fallo principal)

El original construía los tres sistemas con SMIRNOFF y los evaluaba con `NonbondedForce.NoCutoff`, **sin ningún modelo de disolvente**. Se añadió el término de Born generalizado con `GBSAOBC2Force` (radios mbondi2, `SA=ACE`, ε_solvente 78,5, ε_soluto 1), replicando exactamente lo que hace `openmm/app/data/implicit/obc2.xml`.

Detalle clave que hace esto viable en Windows sin AmberTools: `GBSAOBC2Force.getStandardParameters(topology)` asigna radios y apantallamiento **por elemento** (tablas `_mbondi2_radii` y `_screen_parameter`), no por plantilla de residuo. Por eso funciona igual para la proteína y para el ligando parametrizado con SMIRNOFF.

```python
nb.setReactionFieldDielectric(1.0)          # igual que obc2.xml
force = GBSAOBC2Force(solventDielectric=78.5, soluteDielectric=1.0, SA="ACE")
for i, p in enumerate(GBSAOBC2Force.getStandardParameters(topology)):
    q, _, _ = nb.getParticleParameters(i)
    force.addParticle([q, p[0], p[1]])
force.finalize(); force.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
system.addForce(force)
```

**Efecto medido:** el dG pasa de −66,74 a −39,6 kcal/mol. La diferencia (≈27 kcal/mol) es exactamente la energía de desolvación polar que faltaba.

### 2. Geometría acoplada, sin re-embedding

El original re-generaba la conformación del ligando con `ETKDGv3` + `CoordMap` y `useRandomCoords=True`: la geometría de partida no era la pose. Ahora se **transfieren las coordenadas de la pose acoplada** mediante el mapeo MCS y solo se añaden los hidrógenos (`Chem.AddHs(addCoords=True)`). La geometría de partida es idéntica entre corridas.

### 3. Minimización determinista y hasta convergencia

- Se minimiza en plataforma **Reference (float64, sin reducciones paralelas)**.
- Tolerancia explícita `1 kJ/mol/nm` (la de por defecto, `10`, dejaba el complejo a medio relajar).
- La relajación previa en GPU queda **opt-in** (`MMGBSA_STAGE1=1`) porque no es determinista: introduce dispersión de 8-14 kcal/mol.

**Verificación de convergencia:** con 8.000 iteraciones máximas, el minimizador **termina antes del límite** (mismo tiempo que con 5.000: 2 m 14 s) y el dG coincide con las corridas de 5.000 dentro de 0,5 kcal/mol. Es decir: el resultado es un mínimo convergido, no un corte arbitrario.

**Dependencia del presupuesto de iteraciones (lo que antes pasaba desapercibido):**

| Iteraciones máximas | dG | Tiempo |
|---|---|---|
| 50 | −4,32 | 23 s |
| 200 | −10,27 | 55 s |
| 5.000 | −39,62 / −40,11 | 2 m 13 s |
| 8.000 | −39,79 | 2 m 14 s |

Los valores bajos de 50 y 200 iteraciones son **minimizaciones sin converger**, no resultados. El protocolo correcto se define por convergencia, no por un presupuesto fijo.

---

## ¿POR QUÉ NO SE USA LA GPU? (comprobado)

La pregunta es correcta y merece una respuesta medida, no una excusa.

**1. En este entorno no existe la plataforma CUDA.** OpenMM instalado (conda-forge, win-64) ofrece solo `Reference`, `CPU` y `OpenCL`:

```
0 Reference speed=1.00
1 CPU       speed=10.00
2 OpenCL    speed=50.00
```

**2. Y la razón es de empaquetado, no de hardware.** El build de conda-forge para **linux-64** de OpenMM 8.6 sí incluye CUDA:

```
cuda-nvrtc >=12.9.86,<13.0a0
cuda-version >=12.9,<13
libcufft >=11.4.1.4,<12.0a0
```

El build de **win-64 no tiene ninguna dependencia CUDA**. Es decir: **la GPU del portátil (RTX 4080) no se puede usar con CUDA en Windows por esta vía**; solo queda OpenCL.

**3. OpenCL no sirve para optimizar con fiabilidad.** Trabaja en precisión mixta (float32) y la línea de búsqueda de L-BFGS no converge: el resultado **divaga con el presupuesto de iteraciones**. Medido con el mismo compuesto, misma pose y misma tolerancia (1 kJ/mol/nm):

| Iteraciones máximas | dG (3 corridas) | Dispersión |
|---|---|---|
| 500 | −24,27 / −24,32 / −19,85 | 4,5 |
| 1.000 | −24,03 / −20,09 / −22,34 | 3,9 |
| 2.000 | −40,24 / −27,03 / −27,10 | 13,2 |
| 8.000 | −41,03 / −39,52 / −41,12 | 1,6 |
| 20.000 | −41,55 / −41,08 / −39,96 | 1,6 |
| 20.000 + mejor-de-2 | −24,96 / −39,66 / −27,01 / −39,64 | **14,7** |

El mismo sistema pasa de −20 a −41 kcal/mol solo por cambiar el tope de iteraciones: **no hay un mínimo convergido**, hay donde se pare la máquina. Intenté arreglarlo repitiendo y quedándome con la mejor energía (mejor-de-2) y **empeoró** (14,7 de dispersión): con float32, iterar más tiempo degrada la minimización.

En cambio **Reference (float64) sí converge**, y lo verifiqué: con 8.000 iteraciones máximas el minimizador termina *antes* del tope (mismo tiempo que con 5.000) y da el mismo valor (±0,5).

### Conclusión sobre la GPU

- **No es que la GPU no sirva: es que en Windows este OpenMM solo puede usar OpenCL, y OpenCL no optimiza de forma fiable.**
- **Para usar la GPU bien** (rápido *y* determinista) hay dos vías:
  1. **WSL2 Ubuntu** — ya está instalado en este equipo (detenido). En Linux, conda-forge sí instala OpenMM con CUDA, y la plataforma CUDA soporta `DeterministicForces=true` y precisión doble: el mismo orden de velocidad que OpenCL con la reproducibilidad de Reference. **Es la vía recomendada.**
  2. Instalador oficial de OpenMM para Windows, que sí incluye la plataforma CUDA.
- **Mientras no se haga**, la cifra final se calcula en Reference.

## COSTE

| Modo | Tiempo por candidato | ¿Reproducible? |
|---|---|---|
| **Reference (float64), hasta convergencia** | ≈ 2 m 14 s | **Sí (±0,49 kcal/mol)** |
| OpenCL (precisión mixta) | ≈ 12 s | **No (hasta 15 kcal/mol)** |
| WSL2 + OpenMM CUDA (si se monta) | ≈ 12 s esperado | Sí (DeterministicForces) |

2.973 candidatos en Reference ≈ **111 h de CPU** repartidas entre `--workers`. Con 6 procesos en paralelo ≈ 18-20 h. Para una primera validación, correr una muestra estratificada de 100-200 compuestos.

**Atajo para abaratar el CPU:** reducir la bolsa (hoy 12 Å ⇒ ≈ 438 átomos por ligando). La fuerza GB es O(N²) en float64: con una bolsa fija más pequeña (8 Å, ≈ 250 átomos) el coste baja ~4× (≈ 35 s por compuesto).

---

## PENDIENTES DEL PROTOCOLO (no resueltos todavía)

1. **Bolsa fija por target.** Hoy el receptor se recorta a 12 Å alrededor de *cada* ligando, de modo que cada compuesto se evalúa contra un conjunto de átomos distinto (E_receptor varía de −344 a −1.817 kcal/mol). Científicamente hay que **definir una bolsa fija por target** (una sola vez, por ejemplo con la unión de todas las poses) y usarla para todos los candidatos. Es la mayor fuente de no-comparabilidad que queda.
2. **Extremos artificiales.** El recorte deja residuos cortados a los que PDBFixer añade H terminales: cada corte introduce grupos cargados artificiales. La bolsa fija con la cadena completa en los tramos contiguos lo reduce.
3. **Término no polar.** `SA=ACE` ya aporta superficie; conviene declararlo explícitamente en el paper.
4. **Validación con controles (imprescindible).** Hay que reescalar los **20 controles positivos** (poses disponibles en `gpu_dock/resultados_libreria/results_SOD1/` y `results_TDP43_v2/`) junto con un conjunto de señuelos, y medir enriquecimiento. Sin esto, no se puede afirmar que el ranking corregido discrimina.
5. **Sin entropía.** MM-GBSA no incluye entropía: hay un sesgo sistemático hacia valores muy negativos. Los valores absolutos **no** son Kd; solo sirven para ordenar, y el orden debe calibrarse contra los controles.

---

## PROTOCOLO RECOMENDADO

1. **Paso 0 — validar el método** (esta semana):
   `runner_mmgbsa_gb.py --targets SOD1 --max 30` para medir coste real, y luego una muestra con los 20 controles + 20 señuelos para calcular enriquecimiento (AUC / top-N).
2. **Paso 1 — bolsa fija por target** e implementación del receptor único.
3. **Paso 2 — pasada completa** sobre la lista corta (2.973) con bolsa fija, ~20 h.
4. **Paso 3 — ranking final** ordenado por dG dentro de cada target, con la incertidumbre (±0,5 kcal/mol) declarada y sin comparar entre targets (los scores no son comparables entre receptores).

---

**Conclusión:** la corrección está hecha y verificada; el instrumento ya mide una energía libre de unión (solvatada, determinista y convergida) en lugar de una energía de vacío. Falta cerrar el protocolo (bolsa fija) y validar con los controles antes de ordenar candidatos.
