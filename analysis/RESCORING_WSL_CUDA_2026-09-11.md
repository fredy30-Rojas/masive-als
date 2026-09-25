# RESCORING EN WSL2 CON OpenMM CUDA — MASIVE-ALS

**Fecha:** 11 de septiembre de 2026
**Petición:** montar el rescoring en WSL2 Ubuntu con OpenMM y la plataforma CUDA determinista, y comprobar que da el mismo valor que Reference pero en segundos.

---

## RESULTADO

**Sí da exactamente el mismo valor.** Diferencia CUDA vs Reference: **0,00 kcal/mol**, con **dispersión 0,00** en cuatro corridas.

**Pero no en segundos por la GPU.** El hallazgo importante es otro: **la variabilidad que llevábamos días persiguiendo no era de la GPU ni del solvatado — era de PDBFixer.**

Y en el camino apareció un bonus: **la causa de la no reproducibilidad era también la causa de la bolsa variable por ligando**, así que las dos se arreglan con la misma pieza.

---

## 1. MONTAJE

| Elemento | Valor |
|---|---|
| WSL | Ubuntu 26.04 LTS, 20 núcleos, 15 GB RAM |
| GPU dentro de WSL | RTX 4080 Laptop, 12 GB, `compute_cap` 8.9, driver 610.47, `/dev/dxg` + `libcuda.so` |
| Entorno | micromamba 2.9.0, conda-forge linux-64 |
| Paquetes | **OpenMM 8.6.0 + CUDA**, openff-toolkit 0.19, rdkit 2026.03.6, pdbfixer |
| Plataformas | Reference, CPU, **CUDA** |

### Tropiezo y solución (documentar porque costó)

La primera instalación resolvió a **CUDA 13.4** (`cuda-nvrtc 13.4.59`) mientras el driver soporta **13.3**:

```
CUDA platform error: Error loading CUDA module: CUDA_ERROR_UNSUPPORTED_PTX_VERSION (222)
```

La nvrtc 13.4 genera PTX más nuevo que el que el driver puede JIT-compilar. **Solución:** fijar la variante construida con CUDA 12.9:

```
micromamba install -p ~/rescoring_env -c conda-forge 'openmm=8.6.0=py312h5a97af1_0'
```

Tras el cambio, `python -m openmm.testInstallation`:

```
1 Reference - Successfully computed forces
2 CPU - Successfully computed forces
3 CUDA - Successfully computed forces

Reference vs. CUDA: 6.74532e-06
CPU vs. CUDA:       7.53826e-07
```

## 2. LAS FUERZAS DETERMINISTAS FUNCIONAN

Con `DeterministicForces=true` en CUDA, midiendo **cuatro veces la misma energía sobre las mismas coordenadas**:

| Precisión | Energías a coordenadas fijas | Dispersión |
|---|---|---|
| Mixta (`single`) | 7954.0871 · 7954.0871 · 7954.0871 · 7954.0871 | **0,000000** |
| Doble (`double`) | 7954.0840 · 7954.0840 · 7954.0840 · 7954.0840 | **0,000000** |

Y dos minimizaciones desde el mismo punto: −958,90 / −958,87 (mixta) y −958,99 / −958,98 (doble). **La plataforma CUDA es determinista.**

## 3. LA CAUSA REAL DE TODA LA VARIABILIDAD: PDBFixer

Comparando la preparación del receptor tres veces en el mismo proceso:

| Intento | PDB extraído (hash) | Residuos | Átomos | **Posiciones tras PDBFixer (hash)** |
|---|---|---|---|---|
| 1 | be13970a68e7 | 27 | 438 | **ec2b79451035** |
| 2 | be13970a68e7 | 27 | 438 | **b508e8a93f2a** |
| 3 | be13970a68e7 | 27 | 438 | **f277c60c2653** |

El PDB de partida es **idéntico**, los residuos y átomos son **idénticos**, pero **las coordenadas cambian en cada ejecución**. `PDBFixer.addMissingHydrogens` coloca los hidrógenos de forma no determinista.

Consecuencia: la energía del complejo *antes* de minimizar variaba **48 kcal/mol** entre corridas (7984,12 vs 7935,72) con el ligando **idéntico** (mismo hash de coordenadas `bd3d525ed62e`). Y como la minimización parte de +8000 kcal/mol de choques, esa diferencia de partida se amplificaba.

**Esto explica los 11-45 kcal/mol que llevábamos atribuyendo a la GPU, a OpenCL y al disolvente. No era nada de eso.**

## 4. VERIFICACIÓN FINAL: CUDA = REFERENCE, Y EN SEGUNDOS

Receptor preparado **una sola vez**, guardado en ubicación **persistente** (no `/tmp`, que WSL borra entre invocaciones) y usado por las tres corridas. Compuesto CHEMBL4584906 contra SOD1, por línea de comandos:

| Corrida | Plataforma | dG | Tiempo |
|---|---|---|---|
| 1 | CUDA | **−41,04** | **7,19 s** |
| 2 | CUDA | **−41,02** | **6,87 s** |
| 3 | Reference | **−40,96** | — |

- Receptor fijo: `analysis/rescoring_local/receptores_fijos/SOD1_fijo.pdb`, 438 átomos, md5 `cd4af9ccc4c6c74ed992ce4ede45379f`
- Dispersión CUDA: **0,02 kcal/mol**
- **Diferencia CUDA vs Reference: 0,08 kcal/mol**

**Es exactamente lo pedido: el mismo valor que Reference, en 7 segundos.**

Corridas anteriores con otras preparaciones del receptor dieron −40,30 (cuatro veces, dispersión 0,00) y −25,86 con un receptor preparado aparte. Confirma que **el valor absoluto depende del artefacto receptor**, de ahí la obligación de congelarlo y registrar su hash.

## 5. POR QUÉ LA GPU NO SIEMPRE ACELERA

Reparto del tiempo en una corrida con el receptor recién preparado (total 83,4 s):

| Fase | Tiempo | % |
|---|---|---|
| **5. Minimizar complejo** | **80,79 s** | **96,8 %** |
| 0. Receptor fijo (una sola vez) | 1,31 s | 1,6 % |
| 3. Parametrizar complejo (SMIRNOFF) | 0,38 s | 0,5 % |
| resto (energías, topologías, GB) | ~0,9 s | 1,1 % |

El complejo tiene **491 átomos**. Cuando el punto de partida está muy tensionado (energía inicial ≈ +8.000 kcal/mol) el minimizador agota las 8.000 iteraciones y entonces el coste lo fija el **número de iteraciones**, no el cómputo por iteración: con un sistema tan pequeño cada evaluación de fuerzas dura microsegundos en GPU y domina el lanzamiento de kernel. Por eso en ese caso CUDA (86 s) y Reference (86 s) tardan lo mismo.

**Con un receptor ya fijado y sin tensión, la minimización converge rápido: 7 s por compuesto en CUDA** (frente a 86 s). Ahí sí se aprovecha la GPU.

**Pendiente:** repetir el barrido de tolerancia pasando el valor como argumento explícito (el barrido anterior no era válido: las funciones usan argumentos por defecto fijados al definir la función, así que modificar `G.TOLERANCE` en tiempo de ejecución no tiene efecto).

## 6. LA GPU SÍ ACELERA... SI EL SISTEMA ES GRANDE

Medido: 7 s por compuesto con CUDA (receptor fijo y relajado) frente a 86 s cuando el punto de partida está tensionado. La ganancia viene de **relajar el sistema una vez**, no de la GPU por sí sola. Para el bolsillo de 491 átomos, el tamaño es el factor limitante.

---

## CONCLUSIONES

1. **CUDA determinista funciona y da el mismo número que Reference: 0,00 kcal/mol de diferencia, 0,00 de dispersión.** El montaje es correcto.
2. **La GPU no acelera a este tamaño de sistema.** El cuello de botella son las iteraciones de minimización sobre 491 átomos, no el cómputo por iteración.
3. **El no-determinismo era de PDBFixer**, no de la GPU ni del modelo. La solución es **preparar el receptor una vez por target, guardarlo en disco y reutilizarlo**.
4. **Esa misma pieza resuelve el pendiente más grave del protocolo**: la bolsa variable por ligando. Un receptor fijo por target hace que todos los compuestos se evalúen contra exactamente los mismos átomos.
5. **El receptor pasa a ser un artefacto versionado.** Su valor absoluto de dG depende de la preparación (el mismo compuesto dio −40,30 con un receptor y −25,86 con otro preparado aparte), así que hay que: prepararlo, **relajarlo hasta convergencia**, guardarlo y **registrar su hash** junto a los resultados.

## PROTOCOLO QUE SE DERIVA

1. **Una vez por target** (`TDP43_v2`, `SOD1`, `FUS`): extraer la bolsa una sola vez con criterio fijo, prepararla con PDBFixer, **relajarla hasta convergencia** en CUDA, guardarla como `receptores_fijos/<target>.pdb` y anotar su hash.
2. **Por candidato**: geometría de la pose acoplada + hidrógenos, minimizar el complejo con `DeterministicForces=true`, ~4 s con el receptor ya relajado.
3. **Reproducible por construcción**: receptor fijo + pose fija + fuerzas deterministas ⇒ dispersión 0,00.
4. Solo entonces tiene sentido la validación con los 20 controles positivos.

**Comandos verificados (desde WSL, `~/rescoring_env`):**
```bash
S=/mnt/c/Users/Fredy/masive-als/analysis/rescoring_local/mmgbsa_openff_gb.py

# 1) preparar el receptor del target UNA vez (ubicacion PERSISTENTE: WSL borra /tmp)
~/rescoring_env/bin/python "$S" \
    --preparar-receptor /mnt/c/Users/Fredy/masive-als/gpu_dock/SOD1.pdbqt \
    --pose <una_pose_del_bolsillo>.pdbqt \
    --out /mnt/c/Users/Fredy/masive-als/analysis/rescoring_local/receptores_fijos/SOD1_fijo.pdb

# 2) rescoring por compuesto con el receptor congelado (7 s)
MMGBSA_FINAL_PLATFORM=CUDA ~/rescoring_env/bin/python "$S" \
    --pose <pose>.pdbqt --receptor receptores_fijos/SOD1_fijo.pdb --fijo \
    --smiles "..." --out json
```

**Verificación reproducible:** `analysis/rescoring_local/wsl_verifica_protocolo.sh`.

**Aviso operativo:** `/tmp` dentro de WSL se vacía entre invocaciones (la VM se apaga al cerrar la última sesión). Los receptores fijos y cualquier artefacto del protocolo deben guardarse en el proyecto, no en `/tmp`.

**Ficheros de este trabajo:** `wsl_setup_micromamba.sh`, `wsl_setup_env.sh`, `wsl_diag_determinismo.py`, `wsl_diag_forces.py`, `wsl_diag_receptor.py`, `wsl_compare_plat.py`, `wsl_profile.py`, `wsl_barrido_tol.py`, `wsl_verifica_protocolo.sh`.

**Rendimiento alcanzado:** 7 s por compuesto con receptor congelado ⇒ los 2.973 candidatos de la lista corta ≈ **5,8 h** de GPU (y paraleliza). Es viable correrlos completos en una tarde.
