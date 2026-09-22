# MM-GBSA por estratos de tamaño: ordena mejor, pero no quita el sesgo de tamaño

**22 de septiembre de 2026, 16:32.** Números **definitivos**: la repetición de los 90
ligandos con el método único arreglado (v5) terminó en Oracle sin ningún error
(57,2 min, log `runner_estratos_stdout_v6.log`) y este informe ya no mezcla métodos.
La revisión de mediodía con los números provisionales queda resumida en §0 y sus
tablas se conservan en el histórico del repositorio.

**Scripts:** `analysis/preparar_mmgbsa_estratos.py`, `analysis/mmgbsa_runner_estratos.py`
(en Oracle), `analysis/rescoring_estrategia_estratos.py`,
`analysis/rescoring_mmgbsa_robusto.py` (v5) · **Datos:**
`analysis/rescoring_estrategia_estratos.csv` (lista completa, N=513),
`rescoring_estrategia_estratos_muestra.csv`,
`mmgbsa_estratos/resultado_oracle.csv` (**v5 limpio, 90/90**, este es el bueno),
`mmgbsa_estratos/resultado_oracle_v4_v5_mixto.csv` (mezcla vieja, solo para comparar).

Responde a la pregunta que quedó abierta en `INFORME_VALIDACION_SOD1_LIMPIA_2026-09-21.md`:
si la ordenación por energía no sabe ordenar porque su energía es casi una función del
tamaño, **¿el rescoring físico MM-GBSA ordena mejor, y le queda el mismo vicio?**

---

## 0. Lo que cambió desde la revisión de mediodía

1. **La tabla publicada mezclaba dos métodos y ya no.** Los 82 ligandos buenos de la
   víspera se calcularon con la versión que leía la pose con obabel; los 8 que fallaban
   se repitieron con la v5 (lectura de la pose con el mapa `REMARK SMILES IDX` que
   escribe Vina). El mismo ligando (`ACT_isoproterenol`) daba **−15,65** con la vieja y
   **−4,95** con la nueva. La repetición completa de los 90 con la v5 terminó hoy a las
   12:21 UTC (90 ok, 0 errores). El CSV viejo se conserva como
   `mmgbsa_estratos/resultado_oracle_v4_v5_mixto.csv`.
2. **El cálculo no es repetible al kcal, y el ruido real es mayor de lo que se estimó
   al mediodía.** Comparaciones del mismo ligando entre corridas idénticas del mismo
   método v5: `DECM_CHEMBL4435214` **−28,02 / −33,73**; `ACT_isoproterenol`
   **−13,52 / −4,95** (esta segunda pareja procede de una preparación intermedia del
   receptor cuya energía absoluta difiere en cientos de kcal; aun descartándola queda la
   primera, 5,7 kcal). La causa está localizada: **PDBFixer coloca los átomos que faltan
   con azar** y cada proceso prepara su propio receptor. **Consecuencia: en esta tabla,
   diferencias por debajo de ~2 kcal/mol no son señal, son ruido.**
3. **El "MM-GBSA" no lleva cargas en el ligando.** El script las pone a cero a propósito
   (se salta AM1-BCC para que la tanda cupiera en el tiempo). Un informe intermedio dijo
   "GAFF 2.11 + AM1-BCC": **era falso**, queda corregido. El dG se parece mucho a una
   energía de van der Waals más superficie, lo que explica (§3) que su pendiente con el
   tamaño sea la mayor de las cuatro funciones. Un MM-GBSA con cargas de verdad es el
   escalón siguiente, no un detalle.

Los 8 fallos de la víspera siguen arreglados y verificados (4 de mapa de átomos,
2 de bytes NUL en la pose, 2 de residuo THR aislado tras el recorte). **De los 90
intentados, 90 dan número y ninguno error.**

## 1. Lo que se montó

| Pieza | Qué es |
|---|---|
| Receptor | `SOD1_limpio.pdb` (dímero A–H, sin aguas ni copias de 1HL5), el mismo de la validación |
| Poses | las del acoplado original (`validacion_SOD1_v5/out`), **no** se volvieron a acoplar |
| Positivos | los **11 de unión medida** de `verdad_de_referencia.csv` |
| Fondo | **79 señuelos**: 66 emparejados por tamaño (±2 átomos pesados con un positivo) + 13 repartidos por cuantiles de tamaño |
| Método | `rescoring_mmgbsa_robusto.py` v5: GAFF 2.11 **con cargas del ligando a cero**, openmm, OBC2, una pose por ligando, minimización local de 200 pasos |
| Máquina | Oracle, 4 núcleos, 3 procesos en paralelo (~38 s por ligando de media; los 90 en 57,2 min) |

**Cobertura: 11 positivos + 79 señuelos, 90 de 90.** Cada positivo tiene entre 10 y 21
señuelos de su mismo tamaño dentro de la muestra: el criterio emparejado se puede medir
con apoyo, no con anécdotas.

## 2. La comparación, sobre los MISMOS ligandos

Comparar un AUC medido sobre 90 ligandos emparejados con otro medido sobre 513 no vale:
la muestra emparejada es más fácil para cualquier función. Así que las cuatro funciones
se pasaron por los mismos 90 ligandos:

| Función | AUC crudo | AUC residual | Corr. con tamaño | EF5 % | Quimiotipos |
|---|---|---|---|---|---|
| Acoplado (Vina, control) | 0,457 | 0,490 | −0,826 | 1,82 | 3 de 8 |
| Rescoring Vina (`--score_only`) | 0,461 | 0,514 | −0,823 | 1,82 | 3 de 8 |
| **Rescoring Vinardo** | 0,461 | 0,504 | −0,645 | 1,82 | **6 de 8** |
| Rescoring MM-GBSA | 0,453 | 0,495 | −0,487 | 1,82 | 4 de 8 |

*(Con 90 ligandos el EF5 % topa a 4 ligandos y no discrimina; se deja por completitud,
no como medida.)*

Y por estratos de tamaño, que es lo que se quería ver:

| Estrato | Positivos | Señuelos | MM-GBSA | Vinardo | Acoplado |
|---|---|---|---|---|---|
| 0–12 átomos | 2 | 12 | 0,750 | 0,833 | 0,667 |
| 13–17 átomos | 3 | 18 | 0,204 | 0,481 | 0,296 |
| 18–22 átomos | 3 | 19 | 0,368 | 0,211 | 0,263 |
| 23–27 átomos | 1 | 16 | 0,875 | 1,000 | 1,000 |
| 28–99 átomos | 2 | 14 | 0,536 | 0,536 | 0,607 |
| **Estratificado (ponderado)** | | | **0,515** | 0,578 | 0,534 |

## 3. El sesgo de tamaño NO se va: sigue siendo el más fuerte

Dentro del fondo de la muestra, en kcal/mol por átomo pesado:

- Vina (acoplado): **−0,099**;
- Rescoring Vina: −0,101;
- Vinardo: −0,078;
- **MM-GBSA: −0,498.**

El MM-GBSA conserva **el término de tamaño más fuerte de las cuatro funciones**, cinco
veces el de Vinardo en kcal/mol por carbono. Y con el ligando sin cargas (§0.3) ese
−0,498 por carbono es casi lo que se espera de una energía de van der Waals: parte del
vicio de tamaño es consecuencia del atajo. La mejora de mediodía (0,541/0,608) se
desinfló al quitar la mezcla de métodos: **el "ordenaba mejor" era en parte el método
viejo inflando a los positivos, no virtud del MM-GBSA.**

## 4. El veredicto (cambia respecto al informe provisional)

**El MM-GBSA deja de ser la mejor segunda vuelta; ese puesto pasa a Vinardo.** Con el
método único:

| Quimiotipo | Positivo | dG MM-GBSA | Átomos | ¿Gana a los de su tamaño? |
|---|---|---|---|---|
| aminoalcohol naftalénico | 946 | −23,63 | 23 | **sí** |
| anilina | ZZT | −17,56 | 10 | **sí** |
| benzisoxazol-piperidina | K4I | −28,44 | 31 | **sí** |
| catecolamina | isoproterenol / dopamina | −4,95 / −11,38 | 15 / 11 | **no** (antes sí) |
| quinazolina | 12I | −16,30 | 17 | **no** (antes sí) |
| quinazolina-CF₃ | ZO0 | −25,29 | 22 | **sí** |
| fenantridinona | 6B3 | −15,52 | 33 | no |
| nucleósido | 5-fluorouridina | −12,68 | 18 | no |

**El MM-GBSA recupera 4 de 8** (antes parecía 6). **Vinardo recupera 6 de 8** con la
mitad de pendiente de tamaño (−0,078 frente a −0,498): las mismas quinazolinas que el
MM-GBSA de mediodía, y además la catecolamina y la fenantridinona que este pierde.

**El dato que hay que dar entero:** con 90 ligandos, un solo positivo por estrato y un
ruido inter-corrida de ±1 a 6 kcal, la diferencia entre "4 de 8" y "6 de 8" está dentro
de lo que el propio montaje no puede resolver. Lo que sí sostiene la tabla: **ninguna
función ordena en crudo** (AUC 0,45–0,46), **el MM-GBSA es el que más sesgo de tamaño
arrastra**, y la pregunta emparejada por tamaño sigue siendo la única que da señal.

## 5. Lo que NO dice, en el mismo párrafo

1. **No se publica nada con esto.** Once positivos dan para decidir si una idea se sigue
   o se tira, no para reportar un número.
2. **Los estratos son finos:** de 1 a 3 positivos en cada uno. El estratificado pondera
   por tamaño de estrato, pero el orden dentro de los chicos es anécdota (el "1,000" del
   23–27 es **un** positivo).
3. **El MM-GBSA tiene el vicio de tamaño más fuerte de todas** (§3), y **es un MM-GBSA
   de un punto** —una pose por ligando, sin muestreo— **y sin cargas en el ligando**
   (§0.3). No es el protocolo de referencia; los dG no se comparan con la literatura.
4. **El ruido es del orden del efecto o mayor** (§0.2): ±1 a 6 kcal/mol entre corridas
   idénticas. Cualquier conclusión que dependa de menos de ~2 kcal/mol de esta tabla no
   está sostenida.
5. **El EF5 % no es defendible con 90 ligandos** (§2).
6. **La comparación MM-GBSA v5 tiene dos limitaciones añadidas que hay que nombrar:**
   la preparación del receptor varía entre procesos (mismo script, receptor distinto),
   y la lectura v5 de la pose cambia el sistema respecto a la v4 vieja (el
   `ACT_isoproterenol` v4 daba −15,65; ningún v5 lo reproduce).

## 6. Qué significa para el proyecto

- **Vinardo sustituye al MM-GBSA como segunda vuelta del embudo**, dentro de cada
  estrato de tamaño: el ranking de Vina elige la lista, y **Vinardo desempata dentro del
  estrato**. Es gratis (score de Vina 1.2.3, sin openmm ni GAFF), no tiene el ruido de
  preparación del receptor, y su sesgo de tamaño es el menor de las cuatro funciones.
  MM-GBSA sobre la lista completa vuelve a ordenar por tamaño, con cinco veces esa
  pendiente.
- **Ninguna de las cuatro funciones ordena sola.** Las cuatro están por debajo o al
  borde del azar en crudo (0,435–0,463) y todas pasan en la pregunta emparejada.
  El problema del proyecto no es "qué función se usa", es **cómo se pregunta**.
- **Siguiente paso concreto, en orden:**
  1. **Fijar la semilla/preparación del receptor** (prepararlo una vez, guardarlo y
     reutilizarlo) para que el ruido del §0.2 deje de ser del tamaño del efecto. Es el
     cambio más barato y desbloquea las demás comparaciones.
  2. **Probar el MM-GBSA con cargas AM1-BCC de verdad** sobre los mismos ligandos: si el
     término de tamaño baja, parte del sesgo era del atajo; si no baja, es de la medida.
  3. Solo después de eso, **repetir la validación en TDP-43** con los 7 positivos de
     unión medida y el fondo duro de R-BIND 2.0 (`validar_diana_limpia.py` ya está
     preparado), con Vinardo como desempate por estratos.
  4. Y solo si todo lo anterior sale bien: pensar en la lista de candidatos. Hasta hoy,
     **sigue sin publicarse**.
