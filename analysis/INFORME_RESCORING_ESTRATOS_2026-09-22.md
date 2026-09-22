# MM-GBSA por estratos de tamaño: ordena mejor, pero no quita el sesgo de tamaño

**22 de septiembre de 2026, 16:32.** Números **definitivos**: la repetición de los 90
ligandos con el método único arreglado (v5) terminó en Oracle sin ningún error
(57,2 min, log `runner_estratos_stdout_v6.log`) y este informe ya no mezcla métodos.
La revisión de mediodía con los números provisionales queda resumida en §0 y sus
tablas se conservan en el histórico del repositorio.

**Adenda de la tarde (18:15):** verificado con tres filas repetidas que la preparación
semillada del receptor (v6) elimina el ruido de PDBFixer — la tanda de los 90 ya es
reproducible bit a bit (§0.2). No hace falta repetir la tanda.

**Adenda 2 (17:45):** el experimento "MM-GBSA con cargas AM1-BCC de verdad" ya estaba
hecho sin saberlo: el generador GAFF las calcula siempre y la tabla de los 90 las lleva
(§0.3). El sesgo de tamaño es de la medida, no del atajo. No hay que repetir nada.

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
2. **El ruido entre preparaciones del receptor está localizado, sembrado y VERIFICADO.**
   El origen medido paso a paso: `addMissingAtoms` con semilla ya era determinista
   (0,000 Å), pero **`addMissingHydrogens` coloca los hidrógenos con el generador
   aleatorio de Python** (302 de 875 átomos se movían >0,01 Å) y cada proceso preparaba
   su propio receptor: `DECM_CHEMBL4435214` dio **−28,02 / −33,73** y `ACT_isoproterenol`
   **−15,65 / −13,52 / −7,69** entre corridas de las versiones sin sembrar (v4/v5
   temprana). El arreglo (script v6): semilla fija `SEMILLA_RECEPTOR=20260922` en
   `addMissingAtoms`, `random` y `numpy.random` sembrados antes de
   `addMissingHydrogens`, y un solo hilo de OpenMM (con varios, la suma de fuerzas de
   la CPU lleva otro orden cada vez). **Verificación del 22 sep por la tarde:** tres
   filas repetidas fuera de la tanda (`ACT_adrenalina`, `DECM_CHEMBL9250`,
   `DECM_CHEMBL4435214` — este último, el campeón histórico del ruido) reproducen el
   dG **y las tres componentes de energía al mismo decimal** de `resultado_oracle.csv`.
   La tanda de los 90 corrió ya con la preparación fijada: repetirla daría lo mismo.
   **Consecuencia: dentro de esta tabla, los números son reproducibles bit a bit; el
   ruido de ±1–6 kcal solo aplica al comparar con las tandas viejas sin sembrar.**
3. **CORRECCIÓN de la noche: el ligando SÍ lleva cargas AM1-BCC reales.** El script pone
   las cargas a cero creyendo saltarse el AM1-BCC ("Set zero charges to skip slow
   AM1-BCC calculation"), pero **`GAFFTemplateGenerator` siempre calcula AM1-BCC al
   construir la plantilla** (openmmforcefields 0.16.0, `template_generators.py` línea
   594: `assign_partial_charges(partial_charge_method="am1bcc", normalize=True)`) y
   sobrescribe lo que el script hubiera puesto. Verificado por dos caminos: (a) el
   `NonbondedForce` del sistema montado con "cargas a cero" tiene la misma carga máxima
   (0,827 e) que el montado con AM1-BCC explícito; (b) una variante del script con
   AM1-BCC explícito (`rescoring_mmgbsa_bcc.py`, vía `AmberToolsToolkitWrapper`)
   reproduce **bit a bit** tres filas de la tanda (`ACT_adrenalina` −10,65,
   `ACT_naftalenoaminoalcohol_946` −23,63, `DECM_CHEMBL9250` −18,12, incluidas las tres
   componentes de energía). Un informe intermedio leyó esto al revés dos veces: primero
   dijo "GAFF + AM1-BCC" (falso en la intención del script, cierto en la física), luego
   dijo "cargas a cero, el dG es casi van der Waals" (falso en la física). Lo segundo
   era la conclusión equivocada: **la electrostática del ligando estuvo dentro desde el
   principio, y el atajo nunca ahorró el tiempo que decía ahorrar.**

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
| Método | `rescoring_mmgbsa_robusto.py` v6: GAFF 2.11 **con cargas AM1-BCC reales del ligando** (las pone `GAFFTemplateGenerator`, ver §0.3), openmm, OBC2, una pose por ligando, minimización local de 200 pasos, **preparación del receptor semillada** (`SEMILLA_RECEPTOR=20260922`, un hilo de OpenMM por cálculo) |
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
veces el de Vinardo en kcal/mol por carbono. La hipótesis cómoda —que ese vicio viniera
del atajo de las cargas— **cayó con la corrección del §0.3**: las cargas AM1-BCC
estuvieron dentro desde el principio, así que la pendiente es del propio OBC2/GAFF con
cargas, de la medida y no del atajo. La mejora de mediodía (0,541/0,608) se desinfló al
quitar la mezcla de métodos: **el "ordenaba mejor" era en parte el método viejo
inflando a los positivos, no virtud del MM-GBSA.**

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
   de un punto** —una pose por ligando, sin muestreo—. No es el protocolo de
   referencia; los dG no se comparan con la literatura.
4. **Determinista no significa exacto** (§0.2): la tabla ya es reproducible bit a bit,
   pero sigue siendo un MM-GBSA de una pose y sin cargas en el ligando; su error frente
   a la realidad no lo mide la reproducibilidad.
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
  1. ~~Fijar la semilla/preparación del receptor~~ **HECHO y verificado el 22 sep**
     (§0.2): la preparación va semillada, un hilo, y tres filas repetidas reprodujeron
     la tanda a plena precisión. El ruido de preparación ya no ensucia las comparaciones.
  2. ~~Probar el MM-GBSA con cargas AM1-BCC de verdad~~ **RESUELTO la noche del 22 sep
     (§0.3)**: ya las llevaba — el generador GAFF las calcula siempre. Respuesta a la
     pregunta: la pendiente de tamaño no baja (los números son idénticos bit a bit con
     AM1-BCC explícito), así que **el sesgo es de la medida, no del atajo**. No hay
     ninguna variante "con cargas" pendiente: esta tabla ya es ese experimento.
  3. Solo después de eso, **repetir la validación en TDP-43** con los 7 positivos de
     unión medida y el fondo duro de R-BIND 2.0 (`validar_diana_limpia.py` ya está
     preparado), con Vinardo como desempate por estratos.
  4. Y solo si todo lo anterior sale bien: pensar en la lista de candidatos. Hasta hoy,
     **sigue sin publicarse**.
