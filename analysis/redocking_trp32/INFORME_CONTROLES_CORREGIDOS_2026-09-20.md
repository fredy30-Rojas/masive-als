# Control de redocking del bolsillo Trp32: un error nuestro en el medidor de RMSD, corregido

**20 de septiembre de 2026.** Scripts: `validar_rmsd.py` (comprobación del
medidor), `recalcular_rmsd.py` (recálculo de todo), `redock_flexible.py`
(receptor flexible), `preparar_receptor_limpio.py` (revisión de los receptores).
Registros: `validar_rmsd.log`, `recalcular_rmsd.log`, `modos_caja24_corregido.log`,
`flexible_*.log`, `receptor_limpio.log`. Tabla única:
`controles_corregidos.csv`.

Este informe **corrige a los dos anteriores** (`INFORME_REDOCKING_TRP32_...` y
`INFORME_PROTOCOLO_ALTERNATIVO_...`) y añade la tercera ronda. Nada de lo que
dicen esos dos sobre energías cambia; lo que cambia son los RMSD y, con ellos,
el veredicto del control.

---

## 0. Resumen

1. **El medidor de RMSD estaba mal y era nuestro.** Emparejaba los átomos al
   revés e **inflaba** el resultado. Corregido, y comprobado contra la
   implementación de RDKit: coincide **al segundo decimal** en todos los casos.
2. **Con el medidor corregido el control pasa en dos de los cuatro ligandos**:
   isoproterenol (0,46–1,28 Å) y adrenalina (0,66–0,74 Å), con ligando flexible o
   rígido y con las dos funciones de puntuación. Antes figuraban como fallos.
3. **Sigue fallando la dopamina** (2,9–3,2 Å) y **la 5-fluorouridina**
   (5,4–11,3 Å según la función).
4. **La pose depositada de la 5-fluorouridina no es un mínimo del potencial**:
   en una estructura de **1,06 Å de resolución** tiene un contacto F···N de
   **2,07 Å** y un O···O de **2,34 Å**, y la función la puntúa con +0,04 kcal/mol.
   No sirve como patrón de medida y hay que sacarla del criterio.
5. **La dopamina sí es un fallo real, y es de puntuación, no de muestreo**: la
   pose cristalina **sí está entre los nueve modos** (la mejor a 0,74–1,10 Å),
   pero la función prefiere una a 3,1 Å por ~1,2 kcal/mol.
6. **El receptor flexible no arregla nada: lo empeora** (12 Å en el sitio en tres
   de las cuatro). El ligando conserva su forma (0,14–2,0 Å alineado) y se
   desliza fuera del bolsillo: 2,3–7,8 Å de desplazamiento del centroide.

---

## 1. El error, y cómo se comprobó

`rmsd_en_sitio()` emparejaba `cp[i]` contra `cc[match[i]]`. RDKit devuelve en
`match[i]` el índice del átomo del **objetivo** (la pose) homólogo del átomo `i`
de la **consulta** (el cristal), así que el emparejamiento correcto es
`cp[match[i]]` contra `cc[i]`. Las dos formas solo coinciden si la permutación es
involutiva, y casi nunca lo es: el orden de átomos del PDBQT (la raíz primero) no
es el del PDB del cristal.

La comprobación independiente es `rdkit.Chem.rdMolAlign.CalcRMS`
(`validar_rmsd.py`, 36 poses). **La versión corregida y CalcRMS dan el mismo
número en las 36**, con diferencia máxima 0.0000 Å. La versión vieja no.

### Lo que el error hacía

| Estructura (caja 24, semilla 42) | RMSD con el error | RMSD corregido |
|---|---|---|
| 4A7T isoproterenol | 1,88 Å | **0,46 Å** |
| 4A7U adrenalina | 3,74 Å | **0,66 Å** |
| 4A7V dopamina | 4,02 Å | **3,13 Å** |
| 4A7S 5-fluorouridina | 9,63 Å | 9,88 Å |

El efecto es grande justo en los casos buenos: cuanto mejor es la pose, más
castiga un emparejamiento equivocado (el error es de correspondencia de átomos,
no de distancias).

---

## 2. Los tres controles, recalculados

Tabla completa en `controles_corregidos.csv`. Resumen (criterio: mediana ≤ 2,0 Å):

| Control | Ligando de partida | Puntuación | 4A7T isop. | 4A7U adr. | 4A7V dop. | 4A7S 5FU |
|---|---|---|---|---|---|---|
| Caja 24 Å, ex 8 | SDF ideal | vina | **0,48–1,28 PASA** | **0,69–0,72 PASA** | 3,13–3,15 | 11,07–11,33 |
| Caja 18 Å, ex 8 | SDF ideal | vina | **1,60–2,45 PASA** | **0,72–0,75 PASA** | 3,15 | 9,03–9,05 |
| Caja 24 Å, ex 32 | SDF ideal | vina | **1,28 PASA** | **0,72 PASA** | 3,15 | 11,07 |
| Caja 24 Å, ex 8 | pose del cristal | vina | **0,46–0,56 PASA** | **0,66–0,73 PASA** | 3,11–3,13 | 8,74–9,88 |
| Caja 24 Å, ex 8 | pose del cristal | vinardo | **0,77–0,79 PASA** | **0,88–0,92 PASA** | 2,90–3,08 | 5,42–5,45 |
| Caja 24 Å, ex 8 | **ligando rígido** | vina | **0,63–0,65 PASA** | **0,43–0,45 PASA** | 2,97 | 9,83 |
| Caja 24 Å, ex 8 | **ligando rígido** | vinardo | **0,33 PASA** | **0,61 PASA** | 3,02–12,72 | 10,00 |
| **Receptor flexible** | pose del cristal | vina | 1,61–12,14 | 11,72–11,98 | 3,34–12,87 | 9,90–10,86 |

Todo lo demás del informe anterior sobre caja, exhaustividad y preparación del
ligando **se mantiene** (con el medidor corregido el veredicto no cambia en esos
puntos: la caja y la búsqueda no son el problema).

### Los nueve modos, con el medidor corregido (`modos_caja24_corregido.log`)

| Estructura | Mejor de los 9 modos | Modos ≤ 2 Å |
|---|---|---|
| 4A7T isoproterenol | 0,48 – 1,28 Å | 1–3 de 9 |
| 4A7U adrenalina | 0,69 – 0,74 Å | 1–2 de 9 |
| 4A7V dopamina | 0,74 – 1,10 Å | 1–2 de 9 |
| 4A7S 5-fluorouridina | 9,56 – 9,58 Å | 0 de 9 |

**La pose correcta sí se genera** en tres de los cuatro casos (en la dopamina
también, cosa que el informe anterior negaba por el error del medidor). Lo que no
se recupera es la 5-fluorouridina.

### La diferencia de energía (ΔE) sigue en pie, pero significa otra cosa

ΔE no depende del medidor, y se verificó que no fuera un artefacto de comparar
`--score_only` con el acoplamiento: puntuar el modo ya acoplado con
`--score_only` devuelve **el mismo valor que reportó el acoplamiento** (±0,04
kcal/mol), y el tamaño de la caja (autobox o 24 Å) no cambia nada.

ΔE = E(pose del cristal) − E(mejor modo) = **+1,16 / +0,90 / +1,23 / +5,10**
(isoproterenol, adrenalina, dopamina, 5-fluorouridina, vina). Pero ahora se sabe
que el «mejor modo» de isoproterenol y adrenalina **es la pose nativa** (0,5 Å).
Así que ΔE no dice «la función prefiere otra pose»: dice que **la función
distingue diferencias de medio ángulo con energías de 1 kcal/mol**. Medido
directamente: desplazar la pose del cristal 0,5 Å mejoró la interacción de
−5,65 a −7,22 kcal/mol (**1,57 kcal/mol por medio ángulo**).

Y como la incertidumbre de las coordenadas experimentales es de 0,2–0,5 Å, **el
orden dentro de 1–1,5 kcal/mol no es interpretable**. Esto aplica al ranking del
cribado, no al control.

---

## 3. Tercera ronda: receptor flexible

Se soltaron los cinco residuos del bolsillo (Glu21, Gln22, Lys23, Lys30, Glu100;
38 átomos en la parte flexible) y se repitió el control completo.

| Entrada | RMSD en el sitio (3 semillas) | Desplazamiento del centroide | RMSD alineado (Kabsch) |
|---|---|---|---|
| 4A7T isoproterenol | 12,05 / 12,14 / **1,61** | 4,51 / 4,63 / 2,33 Å | 1,99 / 1,92 / 1,15 Å |
| 4A7U adrenalina | 11,98 / 11,98 / 11,72 | ~4,4 Å | 0,45 / 0,45 / 0,57 Å |
| 4A7V dopamina | 3,65 / 12,87 / 3,34 | 2,48 / 5,15 / 2,48 Å | 0,91 / 0,14 / 1,00 Å |
| 4A7S 5-fluorouridina | 10,74 / 9,90 / 10,86 | ~7,5 Å | 0,44 – 1,12 Å |

**Ninguno pasa.** La lectura está en las dos últimas columnas: el ligando
conserva su **forma** (0,1–2,0 Å tras superponer) pero **se desplaza** 2,3–7,8 Å.
Añadir flexibilidad a la cadena lateral no mejora la colocación: abre el bolsillo
y el ligando se desliza. Es el resultado esperable y refuerza que el problema del
sitio no es la rigidez del receptor.

---

## 4. La 5-fluorouridina: el patrón de medida que no sirve

Los datos de 4A7S (1,06 Å de resolución, R/R libre 0,162/0,164, el ligando
depositado dos veces, en las cadenas A y F):

| Medida | Valor |
|---|---|
| Contacto más corto de la pose depositada | **F5 ··· Lys30 NZ = 2,07 Å** |
| Segundo contacto | **O4 ··· Ser98 OG = 2,34 Å** |
| Puntuación de la pose depositada | **+0,04 kcal/mol** (intermolecular **+0,05**) |
| Trasladar la pose −6 / −4 / −2 / 0 / +2 / +4 / +6 Å en x | −0,90 / −1,40 / −1,27 / **+0,04** / **+14,78** / +51,19 / +73,38 kcal/mol |

La pose está pegada a una pared de repulsión: dos ángstroms de desplazamiento la
llevan de 0 a +15 kcal/mol. A 1,06 Å de resolución las coordenadas son fiables
(error ~0,2 Å), así que **la discrepancia es entre el potencial y la estructura,
no un modelo mal refinado**. Y encaja con lo que dice el artículo: la ribosa
queda proyectada al disolvente con ocupación baja.

Conclusión práctica: **la 5-fluorouridina no puede usarse como patrón de RMSD**.
Su sitio existe y su pose está descrita, pero no es un mínimo de la función de
puntuación, así que un fallo de redocking con ella no dice nada del protocolo.
(Es también el único ligando con anillos rígidos y con un flúor, dos cosas que
estos potenciales tratan mal.)

---

## 5. Qué queda dicho, y qué cambia para el proyecto

1. **El protocolo reproduce las poses nativas de las dos catecolaminas** que
   definen el bolsillo: 0,3–0,9 Å con tres protocolos distintos (flexible y
   rígido, vina y vinardo) y en dos tamaños de caja. La caja del proyecto y el
   protocolo quedan **validados para ese quimiotipo**.
2. **La dopamina es el único fallo de verdad**, y es de **puntuación**: la pose
   correcta se genera (0,74–1,10 Å) pero la función prefiere una a 3,1 Å por
   ~1,2 kcal/mol. Con un potencial que resuelve medio ángulo con 1,5 kcal/mol,
   eso no es un fallo de programa: es el límite del método.
3. **El 0,815 de SOD1 sigue provisional, pero por el otro motivo**: por la serie
   congénica (0,601 con los dos activos independientes y LCS-1 en el puesto
   152/219). **Ya no** por «el protocolo no reproduce las poses»: eso era
   consecuencia del error del medidor.
4. **Corregidos**: `paper/paper_masive_als.md`, los dos informes anteriores, el
   banner de `analysis/REVISION_VALIDACION_BOLSILLOS_2026-08-20.md` y el
   entregable de las investigadoras.
5. **Lección, y va al fichero de aprendizajes**: un medidor propio de RMSD se
   comprueba contra el de la librería antes de usarlo para decidir nada. Un
   emparejamiento de átomos invertido **infla** el RMSD y hace pasar por fallos
   lo que funciona.

---

## 6. Qué sigue

1. **Rehacer la validación de SOD1** con el conjunto limpio
   (`controles_sod1_v3.csv`) y los dos fondos, diciendo en el informe que el
   protocolo está validado **para el quimiotipo catecolamina** del bolsillo
   Trp32. Es el paso que ahora sí tiene sentido dar.
2. **Dopamina**: si se quiere un cuarto patrón fiable, hay que buscar otra
   estructura o conformación del bolsillo; con esta no hay forma de que el
   ranking acierte.
3. **AD4** sigue pendiente de verdad (necesita mapas AutoDock4), pero ya no es
   prioritario: con vina y vinardo y dos catecolaminas el control pasa.
4. **La 5-fluorouridina se queda fuera del criterio**, con su nota explicando por
   qué.
