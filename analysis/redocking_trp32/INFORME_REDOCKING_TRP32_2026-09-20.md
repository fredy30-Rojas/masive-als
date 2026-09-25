# El bolsillo Trp32 de SOD1: control de redocking y verdad de los sitios

**20 de septiembre de 2026, versión corregida.** Scripts en
`analysis/redocking_trp32/`: `redock_trp32.py` (control),
`comparar_sitios.py` (geometría de los sitios), `analizar_modos.py` (modos de
Vina), `conformacion_cristal.py` (conformación del ligando). Registros:
`comparar_sitios.log`, `redocking.log`, `redocking_ex32.log`,
`redocking_caja18.log`, `modos_caja24.log`, `conformacion.log`, más los
`resultado_redocking_caja*.json`. No se usó GPU ni se acopló nada nuevo de
la librería: solo ligandos nativos, con Vina en CPU.

> **AVISO — este informe corrige a su primera versión.** La primera versión
> concluía que los cuatro ligandos «de Trp32» no estaban todos en Trp32. **Era
> falso, y el error era mío**: el script cogía «la primera copia» del ligando de
> cada entrada, y las estructuras traen copias en varios sitios de la superficie.
> En 4A7U la adrenalina está dos veces, y la primera copia (A1000) está en un
> sitio secundario. La copia buena es la A1001. Corregido el criterio, los cuatro
> ligandos **sí** están en el bolsillo de Trp32, como dice el artículo original
> (Wright et al., 2013, *Nat Commun* 4:1758). Lo que sigue ya está con el
> criterio corregido.

> **AVISO 2 — LOS RMSD DE ESTE INFORME ESTÁN INFLADOS Y SU VEREDICTO NO VALE.**
> El medidor de RMSD emparejaba los átomos al revés (corregido el 20 de
> septiembre; ver `INFORME_CONTROLES_CORREGIDOS_2026-09-20.md` y
> `validar_rmsd.py`). Con el medidor corregido, **isoproterenol y adrenalina
> redockean bien** (0,48–1,28 Å y 0,66–0,75 Å) y solo fallan dopamina
> (3,1 Å) y 5-fluorouridina (cuya pose depositada, además, no es un mínimo del
> potencial: tiene un contacto F···N de 2,07 Å en una estructura de 1,06 Å).
> Las afirmaciones sobre **energías** de este informe (ΔE, los nueve modos) no
> cambian; las que dependen del RMSD sí.

---

## 0. Por qué se hizo esto

La auditoría del 19 de septiembre dejó dos cabos: la huella de interacción 3D dio
**0,000** con los positivos independientes de SOD1 porque LCS-1 y PRG-A01 no
comparten los residuos que tocan (0,038), y el paso propuesto era averiguar
**dónde se coloca de verdad cada ligando**. Para eso el proyecto quería usar como
patrón de medida los cuatro ligandos co-cristalizados de Trp32.

Pero los cuatro no salen del receptor del proyecto. El receptor del cribado es
**1HL5**; los cuatro ligandos salen de cuatro entradas del mutante **I113T**:
**4A7S** (5-fluorouridina), **4A7T** (isoproterenol), **4A7U** (adrenalina) y
**4A7V** (dopamina). Y el proyecto nunca había hecho el control que dice si el
protocolo es capaz de reproducir esas poses: un **redocking**. Eso es lo que se
hizo.

---

## 1. La caja del proyecto sí está donde dice

El centro de caja del cribado, **(46,5, 80,0, 73,3)** sobre `receptor_SOD1.pdb`
(1HL5), está de verdad sobre el bolsillo de Trp32:

| Residuo (cadena A de 1HL5) | Distancia mínima a la caja |
|---|---|
| Lys30 | 1,1 Å |
| Glu21 | 3,3 Å |
| Val29 | 3,8 Å |
| Val31 | 4,0 Å |
| **Trp32** | **4,1 Å** |
| Gln22 | 5,2 Å |
| Phe20 | 5,2 Å |
| Glu100 | 5,3 Å |

Es el mismo conjunto que rodea al isoproterenol en el cristal. La caja no está
mal colocada; lo que estaba sin comprobar es el protocolo que puntúa dentro.

---

## 2. Los cuatro ligandos sí comparten el bolsillo de Trp32

Se superpusieron las cuatro estructuras por los 153 Cα de la cadena A
(concordancia excelente: RMSD de Cα de 0,17–0,25 Å) y se llevó cada ligando al
sistema común, **eligiendo en cada entrada la copia que ocupa el bolsillo de
Trp32** (la que más residuos del conjunto {Glu21, Gln22, Lys23, Pro28, Val29,
Lys30, Trp32, Glu100} toca a 4,5 Å; el script imprime qué copia elige y qué
copias descarta).

### Copias del ligando en cada entrada

| Entrada | Copias en la cadena A | Elegida |
|---|---|---|
| 4A7S | 1 | A1159 (toca Glu21, Lys30 y Trp32 del conjunto; y Ser98 fuera del conjunto) |
| 4A7T | 3 (A1000, A1001, A1002) | A1000 (toca los 8 residuos) |
| 4A7U | 2 (A1000, A1001) | A1001 (toca los 8); la A1000 está en un sitio secundario |
| 4A7V | 1 | A1000 (toca 7) |

### Residuos contactados (≤ 4,5 Å, en el sistema de 4A7T)

| Estructura | Ligando | Residuos |
|---|---|---|
| 4A7T | isoproterenol | Glu21, Gln22, Lys23, Pro28, Val29, Lys30, **Trp32**, Glu100 |
| 4A7U | adrenalina | Glu21, Gln22, Lys23, Pro28, Val29, Lys30, **Trp32**, Glu100 |
| 4A7V | dopamina | Glu21, Gln22, Lys23, Pro28, Val29, Lys30, Glu100 |
| 4A7S | 5-fluorouridina | **solo** Lys30, Trp32, Ser98 |

### Superposición de esos contactos

| Comparación | Jaccard de residuos contactados |
|---|---|
| isoproterenol vs adrenalina | **1,00** |
| isoproterenol vs dopamina | **0,875** |
| isoproterenol vs 5-fluorouridina | 0,222 |
| ayer: LCS-1 vs PRG-A01 (positivos de SOD1) | 0,038 |

**Este es el dato que rehabilita la métrica de ayer.** Cuando dos ligandos
comparten bolsillo de verdad, la medida de contactos lo ve perfectamente
(0,875–1,00). El 0,038 de ayer no venía de que la métrica fuese mala: venía de
que LCS-1 y PRG-A01 **no ocupan el mismo sitio**, o al menos el acoplamiento no
los coloca en el mismo sitio. La métrica de contactos sirve, con la condición de
aplicarla dentro de un conjunto de sitio verificado.

### La 5-fluorouridina, con matiz

Su centroide queda a 7,2 Å del isoproterenol y toca solo 3 residuos, pero **sí
está en Trp32**: el artículo describe apilamiento aromático del fluorouracilo con
el Trp32 y un puente de hidrógeno con Ser98, y que **la ribosa queda libre,
proyectada al disolvente y con ocupación baja**. O sea: comparte el bolsillo pero
no la extensión. Medida por centroide parece otro sitio; medida por contactos se
ve que es el mismo ancla con una cola al disolvente.

---

## 3. El control de redocking: el protocolo no reproduce las poses

Es el estándar mínimo de cualquier trabajo de acoplamiento y el proyecto nunca lo
había hecho. Se acopló **cada ligando nativo en su propia estructura** (caja
centrada en el centroide del cristal, exhaustividad 8, tres semillas) y se midió
el RMSD **sin superponer** (las dos poses ya están en el sistema del receptor),
probando los emparejamientos por simetría y tomando el mínimo. Criterio:
**RMSD ≤ 2,0 Å**.

Sin semilla fija el resultado cambiaba entre corridas (isoproterenol 1,89 Å una
vez y 2,55 Å la siguiente), así que se fijaron las del proyecto: **42, 2026, 777**.

### Caja de 24 Å

| Estructura | Ligando | RMSD por semilla | Mejor | Mediana | ¿Pasa? |
|---|---|---|---|---|---|
| 4A7T | isoproterenol | 2,16 / 2,54 / **1,87** | 1,87 Å | 2,16 Å | 1 de 3 (al filo) |
| 4A7U | adrenalina | 3,79 / 3,77 / 3,79 | 3,77 Å | 3,79 Å | 0 de 3 |
| 4A7V | dopamina | 4,06 / 4,04 / 4,05 | 4,04 Å | 4,05 Å | 0 de 3 |
| 4A7S | 5-fluorouridina | 10,89 / 11,20 / 11,15 | 10,89 Å | 11,15 Å | 0 de 3 |

### Con una caja más pequeña (18 Å), para descartar el tamaño

| Estructura | Mejor RMSD | ¿Pasa? |
|---|---|---|
| 4A7T | 2,38 Å | 0 de 3 |
| 4A7U | 3,77 Å | 0 de 3 |
| 4A7V | 4,05 Å | 0 de 3 |
| 4A7S | 9,43 Å | 0 de 3 |

La caja no es la causa.

---

## 4. ¿Por qué falla? Tres comprobaciones que lo acotan

### a) No es falta de búsqueda

Con exhaustividad 32 (cuatro veces la del proyecto) los resultados son los
mismos: 4A7S 10,90 Å, 4A7U 3,79 Å, 4A7V 4,05 Å. Más búsqueda no arregla nada.

### b) La pose correcta no está ni entre los nueve modos

Se midió el RMSD de **cada uno** de los nueve modos que Vina guarda
(`analizar_modos.py`). Si la pose correcta estuviera entre ellos y solo mal
puntuada, el fallo sería de puntuación; si no está, el fallo es de colocación.

| Estructura | Mejor de los 9 modos | Modos ≤ 2 Å |
|---|---|---|
| 4A7T | 1,69–1,92 Å | 1–2 de 9 |
| 4A7V | 3,41–3,56 Å | 0 de 9 |
| 4A7U | 3,57–3,62 Å | 0 de 9 |
| 4A7S | 9,36–9,47 Å | 0 de 9 |

La pose cristalográfica **no se genera**. En tres de cuatro casos el programa
encuentra otro mínimo con mejor puntuación y se queda ahí. Y las afinidades
apenas se mueven entre la pose buena y una a 11 Å (todas entre −4,9 y −5,4
kcal/mol): **la energía no distingue** en este bolsillo.

### c) Tampoco es la preparación del ligando

Antes de culpar al acoplamiento había que descartar que el ligando de entrada
estuviera en una conformación imposible de acercar. Se generaron 300 conformeros
por ligando y se buscó el mínimo RMSD alineado con la pose cristalográfica
(`conformacion_cristal.py`):

| Estructura | 300 conformeros | Un solo conformero (el que se le dio a Vina) |
|---|---|---|
| 4A7T isoproterenol | 0,33 Å | 1,34 Å |
| 4A7V dopamina | 0,23 Å | 1,16 Å |
| 4A7U adrenalina | 0,06 Å | 1,33 Å |
| 4A7S 5-fluorouridina | 0,41 Å | 1,87 Å |

La conformación del cristal **es alcanzable** (0,06–0,41 Å). En las tres
catecolaminas la diferencia es torsional, y el acoplamiento sí mueve torsiones.
Así que la preparación no es la causa. (En la 5-fluorouridina queda la duda
razonable del pucker de la ribosa, que Vina no muestrea porque los anillos son
rígidos; ese caso hay que leerlo con más cautela que los otros tres.)

### El diagnóstico, entonces

**El óptimo global de la función de puntuación de Vina en este bolsillo no es la
pose cristalográfica, y la pose correcta no aparece ni entre los nueve modos.**
No es el tamaño de la caja, ni la profundidad de búsqueda, ni la conformación del
ligando. Es el método tal como se está usando: receptor rígido, programa Vina y
esa función de puntuación.

---

## 5. Consecuencias para el proyecto

1. **El «patrón de medida cristalográfico» existe** (los cuatro ligandos están en
   Trp32, y tres de ellos comparten los mismos residuos), pero **el protocolo no
   lo reproduce**. Un patrón que el método no alcanza no valida un ranking.
2. **Ninguna lista de Trp32 se puede presentar como ordenada** mientras esto no se
   arregle. Aplica al AUC 0,815 de SOD1 —que además ya arrastra el problema de la
   serie congénica— y a cualquier candidato salido de esa caja.
3. **El paso 6 de ayer queda explicado**: LCS-1 no cae donde debería porque el
   acoplamiento no coloca bien ni siquiera a los ligandos cuya pose está resuelta
   por cristalografía. El problema es anterior a cualquier métrica.
4. **Queda una medida nueva y útil**: la superposición de residuos contactados
   entre positivos (0,875–1,00 dentro de Trp32). Sirve como test de coherencia del
   conjunto de positivos de cualquier diana, antes de gastar horas en métricas.
5. **La adrenalina aparece además en un sitio secundario** (copia A1000 de 4A7U,
   a 15,5 Å, tocando Asn26, Gly27, Pro28, Ser102, Val103, His110). El artículo no
   lo reporta; probablemente es un contacto cristalino. Se anota, no se usa.

---

## 6. Qué sigue, y en qué orden

1. **Probar el acoplamiento en el receptor cristalográfico del mutante (4A7T) en
   vez de 1HL5**, con el mismo control. Si tampoco reproduce, el problema está en
   el método y no en el receptor del proyecto.
   > **HECHO (20 sep, ver `INFORME_PROTOCOLO_ALTERNATIVO_2026-09-20.md`)**: el
   > receptor de todo el control era ya el cristalográfico (4A7T y compañía), no
   > 1HL5. El receptor del proyecto no es la causa.
2. **Probar un programa distinto** sobre las mismas cuatro poses: `smina.exe` ya
   está en `tools/`, y con anillos flexibles se puede probar `Vina-Carb` o
   aflojar la rigidez del anillo de la ribosa. La pregunta concreta es si algún
   protocolo razonable baja de 2 Å en las cuatro.
   > **HECHO EN PARTE (20 sep, ver `INFORME_PROTOCOLO_ALTERNATIVO_2026-09-20.md`)**.
   > `tools/smina.exe` está corrupto (9 bytes, «Not Found»), así que no se pudo
   > usar; `ad4` no es evaluable sin receptor con mapas AutoDock4. Con `vina` y
   > `vinardo` (las dos de Vina 1.2.3) y con el ligando **rígido** el control pasa
   > en isoproterenol (0,33–0,65 Å) y adrenalina (0,43–0,61 Å), y sigue fallando
   > en dopamina (2,97–3,02 Å) y 5-fluorouridina (9,8–10,0 Å). El siguiente
   > candidato es el receptor flexible.
3. **Solo si un protocolo pasa el control**, se rehace la validación de SOD1 con
   el conjunto limpio (`controles_sod1_v3.csv`) y los dos fondos.
4. **Avisar a Fredy**: el 0,815 y la lista de SOD1 quedan marcados como
   provisionales hasta que el protocolo se valide. Él pidió que ningún número se
   presente como bueno si no lo está.

---

## 7. Lo que este trabajo deja dicho, y sirve

- El control de redocking es **el test que el proyecto debía haberse hecho antes
  de gastar una hora de GPU**, y ahora está hecho, con números, con tres
  semillas, con dos tamaños de caja, con dos exhaustividades y con el análisis de
  los nueve modos.
- Se descartaron tres explicaciones fáciles (caja, búsqueda, preparación) con
  medidas, no con opiniones.
- Quedan dos herramientas reutilizables: `redock_trp32.py` (control de redocking
  para cualquier sistema nuevo antes de usarlo para cribar) y
  `comparar_sitios.py` (¿los ligandos de un conjunto ocupan de verdad el mismo
  sitio?).
- Y queda dicho, para que no se repita: **en estas estructuras hay que elegir la
  copia del ligando con criterio, no la primera que aparezca.**
