# Segundo control de redocking: otras funciones de puntuación y el ligando rígido

**20 de septiembre de 2026.** Script: `protocolo_alternativo.py` (también en esta
carpeta `analizar_rigido.py`). Registros: `protocolo_vina.log`,
`protocolo_vinardo.log`, `protocolo_rigido.log`, `score_alternativo.log`; datos en
`resultado_protocolo_alternativo_*.json`. Receptor y cajas idénticos a los del
primer control (`redock_trp32.py`): las cuatro estructuras cristalográficas del mutante
I113T (4A7S, 4A7T, 4A7U, 4A7V), caja de 24 Å centrada en el centroide del
ligando nativo, exhaustividad 8, semillas 42 / 2026 / 777.

Este informe responde al paso 2 de la sección 6 del documento anterior: **probar
otras funciones de puntuación y separar el problema de muestreo del problema de
puntuación.**

> **AVISO — LOS RMSD DE ESTE INFORME ESTÁN INFLADOS.** El medidor de RMSD
> emparejaba los átomos al revés; corregido y comprobado contra RDKit el mismo
> día (`validar_rmsd.py`). Las tablas corregidas están en
> `INFORME_CONTROLES_CORREGIDOS_2026-09-20.md`: el control **pasa en
> isoproterenol (0,33–0,79 Å) y adrenalina (0,43–0,93 Å)** con ligando flexible o
> rígido, y falla en dopamina (2,9–3,2 Å) y 5-fluorouridina (5,4–11,3 Å). Las
> cifras de energía (ΔE, deriva de la minimización) no cambian: lo que cambia es
> la lectura de ΔE, porque el «mejor modo» de las catecolaminas **es** la pose
> nativa (0,5 Å), así que ΔE mide la sensibilidad del potencial a medio
> ángstrom, no una preferencia por una pose equivocada. Los apartados 3 y 4
> (receptor flexible, y el caso de la 5-fluorouridina) se completaron en el
> informe nuevo.

---

## 0. Dos cosas que se cayeron antes de empezar

1. **`tools/smina.exe` está corrupto.** Ocupa 9 bytes y contiene la palabra
   «Not Found»: es una descarga fallida, no un binario. El informe anterior lo
   daba por disponible y no lo estaba. No se usó smina.
2. **La función `ad4` no se pudo evaluar con este receptor.** Vina 1.2.3 la
   incluye, pero responde
   `ERROR: No receptor allowed, only --flex argument with the AD4 scoring function.`
   Para AD4 hace falta un receptor tipado con mapas AutoDock4 (autogrid), que
   aquí no está preparado. Queda como pendiente explícito, no como resultado.
3. Lo que sí se pudo usar: **`vina` y `vinardo`**, ambas dentro de Vina 1.2.3,
   sobre las mismas estructuras y las mismas cajas.

---

## 1. Puntuar la pose del cristal sin moverla (`--score_only`)

Es la medida que el primer control no tenía. Si la pose cristalina puntuara
mejor que la que encuentra el acoplamiento, el fallo estaría en la búsqueda; si
puntúa peor, el fallo está en el paisaje de puntuación. **ΔE = E(cristal) −
E(mejor modo).** Positivo significa que la función prefiere una pose que no es la
del cristal.

### Acoplamiento del ligando flexible (la pose nativa como conformación de partida)

| Estructura | Ligando | Puntuación | E(cristal) | E(mejor modo) | **ΔE** | RMSD mejor / mediana | ¿Pasa? |
|---|---|---|---|---|---|---|---|
| 4A7T | isoproterenol | vina | −4,28 | −5,44 | **+1,16** | 1,88 / 1,91 Å | sí (medianas) |
| 4A7U | adrenalina | vina | −4,16 | −5,06 | **+0,90** | 3,74 / 3,74 Å | no |
| 4A7V | dopamina | vina | −3,91 | −5,15 | **+1,23** | 3,98 / 4,01 Å | no |
| 4A7S | 5-fluorouridina | vina | **+0,04** | −5,06 | **+5,10** | 8,59 / 9,59 Å | no |
| 4A7T | isoproterenol | vinardo | −3,59 | −5,30 | **+1,72** | 1,96 / 1,96 Å | sí (medianas) |
| 4A7U | adrenalina | vinardo | −3,73 | −4,74 | **+1,01** | 3,55 / 3,55 Å | no |
| 4A7V | dopamina | vinardo | −3,35 | −4,80 | **+1,46** | 3,92 / 3,96 Å | no |
| 4A7S | 5-fluorouridina | vinardo | +0,96 | −4,00 | **+4,95** | 4,21 / 4,25 Å | no |

**ΔE es positivo en los ocho casos.** Ninguna de las dos funciones prefiere la
pose del cristal. Y no es un detalle menor: las diferencias (0,9–1,7 kcal/mol,
salvo la 5-fluorouridina, que se va a 5) están **dentro del error habitual de
estas funciones**. La 5-fluorouridina es el caso extremo: la pose cristalina
puntúa **+0,04 kcal/mol**, es decir, la función no le ve unión ninguna.

## 2. Minimizar desde el cristal (`--local_only`): ¿es un mínimo local?

| Estructura | Puntuación | E(cristal quieta) | E(minimizada) | Deriva | Ligando **rígido**: deriva |
|---|---|---|---|---|---|
| 4A7T | vina | −4,28 | −4,80 | **1,84 Å** | **0,17 Å** |
| 4A7U | vina | −4,16 | −4,54 | **3,72 Å** | **0,20 Å** |
| 4A7V | vina | −3,91 | −4,41 | **3,42 Å** | **0,16 Å** |
| 4A7S | vina | +0,04 | −1,17 | **5,07 Å** | **2,92 Å** |
| 4A7T | vinardo | −3,59 | −4,61 | 1,83 Å | 0,28 Å |
| 4A7U | vinardo | −3,73 | −4,23 | 3,74 Å | 0,27 Å |
| 4A7V | vinardo | −3,35 | −3,91 | 3,42 Å | 0,25 Å |
| 4A7S | vinardo | +0,96 | −2,18 | 4,71 Å | 4,57 Å |

Aquí está el dato nuevo más útil: **con el ligando rígido, en las tres
catecolaminas la pose del cristal es un mínimo local estable** (deriva de 0,16 a
0,28 Å). Todo el desplazamiento de 3,4–3,7 Å que aparece con el ligando flexible
es **torsional**: la función gira los sustituyentes y saca el ligando del sitio,
ganando medio kcal/mol (0,38–0,52). Que 3,7 Å de desplazamiento cuesten 0,4
kcal/mol quiere decir que, en este bolsillo, **la energía no distingue**: el
paisaje es plano.

## 3. Acoplamiento del ligando rígido: muestreo contra puntuación

Se convirtió la pose cristalina en un ligando **rígido** (todos los átomos en el
ROOT, `TORSDOF 0`): el acoplamiento deja de muestrear torsiones y solo busca
traslación y rotación. Si así se reproduce la pose, el problema estaba en el
muestreo de torsiones; si tampoco, está en la puntuación.

| Estructura | Ligando | vina RMSD | vinardo RMSD | ¿Pasa? |
|---|---|---|---|---|
| 4A7T | isoproterenol | **0,63–0,65 Å** | **0,33 Å** | sí, las dos |
| 4A7U | adrenalina | **0,43–0,45 Å** | **0,61 Å** | sí, las dos |
| 4A7V | dopamina | 2,97 Å | 3,02 Å | no |
| 4A7S | 5-fluorouridina | 9,83 Å | 10,00 Å | no |

(Ojo con las energías rígidas: al ser `TORSDOF 0` la energía libre estimada pierde
el término torsional, así que los valores absolutos **no** son comparables con los
del acoplamiento flexible; solo lo es ΔE dentro de la misma corrida. En rígido,
ΔE sigue siendo positivo en los ocho casos: +0,81 y +1,38 en isoproterenol, +0,64
y +0,94 en adrenalina, +1,07 y +1,46 en dopamina, +5,84 y +5,26 en
5-fluorouridina.)

---

## 4. El diagnóstico, ahora sí partido en dos

**Lo que se arregla con el muestreo:** isoproterenol y adrenalina. Su pose
cristalina es un mínimo local estable del paisaje y el acoplamiento la encuentra
si se le da la conformación correcta y se le prohíbe girar torsionos. El fallo
flexible viene de que existen conformaciones vecinas que la función prefiere por
0,4–1,0 kcal/mol, con el ligando desplazado 3,7 Å. No es que el método no vea el
sitio: es que **no distingue la pose buena de la mala**.

**Lo que no se arregla con el muestreo:** dopamina y 5-fluorouridina.
- La dopamina falla incluso rígida (2,97 Å) frente a una pose que puntúa 1,07
  kcal/mol mejor. El paisaje prefiere otra colocación.
- En la 5-fluorouridina el paisaje **repele** la pose cristalina: sin ninguna
  torsión que girar, la minimización se la lleva a 2,9 Å (vina) y 4,6 Å (vinardo).
  Coincide con lo que dice el artículo: la ribosa queda proyectada al disolvente y
  con ocupación baja, y esa exposición es energéticamente desfavorable para la
  función de puntuación. **Es un caso donde la discrepancia puede ser del cristal,
  no del programa**, y así hay que reportarlo.

**Lo que ninguna de las dos funciones cambia:** vinardo no arregla nada respecto
a vina, ni en signo ni en magnitud (ΔE +0,9 a +1,7 en las catecolaminas, RMSD
iguales o peores). Cambiar la función de puntuación dentro de la misma familia de
potenciales no resuelve el problema.

---

## 5. Consecuencias para el proyecto

1. **El paso 1 de la sección 6 anterior ya estaba hecho y ahora está explícito**:
   el receptor de este control es el cristalográfico (4A7T y compañía), no 1HL5.
   No es el receptor del proyecto la causa.
2. **Ninguno de los protocolos probados pasa el control en las cuatro.** El mejor
   resultado es isoproterenol (1,88 Å de mediana con vina, 1,96 con vinardo), y
   aun así su pose nativa puntúa 1,16–1,72 kcal/mol peor que la que se encuentra.
3. **Ninguna lista ordenada con esta caja es defendible hoy.** El AUC 0,815 de
   SOD1 sigue **provisional**, ahora con un motivo más concreto: no solo el
   conjunto de positivos es una serie congénica, es que en el bolsillo de Trp32 la
   función de puntuación no distingue la pose cristalina de una a 4 Å.
4. La diferencia de 0,4–1,7 kcal/mol es del orden del ruido de la función. Por eso
   el punto 3 no se arregla con más muestreo ni con más semillas: con este
   potencial, el orden de la lista es ruido.

---

## 6. Qué sigue, en orden

1. **Receptor flexible.** Es el candidato inmediato y el único que puede mover la
   pose sin tocar la función de puntuación: `mk_prepare_receptor -f` con los
   residuos del bolsillo (Glu21, Gln22, Lys23, Lys30, Glu100) y el mismo control.
   Con isoproterenol y adrenalina el problema medido es de torsiones y de cadena
   lateral, y aquí es donde se puede ganar.
2. **Ensayo de ligando rígido + conformaciones del receptor (ensemble)**: usar los
   cuatro receptores (4A7S/4A7T/4A7U/4A7V) como un conjunto y quedarse con la
   mejor puntuación por ligando. Barato: los receptores ya están preparados.
3. **5-fluorouridina y dopamina no se rescatan con estos potenciales.** Para la
   ribosa haría falta anillo flexible (Vina-Carb) y para la dopamina un
   tratamiento del sitio que no existe en vina/vinardo. Se anotan como límite
   conocido del método, no como tarea pendiente indefinida.
4. **AD4 queda pendiente de verdad**: hace falta preparar el receptor con mapas
   AutoDock4 (`autogrid`/`prepare_receptor4`) si se quiere una tercera opinión
   con función de puntuación de otra familia. Sin eso, no se puede afirmar nada
   sobre AD4.
5. **Solo si alguno de los pasos 1 o 2 pasa el control**, se rehace la validación
   de SOD1 con `controles_sod1_v3.csv` y los dos fondos.
