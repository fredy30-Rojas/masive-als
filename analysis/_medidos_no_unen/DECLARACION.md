# El tercer fondo: compuestos medidos que no unen (declarado el 25 sep 2026)

**Estado hoy: el fondo no existe.** Existen los compuestos que lo formarían —los que el
propio cribado propone medir, en `peticion_ensayos.csv`— y existen sus números, ya
acoplados con el mismo receptor, la misma caja y la misma exhaustividad que la validación
de su diana (`manifiesto.csv`, con la columna `estado` en `pendiente_medicion`). Lo que
falta es el resultado del laboratorio, y eso no se puede inventar.

Esta declaración se escribe **antes** de tener el primer dato, que es la única forma de
que valga: si la forma de leer el fondo se elige después de ver el número, lo que se mide
es la elección.

---

## 1. Por qué hace falta un tercer fondo

Los dos fondos que ya están tienen virtudes y un hueco:

- los **señuelos emparejados** (122 en TDP-43, 339 en SOD1) se parecen a los positivos en
  tamaño, logP, TPSA, donantes, aceptores, rotables y anillos, pero **de ninguno se sabe si
  une o no**;
- el **fondo duro** (152 unidores de ARN de R-BIND 2.0; 143 quelantes y redox en SOD1) son
  moléculas con unión medida, pero **a otra cosa**.

Ninguno de los dos es lo que un benchmark pediría de verdad: compuestos **medidos y
negativos en ese mismo sitio**. Ese hueco está escrito como limitación en el informe del
CR y en el paper, y este bloque es la forma de cerrarlo cuando el dato llegue.

## 2. Cómo se lee

1. **Es un bloque aparte.** Se puntúa contra los positivos igual que los otros dos, y
   **no se mezcla** con ninguno: no entra en `fondo` ni en `fondo2` ni en el bloque
   `juntos`. Si algún día se mezclara, los números de los bloques ya publicados cambiarían
   sin que nadie hubiera medido nada nuevo, que es exactamente lo que hay que evitar.
2. **Un compuesto entra con su ensayo y su condición.** No basta «no une»: se guarda el
   constructo, el método (SPR, ITC, MST, RMN, CETSA), el rango de concentración y la
   fecha. Un negativo es negativo **frente a ese constructo y ese método**, y así se
   escribe.
3. **El sitio manda.** Un compuesto que no une en otra caja no entra aquí: entra en el
   bloque de su caja, o en ninguno. Es la misma regla que se le aplica a los positivos.
4. **Piso declarado: 30 compuestos.** Con menos de 30 medidos, el bloque **se lee en
   descriptivo** —cuántos son, y en qué puesto queda cada uno— y **no** como fondo con el
   que decidir: el intervalo de un fondo tan corto es tan ancho que no decide nada. El
   piso es una convención declarada hoy, no un cálculo; y no se baja después para llegar a
   un veredicto.
5. **Las dos lecturas, y la que decide.** Contra este fondo se leen el AUC por átomo
   pesado (la que decide en todo el proyecto, porque quita el premio al tamaño) y el
   puesto de cada inactivo medido. El AUC **crudo no se usa para decidir**, como en los
   otros bloques: ordena por tamaño y está por debajo del azar.

## 3. Qué contesta y qué no

**Contesta** si el embudo coloca los positivos medidos por delante de compuestos que se
midieron y no unieron en ese sitio. Es la pregunta más exigente que se le puede hacer con
datos reales, más que contra señuelos parecidos.

**No contesta** —y conviene tenerlo escrito— si un candidato es bueno. No valida
candidatos, no cambia la decisión sobre los otros dos fondos y **no permite decir «estos
compuestos no se unen»** en general: permite decir que no unieron en ese ensayo, en ese
constructo y hasta esa concentración.

## 4. Cómo está montado, y la comprobación de que no estorba

- `analysis/fondo_inactivos.py` prepara y acopla los candidatos con el mismo receptor,
  caja y exhaustividad de la validación de su diana (Vina de CPU, `--cpu 1`), y escribe
  `manifiesto.csv`. Las poses quedan en `SOD1/out/` y `TDP-43/out/` con el prefijo `MED_`,
  para no confundirlas nunca con las `ACT_` de los positivos ni con las de los fondos.
- El bloque existe en la regla **solo cuando aparece una fila con papel `inactivo`** en el
  CSV de una corrida; mientras no haya ninguna, no sale en la tabla de la regla.
- Comprobado el mismo día: con el bloque añadido a `regla_decision_bootstrap.py`, la salida
  de hoy (`regla_decision/regla_decision.csv`) es **idéntica byte a byte** a la copia
  congelada en `regla_decision/_antes_2026-09-25/`. Añadir la posibilidad no cambia ningún
  veredicto que ya esté publicado.
- Comprobado también que el validador sigue dando **los mismos números**: la corrida de
  SOD1 del 25 sep escribe el mismo CSV que el publicado, en el mismo orden y con los mismos
  valores, **solo** con una columna nueva y vacía (`puesto_fondo3`, el puesto dentro de
  este bloque).
- Y probado el montaje entero con un manifiesto de prueba —clozapina marcada `no_une` en
  su misma caja—: entra con su afinidad (−5,935) y sus 23 átomos pesados; los 11 que siguen
  en `pendiente_medicion` **no entran** (por eso el cargador pide la lista de nombres, y no
  carga la carpeta entera); y el negativo de SOD1 no aparece en la corrida de TDP-43. La
  prueba no dejó nada en el repositorio.

## 5. Reproducción

```
python analysis/fondo_inactivos.py --solo-listar   # el manifiesto, sin acoplar nada
python analysis/fondo_inactivos.py                 # prepara y acopla lo que falte
python analysis/regla_decision_bootstrap.py        # hoy el bloque no sale: no hay inactivos
```
