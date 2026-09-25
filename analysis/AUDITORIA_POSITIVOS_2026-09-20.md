# Auditoría de los conjuntos de positivos — 20 de septiembre de 2026

**La pregunta, y es la que decide todo el proyecto:** los 20 activos de SOD1 y los
8 de TDP-43, ¿son compuestos independientes o son un mismo esqueleto repetido?

Importa por un motivo documentado: con activos analógicamente próximos entre sí y
señuelos fáciles, **cualquier** modelo saca enriquecimientos de fantasía. Es el
sesgo de análogos que la propia literatura del conjunto DUD-E documentó (Chen 2019,
*"Hidden bias in the DUD-E dataset leads to misleading performance"*). Y es
exactamente lo que acabamos de medir en el cribado focalizado de TDP-43.

Los scripts: `auditar_positivos_sod1.py` y `auditar_similitud_tdp43.py`. Ninguno
acopla nada: usan puntuaciones y SMILES que ya estaban calculados.

---

## 1. SOD1 — el pilar del proyecto no aguanta la prueba

SOD1 es la diana que el proyecto da por **validada**: AUC 0,815 frente a 199
señuelos emparejados en propiedades (bolsillo Trp32, `validacion_SOD1_trp32.csv`).
Nadie había mirado quiénes son esos 20 activos.

**Resultado:**

| Medida | Valor |
|---|---|
| Esqueletos de Murcko repetidos | **16 de 20 son el mismo núcleo pirazolona** (9 + 3 + 3 + 1) |
| Similitud máxima de cada activo con otro activo | mediana **0,71**, media 0,64 |
| Activos con otro activo a similitud ≥ 0,5 | **18 de 20** |
| Activos fuera de la serie | **2** (PRG-A01 y LCS-1) |

Es decir: no son veinte medidas independientes. Son **una serie congénérica** con
dos compuestos de otra química. Los identificadores van seguidos
(`CHEMBL2165601`…`CHEMBL2165614`, más `CHEMBL1643541/56/57`), lo que apunta a **un
solo trabajo de química medicinal** — la contigüidad y el esqueleto compartido lo
indican; el trabajo en concreto no se ha verificado todavía.

**Y ahora el AUC con y sin la serie:**

| Conjunto | AUC |
|---|---|
| Los 20 activos | 0,815 |
| Sólo la serie pirazolona (18) | **0,839** |
| **Sólo los de fuera de la serie** (2) | **0,601** |

El 0,815 lo sostiene la serie. Quitando la serie queda **0,601 con dos compuestos**,
y de esos dos uno va bien y el otro no:

```
PRG-A01    puesto   9/219   -6,196   mejor que el  4 % del fondo
LCS-1      puesto 152/219   -4,906   mejor que el 76 % del fondo
```

**LCS-1 es el más incómodo de los dos:** es el ligando de SOD1 con más respaldo
independiente y bien documentado del conjunto, y el acoplamiento **no lo encuentra**
—cae en el cuarto peor de la tabla. Una validación cuya señal viene de recuperar una
serie y que pierde al ligando de referencia real es **una validación específica de la
serie**, no una validación del método.

**Lo que esto no dice:** no dice que los pirazolonas no se unan a SOD1. Dice que 18
de los 20 positivos miden **una sola química**, y que el AUC se debe a recuperarla.

**Lo que sí dice:** antes de seguir usando SOD1 como el pilar «validado», hay que
(a) colapsar la serie a **un representante** para no contarla dieciocho veces, y
(b) buscar positivos de **otras quimias** con respaldo experimental.

---

## 2. TDP-43 — el mismo defecto, y peor

Detalle completo en `INFORME_RESCORING_FOCALIZADA_2026-09-19.md`, secciones 3 y 4.
Resumen:

- De los 8 controles, **sólo 2 tienen unión medida** (PE859 y berberrubina); los
  otros 6 entraron **por parecerse a la berberrubina**.
- La métrica de similitud da **AUC 1,000 dentro de su propia familia**, **0,538**
  para los que no son de ella (azar), y **una sola molécula** de referencia bate al
  conjunto entero (nitidina 0,952 contra 0,885 de los ocho).
- Una regla de dos líneas (anillos aromáticos + catión) da **0,921**, igual que la
  huella de 2048 bits.

---

## 3. La conclusión, dicha sin adornos

**Ninguna de las dos validaciones del proyecto mide generalidad química.** Las dos
miden, en el fondo, lo mismo: **la capacidad de recuperar la química con la que se
sembró**. Y eso —recuperar una serie congénérica cuando ya tienes un miembro de
ella— es una función legítima y útil, pero **no es descubrir**.

Consecuencia práctica, y es la que ordena el trabajo:

1. **Mientras la verdad de referencia esté sesgada, ninguna métrica es evaluable.**
   Ni el acoplamiento, ni el MM-GBSA, ni la similitud. Ni en SOD1 ni en TDP-43.
2. **El criterio de aceptación se fija por adelantado y es uno solo:** el ranking
   debe colocar a positivos de **al menos dos quimiotipos distintos** por delante de
   los señuelos. Con un solo quimiotipo de positivos, el resultado no se reporta.
3. **El fondo tiene que ser difícil** (señuelos emparejados + unidores de ARN de
   R-BIND para TDP-43), no 200 moléculas al azar.
4. **Nada de GPU** hasta que 1, 2 y 3 estén resueltos.

---

## 4. Lo que hay que hacer, y no necesita tarjeta

1. **Colapsar las series:** SOD1 pasa de 20 «activos» a ~3 independientes
   (1 pirazolona representativa + PRG-A01 + LCS-1). TDP-43 pasa de 8 a 7 de cinco
   quimias con los que ya están localizados (`PLAN_TDP43_2026-09-20.md`).
2. **Buscar positivos de otras quimias con cita y tipo de ensayo** para SOD1
   (candidatos conocidos fuera de la serie: LCS-1 y su serie propia, y lo que dé la
   literatura de ligandos de SOD1) y para TDP-43 (rTRD01, nTRD22, los tres
   fragmentos de RRM2, AIM4, bis-ANS).
3. **Re-testear las dos dianas con el criterio del punto 2 de arriba** y con los dos
   fondos. Vina CPU, horas — no GPU.
4. **Escribir el resultado**, sea el que sea. Si SOD1 y TDP-43 no pasan el criterio
   con positivos independientes, eso es un resultado publicable y honesto, y es más
   valioso que cualquier lista ordenada por un número que no significa nada.
