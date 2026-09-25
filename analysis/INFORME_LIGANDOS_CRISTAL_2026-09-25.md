# Los positivos que faltan: barrido sistemático del PDB y de la literatura (25 sep 2026)

**De dónde viene esto.** La regla de decisión del 24 de septiembre
(`REGLA_DECISION_2026-09-24.md`) dejó los veredictos claros y dejó escrito lo único que
los puede cambiar: **faltan positivos de unión medida** — entre 3 y 6 en TDP-43 para el
bloque de señuelos emparejados, unos 11 en SOD1 para su fondo duro. Ni más exhaustividad
ni otra GPU mueven eso. Hoy se ha ido a buscar esos positivos por los dos caminos que no
gastan cómputo: **el PDB entero** (un ligando co-cristalizado es evidencia estructural de
unión, con cita y con quimiotipo propio) y **la literatura de ensayos de unión**.

**El resultado, en una línea:** el PDB está **agotado** para las dos dianas y para el
sitio que usamos —los 12 ligandos nuevos de SOD1 se asientan todos fuera del Trp32—, así
que los positivos que faltan solo pueden venir de ensayos en disolución; y de la
literatura sí han salido dos citas que faltaban, una de ellas un unido medido nuevo
(XL20, en el C-terminal, no en el bolsillo de ARN).

---

## 1. El barrido del PDB, y por qué es sistemático y no a ojo

`analysis/buscar_ligandos_cristal.py` (nuevo) hace tres cosas, todas de solo lectura:

1. pide a la RCSB las entradas del PDB por **accesión UniProt** (SOD1 `P00441`,
   TDP-43 `Q13148`, FUS `P35637`), no por búsqueda de texto: así entra todo y no lo que
   uno se acuerde de buscar;
2. saca de cada entrada sus **entidades no poliméricas** (GraphQL, por lotes);
3. **filtra aditivos** de cristalización (iones, sulfato, glicerol, PEG, DMSO, TFA,
   cisteína libre…) y pide el **SMILES depositado** de cada compuesto, agrupando por
   **esqueleto de Murcko**.

| Diana | Entradas en el PDB | Ligandos útiles | Ya en la verdad de referencia | Nuevos |
|---|---|---|---|---|
| SOD1 (P00441) | 156 | 22 | 10 | **12** |
| TDP-43 (Q13148) | 44 | **0** | 0 | 0 |
| FUS (P35637) | 23 | **0** | 0 | 0 |

**TDP-43 no tiene un solo ligando pequeño co-cristalizado en todo el PDB, y FUS tampoco.**
Todo lo que el proyecto usa como positivo de esas dos dianas viene de ensayos en
disolución (MST, RMN, SPR, celulares), nunca de una estructura con ligando. Eso no se
sabía con esta certeza: ahora sí, y está medido.

## 2. Y en SOD1, los nuevos no están donde los necesitamos

`analysis/sitio_ligandos_cristal.py` (nuevo) lee, **de las coordenadas depositadas**, en
qué sitio se asienta cada ligando: qué residuos de proteína quedan a menos de 4,0 Å, si el
**Trp32** está entre ellos y a qué distancia. No acopla nada y no usa GPU: es medición
directa.

**El control interno del método pasa:** los ocho ligandos que la verdad de referencia ya
tenía por Trp32 (las tres quinazolinas, isoproterenol, adrenalina, 5-fluorouridina, el
Lig9 y el 946) se leen **todos en Trp32**, a 3,3–3,9 Å. El método lee lo que hay.

Y los 12 nuevos:

| Compuesto | PDB | Quimiotipo | Sitio medido | Cita del depósito |
|---|---|---|---|---|
| 9JT (Ebselen) | 5O40, 6SPH, 6Z4G, 7T8F | benzoisoselenazolona | **Cys111** | Nat Commun 2018 (5O3Y, hermano) |
| 9JK (Ebsulfur) | 5O3Y | benzamida-S | **Cys111** | Nat Commun 2018 |
| Q6H, Q7T, Q7W, Q8H | 6Z3V, 6Z4O, 6Z4J, 6Z4L | benzoisoselenazolona | **Cys111** (familia) | EBioMedicine 2020 |
| LQN, LQW, LR2 | 6SPK, 6SPI, 6SPJ | benzoisoselenazolona | **Cys111** | Commun Biol 2020 |
| D7Z (PtCl₂-DACH) | 6FFK | complejo de platino | **Cys111** | ACS Med Chem Lett 2018 |
| ORA | 6B79 | colorante azo | residuos 3–7 de un **péptido**, no de la proteína | Protein Sci 2018 |
| I2G | 7XX3 | ácido graso | superficie (Lys9, Gly10, Gln15, Leu42…) | 7XX3 |

**Ninguno en Trp32.** Los once de la familia del ebselen y el de platino se unen a la
**Cys111**, que es un sitio real y muy poblado (es la cisteína reactiva: ebselen se une
covalentemente ahí), pero **no es la caja del proyecto**. El ORA no es un ligando de SOD1:
6B79 es un cristal del **péptido** 28-38 y el colorante se apoya en el péptido (residuos
3-7 de esa construcción), no en la proteína. El I2G es un ácido graso de superficie.

Consecuencias, sin adornos:

- **El PDB está agotado para nuestra pregunta.** No hay más positivos del bolsillo del
  Trp32 esperando a que alguien los busque: el espacio estructural de SOD1 tiene 22
  ligandos y todos los del Trp32 ya estaban.
- **El único sitio con muchos ligandos es Cys111**, y no sirve como sustituto: son
  compuestos reactivos (selenio, platino) que forman enlace covalente con la cisteína, y
  validar ahí sería validar otra cosa. Queda escrito por si algún día se quiere, pero no
  es lo que hace falta hoy.

## 3. Lo que sí se ha ganado: dos citas que faltaban

### 3.1 El 946 ya tiene su fuente y su afinidad medida

La verdad de referencia tenía el `naftalenoaminoalcohol_946` (PDB 5YTO) citado como
«PDB 5YTO (Wright/Hasnain)», sin artículo. La cita primaria del depósito es
**Manjula et al. 2018, FEBS Lett 592:1725**, y ese trabajo mide la unión por **SPR:
Kd 33,1 µM**, con la oxidación del Trp32 inhibida por el ligando. Es el mismo compuesto
(la estructura es la de ese trabajo). El 946 pasa de «co-cristal sin afinidad ni fuente»
a **positivo con unión medida y con cita**, que es lo que el criterio de entrada pide:
SOD1 sigue con 11 aptos, pero ahora con mejor respaldo.

### 3.2 XL20: el primer unido medido del C-terminal de TDP-43

Buscando positivos nuevos apareció **XL20** (Gao et al. 2026, *Nature Aging* 6:1667),
encontrado por cribado virtual sobre la **región conservada (CR) 320-340** del dominio de
baja complejidad. Lo que tiene medido, leído del artículo:

- **SPR sobre TDP-43 completa y sobre el péptido CR**: unión micromolar (el propio
  trabajo la llama **aparente**, no una Kd formal: desvía del modelo 1:1),
- **CETSA** (estabilización térmica en lisado celular),
- **dependencia de sitio medida**: borrando el CR la unión desaparece, y mutando el
  **Trp334** baja mucho; con eso el sitio no es una suposición.

Además, XL21 y XL23 inhiben la agregación in vitro pero **no tienen unión medida**: no
entran. XL20 entra en `verdad_de_referencia.csv` **con su cita, su tipo de ensayo y su
estructura, y marcado como no apto**: se une a otro sitio (el CR del C-terminal), no al
bolsillo de RRM1-RRM2 que acoplamos. Contarlo como positivo de RRM sería mentir sobre el
sitio.

### 3.3 Y su estructura, leída del dibujo

Su SMILES no está en el texto del artículo ni en PubChem: aparece **dibujada** en la
Figura Suplementaria 1(a). La tabla del panel (b) da su código de catálogo de Asinex,
**BDF34019555** (el texto del artículo lo confirma: «XL20 (BDF34019555) was identified
from the Asinex compound library and synthesized by WuXi AppTec»), y ese código tampoco
resuelve a una estructura en PubChem.

Así que se leyó el dibujo, con el camino que ya se usó para los tres fragmentos de
Nshogoza 2019:

- `analysis/leer_figura_xl20.py` (nuevo) corta la fila de ocho dibujos por los huecos de
tinta, asigna los bloques **en orden** a las etiquetas XL20–XL27 (comprobadas leyendo la
tira superior) y pasa cada celda por **DECIMER 2.7.2**, el OCR químico.
- **XL20 = `CN(Cc1ccccc1)[C@H]1CCC[C@H](n2cnc3c(N)ncnc32)[C@H]1O`**: `C19H24N6O`,
352,4 Da, logP 2,0, TPSA 93 Å², cuatro anillos. Es una **adenina unida a un
aminociclohexanol bencilado**, un quimiotipo nuevo y coherente con lo que el artículo
describe: un compuesto que entra al cerebro y que se apoya en el **Trp334** por
contactos hidrófobos y aromáticos.
- **Comprobaciones hechas:** válido en RDKit; **los dos recortes distintos dan el mismo
canónico** —y no solo en XL20: las ocho estructuras coinciden en las dos pasadas,
comparando el fragmento mayor, porque el recorte más ancho coge ruido de las etiquetas—;
y el código de Asinex leído por OCR coincide letra por letra con el del texto.
- **Segunda lectura, con visión.** `analysis/verificar_xl20_vision.py` (nuevo) le pasa el
recorte al modelo de visión local **`qwen3-vl:8b`** y le pide que describa lo que ve (a
propósito **no** se le pide el SMILES: un modelo de visión escribe SMILES mal y eso no se
nota). El modelo nombra **los cinco fragmentos** que dice el SMILES —benceno, ciclohexano
saturado, purina con amino, hidroxilo y un nitrógeno con metilo y bencilo—, y en la
pregunta de comparación lee lo mismo en el recorte original y en el dibujo de RDKit. La
traza entera queda en `_xl20_fig/verificacion_vision.txt`.
- **Trampa medida, y conviene no repetirla:** `qwen3-vl:8b` se pasa el presupuesto de
tokens **pensando** y devuelve el campo de la respuesta **vacío** (`done_reason: length`).
Con 500 y con 1.500 tokens la respuesta salió en blanco y el contenido entero estaba en
`message.thinking`. El script pide 4.000 y guarda las dos cosas; aun así, el formato de
lista fija de cinco líneas se lo salta y hay que buscar los fragmentos en el razonamiento.
- **Lo que sigue faltando:** el contraste con una base de datos. **PubChem no reconoce ni
el código ni la estructura**, así que la identidad descansa en el dibujo del artículo y en
estas dos lecturas, no en un registro. La procedencia completa está en
`analysis/controles_tdp43_xl20.csv`, y para poder mirarlo queda la imagen de comparación
—recorte original a la izquierda, dibujo de RDKit a la derecha— en
`_xl20_fig/XL20_comparacion.png`.

## 4. Qué significa para el plan

1. **Los positivos que faltan no van a salir del PDB.** Ya está barrido entero para las
   tres dianas. La única fuente es la literatura de ensayos de unión (SPR, ITC, MST, RMN,
   CETSA) y los repositorios de actividad —y en ChEMBL, para SOD1, lo que hay con
   `assay_type=B` **no es unión a proteína**: son IC50 de inhibición de la reducción de
   NADPH y EC50 celulares, que es justo lo que el criterio de entrada excluye.
2. **Los veredictos no cambian hoy.** TDP-43 sigue PASA solo contra el fondo duro; SOD1
   sigue sin evidencia; el bloque de señuelos emparejados de TDP-43 sigue necesitando
   10-13 positivos y hoy hay 7. No se ha añadido ningún positivo apto nuevo: se ha añadido
   una cita que faltaba y un positivo de otro sitio.
3. **La lectura honesta de §5:** el bloque que decide (fondo duro de TDP-43) ya pasa y no
   depende de cazar más positivos; el bloque blando no se va a cerrar por esta vía a corto
   plazo. Lo que queda por delante es la línea del C-terminal, y ahí XL20 es el primer
   positivo con unión medida que se puede nombrar.
4. **Lo que NO se hace:** no se toca la lista de candidatos, no se relanza el cribado y no
   se cuenta el 946 dos veces ni el ebselen como positivo del Trp32.

## 5. Cómo se reproduce

```
python C:/Users/Fredy/masive-als/analysis/buscar_ligandos_cristal.py   # barrido del PDB
python C:/Users/Fredy/masive-als/analysis/sitio_ligandos_cristal.py    # sitio real de cada ligando
python C:/Users/Fredy/masive-als/analysis/verdad_de_referencia.py      # regenera el CSV

# estructura de XL20 (figura suplementaria del articulo, OCR quimico)
ocsr_env/Scripts/python.exe analysis/leer_figura_xl20.py analysis/_xl20_fig/suppfig1_full.png

# segunda lectura con el modelo de vision local (hace falta `ollama serve`)
python analysis/verificar_xl20_vision.py
```

- Salidas del barrido: `analysis/ligandos_cristal/{sod1,tdp43,fus}.csv` y `resumen.txt`.
- Salidas del sitio: `analysis/ligandos_cristal/sitio_sod1.csv` y `sitio_sod1.txt`
  (contactos a 4,0 Å, distancia mínima al Trp32 y clasificación).
- Salida de la verificación con visión: `_xl20_fig/verificacion_vision.txt`, con las tres
  respuestas y sus razonamientos.
- Salidas de la lectura de XL20: `_xl20_fig/XL20.png` (la celda recortada),
  `XL20_comparacion.png` (recorte frente a dibujo de RDKit) y `smiles_crudos.csv` con las
  ocho lecturas. La figura sale del paquete suplementario del artículo, que se descarga
  entero con `https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13472882/supplementaryFiles`
  (la página de PMC sirve un aviso de descarga en vez del PDF, así que por ahí no baja).
- Descartes que quedan escritos: cuatro entradas de la familia del ebselen
  (6Z3V, 6Z4J, 6Z4L, 6Z4O) son solo mmCIF y no se leyeron sus coordenadas con el lector de
  PDB; su familia ya está medida en Cys111 por otros cinco miembros, así que no cambia la
  conclusión. El filtro de aditivos se corrigió tras la primera pasada (se colaban un PEG,
  un trietilenglicol, un trifluoroacetato, D-serina y cisteína libre), y la comparación con
  la verdad de referencia se hace **sin estereoquímica**, porque el descriptor del PDB no
  la lleva: sin ese arreglo el mismo compuesto se contaba como nuevo.
