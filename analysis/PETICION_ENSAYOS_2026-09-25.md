# Petición de ensayos de unión medida (25 sep 2026)

**Para qué es este documento.** El 24 de septiembre la regla de decisión del proyecto
señaló el cuello de botella con un número delante: no falta cómputo, **faltan positivos de
unión medida**. El barrido sistemático del PDB del 25 de septiembre lo cerró por el otro
lado: por estructura no van a salir —22 ligandos co-cristalizados de SOD1 y **cero de
TDP-43 y cero de FUS**—. Esto es lo que se pide, a quién y para qué, con los compuestos
concretos y el ensayo de cada uno. La tabla que se puede mandar tal cual sale de
`analysis/peticion_ensayos.py`.

---

## 1. Qué cuenta como positivo y qué no

Un compuesto entra en `verdad_de_referencia.csv` solo si hay **unión directa medida**
(Kd, ΔG o CSP de RMN) **y el sitio está identificado o acotado**. No entran, aunque sean
útiles para otras cosas: la actividad celular (EC50, viabilidad, muerte neuronal), la
inhibición de agregación sin dato de unión, y la similitud con otro ligando.

El motivo está medido, no supuesto: en SOD1, los registros de ChEMBL con `assay_type=B`
—los que a primera vista parecerían de unión— son **IC50 de la inhibición de la reducción
de NADPH y EC50 celulares**, no unión a proteína. Y el sitio importa porque los bloques
de la regla se definen por caja: un unido del CTD no es un positivo del bolsillo de RRM1,
y hay siete positivos de TDP-43 que hoy están repartidos entre RRM1 (rTRD01), NTD
(nTRD22), RRM2 (tres fragmentos de Nshogoza) y la interfaz RRM1-RRM2 (PE859 y
berberrubina). **Un positivo de RRM1 vale más aquí que tres de otros sitios.**

## 2. Cuánto falta, con los números de la regla

| Diana | Bloque | Positivos hoy | Para que el límite inferior toque 0,5 | Veredicto |
|---|---|---|---|---|
| TDP-43 | fondo duro (152 unidores de ARN) | 7 | **ya basta con los que hay** | PASA |
| TDP-43 | señuelos emparejados (122) | 7 | ~10–13 en total → **3–6 más** | SIN EVIDENCIA |
| SOD1 | fondo duro (143 quelantes y redox) | 11 | ~22 en total → **~11 más** | SIN EVIDENCIA |
| SOD1 | señuelos emparejados (339) | 11 | ~444 | SIN EVIDENCIA |

Dos lecturas importantes de esa tabla. La primera: el bloque que decide (el duro) **ya
pasa en TDP-43** y lo único que le falta es el bloque blando, que necesita pocos positivos
pero del mismo sitio. La segunda: el bloque de señuelos emparejados de SOD1 pide 444
positivos, o sea que **no es el objetivo**; el que se puede llenar y sí decide es su fondo
duro, y ese pide del orden de once.

## 3. Los compuestos, por prioridad

### Prioridad 1 — cada uno mueve un número ya publicado (5 compuestos)

| Compuesto | Diana y sitio | De dónde sale | Qué se pide | Por qué |
|---|---|---|---|---|
| **XL23** (Asinex LAS51502065) | TDP-43, CR del C-terminal (320–340) | Gao et al. 2026, *Nature Aging* 6:1667 | SPR aparente con el mismo protocolo que dio micromolar para XL20, más CETSA en células; si hay constructo del CR, repetir la dependencia de sitio (sin CR no une, Trp334 mutado baja) | Es el hermano de quimiotipo de XL20 (19 átomos pesados comunes, misma adenina-aminociclohexanol) y tiene actividad funcional sin unión medida. Sin este dato la pareja es una cabeza común con dos colas |
| **LCS-1** (CAS 41931-13-9) | SOD1, sitio por determinar | Wright et al. 2013 / cribado previo | Kd por SPR o ITC **y** determinación de sitio (soaking o RMN) | Es un unido conocido de SOD1 por actividad celular y sin unión a proteína medida. Hoy no vale como positivo justo por eso |
| **PRG-A01** | SOD1, sitio por determinar | Woo et al. 2021 | Kd por SPR o ITC **y** sitio | Igual que LCS-1: es una de las dos únicas quimias independientes que SOD1 tiene y no está medida contra la proteína |
| **CHEMBL2165613** (pirazolona) | SOD1, sitio por determinar | serie ChEMBL2165601–2165614 | Kd por SPR o ITC | Es el representante **único** de 18 análogos colapsados: una medida cuantifica la serie entera, que es la que sostiene los controles de SOD1 |
| **clozapina** (CHEMBL42) | SOD1, Trp32 (caja del proyecto) | candidatos_total_cns.csv (−7,5 kcal/mol, CNS-MPO 4,77) | SPR de cribado contra SOD1 humana; Kd completa si hay señal | Fármaco aprobado, química independiente de todo lo que hay en la referencia (dibenzodiazepina frente a catecolaminas, quinazolinas y nucleósidos) y candidato de reposicionamiento con dato |

Dos detalles de disponibilidad: **LCS-1 se compra por catálogo** (Sigma 567417, Cayman
35231, Selleckchem, MedChemExpress) y **XL23 lleva código de catálogo de Asinex**
(LAS51502065), que es la misma casa de la librería de la que salió XL20. **PRG-A01 no
consta en catálogo comercial**: hay que pedirlo a los autores del trabajo original.

### Prioridad 2 — panel propio: dice si el embudo acierta (8 compuestos)

Los mejores puntuados del propio cribado, ninguno parecido a los controles. Cada uno que
una es un positivo nuevo; cada uno que no una es un error del ranking que también se
apunta. Y hay una ganancia que conviene decir: hoy **ningún señuelo del proyecto tiene
inactividad medida** —son parecidos en propiedades, no inactivos conocidos—, así que un
compuesto que se mida y no una sería el primer señuelo de verdad del proyecto.

| Diana y sitio | Compuestos | Ensayo |
|---|---|---|
| SOD1, Trp32 | CHEMBL3311449 (el mejor MM-GBSA, −33,4 kcal/mol), CHEMBL4553700, CHEMBL584356, CHEMBL8550 | SPR de cribado (una concentración) y Kd si hay señal |
| TDP-43, RRM1 | CHEMBL1362588, CHEMBL27093, CHEMBL19347, CHEMBL23927 (−7,3 a −7,2 kcal/mol por Vina) | SPR, o **RMN de 15N-HSQC con RRM1-RRM2**, que además localiza el sitio por el desplazamiento químico |

Para TDP-43 conviene el HSQC antes que el SPR: los positivos que faltan tienen que ser
del mismo bolsillo, y el CSP lo dice sin ambigüedad.

### Prioridad 3 — informativos, no cierran ningún bloque (2 compuestos)

**bis-ANS** y **Congo Red** tienen evidencia funcional de separación de fases y agregación
y ninguna unión medida. Medirlos daría positivos del CTD, no del RRM, y bis-ANS es una
sonda promiscua que se pega a casi cualquier cara hidrófoba expuesta: si se mide, hay que
distinguir unión específica de pegado inespecífico.

## 4. La tentación que ya está descartada: cambiar a Cys111 no vale

El barrido del PDB encontró para SOD1 doce ligandos co-cristalizados que no estaban en la
verdad de referencia: once de la familia del ebselen (benzoisoselenazolonas y una benzamida con selenio)
y uno de platino (PtCl2-DACH). Los doce están en **Cys111**, repartidos en **9 esqueletos
de Murcko**, y cada uno es evidencia de unión depositada en el PDB.

La tentación es evidente: mover la caja de SOD1 a Cys111 y tener un bloque con evidencia
estructural sin gastar ni un ensayo. **Está descartado, y conviene que quede escrito aquí
para no reproponerlo.** Son compuestos reactivos —selenio y platino— que forman **enlace
covalente** con la cisteína: validar ahí sería validar otra cosa, no la cara
antiagregación que el proyecto eligió a propósito. Y el PDB está agotado para la pregunta
que sí importa: los 22 ligandos depositados de SOD1 están contados y todos los del Trp32
ya estaban. Por eso lo que se pide es medir, y se pide en Trp32.

## 5. Qué haremos con cada resultado

Si **XL23** une: la pareja XL20–XL23 pasa a ser una relación estructura-actividad de dos
puntos y se puede montar el control que hoy no existe —¿reproduce el acoplamiento la
dependencia del Trp334?—; con eso el CR deja de ser una línea cerrada.

Si alguno del **panel** une: se añade como positivo independiente y se vuelve a correr la
regla con él dentro (`regla_decision_bootstrap.py`), que es lo único que cierra intervalos.
Si no une: queda como señuelo de inactividad medida, que es información que hoy no hay.

Si **LCS-1, PRG-A01 o la pirazolona** se miden y resultan débiles o negativas, eso también
decide, y en la dirección incómoda: quedaría dicho que parte de la referencia de SOD1 se
apoyaba en actividad celular heredada y no en unión a proteína, y que al embudo no se le
puede reprochar no recuperar lo que quizá no se une.

## 6. Lo que NO se pide

No se pide ensayar candidatos como apuesta terapéutica: lo que se pide es **calibrar el
embudo**, y por eso la prioridad 1 son compuestos que ya están en la literatura o en el
congelador, no los mejores de la lista. Un resultado positivo **no valida el embudo**: la
regla mide si ordena química por encima del azar en un bloque, nada más. Y no se pide
trabajo sobre FUS: en todo el PDB no hay un solo ligando co-cristalizado suyo y sus dos
entradas de la referencia no tienen ni tipo de ensayo.

Queda además un pendiente de escritorio, no de laboratorio: los **tres fragmentos de
Nshogoza (2019)** están en la verdad de referencia **sin SMILES** ("hay que dibujarlo"), así
que hoy no se pueden acoplar. Dibujarlos del artículo los mete en el juego sin coste de
ensayo.

## 7. Cómo se reproduce

```
python C:/Users/Fredy/masive-als/analysis/peticion_ensayos.py
```

- Lee `verdad_de_referencia.csv`, `controles_sod1_v4.csv`, `suplementario_tdp43_xl21_27.csv`,
  `candidatos_filtrados.csv`, `candidatos_total_cns.csv` y `regla_decision/regla_decision.csv`.
- Escribe `analysis/peticion_ensayos.csv` (una fila por compuesto, con SMILES, ensayo, motivo
  y qué cambia) y `analysis/peticion_ensayos.txt` (el resumen, con los positivos que faltan
  leídos de la regla y no de la memoria).
- Los SMILES y los identificadores vienen de esas tablas: no se escriben a mano aquí.
