# TDP-43 contra el bolsillo RRM2: con los fragmentos dentro, el ranking se hunde

**23 de septiembre de 2026.** Primera corrida de `validar_diana_limpia.py --target TDP43`.
Es la validación que quedó pendiente en `INFORME_RESCORING_ESTRATOS_2026-09-22.md` (§6, paso 3):
diana limpia, positivos de **unión medida** con su cita, fondo emparejado en propiedades y
lectura **dentro de estrato de tamaño**, no por el puesto absoluto.

Se corrió **dos veces**, y la segunda es la que cuenta.

**Scripts:** `analysis/validar_diana_limpia.py` (elige diana) → `analysis/validar_sod1_limpia.py`
(hace el trabajo) · **Datos:** `validar_tdp43_limpia.log`, `validar_tdp43_limpia.csv`,
`validar_tdp43_limpia_resumen.csv`.

---

## 0. Lo que hubo que arreglar antes de que el número significara algo

1. **El envoltorio no cambiaba la diana.** `validar_diana_limpia.py` sobrescribía
   `RECEPTOR`, `CENTRO`, `TAMANO`, `EXHAUSTIVIDAD`, `POSES_FONDO`, `LIGS_FONDO` y `SALIDA`,
   pero `validar_sod1_limpia.py` filtraba los positivos con **`r["target"] == "SOD1"` clavado
   a mano** y titulaba la corrida «VALIDACION DE SOD1». La primera tentativa llegó a acoplar
   **los 11 positivos de SOD1 dentro de la caja de TDP-43**; se vio a los pocos segundos por el
   encabezado del log y se mató el proceso (ninguna vina quedó viva). `DIANA` y `NOTA_FONDO`
   son ahora constantes de módulo que el envoltorio sí puede cambiar.
2. **Dos frases del log estaban escritas a mano y ya eran falsas:** decía «el AUC crudo es X,
   POR DEBAJO DEL AZAR» pasara lo que pasara, y la correlación del fondo se imprimía como un
   `-0.0` literal. Las dos se calculan ahora del dato.
3. **Los tres fragmentos de RRM2 no tenían SMILES** (`smiles` vacío, nota «hay que dibujarlo»).
   Se resolvió leyendo la Figura 1c — ver §3.

---

## 1. Primera corrida: 4 positivos, sin los fragmentos

| Métrica | Valor |
|---|---|
| AUC crudo | 0,589 |
| AUC por átomo pesado | 0,494 |
| AUC residual | 0,584 |
| EF1 % / EF5 % | 0,00 / 0,00 |
| Quimiotipos que baten a los señuelos de su mismo tamaño | **2 de 5** |
| Veredicto | PASA (por el mínimo) |

Positivos, de mejor a peor: nTRD22 (19… puesto 24), berberrubine (46), PE859 (55),
rTRD01 (85). Los dos que ganan a los señuelos de su tamaño: isoxazol-piperidina (nTRD22) y
bencilisoquinolina (berberrubine).

## 2. Segunda corrida: los 3 fragmentos dentro — y el ranking se hunde

| Métrica | Sin fragmentos | **Con fragmentos** |
|---|---|---|
| AUC crudo | 0,589 | **0,313** |
| AUC por átomo pesado | 0,494 | **0,693** |
| AUC residual | 0,584 | **0,313** |
| Quimiotipos que pasan por tamaño | 2 de 5 | **2 de 5** |

**Los tres fragmentos —los únicos compuestos de la tabla que tienen unión medida contra ESE
bolsillo— quedan los últimos de la lista:**

| Ligando | Quimiotipo | Puesto | Afinidad | AUC crudo de su quimiotipo |
|---|---|---|---|---|
| fragmento_1 | fragmento | **125/129** | −5,55 | **0,008** |
| fragmento_3 | fragmento | **127/129** | −5,37 | 0,008 |
| fragmento_2 | fragmento | **129/129** | −4,72 | 0,008 |

El resto de positivos: nTRD22 19/129, berberrubine 46/129, PE859 55/129, rTRD01 114/129.
**AUC crudo por debajo del azar (0,313).** El AUC por átomo pesado sube a 0,693 porque
descuenta el tamaño, que es justo lo que castiga a los fragmentos.

**Y el detalle que hay que decir: los dos únicos que pasan el criterio no se unen a ese
bolsillo.** nTRD22 es un modulador alostérico del **dominio N-terminal** (Mollasalehi 2020) y
la berberrubina se describe en la **interfaz RRM1–RRM2** (Kapsiani 2026). Ninguno de los dos
se une a G245/E246/H256/I257/S258. El «PASA» sale de compuestos que la propia literatura
coloca en otro sitio.

---

## 3. Cómo se resolvió el SMILES de los tres fragmentos

En la literatura, los tres hits de Nshogoza 2019 aparecen **solo como dibujo**, en la Figura
1c, etiquetados «hit 1», «hit 2» y «hit 3». Se revisó el artículo, su material suplementario
(solo tablas de anchura de línea), la revisión de François-Moutal 2021 (PMC8341936), que los
reproduce como imagen, y el disco: no hay nombre, número CAS, código de catálogo ni fórmula.

Se resolvió por **OCR químico** (DECIMER 2.7.2, Steinbeck lab) sobre la Figura 1c, con
`analysis/leer_figura_nshogoza.py`. Dos cosas hubo que arreglar:

- el recorte por bloques de columnas metía también **la etiqueta de texto** («hit 1»…), y
  DECIMER la leía como átomos: la primera pasada devolvió un complejo de paladio, siete boros
  y sodio que eran literalmente las letras. Se aísla el dibujo por su tramo de filas alto
  (~400 px) y se descarta el tramo bajo (~75 px) del texto;
- con el recorte limpio, las tres lecturas convergen.

| | SMILES leído | Fórmula | Masa |
|---|---|---|---|
| hit 1 | `C1CC(CN(C1)C2=NC=CC=N2)N` | C9H14N4 | 178,2 |
| hit 2 | `CC(C1=CN2C(=N1)SC=N2)N` | C6H8N4S | 168,2 |
| hit 3 | `C1=CN=C(C(=C1)C#N)N2CCC(CC2)N` | C11H14N4 | 202,3 |

**Comprobaciones hechas:** dos lecturas con recortes distintos dan el **mismo SMILES
canónico**; las tres masas caen en el rango de fragmento que declara el artículo
(110–250 Da); un control con una molécula dibujada por RDKit se lee sin un solo error; y la
descripción de las tres estructuras por el modelo de visión (qwen3-vl:8b) coincide con la
lectura del OCR — piperidina con NH₂ unida a pirimidina (1), biciclo con S y NH₂ (2),
piridina con nitrilo unida a piperidina con NH₂ (3).

**Lo que sigue pendiente, y hay que decirlo:** ese SMILES **no está publicado**. Es una
lectura automática, confirmada por dos vías, pero **nadie lo ha visto con ojos de químico**.
Queda anotado así en `verdad_de_referencia.csv`, con copia de seguridad
`verdad_de_referencia.csv.bak_pre_fragmentos_20260923`.

---

## 4. Lo que NO dice, en el mismo párrafo

1. **No se publica nada con esto.** Siete positivos deciden si una idea se sigue o se tira,
   no dan un número que reportar.
2. **El fondo duro sigue faltando:** los 122 señuelos son emparejados en propiedades; **no**
   está el fondo de verdaderos unidores de ARN de R-BIND 2.0. Sin él, el criterio del proyecto
   dice que el resultado **no se reporta como validación**. Esto es una decisión interna.
3. **El corte fijo no sirve en absoluto:** EF1 % y EF5 % son cero en las dos corridas.
4. **Es un MM-GBSA/Vina de una pose**, no dinámica.
5. **Los fragmentos son diminutos** (11–15 átomos pesados): su energía es pequeña por
   construcción, así que su último puesto no sorprende. Lo que **no** se puede decir es que
   esto demuestre que no se unan — solo que **la energía acoplada no los reconoce**.
6. **Ninguna métrica cruda ordena**, ni con 4 ni con 7 positivos.

---

## 5. Qué significa para el proyecto

- **El bolsillo de RRM2 parecía la apuesta buena y no lo es.** La comparación con el sitio de
  ARN (0,465 y 0,517) mejoraba con 4 positivos (0,589), pero al meter los tres compuestos que
  definieron ese bolsillo por RMN el AUC cae a 0,313 y ellos quedan en los tres últimos
  puestos. Era el mejor conjunto de positivos posible para esa caja —unión medida, tres
  quimias nuevas— y es el que peor sale.
- **El «PASA» por emparejado de tamaño no se sostiene:** quien lo hace pasar son dos ligandos
  que, según la literatura, se unen en otro sitio.
- **Lo que toca ahora, en orden:**
  1. **Conseguir el fondo duro de R-BIND 2.0.** Sin él nada de esto se reporta.
  2. **Comprobar los tres SMILES a ojo de químico** antes de que entren en cualquier cuenta
     que salga de aquí. El material para mirarlo está en
     `analysis/_nshogoza_fig1c/hit_*_comparacion.png` (original al lado del redibujado).
  3. **Decidir con esto delante** si TDP-43 sigue, y con qué sitio: el de ARN no ordena, el
     de RRM2 tampoco, y los dos que puntúan bien no se unen donde dicen. La alternativa que
     queda escrita en `PLAN_TDP43_2026-09-20.md` es la **agregación del C-terminal**, con la
     plantilla de la cremallera estérica de tau.
- **Y no toca nada más:** no se relanza el cribado, no se pide supercomputadora, no se
  presenta lista de candidatos.

---

## 6. Anexo: cómo se reprodujo

```
python C:/Users/Fredy/masive-als/analysis/validar_diana_limpia.py --target TDP43
```

Lanzador: `analysis/lanzar_validar_tdp43_limpia.bat`. El OCR de la figura:
`masive-als/ocsr_env/Scripts/python.exe analysis/leer_figura_nshogoza.py _nshogoza_fig1.png`,
y la comparación visual con `analysis/comparar_ocr_figura.py`.
