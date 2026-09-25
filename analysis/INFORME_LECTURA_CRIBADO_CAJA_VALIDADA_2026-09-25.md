# La GPU ya cribó, pero en la caja que no está validada (25 sep 2026)

**Para qué es este informe.** Fredy preguntó para qué es todo lo que estamos haciendo y lo
dijo él mismo: *todo esto es para usar las supercomputadoras*. Este informe contesta a eso
con los números del disco: qué acopló ya la GPU, en qué caja, qué dice el propio embudo de
esa caja, y cuál es el trabajo de GPU que queda por hacer y para qué sirve.

---

## 1. Lo que hay ya hecho, contado

| Qué | Cuánto | Dónde | Fecha |
|---|---|---|---|
| Poses acopladas en el disco | **597.917** | `gpu_dock/resultados_libreria` (436.490) y `tanda_z001`–`z011` (~15.000 cada una) | ago 2026 |
| Pares ligando–caja con energía | **151.381** | `gpu_dock/resultados_vinagpu_total.csv` | 20 ago |
| — de ellos, por diana | FUS 51.432, SOD1 48.574, TDP43 51.375 | el mismo CSV | 20 ago |
| Librería de entrada | 66.478 moléculas | `analysis/full_library_solo.smi` | 20 ago |
| Ligandos ya preparados | **113.304** `.pdbqt` | `gpu_dock/libreria_ligands/` | ago 2026 |

Ese trabajo **no se tira**: es la búsqueda. Lo que hay que saber es con qué instrumento se
lee, y ahí es donde el trabajo de estos días cambia la conclusión.

## 2. La caja en la que se hizo el cribado NO es la que está validada

El cribado de las tandas z001–z011 usó estas cajas (`gpu_dock/run_tanda_z001.py`,
`RECEPTORES`):

| Diana | Receptor del cribado | Centro | Tamaño |
|---|---|---|---|
| TDP-43 | `gpu_dock/TDP43.pdbqt` | (28,3, 43,7, 52,5) | 25 Å |
| SOD1 | `gpu_dock/SOD1.pdbqt` | (27,9, 111,8, 64,4) | 25 Å |
| FUS | `gpu_dock/FUS.pdbqt` | (−14,5, 15,1, −7,8) | 25 Å |

Y estas son las que pasaron la regla de decisión del 24 de septiembre:

| Diana | Receptor validado | Centro | Tamaño | Veredicto |
|---|---|---|---|---|
| TDP-43 | `analysis/_tdp43_bolsillo_v2/4BS2_ph74.pdbqt` | (24,23, 16,89, −15,87) | 26 Å | **PASA** (0,731–0,738) |
| SOD1 | `gpu_dock/SOD1_limpio.pdbqt` | (46,5, 80,0, 73,3) | 22 Å | SIN EVIDENCIA (0,563) |

**No son las mismas cajas, y la diferencia está medida.** La caja vieja de TDP-43 se midió
con el embudo el 20 de agosto y dio un AUC de **0,465**, por debajo de la validación
(`analysis/REVISION_VALIDACION_BOLSILLOS_2026-08-20.md`: «TDP-43 | RRM1 (28,3,43,7,52,5) |
0,465 | ❌ NO VALIDADO»). Por eso se corrigió la caja ese mismo día. La caja corregida es
la que el 24 de septiembre pasó la regla, **y también con el motor de la GPU**
(`analysis/_control_gpu_TDP43/validar_gpu.csv`, search_depth 20: AUC 0,733, PASA).

Conclusión, sin adornos: **los 51.375 números de TDP-43 que hay en el disco no dicen nada
del bolsillo validado**, porque se midieron en una caja que el propio embudo mide como
azar. Usarlos para elegir candidatos es usar un instrumento que ya sabemos que a esa caja
no le acierta.

## 3. Lo único que se re-acopló en la caja corregida, y lo que enseñó

El 20 de agosto se re-acoplaron en la caja corregida **solo los de arriba del ranking
viejo** (`analysis/_redock_tdp43_masivo.py`: 2.776 hechos, con el informe de los 162
mejores). El resultado no fue neutro: **todos mejoran entre −0,4 y −1,3 kcal/mol** con la
caja corregida, y eso quiere decir que el ranking viejo **no ordena lo mismo** que la caja
buena. Da igual mirar la cola o la cabeza: cambiar de caja cambia el orden.

## 4. El trabajo de GPU que queda por hacer (y por qué es el bueno)

Cribar la librería entera en las **dos cajas corregidas**, con el **mismo protocolo del
control que la regla validó** (`analysis/control_gpu_tdp43.py`: `search_depth = 20`,
`num_modes = 3`, `thread = 8000`), porque es el único protocolo de GPU del que se sabe que
ordena química por encima del azar.

Está montado y comprobado en `analysis/lanzar_cribado_caja_corregida.py`:

```
python analysis/lanzar_cribado_caja_corregida.py --solo-listar   # no toca la GPU
python analysis/lanzar_cribado_caja_corregida.py --probar 500     # mide el ritmo real
python analysis/lanzar_cribado_caja_corregida.py                  # la librería entera
python analysis/lanzar_cribado_caja_corregida.py --diana TDP43    # solo una diana
```

Comprobado con `--solo-listar`: los dos receptores existen, las cajas son las de la
validación y hay **113.304 ligandos pendientes** por diana. Es resumible (salta lo que ya
tenga pose) y deja un CSV con el mismo formato que `resultados_vinagpu_total.csv`, para que
lo lea el mismo código. **Orden: TDP-43 primero**, que es la caja que pasa el bloque duro;
SOD1 después, que está en SIN EVIDENCIA y necesita ~11 positivos más para moverse.

La `--probar N` existe por una razón práctica: la última vez, la librería entera en una caja
ocupó del orden de días de 4080. Medir el ritmo con 500 ligandos antes de comprometer el
trabajo es una hora bien gastada.

## 5. Cómo se leerá lo que salga

Con la misma medida que pasó la regla, no con la energía cruda: **energía por átomo
pesado**, leída contra los fondos del proyecto (122 señuelos emparejados y 152 unidores de
ARN de R-BIND en TDP-43), con el control de dónde caen los positivos conocidos. El motivo es
el que ya está medido: en las 16 combinaciones probadas, la energía cruda ordenó **por
debajo del azar** (resumen de la regla del 24). Elegir candidatos por la energía más
negativa, que es lo que hace todo el mundo, aquí es tirar el trabajo de la tarjeta.

## 6. Lo que NO se hace

No se leen los 51.375 números viejos como si fueran de la caja validada. No se mezclan
motores ni protocolos sin decirlo. No se gastan horas de GPU en FUS: no tiene ni un positivo
medido, así que no hay con qué validar nada de lo que salga de ahí. Y no se lanza una caja
nueva por lanzarla: cada hora de GPU se justifica con la pregunta *¿contra qué positivos se
lee esto?*
