# El MM-GBSA de FUS terminó, el vigilante murió, y el consenso se quedó con un cero falso

**4 de octubre de 2026, tarde.** Tres cosas que estaban mal y ya no lo están.

---

## 1. El MM-GBSA de FUS está TERMINADO (y nadie lo había recogido)

En Oracle, `runner_mmgbsa_fus.py` acabó **119 de 119** candidatos a la **01:33** del 4 de
octubre (`[01:33:04] FIN FUS` en `runner_mmgbsa_fus.log`).

| Concepto | Valor |
|---|---|
| Candidatos procesados | **119 / 119** |
| Con energía calculada | **103** |
| Timeouts (>1800 s) | 16 |
| Mejor dG | **−25,62** (CHEMBL3311273) |

El fichero remoto es `~/mmgbsa/rescoring_fus_mmgbsa.csv`, 120 líneas.

**Lo que falló:** el vigilante `vigilar_mmgbsa_fus.py` llevaba desde las 02:43 sin poder
hablar con Oracle (`No route to host`) y se detuvo solo, con este mensaje:

```
[2026-10-04 02:43:05] el proceso no corre pero el CSV tiene 0/119 filas.
                      Puede haber muerto; no toco nada.
```

El mensaje es correcto como precaución, pero el «0/119» es el del **parcial local**, que
no se había podido refrescar. El trabajo en Oracle sí había terminado. Consecuencia: el
paso final del vigilante —recalcular `ranking_final.py` en Oracle y regenerar el cruce de
tres motores— **nunca se ejecutó**, y no dejó marcador (`mmgbsa_fus.done` no existía).

**Hecho a mano ahora:** recalculado el ranking en Oracle, traído el CSV (310 candidatos,
antes 202) y regenerado el consenso. Marcador `estado_400/mmgbsa_fus.done` escrito a mano
para que el vigilante no vuelva a intentarlo sobre trabajo ya cerrado.

FUS pasa de **9 a 112** candidatos en el ranking. Era la proteína abandonada del proyecto.

---

## 2. El consenso de tres motores llevaba desde las 02:00 con 0 ligandos en los tres

El fichero `resultados/consenso_tres_motores.csv` decía, desde las 02:00 del 4 de octubre:

```
candidatos: 202
con los tres motores medidos: 0
buenos en los tres (tercio superior): 0
```

**Ese cero no era un resultado: era un fallo de emparejamiento.** El motivo es una
discrepancia de nombres de una sola cifra:

- `gnina_repuntar.py` vuelca sus **10.000 filas por diana DENTRO de `gpu_dock/gnina_top5000/`**
  (no crea un `gnina_top10000/`). El propio `vigilar_cnn_10k.py` lo sabe: su `SALIDA_DIR`
  es `gnina_top5000`.
- Pero al regenerar el consenso, el vigilante le pasaba **`--gnina-top 10000`**.
- `consenso_tres_motores.py` buscaba `gpu_dock/gnina_top10000`, **no existía, devolvía
  vacío sin decir nada**, y escribía un CSV con la columna de CNN toda vacía. Un cero
  silencioso que se lee como «ningún candidato es bueno en los tres» cuando lo que pasa
  es que a nadie se le preguntó.

El marcador `estado_400/cnn_10k.done` (02:00) selló ese estado como bueno, así que el
vigilante no volvió a intentarlo.

---

## 3. Arreglado, con guardas, y probado

### `consenso_tres_motores.py`

- `leer_gnina()` ya no devuelve vacío si el directorio no existe: lanza `FuenteAusente`
  y el script **sale con código 2 sin escribir nada**. El mensaje lista los directorios
  de CNN que sí existen, para que el número correcto salte a la vista.
- Si el directorio existe pero no tiene filas utilizables, también se niega (código 2).
- Si **ningún** candidato acaba con los tres motores medidos, se niega y **no escribe**
  (código 3). Un cruce vacío no es un resultado.
- Las comprobaciones van **antes** de abrir el fichero de salida. En la primera versión del
  arreglo iban después: el aviso salía pero el CSV malo ya estaba escrito. Se detectó con
  una prueba de md5 sobre el fichero real y se movió el orden.

### `gpu_dock/vigilar_cnn_10k.py`

- Nueva constante `TOP_DIRECTORIO = 5000`, con el comentario de por qué no es el número de
  filas. `regenerar()` pasa **ese** número al consenso, no el de filas pedidas.
- La prueba comprueba por `assert` que `regenerar()` usa `TOP_DIRECTORIO` y que ya no usa
  el número de filas.

### Pruebas ejecutadas

| Prueba | Resultado |
|---|---|
| Directorio de CNN inexistente (`--gnina-top 10000`) | código **2**, no escribe |
| Directorio con CSV pero sin ningún ligando en común | código **3**, no escribe |
| Directorio con CSV sin filas utilizables | código **2**, no escribe |
| Directorio real (`--gnina-top 5000`) | código **0**, 310 candidatos, 204 con tres motores, 54 buenos en los tres |
| md5 del CSV antes y después de una corrida fallida | **idéntico** (no tocó el fichero) |
| `ast.parse` de los dos scripts | OK |
| `assert` de que `regenerar()` usa `TOP_DIRECTORIO` | OK |

---

## El resultado que ahora sí es cierto

`resultados/consenso_tres_motores.csv`, con `--gnina-top 5000` y el ranking nuevo:

| Concepto | Antes (falso) | Ahora |
|---|---|---|
| Candidatos | 202 | **310** |
| Con los tres motores medidos | 0 | **204** |
| Buenos en los tres (tercio superior) | 0 | **54** |
| Reparto | — | TDP43 99, SOD1 99, FUS 112 |

Los 12 candidatos con dG físicamente imposible (|dG| > 100 kcal/mol, choques de pose) siguen
marcados `sospechoso=SI` y **ninguno** entra en la lista de buenos en los tres. La criba
funciona.

Los dos ficheros anteriores quedan guardados como `resultados/*.bak_20261004` para poder
comparar.

---

## Lo que este episodio deja claro

Dos vigilancias automáticas fallaron la misma noche, y las dos fallaron **en silencio**:

1. Una murió sin marcador y su último paso no se ejecutó, pero nadie se enteró porque el
   fichero de salida seguía teniendo un contenido plausible (el viejo).
2. La otra escribió un cero que parecía un resultado.

Un vigilante que no distingue «no hay señal» de «no pregunté» es peor que no tenerlo,
porque sella el fallo con un marcador. De ahí las guardas: **una fuente que falta es un
error, no un dato vacío.**

---

## Ficheros

- `resultados/ranking_final_masive_als.csv` — 310 candidatos (antes 202), FUS 9 → 112
- `resultados/consenso_tres_motores.csv` — 310 candidatos, 204 con tres motores, 54 buenos
- `resultados/ranking_final_masive_als.csv.bak_20261004` — el anterior
- `resultados/consenso_tres_motores.csv.bak_20261004` — el anterior (el del cero falso)
- `consenso_tres_motores.py` — guardas de fuente ausente
- `gpu_dock/vigilar_cnn_10k.py` — `TOP_DIRECTORIO` y comentario del fallo
- `gpu_dock/estado_400/mmgbsa_fus.done` — marcador escrito a mano

## Reproducir

```
cd C:/Users/Fredy/masive-als
python consenso_tres_motores.py --gnina-top 5000
```

Con `--gnina-top 10000` (el número equivocado) ahora **se niega** en vez de escribir un cero.