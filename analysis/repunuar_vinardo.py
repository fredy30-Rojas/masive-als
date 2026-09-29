#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""repunuar_vinardo.py — re-puntua las poses YA calculadas con Vinardo y con Vina.

POR QUE
------
Los tres fragmentos de RRM2 caen al fondo de la puntuacion de TDP-43. Ya se midio que
no es la caja (cubre el sitio) ni una penalizacion especifica (se puntuan como sus
pares de tamano). El diagnostico que queda es que la funcion de Vina premia la
superficie enterrada, y los fragmentos tienen 11-15 atomos pesados frente a 25-34 de
los demas positivos.

Vinardo es la alternativa que este mismo vina trae: corrige el peso del enlace
carbono-carbono y el radio de la parte no polar, y esta disenada justamente para que la
puntuacion no dependa tanto del numero de atomos. La prueba es directa: si con Vinardo
los fragmentos suben, el problema era la funcion; si se quedan igual, el problema es el
rango de la diana.

QUE HACE Y QUE NO
-----------------
**No re-acopla nada.** Usa `--score_only` sobre las poses que ya estan calculadas:

  * la POSE es la misma que con Vina, y lo unico que cambia es la puntuacion. Eso es
    justo lo que se quiere medir, y ademas sale en minutos en vez de horas;
  * no toca la GPU ni la tarjeta.

La comparacion es justa porque se re-puntua TODO con la misma funcion: los fragmentos,
los demas positivos y los dos fondos. Puntuarlos solo a ellos daria un numero sin
contraponente, que es justo el error que ya se cometio una vez con GNINA (una sola pose
puntuada no es una segunda opinion).

TRAMPAS DE ESTE VINA (1.2.3), LAS CUATRO ENCONTRADAS
---------------------------------------------------
1. `--score_only` **exige caja** (`--center_*` y `--size_*`), y ademas exige que el
   ligando este DENTRO de ella.
2. No acepta `-r`, ni `-l`, ni `-o`: son opciones DOBLES. Con una letra falla al
   instante con "unrecognised option".
3. **El formato del ligando es la trampa gorda.** Un `_out.pdbqt` de Vina trae
   `MODEL`/`ENDMDL` (tres modelos): con eso da "Unexpected multi-MODEL tag". Si se
   quitan esas lineas y se dejan solo los ATOM, da "Unknown or inappropriate tag in
   flex residue or ligand" porque se quedan los `ROOT`/`BRANCH` colgando. Lo que este
   binario quiere es un unico `ROOT ... ENDROOT` con sus `BRANCH`, o los ATOM y un
   `TORSDOF`. Se construye aqui, sin tocar las poses originales.
4. La afinidad sale como `Estimated Free Energy of Binding`, no como `Affinity`, en
   `--score_only`. Con la busqueda normal si sale como `Affinity`.

DIANAS: DE DONDE SALE CADA COSA
-------------------------------
El receptor, el tamano de caja y la carpeta del fondo blando se leen de
`validar_diana_limpia.DIANAS`, que es la unica copia de "que usa cada diana".
Escribirlo aqui aparte ya lo hizo mal una vez: con CDK2 se puntuo contra el receptor
viejo.

El fondo duro para SOD1 NO esta en `DIANAS` (alla esta a None porque lo creo
`validar_sod1_limpia.py`): se lee de su `LIGS_FONDO2`, no de una ruta escrita a mano.

Y el receptor y el tamano van DENTRO de la tarea, no en variables globales: en Windows
el ProcessPoolExecutor arranca los hijos con `spawn`, que re-importa este modulo y
perderia cualquier valor cambiado desde fuera. Con las constantes de modulo, un
`--target SOD1` habria puntuado contra el receptor de TDP-43 sin decir nada.

Uso:
    python repunuar_vinardo.py --target TDP43
    python repunuar_vinardo.py --target SOD1
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
VINA = os.path.join(RAIZ, "tools", "vina.exe")

TMP = os.path.join(BASE, "_vinardo")
HILOS = 12

H_TIPOS = ("H", "HD", "HS", "D", "DD")


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def un_modelo(ruta):
    """(lineas, coords) del PRIMER modelo de un _out.pdbqt de Vina.

    Solo el primero: el fichero trae `num_modes` modelos y contar todos daria tres
    veces el numero de atomos, que ya paso una vez.
    """
    lineas, dentro = [], False
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            if dentro:
                break
            dentro = True
            continue
        if not dentro:
            continue
        if l.startswith("ENDMDL"):
            break
        lineas.append(l.rstrip("\n"))
    coords = []
    for l in lineas:
        if l.startswith(("ATOM", "HETATM")):
            coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return lineas, np.array(coords) if coords else np.zeros((0, 3))


def a_ligando_unico(lineas, destino):
    """Escribe un unico ligando en el formato que este vina acepta con score_only.

    Se queda con los ATOM del primer modelo y los envuelve en ROOT/ENDROOT, mas el
    TORSDOF del original. Sin las lineas MODEL/ENDMDL, que son las que provocan el
    "multi-MODEL"; sin las ramas, que es lo que dispara el otro error.
    """
    atoms = [l for l in lineas if l.startswith(("ATOM", "HETATM"))]
    tors = [l for l in lineas if l.startswith("TORSDOF")]
    with open(destino, "w", encoding="utf-8") as f:
        f.write("ROOT\n")
        for l in atoms:
            f.write(l + "\n")
        f.write("ENDROOT\n")
        f.write(tors[0] if tors else "TORSDOF 0.0\n")


def puntuar(tarea):
    """Puntua UNA pose con UNA funcion. Devuelve (nombre, scoring, valor, error)."""
    nombre, ruta_pose, scoring, receptor, tamano = tarea
    tmp = os.path.join(TMP, "%s_%s.pdbqt" % (nombre.replace("/", "_"), scoring))
    try:
        lineas, coords = un_modelo(ruta_pose)
        if len(coords) == 0:
            return nombre, scoring, None, "sin modelo"
        centro = coords.mean(axis=0)
        # La caja va centrada en la POSE, no en la caja del docking: score_only exige
        # que el ligando este dentro, y algunos ligandos grandes se salen del lado
        # bueno de una caja de 26 A puesta en el centro del bolsillo.
        #
        # Y si aun asi se sale (3 ligandos del fondo duro de SOD1), se agranda lo
        # justo para que quepa. NO cambia el resultado: la afinidad solo cuenta
        # atomos del receptor a menos de 8 A del ligando, y todos esos ya estaban
        # dentro de una caja de 22 A centrada en la pose. La caja no añade ni quita
        # atomos del calculo, solo evita el error de "ligand is outside the grid box".
        # Comprobado: con este cambio, la corrida entera de TDP-43 sale byte a byte
        # igual que con la caja fija.
        span = coords.max(axis=0) - coords.min(axis=0)
        lado = max(float(tamano), float(span.max()) + 4.0)
        a_ligando_unico(lineas, tmp)
        cmd = [VINA, "--receptor", receptor, "--ligand", tmp,
               "--score_only", "--scoring", scoring,
               "--center_x", "%.3f" % centro[0], "--center_y", "%.3f" % centro[1],
               "--center_z", "%.3f" % centro[2],
               "--size_x", "%.1f" % lado, "--size_y", "%.1f" % lado,
               "--size_z", "%.1f" % lado]
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
        t = r.stdout + r.stderr
        m = (re.search(r"Estimated Free Energy of Binding\s*:\s*(-?[0-9.]+)", t)
             or re.search(r"Affinity\s*[:=]\s*(-?[0-9.]+)", t))
        if not m:
            return nombre, scoring, None, (t[-120:].replace("\n", " ") or "sin salida")
        return nombre, scoring, float(m.group(1)), ""
    except Exception as e:  # noqa: BLE001
        return nombre, scoring, None, repr(e)[:110]
    finally:
        if os.path.exists(tmp):
            try:
                os.remove(tmp)
            except OSError:
                pass


def recoger(patron, etiqueta, excluir_act=False):
    """(nombre, ruta, etiqueta) de las poses que casan con el patron.

    `recoger` quita el prefijo `ACT_` porque en el CSV los positivos van por su nombre
    desnudo. Por eso, para un fondo hay que pedir `excluir_act=True`: si no, las poses
    de los positivos que estan en la MISMA carpeta del fondo se contarian como si
    fueran senuelos. Es la misma regla que aplica `cargar_fondo` en el validador, y sin
    ella el fondo de SOD1 tendria 31 ligandos de mas.
    """
    rutas = sorted(glob.glob(os.path.join(BASE, patron)))
    out = []
    for r in rutas:
        crudo = os.path.basename(r)
        if excluir_act and crudo.startswith("ACT_"):
            continue
        nombre = crudo.replace("_out.pdbqt", "").replace(".pdbqt", "")
        if nombre.startswith("ACT_"):
            nombre = nombre[4:]
        out.append((nombre, r, etiqueta))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="TDP43", choices=["TDP43", "SOD1"])
    ap.add_argument("--salida", default=None)
    ap.add_argument("--hilos", type=int, default=HILOS)
    args = ap.parse_args()

    sys.path.insert(0, BASE)
    import validar_diana_limpia as D
    import validar_sod1_limpia as V
    cfg = D.DIANAS[args.target]
    receptor = cfg["receptor"]
    tamano = cfg["tamano"]
    salida = args.salida or os.path.join(BASE, "vinardo_%s.csv" % args.target.lower())

    os.makedirs(TMP, exist_ok=True)
    if not os.path.exists(receptor):
        raise SystemExit("no esta el receptor %s" % receptor)
    log("diana %s | receptor %s | tamano de caja %d"
        % (args.target, os.path.basename(receptor), tamano))

    poses2 = cfg.get("poses2") or os.path.join(os.path.dirname(V.LIGS_FONDO2), "out")
    log("fondo blando: %s" % cfg["poses"])
    log("fondo duro:   %s" % poses2)

    # --- cuales de los `ACT_` son positivos de verdad --------------------------
    # No basta con el prefijo: en la carpeta hay 14 poses `ACT_` y solo 7 son de
    # UNION MEDIDA (las de `verdad_de_referencia.csv` con `apto=si`). Las otras 7
    # no son ni positivos ni senuelos: el validador las deja fuera de todo
    # (`cargar_fondo` se salta los `ACT_`). Contarlas como fondo era justo el error
    # que se cuela aqui: 9 de los 131 "fondo blando" de la primera corrida de
    # TDP-43 eran positivos. La lista sale de la tabla, no del prefijo.
    truths = [r["ligand"] for r in csv.DictReader(
        open(os.path.join(BASE, "verdad_de_referencia.csv"), encoding="utf-8"))
        if r["target"] == args.target and r["apto"] == "si"]
    sin_unir = set(truths)
    log("positivos de union medida en la tabla: %d -> %s"
        % (len(truths), ", ".join(sorted(truths))))

    # Los positivos se buscan en DOS sitios, en el orden que usa el validador: la
    # carpeta del fondo, de donde salieron los numeros publicados, y la salida de la
    # corrida, que es donde quedan los re-acoplados. Si esta en los dos, gana la
    # salida (es la version con la que se puntuaron las metricas publicadas).
    out_salida = os.path.join(cfg["salida"], "out")
    positivos, act_descartados = {}, set()
    for d in (out_salida, cfg["poses"]):
        for n, r, _e in recoger(os.path.join(d, "ACT_*_out.pdbqt"), "positivo"):
            if n not in sin_unir:
                act_descartados.add(n)
                continue
            if d == out_salida or n not in positivos:
                positivos[n] = (n, r, "positivo")
    positivos = [positivos[n] for n in sorted(positivos)]
    if act_descartados:
        log("poses ACT_ SIN union medida, fuera de todo (ni positivas ni fondo): "
            "%s" % ", ".join(sorted(act_descartados)))

    blando = recoger(os.path.join(cfg["poses"], "*_out.pdbqt"), "fondo blando",
                     excluir_act=True)
    duro = recoger(os.path.join(poses2, "*_out.pdbqt"), "fondo duro",
                   excluir_act=True)
    log("positivos %d | fondo blando %d | fondo duro %d"
        % (len(positivos), len(blando), len(duro)))
    if not positivos:
        raise SystemExit("ningun positivo con union medida en %s ni en %s"
                         % (out_salida, cfg["poses"]))

    tareas = []
    for grupo in (positivos, blando, duro):
        for n, r, _e in grupo:
            tareas.append((n, r, "vina", receptor, tamano))
            tareas.append((n, r, "vinardo", receptor, tamano))
    log("a re-puntuar: %d (dos funciones cada uno)" % len(tareas))

    filas = {}
    fallos = 0
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=args.hilos) as ex:
        for k, (nombre, scoring, val, err) in enumerate(ex.map(puntuar, tareas), 1):
            filas.setdefault(nombre, {})[scoring] = val
            if val is None:
                fallos += 1
                if fallos <= 6:
                    log("   FALLO %-24s %-8s %s" % (nombre, scoring, err))
            if k % 200 == 0 or k == len(tareas):
                log("   %d/%d (%.0f s, %d fallos)"
                    % (k, len(tareas), time.time() - t0, fallos))

    papel = {}
    for n, _, _e in positivos:
        papel[n] = "positivo"
    for n, _, _e in blando:
        papel.setdefault(n, "fondo blando")
    for n, _, _e in duro:
        papel[n] = "fondo duro"

    # Puestos por grupo, que es como los lee el validador: cada positivo contra SU
    # fondo, no contra la mezcla. Si se mezclaran, el puesto de un positivo dependeria
    # de cuantos ligandos hubiera en el otro fondo, y eso es otra pregunta.
    def puesto_en(fondo, campo):
        """{nombre: puesto} de los positivos y de ESE fondo, ordenados por `campo`."""
        universo = [n for n in filas if papel.get(n) in (fondo, "positivo")
                    and filas[n].get(campo) is not None]
        universo.sort(key=lambda n: filas[n][campo])
        return {n: i + 1 for i, n in enumerate(universo)}

    campos = ["ligand", "papel", "vina", "vinardo", "delta_vinardo",
              "puesto_vina_duro", "puesto_vinardo_duro",
              "puesto_vina_blando", "puesto_vinardo_blando"]
    pv_d = puesto_en("fondo duro", "vina")
    pa_d = puesto_en("fondo duro", "vinardo")
    pv_b = puesto_en("fondo blando", "vina")
    pa_b = puesto_en("fondo blando", "vinardo")
    with open(salida, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for n in filas:
            v, a = filas[n].get("vina"), filas[n].get("vinardo")
            w.writerow({"ligand": n, "papel": papel.get(n, "?"),
                        "vina": "" if v is None else "%.3f" % v,
                        "vinardo": "" if a is None else "%.3f" % a,
                        "delta_vinardo": "" if (v is None or a is None)
                        else "%.3f" % (a - v),
                        "puesto_vina_duro": pv_d.get(n, ""),
                        "puesto_vinardo_duro": pa_d.get(n, ""),
                        "puesto_vina_blando": pv_b.get(n, ""),
                        "puesto_vinardo_blando": pa_b.get(n, "")})
    log("csv: %s" % os.path.basename(salida))
    log("positivos con puesto frente al fondo duro: %d"
        % sum(1 for n in filas if papel.get(n) == "positivo" and n in pa_d))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
