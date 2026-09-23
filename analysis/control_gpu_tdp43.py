#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""control_gpu_tdp43.py — CONTROL: el mismo conjunto de TDP-43 acoplado con
Vina-GPU, para ver si el motor cambia el orden.

QUE PREGUNTA CONTESTA
---------------------
La validacion limpia de TDP-43 (positivos + señuelos emparejados + fondo duro de
R-BIND 2.0) va con Vina de CPU (`tools/vina.exe`), y el fondo duro tiene que ir con
el mismo motor que los positivos o la comparacion mezcla el motor con la quimica.
Este control hace lo otro: acopla **todo el conjunto** (los 285 ligandos: positivos,
señuelos y fondo duro) con **Vina-GPU 2.1** de una vez, y pasa despues las mismas
metricas. Si el orden sale distinto, la diferencia es del motor.

LO QUE HAY QUE DECIR AL REPORTARLO
----------------------------------
`search_depth` de Vina-GPU NO es `exhaustiveness` de Vina: son perillas distintas y
no hay una equivalencia exacta. Aqui se usa search_depth = 20 (la recomendacion del
propio proyecto para 2.1) frente a exhaustiveness 16 en CPU. Por eso el control
compara **el orden y los veredictos**, no los valores absolutos de afinidad.

Los ligandos NO se vuelven a preparar: se copian los `.pdbqt` ya preparados con la
receta canonica (`preparar_ligando.py`), los mismos que usa la corrida de CPU. Lo
unico que cambia es el motor que los acopla.

Uso:
    python control_gpu_tdp43.py --montar          # copia los ligandos
    python control_gpu_tdp43.py --acoplar          # Vina-GPU sobre lo que falte
    python control_gpu_tdp43.py --puntuar          # las mismas metricas
    python control_gpu_tdp43.py --todo             # las tres cosas
    python control_gpu_tdp43.py --acoplar --limite 3   # prueba corta

Salida: `_control_gpu_TDP43/` con out/ (positivos y señuelos), out_duro/ (R-BIND),
`acoplar_gpu.log` y `validar_gpu.log` + `validar_gpu.csv`.
"""
import argparse
import csv
import glob
import os
import shutil
import subprocess
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
GPU = os.path.join(RAIZ, "gpu_dock")
VINA_GPU = os.path.join(GPU, "Vina-GPU-2.1-win.exe")

# La misma diana, la misma caja y la misma receta que la corrida de CPU.
RECEPTOR = os.path.join(BASE, "_tdp43_bolsillo_v2", "4BS2_ph74.pdbqt")
CENTRO = (24.23, 16.89, -15.87)
TAMANO = 26
SEARCH_DEPTH = 20     # NO es la exhaustividad de Vina: es la perilla de Vina-GPU
NUM_MODES = 3
THREAD = 8000

SALIDA = os.path.join(BASE, "_control_gpu_TDP43")
LIGS = os.path.join(SALIDA, "ligands")
OUT = os.path.join(SALIDA, "out")
OUT_DURO = os.path.join(SALIDA, "out_duro")
PEND = os.path.join(SALIDA, "pendientes")
LOG = os.path.join(SALIDA, "acoplar_gpu.log")

ORIGENES = [os.path.join(BASE, "_validacion_TDP43", "ligands"),
            os.path.join(BASE, "_validacion_TDP43_rbind", "ligands"),
            # Los positivos que la corrida de CPU re-prepara por su cuenta (los que
            # no estaban en el conjunto viejo: rTRD01, nTRD22 y los tres fragmentos)
            # viven aqui. Se copian tambien, para que el control sea GPU de arriba
            # abajo y no cuele cinco poses de CPU en la tabla.
            os.path.join(BASE, "validar_tdp43_limpia", "ligands")]

# El fondo duro son los de R-BIND; todo lo demas (positivos y señuelos) es el
# fondo blando. Se separan por el prefijo, que es como se llaman ya.
PREFIJO_DURO = "RB_SM_"


def log(m):
    print(m, flush=True)
    os.makedirs(SALIDA, exist_ok=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write("[%s] %s\n" % (time.strftime("%Y-%m-%d %H:%M:%S"), m))


def pose_de(nombre):
    """Donde esta (o deberia estar) la pose de un ligando."""
    for d in (OUT, OUT_DURO):
        p = os.path.join(d, nombre + "_out.pdbqt")
        if os.path.exists(p) and os.path.getsize(p) > 100:
            return p
    return None


def tiene_resultado(ruta):
    if not ruta or not os.path.exists(ruta):
        return False
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            return True
    return False


def montar():
    """Copia los .pdbqt ya preparados (no se re-prepara nada)."""
    os.makedirs(LIGS, exist_ok=True)
    puestos = 0
    for origen in ORIGENES:
        if not os.path.isdir(origen):
            log("   !! falta el origen %s" % origen)
            continue
        for p in sorted(glob.glob(os.path.join(origen, "*.pdbqt"))):
            destino = os.path.join(LIGS, os.path.basename(p))
            if not os.path.exists(destino):
                shutil.copy2(p, destino)
                puestos += 1
    total = len(glob.glob(os.path.join(LIGS, "*.pdbqt")))
    log("ligandos en el montaje: %d (copiados %d)" % (total, puestos))
    return total


def escribir_config(ligdir, outdir):
    cfg = os.path.join(SALIDA, "config_gpu.txt")
    with open(cfg, "w", encoding="utf-8") as f:
        f.write("receptor = %s\n" % RECEPTOR.replace("\\", "/"))
        f.write("ligand_directory = %s\n" % ligdir.replace("\\", "/"))
        f.write("output_directory = %s\n" % outdir.replace("\\", "/"))
        f.write("opencl_binary_path = %s\n" % GPU.replace("\\", "/"))
        f.write("center_x = %s\ncenter_y = %s\ncenter_z = %s\n" % CENTRO)
        f.write("size_x = %d\nsize_y = %d\nsize_z = %d\n" % (TAMANO, TAMANO, TAMANO))
        f.write("search_depth = %d\nnum_modes = %d\nthread = %d\n"
                % (SEARCH_DEPTH, NUM_MODES, THREAD))
    return cfg


def acoplar(limite=0):
    """Vina-GPU sobre lo que aun no tenga pose. Reanudable."""
    if not os.path.exists(VINA_GPU):
        log("!! no esta Vina-GPU en %s" % VINA_GPU)
        return 1
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(OUT_DURO, exist_ok=True)

    todos = sorted(os.path.basename(p)[:-len(".pdbqt")]
                   for p in glob.glob(os.path.join(LIGS, "*.pdbqt")))
    pendientes = [n for n in todos if not tiene_resultado(pose_de(n))]
    con_pose = len(todos) - len(pendientes)
    if limite:
        pendientes = pendientes[:limite]
    log("ligandos: %d   ya con pose: %d   a acoplar ahora: %d"
        % (len(todos), con_pose, len(pendientes)))
    if not pendientes:
        return 0

    # Cada tanda va a su carpeta: asi lo que devuelve el motor no se mezcla con
    # lo que ya estaba, y el reparto se hace despues por el prefijo.
    shutil.rmtree(PEND, ignore_errors=True)
    ligdir = os.path.join(PEND, "ligands")
    tmpout = os.path.join(PEND, "out")
    os.makedirs(ligdir)
    os.makedirs(tmpout)
    for n in pendientes:
        shutil.copy2(os.path.join(LIGS, n + ".pdbqt"),
                     os.path.join(ligdir, n + ".pdbqt"))

    cfg = escribir_config(ligdir, tmpout)
    log("Vina-GPU 2.1 -> search_depth %d, num_modes %d, caja de %d A"
        % (SEARCH_DEPTH, NUM_MODES, TAMANO))
    t0 = time.time()
    with open(os.path.join(SALIDA, "vina_gpu_stdout.log"), "a",
              encoding="utf-8", errors="ignore") as f:
        f.write("\n===== tanda de %d ligandos, %s =====\n"
                % (len(pendientes), time.strftime("%Y-%m-%d %H:%M")))
        r = subprocess.run([VINA_GPU, "--config", cfg], cwd=GPU,
                           stdout=f, stderr=subprocess.STDOUT, timeout=6 * 3600)
    log("Vina-GPU termino con codigo %s en %.1f min" % (r.returncode,
                                                        (time.time() - t0) / 60.0))

    hechas = 0
    for p in glob.glob(os.path.join(tmpout, "*_out.pdbqt")):
        nombre = os.path.basename(p).replace("_out.pdbqt", "")
        destino = OUT_DURO if nombre.startswith(PREFIJO_DURO) else OUT
        os.replace(p, os.path.join(destino, nombre + "_out.pdbqt"))
        hechas += 1
    log("poses guardadas: %d   (fuera: %d)"
        % (hechas, len(pendientes) - hechas))
    sin = [n for n in pendientes if not tiene_resultado(pose_de(n))]
    if sin:
        log("   sin pose (el motor no las devolvio): %s" % ", ".join(sin[:20]))
    return 0


def puntuar():
    """Las MISMAS metricas que la corrida de CPU, sobre las poses de la GPU."""
    n1 = len(glob.glob(os.path.join(OUT, "*_out.pdbqt")))
    n2 = len(glob.glob(os.path.join(OUT_DURO, "*_out.pdbqt")))
    log("puntuando: %d poses de positivos+señuelos y %d del fondo duro" % (n1, n2))
    sys.path.insert(0, BASE)
    import validar_sod1_limpia as V

    V.RECEPTOR = RECEPTOR
    V.CENTRO = CENTRO
    V.TAMANO = TAMANO
    V.EXHAUSTIVIDAD = SEARCH_DEPTH
    V.POSES_FONDO = OUT
    V.LIGS_FONDO = LIGS
    V.SALIDA = os.path.join(SALIDA, "validar_gpu")
    V.DIANA = "TDP43"
    V.NOTA_FONDO = ("122 señuelos emparejados en propiedades, acoplados con "
                    "Vina-GPU 2.1 (search_depth %d)" % SEARCH_DEPTH)
    V.ETIQUETA_FONDO = "señuelos emparejados (Vina-GPU)"
    V.POSES_FONDO2 = OUT_DURO
    V.LIGS_FONDO2 = LIGS
    V.NOTA_FONDO2 = ("unidores de ARN de R-BIND 2.0, acoplados con Vina-GPU 2.1 "
                     "(search_depth %d) en la misma caja" % SEARCH_DEPTH)
    V.ETIQUETA_FONDO2 = "unidores de ARN de R-BIND 2.0 (Vina-GPU)"
    V.log("CONTROL con Vina-GPU: %s" % time.strftime("%Y-%m-%d %H:%M"))
    return V.main()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--montar", action="store_true")
    ap.add_argument("--acoplar", action="store_true")
    ap.add_argument("--puntuar", action="store_true")
    ap.add_argument("--todo", action="store_true")
    ap.add_argument("--limite", type=int, default=0)
    args = ap.parse_args()
    if not any([args.montar, args.acoplar, args.puntuar, args.todo]):
        ap.error("hay que decir que hacer: --montar, --acoplar, --puntuar o --todo")

    os.makedirs(SALIDA, exist_ok=True)
    log("")
    log("CONTROL VINA-GPU — conjunto completo de TDP-43  %s"
        % time.strftime("%Y-%m-%d %H:%M"))
    if args.todo or args.montar:
        montar()
    if args.todo or args.acoplar:
        if acoplar(args.limite):
            return 1
    if args.todo or args.puntuar:
        return puntuar()
    return 0


if __name__ == "__main__":
    sys.exit(main())
