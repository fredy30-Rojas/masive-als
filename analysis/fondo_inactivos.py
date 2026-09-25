#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""fondo_inactivos.py — arma el tercer fondo: compuestos MEDIDOS que no unen.

QUE ES Y QUE NO ES
------------------
Los fondos del proyecto son dos, y los dos son parecidos, no inactivos: los señuelos
emparejados en propiedades y los unidores de ARN de R-BIND 2.0. De ninguno de ellos se
sabe que no unan a esta caja; se sabe que se parecen (al primero) o que unen a otra
cosa (al segundo). Eso deja una pregunta sin poder contestar: ¿el embudo ordena los
positivos medidos por delante de compuestos que se midieron y NO unieron?

Ese tercer fondo no se puede inventar —solo existe cuando alguien mide y el resultado
es negativo—, pero **si se puede dejar montado**: los compuestos que el propio cribado
propone medir ya estan elegidos en `peticion_ensayos.csv`, y lo unico que le falta a
cada uno es el resultado. Este script hace esa parte de ahora:

  1. prepara y acopla esos compuestos **con el mismo receptor, la misma caja y la misma
     exhaustividad que la corrida de validacion de su diana**, para que sus numeros
     sean comparables con los que ya estan publicados (es la misma regla que se siguio
     con el fondo duro en `acoplar_fondo_rbind.py`, que este script reutiliza);
  2. deja el manifiesto en `_medidos_no_unen/manifiesto.csv`, con una fila por compuesto
     y su estado: `pendiente_medicion` hasta que el laboratorio diga algo;
  3. deja escrita la declaracion de como se leera ese fondo, **antes** de tener el
     primer dato (ver `DECLARACION.md` en la misma carpeta).

LO QUE NO HACE, A PROPOSITO
---------------------------
No escribe en `verdad_de_referencia.csv` ni marca a nadie como inactivo: sin un ensayo
que lo diga, todos son `pendiente_medicion`. Y no mezcla este fondo con los otros dos:
un compuesto medido y negativo no entra en la lista de señuelos emparejados ni en el
fondo duro. Son tres bloques y se leen por separado.

Uso:
    python analysis/fondo_inactivos.py                 # prepara y acopla lo que falte
    python analysis/fondo_inactivos.py --solo-listar   # solo el manifiesto, sin acoplar
"""

import argparse
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
VINA = os.path.join(RAIZ, "tools", "vina.exe")

PETICION = os.path.join(BASE, "peticion_ensayos.csv")
SALIDA = os.path.join(BASE, "_medidos_no_unen")
MANIFIESTO = os.path.join(SALIDA, "manifiesto.csv")
LOG = os.path.join(SALIDA, "acoplar.log")

# Los MISMOS cuatro numeros que la corrida de validacion de cada diana
# (`validar_diana_limpia.py`): si cambia uno, los numeros dejan de ser comparables.
CAJAS = {
    "SOD1": {
        "receptor": os.path.join(RAIZ, "gpu_dock", "SOD1_limpio.pdbqt"),
        "centro": (46.5, 80.0, 73.3), "tamano": 22, "exhaustividad": 8,
        "nombre_caja": "Trp32 (PDB 1HL5)",
        "sitios_que_valen": ["trp32", "trp 32", "triptofano 32"],
    },
    "TDP-43": {
        "receptor": os.path.join(BASE, "_tdp43_bolsillo_v2", "4BS2_ph74.pdbqt"),
        "centro": (24.23, 16.89, -15.87), "tamano": 26, "exhaustividad": 16,
        "nombre_caja": "RRM1-RRM2 / interfaz Arg151-Asp247 (PDB 4BS2)",
        "sitios_que_valen": ["rrm1", "rrm", "interfaz"],
    },
}

# Piso declarado por adelantado: por debajo de esto el bloque se lee en descriptivo
# (cuantos y en que puesto), no como fondo con el que decidir. Ver DECLARACION.md.
PISO = 30

HILOS = 6

sys.path.insert(0, RAIZ)
from preparar_ligando import escribir   # noqa: E402


def log(m, limpiar=False):
    print(m, flush=True)
    os.makedirs(SALIDA, exist_ok=True)
    with open(LOG, "a" if not limpiar else "w", encoding="utf-8") as f:
        f.write(m + "\n")


def dirs(diana):
    d = os.path.join(SALIDA, diana)
    return d, os.path.join(d, "ligands"), os.path.join(d, "out")


def nombre_de(lig):
    """Nombre de fichero: MED_<compuesto>, para no confundirlo con ACT_ ni DEC_."""
    return "MED_" + "".join(c if c.isalnum() or c in "-_" else "_" for c in lig)


def tiene_pose(ruta):
    if not os.path.exists(ruta):
        return False
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            return True
    return False


def afinidad(ruta):
    if not tiene_pose(ruta):
        return None
    with open(ruta, encoding="utf-8", errors="ignore") as f:
        for l in f:
            if l.startswith("REMARK VINA RESULT:"):
                return float(l.split()[3])
    return None


def _dock(t):
    vina, rec, lig, out, cx, cy, cz, tam, exh = t
    subprocess.run([vina, "--receptor", rec, "--ligand", lig,
                    "--center_x", str(cx), "--center_y", str(cy),
                    "--center_z", str(cz),
                    "--size_x", str(tam), "--size_y", str(tam), "--size_z", str(tam),
                    "--exhaustiveness", str(exh), "--num_modes", "3",
                    "--out", out, "--cpu", "1"], capture_output=True, timeout=7200)
    return out


def candidatos():
    """Los compuestos de la peticion que alimentarian este fondo.

    Devuelve (del_bloque, condicionales): los del bloque son los que ya se piden
    contra la caja de una diana; los condicionales son los que se piden con el sitio
    por determinar y solo entrarian si el ensayo dice que el sitio es ese.
    """
    if not os.path.exists(PETICION):
        raise SystemExit("falta %s: corra antes peticion_ensayos.py" % PETICION)
    del_bloque, condicionales, fuera = [], [], []
    for r in csv.DictReader(open(PETICION, encoding="utf-8")):
        diana = r["diana"]
        if diana not in CAJAS:
            fuera.append(r)
            continue
        caja = CAJAS[diana]
        sitio = r["sitio"].lower()
        # El orden importa: un "por determinar (Trp32, Cys111 u otro)" lleva dentro el
        # nombre de la caja, y no por eso deja de estar por determinar. Primero se mira
        # si el sitio se sabe, y solo despues si es el de la caja.
        if "por determinar" in sitio or "n/d" in sitio or not sitio.strip():
            condicionales.append(r)
        elif any(k in sitio for k in caja["sitios_que_valen"]):
            del_bloque.append(r)
        else:
            # Medido en otra caja (el CR del C-terminal, el CTD): no es de este bloque.
            fuera.append(r)
    return del_bloque, condicionales, fuera


def main():
    global HILOS
    ap = argparse.ArgumentParser()
    ap.add_argument("--solo-listar", action="store_true")
    ap.add_argument("--hilos", type=int, default=HILOS)
    args = ap.parse_args()
    HILOS = args.hilos

    del_bloque, condicionales, fuera = candidatos()
    log("FONDO DE INACTIVIDAD MEDIDA — manifiesto y acoplamiento   %s"
        % time.strftime("%Y-%m-%d %H:%M"), limpiar=True)
    log("compuestos pedidos contra la caja de una diana: %d" % len(del_bloque))
    log("  con el sitio por determinar (solo entran si el ensayo lo confirma): %d"
        % len(condicionales))
    log("  de otra caja (CR del C-terminal, CTD): %d" % len(fuera))

    # ------------------------------------------------------------- preparacion
    filas, tareas, fallos = [], [], []
    t0 = time.time()
    for r in del_bloque:
        diana, lig = r["diana"], r["compuesto"]
        caja = CAJAS[diana]
        _, ligdir, outdir = dirs(diana)
        os.makedirs(ligdir, exist_ok=True)
        os.makedirs(outdir, exist_ok=True)
        nombre = nombre_de(lig)
        out = os.path.join(outdir, nombre + "_out.pdbqt")
        if not args.solo_listar and not tiene_pose(out):
            p = escribir(nombre, r["smiles"], ligdir, quitar_sales=True)
            if p is None:
                fallos.append((lig, "no se pudo preparar"))
            else:
                tareas.append((VINA, caja["receptor"], p, out, caja["centro"][0],
                               caja["centro"][1], caja["centro"][2], caja["tamano"],
                               caja["exhaustividad"]))
        filas.append({
            "diana": diana, "caja": caja["nombre_caja"], "compuesto": lig,
            "identificador": r["identificador"], "smiles": r["smiles"],
            "pesados": "", "de_donde_sale": r["de_donde_sale"],
            "condicional": "no", "estado": "pendiente_medicion",
            "pose": "si" if tiene_pose(out) else "no", "afinidad": "",
            "en_el_fondo": "si",
        })
    for r in condicionales:
        diana = r["diana"]
        filas.append({
            "diana": diana, "caja": CAJAS[diana]["nombre_caja"] + " (si el ensayo lo dice)",
            "compuesto": r["compuesto"], "identificador": r["identificador"],
            "smiles": r["smiles"], "pesados": "", "de_donde_sale": r["de_donde_sale"],
            "condicional": "si", "estado": "pendiente_medicion", "pose": "no",
            "afinidad": "", "en_el_fondo": "no (sitio por determinar)",
        })

    # ------------------------------------------------------------ acoplamiento
    log("")
    log("a acoplar: %d  (fallos de preparacion: %d)" % (len(tareas), len(fallos)))
    for lig, m in fallos:
        log("   %-30s %s" % (lig, m))
    if tareas and not args.solo_listar:
        t0 = time.time()
        hechos = 0
        with ProcessPoolExecutor(max_workers=HILOS) as ex:
            for _ in ex.map(_dock, tareas):
                hechos += 1
                if hechos % 5 == 0 or hechos == len(tareas):
                    log("   acoplados %d/%d  (%.1f min)"
                        % (hechos, len(tareas), (time.time() - t0) / 60.0))

    # --------------------------------------------------------------- manifiesto
    # Los numeros que ya se pueden poner: atomos pesados y afinidad si hay pose.
    for f in filas:
        if f["condicional"] == "si":
            continue
        _, _, outdir = dirs(f["diana"])
        out = os.path.join(outdir, nombre_de(f["compuesto"]) + "_out.pdbqt")
        a = afinidad(out)
        f["pose"] = "si" if a is not None else "no"
        f["afinidad"] = "" if a is None else "%.3f" % a
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
        for f in filas:
            m = Chem.MolFromSmiles(f["smiles"])
            if m is not None:
                f["pesados"] = m.GetNumHeavyAtoms()
    except Exception:
        pass

    cols = ["diana", "caja", "compuesto", "identificador", "smiles", "pesados",
            "de_donde_sale", "condicional", "estado", "pose", "afinidad", "en_el_fondo"]
    os.makedirs(SALIDA, exist_ok=True)
    with open(MANIFIESTO, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for f in filas:
            w.writerow(f)

    listos = [f for f in filas if f["pose"] == "si" and f["condicional"] == "no"]
    log("")
    log("MANIFIESTO: %d compuestos, de los cuales %d ya tienen pose y afinidad en su caja"
        % (len(filas), len(listos)))
    log("   %s" % MANIFIESTO)
    log("")
    log("PARA QUE SIRVE: el dia que llegue el primer ensayo negativo, ese compuesto pasa de")
    log("   'pendiente_medicion' a 'no_une' con su ensayo y su condicion, y el bloque queda")
    log("   legible. Con menos de %d compuestos medidos se lee en descriptivo, no como fondo"
        % PISO)
    log("   (esta dicho por adelantado en DECLARACION.md, en esta misma carpeta).")
    log("   Los numeros NO se tocan: la afinidad es la del MISMO receptor, caja y")
    log("   exhaustividad que la validacion de esa diana, y por eso son comparables.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
