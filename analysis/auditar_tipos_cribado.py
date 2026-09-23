#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""auditar_tipos_cribado.py — la lupa de los tipos de átomo en los ligandos que
YA se prepararon, para encontrar descartes silenciosos como el del boro.

POR QUE
-------
El 23 sep 2026 apareció que `DEC_CHEMBL4543460`, un señuelo del conjunto de
validación de TDP-43, **nunca se había acoplado**: lleva un ácido borónico y Vina
lo rechaza en el parseo (`Atom type B is not a valid AutoDock type`). No dio
error en ningún sitio: simplemente no había pose. Este script busca cuántos más
hay así en la librería del cribado, y busca también el problema contrario, que
es peor: un átomo de boro **nombrado** B pero **tipado** C (o cualquier otro
tipo), porque entonces la malla se construyó como si fuera carbono y la molécula
se acopló y puntuó como su análogo de carbono.

QUE MIRA, POR CADA PDBQT DE ENTRADA
-----------------------------------
Para cada línea ATOM/HETATM se leen dos cosas:
  * el NOMBRE del átomo (columnas 13-16), que dice qué elemento es de verdad;
  * el TIPO (última columna), que es lo que manda al puntuar.
Y se comparan contra dos listas:

  TIPOS_QUE_VINA_ACEPTA   — si aparece un tipo que no está aquí, ese ligando
                            Vina no lo acopla (descarte silencioso).
  ELEMENTOS_QUE_NO_VAN_A_TIPO — elementos que nunca deben acabar como tipo AD
                            (B, Se, Si, Sn, As, Te, Ge...): si aparecen como
                            tipo, el motor no los entiende.

QUE NO HACE
-----------
No repara nada ni toca los ficheros: solo cuenta y lista. La decisión de qué
hacer con lo que salga es otra conversación.

Uso:
    python auditar_tipos_cribado.py                 # todas las carpetas de la lista
    python auditar_tipos_cribado.py --carpeta X     # solo una
Salida: analysis/_auditoria_tipos_cribado.log y un resumen en pantalla.
"""
import argparse
import collections
import os
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
LOG = os.path.join(BASE, "_auditoria_tipos_cribado.log")

# Lo que Vina (y Vina-GPU) saben leer. Si un tipo no esta aqui, el ligando se
# cae en el parseo, sin avisar al que mira el resultado.
TIPOS_QUE_VINA_ACEPTA = {
    "H", "HD", "C", "A", "N", "NA", "OA", "O", "S", "SA", "P", "F", "Cl", "Br",
    "I", "Mg", "Mn", "Zn", "Ca", "Fe",
}

# Elementos cuyo TIPO no existe en AutoDock: si salen como tipo, descarte seguro.
ELEMENTOS_QUE_NO_VAN_A_TIPO = {
    "B", "SE", "SI", "SN", "AS", "TE", "GE", "SB", "BI", "PT", "PD", "HG", "AU",
    "AG", "CD", "PB", "AL", "TI", "CR", "CO", "NI", "CU", "LI", "NA", "K", "RB",
    "CS", "SR", "BA", "BE", "V", "W", "MO", "NB", "TA", "RE", "IR", "OS", "RU",
    "RH", "ZR", "HF", "SC", "Y", "LA", "CE", "PR", "ND", "SM", "EU", "GD", "TB",
    "DY", "HO", "ER", "TM", "YB", "LU",
}

# Elementos que, si van NOMBRADOS como tales pero tipados de otra cosa, cambian la
# quimica de la molecula sin decirlo.
ELEMENTOS_CRITICOS_EN_EL_NOMBRE = ("B", "SE", "SI", "SN", "AS", "TE", "GE")

CARPETAS = [
    os.path.join(RAIZ, "gpu_dock", "libreria_ligands"),
    os.path.join(RAIZ, "gpu_dock", "libreria_ligands_reparados"),
]


def log(m):
    print(m, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")


def mirar(ruta):
    """(tipos_raros, casos_criticos) de un PDBQT.

    tipos_raros:    tipos que Vina no acepta  -> descarte silencioso
    casos_criticos: (elemento, nombre, tipo) de atomos cuyo NOMBRE dice un
                    elemento critico y cuyo TIPO es otra cosa
    """
    raros = collections.Counter()
    criticos = []
    with open(ruta, encoding="utf-8", errors="replace") as f:
        for l in f:
            if not l.startswith(("ATOM", "HETATM")):
                continue
            nombre = l[12:16].strip().upper()
            tipo = l.split()[-1] if l.split() else ""
            tipo_u = tipo.upper()
            if tipo and tipo not in TIPOS_QUE_VINA_ACEPTA:
                raros[tipo] += 1
            raiz = "".join(c for c in nombre if c.isalpha())
            if raiz in ELEMENTOS_CRITICOS_EN_EL_NOMBRE and tipo_u != raiz:
                criticos.append((raiz, nombre, tipo))
    return raros, criticos


def revisar_carpeta(carpeta, limite=0):
    if not os.path.isdir(carpeta):
        log("   !! no existe %s" % carpeta)
        return None
    ficheros = sorted(n for n in os.listdir(carpeta) if n.endswith(".pdbqt"))
    if limite:
        ficheros = ficheros[:limite]
    log("")
    log("=" * 74)
    log("CARPETA %s   (%d ficheros)" % (carpeta, len(ficheros)))
    log("=" * 74)

    con_tipo_raro, con_critico = [], []
    cuenta_tipos = collections.Counter()
    t0 = time.time()
    for i, n in enumerate(ficheros, 1):
        if i % 20000 == 0:
            log("   ...%d de %d (%.0f min)" % (i, len(ficheros),
                                               (time.time() - t0) / 60.0))
        raros, criticos = mirar(os.path.join(carpeta, n))
        if raros:
            con_tipo_raro.append((n[:-len(".pdbqt")], dict(raros)))
            for t, c in raros.items():
                cuenta_tipos[t] += c
        if criticos:
            con_critico.append((n[:-len(".pdbqt")], criticos[:4]))

    log("   tipos que Vina NO acepta (descartes silenciosos): %d ligandos"
        % len(con_tipo_raro))
    if cuenta_tipos:
        log("      tipos encontrados: %s"
            % ", ".join("%s x%d" % (t, c) for t, c in cuenta_tipos.most_common()))
        for nombre, tipos in con_tipo_raro[:40]:
            log("      %-46s %s" % (nombre[:46], tipos))
        if len(con_tipo_raro) > 40:
            log("      ...y %d mas" % (len(con_tipo_raro) - 40))
    log("   elementos criticos NOMBRADOS pero tipados de otra cosa: %d ligandos"
        % len(con_critico))
    for nombre, detalle in con_critico[:40]:
        log("      %-46s %s" % (nombre[:46], detalle))
    if len(con_critico) > 40:
        log("      ...y %d mas" % (len(con_critico) - 40))

    return {"carpeta": carpeta, "ficheros": len(ficheros),
            "tipo_raro": len(con_tipo_raro), "critico": len(con_critico),
            "tipos": dict(cuenta_tipos)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--carpeta", default=None)
    ap.add_argument("--limite", type=int, default=0)
    args = ap.parse_args()

    open(LOG, "w", encoding="utf-8").close()
    log("AUDITORIA DE TIPOS DE ATOMO EN EL CRIBADO   %s"
        % time.strftime("%Y-%m-%d %H:%M"))
    log("Un tipo fuera de la lista de AutoDock = el ligando no se acopla y no avisa.")
    carpetas = [args.carpeta] if args.carpeta else CARPETAS
    resumen = []
    for c in carpetas:
        r = revisar_carpeta(c, args.limite)
        if r:
            resumen.append(r)
    log("")
    log("=" * 74)
    log("RESUMEN")
    log("=" * 74)
    for r in resumen:
        log("   %-52s %6d ficheros | %d descartes | %d criticos"
            % (os.path.basename(r["carpeta"]), r["ficheros"], r["tipo_raro"],
               r["critico"]))
    log("")
    log("log: %s" % LOG)
    return 0


if __name__ == "__main__":
    sys.exit(main())
