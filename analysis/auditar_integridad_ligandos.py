#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""auditar_integridad_ligandos.py — barrido de integridad de TODOS los ligandos.

POR QUE
-------
Este proyecto ha tenido dos formas de fichero de ligando envenenado, y las dos han
costado una conclusion falsa:

  1. **Pseudo-atomos «glue»** (tipos CG0/G0): meeko no cerraba los anillos de 7
     eslabones en adelante, los abria y pegaba los extremos con atomos de
     pegamento. Vina veia el anillo abierto y regalaba 3-6 kcal/mol. Estaba en los
     primeros puestos de TODAS las validaciones. Ver
     `INFORME_GLUE_VALIDACIONES_2026-09-21.md`.
  2. **Atomos que no cuadran con el SMILES** (sales que se cuelan, anillos mal
     escritos, `sanear_tipos` traduciendo CG0 a carbono): el rescoring no puede
     emparejar la pose con la molecula y falla, o la acopla con dos carbonos de
     mas.

Las dos se detectan con la misma comprobacion: **atomos pesados del fichero = del
SMILES, y cero pseudo-atomos**. Y las dos se pueden comprobar en todas las carpetas
en un par de minutos, que es lo que hace este script antes de que nadie vuelva a
mirar un numero.

QUE COMPRUEBA
-------------
Para cada carpeta `.../ligands` y cada carpeta de poses `.../out` del arbol:
  * cuenta atomos pesados y pseudo-atomos (tipos G0/CG0) de cada fichero;
  * lee el SMILES del propio fichero (`REMARK SMILES`) o de la tabla que se le
    indique, y compara;
  * para las poses, comprueba tambien que la PRIMERA pose cuadra.

Salida: `auditar_integridad_ligandos.log` y `.csv`, con una fila por fichero que
falla. Si no falla nada, lo dice y se acaba.

Uso:
    python auditar_integridad_ligandos.py
    python auditar_integridad_ligandos.py --reparar    # arregla lo que encuentre
"""
import argparse
import csv
import glob
import os
import shutil
import sys

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
LOG = os.path.join(BASE, "auditar_integridad_ligandos")

sys.path.insert(0, RAIZ)
from preparar_ligando import construir, contar_atomos   # noqa: E402

# Carpetas que NO se auditan, con el motivo:
#   * `_roto_glue` / `_antes_glue`: copias a proposito de los ficheros rotos, que
#     se conservan para poder reproducir el fallo.
#   * `tanda*`, `test*`, `_tmp*`, `libreria*`: trozos de la libreria del cribado.
#     Sus ficheros los escribe el conversor de Vina-GPU y no llevan
#     `REMARK SMILES`, asi que aqui no se pueden comprobar; esa parte ya tiene su
#     propia herramienta (`gpu_dock/reparar_libreria_glue.py`, 3.803 afectados).
EXCLUIR = ("_roto_glue", "_antes_glue", "__pycache__",
           "tanda", "test1", "test2", "test5", "_tmp", "libreria",
           "resultados_libreria")


def smi_de_texto(txt):
    """SMILES de un PDBQT ya leido. (Se lee del texto y no del fichero porque las
    poses pueden traer bytes nulos y abrir otra vez es tirar el trabajo.)"""
    for l in txt.splitlines():
        if l.startswith("REMARK SMILES "):
            return l[len("REMARK SMILES "):].strip()
        if l.startswith(("ATOM", "HETATM")):
            break
    return None


def pesados_smiles(smi):
    """Atomos pesados del SMILES. `GetNumHeavyAtoms` y no `GetNumAtoms`: los
    analogos deuterados llevan los [2H] explicitos y no son pesados."""
    try:
        m = Chem.MolFromSmiles(smi) if smi else None
    except Exception:
        return None
    return None if m is None else m.GetNumHeavyAtoms()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reparar", action="store_true",
                    help="re-prepara con la receta canonica los ficheros que fallen")
    ap.add_argument("--max-ficheros", type=int, default=2000,
                    help="por carpeta: mas que esto se salta y se dice (la libreria "
                         "entera son 66.475 ficheros y se audita aparte)")
    ap.add_argument("--solo", default=None,
                    help="subcadena de carpeta a auditar (p.ej. 'validacion')")
    args = ap.parse_args()

    print("=" * 78)
    print("AUDITORIA DE INTEGRIDAD DE LIGANDOS Y POSES")
    print("=" * 78)

    ligandos = sorted(glob.glob(os.path.join(RAIZ, "**", "ligands"), recursive=True))
    poses = sorted(glob.glob(os.path.join(RAIZ, "**", "out"), recursive=True))
    ligandos = [d for d in ligandos if not any(x in d for x in EXCLUIR)]
    poses = [d for d in poses if not any(x in d for x in EXCLUIR)]
    if args.solo:
        ligandos = [d for d in ligandos if args.solo in d]
        poses = [d for d in poses if args.solo in d]
    print("carpetas de ligandos: %d | carpetas de poses: %d"
          % (len(ligandos), len(poses)))

    fallos, revisados, saltadas = [], 0, []

    # Una carpeta con decenas de miles de ficheros (la libreria del cribado) tarda
    # horas en esto; se salta y se dice, en vez de dejar el script colgado sin
    # explicar nada. Para esa hay `gpu_dock/reparar_libreria_glue.py`, que ya se
    # paso el 18 de septiembre.
    def demasiado_grande(d):
        n = len(glob.glob(os.path.join(d, "*.pdbqt")))
        if n > args.max_ficheros:
            saltadas.append((d.replace(RAIZ + os.sep, ""), n))
            return True
        return False

    ligandos = [d for d in ligandos if not demasiado_grande(d)]
    poses = [d for d in poses if not demasiado_grande(d)]
    if saltadas:
        print("")
        for d, n in saltadas:
            print("   SALTADA %-52s %d ficheros (mas de %d)"
                  % (d[:52], n, args.max_ficheros))
        print("")

    sin_smiles_total = []
    for carpeta in ligandos:
        n_sin_smi = 0
        for p in sorted(glob.glob(os.path.join(carpeta, "*.pdbqt"))):
            revisados += 1
            nombre = os.path.basename(p)
            txt = open(p, encoding="utf-8", errors="ignore").read()
            n, glue = contar_atomos(txt)
            smi = smi_de_texto(txt)
            n_smi = pesados_smiles(smi)
            if n_smi is None:
                n_sin_smi += 1
                sin_smiles_total.append((carpeta, nombre, n, glue))
                # los pseudo-atomos se pueden detectar SIN SMILES: eso si se mira
                if glue:
                    fallos.append({"carpeta": carpeta.replace(RAIZ + os.sep, ""),
                                   "fichero": nombre, "tipo": "ligando",
                                   "atomos": n, "glue": glue, "smiles": "",
                                   "motivo": "%d pseudo-atomos (sin SMILES para "
                                             "comprobar atomos)" % glue})
                continue
            if glue or n != n_smi:
                fallos.append({"carpeta": carpeta.replace(RAIZ + os.sep, ""),
                               "fichero": nombre, "tipo": "ligando",
                               "atomos": n, "glue": glue, "smiles": n_smi,
                               "motivo": ("%d pseudo-atomos" % glue) if glue
                                         else "atomos %d != SMILES %d" % (n, n_smi)})
    if sin_smiles_total:
        print("")
        print("ficheros sin `REMARK SMILES` (no se puede comparar atomos, pero SI se")
        print("comprobo que no tuvieran pseudo-atomos): %d, en %d carpetas"
              % (len(sin_smiles_total),
                 len({c for c, _, _, _ in sin_smiles_total})))
        con_glue = [r for r in sin_smiles_total if r[3]]
        print("   de ellos, con pseudo-atomos: %d" % len(con_glue))
        if con_glue:
            for c, nombre, n, glue in con_glue[:20]:
                ruta = c.replace(RAIZ + os.sep, "")
                print("      %-46s %-28s %d pseudo-atomos"
                      % (ruta[:46], nombre[:28], glue))

    for carpeta in poses:
        for p in sorted(glob.glob(os.path.join(carpeta, "*_out.pdbqt"))):
            revisados += 1
            nombre = os.path.basename(p)
            txt = open(p, encoding="utf-8", errors="ignore").read()
            primera = txt.split("ENDMDL")[0] if "ENDMDL" in txt else txt
            n, glue = contar_atomos(primera)
            smi = smi_de_texto(txt)
            n_smi = pesados_smiles(smi)
            if n_smi is None:
                continue
            if glue or n != n_smi:
                fallos.append({"carpeta": carpeta.replace(RAIZ + os.sep, ""),
                               "fichero": nombre, "tipo": "pose",
                               "atomos": n, "glue": glue, "smiles": n_smi,
                               "motivo": ("%d pseudo-atomos" % glue) if glue
                                         else "atomos %d != SMILES %d" % (n, n_smi)})

    print("")
    print("ficheros revisados: %d | con problema: %d" % (revisados, len(fallos)))
    with open(LOG + ".csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["carpeta", "fichero", "tipo", "atomos",
                                          "glue", "smiles", "motivo"])
        w.writeheader()
        for r in fallos:
            w.writerow(r)

    if not fallos:
        print("NADA QUE CORREGIR: ningun fichero con pseudo-atomos ni con un numero")
        print("de atomos distinto al de su SMILES.")
        return 0

    print("")
    for r in fallos[:60]:
        print("   %-46s %-30s %s"
              % (r["carpeta"][:46], r["fichero"][:30], r["motivo"]))
    if len(fallos) > 60:
        print("   ... y %d mas (ver el CSV)" % (len(fallos) - 60))

    if args.reparar:
        print("")
        print("reparando con la receta canonica...")
        arreglados, imposibles = 0, []
        for r in fallos:
            if r["tipo"] != "ligando" or "sin SMILES" in r["motivo"]:
                continue
            p = os.path.join(RAIZ, r["carpeta"], r["fichero"])
            smi = smi_de_texto(open(p, encoding="utf-8", errors="ignore").read())
            txt, n_ok, motivo = construir(smi, quitar_sales=True)
            if txt is None:
                imposibles.append((r["fichero"], motivo))
                continue
            shutil.copy(p, p + ".roto")
            with open(p, "w", encoding="utf-8") as fh:
                fh.write(txt)
            arreglados += 1
            print("   %-34s -> %d atomos, 0 pseudo-atomos"
                  % (r["fichero"][:34], n_ok))
        print("")
        print("reparados: %d | no reparables: %d" % (arreglados, len(imposibles)))
        for nombre, motivo in imposibles:
            print("   %s: %s" % (nombre, motivo))
        print("las poses de los reparados hay que re-acoplarlas (se borran):")
        for r in fallos:
            if r["tipo"] != "ligando":
                continue
            base = os.path.splitext(r["fichero"])[0]
            for c in poses:
                for cand in (os.path.join(c, base + "_out.pdbqt"),):
                    if os.path.exists(cand):
                        os.remove(cand)
                        print("   borrada %s" % cand.replace(RAIZ + os.sep, ""))
    else:
        print("")
        print("vuelve a lanzarlo con --reparar para arreglarlo.")

    with open(LOG + ".log", "w", encoding="utf-8") as f:
        f.write("revisados %d, fallos %d\n" % (revisados, len(fallos)))
        for r in fallos:
            f.write("%s | %s | %s\n" % (r["carpeta"], r["fichero"], r["motivo"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
