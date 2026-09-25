#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""La causa de los -10: pseudo-atomos «glue» en dos controles (21 sep 2026).

EL HALLAZGO
-----------
De los 14 controles de la validacion v4, los UNICOS dos que llevan pseudo-atomos
de pegamento (tipos CG0 y G0 que escribe meeko cuando no sabe cerrar un anillo)
son exactamente los dos que dieron el «positivo» del proyecto:

    ACT_diazepanoquinazolina_4MQ   4 pseudo-atomos   -10,02
    ACT_cf3quinazolina_ZO0         4 pseudo-atomos    -9,77
    (los otros 12 controles:        0                 -4,6 a -6,5)

El anillo de diazepano (7 eslabones) sale escrito ABIERTO: meeko lo parte en dos
ramas y pega los extremos con dos atomos fantasma, duplicando dos carbonos (por
eso el fichero tiene 20 atomos donde la molecula tiene 18). Consecuencias:

  * El arbol de torsiones pasa de 2 ramas a 14: Vina cree que la molecula es
    flexible y le da 13,7 torsiones libres ((3) = 3,98 kcal/mol frente a 0,29).
  * La afinidad que Vina reporta es (1) + (2) + (3) - (4), donde (4) es la
    energia interna del ligando LIBRE, medida en el fichero de entrada. Con el
    anillo abierto y dos atomos de mas, esa energia sale +8,87 kcal/mol en vez
    de ~0, y al restarla le regala al compuesto casi 9 kcal/mol de afinidad.
  * El acoplamiento ni se entera: la pose es la misma, la molecula es la misma.

Por eso el «un CH2 no vale 5 kcal/mol» del informe de la v4 no era quimica: era
el fichero de entrada. Y por eso el efecto no dependia del receptor: aparecia
igual en el receptor sucio, en el limpio y en el de la propia estructura.

QUE HACE ESTE SCRIPT
--------------------
1. Detecta que controles llevan pseudo-atomos.
2. Los vuelve a preparar con la receta corregida (`rigid_macrocycles=True`, la
   misma que el proyecto ya usa para reparar la libreria) y comprueba que el
   numero de atomos coincide con el SMILES: si el anillo vuelve cerrado, el
   fichero tiene los atomos que toca.
3. Acopla las dos versiones (rota y reparada) en el MISMO receptor y la MISMA
   caja, con las mismas semillas, y compara.

Uso: python confirmar_glue_controles.py
Salida: confirmar_glue_controles.log y confirmar_glue_controles.csv
"""
import csv
import glob
import os
import re
import subprocess
import sys

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
GPU = os.path.abspath(os.path.join(BASE, "..", "gpu_dock"))
VINA = r"C:\Users\Fredy\masive-als\tools\vina.exe"
CENTRO = (46.5, 80.0, 73.3)
TAMANO = 22
EXHAUSTIVIDAD = 8
SEMILLAS = (42, 2026, 777)
RECEPTORES = {
    "limpio": os.path.join(GPU, "SOD1_limpio.pdbqt"),
    "sucio": os.path.join(GPU, "SOD1.pdbqt"),
}
SIN_H = ("H", "HD", "HS", "D", "DD")


def log(m):
    print(m, flush=True)


def contar(path):
    """(atomos pesados, cuantos pseudo-atomos glue)"""
    n = glue = 0
    for l in open(path, encoding="utf-8", errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        t = l.rsplit(None, 1)[-1].strip()
        if t in SIN_H:
            continue
        n += 1
        if t in ("G0", "CG0"):
            glue += 1
    return n, glue


def smi_del_fichero(path):
    for l in open(path, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK SMILES "):
            return l[len("REMARK SMILES "):].strip()
        if l.startswith(("ATOM", "HETATM")):
            break
    return None


def reparar(origen, destino):
    """Re-prepara el ligando con los anillos rigidos (receta del proyecto)."""
    smi = smi_del_fichero(origen)
    mol = Chem.MolFromSmiles(smi) if smi else None
    if mol is None:
        return None, "sin SMILES valido"
    n_smi = mol.GetNumAtoms()
    mh = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    params.maxIterations = 1000
    ok = False
    for semilla in (42, 2026, 777, 1234, 9999):
        params.randomSeed = semilla
        if AllChem.EmbedMolecule(mh, params) == 0:
            ok = True
            break
    if not ok:
        return None, "no se pudo generar la geometria 3D"
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=200)
    except Exception:
        pass
    setups = MoleculePreparation(rigid_macrocycles=True).prepare(mh)
    if not setups:
        return None, "meeko no encontro parametros"
    txt, bien, err = PDBQTWriterLegacy.write_string(setups[0])
    if not bien:
        return None, "escritura PDBQT: %s" % err
    n, glue = _contar_texto(txt)
    if glue:
        return None, "sigue con pseudo-atomos"
    if n != n_smi:
        return None, "atomos SMILES %d vs preparado %d" % (n_smi, n)
    with open(destino, "w", encoding="utf-8") as f:
        f.write(txt)
    return n_smi, ""


def _contar_texto(txt):
    n = glue = 0
    for l in txt.splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        t = l.rsplit(None, 1)[-1].strip()
        if t in SIN_H:
            continue
        n += 1
        if t in ("G0", "CG0"):
            glue += 1
    return n, glue


def acoplar(receptor, ligando, salida):
    if os.path.exists(salida) and os.path.getsize(salida) > 100:
        return
    subprocess.run([VINA, "--receptor", receptor, "--ligand", ligando,
                    "--center_x", str(CENTRO[0]), "--center_y", str(CENTRO[1]),
                    "--center_z", str(CENTRO[2]),
                    "--size_x", str(TAMANO), "--size_y", str(TAMANO),
                    "--size_z", str(TAMANO),
                    "--exhaustiveness", str(EXHAUSTIVIDAD),
                    "--num_modes", "9", "--seed", "42", "--out", salida],
                   capture_output=True, timeout=3600)


def afinidad(path):
    if not os.path.exists(path):
        return None
    for l in open(path, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            return float(re.search(r"([-+]?\d+\.\d+)", l).group(1))
    return None


def terminos(receptor, ligando):
    """Desglose (1) intermolecular y (4) energia del ligando libre."""
    src = ligando
    solo = os.path.join(BASE, "tmp_pose", "_uno.pdbqt")
    os.makedirs(os.path.dirname(solo), exist_ok=True)
    sal = []
    for l in open(src, encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            continue
        if l.startswith("ENDMDL"):
            break
        if l.startswith(("ATOM", "HETATM", "ROOT", "ENDROOT", "BRANCH",
                         "ENDBRANCH", "TORSDOF")):
            sal.append(l)
    open(solo, "w", encoding="utf-8").writelines(sal)
    r = subprocess.run([VINA, "--receptor", receptor, "--ligand", solo,
                        "--score_only", "--center_x", str(CENTRO[0]),
                        "--center_y", str(CENTRO[1]), "--center_z", str(CENTRO[2]),
                        "--size_x", str(TAMANO), "--size_y", str(TAMANO),
                        "--size_z", str(TAMANO)],
                       capture_output=True, text=True, timeout=1800)
    txt = (r.stdout or "") + (r.stderr or "")
    d = {}
    for clave, pat in (("intermolecular", r"Final Intermolecular Energy\s*:\s*([-+]?\d+\.\d+)"),
                       ("libre", r"Unbound System's Energy\s*:\s*([-+]?\d+\.\d+)")):
        m = re.search(pat, txt)
        d[clave] = float(m.group(1)) if m else None
    return d


def main():
    origenes = sorted(glob.glob(os.path.join(BASE, "validacion_SOD1_v4", "ligands",
                                             "ACT_*.pdbqt")))
    reparados = os.path.join(BASE, "controles_reparados")
    os.makedirs(reparados, exist_ok=True)
    out = os.path.join(BASE, "out")
    os.makedirs(out, exist_ok=True)

    filas, lineas = [], []
    for origen in origenes:
        nombre = os.path.basename(origen).replace(".pdbqt", "")
        n, glue = contar(origen)
        fila = {"control": nombre, "atomos_fichero": n, "pseudo_atomos": glue}
        if not glue:
            filas.append(fila)
            continue
        destino = os.path.join(reparados, nombre + ".pdbqt")
        n_smi, motivo = reparar(origen, destino)
        fila["atomos_smiles"] = n_smi
        fila["reparado"] = "no" if motivo else "si"
        fila["motivo"] = motivo
        lineas.append("%s: %d atomos en el fichero, %d pseudo-atomos -> %s"
                      % (nombre, n, glue, "reparado (%d atomos)" % n_smi if not motivo
                         else "NO reparado: %s" % motivo))
        if motivo:
            filas.append(fila)
            continue
        for receptor_nombre, receptor in RECEPTORES.items():
            for etiqueta, lig in (("roto(v4)", origen), ("reparado", destino)):
                salida = os.path.join(out, "glue_%s_%s_%s.pdbqt"
                                      % (nombre, receptor_nombre, etiqueta))
                acoplar(receptor, lig, salida)
                a = afinidad(salida)
                fila["%s_%s" % (receptor_nombre, etiqueta)] = a
        t_roto = terminos(os.path.join(GPU, "SOD1_limpio.pdbqt"),
                          os.path.join(out, "glue_%s_limpio_roto(v4).pdbqt" % nombre))
        t_rep = terminos(os.path.join(GPU, "SOD1_limpio.pdbqt"),
                         os.path.join(out, "glue_%s_limpio_reparado.pdbqt" % nombre))
        fila["intermol_roto"] = t_roto.get("intermolecular")
        fila["libre_roto"] = t_roto.get("libre")
        fila["intermol_reparado"] = t_rep.get("intermolecular")
        fila["libre_reparado"] = t_rep.get("libre")
        lineas.append("   E(intermolecular): roto %s | reparado %s"
                      % (t_roto.get("intermolecular"), t_rep.get("intermolecular")))
        lineas.append("   E(ligando libre) : roto %s | reparado %s"
                      % (t_roto.get("libre"), t_rep.get("libre")))
        lineas.append("   afinidad: limpio roto %s -> reparado %s | sucio roto %s -> reparado %s"
                      % (fila.get("limpio_roto(v4)"), fila.get("limpio_reparado"),
                         fila.get("sucio_roto(v4)"), fila.get("sucio_reparado")))
        filas.append(fila)

    texto = "\n".join(lineas)
    print(texto)
    with open(os.path.join(BASE, "confirmar_glue_controles.log"), "w",
              encoding="utf-8") as f:
        f.write(texto + "\n")
    campos = sorted({k for fila in filas for k in fila})
    with open(os.path.join(BASE, "confirmar_glue_controles.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        for fila in filas:
            w.writerow(fila)
    log("guardado: confirmar_glue_controles.csv / .log")
    return 0


if __name__ == "__main__":
    sys.exit(main())
