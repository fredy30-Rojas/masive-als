#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""preparar_ligando.py — LA receta canonica de SMILES a PDBQT del proyecto.

POR QUE EXISTE
--------------
Este fichero nace del fallo mas caro del proyecto: los primeros puestos de TODAS
las validaciones (TDP-43, FUS, SOD1) eran pseudo-atomos de pegamento.

Meeko, con los ajustes por defecto, no cierra los anillos de 7 eslabones en
adelante: los abre en dos ramas y pega los extremos con pseudo-atomos `CG0`/`G0`
(elemento `Du`), DUPLICANDO dos atomos por anillo. El fichero resultante es
veneno silencioso:

  * Vina reporta `score = inter + intra - unbound`. Con el anillo abierto la
    energia interna del ligando LIBRE sube a +7...+16 kcal/mol en vez de ~0, y al
    restarla el compuesto recibe ese regalo como afinidad. Medido: cefarantina
    -10,72 -> -7,47; los 14 primeros señuelos de TDP-43 perdieron 3,2-6,1
    kcal/mol cada uno al repararlos.
  * El ligando tiene 2 atomos de mas, asi que el rescoring no puede emparejar la
    pose con el SMILES ("n atomos no coincide: SMILES 28 vs pose 30").
  * Vina-GPU rechaza los tipos `CG0`, y `sanear_tipos` los traducia a carbono:
    el ligando se acoplaba con dos carbonos fantasma y nadie se enteraba.

El proyecto ya habia arreglado la LIBRERIA (`gpu_dock/reparar_libreria_glue.py`,
3.803 afectados) y los conjuntos de validacion despues
(`reparar_ligandos_glue.py` y compañia), pero cada script tenia su propia copia
de la receta y la de `analysis/validar_senuelos.py` seguia con el ajuste
defectuoso. Este modulo es la unica copia: si hay que cambiar la receta, se
cambia aqui y todos los scripts la heredan.

QUE HACE DISTINTO
-----------------
1. `MoleculePreparation(rigid_macrocycles=True)`: el anillo grande se acopla
   rigido en vez de romperse.
2. Incrustacion robusta: cinco semillas en vez de una. Con una sola semilla,
   `EmbedMolecule` devuelve distinto de cero y el ligando se cae del conjunto en
   silencio.
3. **Comprobacion de atomos contra el SMILES** y **cero pseudo-atomos**, y si
   falla devuelve el motivo en vez de escribir un fichero malo. Es la parte que
   faltaba en todos los scripts.

Uso:
    from preparar_ligando import construir, escribir

    txt, n, motivo = construir("C1CCCCCC1...")     # texto PDBQT o None + motivo
    ruta = escribir("ACT_x", smi, "carpeta/ligands")
"""
import os
import sys

from meeko import MoleculePreparation, PDBQTWriterLegacy
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

# tipos que NO son atomos de la molecula (hidrogenos y pegamento)
SIN_ATOMO_REAL = ("H", "HD", "HS", "D", "DD")
TIPOS_GLUE = ("G0", "CG0")

# varias semillas: con una sola, algunos macrociclos no incrustan
SEMILLAS = (42, 2026, 777, 1234, 9999, 20260921)


def contar_atomos(txt):
    """(numero de atomos pesados, cuantos son pseudo-atomos glue) de un PDBQT."""
    n = glue = 0
    for l in txt.splitlines():
        if not l.startswith(("ATOM", "HETATM")):
            continue
        t = l.rsplit(None, 1)[-1].strip()
        if t in SIN_ATOMO_REAL:
            continue
        n += 1
        if t in TIPOS_GLUE:
            glue += 1
    return n, glue


def contar_fichero(path):
    return contar_atomos(open(path, encoding="utf-8", errors="ignore").read())


def smi_del_fichero(path):
    """SMILES que llevan las poses y los ligandos de meeko en `REMARK SMILES`."""
    for l in open(path, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK SMILES "):
            return l[len("REMARK SMILES "):].strip()
        if l.startswith(("ATOM", "HETATM")):
            break
    return None


def mayor_fragmento(mol):
    """El fragmento organico mas grande de una molecula con sal.

    Media docena de controles llegan con su sal en el SMILES (`CC(C)NCC(O)...Cl`
    para el isoproterenol, el tartrato de la epinefrina, y el clorhidrato de la
    dopamina). Meeko rechaza cualquier molecula de mas de un fragmento ("RDKit
    molecule has 2 fragments. Must have 1"), asi que sin esto esos ligandos se
    caen del conjunto en silencio. Acoplar el cation es lo correcto: es la
    especie que se une.
    """
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    if len(frags) <= 1:
        return mol
    return max(frags, key=lambda m: m.GetNumHeavyAtoms())


def construir(smi, quitar_sales=False):
    """SMILES -> (texto PDBQT, atomos pesados, motivo).

    Devuelve (None, 0, motivo) si no se puede preparar BIEN. Nunca devuelve texto
    con pseudo-atomos ni con un numero de atomos distinto al del SMILES.

    `quitar_sales=True` se queda con el fragmento organico mayor: necesario para
    los SMILES que vienen como sal (catecolaminas). Por defecto NO se toca la
    molecula, para que nadie cambie una quimia sin decirlo.
    """
    mol = Chem.MolFromSmiles(smi) if smi else None
    if mol is None:
        return None, 0, "sin SMILES valido"
    if quitar_sales:
        mol = mayor_fragmento(mol)
    # GetNumHeavyAtoms y no GetNumAtoms: los analogos deuterados (ChEMBL esta
    # lleno) llevan los [2H] EXPLICITOS en el SMILES y GetNumAtoms los cuenta.
    # Un deuterio es un hidrogeno, y en el PDBQT va como H: contarlos daba un
    # falso positivo de "atomos no coincide" en DECH_CHEMBL3137326 (29 vs 23).
    n_smi = mol.GetNumHeavyAtoms()

    mh = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    params.maxIterations = 1000
    for semilla in SEMILLAS:
        params.randomSeed = semilla
        if AllChem.EmbedMolecule(mh, params) == 0:
            break
    else:
        return None, 0, "no se pudo generar la geometria 3D"
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=200)
    except Exception:
        pass

    try:
        setups = MoleculePreparation(rigid_macrocycles=True).prepare(mh)
    except Exception as e:
        return None, 0, "meeko fallo: %s" % e
    if not setups:
        return None, 0, "meeko no encontro parametros"
    txt, bien, err = PDBQTWriterLegacy.write_string(setups[0])
    if not bien:
        return None, 0, "escritura PDBQT: %s" % err

    n, glue = contar_atomos(txt)
    if glue:
        return None, 0, "sigue con %d pseudo-atomos" % glue
    if n != n_smi:
        return None, 0, "atomos SMILES %d vs preparado %d" % (n_smi, n)
    return txt, n, ""


def escribir(nombre, smi, outdir, forzar=False, quitar_sales=False):
    """Escribe `outdir/<nombre>.pdbqt`. Devuelve la ruta o None (+ mensaje)."""
    os.makedirs(outdir, exist_ok=True)
    outp = os.path.join(outdir, nombre + ".pdbqt")
    if not forzar and os.path.exists(outp) and os.path.getsize(outp) > 100:
        return outp
    txt, _, motivo = construir(smi, quitar_sales=quitar_sales)
    if txt is None:
        print("   [preparar_ligando] %s NO preparado: %s" % (nombre, motivo))
        return None
    with open(outp, "w", encoding="utf-8") as f:
        f.write(txt)
    return outp


def verificar_carpeta(carpeta, verbose=False):
    """Cuenta los ficheros con pseudo-atomos de una carpeta. Para auditar.

    Devuelve la lista de nombres afectados. Se usa como red de seguridad: si algo
    vuelve a prepararse mal, esta funcion lo dice.
    """
    import glob
    rotos = []
    for p in sorted(glob.glob(os.path.join(carpeta, "*.pdbqt"))):
        n, glue = contar_fichero(p)
        if glue:
            rotos.append((os.path.basename(p), n, glue))
            if verbose:
                print("   %s: %d atomos, %d pseudo-atomos" % (rotos[-1]))
    return rotos


def _autoprueba():
    """Prueba de regresion con las moleculas que de verdad fallaron.

    Recorre los `_roto_glue/` que dejaron las reparaciones (los ficheros malos,
    que llevan su SMILES en `REMARK SMILES`) y comprueba que la receta canonica
    los vuelve a preparar con cero pseudo-atomos y el numero de atomos del
    SMILES. Son los 32 casos reales, no moleculas inventadas.
    """
    fallos = 0

    # 1) anillos de 7 a 9 eslabones: el caso minimo que rompia meeko
    for n_esl in (7, 8, 9):
        smi = "C1" + "C" * (n_esl - 1) + "1"
        txt, n, motivo = construir(smi)
        if txt is None:
            print("FALLO  anillo_%d   %s" % (n_esl, motivo))
            fallos += 1
        else:
            print("OK     anillo_%-3d %2d atomos, 0 pseudo-atomos" % (n_esl, n))

    # 2) los casos reales del fallo
    import glob
    base = os.path.dirname(os.path.abspath(__file__))
    carpetas = sorted(glob.glob(os.path.join(
        base, "analysis", "*", "ligands", "_roto_glue")))
    if not carpetas:
        print("\n(no hay carpetas _roto_glue: no se pueden probar los casos reales)")
    total_ok = 0
    for carpeta in carpetas:
        print("\n--- %s" % carpeta.replace(base + os.sep, ""))
        for p in sorted(glob.glob(os.path.join(carpeta, "*.pdbqt"))):
            nombre = os.path.basename(p).replace(".pdbqt", "")
            smi = smi_del_fichero(p)
            n_viejo, glue_viejo = contar_fichero(p)
            txt, n, motivo = construir(smi)
            if txt is None:
                print("FALLO  %-26s %s" % (nombre, motivo))
                fallos += 1
                continue
            print("OK     %-26s roto %2d atomos/%d glue -> %2d atomos/0 glue"
                  % (nombre, n_viejo, glue_viejo, n))
            total_ok += 1
    print("\ncasos reales probados: %d | fallos: %d" % (total_ok, fallos))
    return fallos


if __name__ == "__main__":
    sys.exit(1 if _autoprueba() else 0)
