#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""CONTROL: se vuelve a acoplar la naftalensulfonamida de hCA2 (6T4P, MHK).

POR QUE ESTE CONTROL
--------------------
Es la prueba de fuego de hCA2. Con el receptor preparado y el zinc dentro
(`preparar_receptor_hca2.py`), la pregunta es una sola: cuando Vina acopla la
sulfonamida en esta caja, ¿la pone donde la puso el cristal? Si la pone, el
receptor y la caja estan bien, y el AUC que venga despues se puede medir con
confianza. Si no la pone, el problema es el receptor o la caja, y hay que
arreglarlo ANTES de generar un solo señuelo.

Que no sea el control de CDK2
----------------------------
CDK2 resulto ser un control MAL ELEGIDO, no una prueba de que el embudo fallara: su
ligando venia modelado dos veces al 50 %, con B ~47 y un brazo suelto a 6-9 A. Que el
motor no encontrara la geometria no decia nada del motor.

6T4P no tiene ninguno de esos problemas, y esta medido:
  - una sola copia, sin conformaciones alternativas, B medio 13;
  - el N de la sulfonamida a 2,00 A del zinc, que es la coordinacion real de hCA2;
  - los tres histidinas del zinc (His94, His96, His119) a 1,99-2,04 A: es la
    geometria verdadera, leida del propio cristal;
  - ningun polar del farmaco sin pareja en la proteina;
  - resolucion 1,75 A y R-free 0,199;
  - 14 atomos pesados y 1 enlace rotatorio;
  - y, lo mas importante, es un compuesto con POTENCIA MEDIDA (CHEMBL6648, Ki 6,88),
    de modo que el control y el AUC hablan de la misma molecula.

QUE MIDE
--------
Lo mismo que el control de CDK2 y que el de andr, con el MISMO codigo importado, no
reescrito. Si la forma de medir cambia, cambia en un sitio y los tres controles miden
igual.

  - `bolsillo`: RMSD de la pose al cristal SIN alinear. **Es el que decide**, con el
    liston de 2,0 A, porque responde a "esta en el mismo sitio del bolsillo?".
  - `piso`: el mejor RMSD que podria dar cualquier emparejamiento que respete el
    elemento. Si el piso supera el liston, el veredicto no depende ni del orden de
    atomos ni de la simetria.
  - `alineado`: GetBestRMS de RDKit, que superpone antes de medir. Es un
    diagnostico de FORMA, no un criterio de sitio.

Y un criterion extra que solo tiene sentido en hCA2: **la distancia del nitrogeno
acoplado al zinc**. El AUC de hCA2 va a Premiar "tengo una sulfonamida pegada al
metal", asi que ademas del RMSD hay que comprobar que la pose acoplada de verdad se
coordina al Zn, y no se queda a 5 A Taverna en el mismo bolsillo. Sin esta medida
un RMSD de 1,5 A podria ser una pose que esta en el sitio pero sin tocar el metal.

LA CORRESPONDENCIA DE ATOMOS NO SE SUPONE
-----------------------------------------
Se deduce con la asignacion optima exigiendo mismo elemento sobre el PDBQT de
entrada, igual que en CDK2 y andr. Si esa asignacion no da ~0,00 A, el control se
para: significaria que lo que se acopla no es la pose del cristal.

Salida: `_control_hca2/<etiqueta>/`, `barrido_redocking_hca2.txt`, `CONTROL_HCA2_RMD.md`

Uso:
    python control_redocking_hca2.py --caja 24 --depth 20
    python control_redocking_hca2.py --barrido
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(BASE))

# El codigo de medicion es el de CDK2, que es el que se uso y se depuro en andr. Se
# importa de verdad (no se copia) para que los tres controles midan igual.
sys.path.insert(0, os.path.join(RAIZ, "analysis", "calibracion_dude_cdk2"))
import control_redocking_cdk2 as base  # noqa: E402

CAJA = os.path.join(BASE, "caja_hca2.json")
RECEPTOR_PDB = os.path.join(BASE, "receptor_hca2.pdb")
LIG_PDB = os.path.join(BASE, "ligando_cristal.pdb")
CARPETA = os.path.join(BASE, "_control_hca2")
RESUMEN = os.path.join(BASE, "barrido_redocking_hca2.txt")
INFORME = os.path.join(BASE, "CONTROL_HCA2_RMD.md")

PDB_ORIG = "6T4P"
NOMBRE_LIG = "naphthalene-1-sulfonamide (CHEMBL6648, Ki p6.88)"
CORTE_ZN_ACOPLADA = 2.6      # N-Sulfonamida -> Zn de una pose bien puesta


def log(m):
    print(m, flush=True)


def zn_del_receptor_pdbqt(ruta):
    """Coordenadas del Zn que hay en el receptor YA PREPARADO, o None."""
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if not l.startswith(("ATOM", "HETATM")):
            continue
        # El tipo sale de la ULTIMA columna del PDBQT y meeko lo escribe como
        # "Zn" (mayuscula inicial), no "ZN". Buscarlo solo en mayusculas devolvia
        # None y hacia que la medida de la union al zinc saliera como no medida.
        # Es el fallo de verdad aqui: el RMSD pasaba y el "NO COORDINA" que salia
        # era del script, no de la pose.
        if l.rsplit(None, 1)[-1].strip().upper() in ("ZN", "DZN"):
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))
    return None


def distancia_n_al_zinc(modelos, zn, tipos_fichero, mol_ref, ref, perm):
    """Distancia del N del sulfonamida al Zn en cada pose, y de que pose va.

    Devuelve la lista [(indice_pose, distancia)]. Los N se eligen por tipo de
    elemento en el PDBQT de salida, emparejados al mol2 con la misma permutacion
    que se uso para el RMSD, para no estar comparando coordenadas de atomos
    distintos.
    """
    if zn is None:
        return []
    idx_n = [i for i, t in enumerate(tipos_fichero) if t.upper() == "N"]
    if not idx_n:
        return []
    out = []
    for k, m in enumerate(modelos, 1):
        pos = np.asarray(m["pos"])
        d = min(float(np.linalg.norm(np.array(zn) - pos[i])) for i in idx_n)
        out.append((k, d))
    return out


def variante(receptor, mol, ref, centro, tamano, depth, etiqueta, zn):
    carpeta = os.path.join(CARPETA, etiqueta)
    import shutil
    shutil.rmtree(carpeta, ignore_errors=True)
    base.escribir_ligando(mol, os.path.join(carpeta, "ligands"))
    tipos, coords = base.atomos_pesados_pdbqt(
        os.path.join(carpeta, "ligands", "FAP_cristal.pdbqt"))
    perm, rmsd_entrada = base.emparejar_orden(tipos, coords, mol, ref)
    if rmsd_entrada > base.TOL_ENTRADA:
        raise SystemExit("lo que se acoplo no es la pose del cristal (%.3f A):"
                         " la correspondencia no vale y el control no"
                         " significaria nada" % rmsd_entrada)
    log("  entrada: %d pesados, correspondencia con el cristal a %.4f A"
        % (len(coords), rmsd_entrada))
    modelos, ruta_out = base.acoplar(carpeta, receptor, centro, tamano, depth)
    if not modelos:
        log("  el motor no devolvio ninguna pose")
        with open(RESUMEN, "a", encoding="utf-8") as f:
            f.write("caja %2d A | depth %2d | SIN POSES\n" % (tamano, depth))
        return None
    filas = base.medir(mol, ref, modelos, perm, tipos)
    validas = [f for f in filas if f[2] is not None]
    mejor = min(validas, key=lambda f: f[2]) if validas else None
    if mejor is None:
        log("  sin poses comparables")
        with open(RESUMEN, "a", encoding="utf-8") as f:
            f.write("caja %2d A | depth %2d | sin poses comparables\n"
                    % (tamano, depth))
        return None

    # La medida propia de hCA2: ¿la pose acoplada toca el zinc?
    # El PDBQT de AutoDock escribe N a un nitrogeno aromatico y NA a uno aceptor, y
    # los dos son nitrogeno. Para la union al zinc hay que mirar los dos: si solo
    # se mirara "N" y el nitrogeno de la sulfonamida saliera como "NA", no se
    # mediria nada. Se comparan por elemento, que es lo que el PDBQT y el mol2
    # comparten.
    idx_n = [i for i, t in enumerate(tipos)
             if base.GRUPO.get(t.upper(), t.upper()) == "N"]
    dzn_pose = None
    if zn is not None and idx_n:
        pos = np.array(modelos[mejor[0] - 1]["pos"])
        dzn_pose = min(float(np.linalg.norm(np.array(zn) - pos[i]))
                       for i in idx_n)
    dzn_cristal = None
    if zn is not None and idx_n:
        dzn_cristal = min(float(np.linalg.norm(np.array(zn) - ref[perm[i]]))
                          for i in idx_n)

    log("  pose %d: %+.2f kcal/mol | bolsillo %.2f A | piso %.2f A | alineado"
        " %.2f A | solape %.2f A" % (mejor[0], mejor[1], mejor[2], mejor[3],
                                     mejor[4], mejor[5]))
    log("  union al Zn: pose acoplada %.2f A | pose del cristal %.2f A %s"
        % (dzn_pose if dzn_pose is not None else -1,
           dzn_cristal if dzn_cristal is not None else -1,
           "COORDINA" if dzn_pose is not None and dzn_pose <= CORTE_ZN_ACOPLADA
           else "NO COORDINA"))
    if mejor[3] is not None and mejor[3] > base.LISTON_A:
        log("  ni con el mejor emparejamiento de atomos posible baja de"
            " %.2f A: el veredicto no depende del orden" % mejor[3])

    pasa = mejor[2] < base.LISTON_A
    if dzn_pose is not None and dzn_pose > CORTE_ZN_ACOPLADA:
        pasa = False
        log("  aunque el RMSD pase, la pose NO se coordina al zinc: para hCA2 eso")
        log("  no cuenta como control valido")
    linea = ("caja %2d A | depth %2d | mejor pose %d | %+.2f kcal/mol | "
             "bolsillo %.2f A | piso %.2f A | alineado %.2f A | N-Zn pose %.2f A"
             " | N-Zn cristal %.2f A | %s"
             % (tamano, depth, mejor[0], mejor[1], mejor[2], mejor[3], mejor[4],
                dzn_pose if dzn_pose is not None else -1,
                dzn_cristal if dzn_cristal is not None else -1,
                "PASA" if pasa else "NO PASA"))
    log(linea)
    with open(RESUMEN, "a", encoding="utf-8") as f:
        f.write(linea + "\n")
    return {"mejor": mejor, "dzn_pose": dzn_pose, "dzn_cristal": dzn_cristal,
            "pasa": pasa, "etiqueta": etiqueta, "caja": tamano, "depth": depth,
            "rmsd_entrada": rmsd_entrada}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--caja", type=int, default=0)
    ap.add_argument("--depth", type=int, default=20)
    ap.add_argument("--barrido", action="store_true")
    args = ap.parse_args()

    log("CONTROL DE REDOCKING DE hCA2: %s (%s)" % (PDB_ORIG, NOMBRE_LIG))
    log("   %s" % time.strftime("%Y-%m-%d %H:%M"))
    if not os.path.exists(CAJA):
        raise SystemExit("falta caja_hca2.json: ejecuta antes "
                         "preparar_receptor_hca2.py")
    caja = json.load(open(CAJA, encoding="utf-8"))
    receptor = os.path.join(BASE, caja["receptor_pdbqt"])

    # El zinc del receptor YA PREPARADO. Si no esta, el control no se puede hacer:
    # hCA2 sin Zn es un receptor falso.
    zn = zn_del_receptor_pdbqt(receptor)
    if zn is None:
        log("")
        log("ATENCION: el receptor preparado NO tiene Zn. Este control no se puede")
        log("           hacer con el receptor sin el metal, y no se va a hacer con")
        log("           otro receptor a proposito: seria medir otra cosa.")
        return 1
    log("   Zn en el receptor preparado: (%.2f, %.2f, %.2f)" % zn)

    # El farmaco del cristal, leido del PDB de referencia (no de un mol2 de nadie)
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    mol = Chem.MolFromMolFile(os.path.join(BASE, "_lig_tmp.mol"), removeHs=True,
                              sanitize=False)
    if mol is None:
        mol = Chem.MolFromPDBFile(LIG_PDB, removeHs=True, sanitize=False)
    if mol is None:
        raise SystemExit("no se pudo leer el farmaco del cristal")
    mol.UpdatePropertyCache(strict=False)
    Chem.SanitizeMol(mol, Chem.SanitizeFlags.SANITIZE_ALL
                     ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE
                     ^ Chem.SanitizeFlags.SANITIZE_SETAROMATICITY)
    conf = mol.GetConformer()
    ref = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y,
                     conf.GetAtomPosition(i).z] for i in range(mol.GetNumAtoms())])
    log("   farmaco del cristal: %d atomos pesados" % mol.GetNumAtoms())

    # Los grupos se comparan por ELEMENTO: el PDBQT de AutoDock llama C a un
    # carbono aromatico y N a un nitrogeno aceptor. El Zn se annade al grupo C para
    # que el emparejamiento de atomos no lo cuente como si fuera otra cosa.
    base.GRUPO["ZN"] = "ZN"
    base.GRUPO["DZN"] = "ZN"

    centro = np.array(caja["centro_caja"], dtype=float)
    os.makedirs(CARPETA, exist_ok=True)
    with open(RESUMEN, "a", encoding="utf-8") as f:
        f.write("\n=== %s | %s (%s) ===\n" % (PDB_ORIG, NOMBRE_LIG,
                                             time.strftime("%Y-%m-%d %H:%M")))

    res = []
    if args.barrido:
        for tam in (20, 22, 24):
            for depth in (20, 32):
                r = variante(receptor, mol, ref, centro, tam, depth,
                             "caja%d_depth%d" % (tam, depth), zn)
                if r:
                    res.append(r)
    else:
        tam = args.caja or caja["tamano_caja"]
        r = variante(receptor, mol, ref, centro, tam, args.depth,
                     "caja%d_depth%d" % (tam, args.depth), zn)
        if r:
            res.append(r)

    # ---------------- veredicto ----------------
    log("")
    log("=" * 74)
    log("VEREDICTO DEL CONTROL")
    log("=" * 74)
    if not res:
        log("no se pudo hacer el control: no hay poses que medir")
        return 1
    pasan = [r for r in res if r["pasa"]]
    mejor = min(res, key=lambda r: r["mejor"][2])
    m = mejor["mejor"]
    log("   mejor RMSD de bolsillo: %.2f A (liston %.1f A) | piso %.2f A |"
        " alineado %.2f A" % (m[2], base.LISTON_A, m[3], m[4]))
    if mejor["dzn_pose"] is not None:
        log("   union al Zn de la pose: %.2f A (cristal %.2f A)"
            % (mejor["dzn_pose"], mejor["dzn_cristal"]))
    log("   %d de %d variantes pasan" % (len(pasan), len(res)))
    log("")
    if pasan:
        log("   EL CONTROL PASA. El receptor y la caja de hCA2 estan bien: Vina")
        log("   devuelve la pose del cristal. A partir de aqui se puede construir")
        log("   el banco de senuelos y leer el AUC, y ese AUC sera del embudo.")
    else:
        log("   EL CONTROL NO PASA. El problema es el RECEPTOR O LA CAJA, no la")
        log("   serie de activos: la sulfonamida no vuelve a su sitio. Se dice")
        log("   ANTES de generar un solo senuelo, que es lo que se hizo con CDK2.")

    md = []
    md.append("# Control de redocking de hCA2: %s\n" % PDB_ORIG)
    md.append("La naftalensulfonamida del cristal, acoplada de nuevo en la misma "
              "caja, con el zinc dentro del receptor.\n")
    md.append("## Resultado\n")
    md.append("| variante | RMSD de bolsillo | piso | alineado | N-Zn pose |"
              " N-Zn cristal | veredicto |")
    md.append("|---|---|---|---|---|---|---|")
    for r in res:
        m = r["mejor"]
        md.append("| caja %d A, depth %d | **%.2f A** | %.2f A | %.2f A |"
                  " %.2f A | %.2f A | %s |"
                  % (r["caja"], r["depth"], m[2], m[3], m[4],
                     r["dzn_pose"] if r["dzn_pose"] is not None else -1,
                     r["dzn_cristal"] if r["dzn_cristal"] is not None else -1,
                     "PASA" if r["pasa"] else "no pasa"))
    md.append("")
    md.append("Liston de bolsillo: %.1f A. Ademas, para hCA2 hace falta que la "
              "pose acoplada se coordine al zinc a menos de %.1f A: sin eso, un "
              "buen RMSD puede ser una pose en el sitio correcto pero sin tocar "
              "el metal, que es justo lo que el AUC de hCA2 no debe premiar.\n"
              % (base.LISTON_A, CORTE_ZN_ACOPLADA))
    md.append("## Que significa\n")
    if pasan:
        md.append("El control pasa. El receptor con su zinc y la caja de 24 A "
                  "reproducen la pose del cristal, asi que el AUC que venga midera "
                  "el embudo y no una falta de preparacion. Este es el primer "
                  "control de redocking del proyecto que sirve para algo, despues "
                  "del de `andr`.\n")
    else:
        md.append("El control no pasa. Eso no dice que el cribado falle: dice que "
                  "el receptor o la caja no estan bien. Hay que arreglarlo antes de "
                  "generar un solo senuelo, igual que se hizo con CDK2.\n")
    md.append("## Ficheros\n")
    md.append("- `analysis/preparar_receptor_hca2.py` — receptor con el zinc.\n"
              "- `analysis/calibracion_hca2/control_redocking_hca2.py` — este "
              "redocking.\n"
              "- `analysis/calibracion_hca2/barrido_redocking_hca2.txt` — todas "
              "las variantes.\n")
    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("\n".join(md))
    log("")
    log("informe: %s" % os.path.basename(INFORME))
    return 0 if pasan else 2


if __name__ == "__main__":
    raise SystemExit(main())
