#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validación de SOD1 con el conjunto limpio de controles y DOS fondos (20 sep 2026).

POR QUE SE REHACE
-----------------
La validación anterior usaba 20 «activos» de los que 16 son la misma serie
congénica (AUDITORIA_POSITIVOS_2026-09-20.md) y daba AUC 0,815 en la caja
Trp32. Con los dos únicos activos independientes el número caía a 0,601 y
LCS-1 no se recuperaba (puesto 152 de 219). Aquí se rehace con:

ACTIVOS, tal y como los deja `controles_sod1_v3.csv`, en dos clases separadas
por el TIPO DE EVIDENCIA (que es lo que cambia la lectura del resultado):

  * independientes ............ LCS-1 (piridazinona), PRG-A01 (cumarina):
                                quimias distintas entre sí y de la serie.
  * serie colapsada ........... CHEMBL2165613, un solo representante del núcleo
                                pirazolona (contar 18 veces la misma química no
                                son 18 positivos).
  * co-cristalizados .......... isoproterenol, adrenalina, dopamina y
                                5-fluorouridina: evidencia ESTRUCTURAL directa
                                (PDB 4A7T, 4A7U, 4A7V, 4A7S). Sirven para
                                validar el SITIO, no para preguntar si el método
                                ordena quimias distintas.

DOS FONDOS
----------
  1. EMPAREJADO EN PROPIEDADES (el clásico): señuelos de la misma librería con
     peso, logP y enlaces rotables parecidos y Tanimoto < 0,35 frente al activo
     (`generar_decoys`, el mismo generador que la validación anterior). Se
     calculan dos versiones: los nuevos, generados para ESTE conjunto limpio, y
     los 199 antiguos, que ya estaban acoplados en la misma caja.
  2. DURO (nuevo): compuestos de la misma librería que llevan grupos que
     QUELAN METALES O SON REDOX-ACTIVOS (catecol, pirogalol, quinona,
     hidroxamato, tiol, sulfonato, fosfonato, hidrazida, tiosemicarbazona).
     En una metaloenzima como SOD1 eso es un fondo mucho más difícil que
     moléculas emparejadas por tamaño: son quimias que de verdad podrían
     pegarse a la superficie.

Todo se acopla con la misma caja Trp32 del cribado (46,5 80,0 73,3; 22 Å),
exhaustividad 8 y el mismo receptor `gpu_dock/SOD1.pdbqt`, para que los números
sean comparables con los anteriores. Los ligandos se preparan con el mismo
pipeline que los señuelos (ETKDGv3 + MMFF + meeko) para no introducir sesgo.

Uso:
    python validar_sod1_v3.py --preparar      # prepara ligandos y fondos
    python validar_sod1_v3.py --acoplar       # acopla (resumible, por tandas)
    python validar_sod1_v3.py --analizar      # métricas y tablas
"""
import argparse
import csv
import glob
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

import validar_senuelos as VS

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(BASE, "validacion_SOD1_v3")
LIGDIR = os.path.join(WORK, "ligands")
OUTDIR = os.path.join(WORK, "out")
CONTROLES = os.path.join(BASE, "controles_sod1_v3.csv")
LIBRERIA = os.path.abspath(os.path.join(BASE, "..", "compounds", "decoys_library.smi"))
RECEPTOR = os.path.abspath(os.path.join(BASE, "..", "gpu_dock", "SOD1.pdbqt"))
CENTRO = (46.5, 80.0, 73.3)
TAMANO = 22
EXHAUSTIVIDAD = 8
ANTIGUO = os.path.join(BASE, "_validacion_SOD1", "validacion_SOD1_trp32.csv")

# grupos que quelan metales o son redox-activos: el fondo duro de una metaloenzima
SMARTS_DURO = [
    ("catecol", "c1cc(O)c(O)cc1"),
    ("pirogalol", "c1c(O)c(O)c(O)cc1"),
    ("quinona", "O=C1[#6]=[#6][#6](=O)[#6]=[#6]1"),
    ("o-aminofenol", "[NX3]c1ccccc1[OH]"),
    ("hidroxamato", "C(=O)[NX3][OH]"),
    ("tiol", "[SX2H]"),
    ("sulfonato", "S(=O)(=O)[OX2H1,OX1-]"),
    ("fosfonato", "P(=O)([OX2H1,OX1-])[OX2H1,OX1-]"),
    ("hidrazida", "C(=O)[NX3][NX3]"),
    ("tiosemicarbazona", "[NX3][NX3]C(=S)"),
    ("catecol-en-ona", "O=C1C=CC(=O)C=C1"),
]

FAMILIAS = {
    "independientes": ["LCS-1", "PRG-A01"],
    "serie": ["CHEMBL2165613"],
    "co-cristalizados": ["isoproterenol", "adrenalina", "dopamina", "5-fluorouridina"],
}
# nombre en controles -> nombre de fichero PDBQT (los tres primeros ya estaban
# preparados con el mismo pipeline en la validacion anterior)
FICHERO = {
    "LCS-1": "ACT_LCS-1",
    "PRG-A01": "ACT_PRG-A01",
    "CHEMBL2165613": "ACT_CHEMBL2165613",
    "isoproterenol": "ACT_isoproterenol",
    "adrenalina": "ACT_adrenalina",
    "dopamina": "ACT_dopamina",
    "5-fluorouridina": "ACT_5-fluorouridina",
}


def log(m):
    print(m, flush=True)


def leer_controles():
    out = []
    with open(CONTROLES, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            out.append(r)
    return out


def clave_del_fondo(nombre, smi):
    """Nombre de fichero del fondo duro: que conste por que entro."""
    return "DECH_%s" % nombre


def preparar_activos(controles):
    os.makedirs(LIGDIR, exist_ok=True)
    hechos = []
    for r in controles:
        nombre, smi = r["ligand"], r["smiles"]
        fichero = FICHERO.get(nombre)
        # los tres primeros ya existen preparados en la validacion anterior:
        # se copian tal cual para mantener el MISMO pipeline
        viejo = os.path.join(BASE, "_validacion_SOD1", "ligands", fichero + ".pdbqt")
        destino = os.path.join(LIGDIR, fichero + ".pdbqt")
        if not os.path.exists(destino):
            if os.path.exists(viejo):
                open(destino, "w", encoding="utf-8").write(
                    open(viejo, encoding="utf-8").read())
            elif not VS.convertir_pdbqt(fichero, smi, LIGDIR):
                log("  FALLO preparando %s" % nombre)
                continue
        hechos.append((nombre, fichero, r.get("quimia", ""), r.get("tipo_ensayo", "")))
        log("  activo %-16s -> %s.pdbqt (%s)" % (nombre, fichero, r.get("quimia", "")))
    return hechos


def preparar_fondos(controles):
    libreria = VS.leer_libreria(LIBRERIA)
    log("libreria de señuelos: %d compuestos" % len(libreria))

    # --- fondo 1: emparejado en propiedades, generado para ESTE conjunto limpio
    usados, n = set(), 0
    lista_emp = []
    activos_para_emparejar = ["LCS-1", "PRG-A01", "CHEMBL2165613",
                              "isoproterenol", "adrenalina", "dopamina",
                              "5-fluorouridina"]
    smis = {r["ligand"]: r["smiles"] for r in controles}
    for nombre in activos_para_emparejar:
        for dn, ds in VS.generar_decoys(smis[nombre], libreria, 25):
            if dn in usados:
                continue
            usados.add(dn)
            if VS.convertir_pdbqt("DECM_" + dn, ds, LIGDIR):
                n += 1
                lista_emp.append(dn)
    log("fondo EMPAREJADO nuevo: %d señuelos" % n)

    # --- fondo 2: duro (quelantes de metales y redox-activos)
    patrones = [(et, Chem.MolFromSmarts(s)) for et, s in SMARTS_DURO]
    mw_objetivo = np.median([Descriptors.MolWt(Chem.MolFromSmiles(r["smiles"]))
                             for r in controles if Chem.MolFromSmiles(r["smiles"])])
    candidatos, motivos = [], {}
    for dn, ds in libreria:
        m = Chem.MolFromSmiles(ds)
        if m is None:
            continue
        etiquetas = [et for et, pat in patrones if pat is not None and m.HasSubstructMatch(pat)]
        if not etiquetas:
            continue
        candidatos.append((abs(Descriptors.MolWt(m) - mw_objetivo), dn, ds, etiquetas))
    candidatos.sort(key=lambda x: x[0])
    n_duro = 0
    lista_duro = []
    for _, dn, ds, etiquetas in candidatos[:200]:
        if dn in usados:
            continue
        usados.add(dn)
        if VS.convertir_pdbqt(clave_del_fondo(dn, ds), ds, LIGDIR):
            motivos[dn] = "+".join(etiquetas)
            lista_duro.append(dn)
            n_duro += 1
    log("fondo DURO (quelantes/redox): %d señuelos de %d candidatos"
        % (n_duro, len(candidatos)))

    # los ficheros de fondos de corridas anteriores (con otra lista de activos)
    # no deben colarse en el analisis: se borran si no estan en la lista de hoy
    validos = set("DECM_" + d for d in lista_emp) | set(
        "DECH_" + d for d in lista_duro)
    obsoletos = [p for p in glob.glob(os.path.join(LIGDIR, "*.pdbqt"))
                 if os.path.basename(p).replace(".pdbqt", "").split("_")[0]
                 in ("DECM", "DECH") and os.path.basename(p).replace(".pdbqt", "")
                 not in validos]
    for p in obsoletos:
        os.remove(p)
    log("  ficheros de fondos obsoletos borrados: %d" % len(obsoletos))
    with open(os.path.join(WORK, "listas_fondos.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["fondo", "nombre", "motivos"])
        for d in sorted(lista_emp):
            w.writerow(["emparejado", d, ""])
        for d in sorted(lista_duro):
            w.writerow(["duro", d, motivos.get(d, "")])
    with open(os.path.join(WORK, "fondo_duro_motivos.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["nombre", "motivos"])
        for k, v in sorted(motivos.items()):
            w.writerow([k, v])


def _dock_one(t):
    vina, receptor, lig, out, cx, cy, cz, size, exh = t
    if os.path.exists(out) and os.path.getsize(out) > 100:
        return out
    subprocess.run([vina, "--receptor", receptor, "--ligand", lig,
                    "--center_x", str(cx), "--center_y", str(cy),
                    "--center_z", str(cz), "--size_x", str(size),
                    "--size_y", str(size), "--size_z", str(size),
                    "--exhaustiveness", str(exh), "--num_modes", "3",
                    "--out", out, "--cpu", "1"],
                   capture_output=True, timeout=3600)
    return out


def acoplar(workers=12, limite=None):
    """Acopla lo que falte. Resumible: se puede llamar tantas veces como haga falta."""
    os.makedirs(OUTDIR, exist_ok=True)
    ligs = sorted(glob.glob(os.path.join(LIGDIR, "*.pdbqt")))
    tareas = []
    for p in ligs:
        nombre = os.path.basename(p).replace(".pdbqt", "")
        out = os.path.join(OUTDIR, nombre + "_out.pdbqt")
        if os.path.exists(out) and os.path.getsize(out) > 100:
            continue
        tareas.append((VS.VINA_CPU, RECEPTOR, p, out, CENTRO[0], CENTRO[1],
                       CENTRO[2], TAMANO, EXHAUSTIVIDAD))
    log("ligandos preparados: %d | pendientes de acoplar: %d" % (len(ligs), len(tareas)))
    if limite:
        tareas = tareas[:limite]
        log("  esta tanda: %d" % len(tareas))
    if tareas:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            for i, _ in enumerate(ex.map(_dock_one, tareas), 1):
                if i % 25 == 0:
                    log("   ... %d/%d" % (i, len(tareas)))
    hechos = len(glob.glob(os.path.join(OUTDIR, "*_out.pdbqt")))
    log("acoplados en total: %d" % hechos)
    return hechos


def afinidades():
    res = {}
    for p in glob.glob(os.path.join(OUTDIR, "*_out.pdbqt")):
        nombre = os.path.basename(p).replace("_out.pdbqt", "")
        a = VS.parse_affinity(p)
        if a is not None:
            res[nombre] = a
    return res


def roc_auc(act, dec):
    if not act or not dec:
        return None
    act, dec = np.array(act, float), np.array(dec, float)
    return sum(float(np.sum(a <= dec)) for a in act) / (len(act) * len(dec))


def ef(act, dec, pct):
    todo = sorted([(s, 1) for s in act] + [(s, 0) for s in dec])
    k = max(1, int(len(todo) * pct / 100.0))
    return (sum(1 for _, l in todo[:k] if l == 1) / max(1, len(act))) / (pct / 100.0)


def puesto(act, dec, valor):
    """Puesto (1 = mejor) del activo en la lista combinada."""
    todo = sorted([(s, 1) for s in act] + [(s, 0) for s in dec])
    for i, (s, l) in enumerate(todo, 1):
        if l == 1 and abs(s - valor) < 1e-9:
            return i
    return None


def contar_pesados(libreria):
    """nombre -> numero de atomos pesados, de la libreria y de los controles."""
    n = {}
    for nombre, smi in libreria:
        m = Chem.MolFromSmiles(smi)
        if m is not None:
            n[nombre] = m.GetNumHeavyAtoms()
    for r in leer_controles():
        m = Chem.MolFromSmiles(r["smiles"])
        if m is not None:
            n[r["ligand"]] = m.GetNumHeavyAtoms()
    return n


def normalizar_por_tamano(scores, tamanos, fondo):
    """Quita la parte del score que explica el numero de atomos pesados.

    Es el sesgo de tamaño que la auditoria ya habia medido: en una caja pequeña
    las moleculas grandes puntuan mejor por su sola superficie de contacto. Se
    ajusta afinidad = a + b*N sobre el FONDO (nunca sobre los activos) y se
    trabaja con el residuo: lo que queda del score cuando ya se ha descontado el
    tamaño. Si el AUC tampoco sube con el residuo, el problema no es el tamaño,
    es que la puntuacion no tiene señal.
    """
    xs = np.array([tamanos[k] for k in fondo if tamanos.get(k)], dtype=float)
    ys = np.array([scores[k] for k in fondo if tamanos.get(k)], dtype=float)
    if len(xs) < 10:
        return None, None, None
    b, a = np.polyfit(xs, ys, 1)
    resid = {k: scores[k] - (a + b * tamanos[k])
             for k in scores if tamanos.get(k)}
    return resid, a, b


def analizar():
    aff = afinidades()
    # los 20 activos y 199 señuelos antiguos ya estan acoplados con la misma caja
    antiguo = {}
    with open(ANTIGUO, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r["affinity"]:
                antiguo[r["ligand"]] = float(r["affinity"])

    controles = leer_controles()
    resultados = []
    fondo_emp_nuevo = [v for k, v in aff.items() if k.startswith("DECM_")]
    fondo_duro = [v for k, v in aff.items() if k.startswith("DECH_")]
    fondo_emp_antiguo = [v for k, v in antiguo.items() if k.startswith("DEC_")]
    fondo_todo = fondo_emp_nuevo + fondo_duro + fondo_emp_antiguo

    log("fondos: emparejado nuevo %d | duro %d | emparejado antiguo %d"
        % (len(fondo_emp_nuevo), len(fondo_duro), len(fondo_emp_antiguo)))

    # --- normalizacion por tamaño
    tamanos = contar_pesados(VS.leer_libreria(LIBRERIA))

    def tam(k):
        """Tamaño de un ligando: los ficheros llevan prefijo DECM_/DECH_/DEC_/ACT_."""
        base = k.split("_", 1)[1] if k.startswith(
            ("DEC_", "DECM_", "DECH_", "ACT_")) else k
        return tamanos.get(base) or tamanos.get(k)

    scores_todo = dict(antiguo)
    scores_todo.update(aff)
    tamanos_lig = {k: tam(k) for k in scores_todo if tam(k)}
    keys_fondo = [k for k in scores_todo if k.startswith(("DEC_", "DECM_", "DECH_"))]
    resid, a_ajuste, b_ajuste = normalizar_por_tamano(scores_todo, tamanos_lig, keys_fondo)
    if b_ajuste is not None:
        log("sesgo de tamaño medido en el fondo: afinidad = %.2f %+.3f * (atomos pesados)"
            % (a_ajuste, b_ajuste))
        log("   (un carbono mas vale %.2f kcal/mol en este bolsillo)" % b_ajuste)
    fondo_emp_nuevo_r = [resid[k] for k in aff if k.startswith("DECM_") and k in resid]
    fondo_duro_r = [resid[k] for k in aff if k.startswith("DECH_") and k in resid]
    fondo_emp_antiguo_r = [resid[k] for k in scores_todo
                           if k.startswith("DEC_") and k in resid]
    fondo_todo_r = fondo_emp_nuevo_r + fondo_duro_r + fondo_emp_antiguo_r


    familias = dict(FAMILIAS)
    familias["todos"] = FAMILIAS["independientes"] + FAMILIAS["serie"] + \
        FAMILIAS["co-cristalizados"]

    filas = []
    detalle = []
    for nombre in familias["todos"]:
        fichero = FICHERO[nombre]
        valor = aff.get(fichero, antiguo.get(fichero))
        detalle.append({"ligando": nombre, "fichero": fichero, "afinidad": valor})

    for fam, miembros in familias.items():
        act = [aff.get(FICHERO[m], antiguo.get(FICHERO[m])) for m in miembros]
        act = [a for a in act if a is not None]
        act_r = [resid[FICHERO[m]] for m in miembros
                 if resid and FICHERO[m] in resid]
        for etiqueta, fondo, fondo_r in (
                ("emparejado nuevo", fondo_emp_nuevo, fondo_emp_nuevo_r),
                ("emparejado antiguo", fondo_emp_antiguo, fondo_emp_antiguo_r),
                ("duro", fondo_duro, fondo_duro_r),
                ("los tres juntos", fondo_todo, fondo_todo_r)):
            auc = roc_auc(act, fondo)
            auc_r = roc_auc(act_r, fondo_r) if act_r and fondo_r else None
            filas.append({
                "familia": fam, "n_activos": len(act), "fondo": etiqueta,
                "n_fondo": len(fondo),
                "media_activos": round(float(np.mean(act)), 2) if act else None,
                "media_fondo": round(float(np.mean(fondo)), 2) if fondo else None,
                "AUC": round(auc, 3) if auc is not None else None,
                "AUC_sin_tamano": round(auc_r, 3) if auc_r is not None else None,
                "EF1%": round(ef(act, fondo, 1.0), 2) if act and fondo else None,
                "EF5%": round(ef(act, fondo, 5.0), 2) if act and fondo else None,
            })
    # puesto de cada activo dentro de la lista combinada de los tres fondos
    for d in detalle:
        if d["afinidad"] is None:
            continue
        d["puesto_frente_a_los_tres_fondos"] = puesto([d["afinidad"]], fondo_todo,
                                                      d["afinidad"])
        d["n_fondo"] = len(fondo_todo)
        d["cuantil"] = round(100.0 * d["puesto_frente_a_los_tres_fondos"]
                             / (len(fondo_todo) + 1), 1)

    with open(os.path.join(WORK, "analisis_sod1_v3.csv"), "w", newline="",
              encoding="utf-8") as f:
        campos = ["familia", "n_activos", "fondo", "n_fondo", "media_activos",
                  "media_fondo", "AUC", "AUC_sin_tamano", "EF1%", "EF5%"]
        w = csv.DictWriter(f, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        for r in filas:
            w.writerow(r)
    with open(os.path.join(WORK, "detalle_activos_sod1_v3.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["ligando", "fichero", "afinidad",
                                          "puesto_frente_a_los_tres_fondos",
                                          "n_fondo", "cuantil"],
                           extrasaction="ignore")
        w.writeheader()
        for d in detalle:
            w.writerow(d)

    log("\n%-18s %-20s %-6s %-19s %-8s %-14s %-6s %-6s"
        % ("familia", "fondo", "n_act", "n_fondo (media)", "AUC",
           "AUC sin tamaño", "EF1%", "EF5%"))
    for r in filas:
        log("%-18s %-20s %-6s %-19s %-8s %-14s %-6s %-6s"
            % (r["familia"], r["fondo"], r["n_activos"],
               "%d (%.2f kcal/mol)" % (r["n_fondo"], r["media_fondo"] or 0),
               r["AUC"], r["AUC_sin_tamano"], r["EF1%"], r["EF5%"]))
    log("\nafinidades y puesto de los controles (frente a los tres fondos juntos, N=%d):"
        % len(fondo_todo))
    for d in sorted(detalle, key=lambda x: (x["afinidad"] is None, x["afinidad"])):
        log("   %-16s %-8s puesto %s de %d (mejor que el %.1f%% del fondo)"
            % (d["ligando"],
               "%.2f" % d["afinidad"] if d["afinidad"] is not None else "n/d",
               d.get("puesto_frente_a_los_tres_fondos", "n/d"), len(fondo_todo),
               100.0 - (d.get("cuantil") or 100.0)))
    return 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preparar", action="store_true")
    ap.add_argument("--acoplar", action="store_true")
    ap.add_argument("--analizar", action="store_true")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--limite", type=int, default=None,
                    help="maximo de ligandos en esta tanda")
    args = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    controles = leer_controles()
    if args.preparar or not any([args.preparar, args.acoplar, args.analizar]):
        log("=== preparando activos ===")
        preparar_activos(controles)
        log("=== preparando fondos ===")
        preparar_fondos(controles)
    if args.acoplar:
        log("=== acoplando en la caja Trp32 (%.1f, %.1f, %.1f; %d A; ex %d) ==="
            % (CENTRO[0], CENTRO[1], CENTRO[2], TAMANO, EXHAUSTIVIDAD))
        acoplar(workers=args.workers, limite=args.limite)
    if args.analizar:
        log("=== analizando ===")
        analizar()
    return 0


if __name__ == "__main__":
    sys.exit(main())
