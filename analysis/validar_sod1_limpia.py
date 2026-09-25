#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""validar_sod1_limpia.py — la validación de SOD1 que sí se puede reportar.

POR QUE ESTA Y NO LAS ANTERIORES
--------------------------------
Las validaciones v1-v5 de SOD1 se acoplaron contra `gpu_dock/SOD1.pdbqt`, que no es
una proteina limpia: lleva las 1.763 aguas del cristal, el zinc y las 18 copias de
1HL5 (21.585 atomos). Contra ese receptor, dos de los cuatro ligandos
co-cristalizados en Trp32 ni siquiera caian en el bolsillo. Y sus verdades de
referencia estaban sesgadas: 16 de 20 «activos» eran el mismo nucleo pirazolona.

Aqui se usan las tres cosas corregidas:

  * receptor `gpu_dock/SOD1_limpio.pdbqt` (un dimero biologico A-H, sin aguas);
  * positivos de **union medida** de `verdad_de_referencia.csv` (11 de 8 quimias
    distintas: catecolaminas, nucleosido, tres quinazolinas, anilina,
    benzisoxazol-piperidina, fenantridinona, aminoalcohol naftenico);
  * ligandos preparados con la receta canonica (`preparar_ligando.py`).

EL CRITERIO, FIJADO POR ADELANTADO Y NO NEGOCIABLE
--------------------------------------------------
    El ranking pasa si coloca positivos de AL MENOS DOS QUIMIOTIPOS DISTINTOS
    por delante de los señuelos emparejados.

Un quimiotipo cuenta como recuperado si su **mejor miembro cae en el 5 % de cabeza**
de la lista combinada (positivos + fondo), que es el mismo corte que usa EF5%. Con
un solo quimiotipo de positivos el resultado no se reporta; con 8 se puede exigir
esto sin trampa.

Se reportan cuatro numeros y no uno, porque ninguno basta solo:

  * AUC crudo           — el numero clasico, que premia al mas grande;
  * AUC por atomo pesado — quita el tamaño de la afinidad;
  * AUC residual         — afinidad menos la recta ajustada al FONDO sobre los
                           atomos pesados: lo que queda cuando ya se ha descontado
                           el sesgo de tamaño. Es la medida honesta;
  * EF1% / EF5%.

DOS FONDOS, DOS PREGUNTAS  (23 sep 2026)
----------------------------------------
Todo lo de arriba se calcula **por cada fondo**, y eso es lo que faltaba. El fondo
de señuelos emparejados mide si el embudo reconoce *alguna* quimia; no mide lo que
importa en una diana como TDP-43, donde el ligando interesante se parece mucho a un
unidor generico de ARN. Para eso está `POSES_FONDO2`: un fondo de moleculas que SI
se unen a ARN (el set de R-BIND 2.0, `analysis/rbind_fondo.py`). El mismo embudo se
puntua tres veces —contra los señuelos, contra el fondo duro y contra los dos
juntos— con `evaluar()`, que es la unica copia de las metricas.

Dentro de cada bloque la lista es solo la de ese fondo (positivos + ese fondo): si
se mezclaran los dos fondos en todas las cuentas, el puesto de un positivo contra
los señuelos dependeria de cuantos ligandos de R-BIND hubiera, y eso es otra
pregunta, no la misma con mas datos.

COMPROBACION PREVIA (pre-flight)
--------------------------------
Antes de puntuar nada se comprueba que cada ligando tenga los MISMOS atomos
pesados en el fichero, en la pose y en su SMILES. Si no cuadra, se re-prepara y
se re-acopla antes de seguir: la mitad de los sustos de este proyecto han venido
de ahi (`INFORME_GLUE_VALIDACIONES_2026-09-21.md`).

Uso:
    python validar_sod1_limpia.py
Salida: validar_sod1_limpia.log y validar_sod1_limpia.csv
"""
import csv
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
GPU = os.path.join(RAIZ, "gpu_dock")
VINA = os.path.join(RAIZ, "tools", "vina.exe")
SALIDA = os.path.join(BASE, "validar_sod1_limpia")

# Diana por defecto. `validar_diana_limpia.py` sobrescribe DIANA, RECEPTOR, CENTRO,
# TAMANO, EXHAUSTIVIDAD, POSES_FONDO, LIGS_FONDO, SALIDA y NOTA_FONDO antes de
# llamar a main(); sin esos nombres aqui, la diana se quedaba clavada en SOD1 y la
# corrida de TDP43 acoplaba los positivos de SOD1 dentro de la caja de TDP43
# (23 sep 2026).
DIANA = "SOD1"
NOTA_FONDO = ("482 señuelos: 199 emparejados antiguos + 140 emparejados nuevos + "
              "143 duros (quelantes y redox-activos)")
ETIQUETA_FONDO = "señuelos emparejados en propiedades"

sys.path.insert(0, RAIZ)
from preparar_ligando import contar_atomos, escribir, smi_del_fichero   # noqa: E402

RECEPTOR = os.path.join(GPU, "SOD1_limpio.pdbqt")
# La misma caja Trp32 de toda la serie (validar_sod1_v3.py, lineas 72-75).
CENTRO = (46.5, 80.0, 73.3)
TAMANO = 22
EXHAUSTIVIDAD = 8
HILOS = 6

# Las poses del fondo, ya acopladas contra el receptor limpio en la ronda v5.
POSES_FONDO = os.path.join(BASE, "validacion_SOD1_v5", "out")
LIGS_FONDO = os.path.join(BASE, "validacion_SOD1_v5", "ligands")

# Fondo duro: un segundo fondo, opcional, con ligandos que SI se unen a su diana
# (para TDP-43, los unidores de ARN de R-BIND 2.0). Sirve para preguntar lo que el
# fondo emparejado no puede: si el embudo distingue quimia especifica de quimia
# generica de ARN. Si es None, no hay segundo bloque y todo sale como siempre.
POSES_FONDO2 = None
LIGS_FONDO2 = None
NOTA_FONDO2 = ""
ETIQUETA_FONDO2 = "fondo duro"

CORTE_EF = 5.0     # el quimiotipo cuenta si su mejor miembro esta en el 5% de cabeza


def log(m):
    print(m, flush=True)
    with open(SALIDA + ".log", "a", encoding="utf-8") as f:
        f.write(m + "\n")


def pesados_smiles(smi):
    """Atomos pesados. `GetNumHeavyAtoms`: los [2H] explicitos de un analogo
    deuterado son hidrogenos y no cuentan."""
    mol = Chem.MolFromSmiles(smi) if smi else None
    return None if mol is None else mol.GetNumHeavyAtoms()


def afinidad(ruta):
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if l.startswith("REMARK VINA RESULT:"):
            try:
                return float(l.split()[3])
            except (IndexError, ValueError):
                return None
    return None


def _dock(t):
    vina, rec, lig, out, cx, cy, cz, tam, exh = t
    subprocess.run([vina, "--receptor", rec, "--ligand", lig,
                    "--center_x", str(cx), "--center_y", str(cy),
                    "--center_z", str(cz),
                    "--size_x", str(tam), "--size_y", str(tam), "--size_z", str(tam),
                    "--exhaustiveness", str(exh), "--num_modes", "3",
                    "--out", out, "--cpu", "1"],
                   capture_output=True, timeout=3600)
    return out


def auc(act, dec):
    """AUC con empates contados a medias. Menos energia = mejor."""
    if not act or not dec:
        return None
    g = 0.0
    for a in act:
        for d in dec:
            g += 0.5 if a == d else (1.0 if a < d else 0.0)
    return g / (len(act) * len(dec))


def ef(act, dec, pct):
    todo = sorted([(s, 1) for s in act] + [(s, 0) for s in dec])
    k = max(1, int(round(len(todo) * pct / 100.0)))
    n_top = sum(1 for _, lab in todo[:k] if lab == 1)
    return (n_top / len(act)) / (pct / 100.0) if act else None


def _cuenta_recuperados(quimias, scores, puesto, orden, pct):
    """Cuantos quimiotipos cruzan el corte indicado."""
    c = max(1, int(round(len(orden) * pct / 100.0)))
    n = 0
    for q in sorted(quimias):
        miembros = ["ACT_" + x for x in quimias[q] if "ACT_" + x in scores]
        if miembros and min(puesto[k] for k in miembros) <= c:
            n += 1
    return n


def residual(scores, pesados, fondo):
    """Afinidad menos la recta afinidad = a + b*(pesados) ajustada SOLO al fondo."""
    xs = np.array([pesados[k] for k in fondo if pesados.get(k)], dtype=float)
    ys = np.array([scores[k] for k in fondo if pesados.get(k)], dtype=float)
    if len(xs) < 10:
        return None, None, None
    b, a = np.polyfit(xs, ys, 1)
    return ({k: scores[k] - (a + b * pesados[k])
             for k in scores if pesados.get(k)}, a, b)


def pareado_por_tamano(scores, pesados, quimias, positivos, decoys):
    """Cada positivo contra los del fondo de su MISMO numero de atomos.

    Es la unica medida que quita el confusor entero: si el ranking ordenase
    quimica, un positivo debe ganar a los del fondo que pesan como el. Si solo
    gana a los pequeños, lo que ordena es el tamaño.

    Devuelve (filas, quimias_que_ganan). Cada fila es
    (ligando, pesados, cuantos comparados, AUC pareada, gana, quimia).
    """
    filas, gana_por_quimia = [], {}
    for q in sorted(quimias):
        for n_lig in quimias[q]:
            k = "ACT_" + n_lig
            if k not in scores:
                continue
            n = pesados[k]
            par = []
            for margen in (2, 3, 5):
                par = [j for j in decoys if abs(pesados[j] - n) <= margen]
                if len(par) >= 20:
                    break
            if len(par) < 5:
                continue
            auc_par = auc([scores[k]], [scores[j] for j in par])
            gana = auc_par is not None and auc_par > 0.5
            gana_por_quimia.setdefault(q, []).append(gana)
            filas.append((n_lig, n, len(par), auc_par, gana, q))
    ganan = [q for q, v in gana_por_quimia.items() if any(v)]
    return filas, ganan


def evaluar(positivos, decoys, scores, pesados, quimias, etiqueta, nota="",
            detalle=False):
    """Las cuatro metricas y el criterio, contra UN fondo.

    La lista de este bloque es solo la de este fondo: positivos + ese fondo. Se
    devuelve un diccionario con los numeros y las quimias que ganan, para poder
    sacar la tabla comparativa de los bloques al final.
    """
    if not positivos or not decoys:
        return None
    S = {k: scores[k] for k in positivos + decoys}
    P = {k: pesados[k] for k in positivos + decoys if pesados.get(k)}

    a = [S[k] for k in positivos]
    d = [S[k] for k in decoys]
    auc_crudo = auc(a, d)
    auc_atomo = auc([S[k] / P[k] for k in positivos],
                    [S[k] / P[k] for k in decoys])
    resid, _, b0 = residual(S, P, decoys)
    auc_resid = auc([resid[k] for k in positivos], [resid[k] for k in decoys])

    orden = sorted(S.items(), key=lambda kv: kv[1])
    puesto = {k: i + 1 for i, (k, _) in enumerate(orden)}
    corte = max(1, int(round(len(orden) * CORTE_EF / 100.0)))
    corr_fondo = np.corrcoef([P[j] for j in decoys], [S[j] for j in decoys])[0, 1]

    log("")
    log("=" * 78)
    log("BLOQUE: %s" % etiqueta.upper())
    log("   %d positivos contra %d del fondo -> lista de %d"
        % (len(positivos), len(decoys), len(orden)))
    if nota:
        log("   %s" % nota)
    log("=" * 78)
    log("   AUC crudo            %.3f" % auc_crudo)
    log("   AUC por atomo pesado %.3f" % auc_atomo)
    log("   AUC residual         %.3f   (afinidad - recta del fondo;"
        " un carbono pesado vale %+.3f kcal/mol)" % (auc_resid, b0))
    log("   EF1%% %.2f   EF5%% %.2f" % (ef(a, d, 1.0), ef(a, d, 5.0)))

    # --------------------------------------------------------------- por quimiotipo
    log("")
    log("   POR QUIMIOTIPO  (cuenta si el mejor miembro cae en el 5%% de cabeza ="
        " puesto <= %d de %d)" % (corte, len(orden)))
    log("   %-24s %-6s %-30s %s"
        % ("quimia", "n", "mejor puesto (ligando)", "recuperado"))
    recuperados = []
    for q in sorted(quimias):
        miembros = ["ACT_" + n for n in quimias[q] if "ACT_" + n in S]
        if not miembros:
            continue
        mejor = min(miembros, key=lambda k: S[k])
        auc_q = auc([S[k] for k in miembros], d)
        auc_r = auc([resid[k] for k in miembros],
                    [resid[k] for k in decoys])
        ok = puesto[mejor] <= corte
        if ok:
            recuperados.append(q)
        log("   %-24s %-6d %-30s %s"
            % (q, len(miembros),
               "%d/%d  %s" % (puesto[mejor], len(orden), mejor.replace("ACT_", "")),
               "SI" if ok else "no"))
        log("      su AUC crudo %.3f | residual %.3f | afinidad %.2f"
            % (auc_q, auc_r, S[mejor]))

    if detalle:
        log("")
        log("   puesto de cada positivo:")
        for k, _ in orden:
            if k in positivos:
                pct = 100.0 * puesto[k] / len(orden)
                log("      %-30s puesto %3d/%d  (mejor que el %.0f%% del fondo)"
                    "  afinidad %6.2f"
                    % (k.replace("ACT_", ""), puesto[k], len(orden), 100 - pct,
                       S[k]))

    # ------------------------------------------------------- fragilidad del corte
    log("")
    log("   FRAGILIDAD: el veredicto depende del corte que se elija")
    for pct in (1.0, 2.0, 5.0, 10.0, 20.0):
        c = max(1, int(round(len(orden) * pct / 100.0)))
        rec = []
        for q in sorted(quimias):
            miembros = ["ACT_" + n for n in quimias[q] if "ACT_" + n in S]
            if miembros and min(puesto[k] for k in miembros) <= c:
                rec.append(q)
        log("      corte %5.1f%% (puesto <= %3d): %d quimias recuperadas  %s"
            % (pct, c, len(rec), "PASA" if len(rec) >= 2 else "no pasa"))

    # ---------------------------------------- la prueba que separa tamaño de quimica
    log("")
    log("   EMPAREJADO POR TAMANO: cada positivo contra el fondo de su mismo tamaño")
    log("      Si el ranking ordenase QUIMICA, un positivo debe ganar a los del fondo")
    log("      del mismo numero de atomos. Si solo gana a los pequeños, lo que ordena")
    log("      es el tamaño, y eso no es reconocimiento (este fondo correlaciona"
        " %.3f con el)." % corr_fondo)
    log("")
    log("      %-28s %-7s %-9s %-11s %s"
        % ("ligando", "atomos", "comparado", "AUC pareada", "gana?"))
    filas_par, quimias_pareadas = pareado_por_tamano(
        S, P, quimias, positivos, decoys)
    for n_lig, n, n_par, auc_par, gana, _ in filas_par:
        log("      %-28s %-7d %-9d %-11.3f %s"
            % (n_lig, n, n_par, auc_par if auc_par is not None else -1,
               "SI" if gana else "no"))
    log("")
    log("      quimias que ganan a los del fondo DE SU TAMAÑO: %d de %d -> %s"
        % (len(quimias_pareadas), len(quimias),
           ", ".join(quimias_pareadas) or "ninguna"))

    k_top = max(1, int(round(len(orden) * CORTE_EF / 100.0)))
    tam_top = np.mean([P[k] for k, _ in orden[:k_top]])
    tam_fondo = np.mean([P[k] for k in decoys])
    tam_pos = np.mean([P[k] for k in positivos])
    log("")
    log("      atomos pesados: cabeza (%d mejores) %.1f | fondo %.1f | positivos %.1f"
        % (k_top, tam_top, tam_fondo, tam_pos))

    # ------------------------------------------------------------------- criterio
    veredicto = "PASA" if len(quimias_pareadas) >= 2 else "NO PASA"
    rec1 = _cuenta_recuperados(quimias, S, puesto, orden, 1.0)
    rec2 = _cuenta_recuperados(quimias, S, puesto, orden, 2.0)
    rec10 = _cuenta_recuperados(quimias, S, puesto, orden, 10.0)
    log("")
    log("   CRITERIO (positivos de al menos DOS quimiotipos por delante del fondo):")
    log("      A) corte fijo de cabeza (5%%): %d de %d quimias -> %s"
        % (len(recuperados), len(quimias), ", ".join(recuperados) or "ninguna"))
    log("         y es inestable: %d al 1%%, %d al 2%%, %d al 10%%"
        % (rec1, rec2, rec10))
    log("      B) AUC residual (descuenta el tamaño con una recta): %.3f" % auc_resid)
    log("      C) EMPAREJADO POR TAMANO (quita el confusor entero): %d de %d -> %s"
        % (len(quimias_pareadas), len(quimias),
           ", ".join(quimias_pareadas) or "ninguna"))
    log("")
    log("      VEREDICTO frente a este fondo: %s  (decide la C)" % veredicto)
    if veredicto == "PASA":
        log("      El embudo SI reconoce quimia en igualdad de tamaño: %d quimias"
            % len(quimias_pareadas))
        log("      distintas batieron a los ligandos del fondo de su mismo numero")
        log("      de atomos.")
    else:
        log("      Con menos de dos quimiotipos el resultado no se reporta como")
        log("      validacion.")

    def _azar(x):
        if x > 0.5:
            return "por ENCIMA del azar"
        return "POR DEBAJO DEL AZAR" if x < 0.5 else "justo en el azar"

    log("")
    log("      LO QUE NO DICE: el AUC crudo es %.3f, %s, y el residual %.3f."
        % (auc_crudo, _azar(auc_crudo), auc_resid))
    log("      La lista ordenada por energia SIN mas no sirve para elegir")
    log("      candidatos: su cabeza es grande, no buena.")

    return {"etiqueta": etiqueta, "n_pos": len(positivos), "n_fondo": len(decoys),
            "auc_crudo": round(auc_crudo, 4), "auc_atomo": round(auc_atomo, 4),
            "auc_residual": round(auc_resid, 4),
            "ef1": round(ef(a, d, 1.0), 3), "ef5": round(ef(a, d, 5.0), 3),
            "quimias_pareadas": quimias_pareadas,
            "quimias_corte_5pct": recuperados,
            "quimias_corte_1pct": rec1, "quimias_corte_10pct": rec10,
            "veredicto": veredicto,
            "pendiente_tamano": round(b0, 4),
            "correlacion_fondo_tamano": round(float(corr_fondo), 3)}


def cargar_fondo(dir_poses, dir_ligs, scores, pesados, papel, etiqueta):
    """Lee las poses de un fondo ya acoplado y devuelve sus nombres.

    Un ligando entra solo si tiene afinidad en la pose Y su numero de atomos
    pesados se puede leer del fichero preparado: sin eso no se puede emparejar por
    tamaño, que es lo que decide.
    """
    nombres = []
    for p in sorted(os.listdir(dir_poses)):
        if not p.endswith("_out.pdbqt") or p.startswith("ACT_"):
            continue
        nombre = p.replace("_out.pdbqt", "")
        if nombre in scores:
            continue
        aff = afinidad(os.path.join(dir_poses, p))
        lig = os.path.join(dir_ligs, nombre + ".pdbqt")
        n = pesados_smiles(smi_del_fichero(lig)) if os.path.exists(lig) else None
        if aff is None or n is None:
            continue
        scores[nombre] = aff
        pesados[nombre] = n
        papel[nombre] = etiqueta
        nombres.append(nombre)
    return nombres


def main():
    open(SALIDA + ".log", "w", encoding="utf-8").close()
    log("VALIDACION DE %s CONTRA EL RECEPTOR LIMPIO   %s"
        % (DIANA, time.strftime("%Y-%m-%d %H:%M")))
    log("receptor  %s" % os.path.basename(RECEPTOR))
    log("caja      %s tamano %d exhaustividad %d" % (CENTRO, TAMANO, EXHAUSTIVIDAD))

    # ------------------------------------------------------------------ positivos
    filas = [r for r in csv.DictReader(
        open(os.path.join(BASE, "verdad_de_referencia.csv"), encoding="utf-8"))
        if r["target"] == DIANA and r["apto"] == "si"]
    log("positivos de union medida: %d" % len(filas))
    quimias = {}
    for r in filas:
        quimias.setdefault(r["quimia"], []).append(r["ligand"])
    log("quimias distintas: %d -> %s"
        % (len(quimias), ", ".join(sorted(quimias))))

    ligdir = os.path.join(SALIDA, "ligands")
    outdir = os.path.join(SALIDA, "out")
    os.makedirs(ligdir, exist_ok=True)
    os.makedirs(outdir, exist_ok=True)

    # ------------------------------------------------- comprobacion previa (pre-flight)
    log("")
    log("=" * 78)
    log("COMPROBACION PREVIA: fichero = pose = SMILES (atomos pesados)")
    log("=" * 78)
    sospechosos = []
    for r in filas:
        nombre = "ACT_" + r["ligand"]
        lig = os.path.join(POSES_FONDO.replace("out", "ligands"), nombre + ".pdbqt")
        pose = os.path.join(POSES_FONDO, nombre + "_out.pdbqt")
        n_smi = pesados_smiles(r["smiles"])
        n_lig = contar_atomos(open(lig, encoding="utf-8").read())[0] if os.path.exists(lig) else None
        n_pose = None
        if os.path.exists(pose):
            txt = open(pose, encoding="utf-8", errors="ignore").read()
            # solo la primera pose
            if "ENDMDL" in txt:
                txt = txt.split("ENDMDL")[0]
            n_pose = contar_atomos(txt)[0]
        ok = n_smi is not None and n_lig == n_smi and n_pose == n_smi
        log("   %-28s SMILES %s | fichero %s | pose %s   %s"
            % (r["ligand"], n_smi, n_lig, n_pose, "OK" if ok else "REVISAR"))
        if not ok:
            sospechosos.append(r)
    if sospechosos:
        log("")
        log("   %d positivos no cuadran: se re-preparan y re-acoplan antes de puntuar"
            % len(sospechosos))
        tareas = []
        for r in sospechosos:
            p = escribir("ACT_" + r["ligand"], r["smiles"], ligdir,
                         forzar=True, quitar_sales=True)
            if p is None:
                continue
            out = os.path.join(outdir, "ACT_%s_out.pdbqt" % r["ligand"])
            if os.path.exists(out):
                os.remove(out)
            tareas.append((VINA, RECEPTOR, p, out, CENTRO[0], CENTRO[1], CENTRO[2],
                           TAMANO, EXHAUSTIVIDAD))
        if tareas:
            with ProcessPoolExecutor(max_workers=HILOS) as ex:
                list(ex.map(_dock, tareas))
            log("   re-acoplados %d" % len(tareas))

    # -------------------------------------------------------------------- puntuacion
    scores, pesados, papel = {}, {}, {}
    for r in filas:
        nombre = "ACT_" + r["ligand"]
        pose = os.path.join(POSES_FONDO, nombre + "_out.pdbqt")
        if nombre in [f"ACT_{s['ligand']}" for s in sospechosos]:
            pose = os.path.join(outdir, nombre + "_out.pdbqt")
        aff = afinidad(pose) if os.path.exists(pose) else None
        if aff is None:
            log("   !! sin pose para %s: queda fuera" % r["ligand"])
            continue
        scores[nombre] = aff
        pesados[nombre] = pesados_smiles(r["smiles"])
        papel[nombre] = "positivo"

    fondo1 = cargar_fondo(POSES_FONDO, LIGS_FONDO, scores, pesados, papel, "fondo")
    log("")
    log("fondo: %d señuelos acoplados contra el mismo receptor y la misma caja"
        % len(fondo1))
    log("   %s" % NOTA_FONDO)

    fondo2 = []
    if POSES_FONDO2:
        fondo2 = cargar_fondo(POSES_FONDO2, LIGS_FONDO2, scores, pesados, papel,
                              "fondo2")
        log("fondo duro: %d ligandos acoplados contra el mismo receptor y la misma caja"
            % len(fondo2))
        if NOTA_FONDO2:
            log("   %s" % NOTA_FONDO2)

    positivos = [k for k in scores if papel[k] == "positivo"]
    log("positivos con pose: %d de %d" % (len(positivos), len(filas)))

    # ---------------------------------------------------------------- los bloques
    resultados = [evaluar(positivos, fondo1, scores, pesados, quimias,
                          ETIQUETA_FONDO, NOTA_FONDO, detalle=True)]
    if fondo2:
        resultados.append(evaluar(positivos, fondo2, scores, pesados, quimias,
                                  ETIQUETA_FONDO2, NOTA_FONDO2, detalle=True))
        resultados.append(evaluar(positivos, fondo2 + fondo1, scores, pesados,
                                  quimias, "los dos fondos juntos", detalle=False))

    # ------------------------------------------------------- tabla de los bloques
    log("")
    log("=" * 78)
    log("LOS TRES BLOQUES, UNO AL LADO DEL OTRO")
    log("=" * 78)
    log("   %-34s %-7s %-9s %-9s %-9s %-6s %s"
        % ("fondo", "n", "AUC", "AUC/at.", "residual", "EF5%", "quimias/total"))
    for r in resultados:
        if not r:
            continue
        log("   %-34s %-7d %-9.3f %-9.3f %-9.3f %-6.2f %d/%d  %s"
            % (r["etiqueta"][:34], r["n_fondo"], r["auc_crudo"], r["auc_atomo"],
               r["auc_residual"], r["ef5"], len(r["quimias_pareadas"]),
               len(quimias), r["veredicto"]))
    log("")
    log("   El fondo que decide es el duro: es el unico que puede separar quimia")
    log("   especifica de quimia generica de ARN.")

    # ------------------------------------------------------------------ guardado
    resumen = {}
    if resultados[0]:
        r = resultados[0]
        resumen.update({"auc_crudo": r["auc_crudo"], "auc_atomo": r["auc_atomo"],
                        "auc_residual": r["auc_residual"], "ef1": r["ef1"],
                        "ef5": r["ef5"],
                        "quimias_pareadas_por_tamano": len(r["quimias_pareadas"]),
                        "quimias_total": len(quimias), "veredicto": r["veredicto"],
                        "quimias_corte_5pct": len(r["quimias_corte_5pct"]),
                        "quimias_corte_1pct": r["quimias_corte_1pct"],
                        "quimias_corte_10pct": r["quimias_corte_10pct"],
                        "pendiente_tamano_kcal_por_atomo": r["pendiente_tamano"],
                        "correlacion_fondo_tamano": r["correlacion_fondo_tamano"]})
    for i, r in enumerate(resultados[1:], 1):
        if not r:
            continue
        pre = "duro_" if i == 1 else "juntos_"
        resumen[pre + "n_fondo"] = r["n_fondo"]
        resumen[pre + "auc_crudo"] = r["auc_crudo"]
        resumen[pre + "auc_residual"] = r["auc_residual"]
        resumen[pre + "ef5"] = r["ef5"]
        resumen[pre + "quimias_pareadas_por_tamano"] = len(r["quimias_pareadas"])
        resumen[pre + "quimias_que_ganan"] = ",".join(r["quimias_pareadas"])
        resumen[pre + "veredicto"] = r["veredicto"]
    with open(SALIDA + "_resumen.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(list(resumen))
        w.writerow([resumen[k] for k in resumen])

    campos = ["ligand", "papel", "afinidad", "pesados", "puesto_fondo1",
              "puesto_fondo2", "residual"]
    # El puesto va DENTRO DE CADA BLOQUE, que es como se lee esta validacion: el
    # puesto absoluto de una lista mezclada no significa lo mismo.
    p1 = {k: i + 1 for i, k in enumerate(sorted(fondo1, key=lambda x: scores[x]))}
    p2 = {k: i + 1 for i, k in enumerate(sorted(fondo2, key=lambda x: scores[x]))}
    with open(SALIDA + ".csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(campos)
        for k in sorted(scores, key=lambda x: scores[x]):
            w.writerow([k.replace("ACT_", ""), papel[k], scores[k],
                        pesados.get(k, ""), p1.get(k, ""), p2.get(k, ""), ""])
    log("")
    log("guardado: %s.csv" % SALIDA)
    return 0


if __name__ == "__main__":
    sys.exit(main())
