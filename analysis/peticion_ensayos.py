#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""peticion_ensayos.py — que hay que medir en el laboratorio, y para que.

POR QUE EXISTE ESTE SCRIPT
--------------------------
La regla de decision del 24 de septiembre dice donde esta el cuello de botella: no
en el computo, en los positivos de union medida. El bloque duro de TDP-43 pasa con
7; los señuelos emparejados de TDP-43 necesitan del orden de 10-13 en total (o sea
3-6 mas) y el fondo duro de SOD1 unos 22 (11 mas). Subir exhaustividad mueve el
criterio sin acercarlo a la verdad, y el barrido del PDB del 25 de septiembre
demostro que por estructura no van a salir: 22 ligandos co-cristalizados de SOD1,
cero de TDP-43 y cero de FUS.

Este script no inventa nada: lee lo que ya esta en el repositorio (la verdad de
referencia, los controles de SOD1, las tablas de candidatos, el suplementario del
articulo del CR y la regla de decision) y escribe la peticion con los SMILES y los
identificadores que ya estan comprobados, para que la tabla se pueda mandar tal cual.

SALIDAS
-------
    analysis/peticion_ensayos.csv    una fila por compuesto y ensayo propuesto
    analysis/peticion_ensayos.txt    el resumen: que falta, cuanto, y por que cada uno

Uso:
    python analysis/peticion_ensayos.py
"""

import csv
import os
import sys

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
A = os.path.join(BASE, "analysis")

VERDAD = os.path.join(A, "verdad_de_referencia.csv")
CONTROLES_SOD1 = os.path.join(A, "controles_sod1_v4.csv")
SUPLEMENTARIO_CR = os.path.join(A, "suplementario_tdp43_xl21_27.csv")
CANDIDATOS = os.path.join(A, "candidatos_filtrados.csv")
CANDIDATOS_CNS = os.path.join(A, "candidatos_total_cns.csv")
CANDIDATOS_V4 = os.path.join(A, "_candidatos_v4.csv")
REGLA = os.path.join(A, "regla_decision", "regla_decision.csv")

SALIDA_CSV = os.path.join(A, "peticion_ensayos.csv")
SALIDA_TXT = os.path.join(A, "peticion_ensayos.txt")


def leer(ruta):
    """Lee un CSV y devuelve la lista de diccionarios; vacia si no esta."""
    if not os.path.exists(ruta):
        return []
    with open(ruta, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def por_nombre(rows, campo, valor):
    for r in rows:
        if r.get(campo, "").strip().lower() == valor.strip().lower():
            return r
    return {}


def candidato(rows, ligando, target):
    """El candidato de una diana con su SMILES ya calculado en la tabla."""
    for r in rows:
        if r.get("ligand", "").strip() == ligando and r.get("target", "").strip() == target:
            return r
    return {}


# ---------------------------------------------------------------------------
# La peticion, fila a fila. Todo lo que va entre comillas es lo que se pide y
# por que; los SMILES y los identificadores los pone el repositorio.
# ---------------------------------------------------------------------------

# Identificadores de catalogo comprobados en las fichas de proveedor (Sigma, Cayman,
# Selleckchem, MedChemExpress) el 25 de septiembre de 2026. Solo se pone lo comprobado.
CATALOGO = {
    "LCS-1": "CAS 41931-13-9",
}


def filas():
    verdad = leer(VERDAD)
    sod1_controles = leer(CONTROLES_SOD1)
    cr = leer(SUPLEMENTARIO_CR)
    cand = leer(CANDIDATOS)
    cns = leer(CANDIDATOS_CNS)
    v4 = leer(CANDIDATOS_V4)

    f = []

    # ---- Prioridad 1: cada uno mueve un numero que el proyecto ya publico ----

    xl23 = por_nombre(cr, "ligand", "XL23")
    xl20 = por_nombre(verdad, "ligand", "XL20")
    f.append(dict(
        prioridad=1,
        diana="TDP-43",
        sitio="CR del C-terminal (320-340)",
        compuesto="XL23",
        identificador=xl23.get("asinex_id", ""),
        smiles=xl23.get("smiles", ""),
        de_donde_sale="Gao et al. 2026, Nature Aging 6:1667 (Fig. suplementaria 1a)",
        que_medir="SPR de union aparente con el mismo protocolo que dio micromolar para XL20, "
                  "mas CETSA en celulas (engagement de diana). Si el laboratorio tiene el "
                  "constructo del CR, repetir la dependencia de sitio: sin CR no une, y con "
                  "Trp334 mutado deberia bajar.",
        por_que="Es el hermano de quimiotipo de XL20 (19 atomos pesados comunes, adenina + "
                "aminociclohexanol) y el unico compuesto de esa tabla con actividad funcional "
                "(inhibe la agregacion del LCD a 100 uM) sin union medida. Sin este dato la "
                "pareja es una cabeza comun con dos colas.",
        que_cambia="Pasa de 1 a 2 uniones medidas en el CR y convierte la pareja en una "
                   "relacion estructura-actividad de dos puntos: ahi si se puede montar un "
                   "control de verdad (¿reproduce el acoplamiento la dependencia del Trp334?).",
        bloque_que_alimenta="CR (linea cerrada hoy por falta de un segundo dato)",
    ))

    lcs1 = por_nombre(verdad, "ligand", "LCS-1")
    prg = por_nombre(verdad, "ligand", "PRG-A01")
    pyra = por_nombre(verdad, "ligand", "CHEMBL2165613")
    for fila, nota_quim in (
        (lcs1, "piridazinona"),
        (prg, "cumarina"),
        (pyra, "pirazolona (representante de la serie de 18)"),
    ):
        if not fila:
            continue
        f.append(dict(
            prioridad=1,
            diana="SOD1",
            sitio="por determinar (Trp32, Cys111 u otro)",
            compuesto=fila["ligand"],
            identificador=CATALOGO.get(fila["ligand"], ""),
            smiles=fila.get("smiles", ""),
            de_donde_sale=fila.get("cita", "") or "controles de SOD1 del proyecto",
            que_medir="Kd por SPR o ITC contra SOD1 humana, y determinacion del sitio "
                      "(soaking en cristal o RMN). Las dos cosas: la afinidad sola no dice a "
                      "que bolsillo, y el sitio es lo que decide si entra en la validacion.",
            por_que=f"Entra hoy en la verdad de referencia como {nota_quim} pero SIN union a "
                    "proteina medida: su evidencia es actividad celular (o, en la pirazolona, "
                    "ser el representante de una serie congenica colapsada). Es materia prima "
                    "que ya esta en el congelador de la literatura y no exige sintetizar nada.",
            que_cambia="Cada uno que se confirme suma un positivo quimicamente independiente al "
                       "fondo duro de SOD1, que es el bloque que decide (hoy 11 positivos, "
                       "0,643, y hacen falta del orden de 22 para cerrar el intervalo).",
            bloque_que_alimenta="SOD1, fondo duro (143 quelantes y redox)",
        ))

    cloza = por_nombre(cns, "ligand", "CLOZAPINE")
    f.append(dict(
        prioridad=1,
        diana="SOD1",
        sitio="Trp32 (caja del proyecto)",
        compuesto="clozapina",
        identificador="CHEMBL42",
        smiles=cloza.get("smiles", ""),
        de_donde_sale="candidatos_total_cns.csv (%s kcal/mol, aprobado, CNS-MPO %s)"
                      % (cloza.get("energy", "?"), cloza.get("cns_mpo", "?")),
        que_medir="SPR contra SOD1 humana en la caja de Trp32 (una concentracion de cribado y, "
                  "si hay señal, Kd completa).",
        por_que="Es un farmaco aprobado, quimicamente independiente de todo lo que hay en la "
                "verdad de referencia (dibenzodiazepina frente a catecolaminas, quinazolinas y "
                "nucleosidos) y el propio cribado lo puso a -7,5 kcal/mol con CNS-MPO 4,77.",
        que_cambia="Si une, es a la vez un positivo nuevo para el bloque y un candidato de "
                   "reposicionamiento con dato medido; si no une, el cribado pierde ese "
                   "candidato y queda dicho.",
        bloque_que_alimenta="SOD1, fondo duro / candidatos",
    ))

    # ---- Prioridad 2: panel propio, para saber si el embudo acierta ----

    panel_sod1 = ["CHEMBL3311449", "CHEMBL4553700", "CHEMBL584356", "CHEMBL8550"]
    for lig in panel_sod1:
        r = candidato(cand, lig, "SOD1") or por_nombre(v4, "ligand", lig)
        if not r:
            continue
        f.append(dict(
            prioridad=2,
            diana="SOD1",
            sitio="Trp32 (caja del proyecto)",
            compuesto=lig,
            identificador=lig,
            smiles=r.get("smiles", ""),
            de_donde_sale="candidatos_filtrados.csv / _candidatos_v4.csv"
                          + (" (mejor MM-GBSA)" if lig == "CHEMBL3311449" else ""),
            que_medir="SPR de cribado (una concentracion) y Kd si hay señal.",
            por_que="Es de los mejor puntuados por el embudo (%s kcal/mol por Vina%s) y no se "
                    "parece a ningun control."
                    % (r.get("energy", r.get("affinity", "?")),
                       ", dG %s kcal/mol por MM-GBSA" % r["mmgbsa_dG"]
                       if r.get("mmgbsa_dG") else ""),
            que_cambia="Cada uno que una es un positivo nuevo; cada uno que no una es un error "
                       "del ranking que se apunta. En los dos casos el numero mejora, porque "
                       "hoy solo hay 11 positivos y son casi todos del mismo par de campanas "
                       "cristalograficas.",
            bloque_que_alimenta="SOD1, fondo duro / candidatos",
        ))

    panel_tdp43 = ["CHEMBL1362588", "CHEMBL27093", "CHEMBL19347", "CHEMBL23927"]
    for lig in panel_tdp43:
        r = candidato(cand, lig, "TDP43")
        if not r:
            continue
        f.append(dict(
            prioridad=2,
            diana="TDP-43",
            sitio="RRM1 (caja del proyecto)",
            compuesto=lig,
            identificador=lig,
            smiles=r.get("smiles", ""),
            de_donde_sale="candidatos_filtrados.csv (%s kcal/mol)" % r.get("energy", "?"),
            que_medir="SPR o RMN de 15N-HSQC con RRM1-RRM2, que ademas localiza el sitio por "
                      "el desplazamiento quimico (CSP).",
            por_que="Los 7 positivos que sostienen el bloque de TDP-43 no son del mismo sitio: "
                    "rTRD01 es RRM1, nTRD22 es del NTD, tres fragmentos son RRM2 y PE859 y "
                    "berberrubina son de la interfaz RRM1-RRM2. Un positivo de RRM1 mas vale "
                    "aqui mas que tres de otros sitios.",
            que_cambia="Es el bloque que ya PASA (AUC/atomo 0,73 contra el fondo duro): con 3-6 "
                       "positivos mas y del mismo sitio se cierra tambien el intervalo contra "
                       "los señuelos emparejados, que es lo unico que le falta.",
            bloque_que_alimenta="TDP-43, señuelos emparejados (10-13 en total)",
        ))

    # ---- Prioridad 3: positivos del CTD, con su aviso ----

    for lig in ("bis-ANS", "CongoRed"):
        r = por_nombre(verdad, "ligand", lig)
        if not r:
            continue
        f.append(dict(
            prioridad=3,
            diana="TDP-43",
            sitio="CTD / LLPS (320-340)",
            compuesto=lig,
            identificador="",
            smiles=r.get("smiles", ""),
            de_donde_sale=r.get("cita", ""),
            que_medir="Kd contra el CTD o el CR (fluorescencia o ITC).",
            por_que="Tienen evidencia funcional de separacion de fases y agregacion y ninguna "
                    "union medida, pero son sondas promiscuas: bis-ANS se une a casi todo lo "
                    "que tenga una cara hidrofoba expuesta.",
            que_cambia="Darian positivos del CTD, no del RRM. Solo valen si algun dia se ataca "
                       "el CR como diana, y aun asi habria que distinguir union especifica de "
                       "pegado inespecifico.",
            bloque_que_alimenta="CR / CTD (informativo, no cierra ningun bloque)",
        ))

    return f


def resumen_regla():
    """Lo que falta, leido de la regla de decision y no de la memoria."""
    rows = leer(REGLA)
    out = []
    for r in rows:
        if r.get("bloque") not in ("duro", "blando"):
            continue
        n = r.get("n_pos", "")
        falta = r.get("positivos_que_harian_falta", "").strip()
        if not falta:
            # La columna queda vacia cuando el bloque ya pasa: no falta ninguno.
            necesita = "ya basta con los que hay"
        else:
            try:
                necesita = "necesita ~%s (%d mas)" % (falta, int(falta) - int(n))
            except ValueError:
                necesita = "necesita ~%s" % falta
        out.append("%-7s %-30s bloque %-6s hoy %-3s  %-30s -> %s"
                   % (r["target"], r["corrida"], r["bloque"], n, necesita, r["veredicto"]))
    return out


def pintar(f):
    lineas = []
    def log(t=""):
        lineas.append(t)
        print(t)

    log("PETICION DE ENSAYOS DE UNION MEDIDA — %d compuestos" % len(f))
    log("=" * 78)
    log("")
    log("QUE FALTA, SEGUN LA REGLA DE DECISION DEL 24 DE SEPTIEMBRE 2026")
    log("-" * 78)
    for l in resumen_regla():
        log("  " + l)
    log("")
    log("  El bloque duro de TDP-43 ya pasa y no necesita nada. Lo que falta es:")
    log("    TDP-43 contra señuelos emparejados: 3-6 positivos mas, y del MISMO sitio")
    log("    SOD1 contra su fondo duro:           del orden de 11 mas")
    log("  Y no van a salir del PDB: el barrido del 25 de septiembre da 22 ligandos")
    log("  co-cristalizados de SOD1, cero de TDP-43 y cero de FUS.")
    log("")
    log("LOS COMPUESTOS, POR PRIORIDAD")
    log("-" * 78)
    for p in (1, 2, 3):
        sel = [r for r in f if r["prioridad"] == p]
        if not sel:
            continue
        log("")
        log("  Prioridad %d (%d): %s" % (
            p, len(sel),
            {1: "cada uno mueve un numero ya publicado",
             2: "panel propio: dicen si el embudo acierta",
             3: "informativos, no cierran ningun bloque"}[p]))
        for r in sel:
            log("    %-9s %-26s %-26s %s" % (r["diana"], r["compuesto"][:25],
                                             r["identificador"][:25], r["sitio"][:30]))
    log("")
    log("TABLA COMPLETA: %s" % SALIDA_CSV)
    with open(SALIDA_TXT, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lineas) + "\n")


def main():
    f = filas()
    if not f:
        print("No hay filas: faltan los CSV de entrada.")
        return 1
    cols = ["prioridad", "diana", "sitio", "compuesto", "identificador", "smiles",
            "de_donde_sale", "que_medir", "por_que", "que_cambia", "bloque_que_alimenta"]
    with open(SALIDA_CSV, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        for r in f:
            w.writerow(r)
    pintar(f)
    return 0


if __name__ == "__main__":
    sys.exit(main())
