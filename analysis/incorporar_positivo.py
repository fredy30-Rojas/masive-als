#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""incorporar_positivo.py — el portero de los positivos que lleguen del laboratorio.

POR QUE ES UN PORTERO Y NO UN AÑADIDOR
--------------------------------------
La verdad de referencia es lo unico que sostiene los veredictos de la regla del 24 de
septiembre, asi que meter una fila a mano es lo mas caro que se puede hacer mal en este
proyecto. Este script **no escribe** en `verdad_de_referencia.csv`: comprueba la entrada,
dice si vale, y deja la fila preparada en `positivos_pendientes.csv` para que la revise
una persona y la mueva ella. Es la misma regla de siempre: lo que no se deshace se
pregunta.

QUE COMPRUEBA, EN ESTE ORDEN
----------------------------
1. **Que sea union medida.** El campo del ensayo tiene que decir como se midio (Kd, SPR,
   ITC, MST, CSP, RMN, CETSA). Si solo hay actividad celular (EC50, IC50, viabilidad,
   agregacion) no entra: es exactamente el error que dejo a SOD1 sin evidencia.
2. **Que sea del sitio de la caja.** Un unido medido en otro bolsillo (el caso de XL20,
   que es del CR y no del RRM) se guarda pero NO entra en el bloque de esa caja.
3. **Que la estructura cuadre.** SMILES valido en RDKit y numero de atomos pesados, que
   es lo que la regla divide para quitar el sesgo de tamaño.
4. **Que sea independiente.** Tanimoto ECFP4 y esqueleto de Murcko contra los positivos
   que ya estan: mismo Murcko = congenere (suma n pero no independencia), >= 0,4 =
   analogo, y por debajo, quimica nueva. Esto no es un veto, es un aviso con numero: la
   regla crece con positivos independientes, y con analogos crece menos de lo que parece.

Uso:
    python analysis/incorporar_positivo.py --target TDP43 --ligando XL23 \\
        --smiles "Nc1ncnc2c1ncn2[C@H]1CCC[C@H](NC(=O)Cn2c(C(F)(F)F)nc3ccccc32)[C@@H]1O" \\
        --sitio "CR del C-terminal (320-340)" --ensayo "SPR Kd 45 uM + CETSA" \\
        --cita "nuestro ensayo, 2026-10" --quimia "adenina-aminociclohexanol"

Salidas:
    analysis/positivos_pendientes.csv   la fila preparada, con el veredicto y el motivo
"""

import argparse
import csv
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
VERDAD = os.path.join(BASE, "verdad_de_referencia.csv")
PENDIENTES = os.path.join(BASE, "positivos_pendientes.csv")

# La caja de cada diana, tal como esta en validar_diana_limpia.py. El sitio de un
# positivo tiene que caer aqui dentro, o no es un positivo de ese bloque.
CAJAS = {
    "TDP43": ("RRM1-RRM2 / interfaz Arg151-Asp247 (PDB 4BS2, caja 24,23 16,89 -15,87)",
              ["rrm1", "rrm2", "interfaz", "rrm"]),
    "SOD1": ("Trp32 (PDB 1HL5, caja 46,5 80,0 73,3)",
             ["trp32", "triptofano 32", "w32"]),
}

# Palabras que indican union medida de verdad, y palabras que NO bastan.
UNION = ["kd", "ki", "ec50 de union", "spr", "itc", "mst", "csp", "rmn", "nmr", "cetesa",
         "calorimetria", "fluorescencia de union", "tr-fret", "biolayer", "blitz",
         "microscale", "resonancia"]
NO_UNION = ["ec50 celular", "ic50", "viabilidad", "mtt", "agregacion", "agregación",
            "supervivencia", "neuronal", "western", "inmunofluorescencia", "similitud"]

COLS = ["fecha", "target", "ligand", "smiles", "quimia", "sitio", "tipo_ensayo", "cita",
        "pesados", "sitio_ok", "union_medida", "independencia", "tanimoto_max",
        "contra_quien", "murcko_igual", "veredicto", "motivo"]


def leer_verdad():
    if not os.path.exists(VERDAD):
        return []
    with open(VERDAD, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def smiles_valido(smi):
    """Atomos pesados y molecule construida; None si el SMILES no vale."""
    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except Exception:
        return None
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    return m


def murcko(m):
    from rdkit.Chem.Scaffolds import MurckoScaffold
    return Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))


def tanimoto(a, b):
    from rdkit.Chem import AllChem, DataStructs
    fa = AllChem.GetMorganFingerprintAsBitVect(a, 2, nBits=2048)
    fb = AllChem.GetMorganFingerprintAsBitVect(b, 2, nBits=2048)
    return float(DataStructs.TanimotoSimilarity(fa, fb))


def comprobar_sitio(target, sitio):
    """¿El sitio declarado cae dentro de la caja de esa diana?"""
    desc, claves = CAJAS[target]
    s = (sitio or "").lower()
    return any(k in s for k in claves), desc


def comprobar_ensayo(ensayo):
    """(es_union, motivo)."""
    e = (ensayo or "").lower()
    if not e.strip():
        return False, "sin tipo de ensayo declarado"
    tiene_no = [p for p in NO_UNION if p in e]
    tiene_si = [p for p in UNION if p in e]
    if tiene_si:
        return True, "union medida por %s" % ", ".join(sorted(set(tiene_si)))
    if tiene_no:
        return False, ("solo actividad no de union (%s): no entra sin una medida de union a "
                       "proteina" % ", ".join(sorted(set(tiene_no))))
    return False, "el ensayo no dice como se midio la union: hay que concretarlo"


def independencia(m, verdad, target):
    """(etiqueta, tanimoto_max, contra_quien, murcko_igual)."""
    mejor, quien = 0.0, ""
    igual = []
    for r in verdad:
        if r["target"] != target or r["apto"] != "si" or not r["smiles"].strip():
            continue
        otro = smiles_valido(r["smiles"])
        if otro is None:
            continue
        t = tanimoto(m, otro)
        if t > mejor:
            mejor, quien = t, r["ligand"]
        try:
            if murcko(otro) == murcko(m):
                igual.append(r["ligand"])
        except Exception:
            pass
    if igual:
        etiqueta = "congenere (%s)" % ", ".join(igual[:3])
    elif mejor >= 0.4:
        etiqueta = "analogo (Tanimoto %.3f)" % mejor
    else:
        etiqueta = "quimica nueva"
    return etiqueta, mejor, quien, bool(igual)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, choices=sorted(CAJAS))
    ap.add_argument("--ligando", required=True)
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--quimia", default="")
    ap.add_argument("--sitio", required=True)
    ap.add_argument("--ensayo", required=True)
    ap.add_argument("--cita", required=True)
    ap.add_argument("--fecha", default="")
    args = ap.parse_args()

    verdad = leer_verdad()
    m = smiles_valido(args.smiles)
    pesados = m.GetNumHeavyAtoms() if m is not None else None
    sitio_ok, caja = comprobar_sitio(args.target, args.sitio)
    union_ok, motivo_union = comprobar_ensayo(args.ensayo)
    if m is None:
        indep, t, quien, igual = "SMILES invalido", 0.0, "", False
    else:
        indep, t, quien, igual = independencia(m, verdad, args.target)

    motivos = []
    if m is None:
        motivos.append("SMILES que RDKit no construye")
    if not sitio_ok:
        motivos.append("el sitio declarado no es el de la caja de %s (%s)" % (args.target, caja))
    if not union_ok:
        motivos.append(motivo_union)
    if not motivos and igual:
        motivos.append("entra, pero es congenere de %s: suma positivos, no independencia"
                       % ", ".join(igual[:3]))
    veredicto = "no_entra" if (m is None or not sitio_ok or not union_ok) else \
                ("entra_congenere" if igual else "entra")

    fila = dict(fecha=args.fecha, target=args.target, ligand=args.ligando,
                smiles=args.smiles, quimia=args.quimia, sitio=args.sitio,
                tipo_ensayo=args.ensayo, cita=args.cita,
                pesados="" if pesados is None else pesados,
                sitio_ok="si" if sitio_ok else "no",
                union_medida="si" if union_ok else "no",
                independencia=indep, tanimoto_max="%.3f" % t, contra_quien=quien,
                murcko_igual="si" if igual else "no",
                veredicto=veredicto, motivo="; ".join(motivos) if motivos else "cumple el criterio")

    nuevo = not os.path.exists(PENDIENTES)
    with open(PENDIENTES, "a", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=COLS)
        if nuevo:
            w.writeheader()
        w.writerow(fila)

    print("=" * 78)
    print("PORTERO DE POSITIVOS: %s -> %s" % (args.ligando, args.target))
    print("=" * 78)
    print("  caja de la diana      %s" % caja)
    print("  sitio declarado       %s  %s" % (args.sitio, "OK" if sitio_ok else "NO ES DE ESTA CAJA"))
    print("  ensayo                %s  %s" % (args.ensayo, "OK" if union_ok else "NO VALE"))
    print("  atomos pesados        %s" % pesados)
    print("  independencia         %s  (Tanimoto max %.3f contra %s)"
          % (indep, t, quien or "nadie"))
    print("")
    print("  VEREDICTO: %s" % veredicto.upper())
    print("  %s" % fila["motivo"])
    print("")
    print("  Fila preparada en %s (revisar y mover a verdad_de_referencia.csv a mano)."
          % os.path.basename(PENDIENTES))
    if veredicto != "no_entra":
        print("")
        print("  Si se mueve a la verdad, la secuencia es esta (y en este orden):")
        print("")
        if args.target == "TDP43":
            print("    python analysis/validar_diana_limpia.py --target TDP43   # exh 16, los dos fondos")
            print("    python analysis/barrido_exhaustividad_tdp43.py            # exh 8 y 32")
            print("    python analysis/control_gpu_tdp43.py                      # la corrida de GPU")
        else:
            print("    python analysis/validar_diana_limpia.py --target SOD1     # exh 8, los tres fondos")
            print("    python analysis/separar_fondos_sod1.py                    # el fondo, partido por prefijo")
        print("    python analysis/regla_decision_bootstrap.py               # el veredicto nuevo")
        print("")
        print("  Antes de correr la regla, guarde el veredicto de hoy:")
        print("    mkdir -p analysis/regla_decision/_antes_AAAA-MM-DD")
        print("    cp analysis/regla_decision/regla_decision.* analysis/regla_decision/_antes_AAAA-MM-DD/")
    return 0 if veredicto != "no_entra" else 1


if __name__ == "__main__":
    sys.exit(main())
