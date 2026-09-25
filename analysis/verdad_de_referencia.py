#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verdad_de_referencia.py — el conjunto de positivos con el que SI se puede medir.

POR QUE EXISTE
--------------
La auditoria del 20 de septiembre encontro la causa raiz de que ninguna metrica del
proyecto sea evaluable: **la verdad de referencia esta sesgada**.

  * SOD1 tenia 20 «activos» y **16 eran el mismo nucleo pirazolona**; 18 tenian otro
    activo a similitud >= 0,5. Sin la serie, el AUC cae de 0,815 a 0,601 con dos
    compuestos. Eso no mide generalidad quimica: mide recuperar la quimia con la que
    se sembro.
  * TDP-43 tenia 8 controles y **6 entraron por parecerse a la berberrubina**. La
    similitud da AUC 1,000 dentro de su familia y 0,538 fuera.

Mientras esto no se arregle, ninguna metrica —ni el acoplamiento, ni el MM-GBSA, ni
la huella de interaccion— se puede evaluar. Es anterior a la metrica.

QUE HACE
--------
Consolida en un solo fichero, `verdad_de_referencia.csv`, los positivos de las tres
dianas con **cita y tipo de ensayo**, **colapsando las series congenericas a un
representante**. La regla de entrada, fijada por adelantado y sin excepciones:

    sin cita y sin tipo de ensayo, el compuesto no entra.

Y se distingue lo que vale para medir de lo que no:

    union_directa   co-cristalizacion, MST, RMN (CSP), Kd medida  -> APTO
    funcional       reduce agregacion / rescata fenotipo, sin union medida -> NO
    celular         actividad en celula, sin union medida -> NO
    similitud       entro por parecerse a otro control -> NO (fuera del fichero)

De donde sale cada dato, sin inventar nada:

  * SOD1 co-cristalizados: `controles_sod1_v4.csv` (cimiento, cita, PDB) y
    `redocking_trp32/controles_corregidos.csv` (veredicto de redocking).
  * SOD1 independientes: `activos_sod1_v2.csv`.
  * TDP-43: la tabla y las fuentes de `PLAN_TDP43_2026-09-20.md`. Ampliado el 25 de
    septiembre con XL20 (Gao 2026, Nature Aging), unido de verdad al dominio C-terminal
    pero NO al bolsillo de RRM: entra marcado como no apto para no contaminar la
    validacion de RRM, y su cita queda registrada.
  * FUS: `activos_fus.csv` (dos compuestos, sin literatura propia).

Uso:
    python verdad_de_referencia.py            # escribe el CSV e imprime el resumen
"""
import csv
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
SALIDA = os.path.join(BASE, "verdad_de_referencia.csv")

APTO = "si"
NO_APTO = "no"

# ---------------------------------------------------------------------------
# Fuentes de la literatura, tal y como las registro PLAN_TDP43_2026-09-20.md.
# No se anade ninguna cita que el proyecto no haya leido.
# ---------------------------------------------------------------------------
CITA_TDP43 = {
    "rTRD01": "François-Moutal 2019, ACS Chem Biol 14:2006 (Kd 89,4 uM, MST+NMR)",
    "nTRD22": "Mollasalehi 2020, ACS Chem Biol 15:2854 (Kd 145 uM)",
    "fragmento_1": "Nshogoza 2019, IJMS 20:3230 (CSP de 15N-RRM)",
    "fragmento_2": "Nshogoza 2019, IJMS 20:3230 (CSP de 15N-RRM)",
    "fragmento_3": "Nshogoza 2019, IJMS 20:3230 (CSP de 15N-RRM)",
    "PE859": "Kapsiani 2026, PMC12918951 (union medida en HEK + C. elegans)",
    "berberrubine": "Kapsiani 2026, PMC12918951 (union medida en HEK + C. elegans)",
    "AIM4": "Prasad 2016, Sci Rep 6:39490 (funcional, levadura)",
    "bis-ANS": "Babinchak 2020, Nat Commun 11:5574 (funcional, LLPS)",
    "XL20": "Gao 2026, Nature Aging 6:1667 (SPR micromolar aparente + CETSA); Asinex BDF34019555",
}


def leer_csv(ruta):
    if not os.path.exists(ruta):
        return []
    return list(csv.DictReader(open(ruta, encoding="utf-8", errors="ignore")))


def smiles_de(*rutas):
    """name (minusculas) -> SMILES, de los ficheros de activos del proyecto."""
    out = {}
    for ruta in rutas:
        for r in leer_csv(os.path.join(BASE, ruta)):
            n = (r.get("name") or r.get("ligand") or "").strip()
            s = (r.get("smiles") or "").strip()
            if n and s:
                out[n.lower()] = s
    return out


def veredicto_redocking():
    """codigo del ligando -> dict con el control de redocking de su estructura.

    Fuente unica: `redocking_trp32/modos_quimias_nuevas.csv`, que trae los 27 modos
    de cada una de las 11 estructuras con su RMSD a la pose cristalina, todos con
    el mismo protocolo (caja 24 A, exhaustividad 8, tres semillas).

    Se guardan las dos cifras y no una, porque significan cosas distintas:

      * `mejor_modo`     : el mejor RMSD de todos los modos. Mide si el acoplador
                           ES CAPAZ de encontrar la pose correcta.
      * `mejor_puntuado` : el RMSD del modo que gana por energia. Mide si la
                           funcion de puntuacion la prefiere — que es lo que
                           decide un cribado.

    El veredicto (criterio de 2 A) usa `mejor_modo`, igual que el informe del
    21 de septiembre. `mejor_puntuado` se reporta aparte porque es la causa.
    """
    ruta = os.path.join(BASE, "redocking_trp32", "modos_quimias_nuevas.csv")
    por_entrada = {}
    lig_de_entrada = {
        "4A7Q": "4MQ", "4A7G": "12I", "2WZ6": "ZO0", "2WZ0": "ZZT",
        "8GSQ": "K4I", "6A9O": "6B3", "5YTO": "946",
        "4A7T": "5FW", "4A7U": "ALE", "4A7V": "LDP", "4A7S": "5UD",
    }
    for r in leer_csv(ruta):
        try:
            rmsd = float(r["rmsd"])
            aff = float(r["afinidad"])
        except (KeyError, TypeError, ValueError):
            continue
        por_entrada.setdefault(r["entrada"], []).append((rmsd, aff))

    out = {}
    for entrada, modos in por_entrada.items():
        cod = lig_de_entrada.get(entrada)
        if cod is None:
            continue
        mejor_modo = min(m[0] for m in modos)
        mejor_puntuado = min(modos, key=lambda m: m[1])[0]
        out[cod] = {
            "entrada": entrada, "mejor_modo": mejor_modo,
            "mejor_puntuado": mejor_puntuado,
            "n_le2": sum(1 for m in modos if m[0] <= 2.0), "n_modos": len(modos),
            "veredicto": "RECUPERA" if mejor_modo <= 2.0 else "FALLA",
        }
    return out


def positivos_sod1(smis):
    """Los positivos de SOD1: co-cristalizados (con su estructura) + los 3 independientes."""
    reglas = {c[0]: c for c in [
        ("isoproterenol", "catecolamina", "Trp32", "5FW", "union_directa"),
        ("adrenalina", "catecolamina", "Trp32", "ALE", "union_directa"),
        ("dopamina", "catecolamina", "Trp32", "LDP", "union_directa"),
        ("5-fluorouridina", "nucleosido", "Trp32", "5UD", "union_directa"),
        ("quinazolina_12I", "quinazolina", "Trp32", "12I", "union_directa"),
        ("diazepanoquinazolina_4MQ", "quinazolina", "Trp32", "4MQ", "union_directa"),
        ("cf3quinazolina_ZO0", "quinazolina-CF3", "Trp32", "ZO0", "union_directa"),
        ("anilina_ZZT", "anilina", "Trp32", "ZZT", "union_directa"),
        ("paliperidona_K4I", "benzisoxazol-piperidina", "Trp32", "K4I", "union_directa"),
        ("Lig9_6B3", "fenantridinona", "Trp32", "6B3", "union_directa"),
        ("naftalenoaminoalcohol_946", "aminoalcohol naftenico", "Trp32", "946",
         "union_directa (co-cristal + SPR Kd 33,1 uM)"),
    ]}

    redock = veredicto_redocking()
    fuentes = {r["ligand"]: r for r in leer_csv(
        os.path.join(BASE, "controles_sod1_v4.csv"))}

    filas = []
    for nombre, (lig, quimia, sitio, cod, ensayo) in reglas.items():
        src = fuentes.get(lig, {})
        d = redock.get(cod)
        nota = ""
        if cod == "5UD":
            nota = ("su pose depositada NO es un minimo del potencial (F...N de 2,07 A "
                    "en una estructura de 1,06 A): no sirve como patron de medida")
        elif cod == "ZO0":
            nota = "recupera la pose en el modo 9, no en el que gana por energia"
        elif cod == "6B3":
            nota = "el ligando de 6A9O esta a medio resolver (24 de 33 atomos)"
        filas.append({
            "target": "SOD1", "ligand": lig, "smiles": smis.get(lig.lower(), ""),
            "quimia": quimia, "sitio": sitio, "tipo_ensayo": ensayo,
            "cita": src.get("cita", "Wright et al. 2013 Nat Commun 4:1758"),
            "estructura_pdb": src.get("estructura_pdb", ""), "ligando_pdb": cod,
            "redocking_mejor_modo_A": "" if d is None else "%.2f" % d["mejor_modo"],
            "redocking_mejor_puntuado_A": "" if d is None else "%.2f" % d["mejor_puntuado"],
            "redocking_modos_le2A": "" if d is None else "%d/%d" % (d["n_le2"], d["n_modos"]),
            "redocking_veredicto": "sin control" if d is None else d["veredicto"],
            "apto": APTO, "nota": nota,
        })

    # --- los tres independientes: aqui esta el valor real del conjunto ---
    alimento = [
        ("LCS-1", "piridazinona", "actividad celular; cribado previo del proyecto",
         "Wright et al. 2013 / cribado previo del proyecto", NO_APTO,
         "sin union medida: la actividad es celular y heredada"),
        ("PRG-A01", "cumarina", "actividad celular", "Woo et al. 2021 (PRG-A01)", NO_APTO,
         "sin union medida"),
        ("CHEMBL2165613", "pirazolona (serie)", "actividad celular (serie)",
         "serie ChEMBL2165601-2165614 / 1643541-57", NO_APTO,
         "es el representante UNICO de 18 analogos colapsados; no vale como "
         "evidencia de generalidad"),
    ]
    for lig, quimia, ensayo, cita, apto, nota in alimento:
        filas.append({
            "target": "SOD1", "ligand": lig, "smiles": smis.get(lig.lower(), ""),
            "quimia": quimia, "sitio": "n/d", "tipo_ensayo": ensayo, "cita": cita,
            "estructura_pdb": "", "ligando_pdb": "",
            "redocking_mejor_modo_A": "", "redocking_mejor_puntuado_A": "",
            "redocking_modos_le2A": "",
            "redocking_veredicto": "sin estructura", "apto": apto, "nota": nota,
        })
    return filas


def positivos_tdp43(smis):
    """Los de `PLAN_TDP43_2026-09-20.md`: 7 de union medida y 5 quimiotipos."""
    filas = [
        # (ligando, quimia, sitio, tipo_ensayo, apto, nota)
        ("rTRD01", "piperidinil-pirimidina", "RRM1 (y RRM2 por RMN)",
         "union_directa (Kd 89,4 uM, MST + NMR)", APTO, ""),
        ("nTRD22", "isoxazol-piperidina", "NTD, alosterico",
         "union_directa (Kd 145 uM, MST)", APTO, ""),
        # Los tres fragmentos de Nshogoza 2019 (IJMS 20:3230) existen SOLO como dibujo en
        # la Figura 1c: el articulo y su suplementario no dan nombre, CAS ni formula. Su
        # SMILES se leyo por OCR quimico (DECIMER 2.7.2) con `leer_figura_nshogoza.py`, se
        # comprobo en RDKit y se cotejo dibujandolo; esta en `fragmentos_nshogoza.csv`,
        # que es de donde lo coge este generador. Hasta el 25 de septiembre de 2026 ese
        # SMILES vivia SOLO en el informe del 23: una corrida sin el prepara los
        # fragmentos de cero o no los prepara, y con ellos fuera el bloque duro de TDP-43
        # baja de 0,734 (PASA) a 0,595 (SIN EVIDENCIA).
        ("fragmento_1", "fragmento", "RRM2 (G245, E246, H256, I257, S258)",
         "union_directa (CSP de 15N-RRM)", APTO,
         "SMILES leido del dibujo por OCR quimico (fragmentos_nshogoza.csv)"),
        ("fragmento_2", "fragmento", "RRM2 (mismos residuos)",
         "union_directa (CSP de 15N-RRM)", APTO,
         "SMILES leido del dibujo por OCR quimico (fragmentos_nshogoza.csv)"),
        ("fragmento_3", "fragmento", "RRM2 (mismos residuos)",
         "union_directa (CSP de 15N-RRM)", APTO,
         "SMILES leido del dibujo por OCR quimico (fragmentos_nshogoza.csv)"),
        ("PE859", "piridil-pirazol", "interfaz RRM1-RRM2",
         "union_directa (HEK + C. elegans)", APTO, ""),
        ("berberrubine", "bencilisoquinolina", "interfaz RRM1-RRM2",
         "union_directa (HEK + C. elegans)", APTO, ""),
        # los que NO valen para medir, y por que
        ("AIM4", "peptido", "CTD", "funcional (agregacion, levadura)", NO_APTO,
         "sin union medida"),
        ("bis-ANS", "sonda fluorescente", "CTD / LLPS",
         "funcional (separacion de fases)", NO_APTO, "sin union medida"),
        ("ketoconazole", "aznol-imidazol", "n/d",
         "funcional (reduce agregacion)", NO_APTO, "sin union medida"),
        # los que entraron por parecerse a la berberrubina: FUERA del fichero
        ("berberine", "bencilisoquinolina", "n/d", "similitud (clon de berberrubina)",
         NO_APTO, "entra por parecerse a la berberrubina, no por su propia medida"),
        ("sanguinarine", "bencilisoquinolina", "n/d", "similitud", NO_APTO,
         "entra por parecerse a la berberrubina"),
        ("cepharanthine", "bisbencilisoquinolina", "n/d", "similitud", NO_APTO,
         "entra por parecerse a la berberrubina"),
        ("epiberberine", "bencilisoquinolina", "n/d", "similitud", NO_APTO,
         "entra por parecerse a la berberrubina"),
        ("coptisine", "bencilisoquinolina", "n/d", "similitud", NO_APTO,
         "entra por parecerse a la berberrubina"),
        ("nitidine", "benzoquinolizinio", "n/d", "similitud", NO_APTO,
         "entra por parecerse a la berberrubina"),
        ("CongoRed", "diazo", "CTD", "funcional (sin union medida)", NO_APTO,
         "sin union medida"),
        # Encontrado el 25 sep de 2026 barriendo la literatura de union medida: XL20 es
        # un unido de verdad (SPR + CETSA, dependiente del Trp334), pero se une al CR
        # del dominio de baja complejidad (320-340), NO al bolsillo de RRM que acoplamos.
        # Por eso entra con la cita pero no apto: no puede contar como positivo de la
        # validacion de RRM sin mentir sobre el sitio. Sirve para la linea de agregacion.
        # El mismo dia se construyo el receptor del CR (analysis/_cr_receptor/, informe
        # INFORME_CR_XL20_XL23_2026-09-25.md) y se acoplo ahi con XL23: sigue sin ser
        # apto para RRM, y en el CR no hay control de redocking posible porque no existe
        # ningun ligando co-cristalizado con TDP-43 en todo el PDB.
        ("XL20", "adenina-aminociclohexanol",
         "CR del C-terminal (320-340, Trp334)",
         "union_directa (SPR micromolar aparente + CETSA)", NO_APTO,
         "union medida real pero de OTRO sitio: entra cuando se ataque el CR, no en la "
         "validacion de RRM. SMILES leido por OCR quimico (DECIMER 2.7.2) de la Figura "
         "Suplementaria 1a y contrastado con la descripcion del modelo de vision local "
         "(ver controles_tdp43_xl20.csv), sin contraste en bases de datos; XL21 y XL23 "
         "solo tienen actividad funcional, sin union medida"),
    ]
    filas_out = []
    for lig, quimia, sitio, ensayo, apto, nota in filas:
        filas_out.append({
            "target": "TDP43", "ligand": lig, "smiles": smis.get(lig.lower(), ""),
            "quimia": quimia, "sitio": sitio, "tipo_ensayo": ensayo,
            "cita": CITA_TDP43.get(lig, "recogido en PLAN_TDP43_2026-09-20.md"),
            "estructura_pdb": "", "ligando_pdb": "",
            "redocking_mejor_modo_A": "", "redocking_mejor_puntuado_A": "",
            "redocking_modos_le2A": "", "redocking_veredicto": "",
            "apto": apto, "nota": nota,
        })
    return filas_out


def positivos_fus(smis):
    filas = []
    for lig in ("Dehydroxymethylflazine", "CleroindicinC"):
        filas.append({
            "target": "FUS", "ligand": lig, "smiles": smis.get(lig.lower(), ""),
            "quimia": "", "sitio": "n/d", "tipo_ensayo": "sin tipo de ensayo",
            "cita": "sin literatura propia (ver ESTRATEGIA_AVANZAR_2026-09-20.md)",
            "estructura_pdb": "", "ligando_pdb": "",
            "redocking_mejor_modo_A": "", "redocking_mejor_puntuado_A": "",
            "redocking_modos_le2A": "", "redocking_veredicto": "", "apto": NO_APTO,
            "nota": "FUS no tiene ligandos propios en la literatura: no es evaluable",
        })
    return filas


def main():
    # `fragmentos_nshogoza.csv` entra aqui a proposito: son SMILES que no estan en ningun
    # otro fichero de activos porque en la literatura solo existen dibujados, y sin ellos
    # la validacion de TDP-43 pierde tres positivos.
    smis = smiles_de("activos_sod1_v2.csv", "activos_tdp43.csv",
                     "activos_tdp43_v2.csv", "activos_fus.csv", "controles_sod1_v4.csv",
                     "controles_tdp43_xl20.csv", "fragmentos_nshogoza.csv")

    filas = positivos_sod1(smis) + positivos_tdp43(smis) + positivos_fus(smis)

    campos = ["target", "ligand", "smiles", "quimia", "sitio", "tipo_ensayo",
              "cita", "estructura_pdb", "ligando_pdb", "redocking_mejor_modo_A",
              "redocking_mejor_puntuado_A", "redocking_modos_le2A",
              "redocking_veredicto", "apto", "nota"]
    with open(SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for r in filas:
            w.writerow(r)

    print("=" * 78)
    print("VERDAD DE REFERENCIA — %d entradas" % len(filas))
    print("=" * 78)
    for target in ("SOD1", "TDP43", "FUS"):
        grupo = [r for r in filas if r["target"] == target]
        aptos = [r for r in grupo if r["apto"] == APTO]
        quimias = {r["quimia"] for r in aptos if r["quimia"]}
        print("")
        print("%s: %d entradas, %d aptos, %d quimias distintas de positivos"
              % (target, len(grupo), len(aptos), len(quimias)))
        for r in grupo:
            marca = "APTO  " if r["apto"] == APTO else "no    "
            print("   %s%-22s %-26s %s" % (marca, r["ligand"][:22], r["quimia"][:26],
                                           r["nota"][:46]))
        if len(quimias) < 2 and aptos:
            print("   !! con menos de dos quimias distintas no se puede validar nada")
    print("")
    print("guardado: %s" % SALIDA)
    return 0


if __name__ == "__main__":
    sys.exit(main())
