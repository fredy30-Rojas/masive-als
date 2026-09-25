#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Receptor del CR del C-terminal de TDP-43 (320-340, Trp334) para acoplar XL20 y XL23.

POR QUE EXISTE
--------------
El unico positivo con union medida de TDP-43 fuera del bolsillo de ARN es **XL20**,
y su sitio no es una suposicion: en el articulo (Gao et al. 2026, Nature Aging 6:1667)
borrar el **CR (320-340)** hace desaparecer la union y mutar el **Trp334** la baja
mucho. Para acoplar ahi hace falta un receptor del CR, y el proyecto no lo tenia:
la caja que usaba era la del bolsillo RRM1-RRM2 (receptor 4BS2).

QUE HACE
--------
1. **Elige la estructura del CR con numeros, no a ojo.** Recorre las entradas del PDB
   que contienen esa region, pide a PDBe su secuencia por entidad y localiza el
   fragmento dentro de la TDP-43 humana (UniProt Q13148) por un ancla interna
   (`MNFGAFSINP`), que es lo que fija la numeracion real de cada construccion.
   Requisitos: cubrir 320-340 entero y traer el anillo indol del Trp334 resuelto.
   Entre las que cumplen se queda **la mas corta**, porque cuanto menos proteina
   haya fuera del CR, menos se le puede atribuir a otra cosa lo que toque el ligando.
2. **Separa los modelos** del conjunto de RMN: un CR intrinsecamente desordenado no
   es una foto, es una ensenada de conformaciones, y acoplar en una sola seria
   elegir a dedo. Cada modelo se acopla por separado.
3. **Prepara cada modelo con la receta de acoplamiento del proyecto** (meeko,
   `mk_prepare_receptor -p -j -a`, la misma con la que los controles de redocking
   del Trp32 dan RMSD <= 2 A). No se pasa por PDBFixer a proposito: la receta
   validada no lo usa, y ademas PDBFixer coloca hidrogenos de forma no determinista
   (medido en RESCORING_WSL_CUDA_2026-09-11.md).
4. **Define la caja sobre el anillo indol del Trp334** de cada modelo (centroide de
   CG-CD1-CD2-NE1-CE2-CE3-CZ2-CZ3-CH2, caja de 22 A como el resto del proyecto) y
   comprueba que la caja ve el CR: que residuos de 320-340 caen dentro.

Uso:
    python analysis/construir_receptor_cr.py

Salida en analysis/_cr_receptor/:
    pdb/<entrada>.pdb          la estructura descargada (cache)
    modelos/cr_modelo<n>.pdb   la cadena A, solo proteina, un modelo por fichero
    modelos/cr_modelo<n>.pdbqt el receptor de ese modelo, ya preparado
    candidatos_cr.csv          el barrido de estructuras del CR y por que se elige esta
    caja.json                  centro y tamano de la caja de cada modelo
    receptor_cr.txt            el informe legible de lo que se hizo y se midio
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import sys
import urllib.request

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
SALIDA = os.path.join(BASE, "_cr_receptor")
PDBS = os.path.join(SALIDA, "pdb")
MODELOS = os.path.join(SALIDA, "modelos")

sys.path.insert(0, os.path.join(BASE, "redocking_trp32"))
import redock_trp32 as R  # noqa: E402  (PYTHON, MK_RECEPTOR y la receta del proyecto)

UNIPROT = "Q13148"
CR = (320, 340)                       # la region conservada del articulo
TRP334 = 334
ANCLA = "MNFGAFSINP"                  # ancla interna para fijar la numeracion real
ANILLO_TRP = ("CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2")
TAMANO_CAJA = 22                      # el mismo lado que usa el proyecto

# Entradas del PDB que contienen el CR. No es una busqueda ciega: son las que
# aparecen al buscar en el PDB la secuencia del CR; se listan aqui para que el
# barrido sea reproducible y se pueda anadir o quitar una sin tocar el codigo.
CANDIDATAS = ["2N2C", "2N3X", "2N4G", "2N4H", "6N37", "6N3A", "6N3B", "6N3C",
              "7KWZ", "7PY2", "7Q3U", "8CG3", "8CGG", "8CGH", "8QX9", "8QXA",
              "8QXB", "9FOF", "9FOR"]


LINEAS = []   # todo lo que se imprime, para dejar ademas el fichero legible


def log(m):
    print(m, flush=True)
    LINEAS.append(m)


def texto(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": "masive-als/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", errors="replace")


def json_de(url):
    return json.loads(texto(url))


def uniprot_seq():
    fasta = texto("https://rest.uniprot.org/uniprotkb/%s.fasta" % UNIPROT)
    return "".join(l.strip() for l in fasta.splitlines() if not l.startswith(">"))


def entidades(pid):
    """Entidades de una entrada segun PDBe: secuencia y cadenas."""
    try:
        d = json_de("https://www.ebi.ac.uk/pdbe/api/pdb/entry/entities/%s" % pid.lower())
    except Exception as e:  # noqa: BLE001
        log("   %s: PDBe no respondio (%s)" % (pid, e))
        return []
    return [e for e in d.get(pid.lower(), []) if (e.get("sequence") or "")]


def numero_del_constructo(seq, referencia):
    """En que residuo de la TDP-43 humana empieza esta construccion.

    Se localiza el ancla en la referencia y en el constructo; la distancia entre las
    dos posiciones da el desplazamiento. Asi funciona tambien con los mutantes
    (G335D en 2N4G, Q343R en 2N4H), donde una comparacion literal fallaria.
    """
    en_ref = referencia.find(ANCLA)
    en_con = seq.find(ANCLA)
    if en_ref < 0 or en_con < 0:
        return None
    return en_ref - en_con + 1


def descargar(pid):
    os.makedirs(PDBS, exist_ok=True)
    destino = os.path.join(PDBS, pid.lower() + ".pdb")
    if os.path.exists(destino) and os.path.getsize(destino) > 1000:
        return destino
    try:
        t = texto("https://files.rcsb.org/download/%s.pdb" % pid, timeout=180)
    except Exception as e:  # noqa: BLE001
        log("   %s: no se pudo descargar (%s)" % (pid, e))
        return None
    if "ATOM" not in t:
        return None
    with open(destino, "w", encoding="utf-8") as f:
        f.write(t)
    return destino


def trocear_modelos(pdb):
    """[(numero_de_modelo, lineas)] de un PDB, con o sin registros MODEL."""
    modelos, actual, n = [], [], 0
    for l in open(pdb, encoding="utf-8", errors="ignore"):
        if l.startswith("MODEL"):
            actual, n = [], n + 1
        elif l.startswith("ENDMDL"):
            if actual:
                modelos.append((n or 1, actual))
            actual = []
        elif l.startswith("ATOM"):
            actual.append(l)
    if actual:
        modelos.append((n or 1, actual))
    return modelos


def atomos_de(lineas, residuo=None, cadena="A"):
    out = []
    for l in lineas:
        if l[21] != cadena:
            continue
        if residuo is not None and l[22:26].strip() != str(residuo):
            continue
        out.append((l[17:20].strip(), l[12:16].strip(),
                    np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return out


def escribir_modelo(lineas, residuo_inicial, residuo_final, destino):
    """La cadena A del modelo, solo aminoacidos dentro del rango del CR, sin HETATM.

    Se resuelven las conformaciones alternas quedando con la de mayor ocupacion (o
    sin altloc), igual que `redock_trp32.extraer_receptor`, para que meeko no tenga
    que elegir y no reconstruya cadenas laterales del sitio.
    """
    elegidas, orden = {}, []
    for l in lineas:
        if not l.startswith("ATOM") or l[21] != "A":
            continue
        try:
            num = int(l[22:26])
        except ValueError:
            continue
        if num < residuo_inicial or num > residuo_final:
            continue
        clave = (l[22:27], l[12:16])
        try:
            occ = float(l[54:60])
        except ValueError:
            occ = 0.0
        valor = (0 if l[16] == " " else 1, -occ, l[16])
        if clave not in elegidas:
            orden.append(clave)
            elegidas[clave] = (valor, l)
        elif valor < elegidas[clave][0]:
            elegidas[clave] = (valor, l)
    with open(destino, "w", encoding="utf-8") as f:
        for clave in orden:
            f.write(elegidas[clave][1])
        f.write("END\n")
    return sum(1 for _ in orden)


def barrido(referencia):
    """Que entradas del PDB traen el CR entero y con el Trp334 resuelto."""
    filas = []
    for pid in CANDIDATAS:
        entrada = {"pdb": pid}
        ents = entidades(pid)
        if not ents:
            filas.append({**entrada, "cubre_el_cr": False, "nota": "sin entidades en PDBe"})
            continue
        e = max(ents, key=lambda x: len(x.get("sequence") or ""))
        seq = e["sequence"]
        inicio = numero_del_constructo(seq, referencia)
        if inicio is None:
            filas.append({**entrada, "cubre_el_cr": False,
                          "nota": "la secuencia no lleva el ancla %s" % ANCLA})
            continue
        fin = inicio + len(seq) - 1
        cubre = inicio <= CR[0] and fin >= CR[1]
        indice = TRP334 - inicio
        trp = indice < len(seq) and seq[indice] == "W"
        campo = {"pdb": pid, "residuos": "%d-%d" % (inicio, fin), "n_residuos": len(seq),
                 "cadenas": ",".join(e.get("in_chains") or []),
                 "trp334_en_secuencia": trp, "cubre_el_cr": cubre,
                 "nota": "" if cubre else "no cubre 320-340"}
        # ¿y esta resuelto el indol en las coordenadas?
        p = descargar(pid)
        if p:
            modelos = trocear_modelos(p)
            campo["n_modelos"] = len(modelos)
            if modelos:
                primeros = atomos_de(modelos[0][1], TRP334)
                atomos = {a[1] for a in primeros}
                campo["trp334_resuelto"] = all(n in atomos for n in ANILLO_TRP)
                presentes = sorted({int(l[22:26]) for l in modelos[0][1]
                                    if l.startswith("ATOM") and l[21] == "A"
                                    and CR[0] <= int(l[22:26]) <= CR[1]})
                campo["residuos_del_cr_presentes"] = len(presentes)
                campo["cr_completo_en_coordenadas"] = len(presentes) == CR[1] - CR[0] + 1
        filas.append(campo)
    return filas


def centro_anillo(lineas):
    """Centroide del anillo indol del Trp334 y sus atomos, en un modelo."""
    pts = [xyz for (rn, nom, xyz) in atomos_de(lineas, TRP334) if nom in ANILLO_TRP]
    if not pts:
        return None, []
    return np.mean(np.array(pts), axis=0), pts


def helicidad(lineas, primero, ultimo):
    """Distancia CA(i)-CA(i+4) dentro del CR: es la medida de si aquello es helice.

    Se mide en vez de darlo por supuesto porque la construccion se eligio justo por
    ser corta, y una construccion corta puede no plegarse como la region entera.
    Una alpha-helice da 5,4-6,5 A; un tramo extendido, mucho mas.
    """
    ca = {}
    for l in lineas:
        if l.startswith("ATOM") and l[12:16].strip() == "CA" and l[21] == "A":
            res = int(l[22:26])
            if primero - 6 <= res <= ultimo + 6:
                ca[res] = np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    d = [float(np.linalg.norm(ca[r] - ca[r + 4]))
         for r in range(primero, ultimo - 3) if r in ca and r + 4 in ca]
    if not d:
        return 0, None, None, None
    return len(d), float(np.mean(d)), float(min(d)), float(max(d))


def main():
    os.makedirs(MODELOS, exist_ok=True)
    referencia = uniprot_seq()
    log("TDP-43 humana (UniProt %s): %d residuos | residuo %d = %s"
        % (UNIPROT, len(referencia), TRP334, referencia[TRP334 - 1]))
    log("CR del articulo = %d-%d: %s" % (CR[0], CR[1], referencia[CR[0] - 1:CR[1]]))
    log("")

    filas = barrido(referencia)
    campos = ["pdb", "residuos", "n_residuos", "cadenas", "n_modelos", "cubre_el_cr",
              "trp334_en_secuencia", "trp334_resuelto", "residuos_del_cr_presentes",
              "cr_completo_en_coordenadas", "nota"]
    with open(os.path.join(SALIDA, "candidatos_cr.csv"), "w", encoding="utf-8",
              newline="") as f:
        w = csv.DictWriter(f, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        w.writerows(filas)

    log("Barrido de estructuras del CR (analysis/_cr_receptor/candidatos_cr.csv):")
    log("  %-6s %-10s %-5s %-6s %-6s %-9s %s"
        % ("pdb", "residuos", "nres", "mods", "indol", "cr_compl", "nota"))
    aptas = []
    for r in filas:
        log("  %-6s %-10s %-5s %-6s %-6s %-9s %s"
            % (r["pdb"], r.get("residuos", "-"), r.get("n_residuos", "-"),
               r.get("n_modelos", "-"), r.get("trp334_resuelto", "-"),
               r.get("cr_completo_en_coordenadas", "-"), r.get("nota", "")))
        if (r.get("cubre_el_cr") and r.get("trp334_resuelto")
                and r.get("cr_completo_en_coordenadas")):
            aptas.append(r)
    if not aptas:
        log("ninguna estructura cumple los requisitos")
        return 1
    # La mas corta: menos proteina fuera del CR que pueda llevarse los contactos.
    elegida = min(aptas, key=lambda r: r["n_residuos"])
    log("")
    log("ELEGIDA: %s (%s residuos, %s). Criterio: cubre 320-340 entero, con el anillo"
        % (elegida["pdb"], elegida["n_residuos"], elegida["residuos"]))
    log("del Trp334 resuelto, y es la construccion MAS CORTA que lo hace: cuanto menos")
    log("haya fuera del CR, menos se le puede atribuir a otra cosa lo que toque el ligando.")

    pdb = os.path.join(PDBS, elegida["pdb"].lower() + ".pdb")
    modelos = trocear_modelos(pdb)
    inicio, fin = [int(x) for x in elegida["residuos"].split("-")]
    log("")
    log("%d modelos en %s (RMN: el CR es una ensenada de conformaciones, no una foto)."
        % (len(modelos), os.path.basename(pdb)))

    cajas = {}
    helices = []
    for n, lineas in modelos:
        crudo = os.path.join(MODELOS, "cr_modelo%d.pdb" % n)
        base = os.path.join(MODELOS, "cr_modelo%d" % n)
        pdbqt = base + ".pdbqt"
        n_atomos = escribir_modelo(lineas, inicio, fin, crudo)
        centro, pts = centro_anillo(lineas)
        if not os.path.exists(pdbqt):
            # -a: permite residuos que meeko no empareja con plantilla (terminales de
            # una construccion recortada). Misma receta que los controles de redocking.
            r = subprocess.run([R.PYTHON, R.MK_RECEPTOR, "--read_pdb", crudo,
                                "-o", base, "-p", "-j", "-a", "--default_altloc", "A"],
                               capture_output=True, text=True, timeout=900)
            if not os.path.exists(pdbqt):
                log("  modelo %d: FALLO meeko: %s"
                    % (n, ((r.stdout or "") + (r.stderr or ""))[-300:]))
                continue
        pesados = sum(1 for l in open(pdbqt, encoding="utf-8", errors="ignore")
                      if l.startswith("ATOM"))
        # que ve la caja: cuantos residuos del CR caen dentro
        vista = sorted({res for (res, nom, xyz) in
                        [(l[22:26].strip(), l[12:16].strip(),
                          np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
                         for l in lineas if l.startswith("ATOM") and l[21] == "A"]
                        if float(np.linalg.norm(xyz - centro)) <= TAMANO_CAJA / 2.0})
        n_helices, hel_media, hel_min, hel_max = helicidad(lineas, CR[0], CR[1])
        helices.append(hel_media)
        cajas[str(n)] = {"centro": [round(float(x), 3) for x in centro],
                         "tamano": TAMANO_CAJA,
                         "pdbqt": os.path.relpath(pdbqt, BASE).replace("\\", "/"),
                         "residuos_dentro_de_la_caja": len(vista),
                         "primero": vista[0] if vista else None,
                         "ultimo": vista[-1] if vista else None,
                         "helice_ca_i_a_i4_media": (round(hel_media, 2)
                                                    if hel_media is not None else None)}
        log("  modelo %d: %d atomos de proteina -> %d atomos en el receptor | caja en"
            " (%.1f, %.1f, %.1f) ve %d residuos (%s-%s)"
            % (n, n_atomos, pesados, centro[0], centro[1], centro[2], len(vista),
               vista[0] if vista else "-", vista[-1] if vista else "-"))
        log("    helice del CR: CA(i)-CA(i+4) media %.2f A (de %.1f a %.1f, %d medidas)"
            % (hel_media, hel_min, hel_max, n_helices))

    # --- medidas de conjunto que cita el informe, para que salgan de un script ---
    centros = np.array([c["centro"] for c in cajas.values()])
    if len(centros) > 1:
        pares = np.sqrt(((centros[:, None, :] - centros[None, :, :]) ** 2).sum(-1))
        separacion = float(pares.max())
    else:
        separacion = 0.0
    log("")
    log("%d de las %d estructuras del CR cumplen los requisitos (cubren 320-340 entero,"
        " con el anillo del Trp334 resuelto)." % (len(aptas), len(filas)))
    log("El centro del anillo del Trp334 se separa hasta %.2f A entre modelos, asi que cada"
        " modelo lleva su caja." % separacion)
    if helices:
        log("Helice del CR en los %d modelos: CA(i)-CA(i+4) media %.2f A (de %.2f a %.2f)."
            % (len(helices), float(np.mean(helices)), float(np.min(helices)),
               float(np.max(helices))))

    with open(os.path.join(SALIDA, "caja.json"), "w", encoding="utf-8") as f:
        json.dump({"receptor": os.path.basename(pdb), "pdb": elegida["pdb"],
                   "residuos": elegida["residuos"], "trp334": TRP334,
                   "anillo": list(ANILLO_TRP),
                   "candidatas_que_cumplen": len(aptas),
                   "candidatas_probadas": len(filas),
                   "separacion_maxima_del_centro_entre_modelos": round(separacion, 2),
                   "helice_ca_i_a_i4_media": round(float(np.mean(helices)), 2)
                   if helices else None,
                   "modelos": cajas}, f, indent=2)
    log("")
    log("cajas y receptores anotados en analysis/_cr_receptor/caja.json")
    with open(os.path.join(SALIDA, "receptor_cr.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LINEAS) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
