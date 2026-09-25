#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""La prueba decisiva de la anomalia: llevar la pose NATIVA al receptor del cribado.

QUE SE PREGUNTA
---------------
En el receptor del cribado (gpu_dock/SOD1.pdbqt, que es 1HL5 con todo su
ensamblaje y sus aguas) las tres quinazolinas dan:

    4MQ  -10,0    ZO0  -9,8    (diazepano, 7 anillos)
    12I   -5,0                 (piperazina, 6 anillos)

Cinco kcal/mol por un CH2. Ya se ha comprobado que:

  * no es ruido: se reproduce con semillas fijas (dispersion 0,02-0,38 kcal/mol);
  * no son las aguas: quitando las 1.763 aguas cristalograficas del receptor el
    orden no cambia (-10,4 / -10,9 / -5,5).

Lo que queda por saber es si esas poses -10 tienen algo que ver con la pose que
el ligando adopta de verdad en el cristal, o si son una pose distinta que la
funcion de puntuacion premia dentro de esa caja.

COMO SE RESPONDE
----------------
Se coge la pose cristalina de cada ligando (de su propia estructura: 4A7Q, 2WZ6,
4A7G), se superpone SU cadena A sobre la cadena A del receptor del cribado usando
los atomos del propio bolsillo de Trp32 (residuos 21, 22, 23, 28, 29, 30, 32, 99,
100 por numero y nombre de atomo), y se puntua la pose ya colocada dentro del
receptor del cribado con `vina --score_only`. Ese numero se compara con el -10,0
que da el acoplamiento.

  * Si la pose nativa tambien puntua cerca de -10, la funcion prefiere la pose
    nativa y el -10 es real: el ligando SI se une ahi, y el problema es que el
    acoplamiento no la encuentra.
  * Si la pose nativa puntua cerca de -5, el -10 es una pose artificial: la
    funcion prefiere otra colocacion, y el «positivo» de la v4 es un artefacto.

Uso: python pose_nativa_en_cribado.py
Salida: pose_nativa_en_cribado.log y pose_nativa_en_cribado.csv
"""
import csv
import os
import re
import subprocess
import sys

import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
GPU = os.path.abspath(os.path.join(BASE, "..", "gpu_dock"))
CRIBAJE = os.path.join(GPU, "SOD1.pdbqt")
VINA = r"C:\Users\Fredy\masive-als\tools\vina.exe"
CENTRO_CRIBAJE = (46.5, 80.0, 73.3)
TAMANO = 22

RESIDUOS_BOLSILLO = {"19", "20", "21", "22", "23", "28", "29", "30", "32",
                     "33", "96", "97", "98", "99", "100", "101", "135"}

SISTEMAS = {
    "4MQ": ("4A7Q", "4MQ", "A", "1000", "diazepanoquinazolina", -10.07),
    "ZO0": ("2WZ6", "ZO0", "F", "206", "cf3-diazepanoquinazolina", -9.83),
    "12I": ("4A7G", "12I", "A", "1001", "piperazinaquinazolina", -4.89),
}


def log(m):
    print(m, flush=True)


def leer_atomos(ruta, prefijos=("ATOM", "HETATM")):
    out = []
    for l in open(ruta, encoding="utf-8", errors="ignore"):
        if not l.startswith(prefijos):
            continue
        tipo = l.rsplit(None, 1)[-1].upper() if len(l) > 70 else ""
        out.append({
            "res": l[17:20].strip(), "cad": l[21], "num": l[22:27].strip(),
            "atom": l[12:16].strip(),
            "elem": (l[76:78].strip() or _elemento(l[12:16])).upper(),
            "xyz": np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]),
            "hidrogeno": tipo in ("H", "HD", "HS", "D", "DD")
            or (len(l) > 78 and l[76:78].strip() == "H"),
        })
    return out


def _elemento(nombre):
    """Simbolo quimico a partir del nombre de atomo del PDB (C1, N12, OE1...)."""
    letras = "".join(c for c in nombre if c.isalpha())
    for dos in ("CL", "BR", "NA", "MG", "ZN", "CU", "FE", "MN"):
        if letras.upper().startswith(dos):
            return dos
    return (letras[:1] or "C").upper()


def kabsch(P, Q):
    """Rotacion + traslacion que lleva P sobre Q (sin escalar)."""
    cp, cq = P.mean(0), Q.mean(0)
    H = (P - cp).T @ (Q - cq)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    D = np.diag([1.0, 1.0, d])
    R = Vt.T @ D @ U.T
    return R, cq - R @ cp


def marco(atomos, cadena):
    """Atomos del bolsillo de Trp32 en una cadena, indexados por residuo:atomo."""
    out = {}
    for a in atomos:
        if a["cad"] != cadena or a["hidrogeno"]:
            continue
        if a["num"] in RESIDUOS_BOLSILLO:
            out["%s:%s" % (a["num"], a["atom"])] = a["xyz"]
    return out


def superponer(marco_a, marco_b):
    """Transformacion que lleva el bolsillo de A sobre el de B (mismos atomos)."""
    comunes = sorted(set(marco_a) & set(marco_b))
    if len(comunes) < 20:
        return None, None, len(comunes)
    P = np.array([marco_a[k] for k in comunes])
    Q = np.array([marco_b[k] for k in comunes])
    R, t = kabsch(P, Q)
    desvio = np.sqrt((((P @ R.T) + t - Q) ** 2).sum(1).mean())
    return (R, t), float(desvio), len(comunes)


def escribir_pdb(atomos, destino):
    with open(destino, "w", encoding="utf-8") as f:
        for i, (nombre, el, xyz) in enumerate(atomos, 1):
            f.write("HETATM%5d %-4s %3s A   1    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
                    % (i, nombre[:4], "LIG", xyz[0], xyz[1], xyz[2], el))
        f.write("END\n")


def afinidad_stdout(txt):
    vals = re.findall(r"Estimated Free Energy of Binding\s*:\s*([-+]?\d+\.\d+)", txt or "")
    return float(vals[-1]) if vals else None


def puntuar(receptor, ligando_pdbqt, centro, lado=TAMANO):
    r = subprocess.run([VINA, "--receptor", receptor, "--ligand", ligando_pdbqt,
                        "--score_only",
                        "--center_x", "%.3f" % centro[0],
                        "--center_y", "%.3f" % centro[1],
                        "--center_z", "%.3f" % centro[2],
                        "--size_x", str(lado), "--size_y", str(lado),
                        "--size_z", str(lado)],
                       capture_output=True, text=True, timeout=1800)
    return afinidad_stdout((r.stdout or "") + (r.stderr or ""))


def coord_pesados(path, solo_primero=True):
    out = []
    for l in open(path, encoding="utf-8", errors="ignore"):
        if solo_primero and out and l.startswith("ENDMDL"):
            break
        if not l.startswith(("ATOM  ", "HETATM")):
            continue
        tipo = l.rsplit(None, 1)[-1].upper()
        if tipo in ("H", "HD", "HS", "D", "DD"):
            continue
        out.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
    return np.array(out)


def minimizar_local(receptor, ligando_pdbqt, centro, salida, lado=TAMANO):
    """Minimizacion local desde la pose dada: dice si la pose AGUANTA ahi.

    Con la pose nativa ya colocada, si al minimizar se queda casi donde estaba y
    su energia es parecida a la mejor pose acoplada, entonces ese modo de union
    es un minimo real de la funcion dentro de ESE receptor. Si se va a otro
    sitio, la funcion no lo sostiene.
    """
    r = subprocess.run([VINA, "--receptor", receptor, "--ligand", ligando_pdbqt,
                        "--local_only", "--out", salida,
                        "--center_x", "%.3f" % centro[0],
                        "--center_y", "%.3f" % centro[1],
                        "--center_z", "%.3f" % centro[2],
                        "--size_x", str(lado), "--size_y", str(lado),
                        "--size_z", str(lado)],
                       capture_output=True, text=True, timeout=1800)
    e = afinidad_stdout((r.stdout or "") + (r.stderr or ""))
    deriva = None
    if os.path.exists(salida) and os.path.getsize(salida) > 100:
        a = coord_pesados(ligando_pdbqt)
        b = coord_pesados(salida)
        if len(a) == len(b) and len(a):
            deriva = float(np.sqrt(((a - b) ** 2).sum(1).mean()))
    return e, deriva


def main():
    from rdkit import Chem
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    atomos_cribaje = leer_atomos(CRIBAJE)
    marco_b = marco(atomos_cribaje, "A")
    log("receptor del cribado: %d atomos | bolsillo Trp32 de la cadena A: %d atomos"
        % (len(atomos_cribaje), len(marco_b)))
    outdir = os.path.join(BASE, "pose_nativa")
    os.makedirs(outdir, exist_ok=True)

    filas, lineas = [], []
    for lig, (pdb_id, codigo, cadena, resseq, nombre, afin_dock) in SISTEMAS.items():
        pdb = os.path.join(BASE, "trp32_cristal", "pdb", pdb_id + ".pdb")
        atomos = leer_atomos(pdb)
        marco_a = marco(atomos, cadena)
        (R, t), desvio, n = superponer(marco_a, marco_b)
        if R is None:
            log("%s: no hay suficientes atomos comunes para superponer (%d)" % (lig, n))
            continue
        log("\n=== %s (%s, %s%s%s) ===" % (lig, nombre, cadena, resseq, ""))
        log("  superposicion del bolsillo: %d atomos, desvio %.2f A" % (n, desvio))

        # ligando del cristal, en su copia del bolsillo
        lig_at = [a for a in atomos
                  if a["res"] == codigo and a["cad"] == cadena
                  and a["num"] == resseq and not a["hidrogeno"]]
        if not lig_at:
            log("  no se encontro el ligando en el PDB")
            continue
        antes = np.array([a["xyz"] for a in lig_at])
        despues = (antes @ R.T) + t
        log("  ligando: %d atomos | centroide movido %.2f A"
            % (len(despues), float(np.linalg.norm(despues.mean(0) - antes.mean(0)))))
        log("  centroide colocado: %.1f %.1f %.1f  (caja del cribado: %.1f %.1f %.1f)"
            % (*despues.mean(0), *CENTRO_CRIBAJE))
        log("  distancia centroide-caja %.2f A"
            % float(np.linalg.norm(despues.mean(0) - np.array(CENTRO_CRIBAJE))))

        # PDBQT de la pose colocada, con los ordenes de enlace del ligando ideal
        crudo = os.path.join(outdir, "%s_nativa.pdb" % lig)
        escribir_pdb([(a["atom"], a["elem"], xyz)
                      for a, xyz in zip(lig_at, despues)], crudo)
        ideal = Chem.MolFromMolFile(os.path.join(BASE, "redocking_trp32", "ligands",
                                                 "%s_ideal.sdf" % codigo))
        nativa = Chem.MolFromPDBBlock(open(crudo).read(), removeHs=False,
                                      proximityBonding=True)
        pdbqt = os.path.join(outdir, "%s_nativa.pdbqt" % lig)
        ok = False
        if ideal is not None and nativa is not None \
                and ideal.GetNumAtoms() == nativa.GetNumAtoms():
            from rdkit.Chem import AllChem
            m = AllChem.AssignBondOrdersFromTemplate(ideal, Chem.RemoveHs(nativa))
            setup = MoleculePreparation(rigid_macrocycles=True).prepare(
                Chem.AddHs(m, addCoords=True))[0]
            txt, ok, _ = PDBQTWriterLegacy.write_string(setup)
            if ok:
                open(pdbqt, "w", encoding="utf-8").write(txt)
        if not ok:
            log("  no se pudo preparar el PDBQT de la pose nativa")
            continue

        fijo = puntuar(CRIBAJE, pdbqt, CENTRO_CRIBAJE)
        propio = puntuar(CRIBAJE, pdbqt, despues.mean(0))
        log("  E(pose nativa) en el receptor del cribado, caja del proyecto: %s kcal/mol"
            % ("%.2f" % fijo if fijo is not None else "no leida"))
        log("  E(pose nativa) con la caja centrada en la propia pose:       %s kcal/mol"
            % ("%.2f" % propio if propio is not None else "no leida"))
        log("  E(mejor pose acoplada por Vina en el mismo receptor):         %.2f kcal/mol"
            % afin_dock)
        if fijo is not None:
            log("  -> la pose nativa es %.2f kcal/mol PEOR que la pose acoplada"
                % (fijo - afin_dock))

        out_min = os.path.join(outdir, "%s_nativa_min.pdbqt" % lig)
        e_min, deriva = minimizar_local(CRIBAJE, pdbqt, despues.mean(0), out_min)
        log("  E(pose nativa MINIMIZADA en el receptor del cribado):        %s kcal/mol "
            "| deriva %.2f A"
            % ("%.2f" % e_min if e_min is not None else "no leida",
               deriva if deriva is not None else float("nan")))
        filas.append({"ligando": lig, "nombre": nombre, "estructura": pdb_id,
                      "copias": "%s%s" % (cadena, resseq),
                      "desvio_superposicion_A": round(desvio, 2),
                      "atomos_superpuestos": n,
                      "dist_centroide_caja_A": round(float(np.linalg.norm(
                          despues.mean(0) - np.array(CENTRO_CRIBAJE))), 2),
                      "E_nativa_caja_proyecto": fijo, "E_nativa_caja_propia": propio,
                      "E_nativa_minimizada": e_min,
                      "deriva_minimizacion_A": (round(deriva, 2)
                                                if deriva is not None else None),
                      "E_dock": afin_dock,
                      "delta_nativa_menos_dock": (round(fijo - afin_dock, 2)
                                                  if fijo is not None else None)})
        lineas.append("%-4s %-26s desvio %.2f A | E(nativa) %s | E(dock) %.2f | dE %s"
                      % (lig, nombre, desvio,
                         "%.2f" % fijo if fijo is not None else "n/d", afin_dock,
                         "%+.2f" % (fijo - afin_dock) if fijo is not None else "n/d"))

    texto = "\n".join(lineas)
    with open(os.path.join(BASE, "pose_nativa_en_cribado.log"), "w", encoding="utf-8") as f:
        f.write(texto + "\n")
    campos = list(filas[0].keys()) if filas else ["ligando"]
    with open(os.path.join(BASE, "pose_nativa_en_cribado.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        for fila in filas:
            w.writerow(fila)
    print("\n" + texto)
    log("guardado: pose_nativa_en_cribado.csv / .log")
    return 0


if __name__ == "__main__":
    sys.exit(main())
