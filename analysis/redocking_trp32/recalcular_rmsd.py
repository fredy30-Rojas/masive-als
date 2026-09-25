#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recalcula TODOS los controles con el medidor de RMSD corregido (20 sep 2026).

No vuelve a acoplar nada: lee los PDBQT que ya estan en `out/` y remide. El
acoplamiento es determinista con semilla fija y lo unico que estaba mal era la
forma de emparejar los atomos, asi que repetirlo no aportaria nada.

Que recalcula:
  1. Control 1 (`redock_trp32.py`): cajas de 24 y 18 A, exhaustividad 8 y 32.
  2. Control 2 (`protocolo_alternativo.py`): ligando flexible, vina y vinardo.
  3. Control 2 con ligando RIGIDO: aqui el acoplamiento conserva el orden de
     atomos, asi que hay un emparejamiento directo; se anade la correccion de
     simetria (los dos metilos del isoproterenol son equivalentes) copiando la
     cabecera REMARK SMILES del ligando de entrada.
  4. Control 3 (`redock_flexible.py`): receptor flexible, con el RMSD en el
     sitio corregido y el alineado (Kabsch) ya calculado en `resumen_controles.py`.

Salidas: los JSON con sufijo `_corregido`, la tabla por pantalla y
`controles_corregidos.csv`.

Uso: python recalcular_rmsd.py
"""
import csv
import glob
import json
import os
import re

import redock_trp32 as R
import protocolo_alternativo as P
import redock_flexible as F

BASE = R.BASE
SISTEMAS = R.SISTEMAS
SEMILLAS = [42, 2026, 777]
CABECERAS = ["control", "condicion", "entrada", "ligando", "puntuacion",
             "caja", "exhaustividad", "rmsd_por_semilla", "rmsd_mejor",
             "rmsd_mediana", "semillas_que_pasan", "criterio"]


def resumen(rmsds):
    vals = [v for v in rmsds if v is not None]
    if not vals:
        return None, None, 0
    mediana = sorted(vals)[len(vals) // 2]
    return min(vals), mediana, sum(1 for v in vals if v <= 2.0)


def fila(control, condicion, entrada, codigo, puntuacion, caja, exh, rmsds):
    mejor, mediana, pasa = resumen(rmsds)
    return {
        "control": control, "condicion": condicion, "entrada": entrada,
        "ligando": codigo, "puntuacion": puntuacion, "caja": caja,
        "exhaustividad": exh,
        "rmsd_por_semilla": " / ".join(
            "%.2f" % v if v is not None else "n/d" for v in rmsds),
        "rmsd_mejor": round(mejor, 2) if mejor is not None else None,
        "rmsd_mediana": round(mediana, 2) if mediana is not None else None,
        "semillas_que_pasan": pasa,
        "criterio": "PASA" if (mediana is not None and mediana <= 2.0) else "FALLA",
    }


def control1(filas):
    """Receptor rigido, ligando flexible, ligando de partida el ideal de RCSB."""
    for caja, exh in ((24, 8), (18, 8), (24, 32)):
        resultado = []
        for entrada, (codigo, _, nombre) in SISTEMAS.items():
            cry = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
            rmsds, detalle = [], []
            for s in SEMILLAS:
                f = os.path.join(BASE, "out", "%s_caja%d_e%d_s%d.pdbqt"
                                 % (entrada, caja, exh, s))
                if not os.path.exists(f):
                    rmsds.append(None)
                    continue
                v = R.rmsd_en_sitio(cry, R.molecula_dockeada(f))
                rmsds.append(v)
                detalle.append({"semilla": s, "rmsd": round(v, 2) if v else None})
            f = fila("control 1", "receptor rigido, ligando flexible", entrada,
                     codigo, "vina", caja, exh, rmsds)
            filas.append(f)
            f["por_semilla"] = detalle
            resultado.append(f)
        with open(os.path.join(BASE, "resultado_redocking_caja%d_e%d_corregido.json"
                               % (caja, exh)), "w", encoding="utf-8") as fh:
            json.dump(resultado, fh, indent=2, ensure_ascii=False)


def control2(filas):
    """Receptor rigido, ligando flexible, ligando de partida la pose del cristal."""
    for puntuacion in ("vina", "vinardo"):
        resultado = []
        for entrada, (codigo, _, nombre) in SISTEMAS.items():
            cry = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
            rmsds, detalle = [], []
            for s in SEMILLAS:
                f = os.path.join(BASE, "out", "%s_%s_caja24_e8_s%d.pdbqt"
                                 % (entrada, puntuacion, s))
                if not os.path.exists(f):
                    rmsds.append(None)
                    continue
                v = R.rmsd_en_sitio(cry, R.molecula_dockeada(f))
                rmsds.append(v)
                detalle.append({"semilla": s, "rmsd": round(v, 2) if v else None,
                                "afinidad": P.afinidad_de_salida(f)})
            f = fila("control 2", "receptor rigido, ligando flexible", entrada,
                     codigo, puntuacion, 24, 8, rmsds)
            filas.append(f)
            f["por_semilla"] = detalle
            resultado.append(f)
        with open(os.path.join(BASE, "resultado_protocolo_alternativo_completo_%s_corregido.json"
                               % puntuacion), "w", encoding="utf-8") as fh:
            json.dump(resultado, fh, indent=2, ensure_ascii=False)


def con_remark_smiles(salida_rigida, cabecera_pdbqt, destino):
    """Copia las lineas REMARK SMILES del ligando de entrada a la salida rigida.

    El ligando rigido se escribio a mano (todos los atomos en el ROOT) y sin esa
    cabecera meeko no puede reconstruir la molecula, asi que el RMSD se midio por
    orden de atomo. Son la misma molecula, asi que la cabecera vale; con ella se
    puede aplicar la correccion de simetria (los dos metiles del isoproterenol).
    """
    cab = [l for l in open(cabecera_pdbqt, encoding="utf-8", errors="ignore")
           if l.startswith("REMARK")]
    cuerpo = [l for l in open(salida_rigida, encoding="utf-8", errors="ignore")
              if not l.startswith("REMARK")]
    with open(destino, "w", encoding="utf-8") as f:
        f.writelines(cab)
        f.writelines(cuerpo)
    return destino


def control2_rigido(filas):
    """Receptor rigido y ligando RIGIDO en la conformacion del cristal."""
    resultado = []
    for etiqueta in ("vina-rigido", "vinardo-rigido"):
        for entrada, (codigo, _, nombre) in SISTEMAS.items():
            cry = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
            entrada_pdbqt = os.path.join(BASE, "ligands", "%s_cristal_%s.pdbqt"
                                         % (entrada, etiqueta.split("-")[0]))
            rmsds, detalle = [], []
            for s in SEMILLAS:
                bruto = os.path.join(BASE, "out", "%s_%s_caja24_e8_s%d.pdbqt"
                                     % (entrada, etiqueta, s))
                if not os.path.exists(bruto):
                    rmsds.append(None)
                    continue
                con = os.path.join(BASE, "out", "_con_remark_%s_%s_s%d.pdbqt"
                                   % (entrada, etiqueta, s))
                con_remark_smiles(bruto, entrada_pdbqt, con)
                v = None
                try:
                    v = R.rmsd_en_sitio(cry, R.molecula_dockeada(con))
                except Exception:
                    v = F.rmsd_posicional(entrada_pdbqt, bruto)
                rmsds.append(v)
                detalle.append({"semilla": s, "rmsd": round(v, 2) if v else None,
                                "afinidad": P.afinidad_de_salida(bruto)})
            f = fila("control 2", "receptor rigido, ligando RIGIDO", entrada,
                     codigo, etiqueta, 24, 8, rmsds)
            filas.append(f)
            f["por_semilla"] = detalle
            resultado.append(f)
    with open(os.path.join(BASE, "resultado_protocolo_alternativo_rigido_corregido.json"),
              "w", encoding="utf-8") as fh:
        json.dump(resultado, fh, indent=2, ensure_ascii=False)


def control3(filas):
    """Receptor flexible en los residuos del bolsillo."""
    for archivo in sorted(glob.glob(os.path.join(
            BASE, "resultado_redocking_flexible_*.json"))):
        if "corregido" in archivo:
            continue
        datos = json.load(open(archivo, encoding="utf-8"))
        for r in datos:
            entrada = r["entrada"]
            cry = os.path.join(BASE, "ligands", "%s_cristal.pdb" % entrada)
            entrada_pdbqt = os.path.join(BASE, "ligands",
                                         "%s_cristal_flexp.pdbqt" % entrada)
            exhaustividad = r.get("exhaustividad", 8)
            rmsds, detalle = [], []
            for s in SEMILLAS:
                f = os.path.join(BASE, "out", "%s_flex%s_%s_caja%d_e%d_s%d.pdbqt"
                                 % (entrada, r["residuos_flexibles"].replace(",", "-"),
                                    r["puntuacion"], r["caja"], exhaustividad, s))
                if not os.path.exists(f):
                    rmsds.append(None)
                    continue
                v = R.rmsd_en_sitio(cry, R.molecula_dockeada(f))
                rmsds.append(v)
                detalle.append({"semilla": s, "rmsd": round(v, 2) if v else None,
                                "afinidad": P.afinidad_de_salida(f),
                                "dist_centroide": round(F.dist_centroides(entrada_pdbqt, f), 2),
                                "rmsd_alineado": round(F.rmsd_alineado(cry, f), 2)})
            f = fila("control 3", "receptor flexible (A:%s), ligando flexible"
                     % r["residuos_flexibles"], entrada, r["ligando"],
                     r["puntuacion"], r["caja"], exhaustividad, rmsds)
            filas.append(f)
            f["por_semilla"] = detalle
        with open(archivo.replace(".json", "_corregido.json"), "w",
                  encoding="utf-8") as fh:
            json.dump(datos, fh, indent=2, ensure_ascii=False)


def main():
    filas = []
    control1(filas)
    control2(filas)
    control2_rigido(filas)
    control3(filas)

    with open(os.path.join(BASE, "controles_corregidos.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=CABECERAS, extrasaction="ignore")
        w.writeheader()
        for r in filas:
            w.writerow(r)

    print("%-10s %-46s %-6s %-14s %-4s %-4s %-22s %-7s %-7s %s"
          % ("control", "condicion", "entrada", "puntuacion", "caja", "exh",
             "rmsd por semilla", "mejor", "mediana", "criterio"))
    for r in filas:
        print("%-10s %-46s %-6s %-14s %-4s %-4s %-22s %-7s %-7s %s"
              % (r["control"], r["condicion"], r["entrada"], r["puntuacion"],
                 r["caja"], r["exhaustividad"], r["rmsd_por_semilla"],
                 r["rmsd_mejor"], r["rmsd_mediana"], r["criterio"]))

    print("\nCaracterizacion del fallo con receptor flexible:")
    for r in filas:
        if r["control"] != "control 3":
            continue
        for s in r.get("por_semilla", []) or []:
            print("   %-5s semilla %-5s RMSD %-7s alineado %-7s centroide %-7s afinidad %s"
                  % (r["entrada"], s["semilla"], s["rmsd"], s["rmsd_alineado"],
                     s["dist_centroide"], s["afinidad"]))
    print("\nescrito: controles_corregidos.csv y los *_corregido.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
