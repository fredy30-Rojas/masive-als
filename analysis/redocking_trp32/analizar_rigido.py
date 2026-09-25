#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Recalcula el RMSD de las corridas rigidas ya hechas (no vuelve a acoplar).

La primera version del analisis no pudo medir el RMSD del ligando rigido porque
la salida de Vina trae los 9 modos concatenados en el mismo fichero y hay que
leer solo el primero. Los acoplamientos ya estan hechos en `out/`, asi que aqui
solo se remiden (serian horas de CPU repetirlos).

Uso: python analizar_rigido.py
"""
import json
import os

import protocolo_alternativo as P

BASE = os.path.dirname(os.path.abspath(__file__))
CAJA, EXH = 24, 8
SEMILLAS = [42, 2026, 777]
JSON = os.path.join(BASE, "resultado_protocolo_alternativo_rigido_vina_vinardo.json")


def main():
    rig = json.load(open(JSON, encoding="utf-8"))
    # E de la pose del cristal quieta y minimizada ya estan en el JSON
    completos = {}
    for f in ("resultado_protocolo_alternativo_completo_vina.json",
              "resultado_protocolo_alternativo_completo_vinardo.json"):
        ruta = os.path.join(BASE, f)
        if os.path.exists(ruta):
            for r in json.load(open(ruta, encoding="utf-8")):
                completos[(r["entrada"], r["puntuacion"])] = r

    for r in rig:
        if "por_semilla" not in r:
            continue
        entrada, etiqueta = r["entrada"], r["puntuacion"]
        inp = os.path.join(BASE, "ligands", "%s_cristal_rigido.pdbqt" % entrada)
        for s in r["por_semilla"]:
            out = os.path.join(BASE, "out", "%s_%s_caja%d_e%d_s%d.pdbqt"
                               % (entrada, etiqueta, CAJA, EXH, s["semilla"]))
            s["rmsd"] = None
            if os.path.exists(out):
                v = P.rmsd_posicional(inp, out)
                s["rmsd"] = round(v, 2) if v is not None else None
        rmsds = [s["rmsd"] for s in r["por_semilla"] if s["rmsd"] is not None]
        r["rmsd_mejor"] = min(rmsds) if rmsds else None
        r["rmsd_mediana"] = sorted(rmsds)[len(rmsds) // 2] if rmsds else None
        r["semillas_que_pasan"] = sum(1 for v in rmsds if v <= 2.0)
        r["criterio"] = ("PASA" if (r["rmsd_mediana"] is not None
                                   and r["rmsd_mediana"] <= 2.0) else "FALLA")
        print("%-5s %-14s cristal %-7s minim %-7s deriva %-5s | mejor %-6s mediana %-6s %s"
              % (entrada, etiqueta, r["E_cristal_quieta"],
                 r["E_cristal_minimizado"], r["deriva_minimizacion"],
                 r["rmsd_mejor"], r["rmsd_mediana"], r["criterio"]))
        for s in r["por_semilla"]:
            print("        semilla %-5d afinidad %-7s RMSD %s"
                  % (s["semilla"], s["afinidad"], s["rmsd"]))

    json.dump(rig, open(JSON, "w", encoding="utf-8"), indent=2, ensure_ascii=False)
    print("actualizado: %s" % JSON)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
