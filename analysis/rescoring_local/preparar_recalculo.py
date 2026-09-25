# -*- coding: utf-8 -*-
"""Arma la lista para recalcular con el PROTOCOLO VALIDADO (receptor fijo).

Por que: los 2.731 resultados de la lista corta se calcularon con los
receptores viejos de gpu_dock/ (SOD1 con ~18 copias del monomero), no con los
congelados del snapshot. Sus energias (-80 a -137 kcal/mol) no son comparables
con las de la validacion (-8 a -16) y no se sabe si ese score discrimina.

Que selecciona:
  * por diana, la union del top-60 por RESIDUAL y el top-60 por EFICIENCIA
    de ligando, con tope de 100 por diana (los mejores por cualquiera de las
    dos metricas normalizadas);
  * los CONTROLES conocidos con pose (para verificar dentro de la propia
    corrida, con el mismo receptor y la misma plataforma).

Salida: lista_recalculo.csv  (target, ligand, smiles, vina_affinity, tipo)
Uso:    python preparar_recalculo.py [tope_por_diana]
"""
import csv
import os
import sys

LOCAL = r"C:\Users\Fredy\masive-als\analysis\rescoring_local"
RANKING = os.path.join(LOCAL, "rescoring_ranking_normalizado.csv")
CORTAS = os.path.join(LOCAL, "rescoring_lista_corta.csv")
VAL = os.path.join(LOCAL, "validacion_controles.csv")
CONT_TDP43 = os.path.join(LOCAL, "lista_controles_tdp43.csv")
SALIDA = os.path.join(LOCAL, "lista_recalculo.csv")

TOPE = 100          # por diana
N_POR_METRICA = 60  # cuantos mira de cada metrica antes de recortar


def leer_smiles():
    """SMILES por (diana, ligando), de la lista corta y de los controles."""
    m = {}
    for path in (CORTAS, VAL, CONT_TDP43):
        if not os.path.exists(path):
            continue
        for r in csv.DictReader(open(path, encoding="utf-8")):
            t, l, s = r.get("target", ""), r.get("ligand", ""), r.get("smiles", "")
            if t and l and s:
                m.setdefault((t, l), s)
    return m


def main():
    tope = int(sys.argv[1]) if len(sys.argv) > 1 else TOPE
    smi = leer_smiles()

    filas = list(csv.DictReader(open(RANKING, encoding="utf-8")))
    por_diana = {}
    for r in filas:
        por_diana.setdefault(r["target"], []).append(r)

    elegidos = {}
    for t, sub in por_diana.items():
        por_res = sorted(sub, key=lambda r: int(r["puesto_residual"]))[:N_POR_METRICA]
        por_eff = sorted(sub, key=lambda r: int(r["puesto_eficiencia"]))[:N_POR_METRICA]
        vistos = {}
        for r in por_res:
            vistos.setdefault(r["ligand"], set()).add("residual")
        for r in por_eff:
            vistos.setdefault(r["ligand"], set()).add("eficiencia")
        orden = sorted(vistos.items(),
                       key=lambda kv: min(int(next(x["puesto_residual"] for x in sub
                                                   if x["ligand"] == kv[0])),
                                          int(next(x["puesto_eficiencia"] for x in sub
                                                   if x["ligand"] == kv[0]))))
        for lig, metricas in orden[:tope]:
            fila = next(x for x in sub if x["ligand"] == lig)
            tipo = ("candidato_ambos" if len(metricas) == 2
                    else "candidato_" + list(metricas)[0])
            elegidos[(t, lig)] = {
                "target": t, "ligand": lig,
                "smiles": smi.get((t, lig), ""),
                "vina_affinity": fila.get("vina_affinity", ""),
                "tipo": tipo,
            }

    # Controles conocidos con pose: verificacion dentro de la propia corrida
    n_ctrl = 0
    for path in (VAL, CONT_TDP43):
        if not os.path.exists(path):
            continue
        for r in csv.DictReader(open(path, encoding="utf-8")):
            if r.get("tipo") != "control":
                continue
            k = (r["target"], r["ligand"])
            if k in elegidos:
                elegidos[k]["tipo"] += "+control"
                continue
            elegidos[k] = {"target": r["target"], "ligand": r["ligand"],
                           "smiles": r.get("smiles", ""),
                           "vina_affinity": r.get("afinidad", ""),
                           "tipo": "control"}
            n_ctrl += 1

    fuera = [k for k, v in elegidos.items() if not v["smiles"]]
    for k in fuera:
        del elegidos[k]

    with open(SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["target", "ligand", "smiles",
                                          "vina_affinity", "tipo"])
        w.writeheader()
        w.writerows(sorted(elegidos.values(), key=lambda v: (v["target"], v["tipo"], v["ligand"])))

    print("lista escrita: %s" % SALIDA)
    print("total filas: %d  (controles: %d | sin smiles descartados: %d)"
          % (len(elegidos), n_ctrl, len(fuera)))
    por_t = {}
    for v in elegidos.values():
        por_t.setdefault(v["target"], []).append(v["tipo"])
    for t in sorted(por_t):
        c = {}
        for x in por_t[t]:
            c[x] = c.get(x, 0) + 1
        print("  %-9s %d  %s" % (t, len(por_t[t]), c))
    return 0


if __name__ == "__main__":
    sys.exit(main())
