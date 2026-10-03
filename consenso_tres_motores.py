#!/usr/bin/env python3
"""
Cruza los tres motores en un solo ranking de consenso.

Los tres miden cosas distintas de la misma molecula:
  - Vina-GPU : afinidad empirica del docking (kcal/mol, mas negativo = mejor).
  - GNINA    : re-puntaje con red neuronal (CNN affinity, mas alto = mejor).
  - MM-GBSA  : energia libre de la pose (kcal/mol, mas negativo = mejor).

Una molecula que sale buena en los tres es mucho mas creible que una que solo
brilla en uno. Este script junta las tres fuentes por (diana, ligando), calcula
un puesto normalizado en cada motor y un consenso (media de los tres puestos,
mas bajo = mejor), y marca los candidatos que aguantan en pie en los tres.

Entradas:
  resultados/ranking_final_masive_als.csv         (MM-GBSA + metadata)
  gpu_dock/resultados_libreria/resultados_libreria_total.csv  (Vina)
  gpu_dock/gnina_top<N>/rescore_<target>.csv      (GNINA, puede estar a medias)

Salida:
  resultados/consenso_tres_motores.csv

Uso:  python3 consenso_tres_motores.py [--gnina-top 5000]
"""
import argparse
import csv
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
RANKING = os.path.join(BASE, "resultados", "ranking_final_masive_als.csv")
VINA = os.path.join(BASE, "gpu_dock", "resultados_libreria",
                    "resultados_libreria_total.csv")
SALIDA = os.path.join(BASE, "resultados", "consenso_tres_motores.csv")

# La diana TDP43_v2 se trata como TDP43 (es la misma proteina, otra conformacion)
ALIAS_TARGET = {"TDP43_v2": "TDP43"}


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def leer_mmgbsa():
    """{ (target, ligand): {mmgbsa, sospechoso, fuente} }"""
    datos = {}
    if not os.path.isfile(RANKING):
        return datos
    with open(RANKING, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            t = ALIAS_TARGET.get(row["target"], row["target"])
            datos[(t, row["ligand"])] = {
                "mmgbsa": num(row.get("mmgbsa_dG")),
                "sospechoso": row.get("sospechoso", ""),
                "fuente": row.get("fuente", ""),
            }
    return datos


def leer_vina():
    """{ (target, ligand): afinidad } con la mejor (mas negativa) por pareja."""
    datos = {}
    if not os.path.isfile(VINA):
        return datos
    with open(VINA, newline="", encoding="utf-8-sig") as fh:
        for row in csv.DictReader(fh):
            t = ALIAS_TARGET.get(row["target"], row["target"])
            af = num(row.get("affinity"))
            if af is None:
                continue
            clave = (t, row["ligand"])
            if clave not in datos or af < datos[clave]:
                datos[clave] = af
    return datos


def leer_gnina(gnina_top):
    """{ (target, ligand): {cnn_score, cnn_affinity, vina} }"""
    datos = {}
    dir_gnina = os.path.join(BASE, "gpu_dock", f"gnina_top{gnina_top}")
    if not os.path.isdir(dir_gnina):
        return datos
    for nombre in os.listdir(dir_gnina):
        if not (nombre.startswith("rescore_") and nombre.endswith(".csv")):
            continue
        ruta = os.path.join(dir_gnina, nombre)
        with open(ruta, newline="", encoding="utf-8-sig") as fh:
            for row in csv.DictReader(fh):
                t = ALIAS_TARGET.get(row["target"], row["target"])
                datos[(t, row["ligand"])] = {
                    "cnn_score": num(row.get("cnn_score")),
                    "cnn_affinity": num(row.get("cnn_affinity")),
                    "error": row.get("error", ""),
                }
    return datos


def puesto(valores, mejor="menor"):
    """Devuelve {clave: puesto} (1 = mejor), ignorando los que no tienen valor."""
    con_valor = [(k, v) for k, v in valores.items() if v is not None]
    con_valor.sort(key=lambda kv: kv[1], reverse=(mejor == "mayor"))
    return {k: i + 1 for i, (k, _) in enumerate(con_valor)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gnina-top", type=int, default=5000)
    args = ap.parse_args()

    mmgbsa = leer_mmgbsa()
    vina = leer_vina()
    gnina = leer_gnina(args.gnina_top)

    if not mmgbsa:
        print("No hay ranking de MM-GBSA; nada que cruzar.")
        return 1

    # universo de candidatos: los del ranking de MM-GBSA
    claves = set(mmgbsa)
    # metricas por motor
    m_mmgbsa = {k: mmgbsa[k]["mmgbsa"] for k in claves}
    m_vina = {k: vina.get(k) for k in claves}
    m_cnn = {k: (gnina.get(k, {}) or {}).get("cnn_affinity") for k in claves}
    m_cnns = {k: (gnina.get(k, {}) or {}).get("cnn_score") for k in claves}

    p_mmgbsa = puesto(m_mmgbsa, "menor")
    p_vina = puesto(m_vina, "menor")
    p_cnn = puesto(m_cnn, "mayor")

    filas = []
    for k in claves:
        target, ligand = k
        puestos = [p_mmgbsa.get(k), p_vina.get(k), p_cnn.get(k)]
        presentes = [p for p in puestos if p is not None]
        n_motores = len(presentes)
        consenso = (sum(presentes) / n_motores) if presentes else None
        # "los tres": esta bien en los tres motores (tercio superior de cada uno)
        umbral = max(1, len(claves) // 3)
        en_los_tres = (
            p_mmgbsa.get(k, 10**9) <= umbral
            and p_vina.get(k, 10**9) <= umbral
            and p_cnn.get(k, 10**9) <= umbral
        )
        filas.append({
            "target": target,
            "ligand": ligand,
            "consenso": round(consenso, 1) if consenso is not None else "",
            "n_motores": n_motores,
            "los_tres": "SI" if en_los_tres else "",
            "puesto_consenso": "",
            "mmgbsa_dG": mmgbsa[k]["mmgbsa"] if mmgbsa[k]["mmgbsa"] is not None else "",
            "puesto_mmgbsa": p_mmgbsa.get(k, ""),
            "vina_affinity": vina.get(k, ""),
            "puesto_vina": p_vina.get(k, ""),
            "cnn_affinity": m_cnn[k] if m_cnn[k] is not None else "",
            "cnn_score": m_cnns[k] if m_cnns[k] is not None else "",
            "puesto_cnn": p_cnn.get(k, ""),
            "sospechoso": mmgbsa[k]["sospechoso"],
        })

    filas.sort(key=lambda f: (f["n_motores"] != 3, f["consenso"] if f["consenso"] != "" else 10**9))
    for i, f in enumerate(filas, 1):
        f["puesto_consenso"] = i

    campos = ["puesto_consenso", "target", "ligand", "consenso", "n_motores",
              "los_tres", "mmgbsa_dG", "puesto_mmgbsa", "vina_affinity",
              "puesto_vina", "cnn_affinity", "cnn_score", "puesto_cnn",
              "sospechoso"]
    with open(SALIDA, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=campos)
        w.writeheader()
        for f in filas:
            w.writerow({c: f[c] for c in campos})

    con_tres = [f for f in filas if f["n_motores"] == 3]
    los_tres = [f for f in filas if f["los_tres"] == "SI"]
    print(f"candidatos: {len(filas)}")
    print(f"con los tres motores medidos: {len(con_tres)}")
    print(f"buenos en los tres (tercio superior): {len(los_tres)}")
    print(f"escrito: {SALIDA}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
