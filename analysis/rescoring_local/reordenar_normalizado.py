# -*- coding: utf-8 -*-
"""Reordena los resultados del rescoring quitando el sesgo de tamano.

CONTEXTO
La validacion del 11-12 sep 2026 (VALIDACION_RESCORING_2026-09-11.md) mostro
que el dG bruto del MM-GBSA esta dominado por el tamano molecular, aprox.
-1 kcal/mol por atomo pesado. Ordenar por dG bruto es, en la practica, ordenar
por tamano, y entierra a los activos pequenos conocidos: AUC 0.253 en SOD1 y
0.140 en TDP43_v2, los dos POR DEBAJO del azar (p < 0.002).

Aqui se reordena la lista corta con las metricas normalizadas por tamano que
la validacion senalo:

    residual = dG - (a + b * atomos_pesados)   con a,b ajustados en la diana
    LE       = dG / atomos_pesados             (eficiencia de ligando)

Residual NEGATIVO = mejor de lo que su tamano predice. LE negativa mas grande
en valor absoluto = mas eficiente por atomo.

ESTADO DE VALIDACION POR DIANA (no inventar: sale del informe del 11 sep):
    SOD1      -> senal moderada (AUC residual 0.754)  -> el reordenamiento aporta
    TDP43_v2  -> NEGATIVO FIRME (AUC residual 0.073)  -> NO usar para ordenar
    FUS       -> sin controles con pose, no evaluable -> NO usar para ordenar

El script calcula las tres, pero el informe debe decir cual es usable.

RUIDO: la corrida se hizo en CUDA float32, que lleva +-1,5 kcal/mol de
dispersion y ademas SESGA el valor (mismo compuesto: -8 en float32 vs -14,5 en
float64). Por eso se marca aparte lo que supera ese margen.

Uso:  python reordenar_normalizado.py
"""
import csv
import os
import sys

import numpy as np

BASE = r"C:\Users\Fredy\masive-als"
LOCAL = os.path.join(BASE, "analysis", "rescoring_local")
POSES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
ENTRADA = os.path.join(LOCAL, "rescoring_lista_corta.csv")
SALIDA = os.path.join(LOCAL, "rescoring_ranking_normalizado.csv")
INFORME = os.path.join(LOCAL, "ranking_normalizado_resumen.txt")

RUIDO = 1.5  # kcal/mol, dispersion medida en float32 sobre CUDA

# Lo que la validacion dejo claro sobre cada diana
VALIDACION = {
    "SOD1": ("senal moderada (AUC residual 0.754) - reordenamiento utilizable",
             True),
    "TDP43_v2": ("NEGATIVO FIRME (AUC residual 0.073) - NO usar para ordenar",
                 False),
    "FUS": ("sin controles con pose - NO evaluable, NO usar para ordenar", False),
}


def atomos_pesados(pose):
    """Cuenta atomos pesados de la primera pose (mismo criterio que analizar_validacion.py)."""
    n = 0
    try:
        with open(pose, encoding="utf-8", errors="ignore") as f:
            for l in f:
                if l.startswith("ENDMDL"):
                    break
                if l.startswith(("ATOM", "HETATM")) and l[76:79].strip().upper() != "H":
                    n += 1
    except OSError:
        return 0
    return n


def ajuste_robusto(ha, dg):
    """Recta dG ~ atomos pesados, quitando atipicos a 3 MAD y reajustando.

    Devuelve (b, a, n_usados). El ajuste OLS simple se deja llevat por las colas
    de dG (poses malas, compuestos raros), asi que se recorta una vez.
    """
    ha = np.asarray(ha, float)
    dg = np.asarray(dg, float)
    b, a = np.polyfit(ha, dg, 1)
    res = dg - (a + b * ha)
    mad = np.median(np.abs(res - np.median(res)))
    if mad > 0:
        keep = np.abs(res - np.median(res)) <= 3 * 1.4826 * mad
        if keep.sum() >= 30:
            b, a = np.polyfit(ha[keep], dg[keep], 1)
            return b, a, int(keep.sum())
    return b, a, int(len(ha))


def spearman(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    rx = np.argsort(np.argsort(x)).astype(float)
    ry = np.argsort(np.argsort(y)).astype(float)
    rx -= rx.mean()
    ry -= ry.mean()
    d = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / d) if d else 0.0


def main():
    if not os.path.exists(ENTRADA):
        print("no existe %s" % ENTRADA)
        return 1

    filas = []
    sin_pose = 0
    for r in csv.DictReader(open(ENTRADA, encoding="utf-8")):
        if not r.get("mmgbsa_dG"):
            continue
        pose = os.path.join(POSES, "results_" + r["target"], r["ligand"] + "_out.pdbqt")
        ha = atomos_pesados(pose)
        if ha == 0:
            sin_pose += 1
            continue
        filas.append({
            "target": r["target"],
            "ligand": r["ligand"],
            "vina": r["vina_affinity"],
            "dG": float(r["mmgbsa_dG"]),
            "ha": ha,
        })

    lineas = []
    def di(msg=""):
        print(msg)
        lineas.append(msg)

    di("REORDENAMIENTO POR TAMANO - MASIVE-ALS")
    di("Fecha: 2026-09-13   |   filas con dG y pose: %d   (sin pose: %d)"
       % (len(filas), sin_pose))
    di("Metrica: residual = dG - (a + b*atomos_pesados), ajustado por diana.")
    di("Residual negativo = mejor de lo que su tamano predice.")
    di("")

    salida = []
    for t in sorted({f["target"] for f in filas}):
        sub = [f for f in filas if f["target"] == t]
        ha = np.array([f["ha"] for f in sub], float)
        dg = np.array([f["dG"] for f in sub], float)
        b, a, n_aj = ajuste_robusto(ha, dg)
        res = dg - (a + b * ha)
        le = dg / ha

        # puestos: menor dG / menor residual / LE mas negativa = mejor
        orden_raw = np.argsort(dg)
        orden_res = np.argsort(res)
        orden_le = np.argsort(le)
        puesto_raw = np.empty(len(sub), int)
        puesto_res = np.empty(len(sub), int)
        puesto_le = np.empty(len(sub), int)
        puesto_raw[orden_raw] = np.arange(1, len(sub) + 1)
        puesto_res[orden_res] = np.arange(1, len(sub) + 1)
        puesto_le[orden_le] = np.arange(1, len(sub) + 1)

        for i, f in enumerate(sub):
            salida.append({
                "target": t,
                "ligand": f["ligand"],
                "vina_affinity": f["vina"],
                "mmgbsa_dG": "%.2f" % f["dG"],
                "atomos_pesados": f["ha"],
                "eficiencia_ligando": "%.4f" % le[i],
                "residual": "%.2f" % res[i],
                "puesto_dG_bruto": int(puesto_raw[i]),
                "puesto_residual": int(puesto_res[i]),
                "puesto_eficiencia": int(puesto_le[i]),
                "sube_puestos": int(puesto_raw[i] - puesto_res[i]),
            })

        rho_raw_res = spearman(dg, res)
        rho_raw_le = spearman(dg, le)
        rho_tam = spearman(ha, dg)
        n = len(sub)
        k100 = min(100, n)
        top_raw = set(orden_raw[:k100])
        top_res = set(orden_res[:k100])
        solapan = len(top_raw & top_res)
        grandes = int((res < -RUIDO).sum())   # mejor que su tamano por mas del ruido

        estado, usable = VALIDACION.get(t, ("sin dato de validacion", False))
        di("=" * 72)
        di("%s   n=%d" % (t, n))
        di("  VALIDACION: %s" % estado)
        di("  recta ajustada: dG = %.2f + %.3f * atomos_pesados   (n usado %d)"
           % (a, b, n_aj))
        di("  atomos pesados: media %.1f  (min %d, max %d)"
           % (ha.mean(), ha.min(), ha.max()))
        di("  Spearman dG vs tamano = %+.3f  <- cuanto premia el tamano" % rho_tam)
        di("  Spearman dG bruto vs residual = %+.3f" % rho_raw_res)
        di("  Candidatos mejores que su tamano por mas del ruido (%.1f): %d"
           % (RUIDO, grandes))
        di("  Solapamiento del top-100 bruto con el top-100 por residual: %d de %d"
           % (solapan, k100))

        di("")
        di("  --- TOP 15 por RESIDUAL (mejor de lo que su tamano predice) ---")
        di("  %-22s %8s %5s %10s %10s %8s" %
           ("ligando", "dG", "atomos", "residual", "puesto_dG", "sube"))
        for idx in orden_res[:15]:
            di("  %-22s %8.1f %5d %10.2f %10d %8d"
               % (sub[idx]["ligand"], sub[idx]["dG"], ha[idx], res[idx],
                  puesto_raw[idx], puesto_raw[idx] - puesto_res[idx]))

        di("")
        di("  --- TOP 5 bruto (para comparar: los que el dG bruto ponia arriba) ---")
        di("  %-22s %8s %5s %10s %10s" %
           ("ligando", "dG", "atomos", "residual", "puesto_res"))
        for idx in orden_raw[:5]:
            di("  %-22s %8.1f %5d %10.2f %10d"
               % (sub[idx]["ligand"], sub[idx]["dG"], ha[idx], res[idx],
                  puesto_res[idx]))
        di("")

    with open(SALIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(salida[0].keys()))
        w.writeheader()
        w.writerows(salida)

    di("CSV: %s  (%d filas)" % (SALIDA, len(salida)))
    di("")
    di("RECORDATORIO CIENTIFICO: solo en SOD1 el residual tiene senal medida.")
    di("En TDP43_v2 y FUS el rescoring esta descartado; su reordenamiento NO")
    di("debe usarse para decidir candidatos. Y todo esto son ±1,5 kcal/mol de ruido.")
    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
