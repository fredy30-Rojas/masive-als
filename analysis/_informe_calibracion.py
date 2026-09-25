# -*- coding: utf-8 -*-
"""Informe cientifico: calibracion del corte con controles positivos + CNS."""
import os
import csv
import json

BASE = r"C:\Users\Fredy\masive-als"
ANALYSIS = os.path.join(BASE, "analysis")


def main():
    with open(os.path.join(ANALYSIS, "ranking_consolidado_resumen.json"), encoding="utf-8") as f:
        resumen = json.load(f)

    cortes = resumen["top5_por_proteina"]

    # controles
    controles = {}
    with open(os.path.join(ANALYSIS, "controles_calibracion.csv"), encoding="utf-8") as f:
        for r in csv.DictReader(f):
            controles.setdefault(r["target"], []).append(r)

    lineas = []
    lineas.append("# CALIBRACION DEL CORTE CON CONTROLES POSITIVOS — MASIVE-ALS")
    lineas.append("")
    lineas.append("Generado: 2026-09-09 desde resultados en disco del GPU (RTX 4080).")
    lineas.append("")
    lineas.append("## 1. Estado del cribado (outputs validos por proteina)")
    lineas.append("")
    lineas.append("| Proteina | Outputs | Top-5%% (n) | Corte Vina (kcal/mol) | Mejor score |")
    lineas.append("|---|---|---|---|---|")
    for t in ["TDP43_v2", "SOD1", "FUS"]:
        c = cortes[t]
        lineas.append("| %s | %d | %d | %.1f | %.1f |" % (t, c["total"], c["n_top5"], c["corte_afinidad"], c["mejor"]))
    lineas.append("")
    lineas.append("Excluidos documentados (ligandos imposibles, sin relanzar): **%d**." % resumen["excluidos_documentados"])
    lineas.append("")

    lineas.append("## 2. Controles positivos conocidos vs corte")
    lineas.append("")
    for t in ["TDP43", "SOD1", "FUS"]:
        key = "TDP43_v2" if t == "TDP43" else t
        corte = cortes.get(key, {}).get("corte_afinidad", None)
        ctr = controles.get(t, [])
        con_aff = [(float(c["afinidad_vina"]), c["control"]) for c in ctr if c["afinidad_vina"]]
        en_lib = sum(1 for c in ctr if c["en_libreria"] == "si")
        if not con_aff:
            lineas.append("**%s**: %d controles, %d en libreria; ninguno con afinidad medible." % (t, len(ctr), en_lib))
            lineas.append("")
            continue
        peor = max(a for a, _ in con_aff)
        mejor = min(a for a, _ in con_aff)
        pasan = sum(1 for a, _ in con_aff if corte is not None and a <= corte)
        lineas.append("**%s**: %d controles (%d en libreria cribada). Afinidades Vina: %.1f a %.1f kcal/mol."
                      % (t, len(ctr), en_lib, peor, mejor))
        if corte is not None:
            lineas.append("- Corte top-5%%: **%.1f** → pasarian el corte **%d de %d** controles."
                          % (corte, pasan, len(con_aff)))
        for a, nombre in sorted(con_aff, reverse=True):
            marca = " ✔" if corte is not None and a <= corte else ""
            lineas.append("  - %-24s %.1f%s" % (nombre, a, marca))
        lineas.append("")

    lineas.append("## 3. Lectura cientifica")
    lineas.append("")
    lineas.append("- **SOD1** es la unica con calibracion util: 18/20 controles en libreria,"
                  " afinidades **-4.7 a -6.0** kcal/mol. El corte top-5%% en SOD1 es **%.1f**,"
                  " por debajo (mas estricto) que el mejor control conocido. Ningun control"
                  " conocido pasa el corte: el docking por si solo NO discrimina a los activos"
                  " conocidos; el ranking necesita rescoring (MM-GBSA) + filtros de quimica"
                  " medicinal (CNS, PAINS) antes de elegir candidatos." % cortes["SOD1"]["corte_afinidad"])
    lineas.append("")
    lineas.append("- **TDP-43** (v2, bolsillo corregido): solo 2/9 controles en libreria"
                  " (berberine -5.7, ketoconazole -6.3); los demas son sondas publicadas"
                  " (PE859, sanguinarine...) ausentes de la libreria de farmacos. Corte top-5%%: %.1f."
                  % cortes["TDP43_v2"]["corte_afinidad"])
    lineas.append("")
    lineas.append("- **FUS**: 0/2 controles en libreria; calibracion pendiente de incluir"
                  " ligandos FUS conocidos en futuras tandas. Corte top-5%%: %.1f."
                  % cortes["FUS"]["corte_afinidad"])
    lineas.append("")
    lineas.append("## 4. Recomendacion de corte (conservadora)")
    lineas.append("")
    lineas.append("Con los datos actuales, el corte por percentil NO es comparable entre"
                  " proteinas (los scores Vina dependen del receptor). Se recomienda:")
    lineas.append("1. Usar el ranking solo como filtro de enriquecimiento, no como valor absoluto.")
    lineas.append("2. Aplicar rescoring MM-GBSA sobre el top-5%% (14.303 compuestos).")
    lineas.append("3. Filtrar por CNS MPO >= 4 (o BBB clasica) — ya calculado en el CSV.")
    lineas.append("4. Validar los candidatos finales con controles positivos por proteina"
                  " (ampliar la lista FUS y TDP-43 en libreria).")
    lineas.append("")
    lineas.append("## 5. Archivos")
    lineas.append("")
    lineas.append("- `ranking_afinidades.csv` — 286.079 afinidades Vina (las 3 proteinas)")
    lineas.append("- `ranking_top5_consolidado.csv` — top-5%% con SMILES, propiedades y filtros CNS")
    lineas.append("- `controles_calibracion.csv` — controles con afinidad real en libreria")
    lineas.append("- `_excluidos_clasificados.csv` — los 5.696 ligandos imposibles con motivo")
    lineas.append("")

    out = os.path.join(ANALYSIS, "CALIBRACION_CORTE_2026-09-09.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas))
    print(out)
    print("\n".join(lineas))


if __name__ == "__main__":
    main()