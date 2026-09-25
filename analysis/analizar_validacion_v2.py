# -*- coding: utf-8 -*-
"""Analiza validaciones v2 (TDP43/FUS) y propone corte calibrado por señuelos."""
import csv
import sys
import io
import numpy as np

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


def leer(path):
    act, dec = [], []
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            try:
                v = float(r["affinity"])
            except Exception:
                continue
            if r["rol"] == "activo":
                act.append(v)
            else:
                dec.append(v)
    return np.array(act), np.array(dec)


def roc_auc(act, dec):
    if len(act) == 0 or len(dec) == 0:
        return None
    auc = 0.0
    for a in act:
        auc += float(np.sum(a <= dec))
    return auc / (len(act) * len(dec))


def ef(act, dec, pct=1.0):
    todo = [(s, 1) for s in act] + [(s, 0) for s in dec]
    todo.sort(key=lambda x: x[0])
    k = max(1, int(len(todo) * pct / 100.0))
    n_act = sum(1 for _, lab in todo[:k] if lab == 1)
    return (n_act / max(1, len(act))) / (pct / 100.0)


def corte_fp5(dec):
    """Score mas negativo que el 95% de los decoys (5% falsos positivos)."""
    if len(dec) == 0:
        return None
    return float(np.percentile(dec, 95))


def corte_fp1(dec):
    return float(np.percentile(dec, 99)) if len(dec) else None


def main():
    for name, path in [("SOD1 (referencia)", "validacion_SOD1.csv"),
                       ("TDP43 v2", "validacion_TDP43_v2.csv"),
                       ("FUS v2", "validacion_FUS_v2.csv")]:
        print("=" * 50)
        print(name, "—", path)
        try:
            act, dec = leer(path)
        except FileNotFoundError:
            print("  (aun no existe)")
            continue
        print("  activos: %d | decoys: %d" % (len(act), len(dec)))
        if len(act) == 0:
            print("  sin activos con score")
            continue
        print("  media activos: %.2f | media decoys: %.2f" % (act.mean(), dec.mean()))
        print("  ROC-AUC: %.3f" % roc_auc(act, dec))
        print("  EF1%%: %.2f | EF5%%: %.2f" % (ef(act, dec, 1.0), ef(act, dec, 5.0)))
        fp5 = corte_fp5(dec)
        fp1 = corte_fp1(dec)
        print("  CORTE CALIBRADO: 5%% FP -> %.2f | 1%% FP -> %.2f" % (fp5, fp1))
        nact_bajo_fp5 = int(np.sum(act <= fp5))
        print("  activos que pasan el corte 5%% FP: %d/%d" % (nact_bajo_fp5, len(act)))


if __name__ == "__main__":
    main()
