#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""diagnostico_score_only.py — ¿por que la afinidad de `--score_only` NO es la del docking?

QUE ESTA PASANDO
----------------
`comprobar_score_only.py` medio una diferencia de 2 kcal de mediana entre lo que Vina
escribio al acoplar y lo que devuelve `--score_only` sobre la misma pose. Eso es
demasiado para una simple relajacion local, asi que antes de usar estas columnas hay
que saber de donde sale. La unica explicacion que se prueba aqui es la CAJA:

Vina solo cuenta los atomos del receptor que caen dentro de la caja de busqueda, y al
acoplar la caja es la del docking (fija, en el centro del bolsillo). Aqui la caja va
centrada en la POSE, porque `--score_only` exige que el ligando este dentro. Con la
caja del docking, la afinidad deberia volver a salir igual; si sale, la diferencia era
la caja y no la funcion.

Uso:
    python diagnostico_score_only.py --target TDP43
"""
import argparse
import os
import subprocess
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
VINA = os.path.join(RAIZ, "tools", "vina.exe")
sys.path.insert(0, BASE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", default="TDP43", choices=["TDP43", "SOD1"])
    ap.add_argument("--n", type=int, default=8)
    args = ap.parse_args()

    import validar_sod1_limpia as V
    import validar_diana_limpia as D
    import repunuar_vinardo as R
    cfg = D.DIANAS[args.target]
    centro, tamano = cfg["centro"], cfg["tamano"]

    poses2 = cfg.get("poses2") or os.path.join(os.path.dirname(V.LIGS_FONDO2), "out")
    d = poses2
    nombres = sorted(p for p in os.listdir(d) if p.endswith("_out.pdbqt"))[:args.n]

    print("DIANA %s | %d poses de %s" % (args.target, len(nombres), os.path.basename(d)))
    print("%-18s %10s %10s %10s %10s"
          % ("ligando", "docking", "caja docking", "caja pose", "d - pose"))
    for p in nombres:
        ruta = os.path.join(d, p)
        aff = V.afinidad(ruta)
        lineas, coords = R.un_modelo(ruta)
        tmp = os.path.join(BASE, "_vinardo", "diag.pdbqt")
        os.makedirs(os.path.dirname(tmp), exist_ok=True)
        R.a_ligando_unico(lineas, tmp)
        out = {}
        for etiqueta, cx, cy, cz, lado in (
                ("docking", centro[0], centro[1], centro[2], tamano),
                ("pose", coords.mean(axis=0)[0], coords.mean(axis=0)[1],
                 coords.mean(axis=0)[2],
                 max(float(tamano), float((coords.max(axis=0)
                                          - coords.min(axis=0)).max()) + 4.0))):
            r = subprocess.run(
                [VINA, "--receptor", cfg["receptor"], "--ligand", tmp,
                 "--score_only", "--scoring", "vina",
                 "--center_x", "%.3f" % cx, "--center_y", "%.3f" % cy,
                 "--center_z", "%.3f" % cz,
                 "--size_x", "%.1f" % lado, "--size_y", "%.1f" % lado,
                 "--size_z", "%.1f" % lado],
                capture_output=True, text=True, timeout=600)
            t = r.stdout + r.stderr
            import re
            m = re.search(r"Estimated Free Energy of Binding\s*:\s*(-?[0-9.]+)", t)
            out[etiqueta] = float(m.group(1)) if m else None
        if os.path.exists(tmp):
            os.remove(tmp)
        print("%-18s %10s %10s %10s %10s"
              % (p.replace("_out.pdbqt", ""), aff,
                 out["docking"], out["pose"],
                 None if None in (out["docking"], out["pose"])
                 else "%.3f" % (out["docking"] - out["pose"])))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
