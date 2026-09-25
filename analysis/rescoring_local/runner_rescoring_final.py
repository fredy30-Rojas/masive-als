# -*- coding: utf-8 -*-
"""Rescoring MM-GBSA de PRODUCCIÓN — MASIVE-ALS (lista corta + controles).

Protocolo (ver RESCORING_WSL_CUDA_2026-09-11.md y AUDITORIA_RESCORING_2026-09-11.md):

  * Disolvente implícito OBC2 real (GBSAOBC2Force, mbondi2, SA=ACE).
  * Geometría del ligando = la de la pose acoplada (sin re-embedding).
  * Receptor CONGELADO por diana (receptores_fijos/<target>_fijo.pdb),
    preparado con criterio de caja de docking por preparar_receptores.py.
    Todos los compuestos se evalúan contra exactamente los mismos átomos.
  * Plataforma determinista y de un solo hilo (CPU Threads=1 o Reference).
    CUDA se descarta: minimiza en float32 y devuelve mínimos distintos con
    entrada idéntica (medido: hasta 9 kcal/mol de dispersión con la energía
    inicial idéntica a 3 decimales).
  * Cada candidato es un PROCESO independiente, así que se paraleliza sin
    compartir estado.

Escribe a rescoring_final.csv (no mezcla con los datos invalidados).
Reanudable: al relanzar solo procesa lo que no tenga mmgbsa_dG.

Uso (desde WSL):
  ~/rescoring_env/bin/python runner_rescoring_final.py --workers 16
  ~/rescoring_env/bin/python runner_rescoring_final.py --controles
  ~/rescoring_env/bin/python runner_rescoring_final.py --workers 16 --max 40
"""
import argparse
import csv
import hashlib
import json
import os
import subprocess
import sys
import time
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

BASE = "/mnt/c/Users/Fredy/masive-als"
if not os.path.isdir(BASE):
    BASE = r"C:\Users\Fredy\masive-als"
LOCAL = os.path.join(BASE, "analysis", "rescoring_local")
HIJO = os.path.join(LOCAL, "mmgbsa_openff_gb.py")
FIJO = os.path.join(LOCAL, "receptores_fijos")
POSES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
LISTA = os.path.join(BASE, "analysis", "lista_corta_candidatos.csv")
CONTROLES = os.path.join(BASE, "analysis", "controles_calibracion.csv")

PLATFORM = os.environ.get("MMGBSA_FINAL_PLATFORM", "CPU")
THREADS = os.environ.get("MMGBSA_CPU_THREADS", "1")
TIMEOUT = int(os.environ.get("MMGBSA_TIMEOUT", "5400"))  # 90 min por candidato

COLS = ["target", "ligand", "vina_affinity", "mmgbsa_dG", "e_complex",
        "e_receptor", "e_ligand", "n_rec", "plataforma", "receptor_md5",
        "segundos", "error"]


def log(path, msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    print(line, flush=True)
    with open(path, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def receptor_de(target):
    return os.path.join(FIJO, "%s_fijo.pdb" % target)


def cargar_lista(targets=None, max_n=None):
    filas = []
    with open(LISTA, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if targets and r["target"] not in targets:
                continue
            filas.append({"target": r["target"], "ligand": r["ligand"],
                          "affinity": r.get("afinidad"), "smiles": r["smiles"]})
    if max_n:
        filas = filas[:max_n]
    return filas


def cargar_controles(targets=None):
    """Controles de calibración que SÍ están en la librería (tienen pose)."""
    filas = []
    with open(CONTROLES, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if (r.get("en_libreria", "").strip().lower() != "si"):
                continue
            t = r["target"]
            if targets and t not in targets:
                continue
            lig = (r.get("ligand_libreria") or "").strip()
            if not lig:
                continue
            filas.append({"target": t, "ligand": lig,
                          "affinity": r.get("afinidad_vina"),
                          "smiles": r["smiles"], "control": r["control"]})
    return filas


def cargar_hechos(out):
    hechos = {}
    if os.path.exists(out):
        with open(out, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("mmgbsa_dG"):
                    hechos[(r["target"], r["ligand"])] = r
    return hechos


def calcular(r):
    t, lig, smi = r["target"], r["ligand"], r["smiles"]
    pose = os.path.join(POSES, "results_" + t, lig + "_out.pdbqt")
    rec = receptor_de(t)
    base = {"target": t, "ligand": lig, "vina_affinity": r.get("affinity"),
            "plataforma": PLATFORM}
    if not os.path.exists(rec):
        return dict(base, mmgbsa_dG=None, error="sin receptor fijo")
    if not os.path.exists(pose):
        return dict(base, mmgbsa_dG=None, error="sin pose")
    base["receptor_md5"] = r.get("_rec_md5")
    cmd = [sys.executable, HIJO, "--pose", pose, "--receptor", rec, "--fijo",
           "--smiles", smi, "--out", "json"]
    env = dict(os.environ, MMGBSA_FINAL_PLATFORM=PLATFORM,
               MMGBSA_CPU_THREADS=THREADS)
    t0 = time.time()
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT,
                           env=env)
        dt = round(time.time() - t0, 1)
        try:
            row = json.loads(p.stdout.strip().splitlines()[-1])
        except Exception:
            row = None
        if row and row.get("mmgbsa_dG") is not None:
            return dict(base, mmgbsa_dG=row["mmgbsa_dG"],
                        e_complex=row.get("e_complex"),
                        e_receptor=row.get("e_receptor"),
                        e_ligand=row.get("e_ligand"), n_rec=row.get("n_rec"),
                        segundos=dt, error="")
        err = (p.stderr or p.stdout or "salida vacia")[-250:]
        return dict(base, mmgbsa_dG=None, segundos=dt, error=err.strip())
    except subprocess.TimeoutExpired:
        return dict(base, mmgbsa_dG=None, error="TIMEOUT>%ds" % TIMEOUT)
    except Exception as e:
        return dict(base, mmgbsa_dG=None, error=str(e)[:200])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--max", type=int, default=None)
    ap.add_argument("--targets", default=None, help="SOD1,TDP43_v2,FUS")
    ap.add_argument("--tag", default="final")
    ap.add_argument("--controles", action="store_true",
                    help="procesa los controles positivos en vez de la lista")
    a = ap.parse_args()

    targets = [t.strip() for t in a.targets.split(",")] if a.targets else None
    sufijo = "_controles" if a.controles else ""
    out = os.path.join(LOCAL, "rescoring_%s%s.csv" % (a.tag, sufijo))
    logf = os.path.join(LOCAL, "runner_%s%s.log" % (a.tag, sufijo))

    # md5 de los receptores congelados: queda registrado en cada fila.
    rec_md5 = {}
    for t in ["SOD1", "TDP43_v2", "FUS"]:
        p = receptor_de(t)
        rec_md5[t] = md5(p) if os.path.exists(p) else None

    pend_all = (cargar_controles(targets) if a.controles
                else cargar_lista(targets, a.max))
    hechos = cargar_hechos(out)
    pend = [r for r in pend_all if (r["target"], r["ligand"]) not in hechos]
    for r in pend:
        r["_rec_md5"] = rec_md5.get(r["target"])

    log(logf, "PROTOCOLO: plataforma=%s threads=%s tol=%s maxiter=%s"
        % (PLATFORM, THREADS, os.environ.get("MMGBSA_TOL", "1.0"),
           os.environ.get("MMGBSA_MAXITER", "8000")))
    log(logf, "receptores fijos: %s" % json.dumps(rec_md5))
    log(logf, "tarea=%s | total=%d | hechos=%d | pendientes=%d %s"
        % ("controles" if a.controles else "lista", len(pend_all), len(hechos),
           len(pend), dict(Counter(r["target"] for r in pend))))
    if not pend:
        log(logf, "Nada que hacer. DONE.")
        _escribir_info(out, rec_md5, a)
        return

    nuevo = not os.path.exists(out) or os.path.getsize(out) == 0
    outf = open(out, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(outf, fieldnames=COLS, extrasaction="ignore")
    if nuevo:
        w.writeheader()

    ok = err = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        fut = {ex.submit(calcular, r): r for r in pend}
        for i, f in enumerate(fut):
            r = fut[f]
            row = f.result()
            w.writerow(row)
            outf.flush()
            if row.get("mmgbsa_dG") is not None:
                ok += 1
                log(logf, "[%d/%d] %s %s dG=%.2f (%.0fs)"
                    % (i + 1, len(pend), r["target"], r["ligand"],
                       row["mmgbsa_dG"], row.get("segundos") or 0))
            else:
                err += 1
                log(logf, "[%d/%d] %s %s ERROR %s"
                    % (i + 1, len(pend), r["target"], r["ligand"],
                       str(row.get("error", ""))[:120]))
    outf.close()
    log(logf, "FIN: ok=%d err=%d | %.1f min | CSV: %s"
        % (ok, err, (time.time() - t0) / 60.0, out))
    _escribir_info(out, rec_md5, a)


def _escribir_info(out, rec_md5, a):
    info = {"csv": out, "plataforma": PLATFORM, "cpu_threads": THREADS,
            "tolerancia_kj_mol_nm": os.environ.get("MMGBSA_TOL", "1.0"),
            "max_iter": os.environ.get("MMGBSA_MAXITER", "8000"),
            "receptores_md5": rec_md5, "workers": a.workers,
            "controles": bool(a.controles),
            "fin": time.strftime("%Y-%m-%d %H:%M:%S")}
    with open(out + ".info.json", "w", encoding="utf-8") as f:
        json.dump(info, f, ensure_ascii=False, indent=2)


if __name__ == "__main__":
    main()
