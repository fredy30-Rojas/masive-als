# -*- coding: utf-8 -*-
"""Lanzador del rescoring MM-GBSA corregido — receptor FIJO por diana.

Correcciones aplicadas en la revisión del 11/09/2026 (Claude):

  * RECEPTOR CONGELADO: usa receptores_fijos/<diana>_fijo.pdb con --fijo.
    La versión anterior apuntaba a los .pdbqt crudos y sin --fijo, de modo
    que cada ligando volvía a preparar el receptor con PDBFixer (colocación
    de H no determinista) y recortaba la bolsa a su alrededor: las dos
    fuentes de ruido que la auditoría identificó seguían activas.
  * NOMBRE DE DIANA: en la lista corta la diana es TDP43_v2, no TDP43. Con
    la clave antigua los 1.237 candidatos de TDP-43 devolvían "sin
    receptor" y se perdían en silencio.
  * TRAZABILIDAD: cada fila registra el md5 del receptor fijo y la
    plataforma. El valor absoluto de dG depende del artefacto receptor
    (mismo compuesto: -40,30 con uno y -25,86 con otro preparado aparte),
    así que sin el hash el resultado no es reproducible.
  * GPU por defecto: MMGBSA_FINAL_PLATFORM=CUDA. En WSL la plataforma CUDA
    es determinista (DeterministicForces, dispersión 0,00 medida) y con el
    receptor ya relajado tarda ~7 s por compuesto, sin cargar el CPU.

Reanudable: salta los pares (diana, ligando) ya presentes en el CSV.

Uso (desde WSL, con el entorno del rescoring en /root/rescoring_env):
    ~/rescoring_env/bin/python runner_mmgbsa_gb.py --workers 4
    ~/rescoring_env/bin/python runner_mmgbsa_gb.py --targets SOD1 --max 30
    ~/rescoring_env/bin/python runner_mmgbsa_gb.py --lista validacion_controles.csv --tag validacion
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

BASE = r"C:\Users\Fredy\masive-als"
if not os.path.isdir(BASE):                      # dentro de WSL
    BASE = "/mnt/c/Users/Fredy/masive-als"
ANALYSIS = os.path.join(BASE, "analysis")
OUTDIR = os.path.join(ANALYSIS, "rescoring_local")
CAND = os.path.join(ANALYSIS, "lista_corta_candidatos.csv")
RES = os.path.join(BASE, "gpu_dock", "resultados_libreria")
CHILD = os.path.join(OUTDIR, "mmgbsa_openff_gb.py")
FIJO = os.path.join(OUTDIR, "receptores_fijos")
# Una sola bolsa por diana, congelada y con hash registrado (determinista).
TARGETS = ("TDP43_v2", "SOD1", "FUS")
RECEPTORES_FIJOS = {t: os.path.join(FIJO, "%s_fijo.pdb" % t) for t in TARGETS}
COLS = ["target", "ligand", "vina_affinity", "smiles", "mmgbsa_dG",
        "dG_sd", "n_rep", "e_complex", "e_receptor", "e_ligand", "segundos",
        "receptor_md5", "plataforma", "error"]
TIMEOUT = 3600  # 1 h por candidato (con receptor fijo: ~7 s en CUDA)
PLATFORM = os.environ.get("MMGBSA_FINAL_PLATFORM", "CUDA")
# De dónde salen las poses: por defecto las de la librería (results_<diana>/).
# Se puede apuntar a otra carpeta con --poses, con la MISMA estructura dentro.
# Existe (17 sep 2026) para rescoring de poses acopladas en otra caja sin pisar
# las de la librería, que son la prueba de la comparación de cajas.
POSES = RES


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def log(path, msg):
    line = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    print(line, flush=True)
    with open(path, "a", encoding="utf-8") as f:
        f.write(line + "\n")


def load_candidates(targets=None, lista=None):
    rows = []
    with open(lista or CAND, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            t = r.get("target", "")
            if targets and t not in targets:
                continue
            rows.append({"target": t, "ligand": r.get("ligand", ""),
                         "affinity": r.get("afinidad") or r.get("vina_affinity") or "",
                         "smiles": r.get("smiles", "")})
    return rows


def load_done(out):
    done = {}
    if os.path.exists(out):
        with open(out, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("mmgbsa_dG"):
                    done[(r["target"], r["ligand"])] = r
    return done


def calcular_fila(r, hashes, repeats=1, doble=False):
    tgt, lig, smi = r["target"], r["ligand"], r["smiles"]
    pose = os.path.join(POSES, "results_" + tgt, lig + "_out.pdbqt")
    receptor = RECEPTORES_FIJOS.get(tgt)
    base = {"target": tgt, "ligand": lig, "vina_affinity": r.get("affinity"),
            "smiles": smi, "plataforma": PLATFORM, "n_rep": repeats}
    if not receptor:
        return dict(base, mmgbsa_dG=None, error="diana desconocida")
    if not os.path.exists(receptor):
        return dict(base, mmgbsa_dG=None,
                    error="falta receptor fijo %s (correr preparar_receptores.py)"
                          % os.path.basename(receptor))
    base["receptor_md5"] = hashes.get(receptor, "")
    if not os.path.exists(pose):
        return dict(base, mmgbsa_dG=None, error="sin pose")
    # --fijo: el receptor NO se vuelve a preparar (PDBFixer no es determinista)
    cmd = [sys.executable, CHILD, "--pose", pose, "--receptor", receptor,
           "--fijo", "--smiles", smi, "--out", "json"]
    # La plataforma hay que pasarla EXPLÍCITAMENTE al hijo: sin esta variable
    # el hijo usa su propio valor por defecto (CPU) y minimizaría en CPU,
    # aunque el CSV dijera otra cosa.
    env = dict(os.environ, MMGBSA_FINAL_PLATFORM=PLATFORM)
    if doble:
        env["MMGBSA_CUDA_DOUBLE"] = "1"
    # Repeticiones independientes del mismo cálculo: la minimización no es
    # perfectamente reproducible en GPU, así que se reporta la media y su
    # dispersión en vez de un único valor que parece exacto y no lo es.
    dgs, e_comp, e_rec, e_lig, dt_total, plataforma = [], [], [], [], 0.0, ""
    for _ in range(max(1, repeats)):
        t0 = time.time()
        try:
            p = subprocess.run(cmd, capture_output=True, text=True,
                               timeout=TIMEOUT, env=env)
        except subprocess.TimeoutExpired:
            return dict(base, mmgbsa_dG=None, error="TIMEOUT>%ds" % TIMEOUT)
        except Exception as e:
            return dict(base, mmgbsa_dG=None, error=str(e)[:200])
        dt_total += time.time() - t0
        if p.returncode != 0 or not p.stdout.strip():
            err = (p.stderr or p.stdout or "salida vacia")[-250:]
            return dict(base, mmgbsa_dG=None, segundos=round(dt_total, 1),
                        error=err.strip())
        try:
            row = json.loads(p.stdout.strip().splitlines()[-1])
        except json.JSONDecodeError:
            row = None
        if not row or row.get("mmgbsa_dG") is None:
            return dict(base, mmgbsa_dG=None, segundos=round(dt_total, 1),
                        error=str(row)[:200] if row else "salida no JSON")
        dgs.append(row["mmgbsa_dG"])
        e_comp.append(row.get("e_complex"))
        e_rec.append(row.get("e_receptor"))
        e_lig.append(row.get("e_ligand"))
        plataforma = row.get("plataforma") or plataforma
    n = len(dgs)
    media = sum(dgs) / n
    sd = (sum((d - media) ** 2 for d in dgs) / n) ** 0.5
    prom = lambda v: round(sum(v) / len(v), 2) if v and v[0] is not None else None
    return dict(base, mmgbsa_dG=round(media, 2), dG_sd=round(sd, 2),
                e_complex=prom(e_comp), e_receptor=prom(e_rec),
                e_ligand=prom(e_lig), plataforma=plataforma or PLATFORM,
                segundos=round(dt_total, 1), error="")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--max", type=int, default=None)
    ap.add_argument("--targets", default=None,
                    help="lista separada por comas (SOD1,TDP43_v2,FUS)")
    ap.add_argument("--lista", default=None,
                    help="CSV alternativo (mismas columnas que la lista corta)")
    ap.add_argument("--receptores", default=FIJO,
                    help="carpeta con <diana>_fijo.pdb (usar un snapshot "
                         "congelado evita que otra sesión cambie el receptor "
                         "a mitad de la corrida)")
    ap.add_argument("--poses", default=RES,
                    help="carpeta con results_<diana>/<ligando>_out.pdbqt "
                         "(por defecto, las poses de la libreria)")
    ap.add_argument("--tag", default="corregido")
    ap.add_argument("--repeats", type=int, default=1,
                    help="repeticiones independientes por compuesto (se "
                         "reporta media y desviación)")
    ap.add_argument("--double", action="store_true",
                    help="forzar precisión doble en CUDA")
    a = ap.parse_args()

    # Carpeta de receptores fijada al arrancar: el hash de cada fila queda
    # atado al artefacto realmente usado.
    global RECEPTORES_FIJOS, POSES
    RECEPTORES_FIJOS = {t: os.path.join(a.receptores, "%s_fijo.pdb" % t)
                        for t in TARGETS}
    POSES = a.poses

    targets = [t.strip() for t in a.targets.split(",")] if a.targets else None
    out = os.path.join(OUTDIR, "rescoring_%s.csv" % a.tag)
    logf = os.path.join(OUTDIR, "runner_%s.log" % a.tag)

    # Hashes de los receptores fijos: quedan registrados en cada fila.
    hashes = {p: (md5(p) if os.path.exists(p) else "") for p in RECEPTORES_FIJOS.values()}
    for t, p in RECEPTORES_FIJOS.items():
        log(logf, "receptor %-9s %s md5=%s"
            % (t, os.path.basename(p), hashes.get(p, "")[:12] or "AUSENTE"))

    cands = load_candidates(targets, a.lista)
    done = load_done(out)
    pend = [r for r in cands if (r["target"], r["ligand"]) not in done]
    if a.max:
        pend = pend[:a.max]
    log(logf, "lista=%s | plataforma=%s%s | repeticiones=%d | candidatos: %d | ya hechos: %d | pendientes: %d %s"
        % (a.lista or os.path.basename(CAND), PLATFORM,
           " (doble)" if a.double else "", a.repeats, len(cands), len(done),
           len(pend), dict(Counter(r["target"] for r in pend))))
    if not pend:
        log(logf, "Nada que hacer. DONE.")
        return

    new_file = not os.path.exists(out) or os.path.getsize(out) == 0
    outf = open(out, "a", newline="", encoding="utf-8")
    w = csv.DictWriter(outf, fieldnames=COLS, extrasaction="ignore")
    if new_file:
        w.writeheader()

    ok = err = 0
    t_start = time.time()
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        fut = {ex.submit(calcular_fila, r, hashes, a.repeats, a.double): r
               for r in pend}
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
        % (ok, err, (time.time() - t_start) / 60.0, out))


if __name__ == "__main__":
    main()
