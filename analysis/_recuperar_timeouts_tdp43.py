# -*- coding: utf-8 -*-
"""Recupera los TIMEOUT del re-docking masivo TDP-43 y genera el ranking final.

- Relanza los ligandos que dieron TIMEOUT con un timeout mayor (360 s).
- Con 6-8 workers (para no saturar la CPU) y checkpoint reanudable.
- Al terminar TODOS (éxitos + recuperados), escribe el ranking definitivo
  resultados_tdp43_corregido.csv comparando caja vieja vs corregida.
"""
import csv, multiprocessing, os, subprocess, sys, time

VINA = r"C:\Users\Fredy\masive-als\tools\vina.exe"
REC = r"C:\Users\Fredy\masive-als\gpu_dock\TDP43.pdbqt"
CENTER = (16.3, 41.1, 48.5)
SIZE = 24
EXHAUST = 8
WORKERS = 8
TIMEOUT_S = 360
OUTDIR = r"C:\Users\Fredy\masive-als\analysis\_redock_tdp43_masivo"
CHECKPOINT = os.path.join(OUTDIR, "checkpoint.csv")
RESULTADO = os.path.join(OUTDIR, "resultados_tdp43_corregido.csv")
RECUP_DIR = os.path.join(OUTDIR, "_recuperados")
os.makedirs(RECUP_DIR, exist_ok=True)


def dock_uno(args):
    tanda, lig, aff_vieja = args
    lig_pdbqt = r"C:\Users\Fredy\masive-als\gpu_dock\tanda_%s\ligands\%s.pdbqt" % (tanda, lig)
    if not os.path.exists(lig_pdbqt):
        return (lig, tanda, aff_vieja, None, "NO_PDBQT")
    out = os.path.join(RECUP_DIR, "%s_%s.pdbqt" % (lig, tanda))
    conf = os.path.join(RECUP_DIR, "%s_%s.conf" % (lig, tanda))
    with open(conf, "w") as f:
        f.write("receptor = %s\nligand = %s\ncenter_x = %.1f\ncenter_y = %.1f\n"
                "center_z = %.1f\nsize_x = %d\nsize_y = %d\nsize_z = %d\n"
                "exhaustiveness = %d\nnum_modes = 3\nseed = 42\nout = %s\n"
                % (REC.replace("\\", "/"), lig_pdbqt.replace("\\", "/"),
                   CENTER[0], CENTER[1], CENTER[2], SIZE, SIZE, SIZE,
                   EXHAUST, out.replace("\\", "/")))
    try:
        p = subprocess.run([VINA, "--config", conf], capture_output=True,
                           text=True, timeout=TIMEOUT_S)
    except subprocess.TimeoutExpired:
        return (lig, tanda, aff_vieja, None, "TIMEOUT_360")
    best = None
    for line in p.stdout.splitlines():
        s = line.strip()
        if s and s[0].isdigit() and len(s.split()) >= 2:
            try:
                best = float(s.split()[1])
            except ValueError:
                pass
    return (lig, tanda, aff_vieja, best,
            None if best is not None else (p.stderr or "SIN_SALIDA")[:80])


def main():
    # leer checkpoint actual
    filas = []
    for r in csv.DictReader(open(CHECKPOINT)):
        filas.append({"ligand": r["ligand"], "tanda": r["tanda"],
                      "aff_vieja": float(r["aff_vieja"]),
                      "aff_nueva": (float(r["aff_nueva"]) if r["aff_nueva"] else None),
                      "error": r.get("error", "")})
    hechos = {(r["ligand"], r["tanda"]) for r in filas}
    timeout = [r for r in filas if r["aff_nueva"] is None]
    print("checkpoint total: %d | timeouts a recuperar: %d" % (len(filas), len(timeout)))

    tareas = [(r["tanda"], r["ligand"], r["aff_vieja"]) for r in timeout]
    nuevos = []
    if tareas:
        t0 = time.time()
        with multiprocessing.Pool(WORKERS) as pool:
            for i, res in enumerate(pool.imap_unordered(dock_uno, tareas, chunksize=4)):
                lig, tanda, affv, affn, err = res
                nuevos.append((lig, tanda, affv, affn, err))
                if (i + 1) % 25 == 0:
                    print("[%d/%d] %.0fs" % (i + 1, len(tareas), time.time() - t0))
        # actualizar checkpoint
        with open(CHECKPOINT, "a") as f:
            for lig, tanda, affv, affn, err in nuevos:
                f.write("%s,%s,%.3f,%s,%s\n" % (lig, tanda, affv,
                        "%.3f" % affn if affn is not None else "", err or ""))
        ok_rec = sum(1 for n in nuevos if n[3] is not None)
        print("recuperados: %d/%d" % (ok_rec, len(nuevos)))

    # ranking definitivo con todos los exitos
    ok = [r for r in filas if r["aff_nueva"] is not None]
    ok += [{"ligand": n[0], "tanda": n[1], "aff_vieja": n[2],
            "aff_nueva": n[3], "error": n[4]} for n in nuevos if n[3] is not None]
    ok.sort(key=lambda r: r["aff_nueva"])
    with open(RESULTADO, "w") as f:
        f.write("ligand,tanda,aff_vieja,aff_corregida,diferencia\n")
        for r in ok:
            f.write("%s,%s,%.3f,%.3f,%.3f\n" % (r["ligand"], r["tanda"],
                    r["aff_vieja"], r["aff_nueva"], r["aff_nueva"] - r["aff_vieja"]))
    print("\n=== RANKING TDP-43 CORREGIDO (top 30 de %d) ===" % len(ok))
    for r in ok[:30]:
        print("  %-14s %-5s vieja=%6.2f  corregida=%6.2f  (delta %+.2f)"
              % (r["ligand"], r["tanda"], r["aff_vieja"], r["aff_nueva"],
                 r["aff_nueva"] - r["aff_vieja"]))
    print("\ntotal dockeado OK: %d / %d | sin score: %d"
          % (len(ok), len(filas) + len(nuevos),
             sum(1 for r in filas if r["aff_nueva"] is None) +
             sum(1 for n in nuevos if n[3] is None)))
    print("resultado -> %s" % RESULTADO)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
