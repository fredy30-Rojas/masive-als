# -*- coding: utf-8 -*-
"""Re-docking masivo TDP-43 con la caja corregida (16.3, 41.1, 48.5; 24 A).

Toma los hits TDP-43 <= -6.3 (corte top-5% por tanda) de las 10 tandas ya
cribadas, los re-dockea contra el receptor TDP-43 con la caja corregida
centrada en el sitio de union a RNA nativo, y genera el nuevo ranking.
Paralelo (multiprocessing), con checkpoint por tanda para reanudar.
"""
import csv, glob, multiprocessing, os, subprocess, sys, time

VINA = r"C:\Users\Fredy\masive-als\tools\vina.exe"
REC = r"C:\Users\Fredy\masive-als\gpu_dock\TDP43.pdbqt"
CORTE = -6.3
CENTER = (16.3, 41.1, 48.5)
SIZE = 24
EXHAUST = 8
WORKERS = 10  # mitad del CPU (i9-13900H, 20 logicos) para que el PC siga fluido
OUTDIR = r"C:\Users\Fredy\masive-als\analysis\_redock_tdp43_masivo"
CHECKPOINT = os.path.join(OUTDIR, "checkpoint.csv")
RESULTADO = os.path.join(OUTDIR, "resultados_tdp43_corregido.csv")
os.makedirs(OUTDIR, exist_ok=True)

TANDAS = ["z001", "z002", "z003", "z004", "z005",
          "z006", "z007", "z008", "z009", "z010"]


def hits_de_tanda(tanda):
    ruta = r"C:\Users\Fredy\masive-als\gpu_dock\tanda_%s\resultados_%s.csv" % (tanda, tanda)
    ligs = []
    for r in csv.DictReader(open(ruta)):
        if r["target"] == "TDP43" and float(r["affinity"]) <= CORTE:
            ligs.append((r["ligand"], float(r["affinity"])))
    return ligs


def cargar_hechos():
    hechos = set()
    if os.path.exists(CHECKPOINT):
        for r in csv.DictReader(open(CHECKPOINT)):
            hechos.add(r["ligand"])
    return hechos


def dock_uno(args):
    tanda, lig, aff_vieja = args
    lig_pdbqt = r"C:\Users\Fredy\masive-als\gpu_dock\tanda_%s\ligands\%s.pdbqt" % (tanda, lig)
    if not os.path.exists(lig_pdbqt):
        return (lig, tanda, aff_vieja, None, "NO_PDBQT")
    out = os.path.join(OUTDIR, "%s_%s.pdbqt" % (lig, tanda))
    conf = os.path.join(OUTDIR, "%s_%s.conf" % (lig, tanda))
    with open(conf, "w") as f:
        f.write("receptor = %s\nligand = %s\ncenter_x = %.1f\ncenter_y = %.1f\n"
                "center_z = %.1f\nsize_x = %d\nsize_y = %d\nsize_z = %d\n"
                "exhaustiveness = %d\nnum_modes = 3\nseed = 42\nout = %s\n"
                % (REC.replace("\\", "/"), lig_pdbqt.replace("\\", "/"),
                   CENTER[0], CENTER[1], CENTER[2], SIZE, SIZE, SIZE,
                   EXHAUST, out.replace("\\", "/")))
    try:
        p = subprocess.run([VINA, "--config", conf], capture_output=True,
                           text=True, timeout=180)
    except subprocess.TimeoutExpired:
        return (lig, tanda, aff_vieja, None, "TIMEOUT")
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
    t0 = time.time()
    hechos = cargar_hechos()
    print("hechos ya:", len(hechos))

    tareas = []
    for tanda in TANDAS:
        for lig, aff in hits_de_tanda(tanda):
            if lig not in hechos:
                tareas.append((tanda, lig, aff))
    print("pendientes:", len(tareas))

    if not os.path.exists(CHECKPOINT):
        open(CHECKPOINT, "w").write("ligand,tanda,aff_vieja,aff_nueva,error\n")

    filas = []
    with multiprocessing.Pool(WORKERS) as pool:
        for i, res in enumerate(pool.imap_unordered(dock_uno, tareas, chunksize=4)):
            lig, tanda, affv, affn, err = res
            filas.append(res)
            with open(CHECKPOINT, "a") as f:
                f.write("%s,%s,%.3f,%s,%s\n" % (lig, tanda, affv,
                        "%.3f" % affn if affn is not None else "", err or ""))
            if (i + 1) % 50 == 0:
                print("[%d/%d] %.0fs" % (i + 1, len(tareas), time.time() - t0))

    # ranking final
    ok = [r for r in filas if r[3] is not None]
    ok.sort(key=lambda r: r[3])
    with open(RESULTADO, "w") as f:
        f.write("ligand,tanda,aff_vieja,aff_corregida,diferencia\n")
        for lig, tanda, affv, affn, err in ok:
            f.write("%s,%s,%.3f,%.3f,%.3f\n" % (lig, tanda, affv, affn, affn - affv))
    print("\n=== RANKING TDP-43 CORREGIDO (top 25 de %d) ===" % len(ok))
    for lig, tanda, affv, affn, err in ok[:25]:
        print("  %-14s %-5s vieja=%6.2f  corregida=%6.2f  (delta %+.2f)"
              % (lig, tanda, affv, affn, affn - affv))
    print("\ntotal dockeado OK: %d / %d pendientes en %.0fs"
          % (len(ok), len(tareas), time.time() - t0))
    print("resultado -> %s" % RESULTADO)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
