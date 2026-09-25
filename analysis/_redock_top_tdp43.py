# -*- coding: utf-8 -*-
"""Re-docking de los top candidatos TDP-43 de z001 con la caja corregida
(16.3, 41.1, 48.5; 24 A) para ver si el ranking cambia vs la caja vieja."""
import csv, os, subprocess, sys

VINA = r"C:\Users\Fredy\masive-als\tools\vina.exe"
REC = r"C:\Users\Fredy\masive-als\gpu_dock\TDP43.pdbqt"
LOTES = r"C:\Users\Fredy\masive-als\gpu_dock\tanda_z001\ligands"
OUTDIR = r"C:\Users\Fredy\masive-als\analysis\_redock_top_tdp43"
os.makedirs(OUTDIR, exist_ok=True)

# top TDP43 de z001
tops = []
with open(r"C:\Users\Fredy\masive-als\gpu_dock\tanda_z001\resultados_z001.csv") as f:
    for r in csv.DictReader(f):
        if r["target"] == "TDP43":
            tops.append((float(r["affinity"]), r["ligand"]))
tops.sort()
tops = tops[:10]
print("re-dockeando %d top TDP43:" % len(tops))
for aff, lig in tops:
    print("  %s (vieja: %.1f)" % (lig, aff))

# localizar cada pdbqt
import glob
found = {}
for aff, lig in tops:
    hits = glob.glob(os.path.join(LOTES, lig + ".pdbqt"))
    if hits:
        found[lig] = hits[0]
print("localizados:", len(found), "de", len(tops))

res = []
for aff, lig in tops:
    if lig not in found:
        res.append((lig, "NO PDBQT"))
        continue
    out = os.path.join(OUTDIR, lig + "_corr.pdbqt")
    conf = os.path.join(OUTDIR, lig + ".conf")
    with open(conf, "w") as f:
        f.write("receptor = %s\nligand = %s\ncenter_x = 16.3\ncenter_y = 41.1\n"
                "center_z = 48.5\nsize_x = 24\nsize_y = 24\nsize_z = 24\n"
                "exhaustiveness = 16\nnum_modes = 3\nseed = 42\nout = %s\n"
                % (REC.replace("\\", "/"), found[lig].replace("\\", "/"),
                   out.replace("\\", "/")))
    p = subprocess.run([VINA, "--config", conf], capture_output=True, text=True,
                       timeout=240)
    best = None
    for line in p.stdout.splitlines():
        s = line.strip()
        if s and s[0].isdigit() and len(s.split()) >= 2:
            try:
                best = float(s.split()[1])
            except ValueError:
                pass
    res.append((lig, "%.2f" % best if best is not None else "FALLO"))
    print("  %s: %s" % (lig, res[-1][1]))

print("\n=== COMPARATIVA (vieja vs corregida) ===")
for (aff, lig), nuevo in zip(tops, res):
    print("  %-14s vieja=%6.1f  nueva=%s" % (lig, aff, nuevo[1]))
