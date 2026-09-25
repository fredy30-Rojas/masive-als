#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""¿Repite el MM-GBSA el mismo numero si se le llama dos veces con lo mismo?

Corre 3 veces el script viejo y 3 veces el nuevo sobre ACT_isoproterenol (row 0)
y saca los tres terminos de energia, para ver cual se mueve.
"""
import csv
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

VIEJO = "rescoring_mmgbsa_robusto.py.bak_20260922"
NUEVO = "rescoring_mmgbsa_robusto.py"
CSV = "estratos/candidatos.csv"


def correr(script, rep):
    cmd = [sys.executable, script, "--candidates", CSV,
           "--receptores", "receptores_estratos.json", "--row", "0"]
    t0 = time.time()
    r = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    lineas = [l for l in r.stdout.splitlines() if l.startswith("{")]
    if not lineas:
        return f"{script[-12:]} rep{rep} SIN SALIDA: {(r.stderr or '')[-160:]}"
    d = json.loads(lineas[-1])
    return (f"{script[-12:]} rep{rep} dG={d.get('mmgbsa_dG')} "
            f"E_complejo={d.get('e_complex')} E_recep={d.get('e_receptor')} "
            f"E_lig={d.get('e_ligand')} err={str(d.get('error'))[:40]} "
            f"({dt:.0f} s)")


def plataforma():
    import openmm
    n = openmm.Platform.getNumPlatforms()
    nombres = [openmm.Platform.getPlatform(i).getName() for i in range(n)]
    return f"openmm {openmm.version.version} plataformas={nombres}"


def main():
    print(plataforma(), flush=True)
    trabajos = [(VIEJO, 1), (VIEJO, 2), (VIEJO, 3),
                (NUEVO, 1), (NUEVO, 2), (NUEVO, 3)]
    with ThreadPoolExecutor(max_workers=3) as ex:
        for res in ex.map(lambda t: correr(*t), trabajos):
            print(res, flush=True)


if __name__ == "__main__":
    main()
