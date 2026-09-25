# -*- coding: utf-8 -*-
"""Cuando el runner v5 termine en Oracle, limpia las filas con error y relanza.

El runner v5 marca como "hecho" cualquier par (ligand,target) presente en el
CSV, incluso con error. Al terminar la tanda de 42, este script:
1. Espera a que el runner termine (proceso desaparece).
2. Borra del CSV las filas con error (mmgbsa_dG vacio).
3. Relanza el runner v5 (solo reprocesara las pendientes).
"""
import subprocess, sys, time

HOST = "ubuntu@79.72.57.253"
KEY = "/c/Users/Fredy/.ssh/oracle_cloud"
REMOTE = "/home/ubuntu/mmgbsa"

def sh(cmd, timeout=30):
    return subprocess.run(["ssh", "-i", KEY, "-o", "BatchMode=yes",
                           "-o", "ConnectTimeout=10", HOST, cmd],
                          capture_output=True, text=True, timeout=timeout)

def log(msg):
    line = "[%s] %s" % (time.strftime("%H:%M:%S"), msg)
    print(line, flush=True)

# 1) esperar a que el runner termine
log("Esperando a que el runner v5 termine...")
quiet = 0
while True:
    p = sh("pgrep -f runner_rescoring_v5.py | grep -v grep | wc -l")
    n = p.stdout.strip()
    if n == "0":
        quiet += 1
        if quiet >= 3:
            break
    else:
        quiet = 0
    time.sleep(60)
log("Runner terminado.")

# 2) limpiar filas con error
p = sh("cd %s && python3 -c \"\nimport csv\nrows=list(csv.DictReader(open('rescoring_mmgbsa_v5.csv')))\nok=[r for r in rows if r['mmgbsa_dG']]\nbad=[r for r in rows if not r['mmgbsa_dG']]\nwith open('rescoring_mmgbsa_v5.csv','w',newline='') as f:\n    w=csv.DictWriter(f, fieldnames=rows[0].keys())\n    w.writeheader()\n    w.writerows(ok)\nprint('ok=%d bad=%d' % (len(ok), len(bad)))\nfor b in bad: print(b['ligand'], b['target'], b['error'][:60])\n\"" % REMOTE)
print(p.stdout)
print(p.stderr)

# 3) relanzar el runner
p = sh("cd %s && export PATH=/home/ubuntu/miniforge3/envs/mmgbsa/bin:$PATH && "
       "nohup python runner_rescoring_v5.py > runner_v5_stdout.log 2>&1 & sleep 5; "
       "tail -3 runner_v5_stdout.log" % REMOTE, timeout=60)
print(p.stdout)
print(p.stderr)
log("Runner v5 relanzado con pendientes.")
