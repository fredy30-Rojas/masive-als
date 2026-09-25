# -*- coding: utf-8 -*-
"""Vigilante del rescoring MM-GBSA de la lista focalizada (17 sep 2026).

Que hace:
  Cada 5 minutos mira el CSV del rescoring y avisa UNA sola vez:
    * cuando ya se han intentado las 179 moleculas de la lista -> parte con el
      recuento y los mejores dG;
    * cuando el runner deja de estar vivo y el CSV se queda quieto -> aviso de
      que se ha parado, con por donde iba.

Por que asi:
  Los 7 compuestos que fallan por desajuste de atomos NUNCA van a tener dG, asi
  que el fin no se puede detectar contando dG. Se detecta contando cuantas
  moleculas de la lista tienen ya su fila en el CSV, tengan nota o error.

Uso: pythonw vigilar_rescoring_focalizada.py   (sin ventana, ver el .vbs)
Estado: vigilante_rescoring_estado.json   Registro: vigilante_rescoring.log
"""
import csv
import json
import os
import subprocess
import sys
import time
import urllib.parse
import urllib.request

AQUI = os.path.dirname(os.path.abspath(__file__))
LISTA = os.path.join(AQUI, "lista_focalizada_rescoring.csv")
CSV = os.path.join(AQUI, "rescoring_focalizada_arn.csv")
LOG_RUNNER = os.path.join(AQUI, "runner_focalizada_arn.log")
LOG = os.path.join(AQUI, "vigilante_rescoring.log")
ESTADO = os.path.join(AQUI, "vigilante_rescoring_estado.json")
TELEGRAM_JSON = r"C:/Users/Fredy/.claude/secure/telegram-bot.json"
HABLAR = r"C:/Users/Fredy/hablar.py"
# El TTS se lanza con el python de verdad, no con pythonw: hablar.py escribe
# errores por stderr y con pythonw eso revienta. El vigilante si va sin ventana.
PYTHON_TTS = r"C:/Users/Fredy/AppData/Local/Python/pythoncore-3.14-64/python.exe"
if not os.path.exists(PYTHON_TTS):
    PYTHON_TTS = sys.executable

# En Windows, un hijo de consola (wsl.exe, el TTS) se abre su propia ventana
# cuando el padre no tiene consola. Ese era el fogonazo cada 5 minutos.
SIN_VENTANA = getattr(subprocess, "CREATE_NO_WINDOW", 0)

INTERVALO = 300          # 5 minutos
QUIETO_MIN = 25          # minutos sin escribir en el CSV para darlo por parado


def log(msg):
    linea = "[%s] %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), msg)
    try:
        with open(LOG, "a", encoding="utf-8") as f:
            f.write(linea + "\n")
    except Exception:
        pass


def leer_json(path, defecto):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return defecto


def guardar_json(path, datos):
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(datos, f, ensure_ascii=False, indent=2)
    except Exception as e:
        log("No pude guardar %s: %s" % (os.path.basename(path), e))


def telegram(texto):
    """Devuelve True si el parte salio."""
    cred = leer_json(TELEGRAM_JSON, {})
    token, chat = cred.get("token"), cred.get("chat_id")
    if not token or not chat:
        log("Sin credenciales de Telegram")
        return False
    datos = urllib.parse.urlencode({"chat_id": chat, "text": texto}).encode()
    url = "https://api.telegram.org/bot%s/sendMessage" % token
    try:
        with urllib.request.urlopen(url, data=datos, timeout=20) as r:
            return r.status == 200
    except Exception as e:
        log("Telegram error: %s" % e)
        return False


def hablar(texto):
    try:
        subprocess.Popen([PYTHON_TTS, HABLAR, texto],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                         creationflags=SIN_VENTANA)
    except Exception as e:
        log("TTS error: %s" % e)


def avisar(texto_telegram, texto_voz):
    log("AVISO: %s" % texto_telegram)
    log("Telegram OK=%s" % telegram(texto_telegram))
    hablar(texto_voz)


def candidatos():
    try:
        with open(LISTA, encoding="utf-8") as f:
            return [r["ligand"] for r in csv.DictReader(f) if r.get("ligand")]
    except Exception as e:
        log("No pude leer la lista: %s" % e)
        return []


def intentados():
    """{ligando: fila} de todo lo que ya tiene linea en el CSV."""
    hechos = {}
    if not os.path.exists(CSV):
        return hechos
    try:
        with open(CSV, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r.get("ligand"):
                    hechos[r["ligand"]] = r
    except Exception as e:
        log("No pude leer el CSV: %s" % e)
    return hechos


def runner_vivo():
    """True si el runner sigue corriendo dentro de WSL.

    Los corchetes del patron son defensivos: evitan que el bash que lanza el
    pgrep (si el shell no hiciera la optimizacion de exec) se encuentre a si
    mismo y conteste siempre que si.

    Comprobado el 19 sep 2026 a mano: con bash -lc en este lanzador el
    resultado es identico con y sin corchetes, porque bash exec'ea el pgrep y
    el proceso del shell desaparece. Es decir: aquí NO habia tal fallo. Se
    dejan puestos por seguridad, no porque arreglen nada.
    """
    for cmd in (["wsl.exe", "-d", "Ubuntu", "-u", "root", "bash", "-lc",
                 "pgrep -f '[r]unner_mmgbsa_gb.py'"],
                ["wsl.exe", "-d", "Ubuntu", "bash", "-lc",
                 "pgrep -f '[r]unner_mmgbsa_gb.py'"]):
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=60,
                               creationflags=SIN_VENTANA)
            if r.returncode == 0 and r.stdout.strip():
                return True
            return False
        except Exception:
            continue
    return False


def segundos_quieto():
    try:
        return time.time() - os.path.getmtime(CSV)
    except Exception:
        return 0


def mejores(hechos, n=5):
    filas = [r for r in hechos.values() if (r.get("mmgbsa_dG") or "").strip()]
    def val(r):
        try:
            return float(r["mmgbsa_dG"])
        except Exception:
            return 0.0
    filas.sort(key=val)
    return [(r["ligand"], val(r)) for r in filas[:n]]


def bucle():
    total = len(candidatos())
    log("Vigilante arrancado | lista: %d moleculas | intervalo %ds" % (total, INTERVALO))
    if not total:
        return 1
    while True:
        est = leer_json(ESTADO, {})
        hechos = intentados()
        n = len(hechos)
        con_dg = sum(1 for r in hechos.values() if (r.get("mmgbsa_dG") or "").strip())
        errores = n - con_dg

        if n >= total:
            if not est.get("avisado_fin"):
                lista_mejores = ", ".join("%s %.2f" % (l, d) for l, d in mejores(hechos))
                avisar(
                    "MASIVE-ALS rescoring focalizada TERMINADO\n"
                    "%d de %d moleculas intentadas\n"
                    "Con dG: %d | Fallos: %d\n"
                    "Mejores: %s" % (n, total, con_dg, errores, lista_mejores),
                    "Fredy, ya termino el rescoring de la lista focalizada. "
                    "Salieron %d con nota buena de las %d, y %d se quedaron fuera por "
                    "el desajuste de atomos. Los mejores son %s. "
                    "Diga usted y le preparo el ranking completo." % (con_dg, total, errores, lista_mejores))
                est["avisado_fin"] = time.strftime("%Y-%m-%d %H:%M:%S")
                guardar_json(ESTADO, est)
            return 0

        vivo = runner_vivo()
        quieto_min = segundos_quieto() / 60.0
        if not vivo and quieto_min >= QUIETO_MIN and not est.get("avisado_parado"):
            avisar(
                "MASIVE-ALS rescoring focalizada PARADO\n"
                "%d de %d moleculas intentadas (%d con dG, %d fallos)\n"
                "Sin proceso vivo y el CSV lleva %.0f min quieto" % (n, total, con_dg, errores, quieto_min),
                "Fredy, se paro el rescoring de la lista focalizada. Iba por %d de %d, "
                "y no hay ningun proceso vivo, asi que se corto. Digame y lo relanzo "
                "donde lo dejo, que no se pierde nada." % (n, total))
            est["avisado_parado"] = time.strftime("%Y-%m-%d %H:%M:%S")
            guardar_json(ESTADO, est)

        est["ultimo_chequeo"] = time.strftime("%Y-%m-%d %H:%M:%S")
        est["intentados"] = n
        est["con_dg"] = con_dg
        est["runner_vivo"] = vivo
        guardar_json(ESTADO, est)
        time.sleep(INTERVALO)


if __name__ == "__main__":
    raise SystemExit(bucle())
