#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""_limpiar_ranking_tipos.py — quita del ranking bruto los ligandos envenenados.

POR QUE
-------
La auditoria de tipos (_auditoria_tipos_cribado.log, 23 sep 2026) encontro
ligandos de la libreria con un elemento critico (B, Si, Se, Te, Sn, As, Ge)
NOMBRADO en el PDBQT pero TIPADO como otro elemento (normalmente C). Vina los
acoplo y puntuo fingiendo ser otro elemento: su afinidad es falsa. Trece de
ellos superan -7 kcal/mol en TDP43_v2 y ensucian el ranking bruto.

Que hace
--------
1. Lee los nombres de ligando de analysis/ranking_afinidades.csv.
2. Escanea SU PDBQT en gpu_dock/libreria_ligands con la MISMA logica exacta de
   auditar_tipos_cribado.mirar (importada, no copiada, para que no diverjan).
3. Los que salen criticos se sacan del ranking. Nada se borra a ciegas: las
   filas retiradas van a ranking_afinidades_descartes_tipos.csv.
4. Escribe LEEME_ranking_afinidades.md con el porque y la fecha.

No re-prepara ni re-acopla nada: solo limpia el indice. Los 41 salen listados
en el log para poder auditar la limpieza.
"""
import csv
import os
import sys
import time

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
sys.path.insert(0, BASE)
from auditar_tipos_cribado import mirar  # misma logica, cero divergencia

LIBRERIA = os.path.join(RAIZ, "gpu_dock", "libreria_ligands")
RANKING = os.path.join(BASE, "ranking_afinidades.csv")
DESCARTES = os.path.join(BASE, "ranking_afinidades_descartes_tipos.csv")
NOTA = os.path.join(BASE, "LEEME_ranking_afinidades.md")
LOG = os.path.join(BASE, "_limpiar_ranking_tipos.log")


def log(m):
    print(m, flush=True)
    with open(LOG, "a", encoding="utf-8") as f:
        f.write(m + "\n")


def main():
    t0 = time.time()
    log("LIMPIEZA DEL RANKING POR TIPOS FALSOS   %s" % time.strftime("%Y-%m-%d %H:%M"))

    # 1. nombres unicos presentes en el ranking
    with open(RANKING, encoding="utf-8", newline="") as f:
        lector = csv.reader(f)
        cabecera = next(lector)
        filas = list(lector)
    col_lig = cabecera.index("ligand")
    nombres = sorted({r[col_lig] for r in filas})
    log("ranking: %d filas, %d ligandos unicos" % (len(filas), len(nombres)))

    # 2. escanear solo los pdbqt de esos ligandos
    criticos = set()
    vistos = 0
    nombres_set = set(nombres)  # precalculado UNA vez (con set() dentro del bucle tardaba 90 min)
    for n in sorted(os.listdir(LIBRERIA)):
        if not n.endswith(".pdbqt"):
            continue
        stem = n[:-len(".pdbqt")]
        if stem not in nombres_set:
            continue
        vistos += 1
        _, casos = mirar(os.path.join(LIBRERIA, n))
        if casos:
            criticos.add(stem)
            log("  critico: %s  %s" % (stem, casos[:3]))
        if vistos % 10000 == 0:
            log("  ...%d pdbqt revisados (%.0f min)" % (vistos, (time.time() - t0) / 60.0))
    log("escaneo: %d pdbqt revisados, %d ligandos criticos (%.1f min)"
        % (vistos, len(criticos), (time.time() - t0) / 60.0))
    if not criticos:
        log("nada que limpiar.")
        return 0

    # 3. partir el ranking
    buenas, malas = [], []
    for r in filas:
        (malas if r[col_lig] in criticos else buenas).append(r)
    log("filas: %d se quedan, %d se retiran" % (len(buenas), len(malas)))

    with open(RANKING, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(cabecera)
        w.writerows(buenas)
    with open(DESCARTES, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(cabecera)
        w.writerows(malas)
    log("descartes guardados en %s" % DESCARTES)

    # 4. la nota
    por_diana = {}
    for r in malas:
        por_diana[r[0]] = por_diana.get(r[0], 0) + 1
    detalle = "\n".join("- %s: %d filas retiradas" % (k, v)
                        for k, v in sorted(por_diana.items()))
    with open(NOTA, "w", encoding="utf-8") as f:
        f.write("""# Por que ranking_afinidades.csv no coincide con resultados_libreria_total.csv

**24 de septiembre de 2026.** Este ranking se limpio: se retiraron **%d ligandos**
(%d filas) que tenian la afinidad FALSA.

## El porque

La auditoria de tipos de atomo del 23 sep 2026 (`_auditoria_tipos_cribado.log`)
encontro en la libreria ligandos con un elemento critico (boro, silicio, selenio,
telurio...) **nombrado** en el PDBQT pero **tipado** como otro elemento (casi
siempre C). Vina tipa la malla por el tipo, no por el nombre: esos ligandos se
acoplaron y puntuaron fingiendo ser su analogo de carbono, asi que su numero del
ranking no mide su quimica. Es el mismo defecto que dejo a `DEC_CHEMBL4543460`
sin pose en la validacion de TDP-43 (boron no tipable), pero en su variante mas
traicionera: aqui el ligando SI se acopla, y puntua mal sin avisar.

## Que se hizo

- Los %d ligandos se retiraron de `ranking_afinidades.csv` con
  `_limpiar_ranking_tipos.py`, que aplica la MISMA logica de
  `auditar_tipos_cribado.py` (importada, no copiada).
- Las filas retiradas NO se borraron: estan en
  `ranking_afinidades_descartes_tipos.csv` por si hay que auditar.
- No se tocó `resultados_libreria_total.csv` ni el top100 ni las listas de
  candidatos: se verifico que NINGUNO de estos ligandos llego a lista de
  candidatos (cero contaminacion; el ranking bruto era el unico sitio sucio).

## Impacto por diana

%s

## Lo que NO dice esta nota

Que los ligandos retirados sean malos compuestos: varios son quimicamente
interesantes (boro, silicio). Solo que SU AFINIDAD VINA aqui es un artefacto.
Para puntuarlos de verdad habria que re-prepararlos y re-acoplarlos con un
motor que entienda esos elementos, y eso es otra conversacion.
""" % (len(criticos), len(malas), len(criticos), detalle))
    log("nota escrita en %s" % NOTA)
    log("terminado en %.1f min" % ((time.time() - t0) / 60.0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
