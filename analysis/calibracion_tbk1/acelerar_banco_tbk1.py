#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Termina el banco de TBK1 repartiendo el trabajo, en vez de una en una.

POR QUE EXISTE
--------------
El 27 sep 2026 `preparar_banco_tbk1.py` llevaba cinco horas y iba por 3.753 de
35.674 ligandos, con un solo nucleo y atascandose minutos en cada molecula
dificil. De los que faltaban, la mayoria ya estaban preparados en
`gpu_dock/libreria_ligands/` con la receta canonica del proyecto: esos no hay
que construirlos, solo copiarlos. Los que de verdad hay que construir de cero
son unos tres mil, y esos si son lentos.

Este script no inventa receta: importa `preparar_ligando.py`, la copia canonica,
exactamente igual que hace `preparar_banco_tbk1.py`. Lo unico que cambia es que
copia lo que ya existe en vez de rehacerlo, y reparte lo que queda entre varios
nucleos.

Uso:
    python acelerar_banco_tbk1.py --comprobar     # solo cuenta, no toca nada
    python acelerar_banco_tbk1.py                 # copia y genera (6 procesos)
    python acelerar_banco_tbk1.py --procesos 8
    python acelerar_banco_tbk1.py --verificar     # audita lo que hay en ligands/
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys
import time
from multiprocessing import Pool

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(os.path.dirname(BASE))
LIGANDS = os.path.join(BASE, "ligands")
ACTIVOS_CSV = os.path.join(BASE, "compuestos_tbk1.csv")
DECOYS_CSV = os.path.join(BASE, "decoys_tbk1.csv")
INFORME = os.path.join(BASE, "banco_tbk1_acelerado.txt")

# La receta canonica, una sola copia en todo el repo.
sys.path.insert(0, RAIZ)
import preparar_ligando as PL  # noqa: E402
from rdkit import Chem  # noqa: E402

LIBRERIAS = [os.path.join(RAIZ, "gpu_dock", "libreria_ligands"),
             os.path.join(RAIZ, "gpu_dock", "libreria_ligands_reparados")]

# Metales: Vina no tiene tipos para ellos y los acoplaria mal, en silencio.
# Los que sigan en el fragmento organico mayor se apartan y se anotan.
METALES = set(("Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr "
               "Nb Mo Tc Ru Rh Pd Ag Cd In Sn Cs Ba La Ce Pt Au Hg Tl Pb Bi "
               "Gd Eu Tb").split())

LEYENDO = 0


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), flush=True)


def esta_bien(ruta):
    """Un PDBQT del banco esta bien si existe, no esta truncado y no lleva
    pseudo-atomos. Lo del truncado importa: la preparacion anterior se paro a
    mitad y pudo dejar algun fichero cortado."""
    if not os.path.exists(ruta) or os.path.getsize(ruta) < 100:
        return False
    try:
        with open(ruta, encoding="utf-8", errors="ignore") as f:
            txt = f.read()
    except OSError:
        return False
    if "TORSDOF" not in txt:
        return False
    return PL.contar_atomos(txt)[1] == 0


def _copiar(tarea):
    cid, prefijo = tarea
    try:
        destino = os.path.join(LIGANDS, prefijo + cid + ".pdbqt")
        for lib in LIBRERIAS:
            origen = os.path.join(lib, cid + ".pdbqt")
            if not (os.path.exists(origen) and os.path.getsize(origen) > 100):
                continue
            n, glue = PL.contar_fichero(origen)
            if glue:
                # Esta copia no vale, pero la otra libreria puede tenerlo
                # bueno: antes se devolvia fallo aqui y la molecula se caia
                # por el hueco entre la lista de copiar y la de generar.
                continue
            shutil.copy2(origen, destino)
            return (cid, "copiado", "")
    except Exception as e:
        return (cid, "no copiable", "%s: %s" % (type(e).__name__, e))
    return (cid, "no copiable", "")


def _smiles_de_acoplamiento(smi):
    """(smiles, era_sal, metal) de lo que de verdad se acopla.

    El 27 sep 2026 se vio por que el banco no avanzaba: de los 3.018 que no
    estaban en la libreria, 2.169 venian como sal (un `.Cl` pegado) y Meeko
    rechaza cualquier molecula de mas de un fragmento, asi que ETKDG se pasaba
    minutos intentando incrustar el cation y el anion a la vez, molecula tras
    molecula, sin sacar ni una. Se acopla el fragmento organico mayor, que es
    la especie que se une; el contraion no. Las 32.638 que ya estaban en el
    banco son de un solo fragmento, asi que el banco queda consistente.
    """
    mol = Chem.MolFromSmiles(smi) if smi else None
    if mol is None:
        return smi, False, ""
    era_sal = len(Chem.GetMolFrags(mol, sanitizeFrags=False)) > 1
    if era_sal:
        mol = PL.mayor_fragmento(mol)
    metal = ""
    for a in mol.GetAtoms():
        if a.GetSymbol() in METALES:
            metal = a.GetSymbol()
            break
    return (smi, era_sal, metal) if not era_sal else \
        (Chem.MolToSmiles(mol), era_sal, metal)


def _generar(tarea):
    cid, smi, prefijo = tarea
    try:
        smi2, era_sal, metal = _smiles_de_acoplamiento(smi)
        if metal:
            return (cid, "apartado (metal)", metal)
        txt, n, motivo = PL.construir(smi2)
        if txt is None:
            return (cid, "fallo", motivo)
        with open(os.path.join(LIGANDS, prefijo + cid + ".pdbqt"), "w",
                  encoding="utf-8") as f:
            f.write(txt)
    except Exception as e:
        # Una molecula rara no puede tumbar el banco entero, que es lo que
        # paso con la barra de AMPHETAMINE/DEXTROAMPHETAMINE: murio el Pool a
        # los 775 de 2.972 y hubo que empezar la pasada otra vez.
        return (cid, "fallo", "%s: %s" % (type(e).__name__, e))
    return (cid, "generado", "sal recortada" if era_sal else "")


def _auditar(ruta):
    return (os.path.basename(ruta), esta_bien(ruta))


def nombre_seguro(cid):
    """Un id que Windows pueda usar como nombre de fichero.

    El 27 sep 2026 el proceso entero se cayo por un solo decoy, el
    `AMPHETAMINE/DEXTROAMPHETAMINE`: la barra hacia que Windows lo leyera como
    una carpeta y el open() fallaba. En 35.674 ids hay uno asi; mejor
    arreglarlo que confiar en que no vuelva a pasar."""
    return re.sub(r'[\\/:*?"<>|]', "_", cid)


def cargar():
    """[(id, smiles, prefijo)] de activos y decoys, en el orden del banco."""
    fuera = []
    for ruta, prefijo, campo in ((ACTIVOS_CSV, "ACT_", "molecule_chembl_id"),
                                 (DECOYS_CSV, "DEC_", "id")):
        for r in csv.DictReader(open(ruta, encoding="utf-8")):
            fuera.append((nombre_seguro(r[campo]), r["smiles"], prefijo))
    return fuera


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--procesos", type=int, default=6)
    ap.add_argument("--comprobar", action="store_true")
    ap.add_argument("--verificar", action="store_true")
    args = ap.parse_args()

    if args.verificar:
        rutas = [os.path.join(LIGANDS, f) for f in os.listdir(LIGANDS)
                 if f.endswith(".pdbqt")]
        log("auditando %d ficheros..." % len(rutas))
        with Pool(args.procesos) as p:
            res = p.map(_auditar, rutas, chunksize=200)
        malos = [n for n, bien in res if not bien]
        log("auditados %d | malos: %d" % (len(res), len(malos)))
        for n in malos[:40]:
            print("   MAL: %s" % n)
        return 0

    t0 = time.time()
    banco = cargar()
    log("banco: %d ligandos (activos + decoys)" % len(banco))

    faltan = [(cid, smi, pref) for cid, smi, pref in banco
              if not esta_bien(os.path.join(LIGANDS, pref + cid + ".pdbqt"))]
    log("ya estaban bien: %d | faltan: %d" % (len(banco) - len(faltan), len(faltan)))

    a_copiar = [(cid, pref) for cid, smi, pref in faltan
                if any(os.path.exists(os.path.join(lib, cid + ".pdbqt")) and
                       os.path.getsize(os.path.join(lib, cid + ".pdbqt")) > 100
                       for lib in LIBRERIAS)]
    copiados_ya = {cid for cid, _ in a_copiar}
    a_generar = [(cid, smi, pref) for cid, smi, pref in faltan
                 if cid not in copiados_ya]
    log("de los que faltan: %d estan en la libreria (se copian) | %d hay que"
        " construirlos de cero" % (len(a_copiar), len(a_generar)))

    if args.comprobar:
        log("(--comprobar: no se toca nada)")
        return 0

    fallos = []
    sal_recortada = 0
    metales = []
    with Pool(args.procesos) as p:
        if a_copiar:
            log("copiando %d..." % len(a_copiar))
            hechos = 0
            recuperables = []
            for cid, estado, motivo in p.imap_unordered(_copiar, a_copiar, chunksize=50):
                hechos += 1
                if estado != "copiado":
                    fallos.append((cid, estado, motivo))
                    recuperables.append(cid)
                if hechos % 2000 == 0:
                    log("   copiados %d/%d (%.0f s)"
                        % (hechos, len(a_copiar), time.time() - t0))
            log("copia terminada: %d en %.1f min"
                % (len(a_copiar), (time.time() - t0) / 60.0))
            # Lo que no se pudo copiar se construye: si no, se quedaba sin
            # hacer y sin aparecer en ninguna lista.
            if recuperables:
                por_id = {cid: (smi, pref) for cid, smi, pref in faltan}
                extra = [(cid,) + por_id[cid] for cid in recuperables
                         if cid in por_id]
                log("de los no copiables, %d pasan a construirse" % len(extra))
                a_generar = a_generar + extra

        # Las pequeñas primero y las gigantes al final. Los triterpenos de
        # sesenta y pico atomos con veinte estereocentros se llevan media hora
        # cada uno, y al ir delante tenian a los seis nucleos ocupados sin
        # dejar avanzar el resto. Detras, el banco se llena casi entero
        # mientras las dificiles se cuecen en paralelo.
        a_generar.sort(key=lambda t: sum(1 for c in t[1] if c.isupper()) if t[1] else 0)

        if a_generar:
            log("generando %d de cero..." % len(a_generar))
            hechos = 0
            t1 = time.time()
            # chunksize 1: una molecula gigante no puede arrastrar a las otras
            # cuatro de su bloque; en cuanto un nucleo acaba, coge la siguiente.
            for cid, estado, motivo in p.imap_unordered(_generar, a_generar, chunksize=1):
                hechos += 1
                if estado == "generado":
                    if motivo:
                        sal_recortada += 1
                else:
                    fallos.append((cid, estado, motivo))
                    if estado == "apartado (metal)":
                        metales.append(cid)
                if hechos % 25 == 0:
                    ritmo = hechos / max(time.time() - t1, 1e-9) * 60
                    quedan = (len(a_generar) - hechos) / max(ritmo, 1e-9)
                    log("   generados %d/%d | %.0f por minuto | faltan ~%.0f min"
                        % (hechos, len(a_generar), ritmo, quedan))
            log("generacion terminada en %.1f min" % ((time.time() - t1) / 60.0))

    n_pdbqt = len([f for f in os.listdir(LIGANDS) if f.endswith(".pdbqt")])
    # La cuenta que importa: quien del banco sigue sin fichero utilizable, y
    # por que. Contar "pedidos" enganya; esto no.
    sin_fichero = [pref + cid for cid, _, pref in faltan
                   if not esta_bien(os.path.join(LIGANDS, pref + cid + ".pdbqt"))]
    por_que = {}
    for cid, estado, motivo in fallos:
        por_que.setdefault(estado, []).append(cid)
    L = ["BANCO DE TBK1 ACELERADO - %s" % time.strftime("%Y-%m-%d %H:%M"),
         "ligandos en el banco: %d | PDBQT en ligands/: %d"
         % (len(banco), n_pdbqt),
         "ya estaban bien al empezar: %d" % (len(banco) - len(faltan)),
         "copiados de la libreria del proyecto: %d pedidos, %d caidos"
         % (len(a_copiar), len(por_que.get("no copiable", []))),
         "construidos de cero: %d pedidos" % len(a_generar),
         "   de esos, con la sal recortada al fragmento mayor: %d" % sal_recortada,
         "   apartados por llevar metal: %d" % len(metales),
         "",
         "SIGUEN SIN FICHERO: %d" % len(sin_fichero),
         "receta: preparar_ligando.py (la canonica del proyecto)"]
    if metales:
        L.append("")
        L.append("apartados por metal (Vina no los acopla bien):")
        L.append("   " + ", ".join(metales))
    if sin_fichero:
        L.append("")
        L.append("sin fichero, primeros 30:")
        L.append("   " + ", ".join(sin_fichero[:30]))
    with open(INFORME, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print("\n".join(L))
    log("total: %.1f min" % ((time.time() - t0) / 60.0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
