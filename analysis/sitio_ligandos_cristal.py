#!/usr/bin/env python3
"""Mide en que sitio de SOD1 se asienta cada ligando co-cristalizado.

Por que existe: un ligando co-cristalizado solo vale como positivo si se asienta en
el MISMO sitio que la caja del acoplamiento del proyecto (el bolsillo del Trp32). Y eso
no hay que suponerlo ni acoplarlo: esta escrito en las coordenadas depositadas. Aqui
se leen, sin acoplar nada y sin GPU.

Lo que se mide, por ligando y estructura:
  - residuos de proteina a menos de 4,0 A del ligando,
  - si el Trp32 esta entre ellos y a que distancia minima queda,
  - si toca el sitio metalico (His46, His48, His63, His71, His80, Asp83) o la Cys111.

Clasificacion: Trp32 / metalico / Cys111 / otro (y las mezclas).

Entradas:  analysis/ligandos_cristal/sod1.csv (lo deja buscar_ligandos_cristal.py)
Salidas:   analysis/ligandos_cristal/sitio_sod1.csv y sitio_sod1.txt
Cache:     analysis/ligandos_cristal/_pdb/{PDB}.pdb

Uso:
  python analysis/sitio_ligandos_cristal.py
"""
from __future__ import annotations

import csv
import math
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

SALIDA = Path(__file__).resolve().parent / "ligandos_cristal"
PDFILES = "https://files.rcsb.org/download/{}.pdb"
CORTE = 4.0  # angstrom

# Numeracion de la SOD1 humana (1-153), la que usan las estructuras del PDB.
TRP32 = {("TRP", 32)}
METAL = {("HIS", 46), ("HIS", 48), ("HIS", 63), ("HIS", 71), ("HIS", 80), ("ASP", 83)}
CYS111 = {("CYS", 111)}


def descargar(pdb: str) -> Path | None:
    SALIDA.joinpath("_pdb").mkdir(parents=True, exist_ok=True)
    destino = SALIDA / "_pdb" / f"{pdb}.pdb"
    if destino.exists() and destino.stat().st_size > 1000:
        return destino
    for i in range(3):
        try:
            with urllib.request.urlopen(PDFILES.format(pdb), timeout=60) as r:
                datos = r.read()
            destino.write_bytes(datos)
            return destino
        except (urllib.error.URLError, TimeoutError) as e:
            if i == 2:
                print(f"    {pdb}: no se pudo descargar ({e})")
                return None
            time.sleep(3)
    return None


def leer(fichero: Path, comp: str):
    """Atomos de proteina y del ligando (HETATM de ese compuesto)."""
    proteina, ligando = [], []
    with open(fichero, encoding="ascii", errors="ignore") as fh:
        for linea in fh:
            if linea.startswith("HETATM"):
                if linea[17:20].strip() == comp or linea[17:21].strip() == comp:
                    ligando.append(
                        (float(linea[30:38]), float(linea[38:46]), float(linea[46:54]))
                    )
            elif linea.startswith("ATOM"):
                proteina.append(
                    (
                        linea[17:20].strip(),
                        int(linea[22:26]),
                        float(linea[30:38]),
                        float(linea[38:46]),
                        float(linea[46:54]),
                    )
                )
    return proteina, ligando


def contactos(proteina, ligando):
    """Residuos a menos de CORTE angstrom de cualquier atomo del ligando."""
    cerca, minimos = set(), {}
    for res, num, x, y, z in proteina:
        mejor = min(
            math.dist((x, y, z), p) for p in ligando
        )
        if mejor <= CORTE:
            cerca.add((res, num))
            minimos[(res, num)] = round(mejor, 2)
    return cerca, minimos


def clasifica(cerca: set) -> str:
    toca = []
    if cerca & TRP32:
        toca.append("Trp32")
    if cerca & METAL:
        toca.append("metalico")
    if cerca & CYS111:
        toca.append("Cys111")
    return "+".join(toca) if toca else "otro"


def main() -> int:
    entrada = SALIDA / "sod1.csv"
    if not entrada.exists():
        print(f"Falta {entrada}: corre antes buscar_ligandos_cristal.py")
        return 1

    with open(entrada, encoding="utf-8") as fh:
        filas = list(csv.DictReader(fh))

    out = []
    for fila in filas:
        comp = fila["comp_id"]
        pdb = fila["pdbs"].split(";")[0]
        fichero = descargar(pdb)
        if fichero is None:
            continue
        proteina, ligando = leer(fichero, comp)
        if not ligando:
            print(f"  {comp:<5} {pdb}: no hay HETATM de ese compuesto (¿4 caracteres?)")
            continue
        cerca, minimos = contactos(proteina, ligando)
        sitio = clasifica(cerca)
        d_trp = min(
            (minimos[k] for k in minimos if k in TRP32), default=""
        )
        residuos = ";".join(f"{r}{n}" for r, n in sorted(cerca, key=lambda x: x[1]))
        out.append(
            {
                "diana": fila["diana"],
                "comp_id": comp,
                "nombre": fila["nombre"],
                "pdbs": fila["pdbs"],
                "pesados": fila["pesados"],
                "ya_en_verdad": fila["ya_en_verdad"],
                "n_contactos": len(cerca),
                "sitio": sitio,
                "dist_min_trp32": d_trp,
                "residuos_4A": residuos[:400],
            }
        )
        print(f"  {comp:<5} {fila['ya_en_verdad'][:18]:<18} {sitio:<18} {len(cerca):>3} residuos  d(Trp32)={d_trp}")

    out.sort(key=lambda f: (f["sitio"].startswith("otro"), f["ya_en_verdad"] != "no"))
    with open(SALIDA / "sitio_sod1.csv", "w", newline="", encoding="utf-8") as fh:
        if out:
            w = csv.DictWriter(fh, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)

    texto = [
        "SITIO REAL DE CADA LIGANDO CO-CRISTALIZADO EN SOD1 (leido de las coordenadas, 25 sep 2026)",
        "Contacto = residuo de proteina a menos de 4,0 A del ligando. Numeracion humana 1-153.",
        "",
    ]
    for f in out:
        texto.append(
            f"{f['comp_id']:<5} {f['sitio']:<18} {f['pesados']:>3} pesados  "
            f"{f['pdbs'][:20]:<20} {f['nombre'][:44]:<44} {'EN VERDAD: ' + f['ya_en_verdad'] if f['ya_en_verdad'] != 'no' else 'nuevo'}"
        )
    (SALIDA / "sitio_sod1.txt").write_text("\n".join(texto) + "\n", encoding="utf-8")
    print("\n" + "\n".join(texto))
    return 0


if __name__ == "__main__":
    sys.exit(main())
