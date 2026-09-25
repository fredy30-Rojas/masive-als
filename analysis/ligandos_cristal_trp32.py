#!/usr/bin/env python3
"""Busca en el PDB TODOS los ligandos cristalográficos que tocan Trp32 de SOD1.

Motivo: el único control positivo con evidencia directa en el bolsillo Trp32 era
el puñado de estructuras ya conocidas (4A7S/T/U/V). La forma sistemática de
encontrar más —y de encontrar QUIMIAS independientes— es recorrer todas las
estructuras de SOD1 humana depositadas y medir qué moléculas pequeñas (HETATM)
están a menos de N Å de los átomos del anillo de Trp32.

No hace falta superponer nada: cada fichero trae su propia Trp32 y su propio
ligando en el mismo sistema de coordenadas.

Salida en analysis/trp32_cristal/:
  ligandos_trp32.csv   un registro por (estructura, ligando) con contactos
  resumen.txt          lo legible, agrupado por química

Uso:
  python analysis/ligandos_cristal_trp32.py [--corte 5.0]
"""
from __future__ import annotations

import argparse
import csv
import io
import json
import math
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

PDBE = "https://www.ebi.ac.uk/pdbe/search/pdb/select"
RCSB_FILES = "https://files.rcsb.org/download/{}.pdb"
RCSB_CHEM = "https://data.rcsb.org/rest/v1/core/chemcomp/{}"
OUT = Path(__file__).resolve().parent / "trp32_cristal"
CACHE = OUT / "pdb"

# Disolvente, iones y aditivos de cristalización típicos: no son ligandos
IGNORAR = {
    "HOH", "DOD", "ZN", "CU", "NA", "K", "CL", "MG", "CA", "MN", "NI", "CD",
    "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "MPD", "ACT", "ACY", "TRS", "MES",
    "EPE", "IMD", "DMS", "FMT", "CIT", "TLA", "NO3", "AZI", "SCN", "BR", "IOD",
    "F", "NH4", "BME", "DTT", "TCE", "2PE", "1PE", "P6G", "PG4", "PE4", "SIN",
    "BCT", "CO3", "GLY", "SER", "EOH", "MOH", "IPA", "ACE", "PCA", "UNX", "UNL",
}


def texto(url: str, intentos: int = 3, timeout: int = 60) -> str:
    for i in range(intentos):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "masive-als/1.0"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:  # noqa: BLE001
            if i == intentos - 1:
                raise
            import time
            time.sleep(2 + i)
    return ""


def lista_estructuras() -> list[str]:
    url = f"{PDBE}?q=uniprot_accession%3AP00441&wt=json&rows=500&fl=pdb_id"
    d = json.loads(texto(url))
    ids = sorted({doc["pdb_id"].lower() for doc in d["response"]["docs"]})
    return ids


def descargar(pid: str) -> tuple[str, str]:
    CACHE.mkdir(parents=True, exist_ok=True)
    f = CACHE / f"{pid}.pdb"
    if f.exists() and f.stat().st_size > 1000:
        return pid, f.read_text(encoding="utf-8", errors="replace")
    try:
        t = texto(RCSB_FILES.format(pid))
    except Exception:  # noqa: BLE001
        return pid, ""
    if t and ("ATOM" in t or "HETATM" in t):
        f.write_text(t, encoding="utf-8")
        return pid, t
    return pid, ""


def parsear(txt: str):
    """Devuelve (atomos_trp32, hetatms) con coordenadas."""
    trp, het = [], []
    for linea in txt.splitlines():
        if linea.startswith("ATOM") and linea[17:20].strip() == "TRP":
            try:
                if int(linea[22:26]) != 32:
                    continue
            except ValueError:
                continue
            elemento = (linea[76:78].strip() or linea[12:16].strip()[0]).upper()
            if elemento in ("C", "N"):
                trp.append(
                    (linea[12:16].strip(), float(linea[30:38]), float(linea[38:46]), float(linea[46:54]))
                )
        elif linea.startswith("HETATM"):
            res = linea[17:20].strip().upper()
            if res in IGNORAR:
                continue
            elemento = (linea[76:78].strip() or linea[12:16].strip()[0]).upper()
            if elemento in ("H", "D"):
                continue
            try:
                het.append(
                    (
                        res,
                        linea[16].strip(),
                        int(linea[22:26]),
                        linea[12:16].strip(),
                        float(linea[30:38]),
                        float(linea[38:46]),
                        float(linea[46:54]),
                    )
                )
            except ValueError:
                continue
    return trp, het


def analizar(pid: str, txt: str, corte: float) -> list[dict]:
    trp, het = parsear(txt)
    if not trp or not het:
        return []
    out = []
    por_lig: dict[tuple, dict] = {}
    for res, cade, num, atomo, x, y, z in het:
        dmin = min(
            math.dist((x, y, z), (tx, ty, tz)) for _, tx, ty, tz in trp
        )
        if dmin > corte:
            continue
        clave = (res, cade, num)
        r = por_lig.setdefault(clave, {"res": res, "cadena": cade, "num": num, "dmin": 99.0, "n": 0})
        r["n"] += 1
        r["dmin"] = min(r["dmin"], dmin)
    for clave, r in por_lig.items():
        out.append({"pdb": pid, "ligando": r["res"], "cadena": r["cadena"], "num": r["num"],
                    "atomos_cerca": r["n"], "dist_min_A": round(r["dmin"], 2)})
    return out


def quimica(codigos: list[str]) -> dict:
    """SMILES, fórmula y descriptor RCSB de cada código químico."""
    out = {}
    for i, cc in enumerate(codigos):
        try:
            d = json.loads(texto(RCSB_CHEM.format(cc), timeout=30))
            desc = d.get("rcsb_chem_comp_descriptor") or {}
            out[cc] = {
                "nombre": (d.get("chem_comp") or {}).get("name", ""),
                "formula": (d.get("chem_comp") or {}).get("formula", ""),
                "peso": (d.get("chem_comp") or {}).get("formula_weight", ""),
                "smiles": desc.get("SMILES_stereo") or desc.get("SMILES") or "",
                "inchi_key": desc.get("InChIKey", ""),
            }
        except Exception:  # noqa: BLE001
            out[cc] = {"nombre": "", "formula": "", "peso": "", "smiles": "", "inchi_key": ""}
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--corte", type=float, default=5.0, help="distancia maxima a Trp32 (A)")
    args = ap.parse_args(argv)

    OUT.mkdir(parents=True, exist_ok=True)
    ids = lista_estructuras()
    print(f"Estructuras de SOD1 humana en el PDB: {len(ids)}", flush=True)

    registros = []
    with ThreadPoolExecutor(max_workers=8) as ex:
        for pid, txt in ex.map(descargar, ids):
            if not txt:
                print(f"  {pid}: sin fichero PDB (¿solo mmCIF?)", flush=True)
                continue
            registros.extend(analizar(pid, txt, args.corte))
    print(f"Contactos HETATM-Trp32 (<= {args.corte} A): {len(registros)}", flush=True)

    codigos = sorted({r["ligando"] for r in registros})
    info = quimica(codigos)

    with open(OUT / "ligandos_trp32.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(
            fh,
            fieldnames=["pdb", "ligando", "cadena", "num", "atomos_cerca", "dist_min_A",
                        "nombre", "formula", "peso", "smiles"],
        )
        w.writeheader()
        for r in sorted(registros, key=lambda r: (r["ligando"], r["pdb"])):
            q = info.get(r["ligando"], {})
            w.writerow({**r, "nombre": q.get("nombre", ""), "formula": q.get("formula", ""),
                        "peso": q.get("peso", ""), "smiles": q.get("smiles", "")})

    por_lig: dict[str, list[dict]] = {}
    for r in registros:
        por_lig.setdefault(r["ligando"], []).append(r)

    lineas = [
        "LIGANDOS CRISTALOGRAFICOS QUE TOCAN Trp32 DE SOD1 (todas las estructuras P00441)",
        f"Estructuras recorridas: {len(ids)}   Contactos: {len(registros)}   "
        f"Codigos distintos: {len(por_lig)}   Corte: {args.corte} A",
        "",
        f"{'cod':<5} {'n_estr':>6} {'dmin':>6}  {'peso':>7}  {'nombre':<44} formula",
    ]
    for cc, lista in sorted(por_lig.items(), key=lambda kv: -len(kv[1])):
        q = info.get(cc, {})
        dmin = min(r["dist_min_A"] for r in lista)
        lineas.append(
            f"{cc:<5} {len({r['pdb'] for r in lista}):>6} {dmin:>6.2f}  "
            f"{str(q.get('peso', ''))[:7]:>7}  {(q.get('nombre') or '')[:44]:<44} {q.get('formula','')}"
        )
    texto_salida = "\n".join(lineas) + "\n"
    (OUT / "resumen.txt").write_text(texto_salida, encoding="utf-8")
    print(texto_salida)
    return 0


if __name__ == "__main__":
    sys.exit(main())
