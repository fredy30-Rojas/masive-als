#!/usr/bin/env python3
"""Barre el PDB buscando ligandos pequenos co-cristalizados con las dianas del proyecto.

Por que existe: el cuello de botella del rescoring no es el computo, es que faltan
positivos de union medida. Un ligando co-cristalizado es evidencia estructural de
union, con cita (el PDB y su articulo) y con quimiotipo propio: es la fuente mas
barata de positivos nuevos, y es sistematica (no depende de que uno se acuerde).

Criterio de entrada, el del proyecto:
  - el ligando es no polimerico y tiene SMILES depositado,
  - no es aditivo de cristalizacion (iones, sulfato, glicerol, PEG, DMSO, crioprotector),
  - tiene al menos 6 atomos pesados.

Lo que NO dice: que el ligando sea afin. Una co-cristalizacion es evidencia de union,
no una Kd. El sitio de cada ligando (Trp32 frente a sitio metalico, por ejemplo) hay
que leerlo de la propia estructura antes de usarlo como positivo; por eso se guarda
el PDB de cada uno.

Fuentes (todas de solo lectura):
  - RCSB search API    -> entradas del PDB por accession UniProt
  - RCSB GraphQL       -> entidades no polimericas de cada entrada
  - RCSB GraphQL       -> descriptor del componente quimico (SMILES, formula)

Salidas:
  analysis/ligandos_cristal/{sod1,tdp43,fus}.csv
  analysis/ligandos_cristal/resumen.txt

Uso:
  python analysis/buscar_ligandos_cristal.py
"""
from __future__ import annotations

import csv
import json
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

SALIDA = Path(__file__).resolve().parent / "ligandos_cristal"
BIOS = Path(__file__).resolve().parent / "verdad_de_referencia.csv"

# Diana -> accession UniProt. Las tres dianas del proyecto.
DIANAS = {
    "sod1": ("P00441", "SOD1 humana"),
    "tdp43": ("Q13148", "TDP-43 humana"),
    "fus": ("P35637", "FUS humana"),
}

SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
GQL = "https://data.rcsb.org/graphql"
LOTE = 50

# Aditivos de cristalizacion y iones: un ligando de estos no es un positivo.
ADITIVOS = {
    "HOH", "DOD", "SO4", "PO4", "ACT", "ACY", "EDO", "PEG", "PG4", "PGE", "1PE",
    "P6G", "GOL", "MPD", "DMS", "DMF", "TRS", "MES", "EPE", "IMD", "BME", "DTT",
    "CIT", "FMT", "NO3", "NCO", "NH4", "K", "NA", "CL", "BR", "IOD", "MG", "MN",
    "CA", "ZN", "CU", "CU1", "FE", "FE2", "NI", "CO", "CD", "HG", "SF4", "FES",
    "TAM", "BOG", "LDA", "OCT", "OLC", "HEX", "DIO", "C8E", "XPE", "IPA", "MOH",
    "EOH", "ACE", "SIN", "MLI", "MLT", "MLA", "TLA", "ETA", "PE4", "PEF", "MYR",
    "PLM", "STE", "FLC", "TAR", "FUM", "SUC", "MAL", "LAC", "PYR", "CO3", "BCT",
    "GLC", "NAG", "FUC", "MAN", "BMA", "TRE", "SER", "GLY", "ALA", "PGO", "BEN",
    "PEO", "CRY", "TAM", "PCA", "UNX", "UNL", "UNK", "2PE", "12P", "15P", "AZI",
    # Anadidos tras la primera pasada del 25 sep: aditivos que se colaron
    "TFA", "DSN", "CYS", "PE8", "TOE", "A1LYN", "UDI", "PE3", "PE4", "PG5",
    "PG6", "1PG", "3PG", "P33", "PEG4", "PG4", "XE", "XE2",
}
PALABRAS_ADITIVO = (
    "ION", "SULFATE", "PHOSPHATE", "ACETATE", "GLYCEROL", "ETHYLENE GLYCOL",
    "DIMETHYL SULFOXIDE", "WATER", "CHLORIDE", "BROMIDE", "SODIUM", "POTASSIUM",
    "CALCIUM", "MAGNESIUM", "ZINC", "COPPER", "NICKEL", "CADMIUM", "MANGANESE",
    "CRYOPROTECTANT", "BUFFER", "PEG", "NITRATE", "FORMATE", "CITRATE",
    # Anadidos tras la primera pasada del 25 sep
    "OXATRICOSANE", "ETHOXY", "TRIETHYLENE", "POLYETHYLENE", "TRIFLUOROACET",
    "DITHIOL", "TRIOL", "D-AMINO", "D-SERINE",
)


def pedir(url: str, datos: dict | None = None, intentos: int = 3):
    for i in range(intentos):
        try:
            if datos is None:
                with urllib.request.urlopen(url, timeout=60) as r:
                    return json.loads(r.read().decode("utf-8"))
            req = urllib.request.Request(
                url, data=json.dumps(datos).encode("utf-8"),
                headers={"Content-Type": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=90) as r:
                return json.loads(r.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            if i == intentos - 1:
                raise
            print(f"    reintento {i + 1} ({e})", flush=True)
            time.sleep(4)


def entradas(accession: str) -> list[str]:
    """Entradas del PDB cuya secuencia de proteina es esta accession."""
    consulta = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": (
                    "rcsb_polymer_entity_container_identifiers."
                    "reference_sequence_identifiers.database_accession"
                ),
                "operator": "exact_match",
                "value": accession,
            },
        },
        "return_type": "entry",
        "request_options": {"return_all_hits": True},
    }
    d = pedir(SEARCH, consulta)
    ids = [r["identifier"] for r in d.get("result_set", [])]
    print(f"  {accession}: {d.get('total_count')} entradas", flush=True)
    return sorted(ids)


def ligandos_de(ids: list[str]):
    """Compuestos no polimericos de cada entrada, por lotes de GraphQL."""
    por_entrada: dict[str, list[dict]] = {}
    titulos: dict[str, str] = {}
    campos = "rcsb_id struct { title } nonpolymer_entities { pdbx_entity_nonpoly { comp_id name } }"
    for i in range(0, len(ids), LOTE):
        lote = ids[i : i + LOTE]
        q = "{ entries(entry_ids: %s) { %s } }" % (json.dumps(lote), campos)
        d = pedir(GQL, {"query": q})
        if "errors" in d:  # campo no disponible: se repite sin el titulo
            q = "{ entries(entry_ids: %s) { rcsb_id nonpolymer_entities { pdbx_entity_nonpoly { comp_id name } } } }" % json.dumps(lote)
            d = pedir(GQL, {"query": q})
        for e in d.get("data", {}).get("entries", []) or []:
            nid = e["rcsb_id"]
            titulos[nid] = ((e.get("struct") or {}).get("title") or "")[:120]
            por_entrada[nid] = [
                x["pdbx_entity_nonpoly"]
                for x in (e.get("nonpolymer_entities") or [])
                if x and x.get("pdbx_entity_nonpoly")
            ]
        print(f"  entradas {i + len(lote)}/{len(ids)}", flush=True)
    return por_entrada, titulos


def descriptores(comps: list[str]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for i in range(0, len(comps), LOTE):
        lote = comps[i : i + LOTE]
        q = (
            "{ chem_comps(comp_ids: %s) { chem_comp { id name formula formula_weight type } "
            "rcsb_chem_comp_descriptor { SMILES } } }" % json.dumps(lote)
        )
        d = pedir(GQL, {"query": q})
        for c in d.get("data", {}).get("chem_comps", []) or []:
            base = c.get("chem_comp") or {}
            desc = c.get("rcsb_chem_comp_descriptor") or {}
            out[base.get("id", "")] = {
                "nombre": base.get("name") or "",
                "formula": base.get("formula") or "",
                "peso": base.get("formula_weight") or "",
                "tipo": base.get("type") or "",
                "smiles": (desc.get("SMILES") or "").replace("\\", ""),
            }
        print(f"  compuestos {i + len(lote)}/{len(comps)}", flush=True)
    return out


def es_aditivo(comp: str, nombre: str) -> bool:
    if comp.upper() in ADITIVOS:
        return True
    mayus = (nombre or "").upper()
    return any(p in mayus for p in PALABRAS_ADITIVO)


def esqueleto(smiles: str) -> tuple[str, int]:
    """Nucleo de Murcko y numero de atomos pesados. Vacio si el SMILES no parsea."""
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold

    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return "", 0
    return MurckoScaffold.MurckoScaffoldSmiles(mol=m), m.GetNumHeavyAtoms()


def canonico(smiles: str) -> str:
    from rdkit import Chem

    m = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(m) if m else ""


def ya_en_verdad(objetivo: str) -> dict[str, str]:
    """SMILES canonicos SIN estereoquimica de los ligandos ya en la verdad de referencia.

    Sin estereoquimica porque el descriptor del PDB no lleva la del ligando depositado
    (el isoproterenol de 4A7T, por ejemplo, aparece sin ella): comparando con estereo el
    mismo compuesto no se reconoceria y se contaria dos veces como positivo nuevo.
    """
    if not BIOS.exists():
        return {}
    out: dict[str, str] = {}
    with open(BIOS, encoding="utf-8") as fh:
        for fila in csv.DictReader(fh):
            if not fila.get("smiles"):
                continue
            from rdkit import Chem

            m = Chem.MolFromSmiles(fila["smiles"])
            if m is None:
                continue
            c = Chem.MolToSmiles(m, isomericSmiles=False)
            out[c] = fila.get("ligand", "")
    return out


def barre(clave: str, accession: str, descripcion: str, conocidos: dict[str, str]):
    print(f"\n{descripcion} ({accession})")
    ids = entradas(accession)
    if not ids:
        return []
    por_entrada, titulos = ligandos_de(ids)

    usos: dict[str, list[str]] = {}
    nombres: dict[str, str] = {}
    for pid, ligs in por_entrada.items():
        for lig in ligs:
            comp, nombre = lig.get("comp_id", ""), lig.get("name", "")
            if not comp or es_aditivo(comp, nombre):
                continue
            usos.setdefault(comp, []).append(pid)
            nombres[comp] = nombre

    comps = sorted(usos)
    desc = descriptores(comps) if comps else {}

    filas = []
    for comp in comps:
        info = desc.get(comp, {})
        smi = info.get("smiles", "")
        esc, pesados = esqueleto(smi) if smi else ("", 0)
        if not smi or pesados < 6:
            continue
        pids = sorted(usos[comp])
        from rdkit import Chem

        _m = Chem.MolFromSmiles(smi)
        plano = Chem.MolToSmiles(_m, isomericSmiles=False) if _m else ""
        filas.append(
            {
                "diana": clave,
                "comp_id": comp,
                "nombre": info.get("nombre") or nombres.get(comp, ""),
                "formula": info.get("formula", ""),
                "pesados": pesados,
                "smiles": smi,
                "esqueleto": esc,
                "n_estructuras": len(pids),
                "pdbs": ";".join(pids[:10]) + ("..." if len(pids) > 10 else ""),
                "titulo_pdb": titulos.get(pids[0], ""),
                "ya_en_verdad": conocidos.get(plano, "no"),
            }
        )
    return filas


def main() -> int:
    SALIDA.mkdir(parents=True, exist_ok=True)
    conocidos = ya_en_verdad(str(BIOS))
    print(f"Ligandos ya en la verdad de referencia: {len(conocidos)}")

    todo, resumen = [], []
    for clave, (accession, descripcion) in DIANAS.items():
        filas = barre(clave, accession, descripcion, conocidos)
        todo += filas
        with open(SALIDA / f"{clave}.csv", "w", newline="", encoding="utf-8") as fh:
            if filas:
                w = csv.DictWriter(fh, fieldnames=list(filas[0].keys()))
                w.writeheader()
                w.writerows(filas)

        nuevos = [f for f in filas if f["ya_en_verdad"] == "no"]
        esqs = {}
        for f in nuevos:
            esqs.setdefault(f["esqueleto"], []).append(f)
        resumen.append(
            f"\n### {descripcion} ({accession})\n"
            f"Ligandos utiles (no aditivos, >=6 atomos pesados): {len(filas)}\n"
            f"Ya en la verdad de referencia: {len(filas) - len(nuevos)}\n"
            f"Nuevos: {len(nuevos)}   Quimias nuevas (Murcko): {len(esqs)}\n"
        )
        for esc, fs in sorted(esqs.items(), key=lambda kv: -len(kv[1]))[:25]:
            fs.sort(key=lambda f: -f["pesados"])
            ej = fs[0]
            resumen.append(
                f"  {len(fs):>3} moleculas  {ej['pesados']:>3} pesados  "
                f"{ej['comp_id']:<4} {ej['nombre'][:38]:<38} {ej['pdbs'][:34]:<34} {ej['esqueleto'][:40]}"
            )
        print(f"  {descripcion}: {len(filas)} ligandos utiles, {len(nuevos)} nuevos, {len(esqs)} quimias nuevas")

    with open(SALIDA / "ligandos_cristal.csv", "w", newline="", encoding="utf-8") as fh:
        if todo:
            w = csv.DictWriter(fh, fieldnames=list(todo[0].keys()))
            w.writeheader()
            w.writerows(todo)

    texto = (
        "LIGANDOS CO-CRISTALIZADOS CON LAS DIANAS DEL PROYECTO (barrido sistematico del PDB)\n"
        "Criterio: no polimerico, con SMILES, no aditivo de cristalizacion, >=6 atomos pesados.\n"
        "Una co-cristalizacion es evidencia de union, NO una Kd: el sitio hay que leerlo\n"
        "de la estructura antes de usarlo como positivo.\n"
        + "".join(resumen)
    )
    (SALIDA / "resumen.txt").write_text(texto, encoding="utf-8")
    print("\n" + texto)
    return 0


if __name__ == "__main__":
    sys.exit(main())
