#!/usr/bin/env python3
"""Busca positivos de SOD1 humana en ChEMBL con ensayo de union DIRECTA.

El criterio de entrada del proyecto es: sin cita y sin tipo de ensayo, el
compuesto no entra. Aqui se aplica tal cual:

  - diana: CHEMBL2354 (superoxido dismutasa [Cu-Zn] humana, UniProt P00441)
  - solo ensayos con assay_type = B (binding, union directa)
  - tipos de medida: Kd, Ki, IC50, EC50 (todas se guardan; el tipo se anota)
  - se descartan los valores sin unidades M y los no numericos

Salida:
  analysis/positivos_sod1_chembl/actividades_binding.csv  (crudo)
  analysis/positivos_sod1_chembl/por_esqueleto.csv        (agrupado por Murcko)
  analysis/positivos_sod1_chembl/resumen.txt              (lo legible)

Uso:
  python analysis/buscar_positivos_sod1.py
"""
from __future__ import annotations

import csv
import json
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

BASE = "https://www.ebi.ac.uk/chembl/api/data"
TARGET = "CHEMBL2354"  # SOD1 humana
DATA = Path(__file__).resolve().parent / "positivos_sod1_chembl"
TIPOS = {"Kd", "Ki", "IC50", "EC50", "KD", "KI"}
MAX_PAGINAS = 40


def get(url: str, intentos: int = 3):
    for i in range(intentos):
        try:
            with urllib.request.urlopen(url, timeout=45) as r:
                return json.loads(r.read().decode("utf-8"))
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as e:
            if i == intentos - 1:
                raise
            print(f"  reintento {i + 1} ({e})", flush=True)
            time.sleep(3)


def actividades():
    """Todas las actividades de union directa contra SOD1."""
    filas, offset = [], 0
    for _ in range(MAX_PAGINAS):
        url = (
            f"{BASE}/activity.json?target_chembl_id={TARGET}"
            f"&assay_type=B&limit=1000&offset={offset}"
        )
        d = get(url)
        acts = d.get("activities", [])
        for a in acts:
            tipo = a.get("standard_type")
            val = a.get("standard_value")
            uni = a.get("standard_units")
            if not tipo or tipo.upper() not in TIPOS:
                continue
            if val is None or uni not in ("nM", "uM", "pM", "M"):
                continue
            filas.append(
                {
                    "molecula": a.get("molecule_chembl_id"),
                    "nombre": (a.get("molecule_pref_name") or ""),
                    "tipo": tipo,
                    "valor": val,
                    "unidades": uni,
                    "pchembl": a.get("pchembl_value") or "",
                    "ensayo": a.get("assay_chembl_id"),
                    "documento": a.get("document_chembl_id"),
                    "ano": a.get("document_year") or "",
                    "descripcion_ensayo": (a.get("assay_description") or "")[:160],
                }
            )
        total = d.get("page_meta", {}).get("total_count")
        print(f"  offset {offset}: {len(acts)} actividades (total {total})", flush=True)
        if not d.get("page_meta", {}).get("next"):
            break
        offset += 1000
    return filas


def smiles(lote, tam=50):
    """SMILES canonico y peso de cada molecula, por lotes."""
    out = {}
    for i in range(0, len(lote), tam):
        ids = ",".join(lote[i : i + tam])
        d = get(f"{BASE}/molecule.json?molecule_chembl_id__in={ids}&limit={tam}")
        for m in d.get("molecules", []):
            props = m.get("molecule_properties") or {}
            # El SMILES vive en molecule_structures, no en la raiz del objeto
            struct = m.get("molecule_structures") or {}
            out[m["molecule_chembl_id"]] = {
                "smiles": struct.get("canonical_smiles") or m.get("canonical_smiles") or "",
                "peso": props.get("full_mwt") or "",
                "pesados": props.get("heavy_atoms") or "",
            }
        print(f"  smiles {i + len(lote[i:i+tam])}/{len(lote)}", flush=True)
    return out


def esqueletos(filas, smi):
    """Nucleo de Murcko con RDKit para agrupar quimias."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError:
        print("RDKit no disponible: se agrupa por SMILES completo.")
        return {f["molecula"]: smi.get(f["molecula"], {}).get("smiles", "") for f in filas}

    cache = {}
    for f in filas:
        mid = f["molecula"]
        if mid in cache:
            continue
        s = smi.get(mid, {}).get("smiles", "")
        if not s:
            cache[mid] = ""
            continue
        m = Chem.MolFromSmiles(s)
        cache[mid] = MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m else ""
    return cache


def main():
    DATA.mkdir(parents=True, exist_ok=True)
    print(f"Actividades de union directa contra {TARGET} (SOD1 humana)...")
    filas = actividades()
    print(f"Total actividades utiles: {len(filas)}")

    # Un registro por molecula: la mejor medida de cada tipo
    por_mol = {}
    for f in filas:
        mid = f["molecula"]
        if not mid:
            continue
        prev = por_mol.get(mid)
        if prev is None or (f["valor"] and prev["valor"] and f["valor"] < prev["valor"]):
            por_mol[mid] = f
    print(f"Moleculas distintas: {len(por_mol)}")

    ids = sorted(por_mol)
    info = smiles(ids) if ids else {}
    murcko = esqueletos([por_mol[i] for i in ids], info)

    with open(DATA / "actividades_binding.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(filas[0].keys()) if filas else ["molecula"])
        w.writeheader()
        w.writerows(filas)

    grupos = defaultdict(list)
    for mid, f in por_mol.items():
        grupos[murcko.get(mid, "")].append(
            {
                "molecula": mid,
                "nombre": f["nombre"],
                "tipo": f["tipo"],
                "valor_nM": round(float(f["valor"]) * {"nM": 1, "uM": 1e3, "pM": 1e-3, "M": 1e9}[f["unidades"]], 1),
                "pchembl": f["pchembl"],
                "ano": f["ano"],
                "documento": f["documento"],
                "smiles": info.get(mid, {}).get("smiles", ""),
                "peso": info.get(mid, {}).get("peso", ""),
                "pesados": info.get(mid, {}).get("pesados", ""),
                "ensayo": f["descripcion_ensayo"],
            }
        )

    filas_g = []
    for esc, miembros in sorted(grupos.items(), key=lambda kv: -len(kv[1])):
        miembros.sort(key=lambda m: m["valor_nM"])
        mejor = miembros[0]
        filas_g.append(
            {
                "esqueleto": esc,
                "n_compuestos": len(miembros),
                "mejor_valor_nM": mejor["valor_nM"],
                "mejor_molecula": mejor["molecula"],
                "mejor_nombre": mejor["nombre"],
                "tipo_medida": mejor["tipo"],
                "ano": mejor["ano"],
                "documento": mejor["documento"],
                "smiles_representante": mejor["smiles"],
                "pesados": mejor["pesados"],
                "todos": "; ".join(f"{m['molecula']}({m['valor_nM']}nM)" for m in miembros[:12]),
            }
        )
    with open(DATA / "por_esqueleto.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(filas_g[0].keys()) if filas_g else ["esqueleto"])
        w.writeheader()
        w.writerows(filas_g)

    lineas = [
        "POSITIVOS DE SOD1 CON ENSAYO DE UNION DIRECTA (ChEMBL, CHEMBL2354)",
        f"Actividades: {len(filas)}   Moleculas: {len(por_mol)}   Quimias (Murcko): {len(grupos)}",
        "",
        f"{'n':>4}  {'mejor_nM':>10}  {'tipo':<5} {'año':<5} {'molecula':<16} {'pesados':>7}  nombre",
    ]
    for g in filas_g[:60]:
        lineas.append(
            f"{g['n_compuestos']:>4}  {g['mejor_valor_nM']:>10}  {g['tipo_medida']:<5} "
            f"{g['ano']:<5} {g['mejor_molecula']:<16} {str(g['pesados']):>7}  {g['mejor_nombre']}"
        )
    texto = "\n".join(lineas) + "\n"
    (DATA / "resumen.txt").write_text(texto, encoding="utf-8")
    print(texto)


if __name__ == "__main__":
    sys.exit(main())
