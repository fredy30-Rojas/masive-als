#!/usr/bin/env python3
"""
Descarga estructuras 3D reales de TDP-43, SOD1 y FUS desde:
- AlphaFold DB (https://alphafold.ebi.ac.uk/)
- RCSB PDB (https://www.rcsb.org/)

UniProt IDs:
- TDP-43: Q13148 (TARDBP_HUMAN)
- SOD1:  P00441 (SODC_HUMAN)
- FUS:   P35637 (FUS_HUMAN)
"""
import os, sys, urllib.request, json, time

PROTEINS_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "proteins")

PROTEINS = {
    "TDP43": {
        "uniprot": "Q13148",
        "gene": "TARDBP",
        "pdb_best": "6b1n",  # TDP-43 RRM domains (crystal structure)
        "description": "TAR DNA-binding protein 43"
    },
    "SOD1": {
        "uniprot": "P00441",
        "gene": "SOD1",
        "pdb_best": "1hl5",  # Human SOD1 (crystal structure)
        "description": "Superoxide dismutase [Cu-Zn]"
    },
    "FUS": {
        "uniprot": "P35637",
        "gene": "FUS",
        "pdb_best": "6g99",  # FUS RRM domain
        "description": "RNA-binding protein FUS"
    }
}

def download_alphafold(uniprot_id, target_name):
    """Descarga estructura AlphaFold desde EBI."""
    out_dir = os.path.join(PROTEINS_DIR, target_name)
    os.makedirs(out_dir, exist_ok=True)
    
    # AlphaFold DB URL
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    out_pdb = os.path.join(out_dir, f"AF-{target_name}.pdb")
    
    if os.path.exists(out_pdb) and os.path.getsize(out_pdb) > 1000:
        print(f"  [{target_name}] AlphaFold ya descargado: {out_pdb}")
        return out_pdb
    
    print(f"  [{target_name}] Descargando AlphaFold de {url}...")
    try:
        urllib.request.urlretrieve(url, out_pdb)
        size_kb = os.path.getsize(out_pdb) / 1024
        print(f"  [{target_name}] OK: {size_kb:.0f} KB")
        return out_pdb
    except Exception as e:
        print(f"  [{target_name}] ERROR AlphaFold: {e}")
        return None


def download_pdb(pdb_id, target_name):
    """Descarga estructura cristalográfica desde RCSB PDB."""
    out_dir = os.path.join(PROTEINS_DIR, target_name)
    os.makedirs(out_dir, exist_ok=True)
    
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    out_pdb = os.path.join(out_dir, f"PDB-{pdb_id}.pdb")
    
    if os.path.exists(out_pdb) and os.path.getsize(out_pdb) > 1000:
        print(f"  [{target_name}] PDB {pdb_id} ya descargado")
        return out_pdb
    
    print(f"  [{target_name}] Descargando PDB {pdb_id}...")
    try:
        urllib.request.urlretrieve(url, out_pdb)
        size_kb = os.path.getsize(out_pdb) / 1024
        print(f"  [{target_name}] OK: {size_kb:.0f} KB")
        return out_pdb
    except Exception as e:
        print(f"  [{target_name}] ERROR PDB: {e}")
        return None


def main():
    print("=" * 60)
    print(" DESCARGA DE ESTRUCTURAS PROTEICAS REALES")
    print("=" * 60)
    
    results = {}
    for name, info in PROTEINS.items():
        print(f"\n>>> {name} ({info['description']})")
        
        # Intentar AlphaFold primero
        af_path = download_alphafold(info["uniprot"], name)
        
        # Intentar PDB como respaldo
        pdb_path = download_pdb(info["pdb_best"], name)
        
        results[name] = {
            "alphafold": af_path,
            "pdb": pdb_path,
            "uniprot": info["uniprot"]
        }
        
        # Guardar metadata
        meta_file = os.path.join(PROTEINS_DIR, name, "metadata.json")
        with open(meta_file, "w") as f:
            json.dump(info, f, indent=2)
    
    print("\n" + "=" * 60)
    print(" RESUMEN DE DESCARGAS")
    print("=" * 60)
    for name, paths in results.items():
        af_ok = "SI" if paths["alphafold"] else "NO"
        pdb_ok = "SI" if paths["pdb"] else "NO"
        print(f"  {name:6s} | AlphaFold: {af_ok:3s} | PDB: {pdb_ok:3s}")
    
    return results


if __name__ == "__main__":
    main()
