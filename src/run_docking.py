#!/usr/bin/env python3
"""
Pipeline de docking molecular REAL con AutoDock Vina.
- Carga proteinas reales (AlphaFold/PDB)
- Convierte compuestos a PDBQT
- Ejecuta Vina para cada proteina x ligando
- Genera resultados y ranking de hits
"""
import os, sys, csv, subprocess, glob, json
from pathlib import Path

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PROTEINS_DIR = os.path.join(PROJECT_ROOT, "proteins")
COMPOUNDS_DIR = os.path.join(PROJECT_ROOT, "compounds")
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "local_dock")
TOOLS_DIR = os.path.join(PROJECT_ROOT, "tools")
VINA_EXE = os.path.join(TOOLS_DIR, "vina.exe")

# Encontrar obabel.exe (OpenBabel CLI)
def _find_obabel():
    candidates = [
        "obabel",
        os.path.join(os.path.dirname(sys.executable), "Scripts", "obabel.exe"),
        os.path.join(os.path.dirname(sys.executable), "Lib", "site-packages", "openbabel", "bin", "obabel.exe"),
    ]
    for c in candidates:
        if os.path.exists(c):
            return c
    import shutil as _shutil
    found = _shutil.which("obabel")
    return found if found else "obabel"

OBABEL = _find_obabel()

TARGETS = ["TDP43", "SOD1", "FUS"]

# Coordenadas de los sitios activos (centro y tamaño de la caja de docking)
# Determinadas desde literatura cientifica
BINDING_SITES = {
    "TDP43": {"center": (12.0, 8.0, 15.0), "size": (22, 22, 22)},  # RRM1-RRM2 interface
    "SOD1":  {"center": (2.0, 15.0, -4.0), "size": (20, 20, 20)},  # Dimer interface
    "FUS":   {"center": (5.0, 3.0, 8.0),  "size": (22, 22, 22)},  # RRM domain
}


def prepare_protein(target_name):
    """Prepara la proteina para docking: limpia y convierte a PDBQT."""
    protein_dir = os.path.join(PROTEINS_DIR, target_name)
    out_pdbqt = os.path.join(RESULTS_DIR, f"{target_name}.pdbqt")
    
    os.makedirs(RESULTS_DIR, exist_ok=True)
    
    # Buscar PDB/AF disponible
    pdb_files = list(Path(protein_dir).glob("*.pdb"))
    if not pdb_files:
        print(f"  [{target_name}] ERROR: No se encontro archivo PDB en {protein_dir}")
        return None
    
    pdb_file = str(pdb_files[0])
    print(f"  [{target_name}] Usando: {os.path.basename(pdb_file)}")
    
    # Si ya existe PDBQT, reusar
    if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 100:
        print(f"  [{target_name}] PDBQT ya existe")
        return out_pdbqt
    
    # Preparar con OpenBabel: convertir a PDBQT (sin gen3d para evitar crash)
    result = subprocess.run(
        [OBABEL, "-ipdb", pdb_file, "-opdbqt", "-O", out_pdbqt, "-xr"],
        capture_output=True, text=True, timeout=30
    )
    
    if result.returncode != 0:
        print(f"  [{target_name}] ERROR OpenBabel: {result.stderr[:200]}")
        # Intentar sin quitar residuos
        result2 = subprocess.run(
            [OBABEL, "-ipdb", pdb_file, "-opdbqt", "-O", out_pdbqt],
            capture_output=True, text=True, timeout=30
        )
        if result2.returncode != 0:
            return None
    
    size_kb = os.path.getsize(out_pdbqt) / 1024
    print(f"  [{target_name}] PDBQT listo: {size_kb:.0f} KB")
    return out_pdbqt


def prepare_single_ligand_3d(smiles, output_pdbqt):
    """Genera estructura 3D con OpenBabel Python API y guarda PDBQT."""
    try:
        from openbabel import openbabel as ob
        
        conv = ob.OBConversion()
        conv.SetInFormat("smi")
        mol = ob.OBMol()
        if not conv.ReadString(mol, smiles):
            return False
        
        # Generar 3D
        mol.AddHydrogens()
        builder = ob.OBBuilder()
        builder.Build(mol)
        
        # Optimizar con UFF (mas simple que MMFF94, no necesita parametros)
        ff = ob.OBForceField.FindForceField("UFF")
        if ff:
            ff.Setup(mol)
            ff.SteepestDescent(200)
            ff.ConjugateGradients(100)
            ff.GetCoordinates(mol)
        
        # Escribir PDBQT
        conv.SetOutFormat("pdbqt")
        conv.WriteFile(mol, output_pdbqt)
        return True
    except Exception as e:
        print(f"    ERROR 3D: {e}")
        return False


def prepare_ligands(target_name):
    """Convierte la libreria de compuestos a archivos PDBQT individuales con 3D."""
    ligands_dir = os.path.join(RESULTS_DIR, f"{target_name}_ligands")
    os.makedirs(ligands_dir, exist_ok=True)
    
    csv_file = os.path.join(COMPOUNDS_DIR, "fda_subset.csv")
    
    # Si ya tenemos ligandos preparados, devolver lista
    existing = list(Path(ligands_dir).glob("*.pdbqt"))
    if len(existing) > 5:
        print(f"  [{target_name}] {len(existing)} ligandos PDBQT ya preparados")
        return existing
    
    # Leer SMILES del CSV y generar 3D con OpenBabel Python API
    if not os.path.exists(csv_file):
        print(f"  [{target_name}] ERROR: No hay compuestos descargados")
        return []
    
    print(f"  [{target_name}] Generando estructuras 3D con OpenBabel Python API...")
    
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        rows = list(reader)
    
    ligand_paths = []
    for i, row in enumerate(rows):
        smiles = row.get("smiles", "")
        if not smiles:
            continue
        
        out_path = os.path.join(ligands_dir, f"ligand{i+1}.pdbqt")
        
        if os.path.exists(out_path) and os.path.getsize(out_path) > 200:
            ligand_paths.append(out_path)
            continue
        
        if prepare_single_ligand_3d(smiles, out_path):
            ligand_paths.append(out_path)
    
    print(f"  [{target_name}] {len(ligand_paths)} ligandos 3D generados")
    return ligand_paths


def run_docking(target_name, protein_pdbqt, ligand_paths):
    """Ejecuta AutoDock Vina para cada ligando contra la proteina."""
    site = BINDING_SITES[target_name]
    cx, cy, cz = site["center"]
    sx, sy, sz = site["size"]
    
    dock_results = []
    output_dir = os.path.join(RESULTS_DIR, f"{target_name}_outputs")
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"  [{target_name}] Ejecutando docking: {len(ligand_paths)} ligandos...")
    print(f"  [{target_name}] Sitio activo: center=({cx},{cy},{cz}) size=({sx},{sy},{sz})")
    
    for i, lig_path in enumerate(ligand_paths):
        lig_name = os.path.basename(lig_path).replace(".pdbqt", "")
        out_pdbqt = os.path.join(output_dir, f"{lig_name}_docked.pdbqt")
        log_file = os.path.join(output_dir, f"{lig_name}_log.txt")
        
        # Saltar si ya procesado
        if os.path.exists(out_pdbqt) and os.path.getsize(out_pdbqt) > 100:
            # Leer energia del log
            energy = None
            if os.path.exists(log_file):
                with open(log_file) as f:
                    for line in f:
                        if "RESULT:" in line or "Affinity:" in line:
                            try:
                                energy = float(line.split(":")[-1].strip().split()[0])
                            except:
                                pass
            if energy:
                dock_results.append({
                    "ligand": lig_name.replace("ligand", "Compuesto"),
                    "target": target_name,
                    "energy": energy,
                    "file": out_pdbqt
                })
            continue
        
        cmd = [
            VINA_EXE,
            "--receptor", protein_pdbqt,
            "--ligand", lig_path,
            "--out", out_pdbqt,
            "--center_x", str(cx),
            "--center_y", str(cy),
            "--center_z", str(cz),
            "--size_x", str(sx),
            "--size_y", str(sy),
            "--size_z", str(sz),
            "--exhaustiveness", "2",
            "--num_modes", "5",
        ]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
            
            # Extraer energia de la salida - buscar en la tabla de resultados
            energy = None
            in_table = False
            for line in result.stdout.split("\n"):
                stripped = line.strip()
                # Detectar inicio de la tabla de resultados
                if "mode" in stripped and "affinity" in stripped:
                    in_table = True
                    continue
                if not in_table or not stripped:
                    continue
                if stripped.startswith("-----"):
                    continue
                # Parsear lineas de modo: "1       -5.2      0.0      0.0"
                parts = stripped.split()
                if len(parts) >= 2 and parts[0].isdigit():
                    try:
                        mode_num = int(parts[0])
                        e = float(parts[1])
                        # Modo 1 puede ser 0 (pose de entrada), buscar modos 2+
                        if mode_num >= 2:
                            energy = e
                            break
                        # Si solo hay modo 1, usar ese
                        if mode_num == 1 and energy is None:
                            energy = e
                    except (ValueError, IndexError):
                        pass
            
            if energy is None:
                # Leer del log
                if os.path.exists(log_file):
                    with open(log_file) as f:
                        for line in f:
                            if "RESULT:" in line or "Affinity:" in line:
                                try:
                                    energy = float(line.split(":")[-1].strip().split()[0])
                                except:
                                    pass
            
            if energy:
                dock_results.append({
                    "ligand": lig_name.replace("ligand", "Compuesto"),
                    "target": target_name,
                    "energy": energy,
                    "file": out_pdbqt
                })
                
            if (i + 1) % 5 == 0 or i == len(ligand_paths) - 1:
                print(f"  [{target_name}] Progreso: {i+1}/{len(ligand_paths)}")
                
        except subprocess.TimeoutExpired:
            print(f"  [{target_name}] TIMEOUT: {lig_name}")
        except Exception as e:
            print(f"  [{target_name}] ERROR docking {lig_name}: {e}")
    
    return dock_results


def main():
    print("=" * 60)
    print(" DOCKING MOLECULAR REAL - MASIVE-ALS")
    print("=" * 60)
    
    if not os.path.exists(VINA_EXE):
        print(f"ERROR: AutoDock Vina no encontrado en {VINA_EXE}")
        sys.exit(1)
    
    all_results = []
    
    for target in TARGETS:
        print(f"\n{'='*40}")
        print(f" PROCESANDO: {target}")
        print(f"{'='*40}")
        
        # Preparar proteina
        protein_pdbqt = prepare_protein(target)
        if not protein_pdbqt:
            print(f"  [{target}] ERROR: No se pudo preparar la proteina. Saltando.")
            continue
        
        # Preparar ligandos
        ligands = prepare_ligands(target)
        if not ligands:
            print(f"  [{target}] ERROR: No hay ligandos. Saltando.")
            continue
        
        # Ejecutar docking
        results = run_docking(target, protein_pdbqt, ligands)
        all_results.extend(results)
        
        # Top 5 de este target
        results.sort(key=lambda x: x["energy"])
        print(f"\n  [{target}] TOP 5 HITS:")
        for r in results[:5]:
            print(f"    {r['ligand']:20s} | Energia: {r['energy']:8.2f} kcal/mol")
    
    # Guardar resultados completos
    all_results.sort(key=lambda x: x["energy"])
    csv_path = os.path.join(RESULTS_DIR, "docking_results.csv")
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["ligand", "target", "energy", "file"])
        writer.writeheader()
        writer.writerows(all_results)
    
    # Mostrar resumen
    print(f"\n{'='*60}")
    print(" RESUMEN FINAL DE DOCKING")
    print(f"{'='*60}")
    print(f"  Total dockings: {len(all_results)}")
    for target in TARGETS:
        t_results = [r for r in all_results if r["target"] == target]
        if t_results:
            best = min(t_results, key=lambda x: x["energy"])
            print(f"  {target:6s}: {len(t_results):3d} docks | Mejor: {best['energy']:7.2f} kcal/mol ({best['ligand']})")
    
    print(f"\n  Resultados guardados en: {csv_path}")
    return all_results


if __name__ == "__main__":
    main()
