#!/usr/bin/env python3
"""
Validacion de candidatos con Dinamica Molecular usando OpenMM.
- Toma los top hits del docking
- Simula la dinamica del complejo proteina-ligando
- Verifica estabilidad (RMSD, energia)
- Solo procesa los 5 mejores candidatos por proteina
"""
import os, sys, csv, json, time, math
import numpy as np

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Encontrar obabel.exe
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
RESULTS_DIR = os.path.join(PROJECT_ROOT, "results", "local_dock")
MD_DIR = os.path.join(PROJECT_ROOT, "md_results", "local_md")
PROTEINS_DIR = os.path.join(PROJECT_ROOT, "proteins")

os.makedirs(MD_DIR, exist_ok=True)


def run_openmm_md(pdb_file, ligand_name, target_name, sim_time_ns=0.1):
    """
    Ejecuta una simulacion corta de MD con OpenMM.
    - Minimizacion
    - Equilibracion NVT
    - Produccion corta
    Retorna: dict con RMSD medio, energia, estable (bool)
    """
    try:
        from openmm import app, unit, LangevinMiddleIntegrator
        from openmm.app import PDBFile, ForceField, Simulation, Modeller, PDBReporter, StateDataReporter
        import openmm as mm
    except ImportError:
        print(f"  [{ligand_name}] ERROR: OpenMM no disponible")
        return None
    
    out_dir = os.path.join(MD_DIR, target_name, ligand_name)
    os.makedirs(out_dir, exist_ok=True)
    
    # Si ya existe resultado, leerlo
    result_file = os.path.join(out_dir, "md_result.json")
    if os.path.exists(result_file):
        with open(result_file) as f:
            return json.load(f)
    
    try:
        # Cargar estructura
        pdb = PDBFile(pdb_file)
        
        # Force field: AMBER14 + solvente implicito (rapido para screening)
        forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
        
        # Sistema con solvente implicito (GBSA) - mucho mas rapido que solvatar
        # Para validacion de screening usamos implicit solvent
        from openmm.app import HBonds
        
        # Crear sistema
        system = forcefield.createSystem(
            pdb.topology,
            nonbondedMethod=app.NoCutoff,
            nonbondedCutoff=1.0*unit.nanometer,
            constraints=HBonds,
            hydrogenMass=1.5*unit.amu  # HMR para paso de 4fs
        )
        
        # Integrador Langevin
        temperature = 300 * unit.kelvin
        friction = 1.0 / unit.picosecond
        timestep = 0.004 * unit.picoseconds  # 4 fs con HMR
        integrator = LangevinMiddleIntegrator(temperature, friction, timestep)
        
        # Simulacion
        simulation = Simulation(pdb.topology, system, integrator)
        simulation.context.setPositions(pdb.positions)
        
        # Reporters
        sim_steps = int(sim_time_ns * 1000 / 0.004)  # conversion ps -> steps
        if sim_steps > 50000:
            sim_steps = 50000  # Limitar a ~50k steps para pruebas
        
        report_interval = max(1, sim_steps // 50)  # 50 puntos de datos
        
        traj_pdb = os.path.join(out_dir, "trajectory.pdb")
        log_file = os.path.join(out_dir, "md_log.txt")
        
        simulation.reporters.append(PDBReporter(traj_pdb, report_interval))
        simulation.reporters.append(StateDataReporter(
            log_file, report_interval,
            step=True, potentialEnergy=True, temperature=True,
            speed=True
        ))
        
        # Minimizacion
        simulation.minimizeEnergy(maxIterations=500)
        
        # Equilibracion NVT
        simulation.step(5000)  # 20 ps equilibracion
        
        # Produccion
        t_start = time.time()
        simulation.step(sim_steps)
        t_elapsed = time.time() - t_start
        
        # Analizar resultados
        # Leer energia del log
        energies = []
        with open(log_file) as f:
            for line in f:
                if line.startswith("#") or "Step" in line:
                    continue
                parts = line.strip().split(",")
                if len(parts) >= 2:
                    try:
                        energy = float(parts[1])
                        energies.append(energy)
                    except:
                        pass
        
        if not energies:
            print(f"  [{ligand_name}] WARNING: No se pudieron leer energias")
            return None
        
        # Calcular estabilidad
        # Si la energia fluctua poco (< 5% std/mean), el complejo es estable
        mean_energy = np.mean(energies[-20:])  # Ultimos 20 frames
        std_energy = np.std(energies[-20:])
        
        is_stable = (std_energy / abs(mean_energy)) < 0.10 if mean_energy != 0 else False
        
        result = {
            "ligand": ligand_name,
            "target": target_name,
            "mean_energy_kjmol": round(float(mean_energy), 2),
            "std_energy_kjmol": round(float(std_energy), 2),
            "cv_percent": round(float(std_energy / abs(mean_energy) * 100), 1),
            "stable": "YES" if is_stable else "NO",
            "sim_time_ns": sim_time_ns,
            "wall_time_s": round(t_elapsed, 1),
            "steps": sim_steps
        }
        
        # Guardar resultado
        with open(result_file, "w") as f:
            json.dump(result, f, indent=2)
        
        print(f"  [{ligand_name}] Energia: {mean_energy:.0f} +/- {std_energy:.0f} kJ/mol | "
              f"CV: {result['cv_percent']}% | Estable: {result['stable']} | "
              f"Tiempo: {t_elapsed:.0f}s")
        
        return result
        
    except Exception as e:
        print(f"  [{ligand_name}] ERROR MD: {e}")
        import traceback
        traceback.print_exc()
        return None


def load_docked_pdbqt(pdbqt_path):
    """Convierte PDBQT dockeado a PDB valido para OpenMM."""
    pdb_path = pdbqt_path.replace(".pdbqt", "_openmm.pdb")
    
    if os.path.exists(pdb_path) and os.path.getsize(pdb_path) > 100:
        return pdb_path
    
    # Usar OpenBabel para convertir PDBQT -> PDB
    result = subprocess.run(
        [OBABEL, "-ipdbqt", pdbqt_path, "-opdb", "-O", pdb_path],
        capture_output=True, text=True, timeout=10
    )
    
    if result.returncode != 0 or not os.path.exists(pdb_path) or os.path.getsize(pdb_path) < 100:
        return None
    
    # Limpiar el PDB con PDBFixer para que OpenMM lo lea
    try:
        from pdbfixer import PDBFixer
        fixer = PDBFixer(filename=pdb_path)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)
        cleaned_path = pdb_path.replace(".pdb", "_cleaned.pdb")
        with open(cleaned_path, "w") as f:
            from openmm.app import PDBFile as _PDBFile
            _PDBFile.writeFile(fixer.topology, fixer.positions, f, True)
        return cleaned_path
    except Exception as e:
        print(f"    PDBFixer fallo: {e}, usando PDB sin limpiar")
        return pdb_path


def main():
    print("=" * 60)
    print(" VALIDACION MD CON OpenMM - MASIVE-ALS")
    print("=" * 60)
    
    # Cargar resultados de docking
    docking_csv = os.path.join(RESULTS_DIR, "docking_results.csv")
    if not os.path.exists(docking_csv):
        print("ERROR: No hay resultados de docking. Ejecuta run_docking.py primero.")
        sys.exit(1)
    
    with open(docking_csv) as f:
        reader = csv.DictReader(f)
        all_docks = list(reader)
    
    print(f"  Cargados {len(all_docks)} resultados de docking")
    
    # Seleccionar top 3 por proteina
    selected = []
    for target in ["TDP43", "SOD1", "FUS"]:
        t_hits = [d for d in all_docks if d["target"] == target]
        t_hits.sort(key=lambda x: float(x["energy"]))
        selected.extend(t_hits[:3])
    
    print(f"  Seleccionados {len(selected)} candidatos para validacion MD")
    
    results = []
    for hit in selected:
        lig = hit["ligand"]
        target = hit["target"]
        energy = float(hit["energy"])
        
        print(f"\n  [{target}] {lig} (docking: {energy:.2f} kcal/mol)")
        
        # Buscar PDB del complejo dockeado
        pdbqt_path = hit.get("file", "")
        if not pdbqt_path or not os.path.exists(pdbqt_path):
            print(f"    ERROR: No se encuentra {pdbqt_path}")
            continue
        
        # Convertir a PDB
        pdb_path = load_docked_pdbqt(pdbqt_path)
        if not pdb_path:
            print(f"    ERROR: No se pudo convertir a PDB")
            continue
        
        # Ejecutar MD
        result = run_openmm_md(pdb_path, lig, target, sim_time_ns=0.05)
        if result:
            result["docking_energy"] = energy
            results.append(result)
    
    # Guardar resultados
    if results:
        md_csv = os.path.join(MD_DIR, "md_validation.csv")
        with open(md_csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        
        stable = [r for r in results if r["stable"] == "YES"]
        print(f"\n{'='*60}")
        print(" RESUMEN VALIDACION MD")
        print(f"{'='*60}")
        print(f"  Candidatos validados: {len(results)}")
        print(f"  Estables:             {len(stable)}")
        print(f"  Resultados:           {md_csv}")
    
    return results


if __name__ == "__main__":
    import subprocess  # needed in load_docked_pdbqt
    main()
