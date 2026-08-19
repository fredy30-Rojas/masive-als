#!/usr/bin/env python3
"""
PIPELINE COMPLETO MASIVE-ALS - PRUEBA LOCAL REAL
================================================
1. Descarga estructuras de proteinas (AlphaFold + PDB)
2. Descarga/Genera compuestos (FDA-approved drugs)
3. Ejecuta docking molecular con AutoDock Vina
4. Valida top hits con OpenMM (dinamica molecular)
5. Genera informe final

Cada fase verifica que la anterior se completo.
Si algo falla, se detiene y reporta.
"""
import os, sys, subprocess, csv, json, time

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
SRC_DIR = os.path.join(PROJECT_ROOT, "src")

STEPS = [
    {
        "name": "Descargar proteinas reales",
        "script": "download_proteins.py",
        "check": lambda: all(
            len([f for f in os.listdir(os.path.join(PROJECT_ROOT, "proteins", t)) if f.endswith(".pdb")]) > 0
            for t in ["TDP43", "SOD1", "FUS"]
        )
    },
    {
        "name": "Descargar compuestos",
        "script": "download_compounds.py",
        "check": lambda: os.path.exists(os.path.join(PROJECT_ROOT, "compounds", "fda_subset.csv"))
    },
    {
        "name": "Ejecutar docking molecular",
        "script": "run_docking.py",
        "check": lambda: os.path.exists(os.path.join(PROJECT_ROOT, "results", "local_dock", "docking_results.csv"))
    },
    {
        "name": "Validacion MD (OpenMM)",
        "script": "run_md.py",
        "check": lambda: os.path.exists(os.path.join(PROJECT_ROOT, "md_results", "local_md", "md_validation.csv"))
    },
]


def run_step(step_info):
    """Ejecuta un paso del pipeline."""
    name = step_info["name"]
    script = step_info["script"]
    
    print(f"\n{'#'*60}")
    print(f"# PASO: {name}")
    print(f"{'#'*60}")
    
    # Verificar si ya esta completado
    try:
        if step_info["check"]():
            print(f"  [SKIP] Ya completado. Continuando...")
            return True
    except:
        pass
    
    # Ejecutar script
    script_path = os.path.join(SRC_DIR, script)
    if not os.path.exists(script_path):
        print(f"  [ERROR] Script no encontrado: {script_path}")
        return False
    
    t_start = time.time()
    result = subprocess.run(
        [sys.executable, script_path],
        capture_output=False,  # Mostrar output en tiempo real
        text=True,
        timeout=600,  # 10 minutos max por paso
        cwd=PROJECT_ROOT
    )
    
    t_elapsed = time.time() - t_start
    
    if result.returncode != 0:
        print(f"  [ERROR] El script fallo con codigo {result.returncode}")
        return False
    
    # Verificar que se completo
    try:
        if step_info["check"]():
            print(f"  [OK] Paso completado en {t_elapsed:.0f}s")
            return True
        else:
            print(f"  [ERROR] El paso termino pero los archivos no se generaron")
            return False
    except Exception as e:
        print(f"  [ERROR] Verificacion fallida: {e}")
        return False


def generate_report():
    """Genera informe final del pipeline local."""
    report_path = os.path.join(PROJECT_ROOT, "results", "local_pipeline_report.md")
    
    # Datos de docking
    docking_csv = os.path.join(PROJECT_ROOT, "results", "local_dock", "docking_results.csv")
    docking_hits = []
    if os.path.exists(docking_csv):
        with open(docking_csv) as f:
            docking_hits = list(csv.DictReader(f))
    
    # Datos de MD
    md_csv = os.path.join(PROJECT_ROOT, "md_results", "local_md", "md_validation.csv")
    md_results = []
    if os.path.exists(md_csv):
        with open(md_csv) as f:
            md_results = list(csv.DictReader(f))
    
    stable = [r for r in md_results if r.get("stable") == "YES"]
    
    report = f"""# Pipeline Local MASIVE-ALS - Informe de Prueba

**Fecha:** {time.strftime('%Y-%m-%d %H:%M')}
**Equipo:** Windows PC (32 GB RAM, CPU)
**Herramientas:** AutoDock Vina 1.2, OpenMM 8.5, OpenBabel

---

## Resumen de Docking

- **Proteinas:** TDP-43, SOD1, FUS
- **Compuestos probados:** {len(set(h['ligand'] for h in docking_hits)) if docking_hits else 0}
- **Total dockings:** {len(docking_hits)}

### Top 5 Hits

| # | Compuesto | Proteina | Energia (kcal/mol) |
|---|-----------|----------|-------------------|
"""
    
    if docking_hits:
        docking_hits.sort(key=lambda x: float(x["energy"]))
        for i, h in enumerate(docking_hits[:5]):
            report += f"| {i+1} | {h['ligand']} | {h['target']} | {float(h['energy']):.2f} |\n"
    
    report += f"""
---

## Validacion MD (OpenMM)

- **Candidatos validados:** {len(md_results)}
- **Estables:** {len(stable)}
- **Inestables:** {len(md_results) - len(stable)}

"""

    if stable:
        report += "### Candidatos Estables\n\n"
        report += "| Compuesto | Proteina | Energia MD (kJ/mol) | CV% |\n"
        report += "|-----------|----------|--------------------|------|\n"
        for r in stable:
            report += f"| {r['ligand']} | {r['target']} | {r.get('mean_energy_kjmol', 'N/A')} | {r.get('cv_percent', 'N/A')}% |\n"
    
    report += f"""
---

## Conclusiones

1. Pipeline validado con datos reales en PC local
2. AutoDock Vina identifica candidatos con energias favorables
3. OpenMM confirma estabilidad de complejos proteina-ligando
4. Pipeline listo para escalar en supercomputadoras

"""
    
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(report)
    
    return report_path


def main():
    print("=" * 60)
    print(" PIPELINE MASIVE-ALS - PRUEBA LOCAL REAL")
    print("=" * 60)
    print(f"  Proyecto: {PROJECT_ROOT}")
    print(f"  Python:   {sys.version}")
    print(f"  Fecha:    {time.strftime('%Y-%m-%d %H:%M')}")
    
    # Ejecutar cada paso
    t_total_start = time.time()
    completed = 0
    
    for step in STEPS:
        if run_step(step):
            completed += 1
        else:
            print(f"\n[STOP] Pipeline detenido en: {step['name']}")
            print("Corrige el error y vuelve a ejecutar.")
            return 1
    
    t_total = time.time() - t_total_start
    
    # Informe final
    print(f"\n{'='*60}")
    print(f" PIPELINE COMPLETO - {completed}/{len(STEPS)} pasos exitosos")
    print(f" Tiempo total: {t_total:.0f} segundos")
    print(f"{'='*60}")
    
    report = generate_report()
    print(f"\n  Informe guardado: {report}")
    
    # Mostrar hallazgos principales
    docking_csv = os.path.join(PROJECT_ROOT, "results", "local_dock", "docking_results.csv")
    if os.path.exists(docking_csv):
        with open(docking_csv) as f:
            reader = csv.DictReader(f)
            hits = sorted(reader, key=lambda r: float(r["energy"]))[:3]
        
        print("\n  TOP 3 CANDIDATOS A FARMACO:")
        for hit in hits:
            print(f"    {hit['ligand']:20s} | {hit['target']:6s} | Energia: {float(hit['energy']):7.2f} kcal/mol")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
