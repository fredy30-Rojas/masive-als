# -*- coding: utf-8 -*-
"""¿Se aplica DeterministicForces? Mide la misma energía varias veces sobre
las MISMAS coordenadas y dos minimizaciones desde el mismo punto."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M
import mmgbsa_openff_gb as G
import numpy as np

POSE = "/mnt/c/Users/Fredy/masive-als/gpu_dock/resultados_libreria/results_SOD1/CHEMBL4584906_out.pdbqt"
REC = "/mnt/c/Users/Fredy/masive-als/gpu_dock/SOD1.pdbqt"
SM = "N=C1N/C(=N/NC(=O)c2ccc(CN3C(=O)c4cccc5cccc3c45)cc2)c2ccccc21"

from openff.toolkit import ForceField as OFFForceField
from openff.toolkit import Topology as OFFTopology
from openmm.app import Modeller
from openmm import LangevinIntegrator, LocalEnergyMinimizer, unit, Platform
from openmm.app import Simulation
import tempfile

ff = OFFForceField(M.FF_PATH)
lig_top, lig_pos, offmol = G.prepare_ligand_pose(SM, POSE, ff)
lig_center = np.array([[a["x"], a["y"], a["z"]]
                       for a in M.parse_pdbqt_atoms(POSE)]).mean(axis=0)
tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
M.extract_receptor_pdb(REC, lig_center, tmp)
rec_top, rec_pos = M.prepare_receptor(tmp)
os.remove(tmp)
rec_off, rec_mols = M.off_topology_from_pdb(rec_top, rec_pos)
complex_mod = Modeller(rec_top, rec_pos)
complex_mod.add(lig_top, lig_pos)
complex_off = OFFTopology()
for m in rec_mols:
    complex_off.add_molecule(m)
complex_off.add_molecule(offmol)
sysx = M.build_system(complex_off, ff)
G.add_gb_obc2(sysx, complex_mod.topology)

# posiciones congeladas: las mismas para todas las medidas
pos = complex_mod.positions


def energia(plat, props):
    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                               2 * unit.femtoseconds)
    sim = Simulation(complex_mod.topology, sysx, integ, plat, props)
    sim.context.setPositions(pos)
    p = sim.context.getPlatform()
    vals = {}
    for k in ("DeterministicForces", "Precision", "UseCpuPme"):
        try:
            vals[k] = p.getPropertyValue(sim.context, k)
        except Exception:
            vals[k] = "n/d"
    e = sim.context.getState(getEnergy=True).getPotentialEnergy() \
        .value_in_unit(unit.kilocalorie_per_mole)
    return e, vals


for nombre, props in (("precisión mixta, DeterministicForces=true",
                       {"DeterministicForces": "true"}),
                      ("precisión doble, DeterministicForces=true",
                       {"DeterministicForces": "true", "Precision": "double"})):
    print("\n=== %s ===" % nombre)
    plat = Platform.getPlatformByName("CUDA")
    es = []
    for _ in range(4):
        e, vals = energia(plat, props)
        es.append(e)
    print("  propiedades del contexto:", vals)
    print("  energías a coordenadas FIJAS:",
          " ".join("%.4f" % x for x in es))
    print("  dispersión: %.6f kcal/mol" % (max(es) - min(es)))

    # dos minimizaciones desde el mismo punto
    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                               2 * unit.femtoseconds)
    sim = Simulation(complex_mod.topology, sysx, integ, plat, props)
    sim.context.setPositions(pos)
    for v in (1, 2):
        sim.context.setPositions(pos)
        LocalEnergyMinimizer.minimize(
            sim.context, tolerance=1.0 * unit.kilojoule_per_mole / unit.nanometer,
            maxIterations=8000)
        e = sim.context.getState(getEnergy=True).getPotentialEnergy() \
            .value_in_unit(unit.kilocalorie_per_mole)
        print("  minimización %d -> %.2f" % (v, e))
