"""MM-GBSA de una pose de TBK1, en Oracle. Prueba del pipeline.

CAMPO DE FUERZAS
El receptor con Amber14 (amber14-all.xml) y el ligando con GAFF. El ligando
de Vina llega como residuo UNL, que Amber no sabe tipar, asi que se usa
GaffTemplateGenerator, que es la combinacion estandar para esto: proteina
Amber, ligando GAFF, agua TIP3P, solvatacion implicita GBn2.

NOTA DE HIDROGENOS
Al ligando NO se le anaden atomos que falten con pdbfixer: eso lo trunca
(una pose de 33 atomos se queda en 4). Los hidrogenos los pone GAFF, que
tipa el ligando entero, y eso incluye los que faltaban.

El MM-GBSA es una diferencia de tres sistemas:

    dG_union = G(complejo + agua) - G(receptor + agua) - G(ligando, implicito)

Uso:  python3 mmgbsa_pose.py receptor.pdb pose.pdb
"""
import sys
import time

from openmm import (Context, LocalEnergyMinimizer, Platform, VerletIntegrator,
                    app, unit)
from openmm.app import ForceField, Modeller, PDBFile
from pdbfixer import PDBFixer

REC = sys.argv[1]
LIG = sys.argv[2]

CAMPOS = "amber14-all.xml"
IMPLICITO = "implicit/gbn2.xml"
PADDING_NM = 1.0
IONICO_M = 0.15


def sin_model(entrada, unico=False):
    """Limpia el PDB para OpenMM.

    - quita MODEL/ENDMDL: Open Babel escribe 'MODEL 1' en vez de
      'MODEL        1', y OpenMM lee el numero de modelo en las columnas
      11-14, que ahi estan vacias.
    - con unico=True renombra los atomos: Vina nombra los atomos del ligando
      por tipo (C, C, C...), y OpenMM los ve como duplicados dentro del mismo
      residuo y los tira. Una pose de 33 atomos se queda en 4. Por eso el
      ligando necesita nombre unico por atomo.
    """
    limpio = entrada + (".1lig" if unico else ".1model")
    contador = 0
    with open(entrada, errors="replace") as origen, \
            open(limpio, "w", encoding="utf-8") as destino:
        for cruda in origen:
            l = cruda.rstrip("\n").rstrip("\r")
            if l.startswith(("MODEL", "ENDMDL")):
                continue
            if unico and l.startswith(("ATOM", "HETATM")):
                contador += 1
                elemento = l[76:78].strip() or l[12:16].strip()[0]
                nombre = ("%-4s" % (elemento + str(contador)[-3:]))
                l = l[:12] + nombre + l[16:]
            destino.write(l + "\n")
    return limpio


t0 = time.time()
fixer = PDBFixer(filename=sin_model(REC))
fixer.findMissingResidues()
fixer.findNonstandardResidues()
fixer.replaceNonstandardResidues()
fixer.removeHeterogens(keepWater=False)
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.4)
rec_top, rec_pos = fixer.topology, fixer.positions
print("receptor: %d atomos  (%.0f s)" % (sum(1 for _ in rec_top.atoms()),
                                          time.time() - t0), flush=True)

# El ligando: nombres unicos por atomo, y sin que pdbfixer lo toque.
# GAFF lo tipa despues.
lig = PDBFile(sin_model(LIG, unico=True))
lig_top, lig_pos = lig.topology, lig.positions
n_lig = sum(1 for _ in lig_top.atoms())
print("ligando:  %d atomos" % n_lig, flush=True)
if n_lig < 5:
    print("El ligando se quedo sin atomos tras renombrar. Revisar la pose.")
    sys.exit(1)

ff = ForceField(CAMPOS, IMPLICITO)

from openmmforcefields.generators import GAFFTemplateGenerator
from openmmforcefields.utils import Molecule

# Esta API espera objetos Molecule, no un Topology.
mol = Molecule.from_topology(lig_top, assign_charges=True)
gaff = GAFFTemplateGenerator(molecules=[mol])
ff.registerTemplateGenerator(gaff.generator)
print("GAFF %d.%d: %d templates"
      % (gaff.gaff_major_version, gaff.gaff_minor_version,
         len(gaff.generator.templates)), flush=True)


def energia(con_ligando, agua):
    m = Modeller(rec_top, rec_pos)
    if con_ligando:
        m.add(lig_top, lig_pos)
    if agua:
        m.addSolvent(ff, padding=PADDING_NM * unit.nanometer,
                     ionicStrength=IONICO_M * unit.molar, neutralize=True)
    s = ff.createSystem(m.topology, nonbondedMethod=app.PME,
                        nonbondedCutoff=1.0 * unit.nanometer,
                        constraints=app.HBonds, rigidWater=True)
    c = Context(s, VerletIntegrator(1 * unit.femtosecond),
                Platform.getPlatformByName("CPU"))
    c.setPositions(m.positions)
    LocalEnergyMinimizer.minimize(c, maxIterations=200)
    e = c.getState(getEnergy=True).getPotentialEnergy()
    valor = e.value_in_unit(unit.kilocalorie_per_mole)
    c.dispose()
    return valor


t0 = time.time()
gc = energia(True, True)
print("G(complejo con agua) = %12.1f kcal/mol   (%.0f s)" % (gc, time.time() - t0),
      flush=True)

t0 = time.time()
gr = energia(False, True)
print("G(receptor con agua) = %12.1f kcal/mol   (%.0f s)" % (gr, time.time() - t0),
      flush=True)

t0 = time.time()
gl = energia(True, False)
print("G(ligando implicito) = %11.1f kcal/mol   (%.0f s)" % (gl, time.time() - t0),
      flush=True)

print()
print("dG_union MM-GBSA = %.1f kcal/mol" % (gc - gr - gl))
