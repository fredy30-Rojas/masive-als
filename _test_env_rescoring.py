import openff.toolkit
import openmm
import pdbfixer
import rdkit
import numpy
import pandas
print("openff.toolkit", openff.toolkit.__version__)
print("openmm", openmm.version.version)
from openff.toolkit import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
print("Molecule + ForceField import OK")