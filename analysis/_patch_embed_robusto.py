# -*- coding: utf-8 -*-
"""Parche: embedding robusto en prepare_ligand (Oracle)."""
import io

PATH = "/home/ubuntu/mmgbsa/rescoring_mmgbsa_robusto.py"
with io.open(PATH, encoding="utf-8") as f:
    txt = f.read()

old = """    rdmol_h = Chem.AddHs(rdmol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    params.maxIterations = 500
    params.SetCoordMap(coordmap)
    AllChem.EmbedMolecule(rdmol_h, params)

    offmol = Molecule.from_rdkit(rdmol_h, allow_undefined_stereo=True)"""

new = """    rdmol_h = Chem.AddHs(rdmol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.useRandomCoords = True
    params.maxIterations = 500
    params.SetCoordMap(coordmap)
    ret = AllChem.EmbedMolecule(rdmol_h, params)
    # Si el embedding con restricciones de coordenadas falla (ret != 0),
    # reintentar sin coordmap con geometria aleatoria.
    if ret != 0:
        p2 = AllChem.ETKDGv3()
        p2.randomSeed = 43
        p2.useRandomCoords = True
        p2.maxIterations = 2000
        ret = AllChem.EmbedMolecule(rdmol_h, p2)
    if ret != 0:
        raise RuntimeError(f"embedding fallo (ret={ret}) para {smiles[:40]}")

    offmol = Molecule.from_rdkit(rdmol_h, allow_undefined_stereo=True)"""

if old in txt:
    txt = txt.replace(old, new)
    with io.open(PATH, "w", encoding="utf-8") as f:
        f.write(txt)
    print("OK: parche embed aplicado")
else:
    print("AVISO: bloque no encontrado (puede estar ya parcheado)")
