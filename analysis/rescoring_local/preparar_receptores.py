# -*- coding: utf-8 -*-
"""Prepara y CONGELA el receptor de cada diana del rescoring MM-GBSA.

Por que (ver RESCORING_WSL_CUDA_2026-09-11.md):
  * PDBFixer coloca los hidrogenos de forma NO determinista: el receptor
    cambiaba en cada corrida (11-48 kcal/mol de deriva).
  * El runner recortaba la bolsa alrededor de la pose de CADA ligando, de
    modo que dos compuestos se evaluaban contra atomos distintos.

Criterio del bolsillo (ver bolsillo.py): union de esferas de radio `radio`
(10 A por defecto) alrededor de TODOS los atomos de TODAS las poses de la
diana en la lista corta. Se calcula una vez y se congela, de modo que todos
los compuestos se evaluan contra exactamente los mismos atomos y el
bolsillo contiene siempre al ligando completo.

Despues se prepara con PDBFixer y se relaja con restricciones posicionales
sobre los atomos pesados: se quita la tension de los hidrogenos colocados
por PDBFixer sin que las cadenas laterales se hundan en el bolsillo.

Artefacto: analysis/rescoring_local/receptores_fijos/<target>_fijo.pdb
           + receptores_fijos_reporte.json (con md5 y metricas).

Uso (desde WSL):
  ~/rescoring_env/bin/python preparar_receptores.py --radio 10
  ~/rescoring_env/bin/python preparar_receptores.py --targets SOD1
"""
import argparse
import csv
import hashlib
import json
import os
import sys
import tempfile
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M        # noqa: E402
import mmgbsa_openff_gb as G     # noqa: E402
import bolsillo as B            # noqa: E402

BASE = r"C:\Users\Fredy\masive-als"
if not os.path.isdir(BASE):
    BASE = "/mnt/c/Users/Fredy/masive-als"
OUTDIR = os.path.join(BASE, "analysis", "rescoring_local", "receptores_fijos")
LISTA = os.path.join(BASE, "analysis", "lista_corta_candidatos.csv")
CONTROLES = os.path.join(BASE, "analysis", "controles_calibracion.csv")

DIANAS = {
    "SOD1":     os.path.join(BASE, "gpu_dock", "SOD1.pdbqt"),
    "TDP43_v2": os.path.join(BASE, "gpu_dock", "TDP43_v2.pdbqt"),
    "FUS":      os.path.join(BASE, "gpu_dock", "FUS.pdbqt"),
    # Receptor antiguo, solo para puntuar los 2 controles con pose propia.
    "TDP43":    os.path.join(BASE, "gpu_dock", "TDP43.pdbqt"),
}

K_PESADOS = float(os.environ.get("MMGBSA_K_RESTRAINT", "2000.0"))  # kJ/mol/nm^2
MAX_ITER = int(os.environ.get("MMGBSA_MAXITER", "8000"))
TOL = float(os.environ.get("MMGBSA_TOL", "1.0"))
PLATFORM = os.environ.get("MMGBSA_FINAL_PLATFORM", "CPU")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:12]


def ligandos_de(target):
    """Ligandos con pose para esa diana (lista corta; controles para TDP43)."""
    ligs = []
    if target in ("SOD1", "TDP43_v2", "FUS"):
        with open(LISTA, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if r["target"] == target:
                    ligs.append(r["ligand"])
    elif os.path.exists(CONTROLES):
        with open(CONTROLES, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                if (r["target"] == target
                        and r.get("en_libreria", "").strip().lower() == "si"
                        and (r.get("ligand_libreria") or "").strip()):
                    ligs.append(r["ligand_libreria"].strip())
    return ligs


def add_restraint(system, topology, positions, k_kj):
    """CustomExternalForce armonico sobre los atomos pesados (ancla = posicion
    inicial). Limita la deriva de la proteina durante la relajacion."""
    from openmm import CustomExternalForce, unit

    f = CustomExternalForce("0.5*k*periodicdistance(x,y,z,x0,y0,z0)^2")
    f.addGlobalParameter("k", k_kj)
    for p in ("x0", "y0", "z0"):
        f.addPerParticleParameter(p)
    nm = positions.value_in_unit(unit.nanometer)
    for i, at in enumerate(topology.atoms()):
        if at.element is not None and at.element.symbol != "H":
            f.addParticle(i, [float(nm[i][0]), float(nm[i][1]), float(nm[i][2])])
    system.addForce(f)
    return f


def minimizar(system, topology, positions, platform_name, max_iter, tol):
    from openmm import LangevinIntegrator, LocalEnergyMinimizer, unit
    from openmm.app import Simulation

    plat = G._platform_by_name(platform_name) or M.get_platform()
    G.configurar_plataforma(plat)
    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                               2 * unit.femtoseconds)
    sim = Simulation(topology, system, integ, plat)
    sim.context.setPositions(positions)
    LocalEnergyMinimizer.minimize(
        sim.context, tolerance=tol * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=max_iter)
    return sim.context.getState(getPositions=True).getPositions()


def preparar(target, radio, relajar=True, max_poses=None):
    from openff.toolkit import ForceField as OFFForceField
    from openmm.app import PDBFile

    info = {"target": target, "pdbqt": DIANAS[target], "radio_bolsillo": radio}
    os.makedirs(OUTDIR, exist_ok=True)
    out_pdb = os.path.join(OUTDIR, "%s_fijo.pdb" % target)

    t0 = time.time()
    ligs = ligandos_de(target)
    centros, faltan = B.atomos_de_poses(target, ligs, max_poses)
    info["poses_usadas"] = len(ligs) - faltan
    info["poses_sin_archivo"] = faltan
    info["atomos_ligando_centros"] = len(centros)
    if len(centros) == 0:
        raise RuntimeError("sin poses para %s" % target)

    crudo = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False).name
    info.update(B.seleccionar_bolsillo(DIANAS[target], centros, radio, crudo))
    info["t_bolsillo_s"] = round(time.time() - t0, 2)

    rec_top, rec_pos = M.prepare_receptor(crudo)
    info["atomos_con_H"] = rec_top.getNumAtoms()
    info["residuos_finales"] = sum(1 for _ in rec_top.residues())
    try:
        os.remove(crudo)
    except OSError:
        pass

    ff = OFFForceField(M.FF_PATH)
    rec_off, _ = M.off_topology_from_pdb(rec_top, rec_pos)
    sys_rec = M.build_system(rec_off, ff)
    G.add_gb_obc2(sys_rec, rec_top)
    info["e_inicial"] = round(M.potential_energy(sys_rec, rec_top, rec_pos), 2)

    pos_final = rec_pos
    if relajar:
        add_restraint(sys_rec, rec_top, rec_pos, K_PESADOS)
        t1 = time.time()
        pos_final = minimizar(sys_rec, rec_top, rec_pos, PLATFORM, MAX_ITER, TOL)
        info["t_relajar_s"] = round(time.time() - t1, 2)
        info["e_relajada"] = round(
            M.potential_energy(sys_rec, rec_top, pos_final), 2)
        info["k_restraint"] = K_PESADOS
        a0 = np.asarray(rec_pos.value_in_unit(rec_pos.unit), dtype=float)
        a1 = np.asarray(pos_final.value_in_unit(pos_final.unit), dtype=float)
        pesados = [i for i, at in enumerate(rec_top.atoms())
                   if at.element is not None and at.element.symbol != "H"]
        d = np.linalg.norm(a0[pesados] - a1[pesados], axis=1)
        info["deriva_pesados_rms_A"] = round(float(np.sqrt((d ** 2).mean())), 4)
        info["deriva_pesados_max_A"] = round(float(d.max()), 4)

    with open(out_pdb, "w") as f:
        PDBFile.writeFile(rec_top, pos_final, f)
    info["receptor_fijo"] = out_pdb
    info["md5"] = md5(out_pdb)
    info["plataforma"] = PLATFORM
    info["t_total_s"] = round(time.time() - t0, 2)
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="SOD1,TDP43_v2,FUS,TDP43")
    ap.add_argument("--radio", type=float, default=10.0,
                    help="radio en A de la esfera alrededor de cada atomo de ligando")
    ap.add_argument("--max-poses", type=int, default=None)
    ap.add_argument("--sin-relajar", action="store_true")
    a = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    rep_path = os.path.join(OUTDIR, "receptores_fijos_reporte.json")
    res = {}
    if os.path.exists(rep_path):
        try:
            with open(rep_path, encoding="utf-8") as f:
                res = json.load(f)
        except Exception:
            res = {}

    for t in [x.strip() for x in a.targets.split(",") if x.strip()]:
        print("=== %s ===" % t, flush=True)
        try:
            info = preparar(t, a.radio, relajar=not a.sin_relajar,
                            max_poses=a.max_poses)
        except Exception as e:
            info = {"target": t, "error": "%s: %s" % (type(e).__name__, e)}
        res[t] = info
        print(json.dumps(info, ensure_ascii=False, indent=2), flush=True)
        with open(rep_path, "w", encoding="utf-8") as f:
            json.dump(res, f, ensure_ascii=False, indent=2)
    print("reporte: %s" % rep_path)


if __name__ == "__main__":
    main()
