# -*- coding: utf-8 -*-
"""Seleccion del bolsillo del receptor: union de esferas sobre las poses.

Por que no un cutoff esferico desde el centro de la caja: medido sobre la
lista corta, los atomos de ligando llegan a 17,7-18,2 A del centro de la
caja (las cajas van de 22 a 26 A de lado, asi que esos atomos estan dentro
de las esquinas). Un bolsillo esferico de radio 16 A dejaba al ligando mas
grande sobresaliendo del receptor, y su energia de interaccion con los
atomos que faltaban se perdia.

Criterio adoptado (determinista y fijo por diana):
  bolsillo = residuos del receptor que tienen algun atomo a menos de
             `radio` A de ALGUN atomo de ALGUNA de las poses acopladas
             de la diana en la lista corta.
    - radio = 10 A: mas alla de eso la contribucion al dG es despreciable
      (las interacciones de corto alcance dominan) y los atomos lejanos
      cancelan entre E_complejo y E_receptor.
    - Se calcula UNA vez por diana y se congela en disco; todos los
      compuestos se evaluan contra exactamente los mismos atomos.

La busqueda de vecinos usa un KD-tree (scipy); el barrido directo
N_receptor x N_centros tardaba minutos por diana.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mmgbsa_openff as M  # noqa: E402

BASE = r"C:\Users\Fredy\masive-als"
if not os.path.isdir(BASE):
    BASE = "/mnt/c/Users/Fredy/masive-als"
POSES = os.path.join(BASE, "gpu_dock", "resultados_libreria")


def atomos_de_poses(target, ligandos, max_poses=None):
    """Coordenadas de todos los atomos de las poses indicadas (N,3)."""
    if max_poses:
        ligandos = ligandos[:max_poses]
    centros = []
    faltan = 0
    for lig in ligandos:
        p = os.path.join(POSES, "results_%s" % target, lig + "_out.pdbqt")
        if not os.path.exists(p):
            faltan += 1
            continue
        try:
            at = M.parse_pdbqt_atoms(p, primer_model=True)
        except Exception:
            faltan += 1
            continue
        for a in at:
            centros.append((a["x"], a["y"], a["z"]))
    return np.asarray(centros, dtype=float), faltan


def agrupar_residuos(atoms):
    """Agrupa atomos en residuos, detectando copias (reinicios numericos)."""
    residues = []
    copy_id, prev_rn = 0, None
    for a in atoms:
        try:
            rn = int(a["resnum"])
        except (TypeError, ValueError):
            continue
        if prev_rn is not None and rn < prev_rn:
            copy_id += 1
        if (not residues or residues[-1]["resnum"] != rn
                or residues[-1]["copy"] != copy_id):
            residues.append({"resnum": rn, "resname": a["resname"],
                             "copy": copy_id, "atoms": []})
        residues[-1]["atoms"].append(a)
        prev_rn = rn
    return residues


def dist_min_por_residuo(residuos, centros):
    """Minima distancia de cada residuo al conjunto de centros (KD-tree)."""
    pts, idx = [], []
    for i, r in enumerate(residuos):
        for a in r["atoms"]:
            pts.append((a["x"], a["y"], a["z"]))
            idx.append(i)
    pts = np.asarray(pts, dtype=float)
    idx = np.asarray(idx, dtype=np.int64)
    try:
        from scipy.spatial import cKDTree
        d = cKDTree(centros).query(pts, k=1, workers=-1)[0]
    except Exception:  # respaldo sin scipy
        d = np.full(len(pts), np.inf)
        for i in range(0, len(centros), 4000):
            c = centros[i:i + 4000]
            dd = np.linalg.norm(pts[:, None, :] - c[None, :, :], axis=2)
            d = np.minimum(d, dd.min(axis=1))
    out = np.full(len(residuos), np.inf)
    np.minimum.at(out, idx, np.asarray(d, dtype=float))
    return out


def escribir_pdb(residuos, out_pdb):
    lines, serial, chain_idx, res_serial = [], 0, 0, 0
    prev = None
    for r in residuos:
        nueva = (prev is not None and
                 (r["copy"] != prev["copy"] or
                  r["resnum"] != prev["resnum"] + 1))
        if nueva:
            lines.append("TER")
            chain_idx += 1
            res_serial = 0
        res_serial += 1
        for a in r["atoms"]:
            serial += 1
            lines.append(
                "ATOM  %5d %4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
                % (serial, a["name"][:4], a["resname"][:3],
                   chr(ord("A") + chain_idx % 26), res_serial,
                   a["x"], a["y"], a["z"], a["elem"][:2]))
        prev = r
    lines.append("TER")
    lines.append("END")
    with open(out_pdb, "w") as f:
        f.write("\n".join(lines) + "\n")
    return out_pdb


def seleccionar_bolsillo(receptor_pdbqt, centros, radio, out_pdb):
    """Escribe el bolsillo (residuos a < radio de algun centro) como PDB."""
    atoms = M.parse_pdbqt_atoms(receptor_pdbqt)
    if not atoms:
        raise RuntimeError("receptor sin atomos: %s" % receptor_pdbqt)
    atoms = [a for a in atoms if a["resname"] in M.STANDARD_AA]
    if not atoms:
        raise RuntimeError("receptor sin aminoacidos estandar")
    residuos = agrupar_residuos(atoms)
    if not residuos:
        raise RuntimeError("receptor sin residuos validos")

    d_res = dist_min_por_residuo(residuos, centros)

    # Por cada numero de residuo se conserva la copia mas cercana a las poses
    # (SOD1 concatena ~18 copias; las lejanas cancelarian en dG).
    info = {}
    for r, d in zip(residuos, d_res):
        k = r["resnum"]
        if k not in info or d < info[k][0]:
            info[k] = (float(d), r)

    kept = [(d, r) for d, r in info.values() if d <= radio]
    kept.sort(key=lambda t: (t[1]["copy"], t[1]["resnum"]))
    escribir_pdb([r for _, r in kept], out_pdb)
    return {"residuos_receptor": len(residuos), "residuos_unicos": len(info),
            "residuos_bolsillo": len(kept), "radio": radio,
            "d_min_global": round(min(d for d, _ in info.values()), 2),
            "d_max_kept": round(max((d for d, _ in kept), default=0), 2)}
