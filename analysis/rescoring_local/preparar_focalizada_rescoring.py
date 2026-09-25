# -*- coding: utf-8 -*-
"""Prepara la lista focalizada para el rescoring MM-GBSA con el protocolo validado.

Qué hace y por qué:
  1. Escribe la lista con las columnas que espera runner_mmgbsa_gb.py.
  2. **Usa el SMILES sin sal** cuando lo hay. Importa: a 12 de las 179 se les
     quitó el contraión para poder acoplarlas, y el hijo del rescoring compara el
     número de átomos del SMILES con el de la pose. Si se le pasa el SMILES con el
     cloruro, esa molécula falla en silencio con «n atomos no coincide» (es el
     error que dejó fuera a la cefarantina).
  3. Deja las poses en una carpeta propia con la estructura que espera el runner
     (results_TDP43_v2/<ligando>_out.pdbqt). NO pisa las de la librería: esas son
     la prueba de la comparación de cajas.

Uso: python preparar_focalizada_rescoring.py
"""
import csv
import os
import shutil

LOCAL = r"C:\Users\Fredy\masive-als\analysis\rescoring_local"
ENTRADA = os.path.join(LOCAL, "lista_focalizada_tdp43v2.csv")
SALIDA = os.path.join(LOCAL, "lista_focalizada_rescoring.csv")
POSES_ORIGEN = r"C:\Users\Fredy\masive-als\gpu_dock\resultados_caja_arn\results_focalizada"
POSES_DESTINO = r"C:\Users\Fredy\masive-als\gpu_dock\rescoring_caja_arn\results_TDP43_v2"
TARGET = "TDP43_v2"


def main():
    filas = list(csv.DictReader(open(ENTRADA, encoding="utf-8")))
    os.makedirs(POSES_DESTINO, exist_ok=True)

    listos, sin_pose, con_sal = [], [], 0
    for f in filas:
        lig = f["ligand"]
        smi = (f.get("smiles_sin_sal") or "").strip() or (f.get("smiles") or "").strip()
        if (f.get("smiles_sin_sal") or "").strip():
            con_sal += 1
        pose = os.path.join(POSES_ORIGEN, lig + "_out.pdbqt")
        if not os.path.exists(pose):
            sin_pose.append(lig)
            continue
        shutil.copy2(pose, os.path.join(POSES_DESTINO, lig + "_out.pdbqt"))
        listos.append({"target": TARGET, "ligand": lig, "smiles": smi,
                       "vina_affinity": f.get("afinidad_tdp43v2") or "",
                       "tipo": "focalizada_caja_arn"})

    with open(SALIDA, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["target", "ligand", "smiles", "vina_affinity", "tipo"])
        w.writeheader()
        w.writerows(listos)

    print("Lista para el rescoring: %d moleculas -> %s" % (len(listos), SALIDA))
    print("Con SMILES sin sal (el que acople): %d" % con_sal)
    print("Poses copiadas a: %s" % POSES_DESTINO)
    if sin_pose:
        print("SIN POSE (se quedan fuera): %s" % ", ".join(sin_pose[:10]))
    return 0 if listos else 1


if __name__ == "__main__":
    raise SystemExit(main())
