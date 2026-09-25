#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""validar_diana_limpia.py — la validación limpia, para SOD1 y para TDP-43.

POR QUE ES UN SOLO SCRIPT
-------------------------
`validar_sod1_limpia.py` fue la primera (21 sep 2026) y este no la sustituye: la
importa. La validación de una diana y la de otra solo se diferencian en tres cosas
—receptor, caja y carpeta de poses—, así que tener dos copias era pedir que una se
quedara atrás (que es exactamente lo que pasó con la receta de preparar ligandos, y
costó la mitad de las conclusiones falsas de este proyecto).

LO QUE MIDE, IGUAL PARA LAS DOS
-------------------------------
  * receptor limpio y caja las de la validación de esa diana;
  * positivos de **unión medida** de `verdad_de_referencia.csv`;
  * el fondo que ya esté acoplado contra ESE receptor;
  * **comprobación previa** de que fichero = pose = SMILES en átomos pesados, y
    re-preparado + re-acoplado de lo que no cuadre;
  * cuatro números (AUC crudo, por átomo, residual, EF) y **el criterio**: cuántos
    quimiotipos distintos baten a los señuelos **de su mismo tamaño**.

El criterio emparejado por tamaño es el único que decide: el crudo mide tamaño
(correlación −0,72 en SOD1) y el corte fijo de cabeza es inestable.

Uso:
  python validar_diana_limpia.py --target TDP43
  python validar_diana_limpia.py --target SOD1
"""
import argparse
import os
import sys

BASE = os.path.dirname(os.path.abspath(__file__))
RAIZ = os.path.dirname(BASE)
GPU = os.path.join(RAIZ, "gpu_dock")

# ---------------------------------------------------------------------------
# Lo unico que cambia entre dianas.
#
# SOD1: receptor corregido (dimero biologico A-H, sin las aguas ni las 18 copias
#       del cristal) y la caja Trp32 de toda la serie.
# TDP43: 4BS2 con protonacion pH 7.4 y la caja del bolsillo RRM1-RRM2 / puente
#       salino Arg151-Asp247 (la que reprodujo las energias publicadas de PE859 y
#       berberrubina). Su exhaustividad original (16) y su fondo de 122 señuelos
#       emparejados en propiedades. Ademas lleva el FONDO DURO (23 sep 2026): los
#       152 unidores de ARN de R-BIND 2.0 acoplados en la misma caja
#       (`acoplar_fondo_rbind.py`), que es el unico fondo que puede separar
#       quimia especifica de quimia generica de ARN.
# ---------------------------------------------------------------------------
DIANAS = {
    "SOD1": {
        "receptor": os.path.join(GPU, "SOD1_limpio.pdbqt"),
        "centro": (46.5, 80.0, 73.3),
        "tamano": 22,
        "exhaustividad": 8,
        "poses": os.path.join(BASE, "validacion_SOD1_v5", "out"),
        "ligands": os.path.join(BASE, "validacion_SOD1_v5", "ligands"),
        "salida": os.path.join(BASE, "validar_sod1_limpia"),
        "nota_fondo": "482 señuelos: 199 emparejados antiguos + 140 emparejados "
                      "nuevos + 143 duros (quelantes y redox-activos)",
        "poses2": None,
        "ligands2": None,
        "nota_fondo2": "",
        "etiqueta_fondo2": "fondo duro",
    },
    "TDP43": {
        "receptor": os.path.join(BASE, "_tdp43_bolsillo_v2", "4BS2_ph74.pdbqt"),
        "centro": (24.23, 16.89, -15.87),
        "tamano": 26,
        "exhaustividad": 16,
        "poses": os.path.join(BASE, "_validacion_TDP43", "out"),
        "ligands": os.path.join(BASE, "_validacion_TDP43", "ligands"),
        "salida": os.path.join(BASE, "validar_tdp43_limpia"),
        "nota_fondo": "122 señuelos emparejados en propiedades",
        "poses2": os.path.join(BASE, "_validacion_TDP43_rbind", "out"),
        "ligands2": os.path.join(BASE, "_validacion_TDP43_rbind", "ligands"),
        "nota_fondo2": "unidores de ARN de R-BIND 2.0 (Donlic 2022), acoplados "
                       "en la misma caja; ver analysis/rbind_fondo.py",
        "etiqueta_fondo2": "unidores de ARN de R-BIND 2.0 (fondo duro)",
    },
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target", required=True, choices=sorted(DIANAS))
    args = ap.parse_args()

    cfg = DIANAS[args.target].copy()
    # `validar_sod1_limpia` es el modulo que ya existe y hace todo el trabajo; aqui
    # solo se le dice para que diana. Se importa con su nombre ya cargado.
    sys.path.insert(0, BASE)
    import validar_sod1_limpia as V

    # sobrescribir la configuracion de la diana
    V.RECEPTOR = cfg["receptor"]
    V.CENTRO = cfg["centro"]
    V.TAMANO = cfg["tamano"]
    V.EXHAUSTIVIDAD = cfg["exhaustividad"]
    V.POSES_FONDO = cfg["poses"]
    V.LIGS_FONDO = cfg["ligands"]
    V.SALIDA = cfg["salida"]
    V.DIANA = args.target
    V.NOTA_FONDO = cfg["nota_fondo"]
    V.ETIQUETA_FONDO = ("señuelos emparejados en propiedades" if args.target == "SOD1"
                        else "señuelos emparejados en propiedades (fondo blando)")
    V.POSES_FONDO2 = cfg["poses2"]
    V.LIGS_FONDO2 = cfg["ligands2"]
    V.NOTA_FONDO2 = cfg["nota_fondo2"]
    V.ETIQUETA_FONDO2 = cfg["etiqueta_fondo2"]
    return V.main()


if __name__ == "__main__":
    sys.exit(main())
