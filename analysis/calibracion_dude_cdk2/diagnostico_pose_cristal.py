#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Por que el motor no devuelve la pose del cristal de CDK2. Sin usar la GPU.

LA PREGUNTA
-----------
Con la medida arreglada (29 sep 2026) el control queda claro: la entrada es la pose
del cristal a 0,0000 A, y las poses que devuelve Vina-GPU caen a 3,67-6,87 A en el
bolsillo. No es la metrica. Entonces es una de estas tres cosas, y hay que saber
cual antes de tocar nada:

  (a) **BUSQUEDA** — la pose del cristal esta en el paisaje pero el motor no la
      encuentra (caja grande, pocas poses, `search_depth` corto);
  (b) **PUNTUACION** — la pose esta y el motor la descarta porque le da peor
      energia que otra (el sitio de ATP de CDK2 es plano y simetrico: un ligando
      con dos brazos puede darse la vuelta sin salirse);
  (c) **PREPARACION** — la pose del cristal NO CABE en el receptor que hemos
      preparado (un choque, una protonacion, un atomo de mas), asi que el motor
      hace lo unico que puede: apartarla.

QUE MIRA ESTE FICHERO, TODO SIN GPU
-----------------------------------
1. **Cabia la pose del cristal?** Contactos de la pose del cristal contra el
   receptor preparado, con la distancia minima de atomos pesados. Si hubiera un
   choque de verdad, (c) es la respuesta y no hay mas que hablar. Ya se midio un
   minimo de 2,12 A (un enlace de hidrogeno), pero se vuelve a medir aqui y al
   lado de la pose acoplada, que es como se ve si el motor esta evitando algo.
2. **Se conserva el anclaje?** El ligando de 1h00 se ancla por la
   diaminopirimidina al enganche de LEU83 y GLU81. Si la pose acoplada mantiene
   esos contactos y se ha dado la vuelta, es (b): el motor ancla bien y elige mal
   el brazo. Si los pierde, es que se ha ido a otra parte del bolsillo.
3. **Que se movio?** RMSD por fragmento (nucleo de pirimidina, fenilo difluorado,
   resto) de la mejor pose contra el cristal. Un brazo girado 180 grados da las
   tres cuentas muy distintas, y eso es la firma de (b) en un bolsillo plano.
4. **Estaba tensa la pose del cristal?** Energia MMFF94 del conformero del cristal
   frente a la de las poses acopladas, con la MISMA molecula y los mismos
   hidrogenos relajados. Si el conformero del cristal esta varios kcal/mol por
   encima de lo que el motor encuentra, el motor esta acertando al no elegirlo, y
   el problema es la estructura de partida (o el refinamiento del cristal).

LO QUE NO PUEDE RESPONDER, Y POR QUE
------------------------------------
Separar (a) de (b) del todo exige volver a acoplar con una caja diminuta, pegada a
la pose del cristal: si con 6 A la reproduce, era busqueda; si ni asi, era
puntuacion. Eso necesita la tarjeta y esta ocupada con TBK1, asi que queda
escrito en `siguiente_prueba_caja_estrecha()` para lanzarlo cuando este libre.

Uso:
    python diagnostico_pose_cristal.py
Salida: diagnostico_pose_cristal.txt
"""
from __future__ import annotations

import glob
import os
import time

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

import control_redocking_cdk2 as C   # noqa: E402  (la receta, una sola copia)

BASE = C.BASE
SALIDA = os.path.join(BASE, "diagnostico_pose_cristal.txt")

# Los fragmentos del ligando de 1h00, por indice del mol2 (0-based). El dibujo es
# una diaminopirimidina que se ancla al enganche, con un fenilo difluorado por un
# lado y un resto saturado con el amonio por el otro.
FRAGMENTOS = {
    "nucleo pirimidina": [9, 10, 17, 18, 21, 22],
    "fenilo difluorado": [2, 3, 4, 15, 16, 19, 28, 29],
    "anillo fusionado": [5, 6, 7, 8, 13, 14],
    "resto saturado + amonio": [0, 1, 11, 12, 20, 23, 24, 25, 26, 27],
}
# El anclaje canonico de CDK2: el NH de LEU83 y el CO de GLU81.
ANCLAJE = {("LEU", 83), ("GLU", 81), ("HIE", 84), ("LEU", 134)}


def log(lineas, m):
    print(m, flush=True)
    lineas.append(m)


def proteinas():
    """(coords, residuo, elemento) de los atomos pesados del receptor preparado.

    Las dos cosas van separadas a proposito: `residuo` (nombre y numero) sirve
    para decir QUE zona del bolsillo se toca, y `elemento` para saber cuales son
    polares. Mezclarlas fue un error el 29 de septiembre: se filtro por el nombre
    del residuo buscando "N" o "O" y no salia ninguno, porque `l[17:20]` es LEU,
    GLU, etc. La columna del elemento es la 77-78, no la 18-20.
    """
    coords, res, elem = [], [], []
    for l in open(os.path.join(BASE, "receptor_cdk2.pdb"), encoding="utf-8",
                  errors="ignore"):
        if not l.startswith("ATOM"):
            continue
        if l[76:78].strip().upper() == "H":
            continue
        coords.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        res.append((l[17:20].strip().upper(), int(l[22:26])))
        elem.append(l[76:78].strip().upper())
    return np.array(coords), res, elem


def polares(lig, idx_lig, prot, idx_prot, corte=3.6):
    """Contactos de enlace de hidrogeno: N/O/F del ligando contra N/O del receptor.

    Es la cuenta que importa para el brazo del amonio: si la pose del cristal
    COLOCA ese cation haciendo mas enlaces de los que hace la pose del motor, y
    aun asi el motor la descarta, lo que falla no es que falten contactos, es lo
    que la funcion de puntuacion hace con ellos.

    Los indices vienen en el MISMO orden que las coordenadas que se le pasan
    (los dos, ligando y receptor, en orden del mol2 y del pdb respectivamente).
    """
    if not idx_lig or not idx_prot:
        return 0, 99.0
    d = np.linalg.norm(np.asarray(lig)[idx_lig][:, None, :]
                       - prot[idx_prot][None, :, :], axis=2)
    n = int(np.count_nonzero(d < corte))
    return n, float(d.min())


def contactos(lig, prot, res, corte=4.0):
    """{residuo: distancia minima} de los atomos del receptor a menos de `corte`."""
    d = np.linalg.norm(prot[:, None, :] - lig[None, :, :], axis=2)
    salida = {}
    for i, r in enumerate(res):
        m = float(d[i].min())
        if m < corte:
            salida[r] = min(salida.get(r, 99.0), m)
    return salida


def energia_mmff(mol_pesado, coords):
    """(energia del esqueleto, minimo local) de una geometria de atomos pesados.

    DOS CUIDADOS QUE NO SON OPCIONALES:

    * los hidrogenos se generan DESDE LA GEOMETRIA QUE SE MIDE, no se reutilizan
      los del cristal. Al mover los atomos pesados y dejar los hidrogenos donde
      estaban, la molecula sale con choques que no existen y la energia se dispara
      a cientos de kcal/mol (paso el 29 de septiembre: daba 137-219 donde tenia
      que dar decenas).

    * la molecula tiene que venir AROMATIZADA, que es como MMFF sabe tipar los
      anillos. La de `leer_ligando()` no lo esta (se salta la aromatizacion a
      proposito, para que las simetrias del RMSD sean las del cristal), asi que
      aqui se usa una copia sanitizada.

    Con los atomos pesados FIJOS y los hidrogenos sueltos se mide la geometria
    del esqueleto, que es lo que se quiere comparar. El minimo completo de esa
    misma molecula da la referencia de tension cero.
    """
    m = Chem.Mol(mol_pesado)
    cf = m.GetConformer()
    for k, p in enumerate(coords):
        cf.SetAtomPosition(k, (float(p[0]), float(p[1]), float(p[2])))
    mh = Chem.AddHs(m, addCoords=True)
    props = AllChem.MMFFGetMoleculeProperties(mh)
    if props is None:
        return None, None
    ff = AllChem.MMFFGetMoleculeForceField(mh, props)
    for a in mh.GetAtoms():
        if a.GetAtomicNum() > 1:
            ff.AddFixedPoint(a.GetIdx())
    ff.Minimize(maxIts=5000)
    e_esqueleto = float(ff.CalcEnergy())
    ff2 = AllChem.MMFFGetMoleculeForceField(mh, props)
    ff2.Minimize(maxIts=5000)
    return e_esqueleto, float(ff2.CalcEnergy())


def main():
    lineas = []
    log(lineas, "POR QUE EL MOTOR NO DEVUELVE LA POSE DEL CRISTAL DE CDK2")
    log(lineas, "   %s | todo sin GPU" % time.strftime("%Y-%m-%d %H:%M"))
    log(lineas, "")

    mol, ref = C.leer_ligando()
    prot, res, elem = proteinas()
    # Los atomos polares del ligando se sacan del mol2 una sola vez: es la MISMA
    # molecula en todas las poses, asi que los indices valen para todas.
    idx_lig_pol = [a.GetIdx() for a in mol.GetAtoms()
                   if a.GetSymbol().upper() in ("N", "O", "F")]
    idx_prot_pol = [i for i, e in enumerate(elem) if e in ("N", "O")]
    log(lineas, "receptor preparado: %d atomos pesados | ligando: %d (%d polares)"
        % (len(prot), mol.GetNumAtoms(), len(idx_lig_pol)))

    # --- 1. cabia la pose del cristal en el receptor preparado? ---
    con_cristal = contactos(ref, prot, res)
    n_pol_cristal, dmin_pol_cristal = polares(ref, idx_lig_pol, prot, idx_prot_pol)
    dmin = min(con_cristal.values()) if con_cristal else 99.0
    log(lineas, "")
    log(lineas, "1. CABIA LA POSE DEL CRISTAL?")
    log(lineas, "   contactos a menos de 4 A: %d residuos | distancia minima de"
        " atomos pesados %.2f A" % (len(con_cristal), dmin))
    log(lineas, "   anclaje (LEU83, GLU81, HIE84, LEU134): %s"
        % (", ".join("%s%d %.2f" % (r[0], r[1], d)
                     for r, d in sorted(con_cristal.items()) if r in ANCLAJE)
           or "NINGUNO"))
    log(lineas, "   contactos polares (N/O/F del ligando contra N/O del receptor"
        " a menos de 3,6 A): %d | el mas corto %.2f A"
        % (n_pol_cristal, dmin_pol_cristal))
    if dmin < 2.0:
        log(lineas, "   AVISO: hay un contacto por debajo de 2 A. Eso NO es un"
            " choque (es un enlace de")
        log(lineas, "   hidrogeno, 2,1 A es lo normal) pero conviene mirarlo si el"
            " motor rehuye esa zona.")

    # --- 2 y 3. que hace el motor con la mejor pose de cada variante ---
    log(lineas, "")
    log(lineas, "2 y 3. QUE DEVUELVE EL MOTOR")
    carpetas = sorted(d for d in glob.glob(os.path.join(BASE, "_control_cdk2", "*"))
                      if os.path.isdir(d))
    resumen = []
    for carpeta in carpetas:
        etiqueta = os.path.basename(carpeta)
        entrada = os.path.join(carpeta, "ligands", "FAP_cristal.pdbqt")
        salidas = sorted(glob.glob(os.path.join(carpeta, "out", "*_out.pdbqt")))
        if not os.path.exists(entrada) or not salidas:
            continue
        tipos, coords = C.atomos_pesados_pdbqt(entrada)
        perm, rmsd_entrada = C.emparejar_orden(tipos, coords, mol, ref)
        if rmsd_entrada > C.TOL_ENTRADA:
            continue
        poses = []
        for l in open(salidas[0], encoding="utf-8", errors="ignore"):
            if l.startswith("MODEL"):
                poses.append({"energia": None, "pos": []})
            elif l.startswith("REMARK VINA RESULT") and poses and \
                    poses[-1]["energia"] is None:
                poses[-1]["energia"] = float(l.split()[3])
            elif l.startswith(("ATOM", "HETATM")) and poses:
                if l.rsplit(None, 1)[-1].strip().upper() in C.H_TIPOS:
                    continue
                poses[-1]["pos"].append([float(l[30:38]), float(l[38:46]),
                                         float(l[46:54])])
        validas = [p for p in poses if len(p["pos"]) == len(ref)]
        if not validas:
            continue
        for p in validas:
            pos = np.array(p["pos"])
            # a orden del mol2 para poder hablar de atomos con nombre
            inv = np.empty(len(perm), dtype=int)
            inv[perm] = np.arange(len(perm))
            p["en_orden_mol2"] = pos[inv]
            p["rmsd"] = float(np.sqrt(((p["en_orden_mol2"] - ref) ** 2)
                                      .sum(1).mean()))
        mejor = min(validas, key=lambda p: p["rmsd"])
        pos = mejor["en_orden_mol2"]

        log(lineas, "")
        log(lineas, "   %s | mejor pose %d (%+.2f kcal/mol) a %.2f A"
            % (etiqueta, validas.index(mejor) + 1, mejor["energia"],
               mejor["rmsd"]))
        # por fragmento
        for nombre, idx in FRAGMENTOS.items():
            r = float(np.sqrt(((pos[idx] - ref[idx]) ** 2).sum(1).mean()))
            log(lineas, "      %-26s %5.2f A" % (nombre, r))
        # contactos
        con_pose = contactos(pos, prot, res)
        n_pol, dmin_pol = polares(pos, idx_lig_pol, prot, idx_prot_pol)
        log(lineas, "      contactos polares %d (el cristal tiene %d) | el mas"
            " corto %.2f A" % (n_pol, n_pol_cristal, dmin_pol))
        comunes = set(con_cristal) & set(con_pose)
        perdidos = sorted(set(con_cristal) - set(con_pose))
        ganados = sorted(set(con_pose) - set(con_cristal))
        log(lineas, "      contactos comunes %d de %d | perdidos %d | nuevos %d"
            % (len(comunes), len(con_cristal), len(perdidos), len(ganados)))
        ancla_sigue = sorted(r for r in comunes if r in ANCLAJE)
        log(lineas, "      el anclaje que se mantiene: %s"
            % (", ".join("%s%d" % r for r in ancla_sigue) or "NINGUNO"))
        if perdidos:
            log(lineas, "      contactos que se pierden: %s"
                % ", ".join("%s%d" % r for r in perdidos[:8]))
        if ganados:
            log(lineas, "      contactos nuevos: %s"
                % ", ".join("%s%d" % r for r in ganados[:8]))
        resumen.append((etiqueta, mejor, pos, ancla_sigue))

    # --- 4. estaba tensa la pose del cristal? ---
    log(lineas, "")
    log(lineas, "4. ESTABA TENSA LA POSE DEL CRISTAL? (MMFF94, hidrogenos"
        " relajados, mismas condiciones)")
    # Copia AROMATIZADA, en el mismo orden de atomos que el mol2: es la que MMFF
    # sabe tipar. La de `leer_ligando()` se salta la aromatizacion a proposito.
    mol_mmff = Chem.MolFromMol2File(C.MOL2, removeHs=True)
    if mol_mmff is None or mol_mmff.GetNumAtoms() != mol.GetNumAtoms():
        log(lineas, "   no se pudo leer el mol2 aromatizado: sin medida de tension")
        mol_mmff = None
    e_cristal = e_min_cristal = None
    if mol_mmff is not None:
        e_cristal, e_min_cristal = energia_mmff(mol_mmff, ref)
    if e_cristal is None:
        log(lineas, "   MMFF no tiene parametros para esta molecula")
    else:
        log(lineas, "   conformero del cristal: %.2f kcal/mol (su minimo local"
            " %.2f -> tension %.2f)"
            % (e_cristal, e_min_cristal, e_cristal - e_min_cristal))
        log(lineas, "   (la tension del cristal se mide con el MISMO protocolo que"
            " la de las poses:")
        log(lineas, "    atomos pesados fijos, hidrogenos regenerados y sueltos)")
        for etiqueta, mejor, pos, _ in resumen:
            e_pose, e_min = energia_mmff(mol_mmff, pos)
            if e_pose is None:
                continue
            log(lineas, "   %-16s mejor pose: %.2f kcal/mol (tension %.2f)"
                % (etiqueta, e_pose, e_pose - e_min))
            if e_cristal <= e_pose:
                log(lineas, "      el cristal NO esta mas tenso que la pose del"
                    " motor (%.2f contra %.2f): lo" % (e_cristal, e_pose))
                log(lineas, "      descarto la PUNTUACION, no la tension.")
            else:
                log(lineas, "      el cristal esta %+.2f kcal/mol MAS tenso que lo"
                    " que encuentra" % (e_cristal - e_pose))
                log(lineas, "      el motor: puede estar evitando una conformacion"
                    " forzada, y el problema")
                log(lineas, "      seria la estructura de partida.")

    # --- lo que queda por probar, que necesita la tarjeta ---
    log(lineas, "")
    log(lineas, "SIGUIENTE PRUEBA (necesita la GPU, que esta con TBK1): caja"
        " estrecha")
    log(lineas, "   Si se repite el control con una caja de 6-8 A pegada a la pose"
        " del cristal y el")
    log(lineas, "   motor la reproduce, lo que fallaba era la BUSQUEDA. Si ni asi"
        " la saca, es la")
    log(lineas, "   PUNTUACION, y lo que hay que cambiar es el motor (Vinardo, o"
        " rescoring por")
    log(lineas, "   interacciones encima de las poses). Se lanza con:")
    log(lineas, "      python control_redocking_cdk2.py --caja 6 --depth 64")
    log(lineas, "      python control_redocking_cdk2.py --caja 8 --depth 64")

    with open(SALIDA, "w", encoding="utf-8") as f:
        f.write("\n".join(lineas) + "\n")
    print("\nguardado: %s" % SALIDA)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
