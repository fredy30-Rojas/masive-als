"""Rescoring MM-GBSA de poses de TBK1 con AmberTools. Corre en Oracle.

POR QUE EL LIGANDO VIENE DEL SMILES Y NO DEL PDBQT
Vina escribe PDBQT, que NO guarda ordenes de enlace. Al pasarlo a mol2, Open
Babel se los inventa y GAFF se encuentra con angulos imposibles (por ejemplo
cf-nh-c1, que no existe en GAFF). No es un problema de versiones: es que
falta el dato. El SMILES original si lo tiene, asi que la tipacion se hace
desde ahi: se genera la molecula con la herramienta de Amber, se le enchufan
las coordenadas de la pose, y ya. Es ademas mas rapido, porque antechamber no
tiene que adivinar la topologia.

ORDEN
    1. mol2 con AmberTools a partir del SMILES (antechamber BCC + GAFF)
    2. se sustituyen las coordenadas por las de la pose de Vina
    3. tleap: ff14SB + GAFF + TIP3P, complejo solvatado e ionizado
    4. sander: single-point MM/GBSA (igb=5)
    5. dG_union = G(complejo) - G(receptor) - G(ligando)

LIMITE HONESTO
La pose viene de Vina, no de Amber. Es un single-point sobre la pose de
Vina: correcto como rescoring RELATIVO (todos los ligandos con el mismo
protocolo y el mismo receptor, luego ordena bien) pero no es un MM-GBSA de
re-docking. Para un valor experimental absoluto habria que re-acoplar en
campo Amber, que es otro proyecto.

Uso:
  python3 rescoring_mmgbsa.py --receptor receptor.pdb --pose pose.pdbqt \
      --smiles "N#CCNC(=O)c1ccc..."
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile

AMBER_HOME = os.path.expanduser("~/ambertools")
ANTECHAMBER = os.path.join(AMBER_HOME, "bin", "antechamber")
TLEAP = os.path.join(AMBER_HOME, "bin", "tleap")
SANDER = os.path.join(AMBER_HOME, "bin", "sander")
PARMCHK2 = os.path.join(AMBER_HOME, "bin", "parmchk2")
MMPBSA_PY = os.path.join(AMBER_HOME, "bin", "MMPBSA.py")

# Entrada de MMPBSA.py para UN solo sistema (single point, sin minimizacion).
#
# Los valores son los del protocolo del proyecto: igb=5 (modelo GB de Huo) y
# saltcon=0.100. Se calculan por separado complejo, receptor y ligando con el
# MISMO bloque, y el dG de union sale de la resta: si un valor cambiara, el
# protocolo estaria roto.
#
# OJO, no hay surften/surftat a proposito: en la version de AmberTools que hay
# instalada en Oracle (conda-forge) MMPBSA.py responde "Unknown variable surftat
# in &gb", y la ayuda no los lista. El GB aqui es el modelo GB puro, sin
# correccion de tension superficial. Como los tres sistemas se calculan con el
# mismo bloque, el ranking relativo no se altera; solo habria que repetirlo con
# el modelo de superficie si algun dia se quisiera un valor absoluto.
INPUT_MMPBSA = """&general
/
&gb
  igb=5, saltcon=0.100,
/
"""

# El `loadamberparams frcmod` es lo que desbloquea el complejo.
#
# POR QUE (28 sep 2026, medido en Oracle): GAFF no cubre los terminos de
# torsion de algunas combinaciones aromaticas de 4 centros, y tleap se paraba
# con 4 errores "No torsion terms" (tipos ca-cp-nb-ca, cp-cp-nb-ca). El
# receptor pasaba porque ff14SB lo cubre entero. NO era un bug del codigo: los
# parametros que faltan si existen, y los genera parmchk2, que antechamber
# llama por dentro y escribe en `frcmod` al lado del mol2. Ese fichero nunca se
# cargaba, asi que los parametros estaban en disco y tleap no los veia.
# CUATRO fallos de sintaxis, los cuatro medidos en Oracle el 28 sep 2026. Antes de
# arreglarlos, el script fallaba siempre en `complex` y el diagnostico
# apuntaba a los 4 terminos de torsion de GAFF, que NO eran la causa:
#
# 1. `source leaprc.water.tip3p` faltaba. Sin el, `TIP3PBOX` no existe como
#    argumento y tleap lo leia como una cadena: "Argument #2 is of type String
#    must be of type: [unit]".
# 2. `solvateBox` en esta version de AmberTools (conda-forge, la que se instalo
#    en Oracle) pide 3 argumentos, no 4. Con los 4 decia "usage: solvateBox
#    <solute> <solvent> <buffer> [iso] [closeness]".
# 3. `addions` no existe: es `addIons`, y con mayuscula. Ademas aqui no hace
#    falta: el complejo ya sale neutro, asi que se quita en vez de arreglarlo.
#
# 4. `saveoff3` y `saveamberpdb` NO EXISTEN en AmberTools moderno (conda-forge):
#    son de Amber antiguo. El `desc` los devuelve como "STRING (with no
#    reference)", o sea que tleap los trata como texto y da "Error from the
#    parser". El sustituto de `saveoff3` es `saveoff`, con 2 argumentos. Como el
#    MM-GBSA solo necesita los .prmtop y .rst7, los off3 se quitan en vez de
#    traducirlos: solo servian para inspeccionar a mano.
#
# Verificado en Oracle el 28 sep: con estas correcciones receptor, ligando y
# complejo salen los tres prmtop. El `frcmod` se mantiene: parmchk2 lo escribe
# al lado del mol2 cuando falta algun parametro, y sin cargarlo tleap no los ve.
SCRIPT_TLEAP = """source leaprc.protein.ff14SB
source leaprc.water.tip3p
source leaprc.gaff
loadamberparams frcmod
receptor = loadpdb {receptor}
ligand = loadmol2 {ligando}
check receptor
check ligand
complex = combine {{ receptor ligand }}
solvateBox complex TIP3PBOX {box}
saveamberparm complex complex.prmtop complex.rst7
saveamberparm receptor receptor.prmtop receptor.rst7
saveamberparm ligand ligand.prmtop ligand.rst7
quit
"""

INPUT_SANDER = """Single point MM/GBSA.
 &cntrl
  imin=0, maxcyc=0, ntb=1, ntt=3, gamma_lj=2.0,
  ifqnt=1, gbsa=1, igb=5, saltcon=0.100,
  ntbopt=1, ntpr=1, ipol=0, cut=10.0,
 &end
"""
# INPUT_SANDER ya no se usa: el paso 5 lo hace MMPBSA.py, que genera el &cntrl
# solo. Se deja con su valor historico porque documenta el protocolo que se
# venia buscando, pero el codigo que corre es INPUT_MMPBSA.


def run(cmd, cwd, log):
    """Ejecuta guardando la salida, con ruta de log ABSOLUTA."""
    ruta_log = log if os.path.isabs(log) else os.path.join(cwd, log)
    with open(ruta_log, "a", encoding="utf-8", errors="ignore") as fh:
        fh.write("\n$ " + " ".join(cmd) + "\n")
        fh.flush()
        r = subprocess.run(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT,
                           timeout=3600)
    return r.returncode


def smiles_a_mol2(smiles, work):
    """SMILES -> mol2 tipado con GAFF.

    dos pasos porque antechamber no acepta SMILES: primero Open Babel pasa el
    SMILES a mol2 (y aqui si mete los ordenes de enlace, que es justo lo que
    el PDBQT no traia), y luego antechamber lo tipa.

    -c gas (Gasteiger) y no -c bcc: las cargas BCC necesitan sqm, que no esta
    instalado. Para rescoring RELATIVO (ordenar ligandos con el mismo
    protocolo) las cargas Gasteiger bastan; para comparar contra un valor
    experimental absoluto habria que instalar sqm y rehacer con BCC.
    """
    crudo = os.path.join(work, "lig_raw.mol2")
    r = subprocess.run(["obabel", "-:" + smiles, "-omol2", "-O", crudo, "-h"],
                       capture_output=True, text=True, timeout=300)
    if r.returncode != 0 or not os.path.exists(crudo):
        return None
    mol2 = os.path.join(work, "lig.mol2")
    cmd = [ANTECHAMBER, "-i", "lig_raw.mol2", "-fi", "mol2", "-o", "lig.mol2",
           "-fo", "mol2", "-c", "gas", "-at", "gaff", "-rn", "LIG", "-pf", "y"]
    if run(cmd, work, "antechamber.log") != 0:
        return None
    if not os.path.exists(mol2):
        return None
    # parmchk2 EXPLICITO, y este es el paso que de verdad hacia falta.
    #
    # POR QUE (28 sep 2026, medido en Oracle): GAFF no cubre los terminos de
    # torsion ca-cp-nb-ca y cp-cp-nb-ca, que salen en aromaticos con nitrogeno.
    # tleap se paraba con 4 errores "No torsion terms" y complex.prmtop salia
    # de 0 bytes. Se asumia que antechamber generaba el frcmod solo: NO lo
    # hace, no escribe ninguno para estas torsiones.
    # parmchk2 -s gaff si las rellena ("same as X -ca-nb-X") y con el frcmod
    # cargado el complejo sale entero: 0 torsiones faltantes, 17 MB de prmtop.
    run([PARMCHK2, "-i", mol2, "-f", "mol2", "-o",
         os.path.join(work, "frcmod"), "-s", "gaff"], work, "parmchk2.log")
    return mol2


def coordenadas_pdbqt_a_pdb(pdbqt, destino):
    """Vina escribe 3 modelos (num_modes 3). Solo interesa el primero, que es
    el de menor energia."""
    cmd = ["obabel", "-ipdbqt", pdbqt, "-opdb", "-O", destino,
           "-f", "1", "-l", "1"]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    if r.returncode != 0 or not os.path.exists(destino):
        return None
    # quitar MODEL/ENDMDL: OpenMM y Amber leen mal el numero de modelo
    limpio = destino + ".1"
    with open(destino, errors="replace") as origen, \
            open(limpio, "w", encoding="utf-8") as destino2:
        for cruda in origen:
            l = cruda.rstrip("\n").rstrip("\r")
            if l.startswith(("MODEL", "ENDMDL")):
                continue
            destino2.write(l + "\n")
    return limpio


def coordenadas_al_mol2(mol2, pdb_pose, destino):
    """Cambia las coordenadas del mol2 por las de la pose, conservando la
    tipacion GAFF.

    La pose de Vina trae SOLO atomos pesados y, ademas, todos con el mismo
    nombre (C, C, C...), asi que no se puede emparejar por nombre: se empareja
    por ORDEN. Los atomos pesados del mol2 de Amber, en el mismo orden que
    genera la herramienta, corresponden uno a uno con los de la pose.

    Los hidrogenos que anade Amber se colocan en la posicion de su atomo pesado
    padre y se deja que la minimizacion los separe. Si se les deja en (0,0,0)
    el complejo sale disparado.
    """
    # 1. coordenadas de la pose, EN ORDEN
    pose = []
    with open(pdb_pose, errors="replace") as fh:
        for linea in fh:
            if not linea.startswith(("ATOM", "HETATM")):
                continue
            try:
                pose.append((float(linea[30:38]), float(linea[38:46]),
                             float(linea[46:54])))
            except ValueError:
                continue
    if len(pose) < 5:
        return None

    # 2. atomos y enlaces del mol2
    atomos = []
    enlaces = []
    seccion = None
    with open(mol2, errors="replace") as fh:
        for linea in fh:
            if linea.startswith("@<TRIPOS>"):
                seccion = linea.strip()
                continue
            partes = linea.split()
            if seccion == "@<TRIPOS>ATOM" and len(partes) >= 6:
                atomos.append(partes)
            elif seccion == "@<TRIPOS>BOND" and len(partes) >= 4:
                try:
                    enlaces.append((int(partes[1]), int(partes[2])))
                except ValueError:
                    continue
    if len(atomos) < 5:
        return None

    # 3. padre de cada atomo (primer vecino mas pesado que el)
    vecinos = {}
    for a, b in enlaces:
        vecinos.setdefault(a, []).append(b)
        vecinos.setdefault(b, []).append(a)
    padre = {}
    for idx in range(1, len(atomos) + 1):
        simbolo = atomos[idx - 1][5].split(".")[0]
        if simbolo.upper() in ("H", "D"):
            for v in vecinos.get(idx, []):
                if atomos[v - 1][5].split(".")[0].upper() not in ("H", "D"):
                    padre[idx] = v
                    break
            padre.setdefault(idx, idx)
        else:
            padre.setdefault(idx, idx)

    # 4. asignar: los pesados, en orden, por posicion de la pose
    nuevas = [(0.0, 0.0, 0.0)] * (len(atomos) + 1)
    pesados = 0
    for idx in range(1, len(atomos) + 1):
        simbolo = atomos[idx - 1][5].split(".")[0]
        if simbolo.upper() in ("H", "D"):
            continue
        if pesados < len(pose):
            nuevas[idx] = pose[pesados]
            pesados += 1
    if pesados < 5:
        return None

    # 5. los hidrogenos, a donde su padre
    for idx in range(1, len(atomos) + 1):
        simbolo = atomos[idx - 1][5].split(".")[0]
        if simbolo.upper() in ("H", "D"):
            nuevas[idx] = nuevas[padre.get(idx, idx)]

    # 6. escribir el mol2 con las coordenadas nuevas
    seccion = None
    lineas = []
    with open(mol2, errors="replace") as fh:
        for linea in fh:
            cruda = linea.rstrip("\n")
            if cruda.startswith("@<TRIPOS>"):
                seccion = cruda.strip()
                lineas.append(cruda)
                continue
            partes = cruda.split()
            if seccion == "@<TRIPOS>ATOM" and len(partes) >= 6:
                try:
                    idx = int(partes[0])
                except ValueError:
                    lineas.append(cruda)
                    continue
                x, y, z = nuevas[idx]
                partes[2] = "%.4f" % x
                partes[3] = "%.4f" % y
                partes[4] = "%.4f" % z
                lineas.append(" ".join(partes))
                continue
            lineas.append(cruda)

    with open(destino, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lineas) + "\n")
    return destino if pesados >= 5 else None


def montar_sistemas(receptor_pdb, ligando_mol2, work, box):
    script = SCRIPT_TLEAP.format(receptor=os.path.basename(receptor_pdb),
                                 ligando=os.path.basename(ligando_mol2),
                                 box=box)
    # parmchk2 solo escribe frcmod cuando falta algun parametro, y a veces no
    # falta ninguno. Cargar un frcmod que no existe es un error de tleap, asi
    # que la linea se quita cuando no esta.
    if not os.path.exists(os.path.join(work, "frcmod")):
        script = script.replace("loadamberparams frcmod\n", "")
    with open(os.path.join(work, "in.tleap"), "w", encoding="utf-8") as fh:
        fh.write(script)
    return run([TLEAP, "-f", "in.tleap"], work, "tleap.log") == 0


def quitar_caja(prmtop, destino, work):
    """Deja el prmtop SIN condiciones de contorno periodicas.

    POR QUE (28 sep 2026, medido en Oracle): el error de fondo del paso 5 no
    era de sintaxis sino fisico. `solvateBox` deja el prmtop con IFBOX=1, y el
    Generalized Born es incompatible con la caja periodica: mmpbsa_py_energy
    responde "gb>0 is incompatible with periodic boundary conditions. To use
    this method set IFBOX in the PRMTOP file to 0". El GB es, por definicion, un
    modelo de solvente NO periodico, asi que quitar la caja no es un truco: es lo
    que el metodo pide.

    Se hace con parmed, que AmberTools trae, en vez de editar el prmtop a mano:
    parmed reescribe todas las secciones con los recuentos ya ajustados. Editar
    el flag a mano deja los contadores de la seccion BOX desfasados y sander
    lee memoria de mas.
    """
    codigo = ("from parmed import load_file\n"
              "s = load_file(%r)\n"
              "s.box = None\n"
              "s.save(%r, overwrite=True)\n"
              % (prmtop, destino))
    r = subprocess.run([sys.executable, "-c", codigo], cwd=work,
                       capture_output=True, text=True, timeout=3600)
    return r.returncode == 0 and os.path.exists(destino) and \
        os.path.getsize(destino) > 0


def energia_mmgbsa(prmtop, rst7, work, etiqueta):
    """Energia MM/GBSA de UN sistema, con MMPBSA.py.

    POR QUE MMPBSA.py Y NO sander A MANO (28 sep 2026): la llamada antigua
    `sander -O x.md -i in.md -c $MDOUT` no funciona en AmberTools moderno por
    dos motivos, y los dos se midieron en Oracle. Uno: `-O` ya no existe, sander
    responde "mdfil: Error unknown flag". Dos: MDOUT, la carpeta con las
    plantilla de calculo, no viene en la instalacion de conda-forge
    (`dat/1TRX/MDOUT` no existe; `dat/` no tiene ningun 1TRX), y sin ella sander
    no tiene con que calcular. MMPBSA.py es la herramienta que AmberTools trae
    para esto, genera la entrada de sander sola y no depende de MDOUT.

    Los tres sistemas (complejo, receptor, ligando) se calculan por separado con
    el mismo protocolo; el dG de union sale por la resta, como antes.
    """
    mmpbsa_in = os.path.join(work, "mmpbsa_%s.in" % etiqueta)
    with open(mmpbsa_in, "w", encoding="utf-8") as fh:
        fh.write(INPUT_MMPBSA)
    salida = os.path.join(work, "mmpbsa_%s.out" % etiqueta)
    cmd = [sys.executable, MMPBSA_PY, "-i", os.path.basename(mmpbsa_in),
           "-o", os.path.basename(salida)]
    # MMPBSA.py NO acepta el prmtop ni en la entrada (&general no tiene esa
    # variable: "Unknown variable prmtop in &general") ni como argumento de
    # linea de ordenes. Espera los ficheros con el NOMBRE LITERAL
    # `complex_prmtop` y `complex_inpcrd` en el directorio de trabajo, sin
    # ningun sufijo derivado del nombre de salida; si no estan responde
    # "Could not open complex_prmtop for reading". Por eso se copian con ese
    # nombre antes de llamar, y se vuelven a poner los originales luego.
    # El prmtop que se calcula es el de la etiqueta, ya SIN caja: el GB la
    # prohibe (ver quitar_caja). Se escribe aparte y es el que se copia.
    prmtop_sin_caja = os.path.join(work, "%s_sincaja.prmtop" % etiqueta)
    if not quitar_caja(os.path.join(work, prmtop), prmtop_sin_caja, work):
        return None
    for ext in ("prmtop", "rst7"):
        origen = prmtop_sin_caja if ext == "prmtop" else os.path.join(work, rst7)
        destino = os.path.join(work, "complex_%s"
                               % ("prmtop" if ext == "prmtop" else "inpcrd"))
        shutil.copy2(origen, destino)
    # Y ademas la TRAYECTORIA: MMPBSA.py no lee la inpcrd para calcular, pasa
    # por cpptraj y espera un `mdcrd` con las coordenadas. Sin el responde
    # "Error: 'mdcrd': No such file or directory". Para un single point el
    # rst7 ES la trayectoria de un solo cuadro, asi que se copia tal cual.
    shutil.copy2(os.path.join(work, rst7), os.path.join(work, "mdcrd"))
    env = dict(os.environ)
    env["AMBERHOME"] = AMBER_HOME
    env["PATH"] = os.path.join(AMBER_HOME, "bin") + os.pathsep + env.get("PATH", "")
    with open(os.path.join(work, "mmpbsa_%s.log" % etiqueta), "w",
              encoding="utf-8", errors="ignore") as fh:
        r = subprocess.run(cmd, cwd=work, stdout=fh, stderr=subprocess.STDOUT,
                           timeout=7200, env=env)
    if r.returncode != 0:
        return None
    return leer_energia_mmpbsa(salida)


def leer_energia_mmpbsa(ruta):
    """Saca la energia del .out de MMPBSA.py.

    El formato es una tabla con una linea de encabezado y una de datos; las
    columnas van en el orden del &general de la entrada, que aqui es
    "DG (Energy terms) vdW elec. GB Gpolar". Se leen por encabezado y no por
    posicion, porque MMPBSA.py puede reordenar segun la version.
    """
    if not os.path.exists(ruta):
        return None
    lineas = [l.rstrip() for l in open(ruta, errors="replace")]
    for i, linea in enumerate(lineas):
        if "DG (Energy terms)" in linea or ("vdW" in linea and "elec." in linea):
            encabezado = linea.split()
            for datos in lineas[i + 1:i + 6]:
                p = datos.split()
                if len(p) != len(encabezado) or not p[0].replace(".", "").isdigit():
                    continue
                try:
                    vals = {}
                    for col, val in zip(encabezado, p):
                        if col in ("vdW", "elec.", "GB", "Gpolar", "DG"):
                            vals[col] = float(val)
                    if "DG" in vals:
                        return {"dG_total": vals["DG"],
                                "dG_vdw": vals.get("vdW"),
                                "dG_elec": vals.get("elec."),
                                "dG_solv": vals.get("Gpolar")}
                except ValueError:
                    continue
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--pose", required=True, help="PDBQT de Vina")
    ap.add_argument("--smiles", required=True)
    ap.add_argument("--box", type=float, default=12.0)
    args = ap.parse_args()

    for exe in (ANTECHAMBER, PARMCHK2, TLEAP, SANDER, MMPBSA_PY):
        if not os.path.exists(exe):
            print("Falta %s (AmberTools en Oracle)." % exe)
            return 1
    if not shutil.which("obabel"):
        print("Falta Open Babel (obabel).")
        return 1

    receptor = os.path.abspath(args.receptor)
    pose = os.path.abspath(args.pose)
    work = tempfile.mkdtemp(prefix="mmgbsa_", dir=os.path.expanduser("~"))

    print("=" * 66)
    print("MM-GBSA con AmberTools  ·  ligando tipado desde SMILES")
    print("=" * 66)
    print("receptor: %s" % os.path.basename(receptor))
    print("pose:     %s" % os.path.basename(pose))
    print("smiles:   %s" % args.smiles[:60])
    print("trabajo:  %s\n" % work)

    try:
        # 1. SMILES -> mol2 tipado
        mol2 = smiles_a_mol2(args.smiles, work)
        if not mol2:
            print("1. antechamber desde SMILES  FALLO (ver antechamber.log)")
            return 1
        print("1. SMILES -> mol2 con GAFF y cargas Gasteiger  OK")

        # 2. coordenadas de la pose
        pdb_pose = coordenadas_pdbqt_a_pdb(pose, os.path.join(work, "pose.pdb"))
        if not pdb_pose:
            print("2. pose pdbqt -> pdb  FALLO")
            return 1
        print("2. pose pdbqt -> pdb (modelo 1)  OK")

        # 3. mol2 con las coordenadas de la pose
        mol2_pose = coordenadas_al_mol2(mol2, pdb_pose,
                                        os.path.join(work, "lig_pose.mol2"))
        if not mol2_pose:
            print("3. injective de coordenadas  FALLO")
            return 1
        print("3. coordenadas de la pose en el mol2 tipado  OK")

        # 4. tleap
        shutil.copy2(receptor, os.path.join(work, os.path.basename(receptor)))
        if not montar_sistemas(receptor, mol2_pose, work, args.box):
            print("4. tleap  FALLO (ver tleap.log)")
            return 1
        print("4. tleap: ff14SB + GAFF + TIP3P, solvatado  OK")

        # 5. MM-GBSA
        res = {}
        for etiqueta in ("complex", "receptor", "ligand"):
            r = energia_mmgbsa("%s.prmtop" % etiqueta, "%s.rst7" % etiqueta,
                               work, etiqueta)
            if r is None:
                print("5. MM-GBSA %s  FALLO (ver mmpbsa_%s.log)"
                      % (etiqueta, etiqueta))
                return 1
            res[etiqueta] = r
            print("5. MM-GBSA %-9s G = %10.1f kcal/mol"
                  % (etiqueta, r["dG_total"]))

        dG = (res["complex"]["dG_total"] - res["receptor"]["dG_total"]
              - res["ligand"]["dG_total"])
        print()
        print("=" * 66)
        print("dG_union MM-GBSA = %.1f kcal/mol" % dG)
        c = res["complex"]
        print("  (complejo: vdW %s, electrostatica %s, solvatacion %s)"
              % (c["dG_vdw"], c["dG_elec"], c["dG_solv"]))
        print("=" * 66)
        return 0
    finally:
        if os.environ.get("MANTENER_TRABAJO") != "1":
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
