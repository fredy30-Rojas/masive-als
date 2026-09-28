"""Rescoring MM-GBSA de poses de TBK1 con AmberTools. Corre en Oracle.

EL LIGANDO SE TIPEA DESDE LA POSE, Y ESO LO ARREGLA TODO (28 sep 2026)
Durante dias se tipo el ligando desde el SMILES con antechamber y luego se
metieron las coordenadas de la pose dentro, Emparejando UNO CON UNO POR
ORDEN. Ese emparejamiento era FALSO, y de ahi salian todos los fallos que
venian uno detras de otro (VDWAALS nan, EEL inf, y al final un H a 0,546 A
del C que tenia al lado).

POR QUE ERA FALSO: antechamber REORDENA los atomos. La pose de Vina viene
en el orden que le dio Open Babel al SMILES; el mol2 tipado sale en el orden
interno de Amber. Los dos son el mismo numero de atomos pero NO el mismo
atomo en la misma posicion. Medido el 28 sep 2026: en el mol2 el C4 y el C9
estan a cuatro enlaces de distancia, y al meter las coordenadas de la pose
salen a 1,489 A, que es un enlace y medio. O sea, el C de la pose se le
ponia la carga y el tipo del C equivocado. Todo lo que se vio despues
(solapamientos, tipos de H que no existian, angulos imposibles) eran
SINTOMAS de este emparejamiento erroneo, no fallos aparte.

LA SOLUCION ES INVERTIR EL ORDEN, no arreglar el emparejamiento: que sea
AmberTools el que reordene, sobre unas coordenadas que ya son las suyas. Se
arma el mol2 A PARTIR de la pose, se le ponen los hidrogenos que falten, y
se le pasa a antechamber, que escribe el mol2 GAFF con las coordenadas que
le ha llegado. Cada mol2 lleva ya sus coordenadas y no hay ningun
emparejamiento que pueda fallar. Ademas el tipado GAFF sale por CONTEXTO
(real) y no por elemento, que era el otro defecto del metodo anterior.

ORDEN
    1. pdbqt de la pose -> mol2 con enlaces (Open Babel, que si los lee)
    2. se le ponen los hidrogenos que falten, con geometria ideal
    3. antechamber + parmchk2: tipa con GAFF y escribe el frcmod
    4. tleap: ff14SB + GAFF, SIN agua (con igb=5 el solvente es
       implicito). Ver el comentario de SCRIPT_TLEAP.
    5. MMPBSA.py: single-point MM/GBSA (igb=5)
    6. dG_union = G(complejo) - G(receptor) - G(ligando)

LIMITE HONESTO
La pose viene de Vina, no de Amber. Es un single-point sobre la pose de
Vina: correcto como rescoring RELATIVO (todos los ligandos con el mismo
protocolo y el mismo receptor, luego ordena bien) pero no es un MM-GBSA de
re-docking. Para un valor experimental absoluto habria que re-acoplar en
campo Amber, que es otro proyecto. Ademas los H que faltan los pone una
regla de geometria ideal, no el optimizador de Amber; sesga el valor
absoluto, no el RANKING, que es lo unico que se usa aqui.

Uso:
  python3 rescoring_mmgbsa.py --receptor receptor.pdb --pose pose_out.pdbqt
"""
import argparse
import math
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
# parmed (lo usa quitar_caja) NO esta en el python3 del sistema: vive en el
# python de AmberTools. Con `sys.executable` daba "ModuleNotFoundError: No
# module named 'parmed'" y el paso 5 se caia sin decir nada, porque el fallo
# pasaba antes de abrir el log de MMPBSA.py. Se usa el de AmberTools.
PYTHON_AMBER = os.path.join(AMBER_HOME, "bin", "python")

# Y MMPBSA.py necesita lo mismo con su propio modulo: sin esto responde
# "ImportError: Could not import Amber Python modules. Please make sure you
# have sourced .../amber.sh" y "No module named 'MMPBSA_mods'". Ese PYTHONPATH
# solo lo pone amber.sh, que es un script de shell, y este es un proceso
# python: se pone a mano desde env_entorno_ambertools(). La version de python
# no se fija a mano, se lee de la propia AmberTools, que aqui es 3.11 pero
# podria no serlo.
def site_packages_ambertools():
    import glob
    encontrados = sorted(glob.glob(os.path.join(
        AMBER_HOME, "lib", "python3.*", "site-packages")))
    return os.pathsep.join(encontrados)


def env_entorno_ambertools():
    """Entorno con AmberTools en PATH, PYTHONPATH y LD_LIBRARY_PATH."""
    env = dict(os.environ)
    env["AMBERHOME"] = AMBER_HOME
    env["PATH"] = os.path.join(AMBER_HOME, "bin") + os.pathsep + env.get("PATH", "")
    libs = site_packages_ambertools()
    if libs:
        env["PYTHONPATH"] = libs + os.pathsep + env.get("PYTHONPATH", "")
    lib = os.path.join(AMBER_HOME, "lib")
    if os.path.isdir(lib):
        env["LD_LIBRARY_PATH"] = lib + os.pathsep + env.get("LD_LIBRARY_PATH", "")
    return env

# Grosor de la capa de agua EXPLICITA, en angstroms. Por defecto NO se pone
# ninguna: con igb=5 el solvente es implicito y el agua explicita solo estorba
# (ver SCRIPT_TLEAP). Se deja elparametro porque --capa sigue existiendo, por
# si algun dia hace falta un single-point con solvente real.

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
#
# OTROS TRES FALLOS DE SINTAXIS, los tres medidos en Oracle el mismo dia:
#
# 1. `addions` no existe: es `addIons`, y con mayuscula. Ademas aqui no hace
#    falta: el complejo ya sale neutro, asi que se quita en vez de arreglarlo.
# 2. `saveoff3` y `saveamberpdb` NO EXISTEN en AmberTools moderno (conda-forge):
#    son de Amber antiguo. El `desc` los devuelve como "STRING (with no
#    reference)", o sea que tleap los trata como texto y da "Error from the
#    parser". El sustituto de `saveoff3` es `saveoff`, con 2 argumentos. Como el
#    MM-GBSA solo necesita los .prmtop y los .rst7, los off3 se quitan en vez
#    de traducirlos: solo servian para inspeccionar a mano.
# 3. `sander -O` no existe: responde "mdfil: Error unknown flag". Y MDOUT, la
#    carpeta con la plantilla de calculo, no viene en la instalacion de
#    conda-forge. Por eso el paso de energia usa MMPBSA.py, que genera la
#    entrada de sander solo.
#
# POR QUE NO HAY AGUA (28 sep 2026, medido en Oracle).
#
# Se probo primero con solvateOct y con solvateBox, y las dos se quitaron por
# motivos distintos y medidos:
#
# 1. solvateBox: 90.570 atomos y el MM-GBSA se colgaba. No era falta de RAM
#    (la maquina tiene 23 GiB), era TAMANO.
#
# 2. solvateOct con 12 A: 15.426 atomos y el MM-GBSA ya corria en 30 s, pero
#    la energia salia con vdW de 1,98e15 kcal/mol. La causa no era el
#    ligando: `solvateOct` pone agua alrededor de lo que le digas, sin
#    preguntar si hay proteina alrededor. Metia 5.685 moleculas de agua
#    DENTRO del receptor, y el par mas cercano quedaba a 0,28 A. Medido:
#    OG 409 contra 14532 a 0,684 A, HD22 462 contra 14190 a 0,385 A, y
#    muchos mas. El motor no protesta: devuelve un numero enorme y sigue.
#
# 3. Y ademas el agua no hacia falta. Con igb=5 el solvente es IMPLICITO:
#    lo calcula el propio Generalized Born a partir de las cargas y los
#    radios atomicos. El agua explicita no aporta nada al dG y, al estar
#    dentro del sistema, GB la trata como si fuera soluto: solo puede
#    ensuciar el resultado. Un MM-GBSA con GB no lleva agua; lo que se
#    compara en la literatura son las energias de un solo cuadro, sin MD.
#
# Asi que el sistema son receptor + ligando, y el GB hace su trabajo. Si
# alguna vez hace falta agua explicita, se pasa --capa 12 y hay que quitar
# despues el agua que ha caído dentro de la proteina, que tleap no hace.
SCRIPT_TLEAP = """source leaprc.protein.ff14SB
source leaprc.gaff
loadamberparams frcmod
receptor = loadpdb {receptor}
ligand = loadmol2 {ligando}
check receptor
check ligand
complex = combine {{ receptor ligand }}
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


def pose_a_mol2(pdbqt, work, receptor_pdb=None):
    """PDBQT de la pose -> mol2 tipado con GAFF.

    EL ORDEN ESTA INVERTIDO RESPECTO AL QUE SE INTENTO ANTES (28 sep 2026).
    Antes: mol2 desde el SMILES con antechamber, y luego se metian las
    coordenadas de la pose Emparejando POR ORDEN. Ese emparejamiento era
    falso, porque antechamber reordena los atomos: el C4 y el C9 estan a
    cuatro enlaces en el mol2, y al poner las coordenadas de la pose salen a
    1,489 A. El C de la pose se quedaba con el tipo y la carga del C
    equivocado, y de ahi salian los nan, los inf y los H a 0,5 A del vecino.
    Ver la cabecera del fichero.

    Ahora: la pose es la que manda. Open Babel lee el pdbqt de Vina, que si
    trae la conectividad (Vina escribe las ramas ROOT/BRANCH), y de ahi sale
    un mol2 CON LA POSE DENTRO. Se le anaden los hidrogenos que falten con
    geometria ideal, se lo pasa a antechamber para que lo tipee, y antechamber
    escribe su mol2 con las coordenadas que le ha llegado. Cada mol2 lleva
    sus propias coordenadas y no hay ningun emparejamiento por orden que
    pueda fallar. De regalo, GAFF se asigna por CONTEXTO real y no por
    elemento.

    -f 1 -l 1: el pdbqt de Vina trae 3 modelos y aqui solo interesa el
    primero, que es el de menor energia.

    -c gas (Gasteiger) y no -c bcc: las cargas BCC necesitan sqm, que no esta
    instalado. Para rescoring RELATIVO (ordenar ligandos con el mismo
    protocolo) las cargas Gasteiger bastan; para comparar contra un valor
    experimental absoluto habria que instalar sqm y rehacer con BCC.
    """
    crudo = os.path.join(work, "pose_h.mol2")
    r = subprocess.run(["obabel", "-ipdbqt", pdbqt, "-omol2", "-O", crudo,
                        "-f", "1", "-l", "1"],
                       capture_output=True, text=True, timeout=300)
    if r.returncode != 0 or not os.path.exists(crudo):
        return None
    if not anadir_hidrogenos(crudo, crudo + ".h", receptor_pdb):
        return None
    mol2 = os.path.join(work, "lig.mol2")
    cmd = [ANTECHAMBER, "-i", os.path.basename(crudo) + ".h", "-fi", "mol2",
           "-o", "lig.mol2", "-fo", "mol2", "-c", "gas", "-at", "gaff",
           "-rn", "LIG", "-pf", "y"]
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


# FORMATO DE COLUMNAS DE UN MOL2 DE AMBERTOOLS (medido el 28 sep 2026 con
# cat -A sobre el mol2 que escribe antechamber, NO al ojo):
#
#   atomos:   id 0-6, espacio, nombre 8-15, x 16-25, y 26-35, z 36-45,
#             espacio 46, tipo 47-52, subst_id 53-55, subst_name 56-61,
#             carga 62-75
#   enlaces:  id 0-5, a1 6-11, a2 12-17, tipo 18-22
#   cabecera: @<TRIPOS>MOLECULE lleva los RECUENTOS en la tercera linea,
#             "%5d%5d%5d%5d%5d" con (atomos, enlaces, ...). Si no se actualizan
#             al quitar hidrogenos, tleap responde "Invalid MOL2 format.
#             Incorrect number of items in atom record" y "Last line read:
#             @<TRIPOS>BOND": el error senala la linea de enlaces pero la causa
#             esta en la cabecera.
#
# tleap lee por COLUMNAS, no por espacios: unir con " ".join descuadra todo.
LINEA_ATOM_MOL2 = ("%7d %-8s%10.4f%10.4f%10.4f %-6s%3s%6s%14.4f")
LINEA_ENLACE_MOL2 = "%6d%6d%6d%5s"
CABECERA_MOL2 = "%5d%5d%5d%5d%5d"

# Longitudes de enlace estandar para colocar los H (angstroms). Son las de
# GAFF/ff14SB; se usan solo para dar geometria a los H que AmberTools no sabe
# colocar en esta build (ver anadir_hidrogenos).
_LONGITUD_ENLADE_H = {"C": 1.090, "N": 1.010, "O": 0.970, "S": 1.340, "P": 1.840}

# Para saber CUANTOS hidrogenos le faltan a cada atomo pesado se usa la
# VALENCIA, no un numero fijo: hidrogenos = valencia - suma de los ordenes de
# enlace - hidrogenos que ya trae. Es lo que hace que funcione con cualquier
# ligando y no solo con los que se probaron a mano. El enlace aromatico de
# Open Babel llega como "am" y cuenta 1,5, que es lo que hace que un carbono
# de anillo con dos vecinos aromaticos pida UN hidrogeno y no dos.
#
# Medido el 28 sep 2026 con el ligando de prueba (CHEMBL1078178, 33 atomos
# en la pose): pide 20 H, uno en cada carbono de anillo, dos en cada metileno
# de la morfolina, y ninguno en los nitrogenos y oxigenos que ya estan
# completos. Suman 53 atomos, que es exactamente el numero del mol2 que
# daba antechamber desde el SMILES. Que cuadre la cuenta desde dos caminos
# distintos es lo que hizo confiar en el metodo.
_VALENCIA = {"C": 4, "N": 3, "O": 2, "S": 2, "P": 3, "F": 1, "CL": 1,
             "BR": 1, "I": 1}
_ORDEN_ENLACE = {"1": 1.0, "2": 2.0, "3": 3.0, "am": 1.5, "ar": 1.5,
                 "du": 0.0, "un": 0.0, "nc": 0.0}


def _normalizar(v):
    m = math.sqrt(sum(c * c for c in v))
    return [c / m for c in v] if m > 1e-9 else None


def _perpendicular(v):
    """Un vector unitario perpendicular a v, tomando el eje con menos componente."""
    referencia = [1.0, 0.0, 0.0] if abs(v[0]) < 0.9 else [0.0, 1.0, 0.0]
    p = [referencia[1] * v[2] - referencia[2] * v[1],
         referencia[2] * v[0] - referencia[0] * v[2],
         referencia[0] * v[1] - referencia[1] * v[0]]
    return _normalizar(p)


def _rotar(v, eje, angulo):
    """Gira v alrededor de eje un angulo en radianes (Rodrigues)."""
    c, s = math.cos(angulo), math.sin(angulo)
    cruz = [eje[1] * v[2] - eje[2] * v[1],
            eje[2] * v[0] - eje[0] * v[2],
            eje[0] * v[1] - eje[1] * v[0]]
    punto = sum(eje[k] * v[k] for k in range(3))
    salida = [v[k] * c + cruz[k] * s + eje[k] * punto * (1.0 - c) for k in range(3)]
    return _normalizar(salida)


def _combinar(vs, factores):
    """Suma vs[0]*factores[0] + vs[1]*factores[1] + ... y normaliza.

    OJO, no usar zip con factores de mas: zip corta por el mas corto y se
    pierde el factor del primer vector, que es el que marca la DIRECCION. Con
    esa version, los tres H de un metilo salian en el mismo sitio y los dos de
    un metileno tambien, y el motor devolvia vdW del orden de 1e17 sin
    quejarse.
    """
    if not vs or any(v is None for v in vs) or len(vs) != len(factores):
        return None
    salida = [0.0, 0.0, 0.0]
    for v, f in zip(vs, factores):
        for k in range(3):
            salida[k] += v[k] * f
    return _normalizar(salida)


def _elemento(nombre):
    """Elemento de un nombre de atomo de mol2: 'C1' -> 'C', 'Cl' -> 'CL'."""
    letras = ""
    for c in nombre:
        if c.isalpha():
            letras += c
        else:
            break
    return letras.upper()


def _es_hidrogeno(nombre, tipo):
    """True si el atomo es un H. Por ELEMENTO, no por el tipo GAFF: el mol2
    que sale del pdbqt lleva 'H' como tipo, y el tipado lleva 'h1', 'hn'...,
    y comparar el tipo con ("H", "D") con upper() jamas coincidio. Era el bug
    que hacia que los H se colaran en la cuenta de los atomos pesados."""
    return _elemento(nombre) in ("H", "D") or tipo.split(".")[0].upper() == "H"


def _leer_mol2(ruta):
    """(atomos, enlaces) de un mol2. Los atomos son listas de campos tal cual,
    y los enlaces triplas (a, b, tipo)."""
    atomos, enlaces, seccion = [], [], None
    with open(ruta, errors="replace") as fh:
        for linea in fh:
            if linea.startswith("@<TRIPOS>"):
                seccion = linea.strip()
                continue
            p = linea.split()
            if seccion == "@<TRIPOS>ATOM" and len(p) >= 6:
                atomos.append(p)
            elif seccion == "@<TRIPOS>BOND" and len(p) >= 4:
                try:
                    enlaces.append((int(p[1]), int(p[2]), p[3]))
                except ValueError:
                    continue
    return atomos, enlaces


def _hidrogenos_que_faltan(atomos, enlaces):
    """(que le falta a cada atomo, vecinos) por VALENCIA.

    Es el corazon del metodo. Sin esto habia que hardcodear cuantos H lleva
    cada tipo de carbono, y en cuanto salia un ligando raro el numero
    fallaba y el prmtop salia con un H de mas (o de menos) pegado al que sea.
    """
    vecinos, orden = {}, {}
    for a, b, t in enlaces:
        vecinos.setdefault(a, []).append(b)
        vecinos.setdefault(b, []).append(a)
        o = _ORDEN_ENLACE.get(t, 1.0)
        orden[(a, b)] = o
        orden[(b, a)] = o
    ya_tiene = {}
    for i, p in enumerate(atomos, 1):
        if not _es_hidrogeno(p[1], p[5]):
            continue
        for v in vecinos.get(i, []):
            if not _es_hidrogeno(atomos[v - 1][1], atomos[v - 1][5]):
                ya_tiene[v] = ya_tiene.get(v, 0) + 1
    faltan = {}
    for i, p in enumerate(atomos, 1):
        if _es_hidrogeno(p[1], p[5]):
            continue
        val = _VALENCIA.get(_elemento(p[1]), 0)
        if not val:
            continue
        gastado = sum(orden.get((i, v), 1.0) for v in vecinos.get(i, []))
        n = int(round(val - gastado - ya_tiene.get(i, 0)))
        if n > 0:
            faltan[i] = n
    return faltan, vecinos


def leer_receptor(pdb, radio):
    """Atomos del receptor cerca de una zona, para que los H del ligando los
    esquiven.

    POR QUE ESTO HACE FALTA (medido el 28 sep 2026 en Oracle): colocar los H
    mirando solo al ligando funciona hasta que el ligando se apoya en la
    pared del bolsillo. Entonces el H sale por el lado de la pared y se mete
    DENTRO de un atomo del receptor: en la corrida 18 el H5 del ligando
    acabo a 0,269 A de un LEU15, y con eso la energia de vdW del complejo
    salio 1,98e15 kcal/mol. No es un numero feo: es dos atomos en el mismo
    sitio, y el motor no protesta.

    Se pasan solo los atomos a menos de `radio` del atomo padre, que para una
    capa de 6 A son unos cuantos cientos: mirarlos todos (4.866) por cada uno
    de los 72 candidatos de cada H seria un minuto por ligando.
    """
    at = []
    if not pdb or not os.path.exists(pdb):
        return at
    with open(pdb, errors="replace") as fh:
        for linea in fh:
            if not linea.startswith(("ATOM", "HETATM")):
                continue
            try:
                at.append((float(linea[30:38]), float(linea[38:46]),
                           float(linea[46:54])))
            except ValueError:
                continue
    return at


def anadir_hidrogenos(mol2, destino, receptor_pdb=None, radio=6.0):
    """Anade al mol2 los hidrogenos que le falten, con geometria ideal.

    EL MOL2 QUE ENTRA YA TIENE LA POSE DE VINA DENTRO y su conectividad: sale
    de `obabel -ipdbqt` (Vina escribe las ramas ROOT/BRANCH, y Open Babel si las
    lee; medido el 28 sep 2026: 33 atomos y 36 enlaces para el ligando de
    prueba, con 20 H anadidos son 53 y 56, los mismos recuentos que daba
    antechamber desde el SMILES).

    POR QUE SE COLOCAN AQUI Y NO CON AmberTools (medido el 28 sep 2026 en
    Oracle, esta build es conda-forge AmberTools 21): no hay manera de que un
    tleap los ponga.
      - `addHydrogens` NO EXISTE. `desc addHydrogens` sale "STRING (with no
        reference)", que en tleap significa que el parser lo trata como texto.
      - `hadd` TAMPOCO existe, mismo sintoma, mismo mensaje de parser.
      - `addH` SI existe, y coloca los hidrogenos, pero los crea SIN TIPO:
        sale "FATAL: Atom .R<LIG 1>.A<H32 32> does not have a type" para todos,
        y el prmtop le queda de 0 bytes.
      - `pdb4amber` (que si esta) no anade H a un ligando nuevo: 33 atomos
        antes y 33 despues, porque sin H de referencia no hay nada que
        completar.
      - `reduce` (compilado a mano, no venia en AmberTools) arranca y no anade
        nada a un residuo UNK: 33 y 33.
      - `antechamber` desde mol2 NO anade H, y desde pdb tampoco.
      - `obabel -h` sobre un mol2 GAFF anade 1 H de 32 y cree que el resto ya
        esta (interpreta los tipos GAFF como "ya hidrogenado").
    Y antechamber, que si sabe, aqui no sirve: al final lo que se le pasa ya
    viene con los H dentro y sus coordenadas.

    LO QUE SE HACE, en este orden:
      1. leer el mol2 y sus enlaces
      2. contar por VALENCIA cuantos H le faltan a cada atomo pesado, sin
         contar los que ya trae (Vina deja los H POLARES, los del N-H)
      3. colocar cada H con distancia de enlace estandar y en direccion
         opuesta a la suma de los vecinos. Si a un atomo le faltan dos o
         tres, se abren en cono alrededor de esa direccion, que es la
         geometria ideal del metilo y del metileno.
      4. COMPROBAR QUE NO HAYA SOLAPES y devolver None si los hay.

    LO QUE ESTO COSTA (y por que da igual): los H no los pone el algoritmo
    completo de Amber sino una regla de geometria ideal, y el tipo GAFF se lo
    pone despues antechamber, por contexto, que ya es mejor que antes. Eso
    sesga el valor ABSOLUTO del MM-GBSA. NO sesga el RANKING, que es lo
    unico que se usa aqui: todos los ligandos pasan por el mismo codigo, con
    la misma regla y el mismo sesgo, y lo que cambia entre ellos es como se
    acomoda el ligando en el sitio.
    """
    atomos, enlaces = _leer_mol2(mol2)
    if len(atomos) < 5 or not enlaces:
        return None
    receptor = leer_receptor(receptor_pdb, radio)
    faltan, vecinos = _hidrogenos_que_faltan(atomos, enlaces)
    if not faltan:
        return None

    # coordenadas de los atomos que ya estan
    coord = {}
    for i, p in enumerate(atomos, 1):
        try:
            coord[i] = (float(p[2]), float(p[3]), float(p[4]))
        except (ValueError, IndexError):
            return None

    # 3. colocar
    nuevos = []
    for padre in sorted(faltan):
        origen = coord.get(padre)
        if origen is None:
            continue
        vecinos_del_padre = [v for v in vecinos.get(padre, []) if v in coord]
        if not vecinos_del_padre:
            continue
        # los vecinos ya colocados (incluidos los H POLARES que trajo Vina)
        # cuentan para la suma: es lo que hace que el segundo H de un NH salga
        # al lado OPUESTO del primero y no encima.
        suma = [0.0, 0.0, 0.0]
        for v in vecinos_del_padre:
            for k in range(3):
                suma[k] += coord[v][k] - origen[k]
        atras = _normalizar([-c for c in suma])
        if atras is None:
            continue
        d = _LONGITUD_ENLADE_H.get(_elemento(atomos[padre - 1][1]), 1.09)
        n = faltan[padre]
        base = _perpendicular(atras)
        # candidatos geometricos: atras para 1 H, cono de 2 o de 3
        ideales = []
        if n == 1:
            ideales = [atras]
        elif n == 2:
            ideales = [_combinar([atras, base],
                                 [-1.0 / 3.0, s * 2.0 * math.sqrt(2.0) / 3.0])
                       for s in (1.0, -1.0)]
        else:
            seg = _perpendicular(_combinar([atras, base], [-0.5, 0.5 ** 0.5]))
            ideales = [_combinar(
                [atras, base, seg],
                [-1.0 / 3.0,
                 math.sqrt(8.0) / 3.0 * math.cos(2.0 * math.pi * j / 3.0),
                 math.sqrt(8.0) / 3.0 * math.sin(2.0 * math.pi * j / 3.0)])
                for j in range(3)]
        # y cada uno, girado alrededor del eje "atras" en 24 pasos
        candidatos = []
        for ideal in ideales:
            if ideal is None:
                continue
            for giro in range(24):
                candidatos.append(ideal if giro == 0 else _rotar(
                    ideal, atras, 2.0 * math.pi * giro / 24.0))
        if not candidatos:
            continue
        # Atomos contra los que hay que alejarse: TODOS los ya colocados
        # (los del ligando, y los H del mismo padre que se acaben de poner)
        # mas los del receptor que esten cerca. Ver leer_receptor: sin el
        # receptor, un H se mete dentro de la proteina.
        #
        # EL CRITERIO ES UNO SOLO, y es el minimo de la distancia a todos
        # ellos. Antes habia dos: uno para la puntuacion (alejarse de todo lo
        # que no es vecino) y otro que descartaba los que se acercaban a un
        # vecino. Se peleaban: el mejor candidato "lejos" estaba a 1,63 de
        # los demas pero a 0,9 de un vecino, y el filtro lo echaba; el que
        # pasaba el filtro se quedaba a 0,83. Medido el 28 sep 2026 en Oracle,
        # con los dos H de dos metilenos contiguos del anillo de la
        # morfolina. Con un solo criterio el H va, por definicion, al sitio
        # mas despejado que existe, y si ese sitio esta a menos de 0,95 A es
        # que la pose viene ya congestionada y se avisa al final.
        ocupados = [c for i, c in coord.items() if i != padre]
        ocupados.extend(
            c for c in receptor
            if (c[0] - origen[0]) ** 2 + (c[1] - origen[1]) ** 2
            + (c[2] - origen[2]) ** 2 < radio * radio)
        for _ in range(n):
            mejor, mejor_puntuacion = None, -1.0
            # mejor_puntuacion esta en ANGSTROMS, y la distancia que se compara
            # dentro del bucle esta al CUADRADO. La poda tiene que usar el
            # cuadrado del mejor, no el mejor. Con unidades mezcladas pasaba
            # esto: un candidato cuyo minimo parcial era 0,5 al cuadrado se
            # comparaba con un mejor de 0,6 A, no lo echaba, y se quedaba con
            # el, aunque su minimo real fuera de 0,25 A. Medido el 28 sep 2026
            # con CHEMBL5758899: los dos H de un carbono salian a 0,532 A.
            techo = mejor_puntuacion ** 2 if mejor_puntuacion > 0 else -1.0
            for cand in candidatos:
                pos = (origen[0] + d * cand[0], origen[1] + d * cand[1],
                       origen[2] + d * cand[2])
                puntuacion = 1e9
                for c in ocupados:
                    dd = ((pos[0] - c[0]) ** 2 + (pos[1] - c[1]) ** 2
                          + (pos[2] - c[2]) ** 2)
                    if dd < puntuacion:
                        puntuacion = dd
                    if puntuacion <= techo:
                        break
                puntuacion = math.sqrt(puntuacion)
                if puntuacion > mejor_puntuacion:
                    mejor, mejor_puntuacion = pos, puntuacion
                    techo = mejor_puntuacion ** 2
            if mejor is None:
                continue
            idx = len(atomos) + len(nuevos) + 1
            nuevos.append([str(idx), "H%d" % len(nuevos), "%.4f" % mejor[0],
                           "%.4f" % mejor[1], "%.4f" % mejor[2], "H",
                           "1", "UNL1", "0.0000"])
            coord[idx] = mejor
            ocupados.append(mejor)
            enlaces.append((padre, idx, "1"))

    if not nuevos:
        return None

    # 4. escribir, con las columnas que espera tleap y los recuentos al dia
    with open(destino, "w", encoding="utf-8") as fh:
        fh.write("@<TRIPOS>MOLECULE\nLIG\n")
        fh.write(CABECERA_MOL2 % (len(atomos) + len(nuevos), len(enlaces), 0, 0, 0)
                + "\n")
        fh.write("SMALL\nUSER_CHARGES\n\n@<TRIPOS>ATOM\n")
        for i, p in enumerate(atomos + nuevos, 1):
            fh.write(LINEA_ATOM_MOL2
                     % (i, p[1], float(p[2]), float(p[3]), float(p[4]),
                        p[5], p[6], p[7], float(p[8])) + "\n")
        fh.write("@<TRIPOS>BOND\n")
        for k, (a, b, t) in enumerate(enlaces, 1):
            fh.write(LINEA_ENLACE_MOL2 % (k, a, b, t) + "\n")
        fh.write("@<TRIPOS>SUBSTRUCTURE\n"
                 "     1 LIG         1 TEMP              0 ****  ****    0 ROOT\n")

    # 5. COMPROBACION DE SOLAPAMIENTO. Sin esto, un fallo de geometria vuelve
    #    a salir como nan o como un vdW de 1e17 en las energias, que es
    #    exactamente lo que mas tiempo ha costado encontrar en este protocolo.
    #    El umbral es 0,95 A, no 0,8: el enlace mas corto que tiene un H es el
    #    N-H, a 1,01 A, y un H a 0,83 de un C es una geometria rota aunque
    #    tecnicamente no se solapen. Con 0,8 un par a 0,834 pasaba y la energia
    #    salia disparada.
    for idx, c in coord.items():
        if abs(c[0]) + abs(c[1]) + abs(c[2]) < 1e-6:
            print("  solapamiento: atomo %d sin coordenadas" % idx)
            return None
    if receptor:
        for idx, c in coord.items():
            for r in receptor:
                if (c[0] - r[0]) ** 2 + (c[1] - r[1]) ** 2 + (c[2] - r[2]) ** 2 < 0.95 ** 2:
                    print("  solapamiento: atomo %d del ligando dentro del "
                          "receptor a %.3f A" % (idx, 0.95))
                    return None
    items = sorted(coord.items())
    for a in range(len(items)):
        ia, ca = items[a]
        for b in range(a + 1, len(items)):
            ib, cb = items[b]
            d = math.sqrt(sum((ca[k] - cb[k]) ** 2 for k in range(3)))
            if d < 0.95:
                print("  solapamiento: atomos %d y %d a %.3f A" % (ia, ib, d))
                return None
    return destino


def leer_punteros(prmtop):
    """Devuelve la seccion POINTERS de un prmtop como lista de enteros.

    NO se cuentan lineas "ATOM": un prmtop no las tiene, sus atomos van en
    secciones de ancho fijo. POINTERS[3] es NATOM y POINTERS[27] es IFBOX.

    LA ESTRUCTURA (medida el 28 sep 2026 con cat -A sobre complex.prmtop): la
    primera linea es "%VERSION ...", y luego cada seccion es una linea
    "%FLAG <NOMBRE>", seguida de "%FORMAT(...)" y despues los datos. El
    %FORMAT va DETRAS de la etiqueta, no al final de la seccion, asi que no
    sirve para cerrar: hay que saltarselo. El cierre es la siguiente "%FLAG".
    NI "romper" en el %VERSION de la cabecera ni tratar el %FORMAT como cierre:
    las dos cosas se hicieron y las dos dejaban la seccion vacia (0 punteros),
    que es como se manifiesta el fallo: info_prmtop devuelve (None, None) sin
    decir por que.
    """
    seccion, valores = None, []
    with open(prmtop, errors="replace", encoding="utf-8") as fh:
        for linea in fh:
            if linea.startswith("%FORMAT"):
                continue
            if linea.startswith("%FLAG"):
                if seccion == "POINTERS":
                    break
                partes = linea.split()
                seccion = partes[1] if len(partes) > 1 else ""
                continue
            if seccion == "POINTERS":
                valores.extend(linea.split())
    return [int(v) for v in valores]


def info_prmtop(prmtop):
    """(numero de atomos, IFBOX) del prmtop, o (None, None) si no esta."""
    if not os.path.exists(prmtop):
        return None, None
    try:
        p = leer_punteros(prmtop)
    except ValueError:
        return None, None
    if len(p) < 24:
        return None, None
    # POINTERS: 0 NATOM, 1 NTYPES, ... 23 IFBOX (contados en el prmtop, no al ojo)
    return p[0], p[23]


# TERMINOS QUE GAFF NO CUBRE, medidos el 28 sep 2026 en Oracle. No es un fallo
# de codigo ni de version: los parametros no estan en gaff.dat, y tleap lo dice
# una vez que los H estan bien colocados:
#
#   Could not find angle parameter for atom types: nb - ca - ha
#
# Es el angulo N-aromatico-H, tipico de la anilina (N-H sobre un anillo). gaff.dat
# tiene ca-ca-ha 48.2 119.88 y no tiene nb-ca-ha, asi que se copia ese. Son 3
# errores y con esto se queda en 0.
#
# Los 4 TERMINOS DE TORSION que faltaban antes (ca-cp-nb-ca, cp-cp-nb-ca) NO se
# ponen aqui: los pone parmchk2, que escribe el frcmod solo al lado del mol2.
# Lo que se hace es AÑADIR estos encima, sin pisar lo que parmchk2 haya escrito.
TERMINOS_ANADIDOS = [
    # (seccion, tipo, ambiguedad, fuerza, equilibrio, por_que)
    ("ANGLE", "nb-ca-ha", 1, 48.20, 119.88, "same as ca-ca-ha"),
]


def ampliar_frcmod(work):
    """Anade a frcmod los terminos que GAFF no cubre, sin pisar los de parmchk2.

    Si parmchk2 no escribio frcmod (no falta ningun parametro), se crea uno con
    las secciones vacias: tleap necesita verlas, y sin el fichero entero el
    "loadamberparams" se queda sin cargar y los parametros tampoco se ven.
    """
    ruta = os.path.join(work, "frcmod")
    if os.path.exists(ruta):
        with open(ruta, errors="replace", encoding="utf-8") as fh:
            texto = fh.read()
    else:
        texto = ("Remark line goes here\nMASS\n\nBOND\n\nANGLE\n\nDIHE\n\n"
                 "IMPROPER\n\nNONBON\n")
    faltan = [(s, t, a, f, e, q) for (s, t, a, f, e, q) in TERMINOS_ANADIDOS
              if ("%s %s" % (s, t)) not in texto]
    if not faltan:
        return
    lineas = texto.splitlines()
    for seccion, tipo, amb, fuerza, equi, porque in faltan:
        # SIN el nombre de la seccion delante: en frcmod va en su propio
        # encabezado ("ANGLE" en su linea) y el termino va solo en la suya.
        # Con "ANGLE nb-ca-ha" en la misma linea, tlepa toma "ANGLE" como si
        # fuera un tipo de atomo, no encuentra el angulo y sigue con los 3
        # errores, sin decir nada raro.
        nueva = "%-12s %1d    %8.3f    %8.3f      %s" % (
            tipo, amb, fuerza, equi, porque)
        # la seccion va justo despues de su encabezado; si no esta, se anade al
        # final, que es donde tleap lo encuentra igual
        try:
            pos = next(i for i, l in enumerate(lineas)
                       if l.strip() == seccion)
            lineas.insert(pos + 1, nueva)
        except StopIteration:
            lineas.append(nueva)
    with open(ruta, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\n".join(lineas) + "\n")


def montar_sistemas(receptor_pdb, ligando_mol2, work, capa=None):
    ampliar_frcmod(work)
    script = SCRIPT_TLEAP.format(receptor=os.path.basename(receptor_pdb),
                                 ligando=os.path.basename(ligando_mol2),
                                 capa=capa)
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

    POR QUE (28 sep 2026, medido en Oracle): el Generalized Born es
    incompatible con la caja periodica, y mmpbsa_py_energy lo dice claro:
    "gb>0 is incompatible with periodic boundary conditions. To use this
    method set IFBOX in the PRMTOP file to 0". El GB es, por definicion, un
    modelo de solvente NO periodico, asi que quitar la caja no es un truco: es
    lo que el metodo pide.

    Se hace con parmed, que AmberTools trae, en vez de editar el prmtop a mano:
    parmed reescribe todas las secciones con los recuentos ya ajustados. Editar
    el flag a mano deja los contadores de la seccion BOX desfasados y sander
    lee memoria de mas.

    AHORA ES CASO RARO (28 sep 2026, con solvateOct): `solvateOct` no crea caja
    periodica, solo una capa de agua, asi que el prmtop sale con IFBOX=0 y no
    hay nada que quitar. Se comprueba y, si ya es 0, se copia tal cual. Antes
    de esto la funcion se ejecutaba siempre y era un crash silencioso: parmed
    no esta en el python3 del sistema, solo en el de AmberTools, y el paso 5
    se caia antes de abrir el log de MMPBSA.py, sin dejar rastro.
    """
    natom, ifbox = info_prmtop(prmtop)
    if natom is None:
        return False
    if ifbox == 0:
        shutil.copy2(prmtop, destino)
        return True
    codigo = ("from parmed import load_file\n"
              "s = load_file(%r)\n"
              "s.box = None\n"
              "s.save(%r, overwrite=True)\n"
              % (prmtop, destino))
    r = subprocess.run([PYTHON_AMBER, "-c", codigo], cwd=work,
                       capture_output=True, text=True, timeout=3600,
                       env=env_entorno_ambertools())
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
    cmd = [PYTHON_AMBER, MMPBSA_PY, "-i", os.path.basename(mmpbsa_in),
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
    env = env_entorno_ambertools()
    with open(os.path.join(work, "mmpbsa_%s.log" % etiqueta), "w",
              encoding="utf-8", errors="ignore") as fh:
        r = subprocess.run(cmd, cwd=work, stdout=fh, stderr=subprocess.STDOUT,
                           timeout=7200, env=env)
    if r.returncode != 0:
        return None
    return leer_energia_mmpbsa(salida)


def leer_energia_mmpbsa(ruta):
    """Saca la energia del .out de MMPBSA.py.

    EL FORMATO REAL (medido el 28 sep 2026 en Oracle, AmberTools 14.0) NO es el
    que estaba escrito aqui. No hay ninguna linea "DG (Energy terms)": lo que
    sale, con su encabezado, es:

        Energy Component            Average              Std. Dev.   Std. Err. of Mean
        BOND                       462.9680                0.0000              0.0000
        ANGLE                     1702.3439                0.0000              0.0000
        DIHED                     6833.2143                0.0000              0.0000
        VDWAALS                    102.0629                0.0000              0.0000
        EEL                     -42431.9303                0.0000              0.0000
        1-4 VDW                   3018.1442                0.0000              0.0000
        1-4 EEL                  26137.6224                0.0000              0.0000
        EGB                      -8100.8949                0.0000              0.0000
        ESURF                      186.3313                0.0000              0.0000
        G gas                    -4175.5746                0.0000              0.0000
        G solv                   -7914.5636                0.0000              0.0000
        TOTAL                   -12090.1382                0.0000              0.0000

    El lector de antes buscaba "DG (Energy terms)" o una linea con "vdW" y
    "elec." a la vez, y en este formato no hay ninguna de las dos cosas: el
    nombre de la fila es "VDWAALS" y "EEL", en filas separadas. Se quedaba
    (None) sin decir nada, y el paso 3 del script lo reportaba como "FALLO"
    aunque MMPBSA.py hubiera terminado bien y con numeros buenos. Tardo en
    verlo porque el fallo de verdad estaba antes (los atomos encima).

    Se lee POR NOMBRE de la fila y no por posicion, que es lo unico que no
    cambia entre versiones.
    """
    if not os.path.exists(ruta):
        return None
    filas = {}
    for linea in open(ruta, errors="replace"):
        p = linea.split()
        if len(p) < 2 or not p[0][0].isalpha():
            continue
        try:
            filas[p[0]] = float(p[1])
        except ValueError:
            continue
    if "TOTAL" not in filas:
        return None
    # G solv = EGB + ESURF, que es el termino de solvatacion completo
    solv = filas.get("G solv")
    if solv is None:
        solv = (filas.get("EGB") or 0.0) + (filas.get("ESURF") or 0.0)
    return {"dG_total": filas["TOTAL"],
            "dG_vdw": filas.get("VDWAALS"),
            "dG_elec": filas.get("EEL"),
            "dG_solv": solv}


def comprobar_solapes_sistema(rst7, umbral=0.9):
    """Primer par de atomos del sistema a menos de `umbral` A, o None.

    POR QUE (medido el 28 sep 2026 en Oracle): cuando dos atomos caen en el
    mismo sitio, el motor NO protesta. Devuelve un vdW de 1,98e15 kcal/mol y
    un TOTAL que parece un numero cualquiera. Es el fallo mas caro de
    depurar de todo este protocolo, porque no dice ni donde ni por que. Con
    esta comprobacion sale antes del MM-GBSA, con los dos indices y la
    distancia.

    El rst7 de Amber no lleva nombres de atomo, solo coordenadas, y el
    complejo son 15.000 atomos: recorrerlos todos por parejas son 118
    millones de distancias, varios minutos. Se usa una rejilla de 1 A: cada
    atomo solo se compara con los de las 26 celdas vecinas, y sale en menos
    de un segundo.
    """
    if not os.path.exists(rst7):
        return None
    with open(rst7, errors="replace") as fh:
        lineas = [l for l in fh if l.strip()]
    try:
        n = int(lineas[1])
    except (IndexError, ValueError):
        return None
    coords = []
    for linea in lineas[2:]:
        v = linea.split()
        for k in range(0, len(v) - 2, 3):
            coords.append((float(v[k]), float(v[k + 1]), float(v[k + 2])))
            if len(coords) == n:
                break
        if len(coords) == n:
            break
    rejilla = {}
    for i, c in enumerate(coords):
        rejilla.setdefault((int(c[0]), int(c[1]), int(c[2])), []).append(i)
    u2 = umbral * umbral
    for (x, y, z), idxs in rejilla.items():
        vecinos = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    vecinos.extend(rejilla.get((x + dx, y + dy, z + dz), ()))
        for a in range(len(idxs)):
            for b in range(a + 1, len(vecinos)):
                i, j = idxs[a], vecinos[b]
                if j <= i:
                    continue
                ci, cj = coords[i], coords[j]
                d2 = ((ci[0] - cj[0]) ** 2 + (ci[1] - cj[1]) ** 2
                      + (ci[2] - cj[2]) ** 2)
                if d2 < u2:
                    return (i + 1, j + 1, math.sqrt(d2))
    return None


# Minimizacion corta ANTES de calcular el MM-GBSA. Sin esto el numero no
# ordena (medido el 28 sep 2026 en Oracle, ver la nota larga de arriba): el
# dG de union salia en un rango de mas de 200 kcal/mol y el area bajo la curva
# daba 0,57, que es azar. La razon es que se estaba calculando sobre la pose de
# Vina SIN RELAJAR, con los H puestos por una regla de geometria ideal: hay
# enlaces de carbonos aromaticos a 1,38 A, angulos fuera de sitio y un par de
# H a 0,86 A. Con eso el campo no esta cerca de su minimo y el dG lleva el
# ruido de la geometria encima, que es mucho mayor que la diferencia entre un
# buen ligando y uno malo.
#
# 200 pasos de minimos sin restriccion NO deforman la proteina (el receptor de
# 9.688 atomos apenas se mueve) y SI relajan el sitio de union. Es lo que se
# hace en la practica antes de un MM-GBSA de un solo cuadro, y cuesta dos
# segundos por sistema.
#
# PERO NO ES LO QUE PASA AQUI, y por eso va DESACTIVADO por defecto (0). Lo
# que se midio el 28 sep 2026 con el ligando de prueba, minimizando el
# SISTEMA entero:
#   - dG_union pasa de -14,3 a -30,1 kcal/mol. En un sistema de 9.741 atomos
#     un dG de union de -30 es fisicamente inverosimil.
#   - el vdW del complejo cae a -5.420 kcal/mol. El vdW de una proteina
#     minimizada es +30 a +80. Un vdW muy NEGATIVO significa que los atomos se
#     han amontonado: la minimizacion los ha colapsado.
#   - sander avisa "Maximum number of minimization cycles reached": no
#     converge en 200 pasos, asi que se para en un punto que no es un minimo.
#   - el coste pasa de 36 s a 6 min y medio por ligando (200 pasos x 9.741
#     atomos x 3 sistemas).
#   - y el rst7 que escribe sander en Amber 26 es BINARIO (un formato
#     Fortran nuevo, con las etiquetas "spatial" y "atom" dentro). MMPBSA.py y
#     cpptraj lo leen, pero cualquier script que lo abra a pelo revienta con
#     ValueError. Para volver a texto hace falta parmed.
#
# CONCLUSION HONESTA: relajar la geometria es la via correcta en principio, pero
# minimizing sin restriccion y sin converger no relaja, colapsa. La version
# correcta es minimizar SOLO el ligando (53 atomos, converge en segundos) e
# inyectar esas coordenadas en el complejo, o minimizar el complejo con el
# receptor restringido con nmropt. Eso esta pendiente de medir.
#
# LA PRIMERA LINEA ES OBLIGATORIA Y NO ES UN ADORNO. Medido el 28 sep 2026 en
# Oracle, con Amber 26: un mdin que empieza directamente por `&cntrl` hace que
# sander responda "Could not find cntrl namelist" y salga con codigo 1, sin
# calcular nada y sin decir cual es el problema. La razon es que sander se
# come la primera linea como TITULO y luego busca el namelist a partir de la
# segunda; si la primera es el propio `&cntrl`, se lo come y no encuentra
# ninguno. Se probo con `&end` en vez de `/`, con el bloque entero en la linea
# uno, con seis y con dos espacios de sangria, con `mdin` en vez del nombre
# pasado con -i, con AMBERHOME y con amber.sh sourced: las ocho fallan igual.
# Con una linea de texto delante, rc=0 y el rst7 sale. (El mdout lo delata:
# "Here is the input file:" salia VACIO, con la linea de titulo sale el
# contenido.) El MM-GBSA usa MMPBSA.py, que genera su propio mdin, y por eso
# este fallo solo aparece en la minimizacion.
INPUT_MINIMIZAR = """minimizacion corta del sistema, para relajar la geometria
&cntrl
  imin=1, maxcyc={pasos}, ntb=0, ntpr=100, cut=999.0, ifqnt=0,
/
"""


def minimizar(prmtop, rst7, work, etiqueta, pasos=200):
    """Minimiza un sistema y devuelve el rst7 nuevo, o None si no se pudo.

    NO se usa la salida de sander para nada mas: solo las coordenadas. El prmtop
    no cambia al minimizar (las fuerzas son las mismas), asi que el MM-GBSA
    sigue con el prmtop que dio tleap y solo cambia el rst7.
    """
    entrada = os.path.join(work, "min_%s.in" % etiqueta)
    salida = os.path.join(work, "min_%s.out" % etiqueta)
    nuevo = os.path.join(work, "%s_min.rst7" % etiqueta)
    with open(entrada, "w", encoding="utf-8") as fh:
        fh.write(INPUT_MINIMIZAR.format(pasos=pasos))
    cmd = [SANDER, "-O", "-i", os.path.basename(entrada),
           "-o", os.path.basename(salida), "-p", os.path.basename(prmtop),
           "-c", os.path.basename(rst7), "-r", os.path.basename(nuevo)]
    run(cmd, work, "minimizar.log")
    if not os.path.exists(nuevo) or os.path.getsize(nuevo) == 0:
        return None
    return nuevo


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--receptor", required=True)
    ap.add_argument("--pose", required=True,
                    help="PDBQT de salida de Vina (el *_out.pdbqt)")
    ap.add_argument("--capa", type=float, default=0.0,
                    help="capa de agua explicita en A. 0 (por defecto) es lo "
                         "correcto con igb=5; ver el comentario de SCRIPT_TLEAP")
    ap.add_argument("--minimizar", type=int, default=0,
                    help="pasos de minimizacion del SISTEMA antes del "
                         "MM-GBSA. Por defecto 0, y hay un motivo medido: ver "
                         "INPUT_MINIMIZAR. Con 200 el dG_union se va de -14 a "
                         "-30 kcal/mol, el vdW del complejo cae a -5.420 (no "
                         "es fisico), sander avisa 'Maximum number of "
                         "minimization cycles reached' y cada ligando pasa de "
                         "36 s a 6 min y medio")
    args = ap.parse_args()
    capa = args.capa

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
    print("MM-GBSA con AmberTools  ·  ligando tipado desde la pose")
    print("=" * 66)
    print("receptor: %s" % os.path.basename(receptor))
    print("pose:     %s" % os.path.basename(pose))
    print("trabajo:  %s\n" % work)

    try:
        # 1-3. pose -> mol2 con H -> tipado con GAFF
        mol2 = pose_a_mol2(pose, work, receptor)
        if not mol2:
            print("1. pose -> mol2 GAFF  FALLO (ver antechamber.log)")
            return 1
        print("1. pose -> mol2 con H -> GAFF (antechamber + parmchk2)  OK")

        # 2. tleap
        shutil.copy2(receptor, os.path.join(work, os.path.basename(receptor)))
        if not montar_sistemas(receptor, mol2, work, capa):
            print("2. tleap  FALLO (ver tleap.log)")
            return 1
        # Con solvateBox esto eran 90.570 atomos y el MM-GBSA se colgaba. Sin
        # agua (lo correcto con igb=5) son unos 9.700. Si aparece un numero de
        # cinco cifras, se ha vuelto a solvatar.
        for etiqueta in ("receptor", "ligand", "complex"):
            natom, ifbox = info_prmtop(os.path.join(work,
                                                    "%s.prmtop" % etiqueta))
            print("2. %-8s %s atomos, IFBOX %s" % (etiqueta, natom, ifbox))
        # Y antes de gastar un minuto en el MM-GBSA, mirar que no haya dos
        # atomos encima: si los hay, la energia sale absurda y no dice por que.
        choque = comprobar_solapes_sistema(os.path.join(work, "complex.rst7"))
        if choque:
            print("2. FALLO: atomos %d y %d del sistema a %.3f A"
                  % choque)
            return 1
        print("2. sin atomos solapados en el complejo  OK")

        # 3. MM-GBSA, sobre la geometria MINIMIZADA
        res = {}
        for etiqueta in ("complex", "receptor", "ligand"):
            rst7 = "%s.rst7" % etiqueta
            if args.minimizar:
                nuevo = minimizar("%s.prmtop" % etiqueta, rst7, work, etiqueta,
                                args.minimizar)
                if nuevo is None:
                    print("3. minimizacion de %s  FALLO" % etiqueta)
                    return 1
                rst7 = "%s_min.rst7" % etiqueta
            r = energia_mmgbsa("%s.prmtop" % etiqueta, rst7, work, etiqueta)
            if r is None:
                print("3. MM-GBSA %s  FALLO (ver mmpbsa_%s.log)"
                      % (etiqueta, etiqueta))
                return 1
            res[etiqueta] = r
            print("3. MM-GBSA %-9s G = %10.1f kcal/mol"
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
