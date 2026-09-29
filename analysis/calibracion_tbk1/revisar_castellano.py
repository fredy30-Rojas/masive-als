"""Busca palabras inglesas que se hayan colado en la documentacion.

POR QUE ESTE SCRIPT (29 sep 2026): se colaron dos veces palabras en ingles
dentro de los textos en castellano, una de ellas en un titulo ("la lista lo
estaba ainingolando", de "engañando"). Las dos se ganglionaron a la vez por
voz, asi que un fallo mio llego a su oido. Escribirlos sin darme cuenta es
facil cuando el modelo esta cansado, y revisarlos a ojo no escala.

QUE HACE
Busca un diccionario de palabras inglesas que en un texto en castellano de
este proyecto son un fallo, y descarta lo legitimo: lo que va entre acentos
grave, lo citado textualmente de un programa, y en los .py SOLO mira los
comentarios y las cadenas de documentacion, porque la sintaxis de Python es en
inglesa por definicion y no tiene sentido marcarla.

Uso:  python revisar_castellano.py            (toda la carpeta)
       python revisar_castellano.py FICHERO...
"""
import glob
import os
import re
import sys

# Palabras inglesas que en un texto en castellano de este proyecto son un fallo.
INGLES = {
    "the", "and", "with", "without", "for", "from", "that", "this", "have",
    "has", "are", "was", "were", "been", "being", "will", "would", "should",
    "could", "there", "their", "they", "them", "then", "than", "when", "what",
    "which", "where", "while", "about", "after", "before", "because", "only",
    "also", "more", "most", "some", "such", "each", "other", "into", "over",
    "under", "between", "both", "same", "very", "much", "many", "now", "here",
    "knows", "correction", "know", "said", "says", "tell", "told", "give",
    "given", "made", "make", "makes", "use", "used", "using", "need", "needs",
    "want", "wants", "good", "better", "best", "bad", "worse", "worst",
    "first", "second", "third", "last", "next", "one", "two", "three",
    "sourced", "builder", "galionadas", "gala", "seek", "wanted", "got",
    "start", "done", "keep", "let", "its", "does", "did", "done", "next",
}

# Palabras que son ingles pero aqui son legitimas: mensaje literal de un
# programa, nombre de fichero, o terminologia que el proyecto ya usa.
PERMITIDAS = {
    "string", "with", "no", "reference", "file", "data", "set", "type", "name",
    "list", "value", "error", "warning", "fatal", "info", "usage", "input",
    "output", "field", "bond", "atom", "atoms", "resid", "residue", "link",
    "min", "max", "use", "used", "format", "version", "index", "total",
    "space", "single", "model", "step", "first", "last", "all", "is", "it",
    "in", "of", "to", "a", "the", "and", "or", "not", "no", "yes", "id",
    "log", "run", "end", "begin", "point", "grid", "cell", "marca", "solo",
}


def solo_comentarios(ruta):
    """En un .py solo interesan los comentarios y las cadenas de documentacion.

    La sintaxis de Python es inglesa (`for`, `from`, `while`, `None`) y marcarla
    daria cientos de falsos positivos: la primera version del script lo hacia
    y salia con 388 Lineas, casi todas `for` y `from`.
    """
    if not ruta.endswith(".py"):
        return None
    dentro = False
    out = []
    with open(ruta, encoding="utf-8", errors="replace") as fh:
        for num, linea in enumerate(fh, 1):
            t = linea.strip()
            es_doc = False
            if dentro:
                out.append((num, t))
                dentro = not t.endswith('"""')
                continue
            if t.startswith('"""') or t.startswith("'''"):
                dentro = not (t.endswith('"""') and len(t) > 3)
                out.append((num, t))
                es_doc = True
                continue
            if t.startswith("#"):
                out.append((num, t))
                es_doc = True
            if es_doc:
                continue
    return out


def revisar(ruta):
    fallos = []
    if ruta.endswith(".py"):
        lineas = solo_comentarios(ruta) or []
    else:
        lineas = []
        with open(ruta, encoding="utf-8", errors="replace") as fh:
            lineas = [(n, l.rstrip()) for n, l in enumerate(fh, 1)]
    for num, linea in lineas:
        limpio = re.sub(r"`[^`]*`", " ", linea)
        limpio = re.sub(r'"[^"]*"', " ", limpio)
        limpio = re.sub(r"\'[^\']*\'", " ", limpio)
        for palabra in re.findall(r"[A-Za-zá-úÁ-Ú']{3,}", limpio):
            b = palabra.lower()
            if b in INGLES and b not in PERMITIDAS:
                fallos.append((num, palabra, linea[:90]))
    return fallos


def main(rucas):
    if not rucas:
        base = os.path.dirname(os.path.abspath(__file__))
        rucas = sorted(glob.glob(os.path.join(base, "*.py"))
                       + glob.glob(os.path.join(base, "*.md")))
    total = 0
    for r in rucas:
        fallos = revisar(r)
        if not fallos:
            continue
        print("%s" % os.path.basename(r))
        for num, palabra, linea in fallos:
            print("   %4d  %-14s  %s" % (num, palabra, linea))
        total += len(fallos)
    print("\n%d palabras inglesas a revisar en %d ficheros" % (total, len(rucas)))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
