#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""verificar_xl20_vision.py — segunda lectura de la estructura de XL20, esta vez con vision.

POR QUE EXISTE
--------------
El SMILES de XL20 (`controles_tdp43_xl20.csv`) salio de un OCR de estructuras (DECIMER
2.7.2) sobre el dibujo de la Figura Suplementaria 1(a) de Gao 2026. Un OCR quimico se
equivoca, y ese compuesto abre una linea de trabajo: conviene una segunda lectura
independiente antes de darlo por bueno.

AQUI NO SE LE PIDE EL SMILES AL MODELO
--------------------------------------
Un modelo de vision escribe SMILES mal con mucha facilidad, y un SMILES mal escrito no se
nota. Lo que si sabe hacer es mirar y contar: cuantos anillos hay, cuales son aromaticos,
cuantos nitrogenos, que grupos funcionales aparecen. Eso se compara despues, a mano, con lo
que dice el SMILES leido por el OCR. Es la misma idea que se uso con los fragmentos de
Nshogoza, alli con la descripcion del modelo de vision de Ollama.

Se hacen dos preguntas:
  1. el recorte del dibujo tal cual -> descripcion de lo que se ve;
  2. la imagen de comparacion (recorte original a la izquierda, dibujo de RDKit a la
     derecha) -> representan la misma molecula, y donde estan las diferencias.

Requiere Ollama en marcha (no arranca solo en el PC de Fredy: `ollama serve`).

TRAMPA MEDIDA DE `qwen3-vl:8b` (25 sep 2026): se pasa el presupuesto de tokens pensando y
devuelve el campo `content` VACIO. Con `num_predict` de 500 y de 1500 el `done_reason`
es `length` y el contenido vacio, mientras `message.thinking` trae la descripcion entera
(y correcta). Por eso aqui se piden 4000 tokens y se guardan las DOS cosas: el razonamiento
es la parte util, y con el presupuesto corto se pierde sin que se note.

Uso:
    python analysis/verificar_xl20_vision.py [modelo]

Salida: analysis/_xl20_fig/verificacion_vision.txt (razonamiento y respuesta de las dos
preguntas, con el tiempo).
"""
import base64
import json
import os
import sys
import time
import urllib.error

BASE = os.path.dirname(os.path.abspath(__file__))
IMGS = os.path.join(BASE, "_xl20_fig")
API = "http://127.0.0.1:11434"
MODELO = sys.argv[1] if len(sys.argv) > 1 else "qwen3-vl:8b"

sys.path.insert(0, r"C:/Users/Fredy/tools")
from probar_ollama import limpiar, peticion  # noqa: E402  (una sola copia de la llamada)

PREGUNTA_DIBUJO = (
    "This image is the drawing of one molecule (a chemical structure). Describe it "
    "carefully and only what you can see: how many rings there are and of what kind "
    "(benzene, purine, saturated ring...), how many nitrogen atoms, and which functional "
    "groups appear (amino, hydroxyl, amide...). Be concrete and short. Do not write SMILES."
)

PREGUNTA_LISTA = (
    "Look at the molecule in this image and answer this fixed checklist, one line per "
    "item, with only the word YES or NO and nothing else:\n"
    "1. Is there an aromatic six-membered carbocycle (a benzene/phenyl ring)?\n"
    "2. Is there a saturated six-membered carbocycle (a cyclohexane-like ring)?\n"
    "3. Is there a fused bicyclic nitrogen heterocycle with an amino group (adenine-like)?\n"
    "4. Is there a hydroxyl group (OH) on the saturated ring?\n"
    "5. Is there a nitrogen that carries both a methyl and a CH2-phenyl?\n"
)

# Lo que el SMILES leido por el OCR dice que tiene que haber. Se busca por vocabulario,
# porque el modelo describe con palabras: si nombra los cinco fragmentos, coincide.
ITEMS = {
    "benceno/phenilo": ("benzen", "phenyl", "bencen", "fenil", "arom"),
    "ciclohexano saturado": ("cyclohexan", "ciclohexan", "saturat", "aliphatic"),
    "purina/adenina con amino": ("purin", "adenin", "NH2", "H\u2082N", "amino"),
    "hidroxilo": ("hydroxyl", "hidroxil", "OH", "HO"),
    "N con metilo y bencilo": ("methyl", "metil", "CH3", "CH\u2083", "benzyl", "bencil"),
}


PREGUNTA_COMPARACION = (
    "This image has two drawings of the same molecule side by side: the left one is a crop "
    "from a scientific figure, the right one was drawn by a program from a SMILES string. "
    "Do they represent the same molecule? Compare the number of rings, which rings are "
    "aromatic, how many nitrogen atoms there are and which functional groups appear, and "
    "point out any difference you see. Do not write SMILES."
)


def base64_de(ruta):
    with open(ruta, "rb") as f:
        return base64.b64encode(f.read()).decode("ascii")


def preguntar(texto, imagenes, num_predict=4000, timeout=1800):
    datos = {
        "model": MODELO,
        "messages": [{"role": "user", "content": texto, "images": imagenes}],
        "stream": False,
        "think": False,
        "keep_alive": "5m",
        "options": {"num_predict": num_predict, "temperature": 0.2},
    }
    t0 = time.time()
    try:
        r = peticion("/api/chat", datos, timeout)
    except urllib.error.HTTPError as e:
        detalle = e.read().decode("utf-8", "replace")[:200]
        if "think" in detalle.lower():
            datos.pop("think", None)
            r = peticion("/api/chat", datos, timeout)
        else:
            raise RuntimeError(detalle)
    mensaje = r.get("message", {})
    return (limpiar(mensaje.get("content", "")),
            limpiar(mensaje.get("thinking", "")),
            r.get("done_reason", ""),
            time.time() - t0)


def comprobar_items(texto):
    """Que fragmentos de los que dice el SMILES nombra el modelo al mirar el dibujo."""
    bajo = texto.lower()
    salida = []
    for nombre, claves in ITEMS.items():
        salida.append((nombre, any(k.lower() in bajo for k in claves)))
    return salida


def main():
    recorte = os.path.join(IMGS, "XL20.png")
    comparacion = os.path.join(IMGS, "XL20_comparacion.png")
    for ruta in (recorte, comparacion):
        if not os.path.exists(ruta):
            sys.exit("falta %s (corre antes leer_figura_xl20.py)" % ruta)

    print("modelo: %s" % MODELO)
    resp = []
    for titulo, pregunta, ruta in (
        ("1) el dibujo de XL20 tal cual", PREGUNTA_DIBUJO, recorte),
        ("2) recorte original frente al dibujo de RDKit", PREGUNTA_COMPARACION, comparacion),
        ("3) lista fija de fragmentos", PREGUNTA_LISTA, recorte),
    ):
        contenido, razonamiento, motivo, seg = preguntar(pregunta, [base64_de(ruta)])
        print("\n--- %s (%.1f s, termino por %s) ---" % (titulo, seg, motivo))
        print("RESPUESTA:\n%s" % (contenido or "(vacia)"))
        print("RAZONAMIENTO:\n%s" % (razonamiento or "(vacio)"))
        if titulo.startswith("3"):
            print("FRAGMENTOS QUE NOMBRA (sobre el razonamiento):")
            for nombre, visto in comprobar_items(razonamiento):
                print("   %-26s %s" % (nombre, "SI" if visto else "NO"))
        resp.append((titulo, contenido, razonamiento, motivo, seg))

    salida = os.path.join(IMGS, "verificacion_vision.txt")
    with open(salida, "w", encoding="utf-8") as f:
        f.write("VERIFICACION DE LA LECTURA DE XL20 CON MODELO DE VISION (%s)\n" % MODELO)
        f.write("imagen 1: %s\nimagen 2: %s\n" % (recorte, comparacion))
        f.write("OJO: qwen3-vl:8b con presupuesto corto devuelve 'content' vacio y deja todo "
                "en el razonamiento. Aqui van los dos.\n")
        for titulo, contenido, razonamiento, motivo, seg in resp:
            f.write("\n=== %s (%.1f s, termino por %s) ===\n" % (titulo, seg, motivo))
            f.write("RESPUESTA:\n%s\n\nRAZONAMIENTO:\n%s\n" % (contenido or "(vacia)", razonamiento or "(vacio)"))
    print("\nguardado: %s" % salida)


if __name__ == "__main__":
    main()
