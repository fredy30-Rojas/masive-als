# -*- coding: utf-8 -*-
"""Sitúa el sitio de unión de ARN de TDP43_v2 y lo compara con la caja del cribado.

Por qué existe (17 sep 2026): el informe del recálculo dejó como único cabo que
bloquea una decisión el sitio de acoplamiento de TDP43_v2. La evidencia era
indirecta: los ligandos conocidos entierran menos que el fondo, sus poses se
dispersan 8,1 Å y ninguna métrica los salva. La hipótesis era que se está
acoplando en una cavidad profunda mientras los ligandos verdaderos se unen a la
superficie plana de RRM1/RRM2.

Aquí se comprueba con datos y no con literatura: el receptor TDP43_v2 sale de la
estructura 4BS2, que es el TDP-43 con ARN unido. Así que el sitio de unión de
ARN no hay que buscarlo: está en la propia estructura. Se calculan los residuos
que tocan al ARN, su centro, y se mide cuánto se parece eso a la caja que usa el
cribado (centro 24.23, 16.89, -15.87; tamaño 26).

Salida: por pantalla, y un JSON con la caja propuesta.
"""
import json
import math
import os

BASE = r"C:\Users\Fredy\masive-als"
PUENTE = os.path.join(BASE, "analysis", "4bs2_complejo_rna.pdb")   # 4BS2 entero (20 modelos)
RECEPTOR = os.path.join(BASE, "gpu_dock", "TDP43_v2.pdbqt")        # el que usa el cribado
SALIDA = os.path.join(BASE, "analysis", "rescoring_local", "caja_tdp43v2_arn.json")
CAJA_VIEJA = ((24.23, 16.89, -15.87), 26.0)
CORTE_CONTACTO = 4.0     # Å: distancia atomo-atomo para decir "toca al ARN"


def leer_atomos(linea):
    """Saca (resnum, resname, nombre, x, y, z) de una linea ATOM/HETATM."""
    return (int(linea[22:26]), linea[17:20].strip(), linea[12:16].strip(),
            float(linea[30:38]), float(linea[38:46]), float(linea[46:54]))


def leer_modelo_1():
    """Devuelve (proteinas, rna) del primer modelo de la estructura."""
    proteinas, rna = [], []
    modelo = 0
    for ln in open(PUENTE, encoding="utf-8", errors="replace"):
        if ln.startswith("MODEL"):
            modelo += 1
            continue
        if ln.startswith("ENDMDL"):
            if modelo >= 1:
                break
            continue
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        if modelo != 1:
            continue
        (proteinas if ln[21] == "A" else rna).append(leer_atomos(ln))
    return proteinas, rna


def leer_receptor():
    atomos = []
    for ln in open(RECEPTOR, encoding="utf-8", errors="replace"):
        if ln.startswith(("ATOM", "HETATM")):
            atomos.append(leer_atomos(ln))
    return atomos


def dist(a, b):
    return math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2)


def main():
    prot, rna = leer_modelo_1()
    recept = leer_receptor()
    print("4BS2 modelo 1: %d atomos de proteina (cadena A), %d de ARN (cadena B)" % (len(prot), len(rna)))
    print("Receptor del cribado: %d atomos" % len(recept))

    # ── 1) Los residuos que tocan al ARN ──────────────────────────
    tocan = set()
    for p in prot:
        for r in rna:
            if abs(p[3] - r[3]) > CORTE_CONTACTO:
                continue
            if dist((p[3], p[4], p[5]), (r[3], r[4], r[5])) <= CORTE_CONTACTO:
                tocan.add((p[0], p[1]))
                break
    print("\nResiduos de la proteina a menos de %.1f A del ARN: %d" % (CORTE_CONTACTO, len(tocan)))
    por_num = sorted(tocan)
    print("  " + ", ".join("%s%d" % (n, r) for r, n in por_num))

    # ── 2) Misma numeracion que el receptor del cribado? ─────────
    # Si las coordenadas fueran de otro modelo o de otra copia, la caja saldria
    # en el sitio equivocado; hay que comprobarlo y no darlo por hecho.
    recept_por_res = {}
    for a in recept:
        recept_por_res.setdefault(a[0], []).append(a)
    comunes = [r for r, n in por_num if r in recept_por_res]
    desvios = []
    for ln in open(PUENTE, encoding="utf-8", errors="replace"):
        pass
    # CA de cada residuo del primer modelo de la estructura
    ca_4bs2 = {p[0]: (p[3], p[4], p[5]) for p in prot if p[2] == "CA"}
    ca_rec = {a[0]: (a[3], a[4], a[5]) for a in recept if a[2] == "CA"}
    for r in sorted(set(ca_4bs2) & set(ca_rec)):
        d = dist(ca_4bs2[r], ca_rec[r])
        desvios.append((r, d))
    if desvios:
        media = sum(d for _, d in desvios) / len(desvios)
        peor = max(desvios, key=lambda x: x[1])
        print("\nMismo sistema de coordenadas que el receptor? %d residuos comparados" % len(desvios))
        print("  desviacion media del CA: %.2f A (peor: residuo %d, %.2f A)" % (media, peor[0], peor[1]))
        if media > 2.0:
            print("  OJO: las coordenadas NO coinciden. Habria que superponer antes de sacar la caja.")
    print("  residuos de la interfaz presentes en el receptor: %d de %d" % (len(comunes), len(por_num)))

    # ── 3) Centro y tamaño del sitio de union de ARN ─────────────
    interfaz = [p for p in prot if (p[0], p[1]) in tocan]
    centro = tuple(sum(p[i] for p in interfaz) / len(interfaz) for i in (3, 4, 5))
    ext = [(max(p[i] for p in interfaz) - min(p[i] for p in interfaz)) for i in (3, 4, 5)]
    tam = max(ext)
    print("\nSitio de union de ARN (de la estructura):")
    print("  centro: %.2f, %.2f, %.2f" % centro)
    print("  extension real: %.1f x %.1f x %.1f A" % tuple(ext))
    lade = math.ceil((tam + 8) / 2) * 2      # 8 A de margen, en numero par
    print("  lado de caja propuesto: %d A" % lade)

    # ── 4) Comparacion con la caja que usa el cribado ────────────
    cv, sv = CAJA_VIEJA
    dist_centros = dist(centro, cv)
    dentro = sum(1 for p in interfaz
                 if all(abs(p[i + 3] - cv[i]) <= sv / 2 for i in range(3)))
    print("\nCaja del cribado: centro %.2f, %.2f, %.2f  tamaño %g" % (cv + (sv,)))
    print("  distancia entre centros: %.1f A" % dist_centros)
    print("  atomos de la interfaz de ARN dentro de esa caja: %d de %d (%.0f%%)"
          % (dentro, len(interfaz), 100.0 * dentro / len(interfaz)))

    # Cuanto hueco se cava en cada caja: una cavidad profunda encierra muchos
    # mas atomos de receptor que una superficie plana del mismo tamaño.
    def atomos_dentro(c, s):
        return sum(1 for a in recept if all(abs(a[i + 3] - c[i]) <= s / 2 for i in range(3)))
    print("\nAtomos de receptor dentro de cada caja (mas = mas enterrado):")
    print("  caja del cribado: %d" % atomos_dentro(cv, sv))
    print("  caja del ARN:     %d" % atomos_dentro(centro, lade))

    # ── 5) Guardar la caja propuesta ─────────────────────────────
    datos = {
        "origen": "4BS2 (TDP-43 tandem RRM + ARN UG-rico), modelo 1, cadena A",
        "centro": [round(c, 3) for c in centro],
        "tamano": lade,
        "residuos_interfaz": ["%s%d" % (n, r) for r, n in por_num],
        "corte_contacto_A": CORTE_CONTACTO,
        "comparacion_caja_vieja": {
            "centro": list(cv), "tamano": sv,
            "distancia_entre_centros_A": round(dist_centros, 2),
            "porcentaje_interfaz_dentro": round(100.0 * dentro / len(interfaz), 1),
        },
        "nota": ("Centro = centro de masas de los atomos de los residuos que tocan "
                 "al ARN en 4BS2; tamaño = extension de esos atomos + 8 A de margen."),
    }
    os.makedirs(os.path.dirname(SALIDA), exist_ok=True)
    with open(SALIDA, "w", encoding="utf-8") as f:
        json.dump(datos, f, indent=2, ensure_ascii=False)
    print("\nGuardado: %s" % SALIDA)


if __name__ == "__main__":
    main()
