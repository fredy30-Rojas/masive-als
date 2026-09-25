# -*- coding: utf-8 -*-
"""Anade el argumento --tamano a validar_senuelos.py (tamano de caja en A)."""
import io

p = r"C:\Users\Fredy\masive-als\analysis\validar_senuelos.py"
with io.open(p, "r", encoding="utf-8") as f:
    txt = f.read()

# 1) variable global SIZE -> parametrizable
txt = txt.replace("SIZE = 25", "SIZE = 18  # por defecto; se puede cambiar con --tamano")

# 2) argparser: anadir --tamano
old_ap = 'ap.add_argument("--decoys-por-activo", type=int, default=30)'
new_ap = 'ap.add_argument("--decoys-por-activo", type=int, default=30)\n    ap.add_argument("--tamano", type=int, default=None, help="tamano de caja en Angstrom")'
assert old_ap in txt, "argparser no encontrado"
txt = txt.replace(old_ap, new_ap)

# 3) usar el tamano si se pasa (despues de leer --centro)
old_centro = '    cx, cy, cz = [float(x) for x in args.centro.split(",")]'
new_centro = ('    cx, cy, cz = [float(x) for x in args.centro.split(",")]\n'
              '    global SIZE\n'
              '    if args.tamano:\n'
              '        SIZE = args.tamano')
assert old_centro in txt, "centro no encontrado"
txt = txt.replace(old_centro, new_centro)

with io.open(p, "w", encoding="utf-8") as f:
    f.write(txt)
print("parche aplicado")
