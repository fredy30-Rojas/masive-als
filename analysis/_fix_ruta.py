import io
p = r"C:\Users\Fredy\masive-als\analysis\_redock_top_tdp43.py"
t = io.open(p, encoding="utf-8").read()
old = r'LOTES = r"C:\Users\Fredy\masive-als\compounds\lotes_100_pdbqt"'
new = r'LOTES = r"C:\Users\Fredy\masive-als\gpu_dock\tanda_z001\ligands"'
assert old in t, "no encontrado"
io.open(p, "w", encoding="utf-8").write(t.replace(old, new, 1))
print("ruta corregida")
