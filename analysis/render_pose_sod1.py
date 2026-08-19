# -*- coding: utf-8 -*-
"""Render 3D interactivo (3Dmol.js) de CHEMBL3311449 en el bolsillo anti-agregacion de SOD1."""
import os

POSE = r"C:/Users/Fredy/masive-als/analysis/_rescoring_stage/poses/CHEMBL3311449_SOD1_out.pdbqt"
REC  = r"C:/Users/Fredy/masive-als/analysis/receptor_SOD1.pdb"
OUT  = r"C:/Users/Fredy/proyectos-ia/preview_pose_sod1.html"

BOX_CENTER = (46.5, 80.0, 73.3)
BOX_SIZE   = 22.0
POCKET_CUT = 6.0  # Angstrom

# ---------- parsear ligando (solo MODEL 1) ----------
def parse_pdbqt_model1(path):
    atoms = []  # (x, y, z, element)
    in_model1 = False
    for line in open(path, "r", errors="ignore"):
        if line.startswith("MODEL"):
            in_model1 = (line.split()[1] == "1")
            continue
        if line.startswith("ENDMDL"):
            break
        if not in_model1:
            continue
        if line.startswith(("ATOM", "HETATM")):
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            atyp = line[76:79].strip() or line[-2:].strip()
            el = atyp[0] if atyp else "C"
            elmap = {"A": "C", "N": "N", "NA": "N", "O": "O", "OA": "O",
                     "S": "S", "SA": "S", "P": "P", "HD": "H", "H": "H"}
            el = elmap.get(el, el)
            atoms.append((x, y, z, el))
    return atoms

# ---------- parsear receptor ----------
def parse_pdb(path):
    atoms = []  # (x,y,z, resnum, resname)
    for line in open(path, "r", errors="ignore"):
        if line.startswith(("ATOM", "HETATM")):
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            except ValueError:
                continue
            resnum = int(line[22:26].strip() or 0)
            resname = line[17:20].strip()
            atoms.append((x, y, z, resnum, resname))
    return atoms

lig = parse_pdbqt_model1(POSE)
rec = parse_pdb(REC)
print("ligandos:", len(lig), "| receptor atoms:", len(rec))

# residuos del bolsillo (cualquier atomo a <= POCKET_CUT del ligando)
pocket = {}
for (x, y, z, _, resname) in rec:
    for (lx, ly, lz, _) in lig:
        d2 = (x-lx)**2 + (y-ly)**2 + (z-lz)**2
        if d2 <= POCKET_CUT**2:
            pocket[resname] = True
            break
pocket_res = sorted(r for r in pocket if r not in ("HOH", "WAT"))
print("residuos del bolsillo (<=6 A):", ", ".join(pocket_res))

# ---------- construir PDB del ligando ----------
pdb_lig = []
for i, (x, y, z, el) in enumerate(lig, 1):
    name = el.ljust(2)
    pdb_lig.append(f"ATOM  {i:5d} {name:>4s} UNL A   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {el:>2s}")
pdb_lig_s = "\n".join(pdb_lig)

rec_s = open(REC, "r", errors="ignore").read()

# ---------- HTML ----------
html = """<!doctype html>
<html lang="es"><head><meta charset="utf-8">
<title>Pose CHEMBL3311449 — SOD1</title>
<script src="https://cdn.jsdelivr.net/npm/3dmol/build/3Dmol-min.js"></script>
<style>
 body{background:#161618;color:#e6e6e6;font-family:Segoe UI,Arial,sans-serif;margin:0;padding:20px;text-align:center}
 h1{font-size:1.35rem;margin:6px 0}
 #viewer{width:100%;height:640px;background:#0f0f11;border-radius:12px;margin:14px auto}
 .info{color:#9aa0a6;font-size:0.95rem;max-width:820px;margin:0 auto 10px}
 table{margin:16px auto;border-collapse:collapse;color:#d0d0d0}
 td,th{border:1px solid #3a3a40;padding:5px 12px;font-size:0.88rem}
 th{background:#2a2a2e}
 .tag{display:inline-block;background:#252529;border:1px solid #3a3a40;border-radius:6px;padding:2px 8px;margin:2px;font-size:0.8rem}
</style></head>
<body>
<h1>CHEMBL3311449 dentro del bolsillo anti-agregación de SOD1</h1>
<div class="info">Proteína en azul (cartoon) · ligando en naranja (sticks) · residuos del bolsillo (≤ 6 Å) en gris ·
 caja de acoplamiento (46.5, 80.0, 73.3 Å, 22³) en amarillo semitransparente.
 <b>Gira con el ratón; zoom con la rueda.</b></div>
<div id="viewer"></div>
<table>
 <tr><th>Parámetro</th><th>Valor</th></tr>
 <tr><td>Ligando</td><td>CHEMBL3311449 (pirrolo[2,3-d]pirimidina)</td></tr>
 <tr><td>Diana</td><td>SOD1 — bolsillo Trp32 anti-agregación</td></tr>
 <tr><td>Vina (pose)</td><td>-6.9 kcal/mol</td></tr>
 <tr><td>MM-GBSA (rescoring)</td><td><b>-15.79 kcal/mol</b> (líder actual)</td></tr>
 <tr><td>Residuos del bolsillo</td><td>__POCKET__</td></tr>
</table>
<script>
const recPDB = `__REC__`;
const ligPDB = `__LIG__`;
const viewer = $3Dmol.createViewer("viewer", {backgroundColor:"0x0f0f11"});
viewer.addModel(recPDB, "pdb");
viewer.setStyle({}, {cartoon:{color:"0x4a90d9"}});
viewer.addModel(ligPDB, "pdb");
viewer.setStyle({model:1}, {stick:{radius:0.18, colorscheme:"orangeCarbon"}, sphere:{scale:0.25, colorscheme:"orangeCarbon"}});
viewer.addBox({center:{x:46.5,y:80.0,z:73.3}, dimensions:{w:22,h:22,d:22},
               color:"yellow", opacity:0.12});
viewer.zoomTo({model:1}, 380);
viewer.render();
</script>
</body></html>"""

html = html.replace("__REC__", rec_s.replace("\\", "\\\\").replace("`", "\\`"))
html = html.replace("__LIG__", pdb_lig_s.replace("\\", "\\\\").replace("`", "\\`"))
html = html.replace("__POCKET__", ", ".join(pocket_res))

with open(OUT, "w", encoding="utf-8") as f:
    f.write(html)
print("HTML ->", OUT, os.path.getsize(OUT), "bytes")
