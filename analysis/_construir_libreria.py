# Construye full_library_solo.smi consolidando todas las fuentes de compounds/
import os, csv, glob

BASE = r"C:\Users\Fredy\masive-als\compounds"
SALIDA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "full_library_solo.smi")

vistos = {}  # smiles -> nombre
def agregar(smiles, nombre):
    s = (smiles or "").strip()
    if not s or s.startswith("#"):
        return
    if s not in vistos:
        vistos[s] = (nombre or s).strip() or s

# TSV: smiles<TAB>nombre
for tsv in ("seleccion_30000.tsv", "seleccion_30000_v2.tsv"):
    p = os.path.join(BASE, tsv)
    if not os.path.exists(p):
        continue
    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            # formato: CHEMBL_ID<TAB>SMILES
            if len(parts) >= 2 and parts[0].strip() and parts[1].strip():
                agregar(parts[1], parts[0])
            elif parts and parts[0].strip():
                agregar(parts[0], "")
    print(tsv, "leido")

# CSV con columna smiles
for csvf in glob.glob(os.path.join(BASE, "*.csv")):
    try:
        with open(csvf, encoding="utf-8", errors="replace") as f:
            rd = csv.DictReader(f)
            cols = rd.fieldnames or []
            col = next((c for c in ("smiles", "canonical_smiles") if c in cols), None)
            if col is None:
                continue
            n = 0
            for r in rd:
                ident = r.get("name") or r.get("chembl_id") or r.get("zinc_id") or ""
                agregar(r[col], ident)
                n += 1
            print(os.path.basename(csvf), n, "filas")
    except Exception as e:
        print(os.path.basename(csvf), "error:", e)

# full_library.smi viejo (si existe)
p = os.path.join(BASE, "full_library.smi")
if os.path.exists(p):
    with open(p, encoding="utf-8", errors="replace") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            agregar(parts[0], parts[1] if len(parts) > 1 else "")
    print("full_library.smi leido")

with open(SALIDA, "w", encoding="utf-8") as f:
    for s, nombre in vistos.items():
        f.write("%s\t%s\n" % (s, nombre))
print("TOTAL compuestos unicos:", len(vistos))
print("->", SALIDA)
