# -*- coding: utf-8 -*-
"""Corrige el paper MASIVE-ALS: fuente TDP-43 = 4IUF (no 6B1N, que es Hsc70)
y caja corregida centrada en el sitio de union a RNA nativo."""
import io

P = r"C:\Users\Fredy\masive-als\paper\paper_masive_als.md"
with io.open(P, "r", encoding="utf-8") as f:
    txt = f.read()

# 1) resumen: TDP-43: 6B1N/4IUF -> 4IUF
old1 = "target structures obtained from the Protein Data Bank (TDP-43: 6B1N/4IUF, clean RRM1 domain; SOD1: 1HL5, anti-aggregation site; FUS: 6G99)"
new1 = "target structures obtained from the Protein Data Bank (TDP-43: 4IUF, clean RRM1 domain; SOD1: 1HL5, anti-aggregation site; FUS: 6G99)"
assert old1 in txt, "old1 no encontrado"
txt = txt.replace(old1, new1, 1)

# 2) metodos: TDP-43 (PDB 6B1N/4IUF, Source PDB) -> correccion
old2 = "Target structures were obtained from the Protein Data Bank: TDP-43 (PDB 6B1N/4IUF, Source PDB), SOD1 (PDB 1HL5), and FUS (PDB 6G99, an NMR ensemble; the first model was retained after repair with Open Babel)."
new2 = ("Target structures were obtained from the Protein Data Bank: TDP-43 (PDB 4IUF, the RRM1 domain in complex with DNA), SOD1 (PDB 1HL5), "
        "and FUS (PDB 6G99, an NMR ensemble; the first model was retained after repair with Open Babel). "
        "**Correction of 20 Aug 2026:** an earlier version of this document cited 6B1N as a TDP-43 structure; "
        "6B1N is in fact Hsc70/HSPA8 (a heat-shock chaperone) and was erroneous. The receptor actually used for "
        "docking is the clean RRM1 domain of TDP-43, which matches chain A of 4IUF in 773/786 heavy atoms "
        "(see `proteins/TDP43/metadata.json`).")
assert old2 in txt, "old2 no encontrado"
txt = txt.replace(old2, new2, 1)

# 3) nota de reproducibilidad: anadir caja corregida
old3 = ("**Important reproducibility note:** the receptor actually used for docking against TDP-43 was the clean "
        "RRM1 domain (`TDP43_RRM1_clean.pdb`, derived from the full PDB entry), not the full-length structure; "
        "this is recorded in `proteins/TDP43/metadata.json` and in the `REMARK Name` header of `gpu_dock/TDP43.pdbqt`.")
new3 = ("**Important reproducibility note:** the receptor actually used for docking against TDP-43 was the clean "
        "RRM1 domain (`TDP43_RRM1_clean.pdb`, derived from the full PDB entry), not the full-length structure; "
        "this is recorded in `proteins/TDP43/metadata.json` and in the `REMARK Name` header of `gpu_dock/TDP43.pdbqt`. "
        "**Box correction of 20 Aug 2026:** the original TDP-43 grid box (center 28.3, 43.7, 52.5 Å; 25×25×25 Å) "
        "was centred ~17 Å away from the co-crystallized DNA of 4IUF (center 16.3, 41.1, 48.5 Å), i.e. outside the "
        "RNA-binding surface defined by residues 109-179 (RNP1/RNP2). Re-docking of the native nucleotide XUA with "
        "the corrected box yields RMSD 1.75 Å (< 2 Å, successful pose reproduction), demonstrating the docking "
        "protocol is valid for TDP-43 when the box is correctly placed. Tandas screened with the old box are being "
        "re-docked with the corrected box; see `analysis/INVESTIGACION_TDP43_2026-08-20.md`.")
assert old3 in txt, "old3 no encontrado"
txt = txt.replace(old3, new3, 1)

# 4) parametros de docking: actualizar la caja TDP43
old4 = ("TDP-43 (TDP43_RRM1 clean domain; center 28.3, 43.7, 52.5 Å; 25×25×25 Å), "
        "SOD1 (center 46.5, 80.0, 73.3 Å; 22×22×22 Å)")
new4 = ("TDP-43 (TDP43_RRM1 clean domain; original box center 28.3, 43.7, 52.5 Å, 25×25×25 Å; corrected box "
        "center 16.3, 41.1, 48.5 Å, 24×24×24 Å centred on the co-crystallized DNA of 4IUF, applied from 20 Aug 2026), "
        "SOD1 (center 46.5, 80.0, 73.3 Å; 22×22×22 Å)")
assert old4 in txt, "old4 no encontrado"
txt = txt.replace(old4, new4, 1)

with io.open(P, "w", encoding="utf-8") as f:
    f.write(txt)
print("paper corregido con 4IUF y caja corregida")
