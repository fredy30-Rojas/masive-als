# -*- coding: utf-8 -*-
"""Parche del paper: validación v4 TDP-43 (caja corregida en 4IUF) + hallazgo 6B1N->4IUF."""
import io

PATH = "C:/Users/Fredy/masive-als/paper/paper_masive_als.md"
with io.open(PATH, encoding="utf-8") as f:
    txt = f.read()

# 1) Actualizar el párrafo de la validación del 20/08 con la corrección 4IUF y la v4
old_parrafo = """On 20 Aug 2026 the benchmark was repeated with a **7× larger decoy set** (6,462 property-matched decoys drawn from the full screening library, versus 73 and 33 before) and with a **re-centered TDP-43 box on the RNA-binding surface** (center 24.4, 44.9, 53.4 Å; 18×18×18 Å), derived from the centroid of the RNP1/RNP2 aromatic residues that define the RNA-binding face of RRM1. The results were: TDP-43 ROC-AUC = 0.506 with 193 decoys (5 actives), FUS ROC-AUC = 0.609 with 87 decoys (2 actives). Neither box discriminates actives from decoys. Two findings stand out. First, the TDP-43 controls are heterogeneous and site-mismatched: rTRD01 and nTRD22 target the RNA-binding region, whereas bis-ANS binds the C-terminal low-complexity domain (outside RRM1), and 5-fluorouridine/isoproterenol are ligands crystallized at the SOD1 Trp32 pocket rather than TDP-43; a fraction of the "actives" were docked against a site they do not target, a mismatch that inherently depresses AUC. Second, the decoy distributions themselves show a **size artefact**: several large decoys score as low as −28 to −25 kcal/mol — far below the best active (−5.2) — a well-known Vina bias in which predicted affinity correlates with heavy-atom count when the box is small and ligand-sized molecules fill it. This biases any fixed-cutoff hit list toward large compounds and reinforces that the Vina score must be used as a **triage filter, not a ranker**, and that the shortlist must be re-scored with an orthogonal physics-based method (MM-GBSA) before prioritization."""

new_parrafo = """On 20 Aug 2026 the benchmark was repeated with a **7× larger decoy set** (6,462 property-matched decoys drawn from the full screening library, versus 73 and 33 before). This round uncovered a **documentation error in the reference structure**: the PDB entry cited for the TDP-43 receptor, 6B1N, is in fact Hsc70/HSPA8 (a heat-shock protein), not TDP-43; the docking receptor actually used is the RRM1 domain of TDP-43 (sequence-identical to PDB 4IUF, 773/786 atoms matched to the co-crystal). The native nucleic-acid ligand co-crystallized in 4IUF (residue XUA + DNA strand) marks the true RNA-binding surface, which lies ~17 Å from the previously used box center. Re-docking the native ligand fragment (XUA) against the re-centered box (center 16.3, 41.1, 48.5 Å; 24×24×24 Å) recovered the crystallographic pose with **RMSD = 1.75 Å** (< 2 Å, the accepted re-docking success threshold), validating the docking protocol itself. Re-running the decoy benchmark on the corrected box gave **TDP-43 ROC-AUC = 0.612** with 193 decoys (5 actives), up from 0.506 on the misplaced box; bis-ANS (a genuine TDP-43 ligand) scored −6.58 kcal/mol, better than 94.8% of decoys. The remaining gap to the SOD1 standard (0.815) is explained by the same two factors as before: the positive controls remain heterogeneous and partly site-mismatched (rTRD01/nTRD22 bind the RNA region while 5-fluorouridine/isoproterenol are SOD1 ligands), and a systematic **Vina size artefact** lets large decoys score as low as −28 to −25 kcal/mol — far below the best active (−5.2) — a well-known bias in which predicted affinity correlates with heavy-atom count when the box is small and ligand-sized molecules fill it. This reinforces that the Vina score must be used as a **triage filter, not a ranker**, and that the shortlist must be re-scored with an orthogonal physics-based method (MM-GBSA) before prioritization. The corrected TDP-43 box is now used to re-dock the previously screened tandas (Section 5.4)."""

if old_parrafo in txt:
    txt = txt.replace(old_parrafo, new_parrafo)
    print("OK: parrafo 20/08 actualizado")
else:
    print("AVISO: parrafo 20/08 no encontrado")

# 2) Actualizar la tabla de validación
old_tabla = "| TDP-43 — RRM1 RNA-binding face (18 Å box) | 5 | 193 | 0.506 | 0.0 | 0.0 |"
new_tabla = "| TDP-43 — RRM1 RNA-binding surface, 4IUF-native box (24 Å) | 5 | 193 | 0.612 | 0.0 | 0.0 |"
if old_tabla in txt:
    txt = txt.replace(old_tabla, new_tabla)
    print("OK: tabla actualizada")
else:
    print("AVISO: tabla no encontrada")

# 3) Actualizar interpretación
old_interp = "TDP-43 (0.506–0.515) and FUS (0.545–0.609), in contrast, do not yet discriminate actives from decoys, even after re-centering the TDP-43 box on the RNA-binding face and enlarging the decoy set 7-fold."
new_interp = "TDP-43 improved from 0.506 to 0.612 once the box was re-centered on the true RNA-binding surface defined by the native ligand of 4IUF (pose re-docking RMSD 1.75 Å), but remains below the SOD1 standard; FUS (0.545–0.609) does not yet discriminate actives from decoys."
if old_interp in txt:
    txt = txt.replace(old_interp, new_interp)
    print("OK: interpretacion actualizada")
else:
    print("AVISO: interpretacion no encontrada")

old_pend = "the TDP-43/FUS boxes will be re-benchmarked against site-matched positive controls or native co-crystal ligands via pose re-docking (RMSD) as soon as such controls become available."
new_pend = "the TDP-43 box has been re-benchmarked against the native co-crystal ligand of 4IUF via pose re-docking (RMSD 1.75 Å) and is now centered on the true RNA-binding surface; FUS still awaits site-matched controls."
if old_pend in txt:
    txt = txt.replace(old_pend, new_pend)
    print("OK: pendiente actualizado")
else:
    print("AVISO: pendiente no encontrado")

with io.open(PATH, "w", encoding="utf-8") as f:
    f.write(txt)
print("Paper guardado.")
