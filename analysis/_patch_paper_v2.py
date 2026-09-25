# -*- coding: utf-8 -*-
"""Aplica al paper la validación v2/v3 de TDP-43 y FUS (20 ago 2026)."""
import io

P = r"C:\Users\Fredy\masive-als\paper\paper_masive_als.md"
with io.open(P, "r", encoding="utf-8") as f:
    txt = f.read()

# 1) párrafo principal de validación TDP43/FUS
old1 = ("were used against 33 decoys: ROC-AUC = 0.545, EF5% = 0.0. These values are "
        "statistically weak (small positive-control sets, particularly FUS, where the "
        "literature currently offers few direct small-molecule binders) and indicate that "
        "the TDP-43 RRM1 box and the FUS model-1 box do **not** currently discriminate "
        "known binders from decoys, unlike the validated SOD1 Trp32 site. The TDP-43 and "
        "FUS rankings should therefore be interpreted with caution until their binding "
        "sites are re-benchmarked (e.g., against RNA-competitive ligands or alternative "
        "pockets), mirroring the SOD1 correction; candidate lists from these two targets "
        "are provisional.")
new1 = ("were used against 33 decoys: ROC-AUC = 0.545, EF5% = 0.0.\n"
        "\n"
        "On 20 Aug 2026 the benchmark was repeated with a **7× larger decoy set** (6,462 "
        "property-matched decoys drawn from the full screening library, versus 73 and 33 "
        "before) and with a **re-centered TDP-43 box on the RNA-binding surface** (center "
        "24.4, 44.9, 53.4 Å; 18×18×18 Å), derived from the centroid of the RNP1/RNP2 "
        "aromatic residues that define the RNA-binding face of RRM1. The results were: "
        "TDP-43 ROC-AUC = 0.506 with 193 decoys (5 actives), FUS ROC-AUC = 0.609 with 87 "
        "decoys (2 actives). Neither box discriminates actives from decoys. Two findings "
        "stand out. First, the TDP-43 controls are heterogeneous and site-mismatched: "
        "rTRD01 and nTRD22 target the RNA-binding region, whereas bis-ANS binds the "
        "C-terminal low-complexity domain (outside RRM1), and 5-fluorouridine/"
        "isoproterenol are ligands crystallized at the SOD1 Trp32 pocket rather than "
        "TDP-43; a fraction of the \"actives\" were docked against a site they do not "
        "target, a mismatch that inherently depresses AUC. Second, the decoy distributions "
        "themselves show a **size artefact**: several large decoys score as low as −28 to "
        "−25 kcal/mol — far below the best active (−5.2) — a well-known Vina bias in which "
        "predicted affinity correlates with heavy-atom count when the box is small and "
        "ligand-sized molecules fill it. This biases any fixed-cutoff hit list toward large "
        "compounds and reinforces that the Vina score must be used as a **triage filter, "
        "not a ranker**, and that the shortlist must be re-scored with an orthogonal "
        "physics-based method (MM-GBSA) before prioritization.\n"
        "\n"
        "These values are statistically weak (small positive-control sets, particularly FUS "
        "and TDP-43, where the literature currently offers few direct small-molecule binders "
        "of the specific domain modelled) and indicate that the TDP-43 RRM1 box and the FUS "
        "model-1 box do **not** currently discriminate known binders from decoys, unlike the "
        "validated SOD1 Trp32 site. A ChEMBL and PubChem query on 20 Aug 2026 confirmed that "
        "**no curated bioactivity records exist** for TDP-43 (ChEMBL target CHEMBL2362981) "
        "or FUS (CHEMBL5724679), so there are no public positive controls with which to "
        "properly calibrate these two targets today. The TDP-43 and FUS rankings should "
        "therefore be interpreted with caution until site-matched positive controls become "
        "available (e.g., RNA-competitive ligands with reported RRM1 Kd values); candidate "
        "lists from these two targets are provisional.")
assert old1 in txt, "old1 no encontrado"
txt = txt.replace(old1, new1, 1)

# 2) filas de la Tabla 1
old2 = ("| TDP-43 — RRM1 domain | 5 | 73 | 0.515 | 0.0 | 0.0 |\n"
        "| FUS — NMR model 1 (near-blind box) | 2 | 33 | 0.545 | 0.0 | 0.0 |")
new2 = ("| TDP-43 — RRM1 domain (25 Å box) | 5 | 73 | 0.515 | 0.0 | 0.0 |\n"
        "| TDP-43 — RRM1 RNA-binding face (18 Å box) | 5 | 193 | 0.506 | 0.0 | 0.0 |\n"
        "| FUS — NMR model 1 (near-blind box) | 2 | 33 | 0.545 | 0.0 | 0.0 |\n"
        "| FUS — NMR model 1 (re-run, 87 decoys) | 2 | 87 | 0.609 | 0.0 | 0.0 |")
assert old2 in txt, "old2 no encontrado"
txt = txt.replace(old2, new2, 1)

# 3) interpretación de la Tabla 1
old3 = ("TDP-43 (0.515) and FUS (0.545), in contrast, do not yet discriminate actives "
        "from decoys. Two caveats temper the latter: (i) the positive-control sets are small "
        "(5 and 2 docked actives; one TDP-43 control, Congo Red, failed PDBQT conversion and "
        "was excluded), so the AUC is noisy; and (ii) the TDP-43 controls are heterogeneous "
        "— rTRD01, nTRD22 and 5-fluorouridine target the RNA-binding region, whereas bis-ANS "
        "binds the C-terminal low-complexity domain, outside the RRM1 box — so a fraction of "
        "the \"actives\" were docked against a site they do not target, a mismatch that "
        "inherently depresses AUC. We therefore treat the Vina score as a **triage filter, "
        "not a ranker**: the per-target top-5% cut reduces the ~2M-compound library to a "
        "manageable candidate set, but the *order* within that set is not claimed to predict "
        "activity. Before any in vitro prioritization, the shortlisted candidates will be "
        "re-scored with an orthogonal physics-based method (MM-GBSA), and the TDP-43/FUS "
        "boxes will be re-benchmarked against site-matched positive controls or native "
        "co-crystal ligands via pose re-docking (RMSD). Until then, the SOD1 candidate list "
        "is the highest-confidence of the three; the TDP-43 and FUS lists remain provisional.")
new3 = ("TDP-43 (0.506–0.515) and FUS (0.545–0.609), in contrast, do not yet discriminate "
        "actives from decoys, even after re-centering the TDP-43 box on the RNA-binding face "
        "and enlarging the decoy set 7-fold. Three factors explain this: (i) the "
        "positive-control sets are tiny and, for TDP-43, site-mismatched (5 and 2 docked "
        "actives; one TDP-43 control, Congo Red, failed PDBQT conversion and was excluded), "
        "so the AUC is noisy; (ii) a systematic Vina size artefact lets large decoys score as "
        "low as −28 kcal/mol, inflating false positives at fixed cutoffs; and (iii) no curated "
        "public bioactivity data exist for these two targets, so a properly calibrated "
        "positive-control set is not currently available. We therefore treat the Vina score as "
        "a **triage filter, not a ranker**: the per-target top-5% cut reduces the ~2M-compound "
        "library to a manageable candidate set, but the *order* within that set is not claimed "
        "to predict activity. Before any in vitro prioritization, the shortlisted candidates "
        "are re-scored with an orthogonal physics-based method (MM-GBSA), and the TDP-43/FUS "
        "boxes will be re-benchmarked against site-matched positive controls or native "
        "co-crystal ligands via pose re-docking (RMSD) as soon as such controls become "
        "available. Until then, the SOD1 candidate list is the highest-confidence of the "
        "three; the TDP-43 and FUS lists remain provisional.")
assert old3 in txt, "old3 no encontrado"
txt = txt.replace(old3, new3, 1)

with io.open(P, "w", encoding="utf-8") as f:
    f.write(txt)
print("paper actualizado OK")
