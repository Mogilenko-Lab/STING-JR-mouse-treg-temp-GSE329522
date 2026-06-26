# 07_synthesis - capstone synthesis (two-arms panel + reply memo)

Stage 07 assembles the cross-arm evidence for the STING-cGAS standard sweep
into one tidy two-arms table and renders the headline figure. NO new
statistics are computed here - it joins/summarises the masters the earlier
arms produced (GSEA, PROGENy, decoupleR-TF, DE; GATOM/CoReSh where present).

Headline result (the publication-relevant payoff): a cGAS-dependence
ASYMMETRY. The IFN/ISG arm is cGAS-dependent (positive, significant
Interaction); the HIF/glycolysis arm rises in BOTH heat arms but is flat in
the Interaction - no detectable cGAS-dependence at n=5. This is an
asymmetry, NOT proven cGAS-independence (the 1-df interaction is the
lowest-powered comparison). HIF1a/HIF2a are NOT crowned as drivers; the arm
names are labels. Figure claims floor at L3 (DE/enrichment statistics);
mechanism (L7: pseudohypoxia / Complex-I) lives only in the reply memo.

PROVISIONAL - inferred sample mapping (Hspa1b/Hsph1 thermometer + Cgas),
pending the collaborator sample sheet; n=5/group.

Artifacts:
- figures/_overview/two_arms_panel.pdf + two_arms_panel.png - the headline panel
- tables/_overview/two_arms_panel.csv - the source table behind the panel
- tables/two_arms_summary.csv - the full cross-arm evidence table (16_synthesis.R)

Reply memo: docs/_internal/reports/2026-06-24_claim-evidence-memo.md

## figures/_overview/two_arms_panel.png

Two-arms cGAS-dependence asymmetry, multi-method: the IFN/ISG arm is
cGAS-dependent (positive, significant Interaction; 18 significant
method-features in the Interaction column), while the HIF/glycolysis
arm rises in BOTH WT_heat and KO_heat yet is flat in the Interaction
(1 significant; no detectable cGAS-dependence at n=5). Convergent
across GSEA, PROGENy, decoupleR-TF, and per-gene DE (plus
GATOM/CoReSh where provisioned). PROVISIONAL; n=5/group; NOT proven
independence.

**How to read:** Two stacked tracks: TOP = IFN/ISG arm (cGAS-dependent), BOTTOM =
HIF/glycolysis arm (no detectable cGAS-dependence at n=5). Rows
within a track = method-feature glyph rows (GSEA gene sets, PROGENy
pathways, decoupleR TFs, DE marker genes; GATOM modules where
present), ordered by Interaction score. Columns = headline contrasts
(WT heat | cGAS-KO heat | Interaction | Temp main). Tile fill =
signed score (orange = up in the numerator condition, blue = down),
clamped to +/-3.5; the printed number is the score. A BLACK RING
means padj < 0.05 (significant). READ THE ASYMMETRY DOWN THE
'Interaction' COLUMN: the IFN/ISG track lights up (positive, ringed)
= cGAS-dependent; the HIF/glycolysis track goes flat / unringed there
while staying lit in BOTH heat columns = no detectable
cGAS-dependence at n=5. The arm NAMES are labels, not claims that any
one TF (e.g. HIF1a/HIF2a) is the driver. Claim tier: L3
(DE/enrichment statistics; provisional, n=5/group). 'Flat
Interaction' is NOT proven cGAS-independence - the 1-df interaction
is the lowest-powered comparison.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_synthesis_viz.R` | `build_two_arms_panel` | `thresholds.gsea_fdr=0.05; figures.nes_cap=3.5; colors.diverging; design.contrasts` | `03_results/07_synthesis/tables/two_arms_summary.csv; 03_results/objects/16_synthesis.rds` |

## figures/_overview/two_arms_panel.pdf

Two-arms cGAS-dependence asymmetry, multi-method (print-resolution vector PDF):
the IFN/ISG arm is cGAS-dependent (positive, significant Interaction; 18
significant method-features in the Interaction column), while the HIF/glycolysis
arm rises in BOTH WT_heat and KO_heat yet is flat in the Interaction (1
significant; no detectable cGAS-dependence at n=5). Convergent across GSEA,
PROGENy, decoupleR-TF, and per-gene DE (plus GATOM/CoReSh where provisioned).
PROVISIONAL; n=5/group; NOT proven independence.

**How to read:** Vector PDF companion to `two_arms_panel.png` — identical data and
geometry, rendered at 18 in × 16 in on a capstone canvas sized to the
40-glyph-row × 4-column layout (expanded from the default 10×8 to prevent column
collapse and x-label overprinting). Two stacked tracks: TOP = IFN/ISG arm
(cGAS-dependent), BOTTOM = HIF/glycolysis arm (no detectable cGAS-dependence at
n=5). Rows within a track = method-feature glyph rows (GSEA gene sets, PROGENy
pathways, decoupleR TFs, DE marker genes; GATOM modules where present), ordered
by descending absolute Interaction score (strongest cGAS-dependence signal at
top); capped to top_n=20 rows per track by |Interaction score|. Columns =
headline contrasts: WT_heat | cGAS-KO heat | Interaction | Temp_main. Tile fill
= signed score (orange = up in the numerator condition, blue = down), clamped to
+/-3.5; the in-tile printed number is the raw score. A BLACK RING (open-square
border, linewidth 0.9) marks padj < 0.05 — the significance indicator is
color-blind and grayscale safe (no hue-only encoding). READ THE ASYMMETRY DOWN
THE 'Interaction' COLUMN: the IFN/ISG track lights up (positive, ringed tiles) =
cGAS-dependent; the HIF/glycolysis track goes flat / unringed there while staying
lit in BOTH heat columns = no detectable cGAS-dependence at n=5. The arm NAMES
(IFN/ISG, HIF/glycolysis) are labels, NOT claims that any one TF (e.g. HIF1a or
HIF2a) is the driver. HSF1 / heatshock_context rows are held out of this headline
panel (reply-memo context only). Sign convention: positive score = feature is
elevated in the numerator condition of each contrast (e.g., WT_heat numerator =
heated WT vs. vehicle WT). Claim tier: L3 (DE/enrichment statistics; provisional,
n=5/group). 'Flat Interaction' is NOT proven cGAS-independence — the 1-df
interaction is the lowest-powered comparison in a 2×2 design at n=5/group.
[PROVISIONAL — inferred sample mapping (Hspa1b/Hsph1 thermometer + Cgas),
pending collaborator sample sheet; n=5/group.]

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_synthesis_viz.R` | `build_two_arms_panel` | `thresholds.gsea_fdr=0.05; figures.nes_cap=3.5; colors.diverging; design.contrasts; width=18; height=16` | `03_results/07_synthesis/tables/two_arms_summary.csv; 03_results/objects/16_synthesis.rds` |

## tables/_overview/two_arms_panel.csv

Source table for the `two_arms_panel` headline figure: one row per
(arm, method, feature, contrast) combination displayed in the heatmap,
carrying all columns needed to reproduce or audit every tile.

**How to read:** Each row represents a single tile in the two-arms heatmap.
Columns: `arm` (IFN_ISG or HIF_glycolysis), `arm_track` (the pretty facet label
shown in the figure), `method` (source arm: GSEA, PROGENy, TF, DE, GATOM),
`feature` (gene set / pathway / TF / gene name), `feature_kind` (gene set
category or analyte type where available), `contrast` (one of: WT_heat |
KO_heat | Interaction | Temp_main), `score` (raw signed statistic before clamping
— GSEA NES, PROGENy/TF MLM-ULM score, or limma t-statistic), `pvalue`, `padj`
(Benjamini-Hochberg corrected), `direction` (up/down/ns character flag),
`significant` (logical; TRUE when padj < 0.05). The `score` column is the
unclamped value; the figure shows it clamped to +/-3.5. The cGAS-dependence test
lives in the `Interaction` contrast rows. Rows are ordered by arm / method /
feature / contrast. Claim tier: L3 (same as the figure it underlies; provisional,
n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_synthesis_viz.R` | `build_two_arms_panel` (panel_table transmute) | `thresholds.gsea_fdr=0.05; figures.nes_cap=3.5; design.contrasts` | `03_results/07_synthesis/tables/two_arms_summary.csv; 03_results/objects/16_synthesis.rds` |

## tables/two_arms_summary.csv

The full cross-arm evidence table: one row per (arm, method, feature,
contrast). Columns: arm, arm_label, arm_hypothesis, method, feature,
feature_kind, contrast, score, pvalue, padj, direction, significant.
Assembled by 16_synthesis.R from the arm masters (no new statistics).

**How to read:** `arm` in {IFN_ISG (cGAS-dependent hypothesis),
HIF_glycolysis (no-detectable-cGAS-dependence hypothesis), heatshock_context
(HSF1; not a cGAS arm)}. `method` names the source arm (GSEA:<db>, PROGENy,
TF:CollecTRI, DE:limma-trend, GATOM). `score` = the method's signed statistic
(GSEA NES / PROGENy-TF MLM-ULM score / limma t / GATOM pseudo-NES). `padj`
drives `significant` (padj < gsea_fdr). The cGAS-dependence test is the
`Interaction` contrast. Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/16_synthesis.R` | `build_*_arm + to_evidence` | `thresholds.gsea_fdr; design.contrasts; schemas.master_gsea_table` | `03_results/master/master_{gsea_table,progeny_activities,tf_activities,de_genes,gatom_modules}.csv` |

