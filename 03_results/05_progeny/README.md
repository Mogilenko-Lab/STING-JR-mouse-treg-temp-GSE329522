# 05_progeny: artifact captions

PROGENy MLM pathway activity inference (14 pathways, all contrasts)
+ CollecTRI ULM TF activity (decoupleR). Key result: Hypoxia (PROGENy) is flat in
the Interaction contrast, so no cGAS-dependence is detectable at n=5;
JAK-STAT/NFkB/TNFa are positive in the Interaction (cGAS-dependent arm).
Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

Figures use the unified style contract: one `.pdf` (print) + one `.png` (screen) per
stem, no `.screen`/`.print` suffix. Heatmap x-axis contrast labels are drawn at
30-degree angle (`hjust = 1`) to prevent overlap in dense multi-column layouts.

---

## figures/_overview/progeny_heatmap.png

PROGENy 14-pathway x 7-contrast activity heatmap (MLM score, clamped
+/-2.5; * raw p<0.05; rows clustered): Hypoxia/glycolysis rise in
BOTH heat arms but are flat in the Interaction; JAK-STAT/NFkB/TNFa
rise in WT_heat and are positive in the Interaction (cGAS-dependent).
n=5/group. Sample-to-condition mapping confirmed against the owner's
sample sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

**How to read:** Rows = PROGENy pathways (hierarchically clustered); columns =
contrasts in design order (left = WT_heat, KO_heat; then Interaction,
Temp_main, Geno*, Geno_main). Fill: orange = pathway activated in
numerator; blue = activated in denominator. Score clamped to +/-2.5;
* = raw p < 0.05 (n=14 pathways; no multi-test correction warranted).
Claim tier: L3 (n=5/group). Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.heat_tile (ggplot geom_tile)` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/_overview/progeny_interaction_split.png

PROGENy Hypoxia vs immune split: Hypoxia activity is similar in
WT_heat and KO_heat with a flat Interaction (no detectable
cGAS-dependence at n=5), whereas JAK-STAT/NFkB/TNFa are positive in
the Interaction (cGAS-dependent arm). This orthogonal pathway
footprint corroborates the two-arms DE/TF result. n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

**How to read:** Grouped bars: x = contrast (3 headline contrasts: WT_heat, KO_heat,
Interaction), fill = PROGENy pathway. Orange (Hypoxia) bars should be
similar in WT_heat and KO_heat and near-zero in the Interaction panel
- meaning no detectable cGAS-dependence. Blue bars
(JAK-STAT/NFkB/TNFa) should be taller in WT_heat than KO_heat and
positive in Interaction - cGAS-dependent. * = raw p < 0.05. CAUTION:
'flat Interaction' is NOT proven cGAS-independence; the study is
powered at n=5/group. Claim tier: L3. Sample-to-condition mapping
confirmed against the owner's sample sheet (2026-07-22): 20 of 20
libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `ggplot/geom_col` | `thresholds.gsea_fdr=0.05; colors.diverging; design.contrasts` | `03_results/objects/09_progeny_activity.rds` |

## figures/_overview/tf_heatmap.png

CollecTRI TF activity heatmap (TFs significant in >=2 contrasts +
watchlist; rows = TFs, cols = contrasts; score clamped +/-2.5; * BH
padj<0.05; row strip = HIF/IFN/other axis): IFN/IRF/STAT TFs cluster
as the cGAS-dependent block (positive in Interaction); HIF-axis TFs
are non-significant in the Interaction. n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

**How to read:** Rows = CollecTRI TFs (hierarchically clustered); columns = contrasts.
Fill: orange = TF activated in numerator; blue = activated in
denominator. Left-edge strip: orange = HIF axis (Hif1a/Epas1), blue =
IFN/NFkB axis, grey = other. Score clamped to +/-2.5. * = BH padj <
0.05. Claim tier: L3 (n=5/group; IFN-arm TFs positive in Interaction
= cGAS-dependent; HIF-axis TFs flat/NS in Interaction = no detectable
cGAS-dependence at n=5; NOT proven independence).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.heat_tile (ggplot geom_tile)` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/_overview/progeny_tf_combined.png

Combined PROGENy (top panel) + CollecTRI TF (bottom panel) activity
for key pathways/TFs across 3 headline contrasts: Hypoxia (PROGENy)
and Hif1a (TF) are active in both heat arms and flat in the
Interaction (no detectable cGAS-dependence at n=5);
JAK-STAT/IFN/Stat1/Irf3 are positive in the Interaction
(cGAS-dependent). Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

**How to read:** Row facets = PROGENy (top) vs CollecTRI TF (bottom); column facets =
headline contrasts. Lollipop color: orange = HIF axis, blue =
IFN/NFkB axis. Shape: diamond = PROGENy pathway, circle = TF.
Direction: rightward = activated in the numerator condition. * = p <
0.05 (raw for PROGENy; BH padj for TF). Key comparison: Hypoxia/Hif1a
bars are near-zero in the Interaction column - no detectable
cGAS-dependence. JAK-STAT/Stat1/Irf3 bars are positive in the
Interaction column - cGAS-dependent. CAUTION: n=5/group; absence of
significant Interaction for HIF arm is NOT proven independence. Claim
tier: L3. Sample-to-condition mapping confirmed against the owner's
sample sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `ggplot/facet_grid` | `thresholds.gsea_fdr=0.05; colors.diverging; design.contrasts` | `03_results/objects/09_progeny_activity.rds; 03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/WT_heat/progeny_barplot.png

PROGENy 14-pathway activity (MLM) in WT_heat: bars sorted by score; *
raw p<0.05; key pathways (Hypoxia/JAK-STAT/NFkB/TNFa) bold-outlined.

**How to read:** Bars = MLM activity score (orange = pathway more active in numerator;
blue = more active in denominator). Glyphs: * = raw p < 0.05; bold
outline = key pathway. Score is not a fold-change; sign tracks
numerator activation direction. Claim tier: L3 (n=5/group).
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/WT_heat/tf_barplot.png

Top CollecTRI TF activity (ULM) in WT_heat: lollipops colored by
HIF/IFN axis; * BH padj<0.05; IFN/HIF watchlist bold-outlined.

**How to read:** Lollipops = ULM activity score (rightward = pathway-level activation
in numerator). Color: orange = HIF axis (Hif1a/Epas1), blue =
IFN/NFkB axis, grey = other. Open circle = watchlist TF. * = BH padj
< 0.05. Claim tier: L3 (n=5/group). Sample-to-condition mapping
confirmed against the owner's sample sheet (2026-07-22): 20 of 20
libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/KO_heat/progeny_barplot.png

PROGENy 14-pathway activity (MLM) in KO_heat: bars sorted by score; *
raw p<0.05; key pathways (Hypoxia/JAK-STAT/NFkB/TNFa) bold-outlined.

**How to read:** Bars = MLM activity score (orange = pathway more active in numerator;
blue = more active in denominator). Glyphs: * = raw p < 0.05; bold
outline = key pathway. Score is not a fold-change; sign tracks
numerator activation direction. Claim tier: L3 (n=5/group).
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/KO_heat/tf_barplot.png

Top CollecTRI TF activity (ULM) in KO_heat: lollipops colored by
HIF/IFN axis; * BH padj<0.05; IFN/HIF watchlist bold-outlined.

**How to read:** Lollipops = ULM activity score (rightward = pathway-level activation
in numerator). Color: orange = HIF axis (Hif1a/Epas1), blue =
IFN/NFkB axis, grey = other. Open circle = watchlist TF. * = BH padj
< 0.05. Claim tier: L3 (n=5/group). Sample-to-condition mapping
confirmed against the owner's sample sheet (2026-07-22): 20 of 20
libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/Interaction/progeny_barplot.png

PROGENy 14-pathway activity (MLM) in Interaction: bars sorted by
score; * raw p<0.05; key pathways (Hypoxia/JAK-STAT/NFkB/TNFa)
bold-outlined.

**How to read:** Bars = MLM activity score (orange = pathway more active in numerator;
blue = more active in denominator). Glyphs: * = raw p < 0.05; bold
outline = key pathway. Score is not a fold-change; sign tracks
numerator activation direction. Claim tier: L3 (n=5/group).
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/Interaction/tf_barplot.png

Top CollecTRI TF activity (ULM) in Interaction: lollipops colored by
HIF/IFN axis; * BH padj<0.05; IFN/HIF watchlist bold-outlined.

**How to read:** Lollipops = ULM activity score (rightward = pathway-level activation
in numerator). Color: orange = HIF axis (Hif1a/Epas1), blue =
IFN/NFkB axis, grey = other. Open circle = watchlist TF. * = BH padj
< 0.05. Claim tier: L3 (n=5/group). Sample-to-condition mapping
confirmed against the owner's sample sheet (2026-07-22): 20 of 20
libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/Temp_main/progeny_barplot.png

PROGENy 14-pathway activity (MLM) in Temp_main: bars sorted by score;
* raw p<0.05; key pathways (Hypoxia/JAK-STAT/NFkB/TNFa)
bold-outlined.

**How to read:** Bars = MLM activity score (orange = pathway more active in numerator;
blue = more active in denominator). Glyphs: * = raw p < 0.05; bold
outline = key pathway. Score is not a fold-change; sign tracks
numerator activation direction. Claim tier: L3 (n=5/group).
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/Temp_main/tf_barplot.png

Top CollecTRI TF activity (ULM) in Temp_main: lollipops colored by
HIF/IFN axis; * BH padj<0.05; IFN/HIF watchlist bold-outlined.

**How to read:** Lollipops = ULM activity score (rightward = pathway-level activation
in numerator). Color: orange = HIF axis (Hif1a/Epas1), blue =
IFN/NFkB axis, grey = other. Open circle = watchlist TF. * = BH padj
< 0.05. Claim tier: L3 (n=5/group). Sample-to-condition mapping
confirmed against the owner's sample sheet (2026-07-22): 20 of 20
libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## tables/progeny_activity.csv

Flat table of raw decoupleR MLM output for all 14 PROGENy pathways across all
contrasts in the study design; used as the primary checkpoint for downstream
visualization (13_activity_viz.R reads the derived master CSV; this root-level
CSV is the serialized output from the compute script 09_activity_progeny.R).

**How to read:** Rows = one record per pathway × contrast combination. Columns include `source` (PROGENy pathway name), `contrast` (model contrast identifier), `score` (decoupleR MLM activity score — NOT a GSEA-NES), `p_value` (raw), and `padj` (BH-adjusted within contrast where available). Sign convention: positive score = pathway gene program more active in numerator condition. This table is the L2/L3 checkpoint object; it is never recomputed by the viz script. Claim tier: L3 (statistical inference from MLM regulon model).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `NOT_TRACED (written by 09_activity_progeny.R)` | `thresholds.gsea_fdr=0.05` | `03_results/objects/09_progeny_activity.rds` |

## tables/_overview/progeny_heatmap.csv

Source table for `figures/_overview/progeny_heatmap.pdf` and
`figures/_overview/progeny_heatmap.png`; long-format activity scores and raw
p-values for all 14 PROGENy pathways × all contrasts, after score clamping to ±2.5.

**How to read:** Columns: `source` (pathway), `contrast`, `score` (clamped MLM activity score — NOT a GSEA-NES), `p_value` (raw), `significant` (logical; raw p < 0.05). One row per pathway × contrast. Same-stem neighbor of the heatmap figure (source-table adjacency contract). Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.heat_tile / save_overview` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5` | `03_results/objects/09_progeny_activity.rds` |

## tables/_overview/progeny_interaction_split.csv

Source table for `figures/_overview/progeny_interaction_split.pdf` and
`figures/_overview/progeny_interaction_split.png`; long-format activity scores for
the four key PROGENy pathways (Hypoxia, JAK-STAT, NFkB, TNFa) across 3 headline
contrasts (WT_heat, KO_heat, Interaction).

**How to read:** Columns: `pathway` (PROGENy pathway name), `contrast`, `score` (MLM activity score — NOT a GSEA-NES), `p_value` (raw), `significant` (logical; raw p < 0.05). One row per pathway × contrast. Same-stem neighbor of the split bar chart figure. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `ggplot/geom_col / save_overview` | `thresholds.gsea_fdr=0.05; colors.diverging; design.contrasts` | `03_results/objects/09_progeny_activity.rds` |

## tables/_overview/tf_heatmap.csv

Source table for `figures/_overview/tf_heatmap.pdf` and
`figures/_overview/tf_heatmap.png`; long-format ULM activity scores, BH padj,
axis-family label, and significance flag for TFs significant in ≥2 contrasts
(plus watchlist) across all contrasts.

**How to read:** Columns: `source` (TF gene symbol), `axis` (HIF/IFN/other axis family), `contrast`, `score` (clamped ULM activity score — NOT a GSEA-NES), `padj` (BH-adjusted p-value within contrast), `significant` (logical; BH padj < 0.05). One row per TF × contrast. Same-stem neighbor of the TF heatmap figure. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.heat_tile / save_overview` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; figures.top_n=20` | `03_results/objects/03_tf_collectri.rds` |

## tables/_overview/progeny_tf_combined.csv

Source table for `figures/_overview/progeny_tf_combined.pdf` and
`figures/_overview/progeny_tf_combined.png`; unified long-format activity scores
for key PROGENy pathways (Hypoxia, JAK-STAT, NFkB, TNFa) and key TFs (Irf3, Irf7,
Irf1, Stat1, Stat2, Nfkb1, Rela, Hif1a, Epas1) across 3 headline contrasts (WT_heat,
KO_heat, Interaction).

**How to read:** Columns: `entity_type` (PROGENy or TF), `source` (pathway or TF gene symbol), `axis` (HIF/IFN/other), `contrast`, `score` (MLM or ULM activity score — NOT a GSEA-NES), `p` (raw p for PROGENy; BH padj for TF), `significant` (logical; p < 0.05). One row per entity × contrast. Same-stem neighbor of the combined lollipop panel. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `ggplot/facet_grid / save_overview` | `thresholds.gsea_fdr=0.05; colors.diverging; design.contrasts` | `03_results/objects/09_progeny_activity.rds; 03_results/objects/03_tf_collectri.rds` |

## tables/by_contrast/WT_heat/progeny_barplot.csv

Source table for `figures/by_contrast/WT_heat/progeny_barplot.pdf` and
`figures/by_contrast/WT_heat/progeny_barplot.png`; all 14 PROGENy pathway activity
scores and p-values for the WT_heat contrast, sorted by ascending raw p-value.

**How to read:** Columns: `source` (pathway), `contrast`, `score` (MLM activity score — NOT a GSEA-NES), `p_value` (raw), `padj` (BH within contrast). One row per pathway. Claim tier: L3 — provisional, n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny / save_overview` | `thresholds.gsea_fdr=0.05` | `03_results/objects/09_progeny_activity.rds` |

## tables/by_contrast/WT_heat/tf_barplot.csv

Source table for `figures/by_contrast/WT_heat/tf_barplot.pdf` and
`figures/by_contrast/WT_heat/tf_barplot.png`; top TFs by BH padj plus watchlist
force-inclusions for the WT_heat contrast.

**How to read:** Columns: `source` (TF gene symbol), `contrast`, `score` (ULM activity score — NOT a GSEA-NES), `p_value` (raw), `padj` (BH within contrast). Rows sorted by ascending padj, then watchlist additions. Claim tier: L3 — provisional, n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf / save_overview` | `thresholds.gsea_fdr=0.05; figures.top_n=20` | `03_results/objects/03_tf_collectri.rds` |

## tables/by_contrast/KO_heat/progeny_barplot.csv

Source table for `figures/by_contrast/KO_heat/progeny_barplot.pdf` and
`figures/by_contrast/KO_heat/progeny_barplot.png`; all 14 PROGENy pathway activity
scores and p-values for the KO_heat contrast.

**How to read:** Same column structure as `WT_heat/progeny_barplot.csv`. Compare score magnitudes with WT_heat to identify cGAS-dependent attenuation of pathway responses. Claim tier: L3 — provisional, n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny / save_overview` | `thresholds.gsea_fdr=0.05` | `03_results/objects/09_progeny_activity.rds` |

## tables/by_contrast/KO_heat/tf_barplot.csv

Source table for `figures/by_contrast/KO_heat/tf_barplot.pdf` and
`figures/by_contrast/KO_heat/tf_barplot.png`; top TFs by BH padj plus watchlist
force-inclusions for the KO_heat contrast.

**How to read:** Same column structure as `WT_heat/tf_barplot.csv`. Compare with WT_heat scores for IFN-axis TFs (Stat1, Irf3, Irf7) to see reduction consistent with cGAS loss. Claim tier: L3 — provisional, n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf / save_overview` | `thresholds.gsea_fdr=0.05; figures.top_n=20` | `03_results/objects/03_tf_collectri.rds` |

## tables/by_contrast/Interaction/progeny_barplot.csv

Source table for `figures/by_contrast/Interaction/progeny_barplot.pdf` and
`figures/by_contrast/Interaction/progeny_barplot.png`; all 14 PROGENy pathway
scores and p-values for the Interaction contrast (WT_heat minus KO_heat
differential; positive = cGAS-dependent).

**How to read:** Same column structure as other `progeny_barplot.csv` files. Positive `score` = cGAS-dependent pathway activation. Key result: Hypoxia score near zero (no detectable cGAS-dependence, n=5); JAK-STAT/NFkB/TNFa positive (cGAS-dependent). Claim tier: L3 — provisional, n=5/group; underpowered for HIF interaction.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny / save_overview` | `thresholds.gsea_fdr=0.05` | `03_results/objects/09_progeny_activity.rds` |

## tables/by_contrast/Interaction/tf_barplot.csv

Source table for `figures/by_contrast/Interaction/tf_barplot.pdf` and
`figures/by_contrast/Interaction/tf_barplot.png`; top TFs by BH padj plus watchlist
for the Interaction contrast (positive score = cGAS-dependent TF activity).

**How to read:** Same column structure as other `tf_barplot.csv` files. Positive `score` = TF activity requires intact cGAS. Key result: Stat1/Irf3/Irf7 positive and significant (cGAS-dependent); Hif1a/Epas1 near-zero (no detectable cGAS-dependence; NOT proven independence — n=5 underpowered). Claim tier: L3 — provisional, n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf / save_overview` | `thresholds.gsea_fdr=0.05; figures.top_n=20` | `03_results/objects/03_tf_collectri.rds` |

## tables/by_contrast/Temp_main/progeny_barplot.csv

Source table for `figures/by_contrast/Temp_main/progeny_barplot.pdf` and
`figures/by_contrast/Temp_main/progeny_barplot.png`; all 14 PROGENy pathway
activity scores and p-values for the Temp_main contrast (average hyperthermia
effect across WT and KO genotypes).

**How to read:** Same column structure as other `progeny_barplot.csv` files. Positive `score` = pathway activated by hyperthermia (39 C vs. 37 C), averaged across both genotypes. Claim tier: L3 — provisional, n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny / save_overview` | `thresholds.gsea_fdr=0.05` | `03_results/objects/09_progeny_activity.rds` |

## tables/by_contrast/Temp_main/tf_barplot.csv

Source table for `figures/by_contrast/Temp_main/tf_barplot.pdf` and
`figures/by_contrast/Temp_main/tf_barplot.png`; top TFs by BH padj plus watchlist
for the Temp_main contrast (genotype-averaged thermal TF response).

**How to read:** Same column structure as other `tf_barplot.csv` files. Positive `score` = TF activated by hyperthermia averaged over both genotypes. Provides context for interpreting the WT_heat vs. KO_heat pair: Temp_main captures the shared thermal TF response while the Interaction contrast isolates the cGAS-dependent component. Claim tier: L3 — provisional, n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf / save_overview` | `thresholds.gsea_fdr=0.05; figures.top_n=20` | `03_results/objects/03_tf_collectri.rds` |

