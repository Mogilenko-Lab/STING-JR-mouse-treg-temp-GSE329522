# 05_progeny — artifact captions

PROGENy MLM pathway activity inference (14 pathways, all contrasts)
+ CollecTRI ULM TF activity (decoupleR). Key result: Hypoxia (PROGENy) is flat in
the Interaction contrast — no detectable cGAS-dependence at n=5 (PROVISIONAL);
JAK-STAT/NFkB/TNFa are positive in the Interaction (cGAS-dependent arm).
PROVISIONAL — inferred sample mapping pending collaborator sample sheet.

Figures use the unified style contract: one `.pdf` (print) + one `.png` (screen) per
stem, no `.screen`/`.print` suffix. Heatmap x-axis contrast labels are drawn at
30-degree angle (`hjust = 1`) to prevent overlap in dense multi-column layouts.

---

## figures/_overview/progeny_heatmap.pdf

PROGENy 14-pathway x 7-contrast activity heatmap (MLM score, clamped +/-2.5;
* raw p<0.05; rows hierarchically clustered): Hypoxia/glycolysis rise in BOTH heat
arms but are flat in the Interaction; JAK-STAT/NFkB/TNFa rise in WT_heat and are
positive in the Interaction (cGAS-dependent). PROVISIONAL; n=5/group.

**How to read:** Rows = PROGENy pathways (hierarchically clustered via `hclust(dist(mat))`); columns = contrasts in design order (left = WT_heat, KO_heat; then Interaction, Temp_main, Geno*, Geno_main). Column labels drawn at 30-degree angle (hjust = 1). Fill: orange = pathway activated in numerator; blue = activated in denominator (diverging scale via `scale_fill_gradient2`, clamped ±2.5 with `scales::squish`). Glyph: `*` inside tile = raw p < 0.05 (via `.sig_star()`); no multi-test correction applied to n=14 pathways. Sign convention: positive score = pathway gene program more active in the numerator condition of each contrast. Claim tier: L3 — statistical inference from decoupleR MLM regulon model; provisional (n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.heat_tile (ggplot geom_tile)` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/_overview/progeny_heatmap.png

PROGENy 14-pathway x 7-contrast activity heatmap (MLM score, clamped
+/-2.5; * raw p<0.05; rows clustered): Hypoxia/glycolysis rise in
BOTH heat arms but are flat in the Interaction; JAK-STAT/NFkB/TNFa
rise in WT_heat and are positive in the Interaction (cGAS-dependent).
PROVISIONAL; n=5/group.

**How to read:** Rows = PROGENy pathways (hierarchically clustered); columns =
contrasts in design order (left = WT_heat, KO_heat; then Interaction,
Temp_main, Geno*, Geno_main). Fill: orange = pathway activated in
numerator; blue = activated in denominator. Score clamped to +/-2.5;
* = raw p < 0.05 (n=14 pathways; no multi-test correction warranted).
Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.heat_tile (ggplot geom_tile)` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/_overview/progeny_interaction_split.pdf

PROGENy Hypoxia vs immune split: Hypoxia activity is similar in WT_heat and
KO_heat with a flat Interaction (no detectable cGAS-dependence at n=5), whereas
JAK-STAT/NFkB/TNFa are positive in the Interaction (cGAS-dependent arm). This
orthogonal pathway footprint corroborates the two-arms DE/TF result. PROVISIONAL;
n=5/group.

**How to read:** Grouped bar chart; x-axis = 3 headline contrasts with two-line labels ("WT heat / 39 vs 37 C", "KO heat / 39 vs 37 C", "Interaction / cGAS-dependence"); fill = PROGENy pathway (Hypoxia = orange via `POS` config color; JAK-STAT = steelblue4; NFkB = cornflowerblue; TNFa = lightsteelblue); bar height = MLM activity score (positive = numerator activation). Figure subtitle: "Hypoxia (orange): flat Interaction (no detectable cGAS-dependence). JAK-STAT/NFkB/TNFa (blue): positive Interaction (cGAS-dependent). PROVISIONAL." Glyph: `*` nudged above bar = raw p < 0.05. Sign convention: positive score = pathway more active in numerator (e.g., for Interaction, numerator = WT_heat, so positive = cGAS-dependent). CAUTION: a flat Interaction for Hypoxia is NOT proven cGAS-independence; the study is underpowered for that comparison at n=5/group. Claim tier: L3 (provisional).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `ggplot/geom_col` | `thresholds.gsea_fdr=0.05; colors.diverging; design.contrasts` | `03_results/objects/09_progeny_activity.rds` |

## figures/_overview/progeny_interaction_split.png

PROGENy Hypoxia vs immune split: Hypoxia activity is similar in
WT_heat and KO_heat with a flat Interaction (no detectable
cGAS-dependence at n=5), whereas JAK-STAT/NFkB/TNFa are positive in
the Interaction (cGAS-dependent arm). This orthogonal pathway
footprint corroborates the two-arms DE/TF result. PROVISIONAL;
n=5/group.

**How to read:** Grouped bars: x = contrast (3 headline contrasts: WT_heat, KO_heat,
Interaction), fill = PROGENy pathway. Orange (Hypoxia) bars should be
similar in WT_heat and KO_heat and near-zero in the Interaction panel
- meaning no detectable cGAS-dependence. Blue bars
(JAK-STAT/NFkB/TNFa) should be taller in WT_heat than KO_heat and
positive in Interaction - cGAS-dependent. * = raw p < 0.05. CAUTION:
'flat Interaction' is NOT proven cGAS-independence; the study is
powered at n=5/group. Claim tier: L3 (provisional).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `ggplot/geom_col` | `thresholds.gsea_fdr=0.05; colors.diverging; design.contrasts` | `03_results/objects/09_progeny_activity.rds` |

## figures/_overview/tf_heatmap.pdf

CollecTRI TF activity heatmap (TFs significant in >=2 contrasts + watchlist;
rows = TFs, cols = contrasts; score clamped +/-2.5; * BH padj<0.05; left-edge
strip = HIF/IFN/other axis): IFN/IRF/STAT TFs cluster as the cGAS-dependent block
(positive in Interaction); HIF-axis TFs are non-significant in the Interaction.
PROVISIONAL; n=5/group.

**How to read:** Rows = CollecTRI TFs (hierarchically clustered); columns = contrasts in design order. Column labels drawn at 30-degree angle (hjust = 1). Fill: orange = TF activated in numerator; blue = activated in denominator (diverging scale clamped ±2.5). Left-edge annotation strip (via `ggnewscale::new_scale_fill`): orange tile = HIF axis (Hif1a/Epas1), blue tile = IFN/NFkB axis (Irf1/Irf3/Irf7/Stat1/Stat2/Nfkb1/Rela), grey = other. Glyph: `*` inside score tile = BH padj < 0.05 (within-contrast BH correction). Row selection: TFs with BH padj < 0.05 in ≥2 contrasts, up to top 40 by minimum padj, union with KEY_TFS watchlist. Sign convention: positive score = TF more active in numerator condition. Claim tier: L3 (provisional, n=5/group; IFN-arm TFs positive in Interaction = cGAS-dependent; HIF-axis TFs flat/NS in Interaction = no detectable cGAS-dependence at n=5; NOT proven independence).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.heat_tile (ggplot geom_tile)` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/_overview/tf_heatmap.png

CollecTRI TF activity heatmap (TFs significant in >=2 contrasts +
watchlist; rows = TFs, cols = contrasts; score clamped +/-2.5; * BH
padj<0.05; row strip = HIF/IFN/other axis): IFN/IRF/STAT TFs cluster
as the cGAS-dependent block (positive in Interaction); HIF-axis TFs
are non-significant in the Interaction. PROVISIONAL; n=5/group.

**How to read:** Rows = CollecTRI TFs (hierarchically clustered); columns = contrasts.
Fill: orange = TF activated in numerator; blue = activated in
denominator. Left-edge strip: orange = HIF axis (Hif1a/Epas1), blue =
IFN/NFkB axis, grey = other. Score clamped to +/-2.5. * = BH padj <
0.05. Claim tier: L3 (provisional, n=5/group; IFN-arm TFs positive in
Interaction = cGAS-dependent; HIF-axis TFs flat/NS in Interaction =
no detectable cGAS-dependence at n=5; NOT proven independence).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.heat_tile (ggplot geom_tile)` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/_overview/progeny_tf_combined.pdf

Combined PROGENy (top panel) + CollecTRI TF (bottom panel) activity for key
pathways/TFs across 3 headline contrasts: Hypoxia (PROGENy) and Hif1a (TF) are
active in both heat arms and flat in the Interaction (no detectable cGAS-dependence
at n=5); JAK-STAT/IFN/Stat1/Irf3 are positive in the Interaction (cGAS-dependent).
PROVISIONAL.

**How to read:** `facet_grid(entity_type ~ contrast)`: row facets = PROGENy (top) vs CollecTRI TF (bottom); column facets = 3 headline contrasts. Figure subtitle: "PROGENy MLM (top) + CollecTRI ULM (bottom). Orange = HIF arm; blue = IFN/NFkB arm. PROVISIONAL." Lollipop color via `scale_color_manual(values = AXIS_COLORS)`: orange = HIF axis (Hypoxia pathway; Hif1a/Epas1 TFs), steelblue4 = IFN/NFkB axis (JAK-STAT/NFkB/TNFa pathways; Irf1/Irf3/Irf7/Stat1/Stat2/Nfkb1/Rela TFs), grey = other. Shape: diamond (`pch=18`) = PROGENy pathway, filled circle (`pch=16`) = TF. Direction: rightward (positive x) = activated in numerator condition. Glyph: `*` = p < 0.05 (raw p for PROGENy pathways; BH padj for TFs). Source rows ordered by mean score across contrasts within each entity-type facet. CAUTION: n=5/group; absence of significant Interaction for HIF arm is NOT proven cGAS-independence. Claim tier: L3 (provisional).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `ggplot/facet_grid` | `thresholds.gsea_fdr=0.05; colors.diverging; design.contrasts` | `03_results/objects/09_progeny_activity.rds; 03_results/objects/03_tf_collectri.rds` |

## figures/_overview/progeny_tf_combined.png

Combined PROGENy (top panel) + CollecTRI TF (bottom panel) activity
for key pathways/TFs across 3 headline contrasts: Hypoxia (PROGENy)
and Hif1a (TF) are active in both heat arms and flat in the
Interaction (no detectable cGAS-dependence at n=5);
JAK-STAT/IFN/Stat1/Irf3 are positive in the Interaction
(cGAS-dependent). PROVISIONAL.

**How to read:** Row facets = PROGENy (top) vs CollecTRI TF (bottom); column facets =
headline contrasts. Lollipop color: orange = HIF axis, blue =
IFN/NFkB axis. Shape: diamond = PROGENy pathway, circle = TF.
Direction: rightward = activated in the numerator condition. * = p <
0.05 (raw for PROGENy; BH padj for TF). Key comparison: Hypoxia/Hif1a
bars are near-zero in the Interaction column - no detectable
cGAS-dependence. JAK-STAT/Stat1/Irf3 bars are positive in the
Interaction column - cGAS-dependent. CAUTION: n=5/group; absence of
significant Interaction for HIF arm is NOT proven independence. Claim
tier: L3 (provisional).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `ggplot/facet_grid` | `thresholds.gsea_fdr=0.05; colors.diverging; design.contrasts` | `03_results/objects/09_progeny_activity.rds; 03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/WT_heat/progeny_barplot.pdf

PROGENy 14-pathway activity (MLM) in WT_heat: bars sorted by score; * raw p<0.05;
key pathways (Hypoxia/JAK-STAT/NFkB/TNFa) marked with bold outline. Hypoxia and
JAK-STAT expected to be positive (heat-induced glycolytic/HIF overlap + STING-IFN
positive control). PROVISIONAL; n=5/group.

**How to read:** Horizontal bar chart (coord_flip); bars sorted ascending by score (left = most suppressed, right = most activated in WT_heat vs. WT_ctrl). Fill: diverging color (orange = positive score; steelblue4 = negative score) via `scale_fill_gradient2`, limits ±CAP (2.5) with squish. Bold outline = key pathway (Hypoxia/JAK-STAT/NFkB/TNFa); overlay via second `geom_col(fill=NA, color="black")`. Glyph: `*` nudged above bar tip = raw p < 0.05 (no BH correction; n=14 pathways). Axis: y = "Activity score (MLM)"; NOT a fold-change; NOT a GSEA-NES. Positive score = pathway more active in WT_heat numerator. Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/WT_heat/progeny_barplot.png

PROGENy 14-pathway activity (MLM) in WT_heat: bars sorted by score; *
raw p<0.05; key pathways (Hypoxia/JAK-STAT/NFkB/TNFa) bold-outlined.

**How to read:** Bars = MLM activity score (orange = pathway more active in numerator;
blue = more active in denominator). Glyphs: * = raw p < 0.05; bold
outline = key pathway. Score is not a fold-change; sign tracks
numerator activation direction. Claim tier: L3 (provisional,
n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/WT_heat/tf_barplot.pdf

Top CollecTRI TF activity (ULM) in WT_heat: lollipops colored by HIF/IFN axis;
* BH padj<0.05; IFN/HIF watchlist TFs marked with open circle outline.
PROVISIONAL; n=5/group.

**How to read:** Horizontal lollipop chart (coord_flip); top-20 TFs by BH padj plus forced KEY_TFS watchlist inclusions (Irf3/Irf7/Irf1/Stat1/Stat2/Nfkb1/Rela/Hif1a/Epas1), sorted ascending by score. Lollipop stem: `geom_segment` colored by axis family. Point: `geom_point(size=3)` colored by axis family — orange = HIF axis (Hif1a/Epas1), steelblue4 = IFN/NFkB axis (Irf*/Stat*/Nfkb1/Rela), grey55 = other. Open circle outline (`geom_point(shape=21, fill=NA)`) = watchlist TF. Glyph: `*` nudged past point = BH padj < 0.05 (within-contrast BH correction). Axis: y = "Activity score (ULM)"; positive = TF more active in WT_heat numerator. Score is NOT a fold-change or NES. Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/WT_heat/tf_barplot.png

Top CollecTRI TF activity (ULM) in WT_heat: lollipops colored by
HIF/IFN axis; * BH padj<0.05; IFN/HIF watchlist bold-outlined.

**How to read:** Lollipops = ULM activity score (rightward = pathway-level activation
in numerator). Color: orange = HIF axis (Hif1a/Epas1), blue =
IFN/NFkB axis, grey = other. Open circle = watchlist TF. * = BH padj
< 0.05. Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/KO_heat/progeny_barplot.pdf

PROGENy 14-pathway activity (MLM) in KO_heat: bars sorted by score; * raw p<0.05;
key pathways (Hypoxia/JAK-STAT/NFkB/TNFa) marked with bold outline. Compare with
WT_heat to identify cGAS-dependent attenuation of pathway responses. PROVISIONAL;
n=5/group.

**How to read:** Same layout as `WT_heat/progeny_barplot`. Fill: diverging (orange = positive/numerator; steelblue4 = negative). Bold outline = key pathways. `*` = raw p < 0.05. Positive score = pathway more active in KO_heat numerator. Hypoxia expected to remain positive (cGAS-independent thermal program); JAK-STAT expected to be reduced vs. WT_heat (cGAS-dependent attenuation). Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/KO_heat/progeny_barplot.png

PROGENy 14-pathway activity (MLM) in KO_heat: bars sorted by score; *
raw p<0.05; key pathways (Hypoxia/JAK-STAT/NFkB/TNFa) bold-outlined.

**How to read:** Bars = MLM activity score (orange = pathway more active in numerator;
blue = more active in denominator). Glyphs: * = raw p < 0.05; bold
outline = key pathway. Score is not a fold-change; sign tracks
numerator activation direction. Claim tier: L3 (provisional,
n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/KO_heat/tf_barplot.pdf

Top CollecTRI TF activity (ULM) in KO_heat: lollipops colored by HIF/IFN axis;
* BH padj<0.05; IFN/HIF watchlist bold-outlined. Compare IFN-axis TF scores
(Stat1/Irf3/Irf7) with WT_heat for cGAS-dependent attenuation. PROVISIONAL;
n=5/group.

**How to read:** Same layout as `WT_heat/tf_barplot`. Orange = HIF axis, steelblue4 = IFN/NFkB axis, grey = other. Open circle = watchlist TF. `*` = BH padj < 0.05. Positive score = TF more active in KO_heat numerator. Attenuation of IFN-axis TFs relative to WT_heat panel is the cGAS-dependent signal. Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/KO_heat/tf_barplot.png

Top CollecTRI TF activity (ULM) in KO_heat: lollipops colored by
HIF/IFN axis; * BH padj<0.05; IFN/HIF watchlist bold-outlined.

**How to read:** Lollipops = ULM activity score (rightward = pathway-level activation
in numerator). Color: orange = HIF axis (Hif1a/Epas1), blue =
IFN/NFkB axis, grey = other. Open circle = watchlist TF. * = BH padj
< 0.05. Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/Interaction/progeny_barplot.pdf

PROGENy 14-pathway activity (MLM) in the Interaction contrast (WT_heat minus
KO_heat differential): positive score = cGAS-dependent pathway activation; key
result is Hypoxia near-zero (no detectable cGAS-dependence at n=5) vs.
JAK-STAT/NFkB/TNFa positive (cGAS-dependent). PROVISIONAL; n=5/group.

**How to read:** Same layout as `WT_heat/progeny_barplot`. Fill: orange = positive score = pathway activation requires intact cGAS; steelblue4 = negative score = more active in KO. Bold outline = key pathways (Hypoxia/JAK-STAT/NFkB/TNFa). `*` = raw p < 0.05. Key comparison: Hypoxia bar near zero (no detectable cGAS-dependence) vs. JAK-STAT/NFkB/TNFa bars positive. CAUTION: near-zero Interaction for Hypoxia is NOT proven cGAS-independence at n=5 (underpowered). Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/Interaction/progeny_barplot.png

PROGENy 14-pathway activity (MLM) in Interaction: bars sorted by
score; * raw p<0.05; key pathways (Hypoxia/JAK-STAT/NFkB/TNFa)
bold-outlined.

**How to read:** Bars = MLM activity score (orange = pathway more active in numerator;
blue = more active in denominator). Glyphs: * = raw p < 0.05; bold
outline = key pathway. Score is not a fold-change; sign tracks
numerator activation direction. Claim tier: L3 (provisional,
n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/Interaction/tf_barplot.pdf

Top CollecTRI TF activity (ULM) in the Interaction contrast (positive score = TF
activity requires intact cGAS): Stat1/Irf3/Irf7 expected positive and significant;
Hif1a/Epas1 expected near-zero. PROVISIONAL; n=5/group.

**How to read:** Same layout as `WT_heat/tf_barplot`. Orange = HIF axis, steelblue4 = IFN/NFkB axis, grey = other. Open circle = watchlist TF. `*` = BH padj < 0.05. Positive score = TF activity is cGAS-dependent (more active in WT_heat than KO_heat at 39 C). Key result: Stat1/Irf3/Irf7 positive and significant = cGAS-dependent IFN program; Hif1a/Epas1 near-zero = no detectable cGAS-dependence for HIF axis at n=5 (NOT proven independence). Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/Interaction/tf_barplot.png

Top CollecTRI TF activity (ULM) in Interaction: lollipops colored by
HIF/IFN axis; * BH padj<0.05; IFN/HIF watchlist bold-outlined.

**How to read:** Lollipops = ULM activity score (rightward = pathway-level activation
in numerator). Color: orange = HIF axis (Hif1a/Epas1), blue =
IFN/NFkB axis, grey = other. Open circle = watchlist TF. * = BH padj
< 0.05. Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/Temp_main/progeny_barplot.pdf

PROGENy 14-pathway activity (MLM) in Temp_main (average hyperthermia effect across
WT and KO genotypes, 39 C vs. 37 C): captures the shared thermal pathway response
irrespective of cGAS status. PROVISIONAL; n=5/group.

**How to read:** Same layout as `WT_heat/progeny_barplot`. Positive score = pathway more active at 39 C vs. 37 C averaged across both genotypes. Provides the genotype-agnostic thermal baseline; compare with WT_heat and KO_heat panels to decompose cGAS-dependent vs. cGAS-independent thermal responses. `*` = raw p < 0.05. Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/Temp_main/progeny_barplot.png

PROGENy 14-pathway activity (MLM) in Temp_main: bars sorted by score;
* raw p<0.05; key pathways (Hypoxia/JAK-STAT/NFkB/TNFa)
bold-outlined.

**How to read:** Bars = MLM activity score (orange = pathway more active in numerator;
blue = more active in denominator). Glyphs: * = raw p < 0.05; bold
outline = key pathway. Score is not a fold-change; sign tracks
numerator activation direction. Claim tier: L3 (provisional,
n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_progeny` | `thresholds.gsea_fdr=0.05; figures.z_clamp=2.5; colors.diverging` | `03_results/objects/09_progeny_activity.rds` |

## figures/by_contrast/Temp_main/tf_barplot.pdf

Top CollecTRI TF activity (ULM) in Temp_main (genotype-averaged thermal TF
response, 39 C vs. 37 C): shows TFs activated by hyperthermia independent of cGAS
status. PROVISIONAL; n=5/group.

**How to read:** Same layout as `WT_heat/tf_barplot`. Orange = HIF axis, steelblue4 = IFN/NFkB axis, grey = other. Open circle = watchlist TF. `*` = BH padj < 0.05. Positive score = TF more active at 39 C vs. 37 C averaged over both genotypes. Context for WT_heat vs. KO_heat comparison: Temp_main captures the shared thermal TF response; the Interaction contrast isolates the cGAS-dependent component. Claim tier: L3 (provisional, n=5/group).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/13_activity_viz.R` | `.barplot_tf` | `thresholds.gsea_fdr=0.05; figures.top_n=20; colors.diverging` | `03_results/objects/03_tf_collectri.rds` |

## figures/by_contrast/Temp_main/tf_barplot.png

Top CollecTRI TF activity (ULM) in Temp_main: lollipops colored by
HIF/IFN axis; * BH padj<0.05; IFN/HIF watchlist bold-outlined.

**How to read:** Lollipops = ULM activity score (rightward = pathway-level activation
in numerator). Color: orange = HIF axis (Hif1a/Epas1), blue =
IFN/NFkB axis, grey = other. Open circle = watchlist TF. * = BH padj
< 0.05. Claim tier: L3 (provisional, n=5/group).

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

