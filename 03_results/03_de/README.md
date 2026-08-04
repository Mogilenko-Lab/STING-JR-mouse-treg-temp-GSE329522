# Differential expression — the 39 °C heat response in WT and cGAS-KO iTregs

**Method:** limma-trend on log2(CPM+0.5) (no voom — input is pre-normalized CPM deposit).
Model: `~0 + group`, trend=TRUE, robust=TRUE. n=5/group. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

## Contrasts (7)

All 7 contrasts are built generically from `analysis_config.yaml:design.contrasts` via `makeContrasts`; no contrast is hard-coded in the script.

| Name | Formula | Biological meaning |
|---|---|---|
| WT_heat | WT\_39 − WT\_37 | Heat response in the cGAS-competent genotype |
| KO_heat | cGASKO\_39 − cGASKO\_37 | Heat response with cGAS removed |
| Interaction | (WT\_39 − cGASKO\_39) − (WT\_37 − cGASKO\_37) | **THE cGAS-dependence test** — 1 df interaction; underpowered at n=5 |
| Geno_at_39 | WT\_39 − cGASKO\_39 | Genotype effect under heat |
| Geno_at_37 | WT\_37 − cGASKO\_37 | Genotype effect at baseline |
| Temp_main | ½(WT\_39+cGASKO\_39) − ½(WT\_37+cGASKO\_37) | Pooled/marginal average 39-vs-37 °C response — NOT the cGAS-dependence test |
| Geno_main | ½(WT\_37+WT\_39) − ½(cGASKO\_37+cGASKO\_39) | Pooled/marginal average genotype effect |

## Artifacts produced by `02_de_limma_trend_viz.R`

- `figures/by_contrast/<contrast>/volcano.*` — per-contrast volcano plots for all 7 contrasts
- `figures/by_contrast/<contrast>/md.*` — per-contrast MD plots for all 7 contrasts
- `figures/_overview/de_counts_summary.*` — cross-contrast DE-count overview (figure-style)
- `tables/by_contrast/<contrast>/volcano.csv` — source table for each volcano panel (x7)
- `tables/by_contrast/<contrast>/md.csv` — source table for each MD panel (x7)
- `tables/_overview/de_counts_summary.csv` — per-contrast signed DE-count summary table

## Artifacts produced by `02b_cgas_dependence_geometry{,_viz}.R`

- `tables/_overview/cgas_dependence_wide.csv` — one row per gene, the four heat-relevant contrasts side by side (joined on `ensembl`)
- `tables/_overview/cgas_dependence_stats.csv` — the scalars the panel prints in-panel
- `figures/_overview/heat_response_wt_vs_ko.*` — WT heat response against cGAS-KO heat response, identity line drawn

## Headline panels (bespoke, owner-curated)

- `tables/marker_cgas_dependence.csv` — ISG / HIF marker statistics across WT_heat, KO_heat, Interaction
- `tables/fig2_marker_means.csv` — plot-ready per-group mean log2(CPM+0.5) for Fig 2 markers; **load-bearing cross-stage dep** (`03_decoupler_tf.R:687` reads this; path+schema frozen)

## House framing

- Never phrase as "cGAS-independent" — use "no detectable cGAS-dependence at n=5".
- Do not crown HIF1α/HIF2α as the effector; figure claims are floored at L3 (association/correlational); mechanism stays in prose.
- Interaction is 1 df and underpowered; treat FDR thresholds cautiously.
- n=5/group; sample→group labels are Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.


## tables/fig2_marker_means.csv

Per-group mean log2(CPM+0.5) tidy table for 15 cGAS-dependence marker
genes (7 ISG-arm: Ifit1, Isg15, Irf7, Oasl2, Mx1, Stat1, Cxcl10; 8
HIF/glycolysis-arm: Slc2a1, Vegfa, Egln3, Bnip3, Pgk1, Ldha, Aldoa,
Hk2) across 4 groups (WT_37, WT_39, cGASKO_37, cGASKO_39). Also
carries Interaction logFC and adj.P.Val from limma-trend (no
recomputation in viz). READ-ONLY in this script; schema-frozen (also
consumed by 03_decoupler_tf.R:687). Also the HIGHLIGHT watchlist
source for the volcano + MD sweeps. Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

**How to read:** Columns: gene, arm (IFN/ISG arm | HIF/glycolysis arm), hif_class
(HIF-specific | shared-glycolytic | NA for ISG), genotype (WT |
cGASKO), temp (37C | 39C), mean_log2cpm, inter_logFC, inter_adjP. One
row per gene x group (4 rows per gene). inter_adjP is constant within
gene. Schema frozen: do NOT rename/relocate; consumed by downstream
scripts. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `read.csv (read-only; produced by 02_de_limma_trend.R)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds + fig2_marker_means.csv` |

## tables/marker_cgas_dependence.csv

Per-marker Interaction contrast statistics table (subset of
02_de_results.rds, Interaction slot) for the cGAS-dependence marker
genes. Lists logFC, P.Value, adj.P.Val, AveExpr for each marker. Used
as supplemental evidence for the two-arms framing. Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

**How to read:** Rows = marker genes; columns = gene_symbol, ensembl, logFC, AveExpr,
t, P.Value, adj.P.Val, contrast. All values from the Interaction
contrast (1 df; logFC > 0 = up in WT relative to cGASKO after
removing shared heat effect). Non-significant adj.P = no detectable
cGAS-dependence at n=5 (NEVER 'cGAS-independent'). Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `produced by 02_de_limma_trend.R (read-only in viz)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/WT_heat/volcano.png

limma-trend volcano for WT_heat (WT: heat (39 vs 37 °C)): 339 genes
pass adj.P < 0.05 & |log2FC| >= 1.0 (213 up in numerator, 126 up in
denominator). cGAS-STING/ISG + HIF watchlist genes always labelled;
sig up/down counts appended (toolkit) as a second line under the
highest-priority populated significance line in the legend; x-tick
step chosen dynamically to avoid label crowding.

**How to read:** x = log2FC (>0 = up in numerator; <0 = up in denominator); y =
-log10(P.Value) (raw p on y; FDR for the color decision). Dashed
lines = adj.P < 0.05 boundary (horizontal) & |log2FC| >= 1.0
(vertical). Colored points = significant by the toolkit
decision-by-FDR rule. Labelled = top 10 genes by significance per
side, PLUS the always-on two-arms watchlist (ISG:
Ifit1/Isg15/Irf7/...; HIF/glyco: Slc2a1/Vegfa/Egln3/...). Data: CPM /
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/WT_heat/volcano.csv

Source table for the WT_heat volcano panel: all 19,679 genes tested by
limma-trend for the WT_heat contrast (WT\_39 − WT\_37), ordered by
adj.P.Val ascending then |logFC| descending. Carries the full per-gene
statistics used to render figures/by_contrast/WT_heat/volcano.png;
same-stem sibling to that figure. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up in WT\_39 vs WT\_37; <0 = down),
AveExpr (average log2CPM across all samples), t (moderated t-statistic), P.Value (raw),
adj.P.Val (BH-corrected FDR), contrast (constant = "WT_heat"). Significance rule for
figures: adj.P.Val < 0.05 AND |logFC| >= 1.0. n=5/group; 19,679 genes tested.
Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/WT_heat/md.png

limma-trend mean-difference (MD) plot for WT_heat (WT: heat (39 vs 37
°C)): logFC vs average expression. LOESS trend (red dashed) +
median-expression guide expose any expression-range normalization
bias; quadrant counts quantify directionality. cGAS-STING/ISG + HIF
watchlist genes always labelled.

**How to read:** x = average expression (AveExpr, log2CPM scale); y = log2FC. Grey
cloud = NS genes; orange = up (FDR < 0.05); blue = down (FDR < 0.05).
Red dashed = LOESS trend (should hug logFC~0; a tilt flags a
normalization confound); dotted vertical = median expression; dashed
horizontals = |log2FC| >= 1.0. Quadrant numbers = significant-gene
counts per quadrant. Labelled = top 5 sig genes per side, PLUS the
always-on two-arms watchlist. Driven from 02_de_results.rds topTable
columns (no MArrayLM fit on disk). Claim tier: L3; n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/WT_heat/md.csv

Source table for the WT_heat MD (mean-difference) plot: all 19,679 genes
tested by limma-trend for the WT_heat contrast (WT\_39 − WT\_37), ordered by
adj.P.Val ascending then |logFC| descending. Same gene-level statistics as
volcano.csv; same-stem sibling to figures/by_contrast/WT_heat/md.png.
Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis in the MD plot; >0 = up in
WT\_39 vs WT\_37), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val,
contrast. The MD panel plots AveExpr (x) vs logFC (y); this table is its
complete underlying data. Significance rule for coloring: adj.P.Val < 0.05.
n=5/group; 19,679 genes tested. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

---

## figures/by_contrast/KO_heat/volcano.png

limma-trend volcano for KO_heat (cGAS-KO: heat (39 vs 37 °C)): 392
genes pass adj.P < 0.05 & |log2FC| >= 1.0 (239 up in numerator, 153
up in denominator). cGAS-STING/ISG + HIF watchlist genes always
labelled; sig up/down counts appended (toolkit) as a second line
under the highest-priority populated significance line in the legend;
x-tick step chosen dynamically to avoid label crowding.

**How to read:** x = log2FC (>0 = up in numerator; <0 = up in denominator); y =
-log10(P.Value) (raw p on y; FDR for the color decision). Dashed
lines = adj.P < 0.05 boundary (horizontal) & |log2FC| >= 1.0
(vertical). Colored points = significant by the toolkit
decision-by-FDR rule. Labelled = top 10 genes by significance per
side, PLUS the always-on two-arms watchlist (ISG:
Ifit1/Isg15/Irf7/...; HIF/glyco: Slc2a1/Vegfa/Egln3/...). Data: CPM /
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/KO_heat/volcano.csv

Source table for the KO_heat volcano panel: all 19,679 genes tested by
limma-trend for the KO_heat contrast (cGASKO\_39 − cGASKO\_37), ordered by
adj.P.Val ascending then |logFC| descending. Carries the full per-gene
statistics used to render figures/by_contrast/KO_heat/volcano.png.
Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up in cGASKO\_39 vs cGASKO\_37),
AveExpr (average log2CPM), t, P.Value, adj.P.Val, contrast (constant = "KO_heat").
Significance rule for figures: adj.P.Val < 0.05 AND |logFC| >= 1.0. n=5/group;
19,679 genes tested. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/KO_heat/md.png

limma-trend mean-difference (MD) plot for KO_heat (cGAS-KO: heat (39
vs 37 °C)): logFC vs average expression. LOESS trend (red dashed) +
median-expression guide expose any expression-range normalization
bias; quadrant counts quantify directionality. cGAS-STING/ISG + HIF
watchlist genes always labelled.

**How to read:** x = average expression (AveExpr, log2CPM scale); y = log2FC. Grey
cloud = NS genes; orange = up (FDR < 0.05); blue = down (FDR < 0.05).
Red dashed = LOESS trend (should hug logFC~0; a tilt flags a
normalization confound); dotted vertical = median expression; dashed
horizontals = |log2FC| >= 1.0. Quadrant numbers = significant-gene
counts per quadrant. Labelled = top 5 sig genes per side, PLUS the
always-on two-arms watchlist. Driven from 02_de_results.rds topTable
columns (no MArrayLM fit on disk). Claim tier: L3; n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/KO_heat/md.csv

Source table for the KO_heat MD plot: all 19,679 genes tested by
limma-trend for the KO_heat contrast (cGASKO\_39 − cGASKO\_37), ordered by
adj.P.Val ascending then |logFC| descending. Same-stem sibling to
figures/by_contrast/KO_heat/md.png. Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's
sample sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis in MD plot; >0 = up in
cGASKO\_39 vs cGASKO\_37), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val,
contrast. Significance rule for coloring: adj.P.Val < 0.05. n=5/group; 19,679 genes
tested. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

---

## figures/by_contrast/Interaction/volcano.png

limma-trend volcano for Interaction (Heat × genotype interaction
(cGAS-dependence of the heat response)): 9 genes pass adj.P < 0.05 &
|log2FC| >= 1.0 (9 up in numerator, 0 up in denominator).
cGAS-STING/ISG + HIF watchlist genes always labelled; sig up/down
counts appended (toolkit) as a second line under the highest-priority
populated significance line in the legend; x-tick step chosen
dynamically to avoid label crowding. Interaction is THE
cGAS-dependence test (1 df, n=5, lowest-powered term); a
non-significant gene = no detectable cGAS-dependence at n=5, NOT
independence.

**How to read:** x = log2FC (>0 = up in numerator; <0 = up in denominator); y =
-log10(P.Value) (raw p on y; FDR for the color decision). Dashed
lines = adj.P < 0.05 boundary (horizontal) & |log2FC| >= 1.0
(vertical). Colored points = significant by the toolkit
decision-by-FDR rule. Labelled = top 10 genes by significance per
side, PLUS the always-on two-arms watchlist (ISG:
Ifit1/Isg15/Irf7/...; HIF/glyco: Slc2a1/Vegfa/Egln3/...). Data: CPM /
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Interaction/volcano.csv

Source table for the Interaction volcano panel: all 19,679 genes tested by
limma-trend for the Interaction contrast ((WT\_39 − cGASKO\_39) − (WT\_37 −
cGASKO\_37); the 1 df cGAS-dependence test), ordered by adj.P.Val ascending
then |logFC| descending. This is the definitive gene-level cGAS-dependence
readout from the model; non-significant adj.P.Val = no detectable
cGAS-dependence at n=5 (NEVER 'cGAS-independent'). Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = the heat response is larger in
WT than cGASKO, i.e., cGAS-dependent induction; <0 = heat response larger in cGASKO),
AveExpr, t, P.Value, adj.P.Val, contrast (constant = "Interaction"). 1 df test;
lowest-powered term in the model at n=5/group. 23 genes pass adj.P.Val < 0.05 (no
logFC gate); 9 pass the combined adj.P.Val < 0.05 AND |logFC| >= 1.0 threshold used
in the volcano figure. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Interaction/md.png

limma-trend mean-difference (MD) plot for Interaction (Heat ×
genotype interaction (cGAS-dependence of the heat response)): logFC
vs average expression. LOESS trend (red dashed) + median-expression
guide expose any expression-range normalization bias; quadrant counts
quantify directionality. cGAS-STING/ISG + HIF watchlist genes always
labelled. Interaction is THE cGAS-dependence test (1 df, n=5,
lowest-powered term); a non-significant gene = no detectable
cGAS-dependence at n=5, NOT independence.

**How to read:** x = average expression (AveExpr, log2CPM scale); y = log2FC. Grey
cloud = NS genes; orange = up (FDR < 0.05); blue = down (FDR < 0.05).
Red dashed = LOESS trend (should hug logFC~0; a tilt flags a
normalization confound); dotted vertical = median expression; dashed
horizontals = |log2FC| >= 1.0. Quadrant numbers = significant-gene
counts per quadrant. Labelled = top 5 sig genes per side, PLUS the
always-on two-arms watchlist. Driven from 02_de_results.rds topTable
columns (no MArrayLM fit on disk). Claim tier: L3; n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Interaction/md.csv

Source table for the Interaction MD plot: all 19,679 genes tested by
limma-trend for the Interaction contrast (1 df cGAS-dependence test),
ordered by adj.P.Val ascending then |logFC| descending. Same-stem
sibling to figures/by_contrast/Interaction/md.png. Non-significant
adj.P.Val = no detectable cGAS-dependence at n=5, NOT independence.
Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = cGAS-dependent heat induction;
<0 = cGAS-dependent heat suppression), AveExpr (average log2CPM), t, P.Value,
adj.P.Val, contrast. 1 df interaction; underpowered at n=5. Significance rule for
coloring in the MD figure: adj.P.Val < 0.05. Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's
sample sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

---

## figures/by_contrast/Geno_at_39/volcano.png

limma-trend volcano for Geno_at_39 (Genotype (WT vs cGAS-KO) @ 39
°C): 17 genes pass adj.P < 0.05 & |log2FC| >= 1.0 (17 up in
numerator, 0 up in denominator). cGAS-STING/ISG + HIF watchlist genes
always labelled; sig up/down counts appended (toolkit) as a second
line under the highest-priority populated significance line in the
legend; x-tick step chosen dynamically to avoid label crowding.

**How to read:** x = log2FC (>0 = up in numerator; <0 = up in denominator); y =
-log10(P.Value) (raw p on y; FDR for the color decision). Dashed
lines = adj.P < 0.05 boundary (horizontal) & |log2FC| >= 1.0
(vertical). Colored points = significant by the toolkit
decision-by-FDR rule. Labelled = top 10 genes by significance per
side, PLUS the always-on two-arms watchlist (ISG:
Ifit1/Isg15/Irf7/...; HIF/glyco: Slc2a1/Vegfa/Egln3/...). Data: CPM /
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_at_39/volcano.csv

Source table for the Geno_at_39 volcano panel: all 19,679 genes tested by
limma-trend for the Geno_at_39 contrast (WT\_39 − cGASKO\_39; genotype effect
under heat stress), ordered by adj.P.Val ascending then |logFC| descending.
Same-stem sibling to figures/by_contrast/Geno_at_39/volcano.png.
Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up in WT vs cGASKO at 39 °C;
<0 = up in cGASKO), AveExpr, t, P.Value, adj.P.Val, contrast (constant = "Geno_at_39").
Significance rule for figures: adj.P.Val < 0.05 AND |logFC| >= 1.0. n=5/group;
64 genes pass adj.P.Val < 0.05 (no logFC gate). Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Geno_at_39/md.png

limma-trend mean-difference (MD) plot for Geno_at_39 (Genotype (WT vs
cGAS-KO) @ 39 °C): logFC vs average expression. LOESS trend (red
dashed) + median-expression guide expose any expression-range
normalization bias; quadrant counts quantify directionality.
cGAS-STING/ISG + HIF watchlist genes always labelled.

**How to read:** x = average expression (AveExpr, log2CPM scale); y = log2FC. Grey
cloud = NS genes; orange = up (FDR < 0.05); blue = down (FDR < 0.05).
Red dashed = LOESS trend (should hug logFC~0; a tilt flags a
normalization confound); dotted vertical = median expression; dashed
horizontals = |log2FC| >= 1.0. Quadrant numbers = significant-gene
counts per quadrant. Labelled = top 5 sig genes per side, PLUS the
always-on two-arms watchlist. Driven from 02_de_results.rds topTable
columns (no MArrayLM fit on disk). Claim tier: L3; n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_at_39/md.csv

Source table for the Geno_at_39 MD plot: all 19,679 genes tested by
limma-trend for the Geno_at_39 contrast (WT\_39 − cGASKO\_39), ordered by
adj.P.Val ascending then |logFC| descending. Same-stem sibling to
figures/by_contrast/Geno_at_39/md.png. Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's
sample sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis; >0 = up in WT vs cGASKO
at 39 °C), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val, contrast.
Significance rule for MD coloring: adj.P.Val < 0.05. n=5/group; 19,679 genes tested.
Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

---

## figures/by_contrast/Geno_at_37/volcano.png

limma-trend volcano for Geno_at_37 (Genotype (WT vs cGAS-KO) @ 37
°C): 2 genes pass adj.P < 0.05 & |log2FC| >= 1.0 (2 up in numerator,
0 up in denominator). cGAS-STING/ISG + HIF watchlist genes always
labelled; sig up/down counts appended (toolkit) as a second line
under the highest-priority populated significance line in the legend;
x-tick step chosen dynamically to avoid label crowding.

**How to read:** x = log2FC (>0 = up in numerator; <0 = up in denominator); y =
-log10(P.Value) (raw p on y; FDR for the color decision). Dashed
lines = adj.P < 0.05 boundary (horizontal) & |log2FC| >= 1.0
(vertical). Colored points = significant by the toolkit
decision-by-FDR rule. Labelled = top 10 genes by significance per
side, PLUS the always-on two-arms watchlist (ISG:
Ifit1/Isg15/Irf7/...; HIF/glyco: Slc2a1/Vegfa/Egln3/...). Data: CPM /
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_at_37/volcano.csv

Source table for the Geno_at_37 volcano panel: all 19,679 genes tested by
limma-trend for the Geno_at_37 contrast (WT\_37 − cGASKO\_37; baseline genotype
effect at 37 °C), ordered by adj.P.Val ascending then |logFC| descending.
Same-stem sibling to figures/by_contrast/Geno_at_37/volcano.png.
Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up in WT vs cGASKO at
37 °C; <0 = up in cGASKO), AveExpr, t, P.Value, adj.P.Val, contrast (constant =
"Geno_at_37"). Significance rule for figures: adj.P.Val < 0.05 AND |logFC| >= 1.0.
n=5/group; 5 genes pass adj.P.Val < 0.05 (no logFC gate). Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's
sample sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Geno_at_37/md.png

limma-trend mean-difference (MD) plot for Geno_at_37 (Genotype (WT vs
cGAS-KO) @ 37 °C): logFC vs average expression. LOESS trend (red
dashed) + median-expression guide expose any expression-range
normalization bias; quadrant counts quantify directionality.
cGAS-STING/ISG + HIF watchlist genes always labelled.

**How to read:** x = average expression (AveExpr, log2CPM scale); y = log2FC. Grey
cloud = NS genes; orange = up (FDR < 0.05); blue = down (FDR < 0.05).
Red dashed = LOESS trend (should hug logFC~0; a tilt flags a
normalization confound); dotted vertical = median expression; dashed
horizontals = |log2FC| >= 1.0. Quadrant numbers = significant-gene
counts per quadrant. Labelled = top 5 sig genes per side, PLUS the
always-on two-arms watchlist. Driven from 02_de_results.rds topTable
columns (no MArrayLM fit on disk). Claim tier: L3; n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_at_37/md.csv

Source table for the Geno_at_37 MD plot: all 19,679 genes tested by
limma-trend for the Geno_at_37 contrast (WT\_37 − cGASKO\_37), ordered by
adj.P.Val ascending then |logFC| descending. Same-stem sibling to
figures/by_contrast/Geno_at_37/md.png. Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's
sample sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis; >0 = up in WT vs cGASKO
at 37 °C), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val, contrast.
Significance rule for MD coloring: adj.P.Val < 0.05. n=5/group; 19,679 genes tested.
Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

---

## figures/by_contrast/Temp_main/volcano.png

limma-trend volcano for Temp_main (Heat main effect (pooled
genotypes)): 346 genes pass adj.P < 0.05 & |log2FC| >= 1.0 (216 up in
numerator, 130 up in denominator). cGAS-STING/ISG + HIF watchlist
genes always labelled; sig up/down counts appended (toolkit) as a
second line under the highest-priority populated significance line in
the legend; x-tick step chosen dynamically to avoid label crowding.
Pooled/marginal heat effect (collapses genotype); NOT the
cGAS-dependence read-out (that is the Interaction).

**How to read:** x = log2FC (>0 = up in numerator; <0 = up in denominator); y =
-log10(P.Value) (raw p on y; FDR for the color decision). Dashed
lines = adj.P < 0.05 boundary (horizontal) & |log2FC| >= 1.0
(vertical). Colored points = significant by the toolkit
decision-by-FDR rule. Labelled = top 10 genes by significance per
side, PLUS the always-on two-arms watchlist (ISG:
Ifit1/Isg15/Irf7/...; HIF/glyco: Slc2a1/Vegfa/Egln3/...). Data: CPM /
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Temp_main/volcano.csv

Source table for the Temp_main volcano panel: all 19,679 genes tested by
limma-trend for the Temp_main contrast (½(WT\_39+cGASKO\_39) − ½(WT\_37+cGASKO\_37);
pooled/marginal heat main effect), ordered by adj.P.Val ascending then |logFC|
descending. Same-stem sibling to figures/by_contrast/Temp_main/volcano.png.
This contrast collapses genotype; it is NOT the cGAS-dependence readout (that is
the Interaction contrast). Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up at 39 °C vs 37 °C pooled
over genotypes; <0 = down under heat), AveExpr, t, P.Value, adj.P.Val, contrast
(constant = "Temp_main"). Significance rule for figures: adj.P.Val < 0.05 AND
|logFC| >= 1.0. n=5/group (10 samples total, pooled); 11,153 genes pass adj.P.Val < 0.05
(no logFC gate). Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Temp_main/md.png

limma-trend mean-difference (MD) plot for Temp_main (Heat main effect
(pooled genotypes)): logFC vs average expression. LOESS trend (red
dashed) + median-expression guide expose any expression-range
normalization bias; quadrant counts quantify directionality.
cGAS-STING/ISG + HIF watchlist genes always labelled. Pooled/marginal
heat effect (collapses genotype); NOT the cGAS-dependence read-out
(that is the Interaction).

**How to read:** x = average expression (AveExpr, log2CPM scale); y = log2FC. Grey
cloud = NS genes; orange = up (FDR < 0.05); blue = down (FDR < 0.05).
Red dashed = LOESS trend (should hug logFC~0; a tilt flags a
normalization confound); dotted vertical = median expression; dashed
horizontals = |log2FC| >= 1.0. Quadrant numbers = significant-gene
counts per quadrant. Labelled = top 5 sig genes per side, PLUS the
always-on two-arms watchlist. Driven from 02_de_results.rds topTable
columns (no MArrayLM fit on disk). Claim tier: L3; n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Temp_main/md.csv

Source table for the Temp_main MD plot: all 19,679 genes tested by
limma-trend for the Temp_main contrast (pooled heat main effect), ordered by
adj.P.Val ascending then |logFC| descending. Same-stem sibling to
figures/by_contrast/Temp_main/md.png. This contrast is NOT the
cGAS-dependence readout. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis; >0 = up under heat pooled
over genotypes), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val, contrast.
Significance rule for MD coloring: adj.P.Val < 0.05. n=5/group; 19,679 genes tested.
Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

---

## figures/by_contrast/Geno_main/volcano.png

limma-trend volcano for Geno_main (Genotype main effect (pooled
temps)): 4 genes pass adj.P < 0.05 & |log2FC| >= 1.0 (4 up in
numerator, 0 up in denominator). cGAS-STING/ISG + HIF watchlist genes
always labelled; sig up/down counts appended (toolkit) as a second
line under the highest-priority populated significance line in the
legend; x-tick step chosen dynamically to avoid label crowding.
Pooled/marginal genotype effect (collapses temperature); NOT the
cGAS-dependence read-out (that is the Interaction).

**How to read:** x = log2FC (>0 = up in numerator; <0 = up in denominator); y =
-log10(P.Value) (raw p on y; FDR for the color decision). Dashed
lines = adj.P < 0.05 boundary (horizontal) & |log2FC| >= 1.0
(vertical). Colored points = significant by the toolkit
decision-by-FDR rule. Labelled = top 10 genes by significance per
side, PLUS the always-on two-arms watchlist (ISG:
Ifit1/Isg15/Irf7/...; HIF/glyco: Slc2a1/Vegfa/Egln3/...). Data: CPM /
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_main/volcano.csv

Source table for the Geno_main volcano panel: all 19,679 genes tested by
limma-trend for the Geno_main contrast (½(WT\_37+WT\_39) − ½(cGASKO\_37+cGASKO\_39);
pooled/marginal genotype main effect), ordered by adj.P.Val ascending then |logFC|
descending. Same-stem sibling to figures/by_contrast/Geno_main/volcano.png.
This contrast collapses temperature; it is NOT the cGAS-dependence readout (that is
the Interaction contrast). Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up in WT vs cGASKO pooled
over temperatures; <0 = up in cGASKO), AveExpr, t, P.Value, adj.P.Val, contrast
(constant = "Geno_main"). Significance rule for figures: adj.P.Val < 0.05 AND
|logFC| >= 1.0. n=5/group (10 samples total, pooled); 60 genes pass adj.P.Val < 0.05
(no logFC gate). Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Geno_main/md.png

limma-trend mean-difference (MD) plot for Geno_main (Genotype main
effect (pooled temps)): logFC vs average expression. LOESS trend (red
dashed) + median-expression guide expose any expression-range
normalization bias; quadrant counts quantify directionality.
cGAS-STING/ISG + HIF watchlist genes always labelled. Pooled/marginal
genotype effect (collapses temperature); NOT the cGAS-dependence
read-out (that is the Interaction).

**How to read:** x = average expression (AveExpr, log2CPM scale); y = log2FC. Grey
cloud = NS genes; orange = up (FDR < 0.05); blue = down (FDR < 0.05).
Red dashed = LOESS trend (should hug logFC~0; a tilt flags a
normalization confound); dotted vertical = median expression; dashed
horizontals = |log2FC| >= 1.0. Quadrant numbers = significant-gene
counts per quadrant. Labelled = top 5 sig genes per side, PLUS the
always-on two-arms watchlist. Driven from 02_de_results.rds topTable
columns (no MArrayLM fit on disk). Claim tier: L3; n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_main/md.csv

Source table for the Geno_main MD plot: all 19,679 genes tested by
limma-trend for the Geno_main contrast (pooled genotype main effect), ordered by
adj.P.Val ascending then |logFC| descending. Same-stem sibling to
figures/by_contrast/Geno_main/md.png. This contrast is NOT the
cGAS-dependence readout. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis; >0 = up in WT vs cGASKO
pooled over temperatures), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val,
contrast. Significance rule for MD coloring: adj.P.Val < 0.05. n=5/group; 19,679 genes
tested. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

---

## figures/_overview/de_counts_summary.png

Heat dominates this design and the cGAS-dependent arm is small and
one-sided: WT heat changes 8,723 genes and cGAS-KO heat 8,901, while
the interaction that tests cGAS-dependence of the heat response
changes 23 -- all 23 up, 0 down. Every bar prints its own gene count,
so that contrast reads without opening the table. Counts pass adj.P <
0.05 with no fold-change cut-off; the frozen mouse-to-human
projection signature adds |log2FC| >= 1.0 and lists 213 up / 126 down
for the same WT heat contrast
(10_signature/tables/_overview/signature_sizes.csv, gate fdr_logfc)
-- a stricter gate on the same statistics, not a different result.
The interaction is 1 df at n=5, the least-powered term: 23 is a
floor. Claim tier: L3.

**How to read:** x = comparison (config order, short display labels); y = number of
genes. Orange above zero = up, blue below zero = down, both relative
to the first-named side of the comparison; for the interaction, up =
heat response larger in WT. The printed label on each bar end gives
the exact count per direction, including a measured 0. Gate for this
panel: adj.P < 0.05, no fold-change cut-off -- the |log2FC| >= 1.0
gate used by the projection export is stricter and keeps far fewer
genes. A non-significant interaction gene = no detectable
cGAS-dependence at n=5, never independence. Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `geom_col + geom_text (signed bars with printed per-direction counts)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/master/master_de_genes.csv` |

## tables/_overview/de_counts_summary.csv

Per-comparison tallies of genes tested, genes significant at adj.P <
0.05 with no fold-change cut-off, and the split of those into up and
down. Source table for figures/_overview/de_counts_summary.png, whose
printed bar labels are exactly n_up and n_down. Claim tier: L3.
Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

**How to read:** Columns: contrast (config key), n_tested (19,679 genes in every
contrast's topTable), n_sig (adj.P.Val < 0.05, any fold change), n_up and n_down
(that split, relative to the first-named side of the comparison; n_up + n_down =
n_sig). The volcano and MD panels apply the stricter adj.P.Val < 0.05 AND |logFC|
>= 1.0 gate, so their counts run lower. Interaction: n_sig = 23, all up. A
non-significant interaction gene = no detectable cGAS-dependence at n=5, never
independence. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `geom_col + geom_text (signed bars with printed per-direction counts) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/master/master_de_genes.csv` |


## figures/_overview/heat_response_wt_vs_ko.png

Warming iTregs to 39 °C changes 8,723 genes in WT and 8,901 in
cGAS-KO, and it changes them in step: among the 10,418
heat-responsive genes the cGAS-KO response is 0.99× the WT response
(r = 0.95, Spearman 0.93), so the dots pile onto the identity line. A
one-sided minority behaves differently. All 23 genes that pass the
genotype comparison of the heat response lie below the line, a weaker
heat response without cGAS, and 0 lie above. 14 keep that direction
in both genotypes and 9 flip sign: up with heat in WT, down without
cGAS. Claim tier: L3.

**How to read:** One dot per gene: x = log2 fold change at 39 vs 37 °C in WT, y = the
same in cGAS-KO, equal scales, so the dashed identity line runs at
45°. Glyphs run pale to dark. Pale dots have no detectable difference
between the genotypes at n=5. Vermillion circles fall with heat in
both genotypes and fall further once cGAS is gone. Black triangles
rise with heat in WT and fall without cGAS. Highlighted genes pass at
adj.P < 0.05, and all 23 of them sit below the dashed line, which is
a weaker heat response without cGAS. The dotted lines lie 1 log2 unit
either side of the identity. The 9 highlighted genes beyond the lower
dotted line, their names in bold, also clear |logFC| >= 1, and 4 of
those are triangles. Bold marks that stringent core. The frozen
03_results/human_projection/ contract exports the arm at both gates,
these 9 at fdr_logfc and all 23 at fdr_only, which is the sensitivity
a 1 df term at n=5 needs (per-gate ortholog counts:
human_projection/manifest.csv). Read one-sidedness as a statement
about the highlighted arm alone. The pale cloud straddles the dashed
line in both directions. Claim tier: L3, n=5/group.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02b_cgas_dependence_geometry_viz.R` | `geom_point + geom_abline (effect-versus-effect scatter, equal axes, gate as parallels)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1; figures.volcano_label_top=10; figures.point_size=2.4; figures.cue_size=4.0; colors.okabe_ito` | `03_results/03_de/tables/_overview/cgas_dependence_wide.csv + cgas_dependence_stats.csv` |

## tables/_overview/cgas_dependence_wide.csv

One row per gene with the four heat-relevant contrasts side by side,
joined on ensembl: the WT and cGAS-KO heat responses, their shared
average, and the genotype comparison of the heat response. 19,679
genes, all tested in every contrast. This is the single source for
the scatter panel and the table in which the panel's geometric claim
is checkable: wt_minus_ko equals logFC_Interaction gene by gene.
Curated composition: 12 of the 23 are Hallmark interferon-alpha or
-gamma members and 11 are unassigned, several of them mouse-specific
paralogs the human-derived sets omit. Claim tier: L3.

**How to read:** Columns: ensembl, gene_symbol, then logFC_/adjP_/sig_ per contrast
(WT_heat, KO_heat, Interaction, Temp_main). sig_ = adj.P.Val < 0.05.
heat_responsive = significant for heat in at least one genotype.
cgas_dependent = passes the genotype comparison of the heat response.
up_with_heat_in_wt and down_with_heat_in_ko are logFC signs,
reverses_without_cgas is both at once, in_stringent_gate adds |logFC|
>= 1, and gate_class and arm_class fold these into the classes the
panel draws. Temp_main is a scored contrast the panel does not plot.
wt_minus_ko = the WT effect minus the cGAS-KO effect.
interaction_rank and label_rank order the arm and are NA off it.
Positive logFC = higher at 39 °C, except the genotype comparison,
where positive = a larger heat response in WT. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02b_cgas_dependence_geometry.R` | `inner_join on ensembl (read-only here; produced by the compute sibling)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1; figures.volcano_label_top=10; figures.point_size=2.4; figures.cue_size=4.0; colors.okabe_ito` | `03_results/objects/02_de_results.rds` |

## tables/_overview/cgas_dependence_stats.csv

The scalars the scatter panel prints, so no number on a figure is
computed at draw time. Correlations and the regression slope between
the two per-genotype heat responses are reported twice, over all
genes and over the heat-responsive ones only, because agreement
measured on the whole universe could be carried by unchanged genes
sitting at the origin. The restricted scope is the one the panel
quotes, and it is the higher of the two. Also carries per-contrast
counts with their up/down split, the anatomy of the cGAS-dependent
arm, typical effect sizes, axis ranges, and the identity residual.
Claim tier: L3.

**How to read:** Long format: metric, scope, subset, value, note, so a lookup is
(metric, scope, subset). subset is all_genes, heat_responsive or
cgas_dependent_arm and must always be read -- pearson_r appears under
two of them. The arm rows count how many of the significant genes
fall with heat in cGAS-KO, how many reverse sign, how many clear the
stringent |logFC| gate, and how far the two nine-gene subsets
overlap. axis_lim gives a symmetric half-range per contrast and
de_logfc the offset the panel draws as dotted lines; y_expansion and
the Temp_main/Interaction axis_lim rows are unused.
max_abs_identity_residual is 0 to numerical precision, which is what
licenses reading distance from the identity line as the genotype
difference. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02b_cgas_dependence_geometry.R` | `cor + lm + tallies (read-only here; produced by the compute sibling)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1; figures.volcano_label_top=10; figures.point_size=2.4; figures.cue_size=4.0; colors.okabe_ito` | `03_results/objects/02_de_results.rds` |

## tables/_overview/heat_response_wt_vs_ko.csv

Plotted values behind figures/_overview/heat_response_wt_vs_ko.png:
one row per gene, the columns the panel draws, sliced from
cgas_dependence_wide.csv. 19,679 genes. Claim tier: L3.

**How to read:** logFC_ / adjP_ columns are the log2 fold change and BH-adjusted
p-value per contrast. cgas_dependent = adj.P < 0.05 on the genotype
comparison of the heat response. in_stringent_gate adds |logFC| >= 1
and is the panel's dotted-line set; reverses_without_cgas is the sign
flip drawn as triangles; arm_class folds the two into the three
plotted classes. interaction_rank orders the arm by evidence and is
NA elsewhere. hallmark_ifn flags membership of a curated interferon
set and is annotation, never a selection rule. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02b_cgas_dependence_geometry_viz.R` | `save_overview (writes the figure's same-stem source table)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1; figures.volcano_label_top=10; figures.point_size=2.4; figures.cue_size=4.0; colors.okabe_ito` | `03_results/03_de/tables/_overview/cgas_dependence_wide.csv` |

