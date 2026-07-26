# Differential expression — the 39 °C heat response in WT and cGAS-KO iTregs

**Method:** limma-trend on log2(CPM+0.5) (no voom — input is pre-normalized CPM deposit).
Model: `~0 + group`, trend=TRUE, robust=TRUE. n=5/group. Sample→group labels are PROVISIONAL (no sample sheet confirmed).

## Contrasts (7)

All 7 contrasts are built generically from `analysis_config.yaml:design.contrasts` via `makeContrasts`; no contrast is hard-coded in the script.

| Name | Formula | Biological meaning |
|---|---|---|
| WT_heat | WT\_39 − WT\_37 | Heat response in the cGAS-competent genotype |
| KO_heat | cGASKO\_39 − cGASKO\_37 | Heat response with cGAS removed |
| Interaction | (WT\_39 − cGASKO\_39) − (WT\_37 − cGASKO\_37) | **THE cGAS-dependence test** — 1 df interaction; underpowered at n=5 |
| Geno_at_39 | WT\_39 − cGASKO\_39 | Genotype effect under heat |
| Geno_at_37 | WT\_37 − cGASKO\_37 | Genotype effect at baseline |
| Temp_main | ½(WT\_39+cGASKO\_39) − ½(WT\_37+cGASKO\_37) | Pooled/marginal average heat program — NOT the cGAS-dependence test |
| Geno_main | ½(WT\_37+WT\_39) − ½(cGASKO\_37+cGASKO\_39) | Pooled/marginal average genotype effect |

## Artifacts produced by `02_de_limma_trend_viz.R`

- `figures/by_contrast/<contrast>/volcano.*` — per-contrast volcano plots for all 7 contrasts
- `figures/by_contrast/<contrast>/md.*` — per-contrast MD plots for all 7 contrasts
- `figures/_overview/de_counts_summary.*` — cross-contrast DE-count overview (figure-style)
- `tables/by_contrast/<contrast>/volcano.csv` — source table for each volcano panel (x7)
- `tables/by_contrast/<contrast>/md.csv` — source table for each MD panel (x7)
- `tables/_overview/de_counts_summary.csv` — per-contrast signed DE-count summary table

## Headline panels (bespoke, owner-curated)

- `tables/marker_cgas_dependence.csv` — ISG / HIF marker statistics across WT_heat, KO_heat, Interaction
- `tables/fig2_marker_means.csv` — plot-ready per-group mean log2(CPM+0.5) for Fig 2 markers; **load-bearing cross-stage dep** (`03_decoupler_tf.R:687` reads this; path+schema frozen)

## House framing

- Never phrase as "cGAS-independent" — use "no detectable cGAS-dependence at n=5".
- Do not crown HIF1α/HIF2α as the effector; figure claims are floored at L3 (association/correlational); mechanism stays in prose.
- Interaction is 1 df and underpowered; treat FDR thresholds cautiously.
- n=5/group; sample→group labels are PROVISIONAL pending sample sheet confirmation.


## tables/fig2_marker_means.csv

Per-group mean log2(CPM+0.5) tidy table for 15 cGAS-dependence marker
genes (7 ISG-arm: Ifit1, Isg15, Irf7, Oasl2, Mx1, Stat1, Cxcl10; 8
HIF/glycolysis-arm: Slc2a1, Vegfa, Egln3, Bnip3, Pgk1, Ldha, Aldoa,
Hk2) across 4 groups (WT_37, WT_39, cGASKO_37, cGASKO_39). Also
carries Interaction logFC and adj.P.Val from limma-trend (no
recomputation in viz). READ-ONLY in this script; schema-frozen (also
consumed by 03_decoupler_tf.R:687). Also the HIGHLIGHT watchlist
source for the volcano + MD sweeps. Claim tier: L3. PROVISIONAL
sample labels.

**How to read:** Columns: gene, arm (IFN/ISG arm | HIF/glycolysis arm), hif_class
(HIF-specific | shared-glycolytic | NA for ISG), genotype (WT |
cGASKO), temp (37C | 39C), mean_log2cpm, inter_logFC, inter_adjP. One
row per gene x group (4 rows per gene). inter_adjP is constant within
gene. Schema frozen: do NOT rename/relocate; consumed by downstream
scripts. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `read.csv (read-only; produced by 02_de_limma_trend.R)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds + fig2_marker_means.csv` |

## tables/marker_cgas_dependence.csv

Per-marker Interaction contrast statistics table (subset of
02_de_results.rds, Interaction slot) for the cGAS-dependence marker
genes. Lists logFC, P.Value, adj.P.Val, AveExpr for each marker. Used
as supplemental evidence for the two-arms framing. Claim tier: L3.
PROVISIONAL.

**How to read:** Rows = marker genes; columns = gene_symbol, ensembl, logFC, AveExpr,
t, P.Value, adj.P.Val, contrast. All values from the Interaction
contrast (1 df; logFC > 0 = up in WT relative to cGASKO after
removing shared heat effect). Non-significant adj.P = no detectable
cGAS-dependence at n=5 (NEVER 'cGAS-independent'). Claim tier: L3.
PROVISIONAL sample labels.

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
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. PROVISIONAL
sample labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/WT_heat/volcano.pdf

Print-format companion to `figures/by_contrast/WT_heat/volcano.png`. Same
limma-trend volcano for WT_heat (WT: heat (39 vs 37 °C)): 339 genes
pass adj.P < 0.05 & |log2FC| >= 1.0 (213 up, 126 down). PDF format
for publication; rendered at 9×8 in to accommodate ggrepel labels and
the toolkit-appended legend count line.

**How to read:** Same content as `volcano.png`; PDF format for print. x = log2FC
(>0 = up in numerator); y = -log10(P.Value); orange = up, blue = down,
grey = NS. Significant up/down gene counts (↑ up in numerator, ↓ down)
are shown as a second line directly beneath the highest-priority populated
significance line in the legend — priority order: `FDR ≤ 0.05 & |log2FC| ≥ 1.0`,
else `FDR ≤ 0.05`, else 0/0 when nothing crosses FDR. Sign convention:
positive logFC = up in WT\_39 vs WT\_37. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/WT_heat/volcano.csv

Source table for the WT_heat volcano panel: all 19,679 genes tested by
limma-trend for the WT_heat contrast (WT\_39 − WT\_37), ordered by
adj.P.Val ascending then |logFC| descending. Carries the full per-gene
statistics used to render figures/by_contrast/WT_heat/volcano.png;
same-stem sibling to that figure. Claim tier: L3. PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up in WT\_39 vs WT\_37; <0 = down),
AveExpr (average log2CPM across all samples), t (moderated t-statistic), P.Value (raw),
adj.P.Val (BH-corrected FDR), contrast (constant = "WT_heat"). Significance rule for
figures: adj.P.Val < 0.05 AND |logFC| >= 1.0. n=5/group; 19,679 genes tested.
Claim tier: L3. PROVISIONAL sample labels.

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
columns (no MArrayLM fit on disk). Claim tier: L3. PROVISIONAL sample
labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/WT_heat/md.pdf

Print-format companion to `figures/by_contrast/WT_heat/md.png`. Same
limma-trend MD plot for WT_heat (WT: heat (39 vs 37 °C)); PDF format
for publication at 9×8 in.

**How to read:** Same content as `md.png`; PDF format for print. x = AveExpr
(log2CPM); y = logFC (positive = up in WT\_39 vs WT\_37); orange = up,
blue = down, grey = NS; red dashed = LOESS trend. Claim tier: L3.
PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/WT_heat/md.csv

Source table for the WT_heat MD (mean-difference) plot: all 19,679 genes
tested by limma-trend for the WT_heat contrast (WT\_39 − WT\_37), ordered by
adj.P.Val ascending then |logFC| descending. Same gene-level statistics as
volcano.csv; same-stem sibling to figures/by_contrast/WT_heat/md.png.
Claim tier: L3. PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis in the MD plot; >0 = up in
WT\_39 vs WT\_37), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val,
contrast. The MD panel plots AveExpr (x) vs logFC (y); this table is its
complete underlying data. Significance rule for coloring: adj.P.Val < 0.05.
n=5/group; 19,679 genes tested. Claim tier: L3. PROVISIONAL sample labels.

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
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. PROVISIONAL
sample labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/KO_heat/volcano.pdf

Print-format companion to `figures/by_contrast/KO_heat/volcano.png`. Same
limma-trend volcano for KO_heat (cGAS-KO: heat (39 vs 37 °C)): 392 genes
pass adj.P < 0.05 & |log2FC| >= 1.0 (239 up, 153 down). PDF format for
publication at 9×8 in.

**How to read:** Same content as `volcano.png`; PDF format for print. x = log2FC
(positive = up in cGASKO\_39 vs cGASKO\_37); y = -log10(P.Value);
orange = up, blue = down. Significant up/down gene counts (↑ up in
numerator, ↓ down) are shown as a second line directly beneath the
highest-priority populated significance line in the legend — priority order:
`FDR ≤ 0.05 & |log2FC| ≥ 1.0`, else `FDR ≤ 0.05`, else 0/0 when nothing
crosses FDR. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/KO_heat/volcano.csv

Source table for the KO_heat volcano panel: all 19,679 genes tested by
limma-trend for the KO_heat contrast (cGASKO\_39 − cGASKO\_37), ordered by
adj.P.Val ascending then |logFC| descending. Carries the full per-gene
statistics used to render figures/by_contrast/KO_heat/volcano.png.
Claim tier: L3. PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up in cGASKO\_39 vs cGASKO\_37),
AveExpr (average log2CPM), t, P.Value, adj.P.Val, contrast (constant = "KO_heat").
Significance rule for figures: adj.P.Val < 0.05 AND |logFC| >= 1.0. n=5/group;
19,679 genes tested. Claim tier: L3. PROVISIONAL sample labels.

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
columns (no MArrayLM fit on disk). Claim tier: L3. PROVISIONAL sample
labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/KO_heat/md.pdf

Print-format companion to `figures/by_contrast/KO_heat/md.png`. Same
limma-trend MD plot for KO_heat; PDF format for publication at 9×8 in.

**How to read:** Same content as `md.png`; PDF format for print. x = AveExpr
(log2CPM); y = logFC (positive = up in cGASKO\_39 vs cGASKO\_37);
orange = up, blue = down; red dashed = LOESS trend. Claim tier: L3.
PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/KO_heat/md.csv

Source table for the KO_heat MD plot: all 19,679 genes tested by
limma-trend for the KO_heat contrast (cGASKO\_39 − cGASKO\_37), ordered by
adj.P.Val ascending then |logFC| descending. Same-stem sibling to
figures/by_contrast/KO_heat/md.png. Claim tier: L3. PROVISIONAL
sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis in MD plot; >0 = up in
cGASKO\_39 vs cGASKO\_37), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val,
contrast. Significance rule for coloring: adj.P.Val < 0.05. n=5/group; 19,679 genes
tested. Claim tier: L3. PROVISIONAL sample labels.

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
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. PROVISIONAL
sample labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Interaction/volcano.pdf

Print-format companion to `figures/by_contrast/Interaction/volcano.png`.
Same limma-trend volcano for the Interaction (cGAS-dependence) term: 9
genes pass adj.P < 0.05 & |log2FC| >= 1.0. PDF format for publication
at 9×8 in.

**How to read:** Same content as `volcano.png`; PDF format for print. x = log2FC
(positive = heat response larger in WT vs cGASKO, i.e., cGAS-dependent);
1 df Interaction term; non-significant = no detectable cGAS-dependence at n=5,
NOT independence. Significant up/down gene counts (↑ up in numerator, ↓ down)
are shown as a second line directly beneath the highest-priority populated
significance line in the legend — priority order: `FDR ≤ 0.05 & |log2FC| ≥ 1.0`,
else `FDR ≤ 0.05`, else 0/0 when nothing crosses FDR. Claim tier: L3.
PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Interaction/volcano.csv

Source table for the Interaction volcano panel: all 19,679 genes tested by
limma-trend for the Interaction contrast ((WT\_39 − cGASKO\_39) − (WT\_37 −
cGASKO\_37); the 1 df cGAS-dependence test), ordered by adj.P.Val ascending
then |logFC| descending. This is the definitive gene-level cGAS-dependence
readout from the model; non-significant adj.P.Val = no detectable
cGAS-dependence at n=5 (NEVER 'cGAS-independent'). Claim tier: L3.
PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = the heat response is larger in
WT than cGASKO, i.e., cGAS-dependent induction; <0 = heat response larger in cGASKO),
AveExpr, t, P.Value, adj.P.Val, contrast (constant = "Interaction"). 1 df test;
lowest-powered term in the model at n=5/group. 23 genes pass adj.P.Val < 0.05 (no
logFC gate); 9 pass the combined adj.P.Val < 0.05 AND |logFC| >= 1.0 threshold used
in the volcano figure. Claim tier: L3. PROVISIONAL sample labels.

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
columns (no MArrayLM fit on disk). Claim tier: L3. PROVISIONAL sample
labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Interaction/md.pdf

Print-format companion to `figures/by_contrast/Interaction/md.png`. Same
limma-trend MD plot for the Interaction (cGAS-dependence) term; PDF
format for publication at 9×8 in.

**How to read:** Same content as `md.png`; PDF format for print. x = AveExpr
(log2CPM); y = logFC (positive = cGAS-dependent heat induction);
orange = up, blue = down; red dashed = LOESS trend. Non-sig y-values =
no detectable cGAS-dependence at n=5. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Interaction/md.csv

Source table for the Interaction MD plot: all 19,679 genes tested by
limma-trend for the Interaction contrast (1 df cGAS-dependence test),
ordered by adj.P.Val ascending then |logFC| descending. Same-stem
sibling to figures/by_contrast/Interaction/md.png. Non-significant
adj.P.Val = no detectable cGAS-dependence at n=5, NOT independence.
Claim tier: L3. PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = cGAS-dependent heat induction;
<0 = cGAS-dependent heat suppression), AveExpr (average log2CPM), t, P.Value,
adj.P.Val, contrast. 1 df interaction; underpowered at n=5. Significance rule for
coloring in the MD figure: adj.P.Val < 0.05. Claim tier: L3. PROVISIONAL sample
labels.

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
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. PROVISIONAL
sample labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Geno_at_39/volcano.pdf

Print-format companion to `figures/by_contrast/Geno_at_39/volcano.png`.
Same limma-trend volcano for Geno_at_39 (genotype effect under heat):
17 genes pass adj.P < 0.05 & |log2FC| >= 1.0. PDF format for
publication at 9×8 in.

**How to read:** Same content as `volcano.png`; PDF format for print. x = log2FC
(positive = up in WT\_39 vs cGASKO\_39); orange = up, blue = down.
Significant up/down gene counts (↑ up in numerator, ↓ down) are shown as a
second line directly beneath the highest-priority populated significance line
in the legend — priority order: `FDR ≤ 0.05 & |log2FC| ≥ 1.0`, else
`FDR ≤ 0.05`, else 0/0 when nothing crosses FDR. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_at_39/volcano.csv

Source table for the Geno_at_39 volcano panel: all 19,679 genes tested by
limma-trend for the Geno_at_39 contrast (WT\_39 − cGASKO\_39; genotype effect
under heat stress), ordered by adj.P.Val ascending then |logFC| descending.
Same-stem sibling to figures/by_contrast/Geno_at_39/volcano.png.
Claim tier: L3. PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up in WT vs cGASKO at 39 °C;
<0 = up in cGASKO), AveExpr, t, P.Value, adj.P.Val, contrast (constant = "Geno_at_39").
Significance rule for figures: adj.P.Val < 0.05 AND |logFC| >= 1.0. n=5/group;
64 genes pass adj.P.Val < 0.05 (no logFC gate). Claim tier: L3. PROVISIONAL sample labels.

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
columns (no MArrayLM fit on disk). Claim tier: L3. PROVISIONAL sample
labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Geno_at_39/md.pdf

Print-format companion to `figures/by_contrast/Geno_at_39/md.png`. Same
limma-trend MD plot for Geno_at_39; PDF format for publication at 9×8 in.

**How to read:** Same content as `md.png`; PDF format for print. x = AveExpr
(log2CPM); y = logFC (positive = up in WT vs cGASKO at 39 °C); orange
= up, blue = down; red dashed = LOESS trend. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_at_39/md.csv

Source table for the Geno_at_39 MD plot: all 19,679 genes tested by
limma-trend for the Geno_at_39 contrast (WT\_39 − cGASKO\_39), ordered by
adj.P.Val ascending then |logFC| descending. Same-stem sibling to
figures/by_contrast/Geno_at_39/md.png. Claim tier: L3. PROVISIONAL
sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis; >0 = up in WT vs cGASKO
at 39 °C), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val, contrast.
Significance rule for MD coloring: adj.P.Val < 0.05. n=5/group; 19,679 genes tested.
Claim tier: L3. PROVISIONAL sample labels.

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
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. PROVISIONAL
sample labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Geno_at_37/volcano.pdf

Print-format companion to `figures/by_contrast/Geno_at_37/volcano.png`.
Same limma-trend volcano for Geno_at_37 (baseline genotype effect at 37 °C):
2 genes pass adj.P < 0.05 & |log2FC| >= 1.0. PDF format for publication at
9×8 in.

**How to read:** Same content as `volcano.png`; PDF format for print. x = log2FC
(positive = up in WT\_37 vs cGASKO\_37); orange = up, blue = down.
Significant up/down gene counts (↑ up in numerator, ↓ down) are shown as a
second line directly beneath the highest-priority populated significance line
in the legend — priority order: `FDR ≤ 0.05 & |log2FC| ≥ 1.0`, else
`FDR ≤ 0.05`, else 0/0 when nothing crosses FDR. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_at_37/volcano.csv

Source table for the Geno_at_37 volcano panel: all 19,679 genes tested by
limma-trend for the Geno_at_37 contrast (WT\_37 − cGASKO\_37; baseline genotype
effect at 37 °C), ordered by adj.P.Val ascending then |logFC| descending.
Same-stem sibling to figures/by_contrast/Geno_at_37/volcano.png.
Claim tier: L3. PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up in WT vs cGASKO at
37 °C; <0 = up in cGASKO), AveExpr, t, P.Value, adj.P.Val, contrast (constant =
"Geno_at_37"). Significance rule for figures: adj.P.Val < 0.05 AND |logFC| >= 1.0.
n=5/group; 5 genes pass adj.P.Val < 0.05 (no logFC gate). Claim tier: L3. PROVISIONAL
sample labels.

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
columns (no MArrayLM fit on disk). Claim tier: L3. PROVISIONAL sample
labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Geno_at_37/md.pdf

Print-format companion to `figures/by_contrast/Geno_at_37/md.png`. Same
limma-trend MD plot for Geno_at_37; PDF format for publication at 9×8 in.

**How to read:** Same content as `md.png`; PDF format for print. x = AveExpr
(log2CPM); y = logFC (positive = up in WT vs cGASKO at 37 °C); orange
= up, blue = down; red dashed = LOESS trend. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_at_37/md.csv

Source table for the Geno_at_37 MD plot: all 19,679 genes tested by
limma-trend for the Geno_at_37 contrast (WT\_37 − cGASKO\_37), ordered by
adj.P.Val ascending then |logFC| descending. Same-stem sibling to
figures/by_contrast/Geno_at_37/md.png. Claim tier: L3. PROVISIONAL
sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis; >0 = up in WT vs cGASKO
at 37 °C), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val, contrast.
Significance rule for MD coloring: adj.P.Val < 0.05. n=5/group; 19,679 genes tested.
Claim tier: L3. PROVISIONAL sample labels.

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
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. PROVISIONAL
sample labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Temp_main/volcano.pdf

Print-format companion to `figures/by_contrast/Temp_main/volcano.png`.
Same limma-trend volcano for Temp_main (pooled heat main effect): 346
genes pass adj.P < 0.05 & |log2FC| >= 1.0. PDF format for publication
at 9×8 in.

**How to read:** Same content as `volcano.png`; PDF format for print. x = log2FC
(positive = up at 39 °C vs 37 °C pooled over genotypes; not the
cGAS-dependence term). Significant up/down gene counts (↑ up in numerator,
↓ down) are shown as a second line directly beneath the highest-priority
populated significance line in the legend — priority order:
`FDR ≤ 0.05 & |log2FC| ≥ 1.0`, else `FDR ≤ 0.05`, else 0/0 when nothing
crosses FDR. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Temp_main/volcano.csv

Source table for the Temp_main volcano panel: all 19,679 genes tested by
limma-trend for the Temp_main contrast (½(WT\_39+cGASKO\_39) − ½(WT\_37+cGASKO\_37);
pooled/marginal heat main effect), ordered by adj.P.Val ascending then |logFC|
descending. Same-stem sibling to figures/by_contrast/Temp_main/volcano.png.
This contrast collapses genotype; it is NOT the cGAS-dependence readout (that is
the Interaction contrast). Claim tier: L3. PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up at 39 °C vs 37 °C pooled
over genotypes; <0 = down under heat), AveExpr, t, P.Value, adj.P.Val, contrast
(constant = "Temp_main"). Significance rule for figures: adj.P.Val < 0.05 AND
|logFC| >= 1.0. n=5/group (10 samples total, pooled); 11,153 genes pass adj.P.Val < 0.05
(no logFC gate). Claim tier: L3. PROVISIONAL sample labels.

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
columns (no MArrayLM fit on disk). Claim tier: L3. PROVISIONAL sample
labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Temp_main/md.pdf

Print-format companion to `figures/by_contrast/Temp_main/md.png`. Same
limma-trend MD plot for Temp_main (pooled heat main effect); PDF format
for publication at 9×8 in.

**How to read:** Same content as `md.png`; PDF format for print. x = AveExpr
(log2CPM); y = logFC (positive = up under heat pooled over genotypes;
not the cGAS-dependence readout); orange = up, blue = down; red dashed
= LOESS trend. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Temp_main/md.csv

Source table for the Temp_main MD plot: all 19,679 genes tested by
limma-trend for the Temp_main contrast (pooled heat main effect), ordered by
adj.P.Val ascending then |logFC| descending. Same-stem sibling to
figures/by_contrast/Temp_main/md.png. This contrast is NOT the
cGAS-dependence readout. Claim tier: L3. PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis; >0 = up under heat pooled
over genotypes), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val, contrast.
Significance rule for MD coloring: adj.P.Val < 0.05. n=5/group; 19,679 genes tested.
Claim tier: L3. PROVISIONAL sample labels.

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
limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. PROVISIONAL
sample labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Geno_main/volcano.pdf

Print-format companion to `figures/by_contrast/Geno_main/volcano.png`.
Same limma-trend volcano for Geno_main (pooled genotype main effect):
4 genes pass adj.P < 0.05 & |log2FC| >= 1.0. PDF format for
publication at 9×8 in.

**How to read:** Same content as `volcano.png`; PDF format for print. x = log2FC
(positive = up in WT vs cGASKO pooled over temperatures; not the
cGAS-dependence term). Significant up/down gene counts (↑ up in numerator,
↓ down) are shown as a second line directly beneath the highest-priority
populated significance line in the legend — priority order:
`FDR ≤ 0.05 & |log2FC| ≥ 1.0`, else `FDR ≤ 0.05`, else 0/0 when nothing
crosses FDR. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_standard_volcano (toolkit) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_main/volcano.csv

Source table for the Geno_main volcano panel: all 19,679 genes tested by
limma-trend for the Geno_main contrast (½(WT\_37+WT\_39) − ½(cGASKO\_37+cGASKO\_39);
pooled/marginal genotype main effect), ordered by adj.P.Val ascending then |logFC|
descending. Same-stem sibling to figures/by_contrast/Geno_main/volcano.png.
This contrast collapses temperature; it is NOT the cGAS-dependence readout (that is
the Interaction contrast). Claim tier: L3. PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (>0 = up in WT vs cGASKO pooled
over temperatures; <0 = up in cGASKO), AveExpr, t, P.Value, adj.P.Val, contrast
(constant = "Geno_main"). Significance rule for figures: adj.P.Val < 0.05 AND
|logFC| >= 1.0. n=5/group (10 samples total, pooled); 60 genes pass adj.P.Val < 0.05
(no logFC gate). Claim tier: L3. PROVISIONAL sample labels.

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
columns (no MArrayLM fit on disk). Claim tier: L3. PROVISIONAL sample
labels; n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## figures/by_contrast/Geno_main/md.pdf

Print-format companion to `figures/by_contrast/Geno_main/md.png`. Same
limma-trend MD plot for Geno_main (pooled genotype main effect); PDF
format for publication at 9×8 in.

**How to read:** Same content as `md.png`; PDF format for print. x = AveExpr
(log2CPM); y = logFC (positive = up in WT vs cGASKO pooled over
temperatures; not the cGAS-dependence readout); orange = up, blue =
down; red dashed = LOESS trend. Claim tier: L3. PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `create_MD_plot (toolkit; driven from topTable AveExpr/logFC) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/objects/02_de_results.rds` |

## tables/by_contrast/Geno_main/md.csv

Source table for the Geno_main MD plot: all 19,679 genes tested by
limma-trend for the Geno_main contrast (pooled genotype main effect), ordered by
adj.P.Val ascending then |logFC| descending. Same-stem sibling to
figures/by_contrast/Geno_main/md.png. This contrast is NOT the
cGAS-dependence readout. Claim tier: L3. PROVISIONAL sample labels.

**How to read:** Columns: gene_symbol, ensembl, logFC (y-axis; >0 = up in WT vs cGASKO
pooled over temperatures), AveExpr (x-axis; average log2CPM), t, P.Value, adj.P.Val,
contrast. Significance rule for MD coloring: adj.P.Val < 0.05. n=5/group; 19,679 genes
tested. Claim tier: L3. PROVISIONAL sample labels.

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
PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `geom_col + geom_text (signed bars with printed per-direction counts)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/master/master_de_genes.csv` |

## figures/_overview/de_counts_summary.pdf

Vector companion to `figures/_overview/de_counts_summary.png`, same bars
and the same printed per-direction gene counts. Wide canvas (12×8 in)
so the seven two-line comparison labels sit side by side without
colliding.

**How to read:** Same content as `de_counts_summary.png`, in vector form for print.
Orange above zero = up, blue below zero = down, relative to the
first-named side of each comparison. Gate: adj.P < 0.05 with no
fold-change cut-off. Claim tier: L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `geom_col + geom_text (signed bars with printed per-direction counts) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/master/master_de_genes.csv` |

## tables/_overview/de_counts_summary.csv

Per-comparison tallies of genes tested, genes significant at adj.P <
0.05 with no fold-change cut-off, and the split of those into up and
down. Source table for figures/_overview/de_counts_summary.png, whose
printed bar labels are exactly n_up and n_down. Claim tier: L3.
PROVISIONAL sample labels.

**How to read:** Columns: contrast (config key), n_tested (19,679 genes in every
contrast's topTable), n_sig (adj.P.Val < 0.05, any fold change), n_up and n_down
(that split, relative to the first-named side of the comparison; n_up + n_down =
n_sig). The volcano and MD panels apply the stricter adj.P.Val < 0.05 AND |logFC|
>= 1.0 gate, so their counts run lower. Interaction: n_sig = 23, all up. A
non-significant interaction gene = no detectable cGAS-dependence at n=5, never
independence. Claim tier: L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/02_de_limma_trend_viz.R` | `geom_col + geom_text (signed bars with printed per-direction counts) via save_overview` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1.0; figures.volcano_label_top=10; colors.diverging` | `03_results/master/master_de_genes.csv` |

