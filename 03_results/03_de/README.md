# 03_de — The 39 °C response in wild-type and cGAS-knockout iTregs

This stage fits the model and reads out the seven contrasts the 2×2 design supports. The input
is the delivered CPM deposit, so the engine is limma-trend on log2(CPM+0.5) under `~0 + group`
with `trend=TRUE` and `robust=TRUE`, five libraries per group, 19,679 genes tested in every
contrast. All seven contrasts are assembled from `analysis_config.yaml:design.contrasts` through
`makeContrasts`.

Two results carry forward. **Heat dominates.** Warming changes 8,723 genes in wild-type and
8,901 in cGAS-knockout cells at FDR 0.05, and the two genotypes warm in step — over the 10,418
heat-responsive genes the knockout response is 0.99× the wild-type response at r = 0.95
(`tables/_overview/cgas_dependence_stats.csv`). **The cGAS-dependent arm is small and
one-sided.** Twenty-three genes clear FDR 0.05 on the interaction, all twenty-three in the same
direction, and every one of them sits below the identity line in
`figures/_overview/heat_response_wt_vs_ko.png` — a weaker heat response once cGAS is gone.

The interaction is a one-degree-of-freedom term at n = 5 and is the least-powered contrast here.
A gene failing it has **no detectable cGAS-dependence at n = 5**, and 23 is a floor.

**Compute** (`02_de_limma_trend.R`, `02b_cgas_dependence_geometry.R`) fits the model, writes the
per-contrast statistics into `03_results/master/` and assembles the four heat-relevant contrasts
into one wide gene-level table. **Visualisation** draws a volcano and an MD plot per contrast,
the cross-contrast count summary, and the effect-versus-effect scatter.

## The seven contrasts

| Name | Formula | What it asks |
|---|---|---|
| `WT_heat` | WT_39 − WT_37 | The heat response with cGAS present. |
| `KO_heat` | cGASKO_39 − cGASKO_37 | The heat response with cGAS removed. |
| `Interaction` | (WT_39 − cGASKO_39) − (WT_37 − cGASKO_37) | **The cGAS-dependence test.** 1 df; least-powered term in the design. |
| `Geno_at_39` | WT_39 − cGASKO_39 | The genotype effect under heat. |
| `Geno_at_37` | WT_37 − cGASKO_37 | The genotype effect at baseline. |
| `Temp_main` | ½(WT_39 + cGASKO_39) − ½(WT_37 + cGASKO_37) | The pooled heat response, genotype collapsed. |
| `Geno_main` | ½(WT_37 + WT_39) − ½(cGASKO_37 + cGASKO_39) | The pooled genotype effect, temperature collapsed. |

Positive `logFC` points to the first-named side. On the interaction, positive means the heat
response is larger in wild-type.

---

## Figures

### `figures/by_contrast/<contrast>/volcano.png` — seven panels, one per contrast

**Per-gene effect against evidence, with the two-arms watchlist always labelled.**
x, log2 fold change, positive toward the contrast numerator; y, −log10 of the raw p-value, which
keeps the per-gene resolution that −log10(FDR) collapses. Colour marks significance decided on
FDR. The dashed horizontal rule is the raw p realising adjusted p 0.05 and the dashed verticals
are |log2FC| = 1. Labels give the top ten genes per side by significance plus a fixed watchlist
of seven interferon-arm and eight HIF/glycolysis-arm genes, so the same genes are findable in
every panel.

Counts at adjusted p < 0.05 **and** |log2FC| ≥ 1: `WT_heat` 339 (213 up), `KO_heat` 392 (239 up),
`Interaction` 9 (9 up), `Geno_at_39` 17 (17 up), `Geno_at_37` 2, `Temp_main` 346 (216 up),
`Geno_main` 4. The 213 up-genes of `WT_heat` become the frozen arm `WT_heat_up`.
*Source* `tables/by_contrast/<contrast>/volcano.csv` · `02_analysis/scripts/02_de_limma_trend_viz.R`.

### `figures/by_contrast/<contrast>/md.png` — seven panels, one per contrast

**Effect against abundance, which exposes any expression-range bias in the normalisation.**
x, average expression (`AveExpr`, log2 CPM scale); y, log2 fold change. Grey draws
non-significant genes, orange up and blue down at FDR 0.05. The red dashed LOESS trend should
hug zero across the abundance range; a tilt would flag a normalisation confound. The dotted
vertical marks median expression, the dashed horizontals |log2FC| = 1, and the quadrant numbers
count significant genes per quadrant. Labels give the top five per side plus the same watchlist
the volcanoes carry.
*Source* `tables/by_contrast/<contrast>/md.csv` · `02_analysis/scripts/02_de_limma_trend_viz.R`.

### `figures/_overview/de_counts_summary.png`

**Heat moves thousands of genes; the cGAS-dependent arm holds twenty-three, all one way.**
x, contrast in config order; y, number of genes at adjusted p < 0.05 with no fold-change cut.
Orange bars rise above zero for genes up in the first-named side and blue bars fall below for
genes down. Each bar end prints its own count, including a measured zero, so the interaction
reads without opening the table: 23 up, 0 down.

This panel's gate is FDR alone. The frozen projection signature adds |log2FC| ≥ 1 and lists
213 up / 126 down for the same `WT_heat` contrast
(`../10_signature/tables/_overview/signature_sizes.csv`, gate `fdr_logfc`) — a stricter gate on
the same statistics.
*Source* `tables/_overview/de_counts_summary.csv` · `02_analysis/scripts/02_de_limma_trend_viz.R`.

### `figures/_overview/heat_response_wt_vs_ko.png`

**The two genotypes warm in step, and the exceptions all fall the same way.**
One point per gene. x, log2 fold change at 39 versus 37 °C in wild-type; y, the same in
cGAS-knockout, on equal scales so the dashed identity line runs at 45°. The dotted lines sit one
log2 unit either side of it.

Glyphs run pale to dark. Pale points carry no detectable genotype difference at n = 5.
Vermillion circles fall with heat in both genotypes and fall further without cGAS; black
triangles rise with heat in wild-type and fall without cGAS. All 23 genes clearing adjusted p
0.05 on the interaction sit **below** the identity line and none above it, so the cGAS-dependent
arm is one-sided. Fourteen keep their direction in both genotypes and nine flip sign. Nine also
clear |log2FC| ≥ 1 and carry bold labels; four of those nine are triangles.

The frozen `../human_projection/` contract exports this arm at both gates — the nine at
`fdr_logfc` and all 23 at `fdr_only` — which is the sensitivity a 1-df term at n = 5 needs.
Read the one-sidedness as a statement about the highlighted arm; the pale cloud straddles the
identity line in both directions.
*Source* `tables/_overview/heat_response_wt_vs_ko.csv` · `02_analysis/scripts/02b_cgas_dependence_geometry_viz.R`.

---

## Tables

### Per contrast — `tables/by_contrast/<contrast>/{volcano,md}.csv`

Fourteen files, the same-stem source of each panel above. One row per gene, all 19,679 tested,
ordered by adjusted p then descending |logFC|. Columns: `gene_symbol`, `ensembl`, `logFC`,
`AveExpr` (average log2 CPM), `t` (moderated), `P.Value`, `adj.P.Val`, `contrast`. The two files
in a directory carry the same statistics and differ in which columns their panel draws.

The volcano gate is adjusted p < 0.05 **and** |logFC| ≥ 1; the MD colouring uses adjusted p
alone. FDR-only counts are `Temp_main` 11,153, `Geno_at_39` 64, `Geno_main` 60, `Interaction`
23, `Geno_at_37` 5.

### Cross-contrast — `tables/_overview/`

| File | What it holds | How to read it |
|---|---|---|
| `de_counts_summary.csv` | One row per contrast: `n_tested` (19,679 everywhere), `n_sig` at adjusted p 0.05 with no fold-change cut, and its `n_up` / `n_down` split. | `n_up + n_down = n_sig`. The bar labels on the summary figure are exactly these two columns. |
| `cgas_dependence_wide.csv` | One row per gene with the four heat-relevant contrasts side by side, joined on `ensembl`. Carries `logFC_`, `adjP_` and `sig_` per contrast, plus the derived classes the scatter draws. | `wt_minus_ko` equals `logFC_Interaction` gene by gene, which is the identity that licenses reading distance from the diagonal as the genotype difference. `heat_responsive`, `cgas_dependent`, `reverses_without_cgas` and `in_stringent_gate` are the flags behind `arm_class`. Of the 23 in the arm, 12 are Hallmark interferon-alpha or -gamma members and 11 are unassigned. |
| `cgas_dependence_stats.csv` | Long format — `metric`, `scope`, `subset`, `value`, `note` — holding every scalar the scatter prints. | A lookup takes all three keys: `pearson_r` appears under two subsets. Correlations and the regression slope are given over all genes and over heat-responsive genes alone, because agreement measured on the whole universe would be carried by unchanged genes at the origin. `max_abs_identity_residual` is 0 to numerical precision. |
| `heat_response_wt_vs_ko.csv` | The plotted slice of `cgas_dependence_wide.csv`, one row per gene, 19,679 rows. | `arm_class` folds `cgas_dependent`, `in_stringent_gate` and `reverses_without_cgas` into the three drawn classes. `hallmark_ifn` is annotation and never a selection rule. |

### Headline marker tables

| File | What it holds | How to read it |
|---|---|---|
| `tables/fig2_marker_means.csv` | Per-group mean log2(CPM+0.5) for fifteen watchlist genes — seven interferon-arm (Ifit1, Isg15, Irf7, Oasl2, Mx1, Stat1, Cxcl10) and eight HIF/glycolysis-arm (Slc2a1, Vegfa, Egln3, Bnip3, Pgk1, Ldha, Aldoa, Hk2) — across the four design cells, with the interaction logFC and adjusted p carried alongside. | Four rows per gene, one per cell; `inter_adjP` is constant within a gene. **Schema-frozen**: `03_decoupler_tf.R` reads this file by path, so keep the location and column names intact. |
| `tables/marker_cgas_dependence.csv` | Interaction-contrast statistics for the same watchlist: `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`. | Positive `logFC` means the induction is larger in wild-type once the shared heat effect is removed. |
