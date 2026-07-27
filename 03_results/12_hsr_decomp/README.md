# HSR Decomposition

I use the mouse `WT_heat` signed-t ranking, the thresholded `WT_heat_up` gene set, and two curated lenses to separate response strength from set membership. The same tables support both statements: heating these cells to 39 °C induces `HSR_core` (`WT_heat` NES +2.049, padj 1.9e-5; `KO_heat` NES +2.065, padj 3.4e-5), while the thresholded `WT_heat_up` gene set is not a proxy for `HSR_core` because it contains 3 of 213 mouse genes and its human projection contains 2 of 199 genes.

## figures/wtheatup_attribution.png

The thresholded `WT_heat_up` set contains 3 `HSR_core` genes, 12 `TCR_activation` genes, and 198 genes in neither lens.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ggplot(geom_col)` + `save_figure()` | `FIG_CFG`; `colors.okabe_ito` | `03_results/12_hsr_decomp/tables/wtheatup_attribution.csv` |

**How to read:** The stacked bar partitions the 213 `WT_heat_up` genes by direct membership: `HSR_core_only`, `TCR_activation_only`, `shared_both`, or `neither`; segment labels give count and percent.

## figures/lens_nes_by_contrast.png

`HSR_core` is enriched toward the 39 °C end of both genotype-stratified heat contrasts (`WT_heat` NES +2.049, padj 1.9e-5; `KO_heat` NES +2.065, padj 3.4e-5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ggplot(geom_col)` + `save_figure()` | `FIG_CFG`; `colors.okabe_ito` | `03_results/12_hsr_decomp/tables/lens_nes_by_contrast.csv` |

**How to read:** Bars show GSEA NES on signed-t rankings; positive NES means the lens is enriched toward genes higher in the numerator / 39 °C direction for heat contrasts, and labels print the NES.

## figures/hsr_rank_position_panel.png

`HSR_core` genes center at 24.4% rank percentile with median t = 2.05, about 2.4-fold deeper than the 10.21% deepest gene admitted by the `WT_heat_up` export gate.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ggplot(geom_density_ridges + geom_point)` + `save_figure()` | `FIG_CFG`; `colors.okabe_ito`; WT_heat gate `FDR < 0.05` and `logFC >= 1` verified in compute | `03_results/12_hsr_decomp/tables/hsr_rank_position_panel.csv` |

**How to read:** Points are individual genes and ridgelines show their density along the WT_heat rank-percentile axis; the shaded empirical gate span runs from rank 1 to rank 2,010 (0.01%-10.21%) and contains all 213 admitted genes, but only 213 of the 2,010 genes inside it passed because the gate also required `logFC >= 1`. The `WT_heat_up` row is labelled as the gate output because its left-shift is circular by construction.

## figures/hsr_lens_membership_euler.png

Drawn to scale, the thresholded `WT_heat_up` set sits almost entirely outside both curated lenses: it shares 3 of the 47 `HSR_core` genes and 12 of the 66 `TCR_activation` genes, and the two lenses share none with each other.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `eulerr::euler(shape = "circle")` + `ggforce::geom_circle()` + `save_figure(void = TRUE)` | `FIG_CFG`; `colors.diverging.up`; `colors.okabe_ito`; `EULER_SEED`; `POI_GRID` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

**How to read:** Areas and overlaps are fitted to the counts, so geometry is evidence here, and the residual is printed in the subtitle: stress 6.1e-19 and diagError 4.2e-10 against a 0 that means exact. Bold numbers sit inside the region they count. The four multi-set regions get named callouts on leader lines, because the two real overlaps are slivers too thin to hold text; open circles anchor the overlaps that exist, open diamonds the two that are empty, each placed where its count would sit if it were not zero.

The human counterpart is 2 of 199 after ortholog projection, and conditioning either curated lens on the other shifts NES by about +0.005 (`HSR_core` +0.00496; `TCR_activation` +0.00498). Note what this does and does not say: 39 °C does induce `HSR_core` here (`WT_heat` NES +2.049, padj 1.9e-5), so the reading is that the thresholded `WT_heat_up` gene set is not a proxy for the heat-shock response — not that the 39 °C response lacks one.

## figures/hsr_lens_membership_venn.png

The same seven counts in a layout where geometry carries nothing: 198 genes are `WT_heat_up` alone, 3 and 12 are shared with the curated lenses, and both the `HSR_core` with `TCR_activation` region and the three-way region read 0.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ggvenn::ggvenn()` + `save_figure(void = TRUE)` | `FIG_CFG`; `colors.diverging.up`; `colors.okabe_ito`; `figures.label_size` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

**How to read:** All three circles are equal and equally overlapping by convention, so only the printed counts vary and the drawing cannot assert a containment the counts contradict — that is exactly why it belongs beside the area-proportional panel. Every one of the seven regions is drawn at conventional size whatever its count, so a printed 0 marks a region that exists and is empty rather than one that was left out. The counts are re-derived by the renderer from one row per gene, and the script asserts that those rows sum back to 213, 47, and 66 before plotting, so a printed number and a table number cannot drift apart.

## figures/hsr_lens_membership_upset.png

Ranked as intersections, the mouse 39 °C-derived up arm is dominated by its 198 genes in neither curated lens, and the two zero intersections sit at the end of the ranking as bars of height 0 rather than as absent columns.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ComplexUpset::upset(intersections = <all seven>)` + `save_figure()` | `FIG_CFG`; `figures.base_size`; `figures.point_size`; `figures.label_size` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

**How to read:** Each bar is one intersection and the dot matrix below it marks which of the three sets that bar belongs to; the horizontal bars on the left are the three whole-set sizes (213, 47, 66). Every combination of the three sets is given its own column whether or not any gene falls in it, so `HSR_core` with `TCR_activation` and the three-way region appear as zero-height bars annotated "empty" — `ComplexUpset` omits the numeric label on a zero bar, so those two are labelled explicitly. This is the view that stays readable if the comparison ever grows past three sets, where circle-based layouts stop being drawable.

## tables/hsr_decomp_lens_nes.csv

The GSEA table shows `HSR_core`, `HSR_sensitivity`, and `TCR_activation` enrichment across `WT_heat`, `Temp_main`, and `KO_heat`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `clusterProfiler::GSEA()` | `GSEA_SEED`; `RANK_METRIC`; `GSEA_MAX_SIZE`; `minGSSize=10` | `03_results/objects/02_de_results.rds`; curated HSR/TCR lens RDS files |

**How to read:** `nes > 0` means enrichment toward genes higher in the numerator / 39 °C direction for heat contrasts; `leading_edge_genes` is semicolon-delimited.

## tables/hsr_decomp_wtheatup_attribution.csv

The gene-level attribution table labels each `WT_heat_up` member by direct `HSR_core` and `TCR_activation` membership.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | direct set membership | `ATTR_LEVELS`; `RANK_METRIC` | `03_results/objects/17_signature_sets.rds`; curated HSR/TCR lens RDS files; HSR taxonomy CSV |

**How to read:** `HSR_core_only` means membership in `HSR_core` but not `TCR_activation`; `TCR_activation_only` means the converse; `shared_both` means both; `wt_heat_t > 0` means higher at 39 °C.

## tables/hsr_decomp_summary.csv

The summary table reports the four attribution buckets over the 213-gene `WT_heat_up` denominator for `HSR_core` and `HSR_sensitivity`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | direct set membership summary | `ATTR_LEVELS` | `03_results/12_hsr_decomp/tables/hsr_decomp_wtheatup_attribution.csv` |

**How to read:** `fraction` is `n / denominator`, so the `HSR_core_only` row gives 3 / 213 = 1.4% for the thresholded `WT_heat_up` set.

## tables/hsr_decomp_overlap.csv

The pairwise overlap table shows that `HSR_core` and `TCR_activation` share 0 genes, while `WT_heat_up` shares 3 with `HSR_core` and 12 with `TCR_activation`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | direct pairwise set overlap | curated lens files | `03_results/objects/17_signature_sets.rds`; curated HSR/TCR lens RDS files |

**How to read:** `n_intersect` is the raw intersection count and `jaccard` is `n_intersect / union`.

## tables/hsr_decomp_rank_concordance.csv

The rank summary places `WT_heat_up`, `HSR_core`, and `TCR_activation` in the full `WT_heat` signed-t ranking.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | rank-position summary | `RANK_METRIC` | `03_results/objects/02_de_results.rds`; curated HSR/TCR lens RDS files |

**How to read:** `median_rank_pct = 0` is the most-up end of the ranking; this table carries a Wilcoxon diagnostic, but the F2 rank-position panel does not use or display that p-value.

## tables/hsr_decomp_conditional.csv

The conditional GSEA table shows that removing the other curated lens changes the target-lens NES by about +0.005.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | conditional `clusterProfiler::GSEA()` | `GSEA_SEED`; `RANK_METRIC`; `minGSSize=10` | `03_results/objects/02_de_results.rds`; curated HSR/TCR lens RDS files |

**How to read:** `delta_nes = nes_cond - nes_uncond`; the small positive deltas support treating `HSR_core` and `TCR_activation` as separate curated lenses for this membership audit.

## tables/wtheatup_attribution.csv

The source table for `figures/wtheatup_attribution.png` contains the HSR_core rows from `hsr_decomp_summary.csv` plus plot labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `read_csv()` + `write_csv()` | `FIG_CFG` | `03_results/12_hsr_decomp/tables/hsr_decomp_summary.csv` |

**How to read:** `label` is the count and percentage printed inside each stacked-bar segment.

## tables/lens_nes_by_contrast.csv

The source table for `figures/lens_nes_by_contrast.png` contains `HSR_core` and `TCR_activation` NES values across the three contrast rankings.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `read_csv()` + `write_csv()` | `FIG_CFG` | `03_results/12_hsr_decomp/tables/hsr_decomp_lens_nes.csv` |

**How to read:** `label` is the rounded NES printed above each bar.

## tables/hsr_rank_position_panel.csv

The source table for `figures/hsr_rank_position_panel.png` contains every plotted gene position plus the gate-span guard metadata repeated on each row.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | WT_heat rank join + gate guard | `de_fdr`; `de_logfc`; `RANK_METRIC` | `03_results/03_de/tables/by_contrast/WT_heat/md.csv`; `03_results/12_hsr_decomp/tables/hsr_decomp_rank_concordance.csv` |

**How to read:** `rank_pct` is `100 * rank / 19679`; `gate_pass_n` must be 213 and `gate_deepest_gene` must be `Cpne6` at rank 2,010.

## tables/hsr_lens_membership.csv

The shared source of all three membership panels: mouse three-set region counts, the external human `WT_heat_up`/`HSR_core` overlap, and conditional NES deltas.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | direct set overlap + external read-only overlap import | curated lens files | `03_results/12_hsr_decomp/tables/hsr_decomp_overlap.csv`; `/workspaces/STING-JR/human_treg_arthritis/03_results/10_hsr_lens/tables/hsr_wtheatup_overlap.csv`; `03_results/12_hsr_decomp/tables/hsr_decomp_conditional.csv` |

**How to read:** `n_mouse` is the plotted mouse count and the seven rows partition the union of the three sets; `human_wt_hsr_intersect` is the human counterpart for `WT_heat_up` with `HSR_core`, carried only when the external table is present and well formed. The Euler, Venn, and UpSet panels all read this one table, so they cannot disagree about a count — only about how a count is drawn.

## tables/hsr_lens_membership_euler.csv

The source table for `figures/hsr_lens_membership_euler.png` records what the area-proportional fit was asked to draw, what it actually drew, and how far apart those are.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `eulerr::euler()` fit diagnostics + label-anchor solve | `EULER_SEED`; `POI_GRID`; `FIG_CFG` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

**How to read:** `original_value` is the count from the shared table, `fitted_value` is the area the fitted circles actually give that region, and `residual` is the difference in gene units; `euler_stress` and `euler_diag_error` are eulerr's own goodness-of-fit measures, repeated on every row, where 0 means the layout reproduces the counts exactly. `anchor_x`/`anchor_y` is the solved label position and `anchor_margin` is its clearance from the nearest region boundary — positive when the region exists (`region_exists = TRUE`), negative for the two empty regions, where the magnitude says how far the anchor is from being inside all the required circles.

## tables/hsr_lens_membership_venn.csv

The source table for `figures/hsr_lens_membership_venn.png` carries the seven counts with the flag that the layout is deliberately non-quantitative.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | region expansion + `read_csv()`/`write_csv()` | `FIG_CFG` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

**How to read:** `WT_heat_up`, `HSR_core`, and `TCR_activation` are the logical membership of each region and `combo` is the same thing as a key; `area_is_quantitative = FALSE` records that nothing in this panel's geometry should be measured, and `drawn_region_count = 7` records that no region was dropped for being empty.

## tables/hsr_lens_membership_upset.csv

The source table for `figures/hsr_lens_membership_upset.png` gives the seven intersections in the order they are plotted, with the three whole-set sizes.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | intersection ordering + `read_csv()`/`write_csv()` | `FIG_CFG` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

**How to read:** `axis_position` is the left-to-right column index, sorted by descending count so the two `is_empty` rows land last; `set_size_*` repeats the marginal size of each set on every row, matching the horizontal bars on the left of the panel.
