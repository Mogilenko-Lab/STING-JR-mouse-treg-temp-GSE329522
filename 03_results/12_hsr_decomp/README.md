# HSR Decomposition

This stage decomposes the mouse `WT_heat` response (WT activated iTregs, 39 C vs 37 C) into three correlative lenses: empirical `WT_heat_up`, curated HSR proteostasis, and TCR/IEG activation. Honest ceiling: even the clean HSR core is proteotoxic-stress-general, not fever-specific; the 37/39 contrast is confirmatory for this temperature perturbation, but it does not make `WT_heat_up` causal for fever.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `build_ranked_vector()` + `clusterProfiler::GSEA()` + direct set attribution | `GSEA_SEED`, `RANK_METRIC`, `GSEA_MAX_SIZE`, `minGSSize=10` | `03_results/objects/02_de_results.rds`; `03_results/objects/17_signature_sets.rds`; HSR/TCR lens RDS files; HSR taxonomy CSV |
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ggplot()` + `save_figure()` | `FIG_CFG`, project theme, Okabe-Ito/diverging colors | `03_results/12_hsr_decomp/tables/hsr_decomp_summary.csv`; `03_results/12_hsr_decomp/tables/hsr_decomp_lens_nes.csv` |

## Tables

### `tables/hsr_decomp_lens_nes.csv`
GSEA scores the HSR core, HSR sensitivity, and TCR/IEG activation lenses across `WT_heat`, `Temp_main`, and `KO_heat`. Read `nes > 0` as enrichment toward genes higher in the numerator / 39 C direction for heat contrasts; leading-edge genes are semicolon-delimited. Claim tier: confirmatory for the frozen contrasts, correlative for fever-cause language.

### `tables/hsr_decomp_wtheatup_attribution.csv`
Each empirical `WT_heat_up` gene is labeled as `thermal_HSR`, `activation`, `shared_both`, or `neither` using HSR_core and TCR_activation membership, with HSR taxonomy category and signed WT_heat t statistic. Read `thermal_HSR` as HSR_core only, `activation` as TCR_activation only, `shared_both` as both, and `neither` as neither; `wt_heat_t > 0` means higher at 39 C. Claim tier: descriptive attribution of a frozen gene set.

### `tables/hsr_decomp_summary.csv`
The headline tally reports counts and fractions of `WT_heat_up` genes in each attribution bucket, first with HSR_core and then with HSR_sensitivity as a robustness lens. Fractions use the empirical `WT_heat_up` size as denominator, and the summary carries the honest-ceiling note. Claim tier: descriptive decomposition, not causal fever attribution.

### `tables/hsr_decomp_overlap.csv`
Pairwise overlaps among `WT_heat_up`, HSR_core, HSR_sensitivity, and TCR_activation show the intersection counts and Jaccard indices, including the HSR_core versus TCR_activation disjointness check. Read larger Jaccard as greater set sharing; this is a set-membership audit, not an enrichment test.

### `tables/hsr_decomp_rank_concordance.csv`
Ranks compare where HSR_core, TCR_activation, and WT_heat_up genes sit in the full WT_heat signed-t ranking, with a Wilcoxon rank-sum p-value for HSR_core versus TCR_activation ranks. `median_rank_pct = 0` is the most-up end of the ranking, and the `higher_set` column names the lens with smaller median rank. Claim tier: competitive rank-position diagnostic.

### `tables/hsr_decomp_conditional.csv`
Conditional GSEA asks whether HSR_core remains enriched after removing all TCR_activation genes from the WT_heat ranking, and symmetrically whether TCR_activation remains after removing HSR_core genes. A `delta_nes` near zero with significant `padj_cond` supports signal not being explained solely by the conditioned lens. Claim tier: competitive/conditional correlative test.

### `tables/wtheatup_attribution.csv`
Source table for `wtheatup_attribution`, copied from the HSR_core rows of the headline summary with plot labels. Read counts and fractions as the four attribution buckets over empirical `WT_heat_up`; claim tier matches `hsr_decomp_summary.csv`.

### `tables/lens_nes_by_contrast.csv`
Source table for `lens_nes_by_contrast`, restricted to HSR_core and TCR_activation NES across `WT_heat`, `Temp_main`, and `KO_heat`. Read positive NES as enrichment toward genes higher in the numerator / 39 C direction for heat contrasts; claim tier is confirmatory for the contrast and correlative for fever/thermal-cause language.

## Figures

### `figures/wtheatup_attribution.pdf` and `figures/wtheatup_attribution.png`
The stacked bar shows how much of empirical `WT_heat_up` falls into thermal HSR, activation, shared, and neither buckets. Segment labels give count and percent; attribution definitions match `hsr_decomp_wtheatup_attribution.csv`, and the source table is `tables/wtheatup_attribution.csv`. Claim tier: descriptive/correlative decomposition.

### `figures/lens_nes_by_contrast.pdf` and `figures/lens_nes_by_contrast.png`
Grouped bars compare signed NES for HSR_core and TCR_activation across the three context contrasts. Bars above zero are enriched toward the up/numerator end of the signed-t ranking, and the source table is `tables/lens_nes_by_contrast.csv`. Claim tier: confirmatory for the WT 39-vs-37 contrast, correlative for fever/thermal-cause language.
