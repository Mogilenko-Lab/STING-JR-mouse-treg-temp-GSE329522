# What the thresholded 39 °C set is made of

This stage answers one question and its three follow-ups, in the order a reader meets them.

I take the mouse `WT_heat` signed-t ranking, the thresholded `WT_heat_up` set the export gate
produced from it, and two curated lenses — a heat-shock core and a TCR/IEG activation panel — and
ask **what the set is made of**. 198 of its 213 genes fall in neither lens.

That raises a second question: is the heat-shock response absent from the 39 °C response? No.
Heating these cells induces the curated `HSR_core` in both genotypes (`WT_heat` NES +2.049, padj
1.9e-5; `KO_heat` NES +2.065, padj 3.4e-5).

A third follows: how can an induced response be missing from the set derived from it? Because the
gate thresholds effect size, not pathway. `WT_heat_up` genes sit at median rank percentile 1.4 %
(median t 10.80) and `HSR_core` at 24.4 % (median t 2.05), against a gate whose deepest admitted
gene sits at 10.21 %. The response moved moderately; the gate kept only extreme movers.

The last panel is the handoff: the mouse gate output is 213 genes, the human set the disease
compartments receive is 199, and it says how much of that 199 each of them can test.

The order to read the panels in: `wtheatup_attribution` → `lens_nes_by_contrast` →
`hsr_rank_position_panel` → `gate_projection_bridge`, with the three `hsr_lens_membership_*` views
corroborating the first of those.

Re-render order matters for the bridge census. Refresh the human compartments that publish the
ranked lists and the JIA `hsr_wtheatup_overlap.csv` first, review/update
`tables/source_hash_manifest.csv`, then re-run `02_analysis/scripts/19_hsr_decomposition.R`. The
stage writes `tables/source_reads_observed.csv` and stops if a present sibling artifact differs from
its pinned hash.

**A standing distinction this whole stage rests on.** A lens enriching in a ranking and that lens
being contained in a thresholded set are two different measurements. Membership here is always
counted over the whole arm, never over a leading edge: classifying only leading-edge genes and
re-testing those subsets would test genes selected because they enriched.

**Honest ceiling.** `HSR_core` is a curated proteotoxic-stress reference and is not fever-specific.
The 37/39 °C contrast is an experimental perturbation in mouse iTregs, so it supports the response
interpretation and licenses nothing about human fever. Read all of this correlatively.

## figures/_overview/wtheatup_attribution.png

198 of the 213 genes in the thresholded `WT_heat_up` set fall in neither curated lens: 12 (5.6 %)
are in `TCR_activation` only, 3 (1.4 %) in `HSR_core` only, and none in both.

**How to read:** One stacked bar over the 213 `WT_heat_up` genes, partitioned by direct set
membership into `HSR_core_only`, `TCR_activation_only`, `shared_both` and `neither`; the segment
labels give count and percent, and the three segments too thin to hold text are called out to the
right of the bar. The subtitle names the question and states the answer. This panel partitions
membership, not rank enrichment — those are separate measurements and the enrichment one is
`lens_nes_by_contrast`. The 93 % remainder is reported as a remainder: it is not named, and it is
evidence for no mechanism. Claim tier: a direct count over frozen sets.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ggplot(geom_col + position_stack)` + `save_figure(overview = TRUE)` | `FIG_CFG`; `colors.okabe_ito`; `figures.label_size` | `03_results/12_hsr_decomp/tables/_overview/wtheatup_attribution.csv` |

## figures/_overview/lens_nes_by_contrast.png

Warming to 39 °C induces the curated heat-shock core in both genotypes and it is the stronger of
the two lenses in all three contrasts (`WT_heat` NES +2.05, FDR 1.9e-05; heat main effect +2.01,
FDR 4.5e-05; `KO_heat` NES +2.07, FDR 3.4e-05).

**How to read:** Horizontal bars are GSEA NES on each contrast's signed-t ranking, grouped by lens;
NES > 0 means the lens is enriched toward genes higher at 39 °C. The NES is printed at the bar end,
the effective set size and FDR in the right-hand column as *n present of n nominal*.

Both lenses are wholly recovered here — 47 of 47 and 66 of 66 — so the three contrasts are
comparable on that count. The human rankings downstream are not, which is what
`gate_projection_bridge` measures. The cGAS-KO row is what makes the induction cGAS-independent at
this design's power. Claim tier: confirmatory for this experimental temperature contrast;
correlative for any fever reading.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ggplot(geom_col + position_dodge, horizontal)` + `save_figure(overview = TRUE)` | `FIG_CFG`; `colors.okabe_ito`; `design.contrast_labels_short` | `03_results/12_hsr_decomp/tables/_overview/lens_nes_by_contrast.csv` |

## figures/_overview/hsr_rank_position_panel.png

The heat-shock core is a moderate responder and the export gate kept only extreme ones: `HSR_core`
genes centre at rank percentile 24.4 % with median t 2.05, roughly 2.4-fold deeper than the 10.21 %
deepest gene the gate admitted, while the gate's own output centres at 1.4 % with median t 10.80.

**How to read:** Each point is one gene; ridgelines show that set's density along the same
`WT_heat` rank-percentile axis, where 0 % is the most-up end. The shaded band is the empirical gate
span, 0.0051 % to 10.21 % — rank 1 to 2,010 of 19,679 — and the dashed line is its deep edge.

Only 213 of the 2,010 genes inside that band passed, because the gate also required logFC ≥ 1; the
deepest gene it admitted is `Cpne6` at t 4.55, logFC 1.03. The `WT_heat_up` row is left-shifted by
construction and is labelled as the gate's output, not as a result. Claim tier: rank geometry over a
frozen ranking and frozen sets.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ggplot(geom_density_ridges + geom_point + gate-span annotate)` + `save_figure(overview = TRUE)` | `FIG_CFG`; `colors.diverging.up`; `colors.okabe_ito`; `figures.label_size` | `03_results/12_hsr_decomp/tables/_overview/hsr_rank_position_panel.csv` |

## figures/_overview/gate_projection_bridge.png

The 213-gene mouse gate output and the 199-gene human set are the same object seen either side of
ortholog projection, and only 17 of that 199 are present in all eighteen human ranked lists that
consume it — so the 199 is not the same 199 in every compartment.

**How to read:** The top half is a funnel of counts on a log axis, because the first and last steps
differ by three orders of magnitude: 19,679 measured genes, 2,010 inside the gate's rank span, 213
through the gate, 199 after mouse-to-human ortholog projection, 17 present everywhere downstream.
Colour marks which species' symbols each row is counted in.

Curated-lens membership is printed beside the two set rows: `HSR_core` 3 of 213 in mouse against a
47-gene lens and 2 of 199 in human against a 56-gene lens, `TCR_activation` 12 of 213 and 12 of 199
against a 66-gene lens in both. The two lenses behave differently under ortholog projection and that
is a real property, not a transcription slip.

`TCR_activation` is a hand-curated panel whose build asserts a strictly 1:1 human-to-mouse map — 66
human symbols to 66 unique mouse symbols, with the script stopping if any gene maps more than once.
The lens is therefore the same size on both sides and the 12 are the same 12 genes
(`ATF3 CD40LG CD83 EGR1 EGR2 EGR3 IER3 IL2 LTA NR4A1 NR4A2 NR4A3`).

`HSR_core` is assembled paralog-complete from Reactome and GO under a 1:many map, so it is 47 mouse
against 56 human. Its mouse three are `Hspa1a`, `Hspa1b`, `Hsph1`, and the `Hspa1a`/`Hspa1b` paralog
pair collapses to one human symbol in the projected set, leaving `HSPA1A` and `HSPH1`. Read the
mouse and human lens counts as counterparts, not as one measurement.

The bottom half spans the fewest to the most of the 199 present in any single ranked list of each
compartment, with the list count and length range printed; the dashed rule at 17 marks the genes
present in every one of them. Recovery is a property of the analysis that produced each ranked list,
not of the biology, and the shallowest compartment — RA synovium at 2,896–5,585 genes per list —
sets the 17. A magnitude comparison across compartments is therefore a different measurement on a
much smaller set. Claim tier: direct counts over frozen sets and published ranked lists.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `patchwork::wrap_plots(lollipop funnel, range bars)` + `save_figure(overview = TRUE)` | `FIG_CFG`; `colors.diverging.up`; `colors.okabe_ito`; `COMPARTMENT_LABELS` | `03_results/12_hsr_decomp/tables/_overview/gate_projection_bridge.csv` |

## figures/_overview/hsr_lens_membership_euler.png

Drawn to scale, the thresholded `WT_heat_up` set sits almost entirely outside both curated lenses:
it shares 3 of the 47 `HSR_core` genes and 12 of the 66 `TCR_activation` genes, and the two lenses
share none with each other.

**How to read:** This corroborates the answer `wtheatup_attribution` gives, and adds the two lens
marginals the bar chart cannot show. Areas and overlaps are fitted to the counts, so geometry is
evidence here, and the residual is on the face: stress 6.1e-19, diagError 4.2e-10, against a 0 that
means exact.

Bold numbers sit inside the region they count. The four multi-set regions get named callouts on
leader lines, because the two real overlaps are slivers too thin to hold text; open circles anchor
the overlaps that exist, open diamonds the two that are empty, each placed where its count would sit
if it were not zero. Every label anchor is solved rather than hand-placed, because a hand-laid
diagram can assert a containment the counts contradict. Claim tier: corroborating.

Two further checks live in the source tables rather than on the face. The human counterpart of the
`WT_heat_up` ∩ `HSR_core` overlap is 2 of 199 after ortholog projection, drawn in
`gate_projection_bridge`. And conditioning either curated lens on the other shifts NES by about
+0.005 (`HSR_core` +0.00496; `TCR_activation` +0.00498), so the two lenses behave as independent
references rather than as one confound wearing two names.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `eulerr::euler(shape = "circle")` + `ggforce::geom_circle()` + `save_figure(overview = TRUE, void = TRUE)` | `FIG_CFG`; `colors.diverging.up`; `colors.okabe_ito`; `EULER_SEED`; `POI_GRID` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

## figures/_overview/hsr_lens_membership_venn.png

The same seven counts in a layout where geometry carries nothing: 198 genes are `WT_heat_up` alone,
3 and 12 are shared with the curated lenses, and both the `HSR_core` with `TCR_activation` region
and the three-way region read 0.

**How to read:** This corroborates the same answer as the Euler panel, under a layout that cannot
mislead. All three circles are equal and equally overlapping by convention, so only the printed
counts vary — which is why it belongs beside the area-proportional panel.

Every one of the seven regions is drawn at conventional size whatever its count, so a printed 0
marks a region that exists and is empty rather than one that was left out. The counts are re-derived
by the renderer from one row per gene, and the script asserts that those rows sum back to 213, 47
and 66 before plotting, so a printed number and a table number cannot drift apart. Claim tier:
corroborating.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ggvenn::ggvenn()` + `save_figure(overview = TRUE, void = TRUE)` | `FIG_CFG`; `colors.diverging.up`; `colors.okabe_ito`; `figures.label_size` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

## figures/_overview/hsr_lens_membership_upset.png

Ranked as intersections, the mouse 39 °C-derived up arm is dominated by its 198 genes in neither
curated lens, and the two zero intersections sit at the end of the ranking as bars of height 0
rather than as absent columns.

**How to read:** This corroborates the same answer again, in the view easiest to audit
count-by-count. Each bar is one intersection and the dot matrix below marks which of the three sets
it belongs to; the horizontal bars on the left are the whole-set sizes (213, 47, 66).

Every combination gets its own column whether or not any gene falls in it, so `HSR_core` with
`TCR_activation` and the three-way region appear as zero-height bars annotated "empty" —
`ComplexUpset` omits the numeric label on a zero bar, so those two are labelled explicitly. This
view stays readable if the comparison ever grows past three sets, where circle layouts stop being
drawable. Claim tier: corroborating.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `ComplexUpset::upset(intersections = <all seven>)` + `save_figure(overview = TRUE)` | `FIG_CFG`; `figures.base_size`; `figures.point_size`; `figures.label_size` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

## tables/hsr_decomp_lens_nes.csv

The GSEA table shows `HSR_core`, `HSR_sensitivity`, and `TCR_activation` enrichment across
`WT_heat`, `Temp_main`, and `KO_heat`.

**How to read:** `nes > 0` means enrichment toward genes higher in the numerator / 39 °C direction
for heat contrasts; `set_size` is the effective size, i.e. the members present in that ranking;
`leading_edge_genes` is semicolon-delimited and is reported, never used to define a subset for
re-testing.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `clusterProfiler::GSEA(by = "fgsea")` | `GSEA_SEED`; `RANK_METRIC`; `GSEA_MAX_SIZE`; `minGSSize=10` | `03_results/objects/02_de_results.rds`; curated HSR/TCR lens RDS files |

## tables/hsr_decomp_wtheatup_attribution.csv

The gene-level attribution table labels each `WT_heat_up` member by direct `HSR_core` and
`TCR_activation` membership.

**How to read:** `HSR_core_only` means membership in `HSR_core` but not `TCR_activation`;
`TCR_activation_only` means the converse; `shared_both` means both; `wt_heat_t > 0` means higher at
39 °C.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `attrib_for()` direct set membership | `ATTR_LEVELS`; `RANK_METRIC` | `03_results/objects/17_signature_sets.rds`; curated HSR/TCR lens RDS files; HSR taxonomy CSV |

## tables/hsr_decomp_summary.csv

The summary table reports the four attribution buckets over the 213-gene `WT_heat_up` denominator
for `HSR_core` and `HSR_sensitivity`.

**How to read:** `fraction` is `n / denominator`, so the `HSR_core_only` row gives 3 / 213 = 1.4 %
for the thresholded `WT_heat_up` set. The wider `HSR_sensitivity` lens moves that to 5 / 213, which
is the sensitivity read on the same count rather than a second result.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `summarise_attrib()` set-membership summary | `ATTR_LEVELS` | `03_results/12_hsr_decomp/tables/hsr_decomp_wtheatup_attribution.csv` |

## tables/hsr_decomp_overlap.csv

The pairwise overlap table shows that `HSR_core` and `TCR_activation` share 0 genes, while
`WT_heat_up` shares 3 with `HSR_core` and 12 with `TCR_activation`.

**How to read:** `n_intersect` is the raw intersection count and `jaccard` is `n_intersect / union`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `utils::combn()` + direct pairwise set overlap | curated lens files | `03_results/objects/17_signature_sets.rds`; curated HSR/TCR lens RDS files |

## tables/hsr_decomp_rank_concordance.csv

The rank summary places `WT_heat_up`, `HSR_core`, and `TCR_activation` in the full `WT_heat`
signed-t ranking.

**How to read:** `median_rank_pct = 0` is the most-up end of the ranking; this table carries a
Wilcoxon diagnostic, but the rank-position panel does not use or display that p-value.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `rank_summary_one()` + `stats::wilcox.test()` | `RANK_METRIC` | `03_results/objects/02_de_results.rds`; curated HSR/TCR lens RDS files |

## tables/hsr_decomp_conditional.csv

The conditional GSEA table shows that removing the other curated lens changes the target-lens NES by
about +0.005.

**How to read:** `delta_nes = nes_cond - nes_uncond`; the small positive deltas support treating
`HSR_core` and `TCR_activation` as separate curated lenses for this membership audit.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | conditional `run_lens_gsea()` on a lens-purged ranking | `GSEA_SEED`; `RANK_METRIC`; `minGSSize=10` | `03_results/objects/02_de_results.rds`; curated HSR/TCR lens RDS files |

## tables/hsr_lens_membership.csv

The shared source of all three membership panels: mouse three-set region counts, the external human
`WT_heat_up`/`HSR_core` overlap, and conditional NES deltas.

**How to read:** `n_mouse` is the plotted mouse count and the seven rows partition the union of the
three sets; `human_wt_hsr_intersect` is the human counterpart for `WT_heat_up` with `HSR_core`,
carried only when the external table is present and well formed. The Euler, Venn, and UpSet panels
all read this one table, so they cannot disagree about a count — only about how a count is drawn.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `overlap_n()` region arithmetic + optional external overlap import | curated lens files | `03_results/12_hsr_decomp/tables/hsr_decomp_overlap.csv`; `../human_treg_arthritis/03_results/10_hsr_lens/tables/hsr_wtheatup_overlap.csv`; `03_results/12_hsr_decomp/tables/hsr_decomp_conditional.csv` |

## tables/_overview/wtheatup_attribution.csv

The source table for `figures/_overview/wtheatup_attribution.png` contains the `HSR_core` rows from
`hsr_decomp_summary.csv` plus the plot labels.

**How to read:** `n` and `fraction` are the plotted quantities over `denominator` 213 and `label` is
the count and percentage printed on or beside each segment.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `dplyr::filter()` + `readr::write_csv()` | `FIG_CFG`; `ATTR_LEVELS` | `03_results/12_hsr_decomp/tables/hsr_decomp_summary.csv` |

## tables/_overview/lens_nes_by_contrast.csv

The source table for `figures/_overview/lens_nes_by_contrast.png` contains `HSR_core` and
`TCR_activation` NES values across the three contrast rankings, with the effective set size the
panel prints.

**How to read:** `n_nominal` is the curated lens size and `n_effective` the members present in that
ranking; they coincide here, and `size_fdr_label` is the string printed in the right-hand column.
`label` is the rounded NES printed at each bar end.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `dplyr::mutate()` + `readr::write_csv()` | `FIG_CFG`; `TERM_LEVELS`; `CONTRAST_LEVELS` | `03_results/12_hsr_decomp/tables/hsr_decomp_lens_nes.csv` |

## tables/_overview/hsr_rank_position_panel.csv

The source table for `figures/_overview/hsr_rank_position_panel.png` contains every plotted gene
position plus the gate-span guard metadata repeated on each row.

**How to read:** `rank_pct` is `100 * rank / 19679`; the guard is that `gate_pass_n` must be 213 and
`gate_deepest_gene` must be `Cpne6` at rank 2,010, which is what proves the plotted gate and the
frozen `WT_heat_up` set are the same object. `gate_min_rank_pct` is 0.00508 % and
`gate_max_rank_pct` 10.2139 %.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `WT_heat` rank join + gate guard (`stop()` on mismatch) | `de_fdr`; `de_logfc`; `RANK_METRIC` | `03_results/03_de/tables/by_contrast/WT_heat/md.csv`; `03_results/12_hsr_decomp/tables/hsr_decomp_rank_concordance.csv` |

## tables/_overview/gate_projection_bridge.csv

The source table for `figures/_overview/gate_projection_bridge.png` reconciles the 213-gene mouse
gate with the 199-gene human set, carries each side's curated-lens membership, and records how much
of the human set is present in each human compartment's ranked lists.

**How to read:** `block` selects the part of the ledger — `funnel` is the five counted steps,
`lens` is curated-lens membership on the mouse and human sides, and `downstream_min` /
`downstream_max` bracket the recovery of the 199 within one compartment. `side` says which species'
symbols a row is counted in, `pct_of_denominator` is `n_genes / denominator`, and `note` names the
provenance of each row. The downstream sweep is a read-only census of every published
`ranked_*.tsv` under the sibling compartments; a ranked list shorter than 1,000 rows is treated as a
partial write and skipped, and a checkout with no siblings writes the funnel and lens blocks alone
and records a census of zero lists.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `bridge_row()` ledger assembly + read-only ranked-list census | `MIN_RANKED_LEN`; `census_roots` | `03_results/human_projection/signatures/WT_heat/WT_heat_up.txt`; curated human HSR/TCR lens RDS files; `../<compartment>/03_results/*/tables/ranked_*.tsv` |

## tables/source_hash_manifest.csv

This manifest pins every sibling artifact the bridge census is allowed to read.

**How to read:** `source_label` names the dependency as used by the stage, `source_path` is relative
to the umbrella checkout, and `sha256` is the required byte hash. If a listed sibling file is present
with a different hash, the compute script stops before updating the bridge tables.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `verify_optional_source_hash()` | pinned SHA-256 | JIA HSR overlap table and human compartment ranked lists |

## tables/source_reads_observed.csv

The bridge census read 19 sibling artifacts, and every present file matched its pinned SHA-256.

**How to read:** One row per optional sibling source considered by the bridge census. `status=read`
means the file existed, matched `source_hash_manifest.csv`, and was read. A standalone clone without
siblings records missing optional sources instead of failing.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition.R` | `record_source()` | `MIN_RANKED_LEN=1000`, pinned SHA-256 | `03_results/12_hsr_decomp/tables/source_hash_manifest.csv` |

## tables/_overview/hsr_lens_membership_euler.csv

The source table for `figures/_overview/hsr_lens_membership_euler.png` records what the
area-proportional fit was asked to draw, what it actually drew, and how far apart those are.

**How to read:** `original_value` is the count from the shared table, `fitted_value` is the area the
fitted circles actually give that region, and `residual` is the difference in gene units;
`euler_stress` and `euler_diag_error` are eulerr's own goodness-of-fit measures, repeated on every
row, where 0 means the layout reproduces the counts exactly. `anchor_x`/`anchor_y` is the solved
label position and `anchor_margin` is its clearance from the nearest region boundary — positive when
the region exists (`region_exists = TRUE`), negative for the two empty regions, where the magnitude
says how far the anchor is from being inside all the required circles.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | `eulerr::euler()` fit diagnostics + `solve_anchors()` label-anchor solve | `EULER_SEED`; `POI_GRID`; `FIG_CFG` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

## tables/_overview/hsr_lens_membership_venn.csv

The source table for `figures/_overview/hsr_lens_membership_venn.png` carries the seven counts with
the flag that the layout is deliberately non-quantitative.

**How to read:** `WT_heat_up`, `HSR_core`, and `TCR_activation` are the logical membership of each
region and `combo` is the same thing as a key; `area_is_quantitative = FALSE` records that nothing
in this panel's geometry should be measured, and `drawn_region_count = 7` records that no region was
dropped for being empty.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | region expansion to one row per gene + `readr::write_csv()` | `FIG_CFG`; `SET_ORDER` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |

## tables/_overview/hsr_lens_membership_upset.csv

The source table for `figures/_overview/hsr_lens_membership_upset.png` gives the seven intersections
in the order they are plotted, with the three whole-set sizes.

**How to read:** `axis_position` is the left-to-right column index, sorted by descending count so
the two `is_empty` rows land last; `set_size_*` repeats the marginal size of each set on every row,
matching the horizontal bars on the left of the panel.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/19_hsr_decomposition_viz.R` | intersection ordering + `readr::write_csv()` | `FIG_CFG`; `SET_ORDER` | `03_results/12_hsr_decomp/tables/hsr_lens_membership.csv` |
