# 13_semantic_decomp — What the 39 °C set looks like to the ontology

[`../12_hsr_decomp/`](../12_hsr_decomp/) asked how much of `WT_heat_up` falls inside curated
lenses. That is a membership question and its answer is 22 of 213 in, the rest out. This stage
asks what the same 213 genes look like to the Gene Ontology itself, using semantic similarity
over the Biological Process graph. Nothing here reads expression, so no grouping it produces can
have come from an enrichment.

The published figure takes the graded question membership cannot ask: for a gene no lens
contains, how close does its annotation sit to that lens, measured against background genes
annotated to a similar number of terms? `TCR_activation` has a real penumbra — 23 of 153
scorable genes above the 95th percentile against 7.7 by chance. `HSR_core` has 11 against the
same 7.7, so beyond `Hspa1a`, `Hspa1b` and `Hsph1` this set resembles the heat-shock core about
as closely as chance does.

A third lens, `Lombardi2022_HIF`, is scored by the stage and left off the figure. It is scorable
on 72 genes where the other two reach 153, so putting it on the same axis would invite a
comparison across two denominators; its 4 against a chance count of 3.6 is thin either way. It
stays in `lens_proximity_per_gene.csv` and `semantic_lens_coverage_loss.csv`, where the
denominator travels with it. The set is the Lombardi 2022 conserved HIF signature, a pan-cancer
consensus derived under hypoxia in cancer cells, used here as a proximity lens over GO
annotation and scored as a gene-set database by no stage in these results.

## What the tables return

**Coverage.** GO BP places 196 of the 213 genes on curated evidence, and `semantic_coverage.csv`
names the 17 it cannot. The gap bears on the biology at issue: `P4ha2` carries no Biological
Process term in mouse or in human and is one of the seven `Lombardi2022_HIF` genes in the set, so
the hypoxia overlap the euler panels draw as 7 is 6 here.

**Coherence.** Mean pairwise similarity is 0.269, +1.38 SD from a null matched on annotation
depth (p 0.095), so the set behaves like a random set of equally well-studied genes. Against a
uniform null it sits +4.39 SD, because the export gate thresholds effect size and so admits
well-studied genes: a median 9 Biological Process terms against the background's 4. `KO_heat_up`
behaves the same way, so this is the shape of a threshold-gated set.

**Partition.** None is reported. The best split leaves 99% of the genes in one group and a
depth-matched null does no better at 95%, so being unsplittable is a property of well-annotated
gene sets.

## Two properties of the measure to hold before reading a number

Semantic proximity says two genes carry similar annotation. Co-regulation is a separate
question, and so is whether a gene takes part in the process a lens is named for.

The measure saturates for sparsely annotated genes: a gene with one or two terms either shares a
term with the lens, scoring 1, or does not. Those genes draw as crosses, sit outside the counted
excess, and their number prints beside each lens.

## Which build produced these numbers

Wang similarity weights each parent edge of the GO graph by its relationship: `is_a` 0.8,
`part_of` 0.6, everything else 0.7. GOSemSim through 2.36.0 keys that weight table on the strings
`is_a` and `part_of` and matches it against its own edge table, which writes the same two
relationships as `isa` and `part of`. Every edge missed both entries and took 0.7. On Biological
Process all 66,947 edges land on that one weight.

GOSemSim 2.39.2 normalises the spellings before the match, and this stage runs on 2.39.2 from
the project-local R library described in `02_analysis/config/env/r_library.md`.

The stage takes the measurement rather than the version string. It recovers the `is_a` weight at
run time from a term whose only parent is an `is_a` edge to the BP root, and stops if the answer
misses 0.8. `semantic_provenance.csv` records what it measured as `wang_isa_weight_measured`:
0.801 here, against 0.703 under 2.36.0.

Every similarity on this page runs about 10% above its 2.36.0 value, because `is_a` carries most
of the graph. The two builds' full 196 × 196 matrices correlate at Pearson 0.994 (Spearman
0.992) with 19,034 of 19,110 pairs moving, so every comparison ranks the same way at higher
absolute values. Every count in the proximity figure came out the same under both builds;
`semantic_engine_validation.csv` carries every published number side by side.

Reproducing this is cheap after the first run. `20_semantic_decomposition.R` caches the term
similarity matrix under `03_results/_scratch/`, keyed on the term set, measure, ontology, GO.db
version and GOSemSim version, and the filename carries the GOSemSim version so the two builds
can be compared without evicting each other.

---

## Figures

### `figures/_overview/wtheatup_lens_proximity.png`

**How close the arm's genes sit to two curated lenses, against a depth-matched background.**
Upper panel: one dot per gene per lens, at the percentile of its similarity to the nearest other
lens member. No gene is scored against itself, so the red dots of the lens members work as a
positive control on the measure. The dashed rule is the 95th percentile, and each axis label
gives the observed count above it beside the count expected by chance. Crosses are genes whose
depth band already saturates: drawn, excluded from the counts, tallied on the axis. Lower panel:
the banding that makes the comparison possible.

Of the 196 annotated genes in the arm, the two drawn lenses contain 15. `HSR_core` holds 3, and
11 of 153 scorable genes sit above its 95th percentile against 7.7 by chance.
`TCR_activation` holds 12, and 23 of the same 153 sit above its 95th percentile against the same
7.7. A lens whose two counts match carries no graded signal beyond its listed members.
*Source* `tables/_overview/wtheatup_lens_proximity.csv` ·
`02_analysis/scripts/21_semantic_decomposition_viz.R`.

---

## Tables

All eight are written by `20_semantic_decomposition.R`, except `semantic_engine_validation.csv`
(`20b_semantic_engine_validation.R`) and the `_overview/` one (the viz script).

| File | What is in it | How to read it |
|---|---|---|
| `semantic_coverage.csv` | The annotation ladder for `WT_heat_up`: nominal 213, 196 with a curated GO BP term, and the 17 shortfall split into three mutually exclusive classes with the gene names spelled out. | `fraction` is out of the nominal 213. Step 6 is the cost of setting IEA evidence aside and sits outside the shortfall sum. |
| `semantic_lens_coverage_loss.csv` | Per curated lens, how many of its members in the set the ontology can place. | `lost_genes` names the ones an overlap drawn elsewhere counts and this stage drops. `Lombardi2022_HIF` loses `P4ha2`. |
| `semantic_coherence_sets.csv` | One row per set: the query, the two sibling gate sets, the three curated lenses, 17 single-GO-term reference draws, 20 uniform nulls, 20 depth-matched nulls. | `mean_sim` is the mean pairwise Wang/BMA similarity within the set. `largest_cluster_frac` is the share the best split leaves in one group, so high means unsplittable. Compare within `family`. |
| `semantic_scale_summary.csv` | The query's mean similarity against each null, as z and as a permutation p, plus the median terms per gene for all three. | The p is over 20 draws and floors at 1/21 = 0.048, so it resolves a small effect no further. The uniform-null row sizes the annotation-depth confound. |
| `lens_proximity_per_gene.csv` | One row per (gene, lens) for all 196 annotated query genes: similarity to the nearest other lens member, its percentile against a uniform and a depth-matched background, the depth band, and whether it was counted. | `counted` is FALSE when the band held under 20 background genes or already reached the ceiling at its own 95th percentile. Only counted rows enter the excess. `is_lens_member` rows are the positive control. |
| `query_per_gene_coherence.csv` | Each query gene's mean similarity to the rest of the set, with its term count. | Sorted descending. Read alongside `n_terms`: a gene annotated to more terms sits closer to everything. |
| `semantic_provenance.csv` | Package versions, ontology, measure, combiner, evidence filter, seed, and every set size the run used. | `GOSemSim` records which build produced the numbers here, and `wang_isa_weight_measured` records what the run measured rather than what the version claims. |
| `semantic_engine_validation.csv` | Every number this stage publishes under GOSemSim 2.36.0 and 2.39.2 side by side, plus probe comparisons on 4,000 GO BP term pairs, 780 `WT_heat_up` gene pairs, and the full 196 × 196 matrix. | `old_value` / `new_value` are the two builds and `changed` is TRUE beyond 1e-9. Rows whose `reading` says so hold one between-build statistic and repeat it in both columns, so `delta` is 0 there by construction. |
| `_overview/wtheatup_lens_proximity.csv` | Per-lens summary of the figure, for the two lenses it draws: members in the set, genes scorable, genes above the depth-matched 95th percentile, the chance count, and the saturated-band exclusions. | `n_above_p95` against `n_expected_by_chance` is the whole read — `TCR_activation` 23 against 7.7, `HSR_core` 11 against 7.7. `Lombardi2022_HIF` is scored by the stage and absent here; it sits in `lens_proximity_per_gene.csv`. |
