# What the 39 °C set looks like to the ontology

The three-lens stage asked how much of `WT_heat_up` falls inside curated lenses. That is a
membership question and it answers 22 of 213 in, the rest out. This stage asks what the same
213 genes look like to GO itself, using semantic similarity over the Biological Process graph.
Nothing here reads expression, so no grouping it produces can have come from an enrichment.

The published figure takes the graded question membership cannot ask: for a gene no lens contains,
how close does its annotation sit to that lens, measured against background genes annotated to a
similar number of terms. `TCR_activation` has a real penumbra, 23 of 153 scorable genes above the
95th percentile against 7.7 by chance. `HSR_core` has none worth the name at 11 against the same
7.7, so beyond `Hspa1a`, `Hspa1b` and `Hsph1` nothing in this set resembles the heat-shock core
more than chance does.

A third lens, `Lombardi2022_HIF`, is scored by the stage and left off that figure. It is scorable
on 72 genes where the other two reach 153, so putting it on the same axis would invite a
comparison across two denominators. Its 4 against a chance count of 3.6 is too thin to read either
way. It stays in `lens_proximity_per_gene.csv` and `semantic_lens_coverage_loss.csv`, where the
denominator travels with it.

The set is the Lombardi 2022 conserved HIF signature, a pan-cancer consensus derived under hypoxia
in cancer cells, used here as a proximity lens over GO annotation. No stage in these results scores
it as a gene-set database; `mouse_anchor/03_results/06_gsea/README.md` records where it was removed
from and why.

The stage's other readings live in its tables. **Coverage:** GO BP places 196 of the 213 genes on
curated evidence, and `semantic_coverage.csv` names the 17 it cannot. The gap bears on the biology
at issue: `P4ha2` carries no Biological Process term in mouse or in human and is one of the seven
`Lombardi2022_HIF` genes in the set, so the hypoxia overlap the euler draws as 7 is 6 here.

**Coherence:** mean pairwise similarity 0.269, +1.38 SD from a null matched on annotation depth
(p 0.095), so the set is not distinguishable from a random set of equally well-studied genes.
Against a uniform null it sits +4.39 SD, because the export gate thresholds effect size and so
admits well-studied genes: a median 9 Biological Process terms against the background's 4.
`KO_heat_up` behaves the same way, so this is what a threshold-gated set looks like.

**Partition:** none is reported anywhere in this stage. The best split leaves 99 % of the genes in
one group and a depth-matched null is no better at 95 %, so being unsplittable is a property of
well-annotated gene sets rather than a finding about this one.

**Two things to know before trusting any number here.** Semantic proximity says two genes carry
similar annotation. It is not evidence of co-regulation and it is not evidence that a gene takes
part in the process a lens is named for. And the measure saturates for sparsely annotated genes:
a gene with one or two terms either shares a term with the lens, scoring 1, or does not. Those
genes are drawn as crosses, excluded from the counted excess, and their number is printed beside
each lens.

**Which GOSemSim build produced these numbers.** Wang similarity weights each parent edge of the
GO graph by the relationship it carries: is_a 0.8, part_of 0.6, everything else 0.7. GOSemSim
through 2.36.0 keys that weight table on the strings `is_a` and `part_of` and matches it against
its own edge table, which writes the same two relationships as `isa` and `part of`. Every edge
missed both entries and took the 0.7 weight. On Biological Process all 66,947 edges land on that
one weight and none receives 0.8 or 0.6.

GOSemSim 2.39.2 normalises the spellings before the match, and this stage now runs on 2.39.2 from
the project-local R library described in `02_analysis/config/env/r_library.md`. That build is the
Bioconductor devel branch, where the fix landed first.

The stage does not take the version string as evidence. It recovers the is_a weight from the
measure at run time, using a term whose only parent is an is_a edge to the BP root, and stops if
the answer is not 0.8. `semantic_provenance.csv` records what it measured as
`wang_isa_weight_measured`, so that table shows which weighting produced the numbers beside it:
0.801 here, against 0.703 under 2.36.0.

Every similarity on this page is about 10 % higher than its 2.36.0 value, because is_a carries
most of the graph and its weight rose from 0.7 to 0.8. The two builds' full 196 × 196 pairwise
matrices correlate at Pearson 0.994 (Spearman 0.992), with 19,034 of 19,110 gene pairs moving,
so every comparison here ranks the same way at higher absolute values.

Every count in the proximity figure came out the same under both builds. Those counts are
percentiles within a depth band, which a near-monotone rescaling barely disturbs, so they had
room to survive. They are reported as an observation.
`tables/semantic_engine_validation.csv` carries every published number under both builds.

Reproducing this is cheap after the first run. `20_semantic_decomposition.R` caches the term
similarity matrix under `03_results/_scratch/`, keyed on the term set, measure, ontology, GO.db
version and GOSemSim version, so only a genuine change recomputes it. The cache filename carries
the GOSemSim version as well, so the two builds can be compared without each evicting the other.

## figures/_overview/wtheatup_lens_proximity.png

Of the 196 annotated genes in the mouse 39 °C-derived up arm, the two
curated lenses drawn here contain 15. HSR_core holds 3, and of the
153 genes scorable against a depth-matched background 11 sit above
its 95th percentile against a chance count of 7.7; TCR_activation
holds 12, and of the 153 genes scorable against a depth-matched
background 23 sit above its 95th percentile against a chance count of
7.7.

**How to read:** Upper panel: one dot per gene per lens, at the percentile of its
similarity to the nearest other lens member. No gene is scored
against itself, so the red dots control that the measure works. The
dashed rule is the 95th percentile; each axis label gives the
observed count above it beside the count expected by chance, and a
lens whose two counts match has no graded signal beyond its listed
members. Crosses are genes whose depth band already saturates: drawn,
excluded from the counts, tallied on the axis. Lower panel: why the
banding is needed. A third scored lens, Lombardi2022_HIF, is left off
over a different denominator and stays in the stage's per-gene and
coverage-loss tables.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/21_semantic_decomposition_viz.R` | `top-level (pD/pE)` | `semantic.depth_bands = 0,3,6,10,15,22,32,50,Inf, semantic.bg_sample_size = 2000, semantic.measure = Wang, lens_not_drawn = Lombardi2022_HIF` | `03_results/objects/20_semantic_decomp.rds` |

## tables/

All eight are written by `02_analysis/scripts/20_semantic_decomposition.R`, except
`semantic_engine_validation.csv`, which is written by
`02_analysis/scripts/20b_semantic_engine_validation.R`. The one under `_overview/` is the
sibling source table of the figure above and is written by the viz script.

| File | What is in it | How to read it |
|---|---|---|
| `semantic_coverage.csv` | The annotation ladder for `WT_heat_up`: nominal 213, 196 with a curated GO BP term, and the 17 shortfall split into three mutually exclusive classes with the gene names spelled out. | `fraction` is out of the nominal 213. Step 6 is the cost of setting IEA evidence aside and is not part of the shortfall sum. |
| `semantic_lens_coverage_loss.csv` | Per curated lens, how many of its members in the set the ontology can place. | `lost_genes` names the ones an overlap drawn elsewhere counts and this stage has to drop. `Lombardi2022_HIF` loses `P4ha2`. |
| `semantic_coherence_sets.csv` | One row per set: the query, the two sibling gate sets, the three curated lenses, 17 single-GO-term reference draws, 20 uniform nulls, 20 depth-matched nulls. | `mean_sim` is the mean pairwise Wang/BMA similarity within the set. `largest_cluster_frac` is the share of genes the best split leaves in one group, so high means unsplittable. Compare within `family`, never across sizes. |
| `semantic_scale_summary.csv` | The query's mean similarity against each null, as z and as a permutation p, plus the median terms per gene for all three. | The p is over 20 draws and floors at 1/21 = 0.048, so it cannot resolve a small effect. The uniform-null row sizes the annotation-depth confound and supports no claim. |
| `lens_proximity_per_gene.csv` | One row per (gene, lens) for all 196 annotated query genes: similarity to the nearest other lens member, its percentile against a uniform and against a depth-matched background, the gene's depth band, and whether it was counted. | `counted` is FALSE when the band held under 20 background genes or already reached the ceiling at its own 95th percentile. Only counted rows enter the excess. No gene is scored against itself, so `is_lens_member` rows work as a positive control. |
| `query_per_gene_coherence.csv` | Each query gene's mean similarity to the rest of the set, with its term count. | Sorted descending. Read alongside `n_terms`: a gene annotated to more terms sits closer to everything. |
| `semantic_provenance.csv` | Package versions, ontology, measure, combiner, evidence filter, seed, and every set size the run used. | `GOSemSim` records which build produced the numbers in this directory. Any change to it invalidates the absolute values on this page. |
| `semantic_engine_validation.csv` | Every number this stage publishes under GOSemSim 2.36.0 and 2.39.2 side by side, plus probe comparisons on 4,000 GO BP term pairs, 780 `WT_heat_up` gene pairs, and the full 196 × 196 pairwise matrix. | `old_value`/`new_value` are the two builds; `changed` is TRUE when they differ beyond 1e-9. Rows whose `reading` says so hold one between-build statistic and repeat it in both columns, so `delta` is 0 there by construction. |
| `_overview/wtheatup_lens_proximity.csv` | Per-lens summary of the proximity figure, for the two lenses it draws: members in the set, genes scorable, genes above the depth-matched 95th percentile, the chance count, and the saturated-band exclusions. | `n_above_p95` against `n_expected_by_chance` is the whole read. `TCR_activation` 23 against 7.7, `HSR_core` 11 against 7.7. `Lombardi2022_HIF` is scored by the stage but not drawn, so it is absent here and present in `lens_proximity_per_gene.csv`. |

