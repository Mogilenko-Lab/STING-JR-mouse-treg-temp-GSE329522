# What the 39 °C set looks like to the ontology

The three-lens stage asked how much of `WT_heat_up` falls inside curated lenses. That is a
membership question and it answers 22 of 213 in, the rest out. This stage asks what the same
213 genes look like to GO itself, using semantic similarity over the Biological Process graph.
Nothing here reads expression, so no grouping it produces can have come from an enrichment.

Three questions, in the order a reader meets them.

**Can the ontology see the set?** It places 196 of 213 on curated evidence. The 17 it cannot are
named here, and the gap bears on the biology at issue: `P4ha2` carries no Biological Process term
in mouse or in human, and `P4ha2` is one of the seven `Lombardi2022_HIF` genes in the set. The
hypoxia overlap the euler draws as 7 is 6 here, for a reason that has nothing to do with the data.

**Is it one process or many?** Its mean pairwise similarity is 0.269, against 0.324 to 0.385 for
single-GO-term references and 0.374 and 0.365 for the curated activation and heat-shock lenses.
Measured against a null matched on annotation depth it sits +1.38 SD (p 0.095), so it is not
distinguishable from a random set of equally well-studied genes. That p is a permutation p over
20 null draws and cannot resolve below 0.048; one draw of the 20 exceeds the query where two did
under the previous GOSemSim build, and that single draw is the whole of the move from 0.14.

Against a uniform null the same set sits +4.39 SD. The export gate thresholds effect size, so it
admits well-studied genes: its members carry a median 9 Biological Process terms where the
expressed background carries 4. Annotation depth is what separates the two nulls.
`KO_heat_up` behaves the same way, so this is what a threshold-gated set looks like.

**Can it be split into modules?** No partition is reported anywhere in this stage, and the reason
is drawn on the panel. The best split available leaves 99 % of the genes in one group. A
depth-matched null is no better at 95 %, and single-GO-term references sit at 88 %, so being
unsplittable is a property of well-annotated gene sets. Only the uniformly drawn null separates,
at 15 %, because genes carrying one or two terms are semantic islands.

The second figure takes the graded question membership cannot ask: for a gene no lens contains,
how close does its annotation sit to that lens, measured against background genes annotated to a
similar number of terms. `TCR_activation` has a real penumbra, 23 of 153 scorable genes above the
95th percentile against 7.7 by chance. `HSR_core` has none worth the name at 11 against the same
7.7, so beyond `Hspa1a`, `Hspa1b` and `Hsph1` nothing in this set resembles the heat-shock core
more than chance does. `Lombardi2022_HIF` is 4 above 3.6, on only 72 scorable genes, which is too
thin to read either way and is reported as such.

**Two things to know before trusting any number here.** Semantic proximity says two genes carry
similar annotation. It is not evidence of co-regulation and it is not evidence that a gene takes
part in the process a lens is named for. And the measure saturates for sparsely annotated genes:
a gene with one or two terms either shares a term with the lens, scoring 1, or does not. Those
genes are drawn as crosses, excluded from the counted excess, and their number is printed beside
each lens, because for `Lombardi2022_HIF` it is 124 of 196.

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

Every count in the second figure came out the same under both builds. Those counts are
percentiles within a depth band, which a near-monotone rescaling barely disturbs, so they had
room to survive. They are reported as an observation.
`tables/semantic_engine_validation.csv` carries every published number under both builds.

Reproducing this is cheap after the first run. `20_semantic_decomposition.R` caches the term
similarity matrix under `03_results/_scratch/`, keyed on the term set, measure, ontology, GO.db
version and GOSemSim version, so only a genuine change recomputes it. The cache filename carries
the GOSemSim version as well, so the two builds can be compared without each evicting the other.

## figures/_overview/wtheatup_semantic_coherence.png

GO BP places 196 of the 213 genes in the mouse 39 °C-derived up arm.
Its mean pairwise semantic similarity is 0.269, which is +1.38 SD
from a null matched on annotation depth (p 0.095) and so not
distinguishable from a random set of equally well-studied genes,
while single-GO-term references reach 0.346 and the uniform null sits
at 0.238. No gene-level partition is reported: the best available
split leaves 99% of the genes in one group, and the depth-matched
null behaves the same way (95%).

**How to read:** Top, coverage: grey is the set as the gate produced it, green how
much carries a GO Biological Process annotation, and the inset names
what it cannot place. Middle, coherence: each dot is one set's mean
pairwise similarity, the tick a family mean, the dashed rule the set
under test. Both nulls are drawn on purpose, and only the
depth-matched one carries the comparison, because the uniform one is
confounded by the gate admitting well-studied genes. Bottom,
separability: best silhouette over three linkages, k-medoids and k to
15, against the share of genes that split puts in one cluster. A set
that decomposes lands right and low. This one lands top-left, so no
partition is drawn.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/21_semantic_decomposition_viz.R` | `top-level (pA/pB/pC)` | `semantic.measure = Wang, semantic.combine = BMA, semantic.ontology = BP, semantic.n_null_matched = 20, semantic.seed = 20260728` | `03_results/objects/20_semantic_decomp.rds` |

## figures/_overview/wtheatup_lens_proximity.png

Of the 196 annotated genes in the mouse 39 °C-derived up arm, the
curated lenses contain 21. HSR_core holds 3, and of the 153 genes
scorable against a depth-matched background 11 sit above its 95th
percentile against a chance count of 7.7; TCR_activation holds 12,
and of the 153 genes scorable against a depth-matched background 23
sit above its 95th percentile against a chance count of 7.7;
Lombardi2022_HIF holds 6, and of the 72 genes scorable against a
depth-matched background 4 sit above its 95th percentile against a
chance count of 3.6.

**How to read:** Upper panel: each dot is one gene against one lens, at the percentile
of its similarity to the nearest OTHER lens member, scored against
background genes carrying a similar number of terms. No gene is
scored against itself, so the red dots act as a control that the
measure works. The dashed rule is the 95th percentile, and each axis
label gives the observed count above it beside the count expected by
chance. A lens whose two counts match has no graded signal beyond its
listed members. Crosses are genes whose depth band already reaches
the ceiling: drawn, excluded from the counts, tallied on the axis.
Lower panel: why the banding is needed.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/21_semantic_decomposition_viz.R` | `top-level (pD/pE)` | `semantic.depth_bands = 0,3,6,10,15,22,32,50,Inf, semantic.bg_sample_size = 2000, semantic.measure = Wang` | `03_results/objects/20_semantic_decomp.rds` |

## tables/

All eight are written by `02_analysis/scripts/20_semantic_decomposition.R`, except
`semantic_engine_validation.csv`, which is written by
`02_analysis/scripts/20b_semantic_engine_validation.R`. The two under `_overview/` are the
sibling source tables of the two figures above and are written by the viz script.

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
| `_overview/wtheatup_semantic_coherence.csv` | `semantic_coherence_sets.csv` as the coherence figure consumed it. | Same columns; kept as the figure's same-stem neighbour. |
| `_overview/wtheatup_lens_proximity.csv` | Per-lens summary of the proximity figure: members in the set, genes scorable, genes above the depth-matched 95th percentile, the chance count, and the saturated-band exclusions. | `n_above_p95` against `n_expected_by_chance` is the whole read. `TCR_activation` 23 against 7.7, `HSR_core` 11 against 7.7, `Lombardi2022_HIF` 4 against 3.6 on 72 scorable genes. |

