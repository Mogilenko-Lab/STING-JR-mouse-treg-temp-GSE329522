# What the 39 °C set looks like to the ontology

The three-lens stage asked how much of `WT_heat_up` falls inside curated lenses. That is a
membership question and it answers 22 of 213 in, the rest out. This stage asks what the same
213 genes look like to GO itself, using semantic similarity over the Biological Process graph.
Nothing here reads expression, so no grouping it produces can have come from an enrichment.

Three questions, in the order a reader meets them.

**Can the ontology see the set?** It places 196 of 213 on curated evidence. The 17 it cannot are
named rather than dropped quietly, and the gap is not neutral with respect to the biology at
issue: `P4ha2` carries no Biological Process term in mouse or in human, and `P4ha2` is one of the
seven `Lombardi2022_HIF` genes in the set. The hypoxia overlap the euler draws as 7 is 6 here,
for a reason that has nothing to do with the data.

**Is it one process or many?** Its mean pairwise similarity is 0.244, against 0.295 to 0.358 for
single-GO-term references and 0.345 and 0.342 for the curated activation and heat-shock lenses.
Measured against a null matched on annotation depth it sits +1.33 SD (p 0.14), which is not
distinguishable from a random set of equally well-studied genes.

Against a uniform null the same set sits +3.27 SD, and that difference is the confound rather
than the biology. The export gate thresholds effect size, so it admits well-studied genes: its
members carry a median 9 Biological Process terms where the expressed background carries 4.
`KO_heat_up` behaves the same way, so this is what a threshold-gated set looks like.

**Can it be split into modules?** No partition is reported anywhere in this stage, and the panel
says why rather than leaving the absence to be inferred. The best split available leaves 93 % of
the genes in one group. A depth-matched null is no better at 96 %, and single-GO-term references
sit at 87 %, so being unsplittable is a property of well-annotated gene sets. Only the uniformly
drawn null separates, at 16 %, because genes carrying one or two terms are semantic islands.

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

Reproducing this is cheap after the first run. `20_semantic_decomposition.R` caches the term
similarity matrix under `03_results/_scratch/`, keyed on the term set, measure, ontology and
GO.db version, so only a genuine change recomputes it.

## figures/_overview/wtheatup_semantic_coherence.png

GO BP places 196 of the 213 genes in the mouse 39 °C-derived up arm.
Its mean pairwise semantic similarity is 0.244, which is +1.33 SD
from a null matched on annotation depth (p 0.14) and so not
distinguishable from a random set of equally well-studied genes,
while single-GO-term references reach 0.318 and the uniform null sits
at 0.220. No gene-level partition is reported: the best available
split leaves 93% of the genes in one group, and the depth-matched
null behaves the same way (96%).

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

