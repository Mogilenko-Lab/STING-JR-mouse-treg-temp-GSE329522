# Which biology the projected up arms organise into, asked of the whole ontology

The three-lens decomposition and the human compartment's arm decomposition both ask how much
of `WT_heat_up` falls inside a hand-picked panel of curated programs. A panel that already
contains hypoxia can answer "is hypoxia in this arm", and two thirds of the arm falls in no
panel member at all. This directory asks the ontology instead: over the whole GO Biological
Process space, and separately the whole Molecular Function space, which terms does each arm's
gene content organise into. The candidate space is fixed by GO rather than by whoever is asking.

**What this directory is, and what it is not.** It describes the biology the arm's gene
content organises into. It is **not** a composition of the arm, and its own count-level null is
why. Over-representation gives a gene as many memberships as it has terms: the 384 enriched GO
BP terms for `WT_heat_up` hand out 3,597 gene-term memberships over a union of 160 genes, and 10
of the 202 genes sit in exactly one enriched term. Collapsed to the 36 similarity blocks it is
974 memberships with 22 genes in one block. Shares taken off that do not add up to the arm. The
partition question is asked off frozen MSigDB collections in `16_arm_composition`.

**Read `arms_observed_vs_null.png` before anything else.** The gate that produced these arms
thresholds effect size, which preferentially admits well-studied genes, and the hypergeometric
conditions on set size and universe only. Every count here is therefore compared against 2000
random gene sets matched to the arm band by band on GO BP terms per gene. For `WT_heat_up` the
384 enriched terms sit against a matched-null median of 230, 95th percentile 370 and maximum
603, permutation p 0.037. The 160 genes those terms cover sit against a null median of 147, 95th
percentile 159 and maximum 167, permutation p **0.040**. A uniform null of the same size reaches
a median of 0 terms.

That coverage p is the number to hold onto. "GO BP places 160 of the 202 genes in an enriched
term" is the headline this directory would otherwise carry, and it is about what a random set of
equally well-studied genes achieves — the null's own maximum, 167, is above the observed 160.
The other three arms behave differently: `KO_heat_up` covers 166 against a null median of 129
(p 0.001), and both `Interaction` arms clear their nulls comfortably because a 7-gene or 19-gene
draw usually reaches nothing at all. The arm this project cares about is the weakest of the four
on both counts.

**Recurrence, and why the head of the list is the wrong place to look.** For each enriched term
`wtheatup_null_recurrence.png` gives the share of matched draws that also take it to
significance. The median term recurs in 10.0% of draws, 126 of 384 terms recur in more than 20%
and 24 in more than 50%. Recurrence is highest among the terms a reader would quote first:
chemotaxis 0.274, angiogenesis 0.367, blood vessel morphogenesis 0.458, positive regulation of
cell migration 0.574, positive regulation of cytokine production 0.612. A term at 0.6 is one
the ontology hands to three of every five well-annotated gene sets of this size.

**The proteostasis probe.** Nineteen distinct terms a heat-shock reading would predict were
looked up by id, plus a name sweep over the enriched list. Fifteen entered the hypergeometric
and two cleared the adjusted-p cutoff: `GO:0031649 heat generation` (IL1A, PTGS2, TNFSF11) and
`GO:0001659 temperature homeostasis` (EGR1, GATM, IGF2BP2, IL1A, KDM6B, NPR3, PTGS2, TNFSF11,
p_matched 0.103). Both are carried by inflammatory mediators. Every chaperone, folding and
unfolded-protein term that entered came back above the cutoff, and the chaperone genes present
in the arm are HSPA1A and HSPH1. The four terms that never entered carry their reason in the
`status` column rather than a bare FALSE.

**Which term is tested, and which is only named.** `go_vs_curated_lenses.csv` pairs each
curated lens with the GO term whose name a reader would take to mean the same thing.
`GO:0006954 inflammatory response` holds 689 genes in the with_iea universe, above the 500-gene
cap, so it never entered the test at all while claiming 43 of the arm's 202 genes, the largest
claim of any term in the panel. Under `no_iea` it falls to 486 genes, enters, and comes back
5th of 225 for `WT_heat_up` and 2nd for `KO_heat_up`. Both variants sit on the same row.

**Two IEA variants, and a mouse-space replicate.** `with_iea` is primary because it is the
canonical one-line clusterProfiler path, and the study-depth confound that motivates dropping
IEA is handled by the depth-matched null instead of by an evidence filter. The two variants
share 182 of the 384 with_iea BP terms (Jaccard 0.43) and only 4 of their top 10, so the
evidence filter moves the headline list. The same question asked in mouse space before the
ortholog map returns 456 BP terms sharing 313 with the human-space 384 (Jaccard 0.59).

**Provenance.** GOSemSim 2.39.2 from the project-local R library, GO.db 3.22.0,
org.Hs.eg.db 3.22.0, clusterProfiler 4.18.4. Two seams are checked and recorded in
`go_provenance.csv`: the hand-built annotation map reproduces `enrichGO()` exactly on the
primary arm (same term ids, max |dp| 0), and the sparse-matrix recomputation the null runs on
reproduces `enricher()` term by term (max |dk| 0, max |dq| 0).

## figures/_overview/wtheatup_term_blocks.png

The 384 GO BP terms enriched for WT_heat_up collapse to 36 blocks.
The largest is blood vessel morphogenesis with 56 terms over 81 of
the arm's 202 genes, and 61% of its terms have a depth-matched
p_matched under 0.05. Across all blocks the union of genes the
enriched terms hold is 160. The shaded row is the response to hypoxia
block: 3 terms over 12 arm genes (ADAM8, AK4, EGR1, FOSL2, HK2,
HSPG2, INHBA, LTA, NPPC, NR4A2, PAK1, PTGS2), best adjusted p 0.0118,
and 0% of its terms with a depth-matched p_matched under 0.05.

**How to read:** The 384 enriched terms collapse to 36 blocks at Wang similarity 0.7,
average linkage cut at 0.7 (GOSemSim 2.39.2), one row each. Bar
length is the number of enriched terms in the block, fill the
smallest adjusted p among them, the bar-end text the arm genes those
terms cover between them, out of the arm's 202. The right panel is
the same block's depth-matched verdict: the share of its terms whose
observed hit count a matched random draw equals or beats in under 5%
of 2000 replicates. A long bar with a low dot is a large block a
random draw reproduces. The shaded row spanning both panels is the
response to hypoxia block, shaded so a reader can find it among the
36. Blocks are a similarity cut over term annotation, so a gene can
sit in several and the gene counts do not sum to the arm. Claim tier:
descriptive, hypothesis-generating.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/25_go_decomposition_viz.R` | `top-level (p1a | p1b)` | `go_decomposition.primary_ontology = BP, go_decomposition.primary_iea = with_iea, go_decomposition.p_cutoff = 0.05, go_decomposition.min_gs_size = 10, go_decomposition.max_gs_size = 500, go_decomposition.n_null = 2000, go_decomposition.seed = 20260730` | `03_results/objects/24_go_decomposition.rds` |

## figures/_overview/arms_observed_vs_null.png

For WT_heat_up the 384 enriched GO BP terms sit against a
depth-matched null median of 230 and 95th percentile of 370,
permutation p 0.037 over 2000 replicates. The 160 arm genes in at
least one enriched term sit against a null median of 147 and 95th
percentile of 159, permutation p 0.04. A uniformly drawn null of the
same size reaches a median of 0 terms, so annotation depth accounts
for most of the term count before any biology is read.

**How to read:** Grey is the null over 2000 depth-matched draws, the red rule the
arm's own value, the blue triangle the median of a uniformly drawn
null of the same size. The permutation p is the share of draws
reaching the observed value or more, floored at 1/2001. Rows are arms
and carry the genes drawn per replicate, which is the arm's annotated
size. Left column: terms clearing the adjusted-p gate. Right column:
drawn genes landing in one of those terms, the null for the coverage
number the rest of this directory rests on. Claim tier: this panel
bounds what the rest of the stage can support.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/25_go_decomposition_viz.R` | `top-level (p2)` | `go_decomposition.primary_ontology = BP, go_decomposition.primary_iea = with_iea, go_decomposition.p_cutoff = 0.05, go_decomposition.min_gs_size = 10, go_decomposition.max_gs_size = 500, go_decomposition.n_null = 2000, go_decomposition.seed = 20260730` | `03_results/objects/24_go_decomposition.rds` |

## figures/_overview/wtheatup_null_recurrence.png

Over the 384 GO BP terms enriched for WT_heat_up, a depth-matched
random draw of the same size reaches q below 0.05 for the median term
in 10.0% of 2000 replicates, for 126 terms in more than 20% and for
24 terms in more than 50%. Recurrence rises towards the head of the
list: angiogenesis 0.367, blood vessel morphogenesis 0.458, positive
regulation of cytokine production 0.612, chemotaxis 0.274, positive
regulation of cell migration 0.574.

**How to read:** One dot per enriched term, at its rank by adjusted p against how
often a depth-matched draw of the same size also takes it to
significance over 2000 replicates. Blue dots are terms a matched draw
equals or beats in under 5% of replicates. The red curve is a loess
fit and the right panel is the same y axis as a histogram over all
terms. The dashed rule is the median recurrence, 0.100, and the
dotted rules are 20% and 50%: 126 terms sit above 20% and 24 above
50%. The 8 most significant terms are named, and so are the three
hypoxia and oxygen-level terms, in bold at ranks 161, 189, 250. All
three carry a p_matched at or above 0.05. A term high on this axis is
one the ontology hands to any well-annotated gene set of this size.
Claim tier: a property of the term and of the background, bounding
how the block panel reads.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/25_go_decomposition_viz.R` | `top-level (p3a | p3b)` | `go_decomposition.primary_ontology = BP, go_decomposition.primary_iea = with_iea, go_decomposition.p_cutoff = 0.05, go_decomposition.min_gs_size = 10, go_decomposition.max_gs_size = 500, go_decomposition.n_null = 2000, go_decomposition.seed = 20260730` | `03_results/objects/24_go_decomposition.rds` |

## figures/_overview/wtheatup_proteostasis_probe.png

This is the direct by-name test of whether proteostasis and
heat-shock terms enrich in WT_heat_up. Every term on the configured
probe list is looked up and reported here, whether or not it returned
anything. Every chaperone, protein-folding, unfolded-protein and
heat-shock term that entered the hypergeometric returned an adjusted
p above the 0.05 cutoff, and the chaperone genes present in the arm
are HSPA1A and HSPH1. Of 19 distinct probe terms, 15 entered and 2
cleared the cutoff: heat generation (3 genes: IL1A, PTGS2, TNFSF11,
p_matched 0.011) and temperature homeostasis (8 genes: EGR1, GATM,
IGF2BP2, IL1A, KDM6B, NPR3, PTGS2, TNFSF11, p_matched 0.103). Both of
those are organism-level thermoregulation terms, and the genes
carrying them are inflammatory mediators.

**How to read:** One row per probe term, split by ontology. The dot is the raw
hypergeometric p on a -log10 axis and the dashed rule is p = 0.05.
Grey is tested and above the adjusted-p cutoff, vermillion tested and
below it, a blue cross at zero a term the test never saw. 15 of the
19 probe terms entered. The right panel names the arm genes each
tested term holds, so a term clearing the cutoff can be read by its
gene content, and for an untested term gives the reason it never
entered: above the 500-gene cap, below the 10-gene floor, or holding
no arm gene. Those reasons are properties of the term's size in the
13128-symbol background or of its gene content and carry no
information about the arm. The `status` and `found_by` columns of the
source table hold the same values. Claim tier: these are
searched-and-absent rows, reported so an absence carries a number.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/25_go_decomposition_viz.R` | `top-level (p4a | p4b)` | `go_decomposition.primary_ontology = BP, go_decomposition.primary_iea = with_iea, go_decomposition.p_cutoff = 0.05, go_decomposition.min_gs_size = 10, go_decomposition.max_gs_size = 500, go_decomposition.n_null = 2000, go_decomposition.seed = 20260730` | `03_results/objects/24_go_decomposition.rds` |

## figures/_overview/arms_coverage_ladder.png

160 of WT_heat_up's 202 genes sit in at least one enriched GO BP
term, 31 carry an annotation and sit in none, and 11 carry no BP
annotation at all (7 of those are visible in MF). A depth-matched
random draw of 191 genes covers a median of 147, so the permutation p
on the coverage number is 0.04.

**How to read:** One stacked bar per arm, in gene counts. Green is genes held by at
least one term that reached significance, orange is genes the
ontology can place that no such term holds, grey is genes with no
annotation in this ontology. Segment counts are printed inside the
bar when the segment is wide enough. The black tick is the median
coverage over 2000 depth-matched draws of the same size, and the text
at the bar end carries the arm total, its coverage, that median and
the permutation p. Read the tick before the green segment: coverage
this high is what the ontology gives any set of equally
well-annotated genes. Claim tier: descriptive.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/25_go_decomposition_viz.R` | `top-level (p5)` | `go_decomposition.primary_ontology = BP, go_decomposition.primary_iea = with_iea, go_decomposition.p_cutoff = 0.05, go_decomposition.min_gs_size = 10, go_decomposition.max_gs_size = 500, go_decomposition.n_null = 2000, go_decomposition.seed = 20260730` | `03_results/objects/24_go_decomposition.rds` |

## tables/

Thirteen tables are written by `02_analysis/scripts/24_go_decomposition.R`. The five under
`_overview/` whose stems match a figure above are that figure's source table and are written by
`02_analysis/scripts/25_go_decomposition_viz.R`. Everything lives in `tables/_overview/`.

Two conventions run through all of them. `p_matched` is the one-sided permutation p from 2000
depth-matched draws: the share of draws whose hit count for that term reaches the observed one,
so small means the arm beats matched-random. `frac_matched_reaching_q` is recurrence: the share
of those same draws that take the term to q below 0.05 on their own, so large means the term is
cheap to enrich. A term wants both a small `p_matched` and a small recurrence.

| File | What is in it | How to read it |
|---|---|---|
| `go_enriched_terms.csv` | Every enriched term for every arm x ontology x IEA variant, 1,515 rows, with the clusterProfiler columns plus `term_size`, `fold_enrichment`, `simplify_kept`, the block assignment, and the null columns. | Filter on `arm`, `ontology` and `iea_variant` before reading anything. The null columns are populated for the primary combination (BP, with_iea) only and are NA elsewhere, because the null is computed for the primary ontology alone. `Count` is arm genes in the term, `term_size` the term's size in the 13,128-symbol background. |
| `go_term_summary.csv` | One row per arm x ontology x IEA variant: arm size, annotated size, background size, terms testable, terms with at least one hit, terms enriched, terms after `simplify()`, blocks, and the null summary for both headline counts. | `n_enriched` against `null_median_n_signif` and `null_q95_n_signif` is the term-count read, `n_genes_covered` against `null_median_covered` and `null_q95_covered` the coverage read, with `p_n_signif_matched` and `p_covered_matched` the permutation p for each. `uniform_median_*` sizes the annotation-depth confound and supports no claim on its own. |
| `go_depth_null.csv` | Per-term null detail for the enriched terms of the primary ontology: observed hit count, recomputed p and q, the matched-null mean and 95th percentile of the hit count, `p_matched`, recurrence, and the same three against a uniform null. | Sorted by adjusted p within arm. `k_obs` and `p_obs_recomputed` are the matrix recomputation, and they agree with `go_enriched_terms.csv` to zero by the seam check in `go_provenance.csv`. |
| `go_null_design.csv` | Two row kinds. `band` rows give the annotation-depth bands: edges, how many query genes fall in each, how big the pool is, and whether the pool was short. `null_summary` rows give both nulls' summary statistics per arm. | A TRUE in `pool_short_of_query` means that band could not supply enough genes and the shortfall was borrowed from neighbouring bands; `mean_genes_borrowed` and `max_genes_borrowed` on the summary row say how much. Borrowing keeps every replicate exactly the query's size, which matters because both statistics scale with n. |
| `go_null_draws.csv` | Every one of the 2000 replicates of both nulls for all four arms: terms reaching significance and genes covered, per replicate, 16,000 rows. | The substrate behind `arms_observed_vs_null.png`. `n_signif_obs` and `n_covered_obs` repeat the arm's own values on every row so a permutation p can be recomputed from this file alone. |
| `go_term_clusters.csv` | Every enriched term's block assignment under average-linkage clustering of the Wang similarity matrix at height 0.7, with the block's size, representative and best adjusted p. | `cluster_representative` is the block's most significant term, chosen by adjusted p. Blocks are numbered per arm x ontology x IEA variant and the ids do not carry across those. |
| `go_iea_comparison.csv` | Per arm x ontology, what the IEA evidence filter changes: annotated query and background sizes, term counts each way, overlap, Jaccard, and both top-10 lists. | `n_top10_shared` is the number that matters: for `WT_heat_up` BP it is 4 of 10, so the headline list is not stable to the evidence filter and both variants are reported everywhere. |
| `go_vs_curated_lenses.csv` | For every arm x curated lens, the GO term whose name means the same thing, how many arm genes each of the two claims, the shared count, the Jaccard, and all three gene lists spelled out. Both IEA variants sit on the same row. | `go_term_status` says why a row has the result it has: `tested`, `above_max_gs_size`, `below_min_gs_size`, `no_arm_gene_in_term`, `term_absent_from_universe`, `ontology_not_run`, `absent_from_GO.db`. `go_term_enriched` is TRUE or FALSE only when `go_term_status` is `tested` and NA otherwise, so an untested term can never be read as a negative result. The `no_iea_*` columns repeat the same set of answers under the other evidence filter. |
| `go_gene_coverage.csv` | One row per arm x gene: whether the background holds it, how many BP and MF terms it carries, how many of the enriched ones hold it, its coverage class, and whether the other ontology can see it. | `coverage_class` is the three-way split drawn in `arms_coverage_ladder.png`. `visible_in_other_ontology` is why MF rides alongside BP: P4HA2 carries no BP term at all and a BP-only read would call it un-annotated. |
| `go_coverage_summary.csv` | That split aggregated per arm, with the coverage null beside it. | Read `n_in_enriched_term` against `null_median_covered` and `p_covered_matched`, never on its own. The null draws sets of the arm's ANNOTATED size (`n_arm_annotated`), while the three-way split is out of the nominal size (`n_genes`), which is why the two denominators differ. |
| `go_proteostasis_probe.csv` | Every configured proteostasis probe id for every arm, plus whatever the name sweep found, with the term's size in the universe, the arm genes it holds, the expected count, raw and adjusted p, `p_matched`, recurrence, and a status. | `status` and `tested` carry the same vocabulary as the lens table, and `enriched` is NA for any row that never entered the test. Rows appearing twice for one arm were found both by the configured id and by the name sweep. `genes_hit` is the whole point: read the gene list before the term name. |
| `go_mouse_replicate.csv` | The same over-representation asked in mouse space before the ortholog map, against the mouse universe, with the overlap against the human-space result. | Secondary by construction, since the arms are defined in human projection space. A term that appears only in human space is attributable to the projection rather than to the biology. |
| `go_provenance.csv` | Script, stage, run date, every package version, every cutoff, the null design, the seed, and the two seam checks. | `seam_enricher_vs_enrichGO_*` checks the hand-built annotation map against the canonical `enrichGO()` call. `seam_null_matrix_*` checks the sparse-matrix recomputation the null runs on against `enricher()`. Both are 0 here, which is what makes a null count and an observed count the same quantity. |
| `_overview/wtheatup_term_blocks.csv` | Per block for `WT_heat_up`: terms, genes, best adjusted p, terms passing `p_matched`, and median recurrence. | The figure's source table. `n_genes` is a union across the block's terms and the column does not sum to the arm. |
| `_overview/arms_observed_vs_null.csv` | Per arm x quantity: genes drawn per replicate, observed, null median, null 95th percentile, uniform null median, permutation p. | The eight numbers the panel prints, in one place. |
| `_overview/wtheatup_null_recurrence.csv` | All 384 enriched terms for `WT_heat_up` in rank order with size, hit count, adjusted p, `p_matched`, recurrence and the `simplify()` flag. | Sorted by adjusted p, so `rank` 1 is the strongest hypergeometric result. Scan `frac_matched_reaching_q` down the first twenty rows before quoting any of them. |
| `_overview/wtheatup_proteostasis_probe.csv` | The probe table for the primary arm with one row per distinct term and a `found_by` column. | Deduplicated against `go_proteostasis_probe.csv`: a term found both by id and by the name sweep appears once, with `found_by` naming both. |
| `_overview/arms_coverage_ladder.csv` | The three-way coverage split per arm with the coverage null and its permutation p. | Same columns as `go_coverage_summary.csv`, kept as the figure's same-stem neighbour. |

