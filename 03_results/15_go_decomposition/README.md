# 15_go_decomposition — Which biology the arms organise into, asked of the whole ontology

The curated decompositions ask how much of an arm falls inside a hand-picked panel of programs.
A panel that already contains hypoxia can answer whether hypoxia is in the arm, and two thirds
of the arm falls in no panel member. This stage asks the ontology instead: over the whole GO
Biological Process space, and separately the whole Molecular Function space, which terms does
each arm's gene content organise into? The candidate space is fixed by GO rather than by whoever
is asking.

**This describes the biology the arm's gene content organises into.** It is a different
measurement from composition, and its own count-level null is why. Over-representation gives a
gene as many memberships as it has terms: the 384 enriched GO BP terms for `WT_heat_up` hand out
3,597 gene-term memberships over a union of 160 genes, and 10 of the 202 genes sit in exactly one
enriched term. Collapsed to 36 similarity blocks it is 974 memberships with 22 genes in one
block. Shares taken off that do not add up to the arm. The partition question is asked off frozen
MSigDB collections in [`../16_arm_composition/`](../16_arm_composition/).

## Read the null before anything else

The gate that produced these arms thresholds effect size, which preferentially admits
well-studied genes, and the hypergeometric conditions on set size and universe alone. Every
count here is therefore compared against 2,000 random gene sets matched to the arm band by band
on GO BP terms per gene (`go_term_summary.csv`).

For `WT_heat_up` the 384 enriched terms sit against a matched-null median of 230, 95th percentile
370, maximum 603, permutation p 0.037. The 160 genes those terms cover sit against a null median
of 147, 95th percentile 159, maximum 167, permutation p **0.040**. A uniform null of the same
size reaches a median of 0 terms.

That coverage p is the number to hold onto. "GO BP places 160 of the 202 genes in an enriched
term" is the headline this directory would otherwise carry, and it describes what a random set of
equally well-studied genes achieves — the null's own maximum, 167, sits above the observed 160.
The other arms behave differently: `KO_heat_up` covers 166 against a null median of 129 (p 0.001),
and both `Interaction` arms clear their nulls comfortably because a 7-gene or 19-gene draw usually
reaches nothing at all. The arm this project cares about is the weakest of the four on both counts.

## Recurrence, and where in the list to look

For each enriched term, `wtheatup_null_recurrence.png` gives the share of matched draws that also
take it to significance. The median term recurs in 10.0% of draws, 126 of 384 terms recur above
20% and 24 above 50%. Recurrence is highest among the terms a reader would quote first:
chemotaxis 0.274, angiogenesis 0.367, blood vessel morphogenesis 0.458, positive regulation of
cell migration 0.574, positive regulation of cytokine production 0.612. A term at 0.6 is one the
ontology hands to three of every five well-annotated gene sets of this size.

## The proteostasis probe

Nineteen terms a heat-shock reading would predict were looked up by id, plus a name sweep over
the enriched list. Fifteen entered the hypergeometric and two cleared the adjusted-p cutoff:
`GO:0031649 heat generation` (IL1A, PTGS2, TNFSF11) and `GO:0001659 temperature homeostasis`
(EGR1, GATM, IGF2BP2, IL1A, KDM6B, NPR3, PTGS2, TNFSF11, p_matched 0.103). Both are
organism-level thermoregulation terms and inflammatory mediators carry them. Every chaperone,
folding and unfolded-protein term that entered came back above the cutoff, and the chaperone
genes present in the arm are HSPA1A and HSPH1. The four terms that never entered carry their
reason in the `status` column.

## Which term is tested, and which is only named

`go_vs_curated_lenses.csv` pairs each curated lens with the GO term whose name a reader would
take to mean the same thing. `GO:0006954 inflammatory response` holds 689 genes in the
`with_iea` universe, above the 500-gene cap, so it never entered the test while claiming 43 of
the arm's 202 genes — the largest claim of any term in the panel. Under `no_iea` it falls to 486
genes, enters, and comes back 5th of 225 for `WT_heat_up` and 2nd for `KO_heat_up`. Both variants
sit on the same row.

## Two evidence variants, and a mouse-space replicate

`with_iea` is primary because it is the canonical one-line `clusterProfiler` path, and the
study-depth confound that motivates dropping IEA is handled by the depth-matched null rather than
by an evidence filter. The two variants share 182 of the 384 `with_iea` BP terms (Jaccard 0.43)
and 4 of their top 10, so the evidence filter moves the headline list. The same question asked in
mouse space before the ortholog map returns 456 BP terms sharing 313 with the human-space 384
(Jaccard 0.59).

**Provenance.** GOSemSim 2.39.2 from the project-local R library, GO.db 3.22.0, org.Hs.eg.db
3.22.0, clusterProfiler 4.18.4, seed 20260730. Two seams are checked and recorded in
`go_provenance.csv`: the hand-built annotation map reproduces `enrichGO()` exactly on the primary
arm (same term ids, max |dp| 0), and the sparse-matrix recomputation the null runs on reproduces
`enricher()` term by term (max |dk| 0, max |dq| 0).

---

## Figures

### `figures/_overview/wtheatup_term_blocks.png`

**384 enriched terms, collapsed to 36 semantic blocks.**
Two panels sharing a row axis of blocks. Left: bar length gives the number of enriched terms in
the block, fill gives the smallest adjusted p among them, and the bar-end text gives the arm
genes those terms cover between them out of the arm's 202. Right: the block's depth-matched
verdict — the share of its terms whose observed hit count a matched random draw equals or beats
in under 5% of 2,000 replicates. A long bar beside a low dot is a large block a random draw
reproduces.

Blocks are an average-linkage cut of the Wang similarity matrix at height 0.7. The largest is
blood vessel morphogenesis: 56 terms over 81 arm genes, with 61% of its terms under
`p_matched` 0.05. Across all blocks the union of covered genes is 160. The shaded row spanning
both panels is the response-to-hypoxia block — 3 terms over 12 arm genes (ADAM8, AK4, EGR1,
FOSL2, HK2, HSPG2, INHBA, LTA, NPPC, NR4A2, PAK1, PTGS2), best adjusted p 0.0118, and 0% of its
terms under the depth-matched cut — shaded so a reader can find it among the 36.

A gene can sit in several blocks, so the gene counts do not sum to the arm.
*Source* `tables/_overview/wtheatup_term_blocks.csv` ·
`02_analysis/scripts/25_go_decomposition_viz.R`.

### `figures/_overview/wtheatup_null_recurrence.png`

**How often a depth-matched random draw finds the same term.**
Left panel: one dot per enriched term, at its rank by adjusted p (x) against how often a
depth-matched draw of the same size also takes it to significance over 2,000 replicates (y).
Blue dots mark terms a matched draw equals or beats in under 5% of replicates. The red curve is a
loess fit. Right panel: the same y axis as a histogram over all terms.

The dashed rule is the median recurrence, 0.100; the dotted rules are 20% and 50%, with 126
terms above the first and 24 above the second. The eight most significant terms are named, and
so are the three hypoxia and oxygen-level terms, in bold at ranks 161, 189 and 250 — all three
carry a `p_matched` at or above 0.05. A term high on this axis is one the ontology hands to any
well-annotated gene set of this size.
*Source* `tables/_overview/wtheatup_null_recurrence.csv` ·
`02_analysis/scripts/25_go_decomposition_viz.R`.

### `figures/_overview/wtheatup_proteostasis_probe.png`

**The direct by-name search for proteostasis and heat-shock terms.**
One row per probe term, split by ontology. Left panel: the dot gives the raw hypergeometric p on
a −log10 axis, with a dashed rule at p 0.05. Grey marks a term tested and above the adjusted-p
cutoff, vermillion tested and below it, and a blue cross at zero a term the test never saw.
Right panel: the arm genes each tested term holds, so a term clearing the cutoff reads by its
gene content; for an untested term it gives the reason — above the 500-gene cap, below the
10-gene floor, or holding no arm gene.

Fifteen of the nineteen probe terms entered and two cleared the cutoff: heat generation (3
genes: IL1A, PTGS2, TNFSF11; `p_matched` 0.011) and temperature homeostasis (8 genes;
`p_matched` 0.103). Every reported term appears whether or not it returned anything, so an
absence carries a number.
*Source* `tables/_overview/wtheatup_proteostasis_probe.csv` ·
`02_analysis/scripts/25_go_decomposition_viz.R`.

### `figures/_overview/arms_coverage_ladder.png`

**How many of each arm's genes an enriched term holds, against the null.**
One stacked bar per arm, in gene counts. Green gives genes held by at least one term that
reached significance, orange genes the ontology can place that no such term holds, grey genes
with no annotation in this ontology. Counts print inside a segment where it is wide enough. The
black tick is the median coverage over 2,000 depth-matched draws of the same size, and the
bar-end text gives the arm total, its coverage, that median and the permutation p.

For `WT_heat_up`: 160 of 202 genes sit in at least one enriched term, 31 carry an annotation and
sit in none, 11 carry no BP annotation at all (7 of those are visible in MF). A depth-matched
draw of 191 genes covers a median of 147, so the permutation p is 0.04. **Read the tick before
the green segment.**
*Source* `tables/_overview/arms_coverage_ladder.csv` ·
`02_analysis/scripts/25_go_decomposition_viz.R`.

---

## Tables

Thirteen are written by `24_go_decomposition.R`; the four under `_overview/` whose stems match a
figure are written by `25_go_decomposition_viz.R`. Everything lives under `tables/_overview/`.

Two conventions run through all of them. `p_matched` is the one-sided permutation p from 2,000
depth-matched draws: the share of draws whose hit count for that term reaches the observed one,
so small means the arm beats matched-random. `frac_matched_reaching_q` is recurrence: the share
of the same draws that take the term to q below 0.05 on their own, so large means the term is
cheap to enrich. A term wants both small.

| File | What is in it | How to read it |
|---|---|---|
| `go_enriched_terms.csv` | Every enriched term for every arm × ontology × IEA variant, 1,515 rows, with the clusterProfiler columns plus `term_size`, `fold_enrichment`, `simplify_kept`, the block assignment and the null columns. | Filter on `arm`, `ontology` and `iea_variant` first. Null columns are populated for the primary combination (BP, `with_iea`) and are NA elsewhere. `Count` is arm genes in the term; `term_size` is the term's size in the 13,128-symbol background. |
| `go_term_summary.csv` | One row per arm × ontology × IEA variant: arm size, annotated size, background size, terms testable, terms with a hit, terms enriched, terms after `simplify()`, blocks, and the null summary for both headline counts. | `n_enriched` against `null_median_n_signif` and `null_q95_n_signif` is the term-count read; `n_genes_covered` against `null_median_covered` and `null_q95_covered` the coverage read, with `p_n_signif_matched` and `p_covered_matched` beside them. `uniform_median_*` sizes the annotation-depth confound. |
| `go_depth_null.csv` | Per-term null detail for the primary ontology: observed hit count, recomputed p and q, the matched-null mean and 95th percentile, `p_matched`, recurrence, and the same against a uniform null. | Sorted by adjusted p within arm. `k_obs` and `p_obs_recomputed` are the matrix recomputation and agree with `go_enriched_terms.csv` to zero by the seam check. |
| `go_null_design.csv` | Two row kinds. `band` rows give the annotation-depth bands — edges, query genes per band, pool size, and whether the pool was short. `null_summary` rows give both nulls' statistics per arm. | TRUE in `pool_short_of_query` means the band borrowed from neighbours; `mean_genes_borrowed` and `max_genes_borrowed` say how much. Borrowing keeps every replicate exactly the query's size, which matters because both statistics scale with n. |
| `go_null_draws.csv` | All 2,000 replicates of both nulls for all four arms, 16,000 rows: terms reaching significance and genes covered, per replicate. | The substrate behind every null number above. `n_signif_obs` and `n_covered_obs` repeat the arm's own values on every row, so a permutation p recomputes from this file alone. |
| `go_term_clusters.csv` | Every enriched term's block assignment, with the block's size, representative and best adjusted p. | `cluster_representative` is the block's most significant term. Block ids are numbered per arm × ontology × IEA variant and do not carry across. |
| `go_iea_comparison.csv` | Per arm × ontology, what the evidence filter changes: annotated sizes, term counts each way, overlap, Jaccard, and both top-10 lists. | `n_top10_shared` is the number that matters: 4 of 10 for `WT_heat_up` BP, which is why both variants are reported everywhere. |
| `go_vs_curated_lenses.csv` | Per arm × curated lens, the GO term whose name means the same thing, how many arm genes each claims, the shared count, the Jaccard and all three gene lists. Both IEA variants sit on one row. | `go_term_status` says why a row reads as it does: `tested`, `above_max_gs_size`, `below_min_gs_size`, `no_arm_gene_in_term`, `term_absent_from_universe`, `ontology_not_run`, `absent_from_GO.db`. `go_term_enriched` is TRUE or FALSE only where `go_term_status` is `tested`, so an untested term never reads as a negative result. |
| `go_gene_coverage.csv` | One row per arm × gene: whether the background holds it, how many BP and MF terms it carries, how many enriched terms hold it, its coverage class, and whether the other ontology can see it. | `coverage_class` is the three-way split the coverage ladder draws. `visible_in_other_ontology` is why MF rides alongside BP: `P4HA2` carries no BP term at all and a BP-only read would call it unannotated. |
| `go_coverage_summary.csv` | That split aggregated per arm, with the coverage null beside it. | Read `n_in_enriched_term` against `null_median_covered` and `p_covered_matched`. The null draws sets of the arm's **annotated** size while the split is out of the **nominal** size, which is why the two denominators differ. |
| `go_proteostasis_probe.csv` | Every configured probe id for every arm plus whatever the name sweep found, with the term's size in the universe, the arm genes it holds, the expected count, raw and adjusted p, `p_matched`, recurrence and a status. | `status` and `tested` carry the lens table's vocabulary; `enriched` is NA for any row that never entered. A row appearing twice for one arm was found both by id and by the name sweep. Read `genes_hit` before the term name. |
| `go_mouse_replicate.csv` | The same over-representation asked in mouse space before the ortholog map, with the overlap against the human-space result. | Secondary by construction, since the arms are defined in human projection space. A term appearing only in human space is attributable to the projection. |
| `go_provenance.csv` | Script, stage, run date, package versions, every cutoff, the null design, the seed and the two seam checks. | `seam_enricher_vs_enrichGO_*` checks the hand-built annotation map against `enrichGO()`; `seam_null_matrix_*` checks the sparse-matrix recomputation against `enricher()`. Both are 0, which makes a null count and an observed count the same quantity. |

The four `_overview/` figure sources — `wtheatup_term_blocks.csv`,
`wtheatup_null_recurrence.csv`, `wtheatup_proteostasis_probe.csv` and
`arms_coverage_ladder.csv` — carry exactly the rows their panel draws. In the block table
`n_genes` is a union across the block's terms and the column does not sum to the arm; in the
recurrence table `rank` 1 is the strongest hypergeometric result, and
`frac_matched_reaching_q` is the column to scan down the first twenty rows before quoting any of
them.
