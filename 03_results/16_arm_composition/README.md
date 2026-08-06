# 16_arm_composition — What each heat-contrast arm is made of

Two earlier decompositions of these arms could not carry a composition figure. A hand-picked
nine-lens panel left 132 of `WT_heat_up` unclaimed against the 199-gene arm as it then stood, and
that panel has not been re-run against the 202-gene arm the ortholog fix produced. A
whole-ontology GO run gave 972 module memberships over 159 genes, a 6.1-fold double count.

This directory holds the third attempt, and it is the one that answers the partition question.
Sets come from five frozen MSigDB collections — Hallmark, GO BP, GO MF, KEGG, Reactome. A set
enters selection only when a depth-matched random gene set reaches its adjusted p in under 5% of
2,000 replicates, and it is dropped again when that significance recurs in more than 25% of them.
Redundant sets are pruned on the arm genes they share rather than on general set similarity.

## Two accountings, both reported

`fractional` gives a gene held by k selected sets 1/k to each. `winner_take_all` gives the gene
whole to its single best set. Both sum to exactly 1.0 over the arm, so the two are directly
comparable, and where they disagree that disagreement is the result.

`HALLMARK_TNFA_SIGNALING_VIA_NFKB` takes 0.044 of `WT_heat_up` fractional and 0.136
winner-take-all, and under `winner_take_all` it ties exactly with `GOMF_CYTOKINE_ACTIVITY`, both
on 27 genes. The two accountings put different categories at the head of the composition, and the
winner-take-all head is a tie between two collections of very different size, so both readings
are drawn on every face and neither is called primary.

## Which set takes a shared gene, and by what rule

Under `winner_take_all` a gene held by several selected sets goes to one of them. The ordering is
`enriched`, then `p_matched`, then the raw hypergeometric p-value, then fold enrichment, and it
is recorded in `composition_provenance.csv`. The raw p is the comparable quantity across
collections; an adjusted p is computed inside each collection's own family, and those families
run from 50 sets in Hallmark to 4,449 in GO BP. `p_matched` reaches its permutation floor for 30
of the 40 sets selected for `WT_heat_up`, so the key after it settles most assignments, which is
why the rule is stated rather than left implicit.

## The anchored variant is anchored

`hypoxia_anchored` pins a configured list of hypoxia sets into the selection whether or not they
reached significance. It sits beside the `unpinned` variant on every figure. For `WT_heat_up`
two of the pins were tested and stayed above the bar, so they carry NA significance, print in
purple wherever they appear, and their share is a gene count.

## How much of `WT_heat_up` carries a hypoxia annotation

Recorded in `composition_hypoxia_sources.csv`. Three sets hold hypoxia-annotated arm genes: 18,
11 and 4. Those three counts are un-addable — their union is 26 genes, 13.1% of the 199, and
Hallmark and GO BP share only 3 of them at a Jaccard of 0.115. Only `HALLMARK_HYPOXIA` reached
significance. GO MF and KEGG carry no hypoxia set at all, and a fourth pinned GO BP set has 7
background genes against a 10-gene floor, so those three cells record an absence of opportunity.

## What the null frame contains

Recorded in `composition_null_summary.csv`. A gene set drawn at random and matched on annotation
depth reaches significance in a median of 59 GO BP sets already, and in a median of zero sets in
Reactome, GO MF and KEGG. Annotation depth does most of the work in GO BP, and the three smaller
collections separate an observed result from a depth artefact far more sharply. Read any term
count in this directory against that frame.

**Claim tier.** Everything here is descriptive composition of a gene set derived elsewhere, and
it is hypothesis-generating. A set reaching significance is not evidence that the process it is
named for is present. The hypoxia counts measure annotation overlap.

**Provenance.** Background 12,986 human symbols from `WT_heat_ranked.rnk`. Gene sets from
msigdbr 26.1.0 / MSigDB 2026.1.Hs, md5 `c6ec92f75aa6062f511309493822f04b`. Selection under
clusterProfiler 4.18.4, seed 20260731, run 2026-07-31.

---

## Figures

### `figures/_overview/arm_composition.png`

**The composition of each arm, at two grains and under both accountings.**
Panel A places every gene at collection resolution: six bands, two bars per arm, each bar summing
to 1.0. Panels B to D go to set resolution, one row per category with two marks — a circle for
`fractional`, a triangle for `winner_take_all` — joined by a grey segment whose length is the
size of the disagreement.

A mark on the zero line is a real zero under that accounting, which is why the set-level panels
use paired points: 11 of `WT_heat_up`'s 40 categories take a share of exactly zero under one of
the two. The pair of numbers on the right repeats each row as fractional / winner-take-all, and
fill gives the source collection. Rows order by the larger of the two shares.

Under `unpinned`, 74 of 202 `WT_heat_up` genes sit in no selected set. The largest named category
under `fractional` is `HALLMARK_TNFA_SIGNALING_VIA_NFKB` at 0.042, and the same set takes 0.134
under `winner_take_all`, where it ties exactly with `GOMF_CYTOKINE_ACTIVITY` — 27 genes each, two
collections of very different size — so that accounting names no single leading pathway.
`KO_heat_up` leaves 100 of 221 unclaimed and `Interaction_fdrOnly_up` 1 of 19.

Only the ten largest sets per arm are drawn by name; the roll-up row states how many further sets
it holds and carries their summed share, so drawn rows plus roll-up plus residue total 1.0. Every
rolled-up set is named in `composition_shares.csv`.
*Source* `tables/_overview/arm_composition.csv` ·
`02_analysis/scripts/27_arm_composition_viz.R`.

### `figures/_overview/arm_composition_variants.png`

**What pinning the hypoxia sets recovers.**
Every row is drawn twice: a hollow grey circle for `unpinned` and a filled purple circle for
`hypoxia_anchored`, joined by a grey segment. The monospaced text on the right of each row gives
that column's pair as `unpinned → anchored` on the first line and the set's status in that arm on
the second, printed purple wherever that status reads `not enriched, NA`. A purple row is a set
that entered because the configuration named it, was measured, and stayed above the bar; its
share is a gene count.

Panel A runs to 3% and panel B to 100%, because the pinned shares and the residue differ by more
than an order of magnitude. Panel B's `all other selected sets` row aggregates every category
anchoring leaves alone.

Pinning takes `WT_heat_up` from 41 to 43 categories fractional and holds it at 30
winner-take-all. Both additions carry NA significance, because
`GOBP_RESPONSE_TO_OXYGEN_LEVELS` and `GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS` were tested and
stayed above the bar. The residue moves by 0.000, from 0.366 to 0.366, which is the whole of what
anchoring recovers. `HALLMARK_HYPOXIA` gives fractional weight to the pins it shares genes with,
0.026 down to 0.022, while its winner-take-all share holds at 0.020.
`REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA` was tested in every arm, holds none of their genes, and
sits at zero under both variants.
*Source* `tables/_overview/arm_composition_variants.csv` ·
`02_analysis/scripts/27_arm_composition_viz.R`.

---

## Tables

Every table sits in `tables/_overview/`. The ten `composition_*.csv` are written by
`26_arm_composition.R`; the two named after a figure are that figure's same-stem sibling and are
written by `27_arm_composition_viz.R`. Several `composition_*.csv` carry readings this directory
draws in no panel, and each is cited by name where its numbers are quoted above.

| File | What is in it | How to read it |
|---|---|---|
| `composition_ora_terms.csv` | 600 rows, one per set reaching significance, over all four arms and five collections, with the permutation columns joined on. | `p_matched` is the depth-matched permutation p and `frac_matched_reaching_q` how often a random depth-matched draw calls the same set significant. A small `p.adjust` beside a recurrence above 0.25 marks a set the ontology finds for anything, and selection drops it. |
| `composition_selected.csv` | 1,204 rows: every candidate set for every arm and both variants, with the flag that dropped it. For `WT_heat_up` unpinned there are 268 candidates — 65 fail the permutation p, 5 the recurrence bound, 85 are pruned as redundant, 40 are selected. | Read the four boolean columns in order. `pruned_by` names the surviving set that claimed nearly the same arm genes with the Jaccard that triggered it; 77 of the 85 prunes are by a set from the same collection. `pinned` rows are exempt from pruning and from the per-collection cap, and are never used as a pruner. |
| `composition_gene_assignment.csv` | 2,255 rows, one per arm, variant, gene and selected set, plus one row per unclaimed gene. The substrate both accountings are computed from. | `weight_fractional` is `1/n_sets_for_gene` and `is_winner` marks the set that takes the gene under `winner_take_all`. Summing `weight_fractional` over an arm and variant gives the arm's gene count, and so does counting `is_winner`. |
| `composition_shares.csv` | 342 rows, the full bar substrate: arm × variant × accounting × category, with the gene count and the share. Every set behind a roll-up row is here. | `share` sums to exactly 1.0 within each arm, variant and accounting. `n_genes` differs between the two accountings for the same set, because winner-take-all keeps only the genes that chose it. |
| `composition_remainder.csv` | 8 rows, the residue of each arm under both variants, split into three mutually exclusive classes, plus the count sitting in a set that reached significance and lost selection. | The three annotation columns partition `n_unclaimed`. `n_unclaimed_in_an_enriched_set` is a subset of `n_unclaimed_annotated_testable` and sits outside the partition. `genes_no_annotation_at_all` names the genes the five collections cannot place. |
| `composition_hypoxia_sources.csv` | 28 rows, every configured hypoxia pin for every arm, with its status, its size in the background, and the arm genes it holds spelled out. | `status` separates `tested, enriched` from `tested, not enriched`, from a set outside the tested size window, from a collection carrying no hypoxia set. The last two record an absence of opportunity. |
| `composition_hypoxia_overlap.csv` | 6 rows, the pairwise overlap of the hypoxia sets holding arm genes, for the two arms that have any. | `n_shared` and `jaccard` are computed on arm genes alone. This is the table that forbids adding the per-set counts. |
| `composition_null_summary.csv` | 39 rows, per arm and collection, for the depth-matched and uniform nulls side by side, plus one row whose `null` reads `not_tested` for a cell the engine could not run. | Compare `n_observed_signif` with `median_n_signif` and `q95_n_signif`, and `covered_obs` with `median_covered` and `q95_covered`. `p_covered` floors at 0.0005 over 2,000 replicates. Only the depth-matched rows carry the comparison. Filter on `null` before aggregating: a `not_tested` row has empty measurements and its `not_tested_reason` states why. |
| `composition_go_concordance.csv` | 4 rows comparing this GO BP result with the propagated `org.Hs.eg.db` GO BP under [`../15_go_decomposition/`](../15_go_decomposition/). | Two different tests of related content, one using MSigDB C5 GO:BP and the other a propagated map with an explicit IEA switch. Label Jaccard runs 0.48 to 0.74 across the arms; read it as how far the mapping choice moves the answer. |
| `composition_permutation_floor.csv` | 8 rows, one per arm and variant: how many selected sets sit exactly on the `p_matched` permutation floor of 1/2001. | `frac_at_floor` is the share for which the permutation p carries no further resolution, 0.75 for `WT_heat_up` unpinned. Where it is high, the key after `p_matched` settles most winner-take-all assignments. Raising `n_null` is what lowers it. |
| `composition_provenance.csv` | Package versions, the gene-set object md5, every threshold, the seed, the tie-break rule, the permutation floor, and the seam check between the matrix recomputation and `clusterProfiler::enricher`. | `seam_max_abs_dk` and `seam_max_abs_dq` are both 0, so the permutation engine and the ORA call agree set by set. Any change to `genesets_md5` invalidates every share here. |
| `arm_composition.csv` | 36 rows, the composition figure as drawn: the ten named sets per arm, the roll-up with its set count, and the residue. | `n_sets_rolled_up` is filled on the roll-up row alone. The full per-set list behind it is `composition_shares.csv`. |
| `arm_composition_variants.csv` | 96 rows, every category of the three arms carrying a composition, under both variants, with the two deltas. | `delta_fractional` and `delta_winner_take_all` are anchored minus unpinned. A row with `pinned` TRUE and `status` `tested, not enriched` gained its share by configuration and carries no significance. |
