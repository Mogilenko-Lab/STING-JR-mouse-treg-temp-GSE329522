# What each mouse heat-contrast arm is made of

Two earlier decompositions of these arms could not carry a composition figure. A hand-picked
nine-lens panel left 132 of `WT_heat_up` unclaimed, against the 199-gene arm as it then stood — the
arm is 202 genes now, three of which arrived once the ortholog step stopped reporting symbols it
could not key as having no human ortholog, and that earlier panel has not been re-run against it. A
whole-ontology GO run gave 972 module-memberships over 159 genes, a 6.1-fold double count.

This directory holds the third attempt. Sets come from five frozen MSigDB collections
(Hallmark, GO BP, GO MF, KEGG, Reactome). A set enters only when a depth-matched random gene
set reaches its adjusted p in under 5% of 2,000 replicates, and it is dropped again when that
significance recurs in more than 25% of them. Redundant sets are pruned on the arm genes they
share, rather than on general set similarity.

**Two accountings, both reported.** `fractional` gives a gene held by k selected sets 1/k to
each. `winner_take_all` gives the gene whole to its single best set. Both sum to exactly 1.0
over the arm, so the two are directly comparable, and where they disagree that disagreement is
the result. `HALLMARK_TNFA_SIGNALING_VIA_NFKB` takes 0.044 of `WT_heat_up` fractional and 0.136
winner-take-all, and under `winner_take_all` it ties exactly with `GOMF_CYTOKINE_ACTIVITY`, both
on 27 genes. The two accountings put different categories at the head of the composition, and
the winner-take-all head is a tie between two collections of very different size, so neither
reading is called primary and both are drawn on every face.

**Which set takes a shared gene is a choice, and the rule is on the record.** Under
`winner_take_all` a gene held by several selected sets goes to one of them. The ordering is
`enriched, then p_matched, then raw hypergeometric pvalue, then fold_enrichment`, recorded in
`composition_provenance.csv`. The raw p is the comparable quantity across collections. An
adjusted p is not, because it is computed inside each collection's own family and those
families run from 50 sets in Hallmark to 4,449 in GO BP. `p_matched` reaches its permutation
floor for 30 of the 40 sets selected for `WT_heat_up`, so the key after it settles most
assignments, which is why the rule is stated rather than left implicit.

**The anchored variant is anchored.** `hypoxia_anchored` pins a configured list of hypoxia sets
into the selection whether or not they reached significance. It sits beside the `unpinned`
variant on every figure and never replaces it. For `WT_heat_up` two of the pins were tested and
did not reach the bar, so they carry NA significance, are printed in purple wherever they
appear, and their share is a gene count.

**How much of `WT_heat_up` carries a hypoxia annotation.** Recorded in
`composition_hypoxia_sources.csv`. Three sets hold hypoxia-annotated arm genes: 18, 11 and 4.
Those three counts cannot be added. Their union is 26 genes, 13.1% of
the 199, and Hallmark and GO BP share only 3 of them at a Jaccard of 0.115. Only
`HALLMARK_HYPOXIA` reached significance. GO MF and KEGG carry no hypoxia set at all, and a
fourth pinned GO BP set has 7 background genes against a 10-gene floor, so those three cells
record an absence of opportunity rather than a zero.

**What the null frame contains.** Recorded in `composition_null_summary.csv`. A gene set drawn
at random and matched on annotation depth
reaches significance in a median of 59 GO BP sets already, while reaching a median of zero in
Reactome, GO MF and KEGG. Annotation depth does most of the work in GO BP, and the three
smaller collections separate an observed result from a depth artefact far more sharply. Read
any term count in this directory against that frame.

**Claim tier.** Everything here is descriptive composition of a gene set derived elsewhere, and
it is hypothesis-generating. A set reaching significance is not evidence that the process it is
named for is present. The hypoxia counts measure annotation overlap and say nothing about
causal structure.

Background: 12,986 human symbols from `WT_heat_ranked.rnk`. Gene sets: msigdbr 26.1.0 / MSigDB
2026.1.Hs, md5 `c6ec92f75aa6062f511309493822f04b`. Selection ran under clusterProfiler 4.18.4
with seed 20260731 on 2026-07-31.

## figures/_overview/arm_composition.png

Under the `unpinned` variant 74 of 202 genes of WT_heat_up sit in no
selected set. The largest named category under `fractional` is
HALLMARK_TNFA_SIGNALING_VIA_NFKB at 0.042 of the arm, and the same
set takes 0.134 under `winner_take_all`. 2 categories tie exactly at
the head of `winner_take_all`, GOMF_CYTOKINE_ACTIVITY and
HALLMARK_TNFA_SIGNALING_VIA_NFKB, each on 27 genes and each taking
0.134 of the arm, and they come from collections of very different
size, so that accounting names no single leading pathway here. Which
pathway leads the composition therefore depends on how a gene held by
several sets is counted, which is why both readings are drawn and
neither is called primary. KO_heat_up leaves 100 of 221 unclaimed and
Interaction_fdrOnly_up 1 of 19.

**How to read:** Panel A is every gene at collection resolution: six bands, two bars
per arm, each bar summing to 1.0. Panels B to D go to set resolution.
Each row is one category with two marks, a circle for `fractional`
where a gene in k selected sets gives 1/k to each, and a triangle for
`winner_take_all` where the gene goes whole to its single best set.
The grey segment between them is the size of the disagreement and
both readings are reported. A mark on the zero line is a real zero
under that accounting, which is why the set-level panels use paired
points and not bars: 11 of WT_heat_up's 40 categories take a share of
exactly zero under one of the two. The pair of numbers on the right
repeats each row as fractional / winner-take-all, and fill is the
collection the set came from. Rows are ordered by the larger of the
two shares, so the residue and the roll-up sort in among the sets.
Only the ten largest sets per arm are drawn by name, and the roll-up
row states how many further sets it holds and carries their summed
share, so the drawn rows plus the roll-up plus the residue still
total 1.0. Every rolled-up set is named in composition_shares.csv.
Claim tier: descriptive composition of a derived gene set,
hypothesis-generating.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/27_arm_composition_viz.R` | `set_panel() / top-level (pA-pD)` | `arm_composition.collections = Hallmark,GO_BP,GO_MF,KEGG,Reactome, p_matched_cutoff = 0.05, max_null_recurrence = 0.25, prune_hit_jaccard = 0.5, top_n_per_collection = 10, n_null = 2000, seed = 20260731, winner_take_all_tiebreak = enriched, then p_matched, then raw hypergeometric pvalue, then fold_enrichment` | `03_results/objects/26_arm_composition.rds (md5 f06575ef9131a0061c39156a3065a0e6)` |

## figures/_overview/arm_composition_variants.png

Pinning the hypoxia sets takes WT_heat_up from 41 to 43 categories
fractional and 30 to 30 winner-take-all. Both additions carry NA
significance, because GOBP_RESPONSE_TO_OXYGEN_LEVELS and
GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS were tested and did not reach
the bar. The residue falls by 0.000, from 0.366 to 0.366, which is
the whole of what anchoring recovers. HALLMARK_HYPOXIA gives
fractional weight to the pins it shares genes with, 0.026 down to
0.022, while its winner-take-all share holds at 0.020.
REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA was tested in every arm and
holds none of their genes, so it sits at zero under both variants.

**How to read:** Every row is drawn twice: a hollow grey circle for the `unpinned`
variant and a filled purple circle for `hypoxia_anchored`, joined by
a grey segment. The monospaced text on the right of each row gives
that column's pair as `unpinned -> anchored` on the first line and
the set's status in that arm on the second, printed purple wherever
that status reads `not enriched, NA`. A purple row is a set that
entered because the configuration named it, was measured, and did not
reach the bar: its share is a gene count and supports no enrichment
claim. Panel A runs to 3% and panel B to 100%, because the pinned
shares and the residue differ by more than an order of magnitude and
one axis would flatten one of them. Panel B's `all other selected
sets` row aggregates every category anchoring leaves alone; the
per-category deltas are in the sibling table. Claim tier:
descriptive, and the purple rows support no enrichment claim at all.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/27_arm_composition_viz.R` | `variant_panel() / top-level (pV1/pV2)` | `arm_composition.collections = Hallmark,GO_BP,GO_MF,KEGG,Reactome, p_matched_cutoff = 0.05, max_null_recurrence = 0.25, prune_hit_jaccard = 0.5, top_n_per_collection = 10, n_null = 2000, seed = 20260731, winner_take_all_tiebreak = enriched, then p_matched, then raw hypergeometric pvalue, then fold_enrichment` | `03_results/objects/26_arm_composition.rds (md5 f06575ef9131a0061c39156a3065a0e6)` |


## tables/_overview/

Every table sits in `tables/_overview/`. The ten `composition_*.csv` are written by
`02_analysis/scripts/26_arm_composition.R`. The two named after a figure are that figure's
same-stem sibling and are written by `02_analysis/scripts/27_arm_composition_viz.R`. Several
`composition_*.csv` carry readings this directory no longer draws as a panel; each is cited by
name where its numbers are quoted above.

| File | What is in it | How to read it |
|---|---|---|
| `composition_ora_terms.csv` | 600 rows, one per set reaching significance, over all four arms and five collections, with the permutation columns joined on. | `p_matched` is the depth-matched permutation p and `frac_matched_reaching_q` is how often a random depth-matched draw calls the same set significant. A set with a small `p.adjust` and a `frac_matched_reaching_q` above 0.25 is a set the ontology finds for anything, and selection drops it. |
| `composition_selected.csv` | 1,204 rows: every candidate set for every arm and both variants, with the flag that dropped it. For `WT_heat_up` under the unpinned variant there are 268 candidates, 65 fail the permutation p, 5 fail the recurrence bound, 85 are pruned as redundant, and 40 are selected. | Read the four boolean columns in order. `pruned_by` names the surviving set that claimed nearly the same arm genes, with the Jaccard that triggered it, and 77 of the 85 prunes are by a set from the same collection. `pinned` rows are exempt from pruning and from the per-collection cap, and are never used as a pruner. |
| `composition_gene_assignment.csv` | 2,255 rows, one per arm, variant, gene and selected set, plus one row per unclaimed gene. This is the substrate both accountings are computed from. | `weight_fractional` is 1/`n_sets_for_gene` and `is_winner` marks the single set that takes the gene under `winner_take_all`. Summing `weight_fractional` over an arm and variant gives the arm's gene count, and so does counting `is_winner`. |
| `composition_shares.csv` | 342 rows, the full bar substrate: arm by variant by accounting by category, with the gene count and the share. Every set behind a roll-up row in the composition figure is here. | `share` sums to exactly 1.0 within each arm, variant and accounting. `n_genes` differs between the two accountings for the same set, because winner-take-all keeps only the genes that chose it. |
| `composition_remainder.csv` | 8 rows, the residue of each arm under both variants, split into three mutually exclusive classes, plus the count that sits in a set that reached significance and lost selection. | The three annotation columns partition `n_unclaimed`. `n_unclaimed_in_an_enriched_set` is a subset of `n_unclaimed_annotated_testable` and sits outside the partition. `genes_no_annotation_at_all` names the genes the five collections cannot place. |
| `composition_hypoxia_sources.csv` | 28 rows, every configured hypoxia pin for every arm, with its status, its size in the background, and the arm genes it holds spelled out. | `status` separates `tested, enriched` from `tested, not enriched`, from a set outside the tested size window, from a collection carrying no hypoxia set. The last two are an absence of opportunity and are not a zero share. |
| `composition_hypoxia_overlap.csv` | 6 rows, the pairwise overlap of the hypoxia sets that hold arm genes, for the two arms that have any. | `n_shared` and `jaccard` are computed on arm genes only, not on the whole sets. This is the table that forbids adding the per-set counts. |
| `composition_null_summary.csv` | 39 rows, per arm and collection, for the depth-matched and the uniform null side by side, plus one row whose `null` reads `not_tested` for a cell the engine could not run. | Compare `n_observed_signif` with `median_n_signif` and `q95_n_signif`, and `covered_obs` with `median_covered` and `q95_covered`. `p_covered` floors at 0.0005 over 2,000 replicates. Only the depth-matched rows carry the comparison, because the uniform draw is confounded by the gate admitting well-studied genes. Filter on `null` before aggregating: a `not_tested` row has empty measurements and its `not_tested_reason` states why, so the absence never has to be inferred. |
| `composition_go_concordance.csv` | 4 rows comparing this GO BP result with the propagated `org.Hs.eg.db` GO BP under `03_results/15_go_decomposition/`. | The two are different tests of related content, one using MSigDB C5 GO:BP and the other a propagated map with an explicit IEA switch. Label Jaccard runs 0.48 to 0.74 across the arms. Read it as how far the mapping choice moves the answer. |
| `composition_permutation_floor.csv` | 8 rows, one for each arm and variant: how many of its selected sets sit exactly on the `p_matched` permutation floor of 1/2001. | `frac_at_floor` is the share of selected sets for which the permutation p carries no further resolution, 0.75 for `WT_heat_up` under `unpinned`. Where it is high, the key after `p_matched` in the winner-take-all ordering is what settles most assignments, so read this table alongside `winner_take_all_tiebreak` in the provenance. Raising `n_null` is what would lower it. |
| `composition_provenance.csv` | Package versions, the gene-set object md5, every threshold, the seed, the winner-take-all tie-break rule, the permutation floor, and the seam check between the matrix recomputation and `clusterProfiler::enricher`. | `seam_max_abs_dk` and `seam_max_abs_dq` are both 0, so the permutation engine and the ORA call agree set by set. Any change to `genesets_md5` invalidates every share in this directory. |
| `arm_composition.csv` | 36 rows, the rows as the composition figure drew them: the ten named sets per arm, the roll-up with its set count, and the residue. | `n_sets_rolled_up` is filled only on the roll-up row. For the full per-set list behind that row, use `composition_shares.csv`. |
| `arm_composition_variants.csv` | 96 rows, every category of the three arms carrying a composition, under both variants, with the two deltas. | `delta_fractional` and `delta_winner_take_all` are anchored minus `unpinned`. A row with `pinned` TRUE and `status` `tested, not enriched` gained its share by configuration and carries no significance. |

