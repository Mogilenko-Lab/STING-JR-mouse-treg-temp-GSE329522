# 11_projection: artifact captions

## figures/_overview/human_signature_sizes.png

Frozen human up/down set sizes per exported contrast and threshold
gate, after mouse->human ortholog mapping.

**How to read:** Grouped bars per exported contrast; orange = up (higher in numerator
/ 39 °C), blue = down; numbers = human genes in the frozen set. A
contrast carried at two threshold gates gets one bar pair per gate,
with the gate named in square brackets under its tick — the looser
gate is the sensitivity read, not a second result. This is the
human-space counterpart of the mouse-side size bars. Claim tier: L3
(DE statistics), n=5/group. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_projection_export_viz.R` | `ggplot(geom_col)` | `decisions.projection.gate; decisions.projection.secondary_gate; colors.diverging` | `03_results/11_projection/tables/_overview/human_signature_sizes.csv` |

## figures/_overview/mapping_loss.png

Mapping did not decimate the primary signature: most mouse genes map
1:1; dropped/many-mapped fractions are small (alarm-flagged if > half
lost).

**How to read:** Per exported contrast (x) and direction (facet), the mouse up/down
set is stacked by fate. Green = mapped 1:1. Orange = one mouse to
many human, each ortholog inheriting the mouse t. Grey = the
orthology source keyed the symbol and had no human counterpart.
Purple = it could not key the symbol at all.

Grey and purple were one segment until now, which reported a
vocabulary loss as a biological one. A contrast carried at two gates
gets one bar per gate, the gate named in brackets under its tick. Bar
height = mouse-set size; '->N' = human-set size. A red * marks > 50%
dropped, a decimated signature to demote. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_projection_export_viz.R` | `ggplot(geom_col)` | `alarm_frac=0.5; babelgene(offline); decisions.projection.ortholog_ambiguity; decisions.projection.secondary_gate` | `03_results/11_projection/tables/_overview/mapping_loss.csv` |

## figures/_overview/projection_overlap_ledger.png

The projected WT_heat and KO_heat up arms share 185 of 238 human
genes, while both Interaction gates remain disjoint from the heat-set
union.

**How to read:** Facets separate up and down arms. The left stacked bar partitions the
WT_heat/KO_heat union into WT_heat-only, shared, and KO_heat-only
slices; the two right bars show Interaction at the primary and
secondary gates. Counts are printed on the bars; open diamond means
structurally empty. The up-arm heat sets share 185 of 238 genes
(Jaccard 0.777), so their enrichment scores from the same ranked list
are not independent; seeing both move together is close to guaranteed
and is not corroboration. The down-arm heat sets share 81 of 130
genes (Jaccard 0.623). The WT-only up-arm slice has 17 projected
human genes; 0 map from mouse genes flagged cgas_dependent in
cgas_dependence_wide.csv. The separate HSR/TCR lens-membership
question is shown by `hsr_lens_membership`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_projection_export_viz.R` | `ggplot(geom_col + structural-empty glyph)` | `decisions.projection.gate; decisions.projection.secondary_gate; colors.diverging; colors.okabe_ito` | `03_results/11_projection/tables/_overview/projection_overlap_ledger.csv` |

## figures/_overview/conversion_ledger.png

Each frozen mouse up arm loses genes three ways: no human
counterpart, a symbol the orthology source could not key, and
paralogs collapsing onto one human symbol. WT_heat_up carries 202 of
213 (8 no ortholog, 1 not keyable, 2 collapsed, 3 recovered by
normalisation); KO_heat_up carries 221 of 239 (11 no ortholog, 5 not
keyable, 2 collapsed, 3 recovered by normalisation); Interaction_up
carries 7 of 9 (1 no ortholog, 0 not keyable, 1 collapsed, 0
recovered by normalisation).

One-to-one arrivals: 197 of 202, 216 of 221, 3 of 7 in the same
order. One collapse lands on the heat-shock genes: mouse Hspa1a and
Hspa1b both map to human HSPA1A, so the 3 curated HSR core genes in
the WT_heat_up and KO_heat_up mouse arms (Hspa1a, Hspa1b, Hsph1)
arrive as 2 in human (HSPA1A, HSPH1).

**How to read:** One horizontal bar per frozen up arm. Bar length is the mouse genes
that passed the fdr_logfc gate; the segments are their fate. Blue is
carried into human, orange is lost when several mouse paralogs
collapsed onto one human symbol, grey is dropped for want of an
accepted ortholog, and purple is dropped because the orthology source
could not key the symbol at all. Those last two used to share one
grey segment, which read every vocabulary loss as a biological one;
purple is the honest half of that bucket. A count sits inside its
segment where it fits, and is otherwise parked to the right of the
bar on a leader line back to the sliver.

The map is babelgene 22.9, queried mouse to human at min_support = 3,
so an edge is accepted when at least three source databases agree; it
holds 13,489 edges over 13,313 mouse and 13,128 human symbols. A
binary set takes the union both ways: several mouse onto one human
shrinks it, one mouse across several human grows it. Every human
denominator downstream is the carried set rather than the mouse set.

Per-arm arrival routes (one-to-one, several mouse collapsed onto it,
one mouse split across several human) and the paralog groups that
collided inside each arm are in the same-stem CSV. Claim tier: a
direct count over frozen sets.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_projection_export_viz.R` | `ggplot(geom_rect + leader-line sliver labels)` | `decisions.projection.gate=fdr_logfc; decisions.projection.contrasts_primary; decisions.projection.ortholog_ambiguity.min_support=3; babelgene=22.9; colors.okabe_ito` | `03_results/human_projection/{manifest.csv,ortholog_map.tsv,signatures/}; 03_results/objects/17_signature_sets.rds` |

## tables/_overview/query_normalisation_ledger.csv

Of 295 matrix symbols the orthology source could not key, 262
resolved one-to-one to a current MGI symbol and 144 of those then
mapped to a human ortholog — genes that had been counted as having no
human ortholog when the truth was that their name had changed. 31
pairs were withheld by a guard and 2 by human decision.

**How to read:** One row per candidate. `matrix_symbol` is the name the data carries,
`current_symbol` what org.Mm.eg.db calls the same Entrez id today,
`resolution` what happened.

`accepted` means the symbol resolved to exactly one Entrez id, the
current symbol is that id's own official symbol and no other gene's,
no second candidate claimed it, and it is not already a separate row
here. `flagged_for_review` is withheld by decision, not by a guard:
Ndufb1-ps because handing a pseudogene-named row the real gene's
ortholog judges the vM25 annotation, H2aw because histone renaming is
many-to-many.

`rejected_current_symbol_already_in_vocabulary` is the guard that
matters most quietly. Normalising there would give one gene two rows
and let the many-mouse-to-one-human rule silently replace an existing
gene's ranked metric.

`mapped_after_normalisation` splits the two outcomes an accepted pair
can have. TRUE is a recovery. FALSE means the symbol was re-asked
correctly and still has no human ortholog at this min_support — a
real orthology result, not a vocabulary loss. Claim tier: a direct
count over an annotation database, no statistics.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_projection_export_viz.R` | `helpers/ortholog_utils.R::normalise_mouse_query (written by 18_projection_export.R)` | `symbol_alias.ortholog_query_flagged_for_review; org.Mm.eg.db` | `03_results/objects/17_signature_sets.rds (the modelled symbol universe)` |

