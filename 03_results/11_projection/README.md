# 11_projection — Crossing the signatures into human symbols

The human compartments consume gene lists in HGNC symbols. This stage maps the frozen mouse
signatures across, and it records the fate of every gene, because that fate decides every
denominator downstream.

The map is babelgene 22.9, queried mouse to human at `min_support = 3`, so an edge is accepted
when at least three source databases agree. It holds 13,489 edges over 13,313 mouse and 13,128
human symbols. A binary gene set takes the union both ways: several mouse genes onto one human
symbol shrinks the set, one mouse gene across several human symbols grows it.

The three up arms arrive as: `WT_heat_up` 202 of 213, `KO_heat_up` 221 of 239, `Interaction_up`
7 of 9 (`tables/_overview/conversion_ledger.csv`). Every human denominator in this project is the
carried set rather than the mouse set.

**A gene lost at the seam is a vocabulary result until proven a biological one.** babelgene keys
on current MGI symbols and this matrix was quantified against GENCODE vM25, so a symbol whose
name has changed fails to key and used to be counted as having no human ortholog. Splitting that
bucket is what `tables/_overview/query_normalisation_ledger.csv` does, and it recovers 144 genes
across the modelled universe and three in `WT_heat_up`.

---

## Figures

### `figures/_overview/human_signature_sizes.png`

**Frozen human set sizes, per exported contrast and gate.**
Grouped bars. x, exported contrast; y, human genes in the frozen set, printed on the bar.
Orange gives the up set and blue the down set. A contrast carried at two gates gets one bar pair
per gate, with the gate named under its tick; the looser gate is the sensitivity read. This is
the human-space counterpart of `../10_signature/figures/_overview/signature_sizes.png`.
*Source* `tables/_overview/human_signature_sizes.csv` ·
`02_analysis/scripts/18_projection_export_viz.R`.

### `figures/_overview/mapping_loss.png`

**Most genes cross one-to-one, and the losses split four ways.**
Stacked bars. x, exported contrast; facets, direction; bar height, the mouse set size, with
`→N` at the bar end giving the human set size. Segments give the fate: green, mapped one-to-one;
orange, one mouse gene to several human, each ortholog inheriting the mouse t; grey, the
orthology source knew the symbol and had no human counterpart; purple, the source could not key
the symbol at all.

Grey and purple were one segment until this run, which reported a vocabulary loss as a
biological one; purple is the honest half of that bucket. A red asterisk marks a contrast losing
more than half its genes.
*Source* `tables/_overview/mapping_loss.csv` · `02_analysis/scripts/18_projection_export_viz.R`.

### `figures/_overview/conversion_ledger.png`

**Where each up arm's genes go on the way to human.**
One horizontal bar per frozen up arm. Bar length is the mouse genes that passed the `fdr_logfc`
gate, and the segments are their fate: blue carried into human, orange lost when several mouse
paralogs collapsed onto one human symbol, grey dropped for want of an accepted ortholog, purple
dropped because the source could not key the symbol. A count sits inside its segment where it
fits and is otherwise parked to the right on a leader line.

`WT_heat_up` carries 202 of 213 (8 no ortholog, 1 unkeyable, 2 collapsed, 3 recovered by
normalisation); `KO_heat_up` 221 of 239 (11, 5, 2, 3); `Interaction_up` 7 of 9 (1, 0, 1, 0).
One-to-one arrivals are 197 of 202, 216 of 221 and 3 of 7 in the same order.

One collapse lands on the heat-shock genes: mouse `Hspa1a` and `Hspa1b` both map to human
`HSPA1A`, so the three curated heat-shock-core genes in the `WT_heat_up` and `KO_heat_up` mouse
arms (`Hspa1a`, `Hspa1b`, `Hsph1`) arrive as two in human (`HSPA1A`, `HSPH1`).
*Source* `tables/_overview/conversion_ledger.csv` ·
`02_analysis/scripts/18_projection_export_viz.R`.

### `figures/_overview/projection_overlap_ledger.png`

**How much the projected heat arms share, and where the interaction arm sits.**
Facets separate up and down arms. The left stacked bar partitions the `WT_heat` / `KO_heat`
union into a wild-type-only slice, a shared slice and a knockout-only slice; the two right bars
give `Interaction` at each gate. Counts print on the bars and an open diamond marks a
structurally empty region.

The up-arm heat sets share 185 of 238 human genes (Jaccard 0.777), so their enrichment scores
computed from the same ranked list are dependent, and seeing both move together is close to
guaranteed. The down-arm heat sets share 81 of 130 (Jaccard 0.623). Both interaction gates stay
disjoint from the heat-set union. The wild-type-only up slice holds 17 human genes, of which 0
map from mouse genes flagged `cgas_dependent` in
`../03_de/tables/_overview/cgas_dependence_wide.csv`.
*Source* `tables/_overview/projection_overlap_ledger.csv` ·
`02_analysis/scripts/18_projection_export_viz.R`.

---

## Tables

### `tables/_overview/query_normalisation_ledger.csv`

The record of the symbol-normalisation step, one row per candidate. Of 295 matrix symbols the
orthology source could not key, 262 resolved one-to-one to a current MGI symbol and 144 of those
then mapped to a human ortholog — genes counted as having no human ortholog when the truth was
that their name had changed. Thirty-one pairs were withheld by a guard and two by decision.

`matrix_symbol` is the name the data carries, `current_symbol` what `org.Mm.eg.db` calls the
same Entrez id today, `resolution` what happened. `accepted` means the symbol resolved to exactly
one Entrez id, the current symbol is that id's own official symbol and no other gene's, no second
candidate claimed it, and it is not already a row here.

`rejected_current_symbol_already_in_vocabulary` is the guard that matters most quietly:
normalising there would give one gene two rows and let the many-mouse-to-one-human rule replace
an existing gene's ranked metric. `flagged_for_review` is withheld by decision rather than by a
guard — `Ndufb1-ps`, because handing a pseudogene-named row the real gene's ortholog judges the
vM25 annotation, and `H2aw`, because histone renaming is many-to-many.

`mapped_after_normalisation` splits the two outcomes an accepted pair can have. TRUE is a
recovery. FALSE means the symbol was re-asked correctly and still has no human ortholog at this
support level, which is an orthology result.

### The other four

| File | What it holds |
|---|---|
| `tables/_overview/human_signature_sizes.csv` | One row per exported contrast × direction × gate with the human set size. |
| `tables/_overview/mapping_loss.csv` | Per contrast, direction and gate, the four mapping-fate counts and the resulting human set size. |
| `tables/_overview/conversion_ledger.csv` | Per up arm: the mouse gate size, the four fates, the human size, the arrival routes (one-to-one, several-mouse-collapsed, one-mouse-split), and the paralog groups that collided inside the arm. |
| `tables/_overview/projection_overlap_ledger.csv` | The union partition per direction and the two interaction gates, with counts and Jaccard. |

The frozen lists themselves, the applied ortholog map and the manifest of sizes and gates are in
[`../human_projection/`](../human_projection/), which is the contract the human compartments
read.
