# 14_signature_expression — Each signature gene in the four design cells

The signatures are frozen. This stage puts their members back on the expression they were
selected from, one gene at a time across the four cells of the 2×2 design, so a reader can see
which cell each gene calls high and which it calls low. It also carries every member's DE
statistics and its fate on the way to human, so a signature's size, its power and its mapping
loss read from one place.

The pattern is clean and it separates the arms. Across every member — plotted and unplotted —
`WT_heat_up` and `KO_heat_up` genes have a 39 °C cell as their high cell and a 37 °C cell as
their low cell in both genotypes. For all 9 `Interaction_up` genes and all 23
`Interaction_up_fdrOnly` genes, the low cell is cGAS-knockout at 39 °C, which is the arithmetic
content of an interaction term: a deficit in one cell.

Up arms only, throughout.

---

## Figures

### `figures/_overview/signature_dotplot.png`

**Every plotted signature gene across the four design cells.**
One row per gene, one column per design cell (genotype × temperature, five libraries each).
Dot fill gives the gene's four cell means z-scored across those cells — orange marks its high
cell, blue its low cell, pale a flat gene — on a scale running to ±1.5. Dot size gives mean log2
CPM. Colour and size are independent channels, so an abundant flat gene draws large and pale.

Rows run by the arm's own source-contrast moderated t, rank 1 at the top, and each strip prints
the on-screen count. The cap is fifteen genes per arm, so a strip reading 15 of 213 is a cap
rather than a set size. The interaction panels' blue column is cGAS-knockout at 39 °C.
`Interaction_up` and `Interaction_up_fdrOnly` are one 1-df contrast at two gates.
*Source* `tables/_overview/signature_dotplot.csv` ·
`02_analysis/scripts/23_signature_expression_viz.R`.

---

## Tables

### `tables/_overview/signature_dotplot.csv`

The 216 rows the dot plot draws: the top fifteen members of each up arm by descending
source-contrast t, or all of them where an arm holds fewer, crossed with the four design cells.
Deliberately self-contained — the figure joins nothing at draw time, so a row here is a dot
there. `z_across_groups` is the colour channel, `mean_logcpm` the size channel,
`rank_within_signature` the row order. `n_genes_in_signature` carries the arm's full member
count so the cap reads as a cap, and `z_bound` carries the ±1.5 limit the colour scale uses.
Downstream notebooks discover this table by its same-stem pairing with the PNG; keep the pair
intact.

### `tables/signature_gene_group_expression.csv`

Every member of all four up arms in each of the four design cells. Long format, one row per
(signature, gene, group), with the group also split into its genotype and temperature factors;
`n_samples` is 5 everywhere. `mean_logcpm` and `sd_logcpm` are the mean and sample standard
deviation over those five libraries and `se_logcpm` is `sd/sqrt(n)`.

`z_across_groups` z-scores a gene's **four cell means** against each other rather than its
samples, so "which cell is this gene's high cell" reads the same for an abundant and a scarce
gene. It is bounded at ±1.5, and a gene flat across all four cells reports 0. One gene can appear
under several signatures, so `signature` is part of the key.

### `tables/signature_gene_stats.csv`

One row per (signature, gene): the DE statistics the gene was selected on, its rank inside its
own arm, and what happens to it on the way to human.

`logFC`, `t` and `adj_P_Val` are read verbatim from `source_contrast`'s limma table; nothing is
re-fitted, and the two DE sources (`../objects/02_de_results.rds` and
`../master/master_de_genes.csv`) are cross-checked against each other before any row is written.
`rank_within_signature = 1` is the arm's most extreme positive t, ties broken on symbol so the
top-N cut is deterministic.

`mapped_to_human` is FALSE where the gene has no ortholog in the frozen applied map, and
`human_symbol` lists all orthologs, semicolon-separated. `mapping_type` is the gene's edge class
over the **whole** applied map, so `many2one` means its ortholog is shared with another mouse
gene somewhere in the map rather than inside this signature.

### `tables/signature_expression_summary.csv`

One row per up arm: its size in mouse symbols, its ortholog fate, and the typical effect size of
its members. Every value reproduces the frozen projection ledger in
[`../11_projection/`](../11_projection/).

**Two senses of "mapped" sit side by side, and both matter downstream.**
`n_mapped_to_human` counts **mouse** genes with at least one human ortholog, so
`n_genes_mouse − n_dropped_unmapped_total = n_mapped_to_human = n_one2one + n_one2many +
n_many2one`; it is the sum of the per-gene `mapped_to_human` flag in `signature_gene_stats.csv`.
`n_distinct_human_symbols` counts the **human** symbols the arm lands on — what
`../human_projection/manifest.csv` calls `n_human` — and runs smaller wherever two mouse genes
collapse onto one ortholog.

The drop total then splits three ways. `n_dropped_no_ortholog` is the orthology source knowing
the symbol and having no human counterpart. `n_dropped_stale_query_symbol` is the source being
unable to key the symbol at all, because this matrix carries an older MGI vintage than the source
expects. `n_query_symbol_normalised` counts genes that arrived only because their symbol was
lifted to its current form before the source was asked.

`median_abs_t` reads 10.8 and 10.3 for the heat arms against 6.7 and 5.8 for the interaction
arms, which is the power gap between a 10-versus-10 contrast and a 1-df interaction term.
`n_genes_plotted` records the figure's cap.
