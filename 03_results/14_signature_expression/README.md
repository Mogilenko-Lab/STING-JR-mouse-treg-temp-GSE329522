# 14_signature_expression: artifact captions

## figures/_overview/signature_dotplot.png

Across every member, not only the plotted ones, WT_heat_up and
KO_heat_up genes have a 39 °C cell as their high cell and a 37 °C
cell as their low cell in BOTH genotypes, whereas cGAS-KO 39 °C is
the LOW cell for all 9 Interaction_up and all 23
Interaction_up_fdrOnly genes.

**How to read:** One row per gene, one column per design cell (genotype x temperature,
n=5). Dot FILL = the gene's four cell means z-scored across those
cells: orange = its high cell, blue = its low cell, pale = flat; the
scale runs to ±1.5. Dot SIZE = mean log2 CPM. Colour and size are
independent channels. Rows run by the arm's own source-contrast t,
rank 1 at top; each strip gives the on-screen count. The interaction
panels' blue column is cGAS-KO 39 °C — a deficit in one cell, the
arithmetic content of an interaction term. Interaction_up and
Interaction_up_fdrOnly are one 1-df contrast at two gates, not two
signatures. Up arms only. Claim tier: L3 (n=5/group).
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/23_signature_expression_viz.R` | `ggplot(geom_point, shape=21)` | `signature_expression.top_n_genes=15; signature_expression.dot_size_range=[1.8, 6.5]; thresholds.de_fdr=0.05; thresholds.de_logfc=1; colors.diverging` | `03_results/14_signature_expression/tables/_overview/signature_dotplot.csv` |

## tables/_overview/signature_dotplot.csv

The 216 rows the dot plot draws: the top 15 members of each up arm by
descending source-contrast t (or all of them where an arm has fewer),
crossed with the four design cells.

**How to read:** Same-stem neighbour of figures/_overview/signature_dotplot.png and
deliberately self-contained: the figure joins nothing at draw time,
so a row here is a dot there. z_across_groups is the colour channel,
mean_logcpm the size channel, rank_within_signature the row order (1
= top). n_genes_in_signature is the arm's FULL member count, so
15-of-213 reads as a cap rather than a set size; z_bound carries the
four-cell z limit (1.5) the colour scale uses. Downstream notebooks
discover this table by its same-stem pairing with the PNG — keep the
pair intact. Claim tier: L3 (n=5/group). Sample-to-condition mapping
confirmed against the owner's sample sheet (2026-07-22): 20 of 20
libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/22_signature_expression.R` | `top-N rank cut + inner_join (rewritten byte-stably by save_overview)` | `signature_expression.top_n_genes=15` | `03_results/14_signature_expression/tables/signature_gene_stats.csv (+ signature_gene_group_expression.csv)` |

## tables/signature_gene_group_expression.csv

Every member of all four up arms in each of the four design cells:
group-mean log2 CPM with its dispersion and across-cell z.

**How to read:** Long format, one row per (signature, gene_symbol, group). group is
the genotype x temperature design cell, also split into genotype and
temperature; n_samples is 5 everywhere. mean_logcpm / sd_logcpm are
the mean and SAMPLE sd over those 5 libraries, se_logcpm is
sd/sqrt(n_samples). z_across_groups z-scores a gene's FOUR CELL MEANS
against each other, not its samples, so 'which cell is this gene's
high cell' reads the same for an abundant and a scarce gene; it is
bounded at ±1.5, and a gene flat across all four cells is reported as
0, not NA. One gene can appear under several signatures, so signature
is part of the key. Up arms only. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/22_signature_expression.R` | `group-mean/sd reshape + z_across()` | `signature_expression.expression_matrix=logcpm_mat; design.samples_per_group=5` | `03_results/objects/01_eda.rds, 03_results/objects/17_signature_sets.rds` |

## tables/signature_gene_stats.csv

One row per (signature, gene): the DE statistics the gene was
selected on, its rank inside its own arm, and what happens to it on
the way to human.

**How to read:** logFC, t and adj_P_Val are read verbatim from source_contrast's limma
table; nothing is re-fitted, and the two DE sources
(02_de_results.rds, master_de_genes.csv) are cross-checked against
each other before any row is written. rank_within_signature = 1 is
the arm's most extreme positive t, ties broken on symbol so the top-N
cut is deterministic. mapped_to_human is FALSE when the gene has no
ortholog in the frozen applied map, and is therefore dropped from the
human projection; human_symbol lists all orthologs, ';'-separated.
mapping_type is the gene's edge class over the WHOLE applied map, so
many2one means its ortholog is shared with another mouse gene
somewhere in the map, not inside this signature. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/22_signature_expression.R` | `de_of() + ortholog_fate left-join` | `gsea.rank_metric=t; thresholds.de_fdr=0.05; thresholds.de_logfc=1` | `03_results/objects/02_de_results.rds, 03_results/master/master_de_genes.csv, 03_results/human_projection/ortholog_map.tsv` |

## tables/signature_expression_summary.csv

One row per up arm: its size in mouse symbols, its ortholog fate, and
the typical effect size of its members — every value reproduces the
frozen projection ledger.

**How to read:** TWO senses of 'mapped' sit side by side because they differ and both
matter downstream. n_mapped_to_human counts MOUSE genes with at least
one human ortholog, so n_genes_mouse - n_dropped_no_ortholog =
n_mapped_to_human = n_one2one + n_one2many + n_many2one; it is the
sum of the per-gene mapped_to_human flag in signature_gene_stats.csv.
n_distinct_human_symbols counts the HUMAN symbols the arm lands on —
what 03_results/human_projection/manifest.csv calls n_human — and is
smaller whenever two mouse genes collapse onto one ortholog: for
WT_heat_up, 201 against 199. median_abs_t is 10.8/10.3 for the heat
arms against 6.7/5.8 for the interaction arm, the power gap between a
10-vs-10 contrast and a 1-df interaction term. n_genes_plotted
records the figure's cap. Tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/22_signature_expression.R` | `per-signature summary + manifest_row() cross-check` | `signature_expression.signatures; signature_expression.top_n_genes=15` | `03_results/objects/17_signature_sets.rds, 03_results/human_projection/{ortholog_map.tsv,manifest.csv}` |

