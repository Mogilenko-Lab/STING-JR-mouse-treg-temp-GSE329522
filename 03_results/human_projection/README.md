# human_projection — The frozen signature contract

This directory is the interface between the mouse anchor and every human compartment. It holds
the mouse 39 °C signatures in human HGNC symbols, the ortholog map that produced them, and a
manifest recording every size and every loss. A human compartment scores against these files and
recomputes nothing.

It is a reshape and an ortholog map of the finished differential expression. The biology is in
[`../03_de/`](../03_de/); the mapping arithmetic is reported in
[`../11_projection/`](../11_projection/).

**Gate.** `fdr_logfc` — adjusted p < 0.05 and |log2FC| ≥ 1 — for every binary set, with the
`Interaction` contrast additionally exported at `fdr_only` as the sensitivity read a 1-df term
at n = 5 needs. Ranked lists carry the signed moderated **t**.

**Ortholog source.** babelgene 22.9, mouse to human, `min_support = 3`, so an edge is accepted
when at least three source databases agree. One mouse gene to many human takes the union for a
binary set and assigns the mouse t to each ortholog in a ranked list. Many mouse to one human
takes the union for a binary set and keeps the maximum |t| for a ranked list, except that an edge
from the original query outranks one recovered by symbol normalisation.

**Symbol normalisation.** babelgene keys on current MGI symbols and this matrix was quantified
against GENCODE vM25, so 2,341 of its 19,679 symbols are no longer current. Those babelgene could
not key are re-asked under their current symbol through `org.Mm.eg.db` 3.22.0 and the edges
mapped back; `mouse_symbol` always stays the symbol the data carries. 144 genes arrived that way,
each previously counted as having no human ortholog. The recovery is strictly additive — 155
edges added, none removed — and every candidate, accepted or withheld, is a row in
[`../11_projection/tables/_overview/query_normalisation_ledger.csv`](../11_projection/tables/_overview/).

---

## Files

### `manifest.csv`

The contract, one row per (contrast, direction, gate). `n_mouse` is the set size the gate
produced and `n_human` the size after mapping. The four drop columns split what happened in
between: `n_dropped_no_ortholog` (the source knew the symbol and returned no human counterpart),
`n_dropped_stale_query_symbol` (the source could not key the symbol at all),
`n_dropped_unmapped_total` (their sum, kept for comparison with earlier builds) and
`n_query_symbol_normalised` (genes that arrived only after normalisation). `n_many_mapped` counts
mouse paralogs collapsing onto one human symbol. `file` is the path relative to this directory.

The eleven rows at a glance:

| Contrast | Role | Gate | Up | Down | Ranked |
|---|---|---|---|---|---|
| `WT_heat` | primary | `fdr_logfc` | 213 → **202** | 126 → 96 | 19,679 → 13,128 |
| `KO_heat` | comparator | `fdr_logfc` | 239 → **221** | 153 → 115 | 19,679 → 13,128 |
| `Interaction` | comparator | `fdr_logfc` | 9 → **7** | 0 → 0 | 19,679 → 13,128 |
| `Interaction` | comparator | `fdr_only` | 23 → **19** | 0 → 0 | — |

The `Interaction` down arm is 0 genes at both gates, which is a structural absence and reads as
one everywhere downstream.

### `signatures/<contrast>/`

Three or four files per contrast. `<contrast>_up.txt` and `<contrast>_down.txt` are plain
newline-delimited HGNC symbols, one per line, sorted, with the direction carried by the filename.
`<contrast>_ranked.rnk` is the two-column signed ranked list — symbol and moderated t, descending
— over the full 13,128-symbol human universe. The `Interaction` directory additionally carries
`Interaction_fdrOnly_{up,down}.txt` at the relaxed gate.

### `ortholog_map.tsv`

Every applied edge: `mouse_symbol`, `human_symbol`, and the support and recovery flags behind it.
`mouse_symbol` is always the symbol the count matrix carries, so a lookup from this compartment's
own tables needs no translation.

### `SIGNATURES.md`

The full provenance record: the git SHA and build timestamp, the thresholds, the ortholog source
and its collision policy in prose, and the per-contrast export notes. Read it when auditing how a
file was made; read `manifest.csv` when checking what is in one.

---

## Two things a consumer should carry

**Every human denominator is the carried set.** `WT_heat_up` is 202 genes downstream, and the
213 belongs to mouse space. A caption quoting one alongside a human statistic should name which.

**One collapse lands on the heat-shock genes.** Mouse `Hspa1a` and `Hspa1b` both map to human
`HSPA1A`, so the three curated heat-shock-core genes in the `WT_heat_up` and `KO_heat_up` mouse
arms arrive as two in human, `HSPA1A` and `HSPH1`.
