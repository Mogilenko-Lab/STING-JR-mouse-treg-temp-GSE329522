# 10_signature — Choosing the gate, and sizing the candidate signatures

The differential expression is done; this stage turns it into named gene lists. Four contrasts
are carried as candidates — `WT_heat`, `Temp_main`, `Geno_at_39` and `Interaction` — each split
into an up set and a down set at two stringency gates: FDR alone (adjusted p < 0.05), and FDR
plus |log2FC| ≥ 1.

Three panels support the gate decision. The first sizes every candidate set at both gates. The
second asks whether the candidates are distinct programs or restatements of one another. The
third previews how much of each set survives the mouse-to-human ortholog step, so a signature
that would be decimated by mapping is visible before it is frozen.

Everything here is in mouse symbols. The frozen human sets are re-mapped from the approved lists
by `18_projection_export.R`, and [`../11_projection/`](../11_projection/) reports that run.

---

## Figures

### `figures/_overview/signature_sizes.png`

**Set size at both gates, per contrast and direction.**
Grouped bars. x, contrast; y, genes in the set, printed above each bar. Orange gives the up set
(higher in the numerator, 39 °C) and blue the down set. The top facet is the FDR-only gate and
the bottom facet is FDR plus |log2FC| ≥ 1. The stringent gate is markedly the smaller of the
two, and the difference is largest for the thin contrasts, which is the reason to look before
freezing.
*Source* `tables/_overview/signature_sizes.csv` · `02_analysis/scripts/17_signature_derive_viz.R`.

### `figures/_overview/updown_overlap.png`

**How much the candidate signatures share.**
Symmetric heatmap, two facets. Each tile gives the Jaccard overlap of two contrasts' up sets
(left facet) or down sets (right facet) at the `fdr_logfc` gate; darker means more shared genes
and the diagonal is 1 by construction. Off-diagonal values near zero mark distinct programs, so
carrying several candidates adds information rather than repeating it. The source table carries
both gates.
*Source* `tables/_overview/updown_overlap.csv` · `02_analysis/scripts/17_signature_derive_viz.R`.

### `figures/_overview/ortholog_coverage_preview.png`

**How much of each set would survive the crossing into human.**
Stacked bars, filled to 1. x, contrast; y, the fraction of that contrast's significant mouse
genes (up and down together) by mapping fate: green, one-to-one to a human ortholog; orange, one
mouse to several human; grey, no human ortholog and therefore dropped. The top facet is the
FDR-only gate and the bottom facet FDR plus |log2FC|.

This is an offline babelgene preview for judging mapping loss ahead of the freeze. The frozen
human sets come from the export run, and [`../11_projection/`](../11_projection/) reports that
one.
*Source* `tables/_overview/ortholog_coverage_preview.csv` ·
`02_analysis/scripts/17_signature_derive_viz.R`.

---

## Tables

| File | What it holds | How to read it |
|---|---|---|
| `tables/_overview/signature_sizes.csv` | One row per contrast × direction × gate, with the gene count. | `gate` takes `fdr_only` or `fdr_logfc`. The `WT_heat` × up × `fdr_logfc` row is 213, which is the arm that becomes `WT_heat_up`. |
| `tables/_overview/updown_overlap.csv` | Pairwise Jaccard between contrasts, per direction and per gate. | Carries both gates; the figure draws `fdr_logfc`. Read the off-diagonal. |
| `tables/_overview/ortholog_coverage_preview.csv` | Per contrast and gate, the counts behind the three mapping fates. | An offline preview at `min_support = 3`. The authoritative per-gene fate lives in [`../11_projection/`](../11_projection/). |

The frozen sets themselves are written by the export step into
[`../human_projection/signatures/`](../human_projection/), and the mouse-side set objects are
cached at `../objects/17_signature_sets.rds`.
