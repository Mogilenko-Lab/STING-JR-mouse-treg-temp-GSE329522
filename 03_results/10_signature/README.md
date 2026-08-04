# 10_signature: artifact captions

## figures/_overview/signature_sizes.png

Per-contrast up/down projection-set sizes in mouse symbols, at both
gates; the FDR+|logFC| gate is markedly more stringent than FDR-only.

**How to read:** Grouped bars per contrast; orange = up (higher in the numerator / 39
C), blue = down. TOP facet = FDR-only gate (adj.P.Val < 0.05), BOTTOM
= FDR + |log2FC| >= 1 (stringent). Numbers above bars = genes in that
set (mouse symbols, pre-ortholog). Use this to pick the gate before
the signatures are frozen: FDR+logFC can decimate small contrasts.
Claim tier: L3 (DE statistics), n=5/group. Sample-to-condition
mapping confirmed against the owner's sample sheet (2026-07-22): 20
of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_signature_derive_viz.R` | `ggplot(geom_col)` | `thresholds.de_fdr=0.05; thresholds.de_logfc=1; gsea.rank_metric=t; colors.diverging` | `03_results/10_signature/tables/_overview/signature_sizes.csv` |

## figures/_overview/updown_overlap.png

Off-diagonal Jaccard shows how distinct the candidate signatures
(WT_heat / Temp_main / Geno_at_39 / Interaction) are.

**How to read:** Symmetric heatmap; each tile = Jaccard overlap of two contrasts' UP
sets (left facet) or DOWN sets (right facet), gate = fdr_logfc.
Darker = more shared genes; diagonal = 1 (self). Off-diagonal near 0
means the signatures are distinct programs (good — comparators add
information). The source table carries BOTH gates. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_signature_derive_viz.R` | `ggplot(geom_tile)` | `gate=fdr_logfc; colors.diverging.down` | `03_results/10_signature/tables/_overview/updown_overlap.csv` |

## figures/_overview/ortholog_coverage_preview.png

Dry-run mouse->human coverage of each contrast's up+down set; a large
'unmapped' fraction flags mapping-loss risk before the freeze.

**How to read:** Stacked bars = fraction of a contrast's significant mouse genes
(up+down) that map 1:1 to a human ortholog (green),
one-mouse->many-human (orange), or have no human ortholog (grey,
dropped downstream). TOP facet = FDR-only, BOTTOM = FDR+logFC. This
is an offline babelgene PREVIEW for judging mapping loss. The frozen
human sets are re-mapped from the approved gene lists by
18_projection_export.R, and that run is what
03_results/11_projection/ reports. Claim tier: L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/17_signature_derive_viz.R` | `ggplot(geom_col, position=fill)` | `decisions.projection.ortholog_ambiguity.min_support; babelgene(offline)` | `03_results/10_signature/tables/_overview/ortholog_coverage_preview.csv` |

