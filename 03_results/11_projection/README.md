# 11_projection — artifact captions

## figures/_overview/human_signature_sizes.png

Frozen human up/down set sizes per exported contrast after
mouse->human ortholog mapping.

**How to read:** Grouped bars per exported contrast; orange = up (higher in numerator
/ 39 C), blue = down; numbers = human genes in the frozen set. This
is the human-space counterpart of the stage-17 (mouse) size bars.
Claim tier: L3 (DE statistics); provisional sample mapping,
n=5/group.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_projection_export_viz.R` | `ggplot(geom_col)` | `decisions.projection.gate; colors.diverging` | `03_results/11_projection/tables/_overview/human_signature_sizes.csv` |


## figures/_overview/mapping_loss.png

Mapping did not decimate the primary signature: most mouse genes map
1:1; dropped/many-mapped fractions are small (alarm-flagged if > half
lost).

**How to read:** Per contrast (x) and direction (facet), the mouse up/down set is
stacked by fate: green = mapped 1:1, orange = one-mouse->many-human
(each human ortholog inherits the mouse t in the ranked list; unioned
in the binary set), grey = dropped (no human ortholog). Bar height =
mouse-set size; '->N' above = resulting human-set size. A red * marks
> 50% dropped — a decimated signature to demote. Read it as the
sanity check that the frozen human sets are not hollow. Claim tier:
L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/18_projection_export_viz.R` | `ggplot(geom_col)` | `alarm_frac=0.5; babelgene(offline); decisions.projection.ortholog_ambiguity` | `03_results/11_projection/tables/_overview/mapping_loss.csv` |

