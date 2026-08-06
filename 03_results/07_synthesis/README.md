# 07_synthesis — The four evidence arms in one table

Four methods have now read the same seven contrasts: gene-set enrichment, PROGENy footprints,
CollecTRI factor activity and per-gene differential expression, with GATOM modules where the
metabolic run returned any. This stage joins them into one tidy evidence table and draws the
single panel that reads across all four. It computes no new statistics.

The panel reports one thing: a **cGAS-dependence asymmetry**. The interferon arm carries the
interaction — 18 significant method-features in that column — while the HIF and glycolysis arm
rises in both heat arms and holds 1 significant feature there. That is an asymmetry between two
arms measured at n = 5. The interaction is the one-degree-of-freedom term and the least-powered
comparison in the design, so a flat column bounds cGAS-dependence at this sample size.

The arm names are labels for how each row was grouped. Naming HIF1α or HIF2α as the effector
would take a separate test, and [`../04_tf/`](../04_tf/) is where that test is made.

---

## Figures

### `figures/_overview/two_arms_panel.png`

**The interferon arm carries the interaction; the HIF and glycolysis arm carries both heat arms
and goes flat there.**
Two stacked tracks: the interferon/ISG arm above, the HIF/glycolysis arm below. Rows within a
track are method-features — GSEA gene sets, PROGENy pathways, CollecTRI factors, differential-
expression marker genes, GATOM modules where present — ordered by their interaction score.
Columns are the four headline contrasts: wild-type heat, cGAS-knockout heat, interaction,
temperature main effect. Tile fill gives the signed score, orange toward the numerator condition
and blue away from it, clamped to ±3.2, with the score printed in the tile. A black ring marks
adjusted p < 0.05.

**Read down the interaction column.** The interferon track lights up positive and ringed. The
HIF and glycolysis track stays lit in both heat columns and goes unringed there.
*Source* `tables/_overview/two_arms_panel.csv` · `02_analysis/scripts/16_synthesis_viz.R`.

---

## Tables

### `tables/two_arms_summary.csv`

The full cross-arm evidence table, one row per (arm, method, feature, contrast). Assembled by
`16_synthesis.R` from the four master tables under [`../master/`](../master/).

`arm` takes `IFN_ISG`, `HIF_glycolysis` or `heatshock_context` (HSF1, which is carried as
context rather than as a cGAS arm), with `arm_label` and `arm_hypothesis` naming what each
grouping was assembled to test. `method` names the source arm — `GSEA:<db>`, `PROGENy`,
`TF:CollecTRI`, `DE:limma-trend`, `GATOM` — and `feature` the gene set, pathway, factor or gene.
`score` is that method's own signed statistic: a normalised enrichment score, a linear-model
activity score, a moderated t, or a GATOM pseudo-NES. `padj` drives `significant` at the
configured FDR. The cGAS-dependence test lives in the `Interaction` rows.

Scores from different methods share no scale. Compare within a `method`, and read the panel for
the pattern across them.

### `tables/_overview/two_arms_panel.csv`

The panel's same-stem source: one row per tile, carrying every column needed to audit it.
`arm`, `arm_track` (the facet label drawn), `method`, `feature`, `feature_kind`, `contrast`,
`score` (raw and unclamped), `pvalue`, `padj`, `direction` and `significant`. Rows order by arm,
method, feature, contrast. The figure clamps `score` for display; this file carries the value the
method returned.
