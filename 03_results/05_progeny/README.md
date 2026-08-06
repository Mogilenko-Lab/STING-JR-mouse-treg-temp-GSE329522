# 05_progeny — Signalling footprints across the seven contrasts

PROGENy infers the activity of fourteen signalling pathways from the expression of their
footprint genes, with continuous weights and no gene-set list. It fails differently from the
gene-set sweep in [`../06_gsea/`](../06_gsea/), which is why both are run on the same rankings.
CollecTRI transcription-factor activity is carried alongside on the same panels so the two
inference layers read together.

The result mirrors the two-arms asymmetry the differential expression shows. Hypoxia and
glycolysis rise in both heat arms and stay flat on the interaction, so those footprints carry
no detectable cGAS-dependence at n = 5. JAK-STAT, NF-κB and TNFα rise in `WT_heat` and stay
positive on the interaction, which is the cGAS-dependent arm.

**Compute** (`09_activity_progeny.R`, `03_decoupler_tf.R`) writes the activity scores into
`03_results/master/`. **Visualisation** (`13_activity_viz.R`) draws four cross-contrast panels
and, for the four headline contrasts, a per-contrast bar pair.

Both inference layers score activity with a linear model over target expression. A footprint
score describes how a pathway's target genes behave; pathway activity itself is untested.

---

## Figures

### `figures/_overview/progeny_heatmap.png`

**Fourteen pathways across all seven contrasts, on one tile grid.**
Rows, PROGENy pathways in hierarchical-clustering order; columns, contrasts in design order.
Fill gives the multivariate-model score, orange where the pathway is activated in the contrast
numerator and blue in the denominator, clamped to ±2.5. An asterisk marks raw p < 0.05; with
fourteen pathways no multiple-testing correction is applied. Hypoxia and glycolysis light both
heat columns and go flat in the interaction column; JAK-STAT, NF-κB and TNFα stay lit there.
*Source* `tables/_overview/progeny_heatmap.csv` · `02_analysis/scripts/13_activity_viz.R`.

### `figures/_overview/progeny_interaction_split.png`

**The hypoxia footprint and the immune footprints separate on the interaction.**
Grouped bars. x, the three headline contrasts (`WT_heat`, `KO_heat`, `Interaction`); y, activity
score; fill, pathway. Orange hypoxia bars stand at comparable height in the two heat contrasts
and near zero in the interaction. Blue JAK-STAT, NF-κB and TNFα bars run taller in `WT_heat`
than `KO_heat` and stay positive in the interaction. An asterisk marks raw p < 0.05. A flat
interaction bounds cGAS-dependence at n = 5 and leaves independence untested.
*Source* `tables/_overview/progeny_interaction_split.csv` · `02_analysis/scripts/13_activity_viz.R`.

### `figures/_overview/tf_heatmap.png`

**CollecTRI factor activity on the same grid.**
Rows, factors significant in at least two contrasts plus the watchlist, in clustering order;
columns, contrasts. Fill gives the univariate score clamped to ±2.5, orange toward the
numerator and blue toward the denominator; an asterisk marks BH-adjusted p < 0.05. The left-edge
strip gives the axis: orange for HIF (Hif1a, Epas1), blue for IFN/NF-κB, grey for other. The
IFN/IRF/STAT factors cluster as the block that is positive on the interaction; the HIF-axis rows
stay above threshold there.
*Source* `tables/_overview/tf_heatmap.csv` · `02_analysis/scripts/13_activity_viz.R`.

### `figures/_overview/progeny_tf_combined.png`

**Both inference layers on one canvas, for the three headline contrasts.**
Row facets separate PROGENy (top) from CollecTRI (bottom); column facets are the three headline
contrasts. Lollipop colour gives the axis, orange for HIF and blue for IFN/NF-κB, and shape
separates the layers, diamond for a pathway and circle for a factor. Rightward is activation in
the numerator condition. An asterisk marks p < 0.05, raw for PROGENy and BH-adjusted for the
factors. The Hypoxia diamond and the Hif1a circle both sit near zero in the interaction column,
while JAK-STAT, Stat1 and Irf3 stay positive there.
*Source* `tables/_overview/progeny_tf_combined.csv` · `02_analysis/scripts/13_activity_viz.R`.

### `figures/by_contrast/<contrast>/{progeny,tf}_barplot.png`

Eight panels over four contrasts — `WT_heat`, `KO_heat`, `Interaction`, `Temp_main` — one bar
plot per inference layer. x, activity score; y, pathway or factor, ordered by score. Orange runs
positive toward the numerator and blue negative. An asterisk marks p < 0.05, raw for PROGENy and
BH-adjusted for the factors. These are the per-contrast expansion of the cross-contrast panels
above and add no statistic.
*Source* `tables/by_contrast/<contrast>/{progeny,tf}_barplot.csv` ·
`02_analysis/scripts/13_activity_viz.R`.

---

## Tables

| File | What it holds | How to read it |
|---|---|---|
| `tables/progeny_activity.csv` | The stage's compute output: one row per pathway × contrast, 98 rows, with the multivariate score, its p-value and the contrast key. | Fourteen pathways × seven contrasts is complete by construction, so a missing row is an error. This table is also appended to `../master/master_gsea_table.csv` under `database = PROGENy`; count PROGENy once when pooling. |
| `tables/_overview/progeny_heatmap.csv` | The tile grid: pathway × contrast with the score, the clamped score and the significance flag. | `score` is unclamped; the figure draws it squished to ±2.5. |
| `tables/_overview/progeny_interaction_split.csv` | The three headline contrasts for the plotted pathway subset. | The hypoxia rows and the immune rows are the comparison; read the interaction column across them. |
| `tables/_overview/tf_heatmap.csv` | Factor × contrast scores with the axis label the left-edge strip draws. | `axis` takes `HIF`, `IFN`, or `other`, assigned from a fixed watchlist and never from the scores. |
| `tables/_overview/progeny_tf_combined.csv` | Both layers stacked with a `layer` column separating them. | Scores from the two layers share no scale; compare within a `layer`. |
| `tables/by_contrast/<contrast>/{progeny,tf}_barplot.csv` | Eight files, the same-stem source of each bar panel. | A per-contrast slice of the tables above, carried so each panel is reproducible from one file. |
