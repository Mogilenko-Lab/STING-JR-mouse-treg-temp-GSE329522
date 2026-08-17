# interactive — Browser views over the same tables

Standalone HTML dashboards built from [`../master/`](../master/). They exist for exploration:
every number a dashboard shows is a number in a committed table, and each view names the table
it reads.

Three files: one bump-chart dashboard with its companion CSV, one reference embedding, and one
badge lookup.

---

## `gsea_bump_interaction.html`

**Two views of the same pathway-level heat response, on one shared filter pipeline.**
A radio selector at the top of the sidebar toggles the views; the database checkboxes, pattern
checkboxes, |NES| range sliders, interaction filter and keyword highlight are shared, and
switching views preserves every active filter.

**View 1 — Interaction, single panel.** x, the two genotype heat contrasts (`WT heat`,
`cGASKO heat`); y, enrichment score or ordinal rank, toggled in the sidebar. Each pathway is a
line joining its two endpoints, and the vertical gap between them is the visual heat-response
difference. Curvature encodes the interaction contrast — magnitude from |interaction NES|,
direction from its sign, drawn only where the interaction is significant. A flat line marks no
detectable cGAS-dependence at n = 5. Colour runs by pathway pattern or by any of the three NES
columns, on a blue-white-orange diverging scale capped at |NES| 3.5.

**View 2 — Heat trajectory, WT beside cGAS-KO.** Two panels sharing a locked y range, wild-type
left and knockout right. x, temperature, with ticks at 37 °C and 39 °C alone; y, enrichment score
labelled "NES (slope from 37 °C reference)". Each pathway is a straight line from (37 °C, 0) to
(39 °C, that genotype's heat-contrast NES), so a positive slope marks a pathway raised at 39 °C.
Comparing a pathway's slope between panels is the read. Steep in wild-type and flat in the
knockout, with a significant interaction, marks a cGAS-dependent heat response. Matched slopes
mark no detectable cGAS-dependence at n = 5.

**The 37 °C anchor sits at zero by construction.** This view plots one heat contrast per
genotype, so every line shares an origin whatever its true 37 °C baseline. It is a slope chart.
The rank toggle is hidden here because rank at the reference is undefined. An absolute
per-temperature trajectory would take per-group single-sample scores, which this dataset has no
run of.

*Built by* `02_analysis/scripts/bump_dashboard.py::DashboardPipeline.run()` *from*
`../master/master_gsea_table.csv`.

## `gsea_bump_interaction.csv`

The dashboard's source table: one flat wide row per pathway, 4,798 rows, carrying enrichment
score, adjusted p and rank for `WT heat`, `cGASKO heat` and the interaction, plus the computed
response-pattern category the colour modes use.

## `reference_embedding.parquet` · `leiden_identity_badges.json`

Supporting substrate: the stored gene-set layout coordinates, and the cluster identity labels
the badges draw.

---

## What a dashboard supports

Every panel here draws enrichment and activity statistics from the master tables. A pathway-level
interaction signal is a property of a ranked list and attributes to no single transcription
factor on its own — [`../04_tf/`](../04_tf/) is where that attribution is tested, and `fig3p`,
`fig3q` and `fig3r` are the three panels that bound it.

