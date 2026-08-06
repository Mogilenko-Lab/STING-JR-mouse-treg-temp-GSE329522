# 06_gsea — The full gene-set sweep

Pre-ranked enrichment of every contrast's moderated-t ranking against thirteen gene-set
databases. Eight come from MSigDB (Hallmark, KEGG, Reactome, WikiPathways, GO_BP, GO_MF, GO_CC,
TF_Targets) and five are curated in this compartment (TransportDB, MitoPathways, MitoXplorer,
HSR_lens, TCR_activation). Seven contrasts × thirteen databases gives 91 cells, and all 91
emitted their four panels.

The enrichment is computed by `05_gsea_msigdb_run.R` and `06_gsea_custom_run.R`, which write
into `../master/master_gsea_table.csv`. Everything in this directory draws that table:
`12_gsea_viz.R` builds the per-cell battery and the asymmetry panel, `12b_gsea_overview_pooled_viz.R`
the three pooled overviews, `28_hypoxia_focus_viz.R` the two hypoxia panels. Running-sum curves
are recomputed deterministically at exponent 1 from the same ranked vector the sweep used, so a
curve and its tabulated NES agree by construction.

## What the sweep returns

Read off `../master/master_gsea_table.csv`, so each of these is checkable:

- **The two interferon-response sets are negative in both heat contrasts and positive on the
  interaction.** INTERFERON_ALPHA_RESPONSE reaches −1.90 (adjusted p 1.7e-04) in `WT_heat` and
  −2.52 (8.0e-12) in `KO_heat`, against +3.14 (2.7e-24) on the interaction. The gamma set runs
  −1.77, −2.20, +2.83 (5.0e-22). The interaction is positive because the decrease runs deeper in
  the knockout.
- **Hypoxia, glycolysis and mTORC1 are positive and significant in both heat contrasts.**
  HALLMARK_HYPOXIA reaches +1.91 (4.2e-06) in `WT_heat` and +1.95 (2.5e-06) in `KO_heat`, with an
  interaction NES of −1.27 at adjusted p 0.17. Oxidative phosphorylation stays negative
  throughout and reaches significance in no contrast.
- **A non-significant interaction NES reads as no detectable cGAS-dependence at n = 5.** The
  one-degree-of-freedom term is the least-powered contrast in this design.

A set enriching locates its gene content in a ranking. Establishing that the program the set is
named for is present is a separate measurement, made in
[`../15_go_decomposition/`](../15_go_decomposition/) and
[`../16_arm_composition/`](../16_arm_composition/).

## Conventions

**Sign.** NES > 0 places the set toward the contrast numerator (39 °C, or wild-type), NES < 0
toward the denominator. On the interaction, positive means the enrichment runs stronger in
wild-type.

**Glyphs.** Orange draws positive and blue negative. A star, a filled point or a black ring
marks adjusted p < 0.05.

**Parameters.** Ranking metric is the limma moderated t. Set-size window 15–500 genes. FDR is
Benjamini-Hochberg within a database, and `06b_gsea_pooled_padj.R` supplies the pooled
alternative where a cross-database comparison is wanted.

**Two ranking metrics govern every absence.** The dotplot and facet panels select their top 20
by adjusted p, and the barplot and running-sum panels select by |NES|. Worked example from the
master table: in the `WT_heat` × Hallmark cell, HALLMARK_HYPOXIA is 6th by |NES| (+1.908) and
4th by adjusted p (4.16e-06), so it falls outside the running sum's top five and sits inside the
dotplot's top twenty. HALLMARK_ANGIOGENESIS runs the other way, 2nd by NES and 10th by adjusted
p. An absence from one panel is a statement about that panel's ranking rule.

---

## The per-cell battery

`figures/by_contrast/<contrast>/<database>/` holds four panels, each as a `.pdf` and `.png`
pair — 728 files over 91 cells. The source table for a cell is
`tables/by_contrast/<contrast>/<database>/gsea_{msigdb,custom}.csv`.

### `dotplot.png`
**Every set in the cell, placed by leading-edge fraction.**
x, GeneRatio, the leading-edge genes divided by the effective set size; y, set name ordered by
GeneRatio descending. Point size gives −log10(adjusted p), fill gives NES clamped to ±3.2, and a
black outline marks adjusted p < 0.05. Selection is the top twenty by adjusted p.

### `facet.png`
**The same dotplot split by direction.**
Two panels, NES > 0 above and NES < 0 below, on the same encodings. Reading the two together
gives the direction balance of the cell at a glance.

### `barplot.png`
**Effect size for the significant sets alone.**
x, NES from zero; y, set name ordered by NES. Bars are drawn for FDR-significant sets only, top
twenty by |NES|.

### `running_sum.png`
**Where each set's genes sit along the ranking.**
Three stacked panels for the top five sets by |NES|. Top, the running enrichment score, which
steps up at each set member and decays between them, with the y range pinned to [−1, 1] so
curves compare across every cell of the sweep. Middle, one tick per member at its rank, one row
per set. Bottom, the ranked moderated t the curves were computed on.

### The contrasts and databases the wildcards range over

| Contrast | Design role |
|---|---|
| `WT_heat` | Heat response, wild-type arm. |
| `KO_heat` | Heat response, cGAS-knockout arm. |
| `Interaction` | The cGAS-dependence test; 1 df, lowest power. |
| `Geno_at_39` · `Geno_at_37` | Genotype simple effects, under heat and at baseline. |
| `Temp_main` | Pooled heat response, genotype collapsed. |
| `Geno_main` | Pooled genotype effect, temperature collapsed. |

| Database | Source | Content |
|---|---|---|
| `Hallmark` | MSigDB H | 50 curated hallmark sets. |
| `KEGG` · `Reactome` · `WikiPathways` | MSigDB C2:CP | Pathway collections. |
| `GO_BP` · `GO_MF` · `GO_CC` | MSigDB C5 | Gene Ontology, three ontologies. |
| `TF_Targets` | MSigDB C3:TFT | Transcription-factor target sets. |
| `TransportDB` · `MitoPathways` · `MitoXplorer` | Curated here | Transporter and mitochondrial metabolic sets, carrying no interferon set. |
| `HSR_lens` | Curated here | Two sets, `HSR_core` and `HSR_sensitivity`, built from msigdbr v2026.1.Hs. |
| `TCR_activation` | Curated here | One curated TCR / immediate-early activation set. |

---

## Cross-contrast figures

`figures/_overview/` holds six panels as `.pdf` / `.png` pairs.

### `gsea_asymmetry_panel.png`

**The two arms side by side, over the four headline contrasts.**
Two-block tile heatmap. Rows are the Hallmark sets whose name matches the block's pattern — top
block `INTERFERON|INNATE_IMMUNE|INFLAMMATORY|TNF|IL6|IL2|ALLOGRAFT` (7 sets), bottom block
`HYPOXIA|GLYCOLYSIS|MTORC1|OXIDATIVE_PHOSPHORYLATION` (4 sets). Columns are `WT_heat`,
`KO_heat`, `Interaction`, `Temp_main`. Fill gives NES squished to ±3.2, orange toward the
numerator and blue toward the denominator. A star marks adjusted p < 0.05. Rows within a block
order by NES. Each block draws its top twelve by max |NES|, and both blocks fall under that cap,
so every matching Hallmark set is on the panel.

The interaction column is the read: the interferon rows light up positive and ringed while the
hypoxia and metabolic rows stay lit in both heat columns and go flat there. The blocks group set
names by pattern. Which regulator drives either block is the question
[`../04_tf/`](../04_tf/) takes up.
*Source* `../master/master_gsea_table.csv` · `02_analysis/scripts/12_gsea_viz.R`.

### `gsea_pooled_overview_WT_heat.png` · `gsea_pooled_overview_Interaction.png` · `gsea_pooled_overview.png`

**What happens when one correction spans every database.**
Lollipop panel, one facet per database. x, NES squished to the configured cap; fill, orange
positive and blue negative. Glyph shape records what a set becomes when the per-database
Benjamini-Hochberg correction is replaced by one pooled correction across all 4,792 tests:
filled circle, significant under both; diamond, per-database only; triangle, pooled only; open
circle, neither. Each facet draws its top five by pooled adjusted p, ties broken on raw p then
|NES|, and the facet header reads (n/N) — how many of that database's sets pass the pooled FDR,
counted over the whole database.

`WT_heat` carries 860 significant sets per-database and 865 pooled (46 lost, 51 gained).
`Interaction` carries 192 and 164 (32 lost, 4 gained). A facet reading (0/23) marks a real null
and keeps its place, so the panel reports every database, including the ones that came back
empty. `gsea_pooled_overview.png` sets the two contrasts abreast on one canvas at double width,
sharing one legend. For reading a single contrast, prefer its own panel.

Pooled adjusted p is a comparability device across databases. It is a calibrated error rate
where the sets are independent, and GO terms and pathway sets share genes.
*Source* `tables/_overview/gsea_pooled_overview.csv` ·
`02_analysis/scripts/12b_gsea_overview_pooled_viz.R`.

### `hypoxia_routes_by_contrast.png`

**Five hypoxia-named sets over four databases, across three contrasts.**
Three stacked rows, one facet per contrast. x, each gene's position in that facet's own ranked
list as a fraction of its length, most numerator-shifted at 0. Top row, five overlaid running
curves keyed by colour **and** dash pattern so the panel holds under any form of colour
blindness and in mono, with y pinned to [−1, 1]. Middle row, one tick per member at its rank in a
labelled row per set. Bottom row, the ranked moderated t as one grey filled area per facet.
NES and adjusted p print inside each facet in the set's own colour.

Three of the five carry a positive NES and clear FDR 0.05 in both heat contrasts, and 7 of the
15 set × contrast cells clear it. `REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA` carries a negative NES
and clears none (adjusted p 0.72, 0.71, 0.94), so shared wording between two set names leaves
shared behaviour an open question.

**The two GO BP curves nest.** All 86 genes of `GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS` present
in the ranking sit inside the 172 of `GOBP_RESPONSE_TO_OXYGEN_LEVELS`. The child is the
cell-intrinsic part, and the parent also pools organism-level branches. Their agreement is partly
built in and their difference is the informative part. The parent clears FDR in both heat
contrasts (5.96e-04 and 9.91e-03) and the child clears it in wild-type alone (0.0234 against
0.119), so the cell-intrinsic half is the weaker half. Every interaction NES here stays above
threshold, which is the bound a 1-df term gives at n = 5.

The five sets were chosen by name, with their statistics playing no part in the choice. The
general panels of the sweep select on those statistics alone.
*Source* `tables/_overview/hypoxia_routes_by_contrast.csv` ·
`02_analysis/scripts/28_hypoxia_focus_viz.R`.

### `runsum_HALLMARK_HYPOXIA.png`

**One set, three contrasts, one axis.**
Three stacked panels sharing an x axis of fractional rank, most 39 °C-shifted at 0. Colour keys
the contrast and runs through all three panels. Top, the running enrichment score with y pinned
to [−1, 1]; middle, one tick row per contrast at each member's rank; bottom, the ranked
moderated t, each contrast filling the area between its curve and zero in transparent colour so
the three rankings stay legible where they coincide. The legend carries each contrast's NES,
adjusted p and genes present.

In `WT_heat` the set reaches NES +1.9240 at adjusted p 5.32e-06, with 172 of its genes in the
19,679-gene ranking, 57 in the leading edge, and a peak enrichment score of +0.458 in the top
12%. Across the three: `WT_heat` +1.92 (5.32e-06), `KO_heat` +1.96 (3.01e-06), `Interaction`
−1.27 (0.171). The set separates warmed from unwarmed iTregs in both genotypes.

This set is 6th by |NES| and 4th by adjusted p in its `WT_heat` Hallmark cell, so the sweep's
general running-sum panels leave it out and this panel draws it by name.
*Source* `tables/_overview/runsum_HALLMARK_HYPOXIA.csv` ·
`02_analysis/scripts/28_hypoxia_focus_viz.R`.

---

## Tables

| Path | What it holds |
|---|---|
| `tables/by_contrast/<contrast>/<database>/gsea_{msigdb,custom}.csv` | The per-cell result: one row per set with NES, p, adjusted p, effective size and leading edge. The source of that cell's four panels. |
| `tables/_overview/gsea_{msigdb,custom}_all.csv` | Every row of the MSigDB and custom sweeps, unfiltered. |
| `tables/_overview/gsea_msigdb_summary.csv` | Per-contrast counts of significant sets and the top set per database. |
| `tables/_overview/gsea_pooled_overview.csv` | The pooled-FDR family: every set with both corrections and the glyph class the overview panels draw. |
| `tables/_overview/gsea_pooled_summary_by_db.csv` | Per database and contrast, how many sets pass each correction. |
| `tables/_overview/hypoxia_routes_by_contrast.csv` · `runsum_HALLMARK_HYPOXIA.csv` | The two hypoxia panels' same-stem sources: NES, adjusted p and effective size per set and contrast, with the plotted curve points. |
| `tables/_overview/gsea_asymmetry_panel.csv` | The tile grid of the asymmetry panel. |

## One gene set removed from this stage, and why

`Lombardi2022_HIF` was scored here and has been withdrawn. It is a pan-cancer consensus derived
under hypoxia in cancer cells, which makes it the wrong reference for a 39 °C in-vitro contrast
in iTregs, and it meets `WT_heat_up` at 7 of that arm's 213 genes, 3.3%. Hypoxia in this stage
is read off the curated versioned collections, and
`figures/_overview/hypoxia_routes_by_contrast.png` is the panel that does it. The set survives
as a **membership lens over gene symbols** in [`../12_hsr_decomp/`](../12_hsr_decomp/) and
[`../13_semantic_decomp/`](../13_semantic_decomp/), where it is read for overlap and scored as a
gene-set database by neither.

The same episode produced a related caution. A bespoke sixteen-gene HIF list circulated as a
hypoxia reference, and its content is roughly 92% heat-shock and glycolytic while its
hypoxia-diagnostic core (Pdk1, Bnip3, Bnip3l, Car9) is repressed. Every gene-set size on a
figure in this stage is read at run time from the master table, so a size printed on a canvas is
the size the test used.
