# Sample-mapping QC of the iTreg libraries

**Generated:** 2026-07-30 21:37 by `02_analysis/scripts/01_mapping_qc.R`

## Verdict (PROVISIONAL -- pending collaborator sample sheet)

The inferred temperature-major mapping (021-025 WT_37, 026-030 cGASKO_37, 031-035 WT_39, 036-040 cGASKO_39) is **SUPPORTED (gate PASS)**: **20/20** libraries have a data-derived label (thermometer-predicted temperature + Cgas-predicted genotype) that matches the inferred label. No discordant libraries.

**Evidence.** The heat-shock thermometer (Hspa1b/Hsph1/Hspa1a/Dnajb1) is monotone with the inferred temperature: mean log2CPM = 4.91 in the hot half (031-040) vs 3.79 in the cool half (021-030), a +1.12 log2 shift. Cgas (Cgas) is higher in WT than cGAS-KO within both temperature halves (37C: WT-KO = 1.02; 39C: WT-KO = 1.24 log2CPM). PCA on all 19657 genes carrying variance places 57.3%/4.4% of variance on PC1/PC2. PC1 is temperature almost exactly (R2 0.98 against temperature, 0.00 against genotype); genotype sits weakly on PC2 (R2 0.07), so the 2x2 is recoverable from expression but its two factors are recoverable to very different degrees.

**Scramble caveat.** The deposited CPM column order is temperature-major; the GEO GSM accessions are genotype-major. A naive positional GSM->column join therefore mislabels **10** libraries (026, 027, 028, 029, 030, 031, 032, 033, 034, 035) -- see `figures/fig1d_scramble.png`. This is why the mapping must be marker-derived, not accession-positional.

**Bottom line:** inferred mapping **SUPPORTED**, PROVISIONAL pending the collaborator sample sheet.

See `tables/mapping_verdict.csv` for the per-library machine-checkable verdict.


## figures/fig1a_thermometer.png

Heat-shock thermometer (label-blind temperature validation):
per-library log2(CPM+0.5) for the 4 heat-shock marker genes (Hspa1b,
Hsph1, Hspa1a, Dnajb1), one facet per gene, all 20 libraries
(12630-RS-021..040) on the x-axis, colored by INFERRED temperature.
Libraries 031-040 (inferred 39C, red) run systematically higher than
021-030 (inferred 37C, blue) in every marker, confirming the
deposited CPM column order is temperature-major. Label-blind: the
temperature is data-derived (marker expression), not read off the
deposited labels. Claim tier: PROVISIONAL (pending collaborator
sample sheet).

**How to read:** x = library (12630-RS-021..040, temperature-major order); y =
log2(CPM+0.5). Fill color = INFERRED temperature: blue = 37C
(negative pole of the diverging map), red = 39C (positive pole). One
facet per heat-shock marker (free y). Read: the hot half (031-040)
bars should tower over the cool half (021-030) in each facet. The
labels are inferred from marker expression (label-blind), NOT the
deposited sample labels. Claim tier: PROVISIONAL, pending the
collaborator sample sheet.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_mapping_qc_viz.R` | `geom_col + facet_wrap (project_theme + save_figure)` | `figures.width_column=3.5; figures.base_size_column=9; figures.base_size=16` | `03_results/01_qc/tables/fig1a_thermometer_data.csv` |

## figures/fig1b_cgas.png

Cgas genotype check (label-blind within each temperature half):
per-library log2(CPM+0.5) of Cgas, all 20 libraries split into the
inferred-37C and inferred-39C facets, colored by INFERRED genotype.
Within EACH temperature half, the WT libraries (green) exceed the
cGAS-KO libraries (purple), the expected signature of a Cgas knockout
-- confirming the genotype axis is data-recoverable. Label-blind:
genotype is inferred from Cgas expression, not the deposited labels.
Claim tier: PROVISIONAL (pending collaborator sample sheet).

**How to read:** x = library (12630-RS-, faceted into inferred-37C | inferred-39C
halves); y = Cgas log2(CPM+0.5). Fill color = INFERRED genotype:
green = WT, purple = cGASKO. Read WITHIN each facet (each temperature
half): WT bars should sit above cGAS-KO bars. The comparison is
intentionally within-temperature so the heat effect does not confound
the genotype call. Genotype is inferred from Cgas expression
(label-blind). Claim tier: PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_mapping_qc_viz.R` | `geom_col + facet_grid (project_theme + save_figure)` | `figures.width_column=3.5; figures.base_size_column=9; figures.base_size=16` | `03_results/01_qc/tables/fig1b_cgas_data.csv` |

## figures/fig1d_scramble.png

Scramble exhibit (mapping-competence demonstration): a 2-row tile per
library comparing the marker-derived condition assignment (top row)
against a naive positional GSM->column join (bottom row), all 20
libraries on the x-axis. 10/20 libraries (026-035) are DISCORDANT
(black outline + x marker): the deposited CPM columns are
temperature-major while the GEO GSM accessions are genotype-major, so
a positional join silently swaps KO-37 <-> WT-39 in the middle block.
This is why the mapping MUST be marker-derived, not
accession-positional. Claim tier: PROVISIONAL (pending collaborator
sample sheet).

**How to read:** x = library (12630-RS-021..040); y = the two label SOURCES (top =
marker-derived, bottom = naive positional GSM->column join). Tile
FILL = assigned condition (WT_37/cGASKO_37/WT_39/cGASKO_39, a
blue->red diverging-by-condition map). A black tile outline + an x
glyph on the strip = the two sources DISAGREE for that library (a
mislabeled column under the naive join). Read: the discordant block
(026-035) is exactly where the temperature-major vs genotype-major
orderings cross. Claim tier: PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_mapping_qc_viz.R` | `geom_tile + geom_point (project_theme + save_figure)` | `figures.width_column=3.5; figures.base_size_column=9; figures.base_size=16` | `03_results/01_qc/tables/fig1d_scramble_data.csv` |

