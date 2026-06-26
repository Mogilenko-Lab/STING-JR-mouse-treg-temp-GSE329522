# 02_eda — Exploratory Analysis

Label-blind PCA confirms that both experimental factors (temperature and genotype) are recoverable from expression data alone, supporting the inferred 2x2 design (WT/cGASKO x 37C/39C). This phase produces no statistics — it is a visual confirmation of the mapping QC verdict from phase 01_qc.

---

## figures/fig1c_pca_2x2.png

PCA on the top 2000 variable genes (label-blind): each library
plotted on PC1 (78.9% var) vs PC2 (3.2% var), colored by INFERRED
temperature and shaped by INFERRED genotype, with library ids
overlaid. The libraries split into the expected 2x2: temperature
dominates PC1 (37C vs 39C separate left/right) and genotype tracks
the secondary axis -- confirming both experimental factors are
recoverable from expression alone. Claim tier: PROVISIONAL (pending
collaborator sample sheet).

**How to read:** x = PC1 (% var in axis title); y = PC2. Point COLOR = inferred
temperature (blue = 37C, red = 39C); point SHAPE = inferred genotype
(circle = WT, triangle = cGASKO). Grey text = library id. Read:
temperature should separate along PC1 (the big axis), genotype along
PC2 -- the canonical 2x2 of a temperature x genotype design. Labels
are inferred (label-blind), NOT deposited. Companion scree panel:
fig1c_pca_scree. Claim tier: PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_mapping_qc_viz.R` | `geom_point + geom_text_repel (project_theme + save_figure)` | `figures.width_column=3.5; figures.base_size_column=9; figures.base_size=16` | `03_results/02_eda/tables/fig1c_pca_data.csv` |

## figures/fig1c_pca_scree.png

PCA scree / variance-explained panel (top 2000 variable genes): %
variance per principal component, PC1/PC2 highlighted. PC1 alone
captures 78.9% of the variance (the temperature axis) and PC2 a
further 3.2% (the secondary genotype axis); PC3+ contribute little,
confirming the 2x2 structure seen in fig1c_pca_2x2 is concentrated in
the leading two axes. Claim tier: PROVISIONAL (pending collaborator
sample sheet).

**How to read:** x = principal component (descending variance); y = % variance
explained. Orange bars = PC1/PC2 (the two axes that carry the
temperature x genotype 2x2); grey = PC3+. Read: a steep drop after
PC2 means the experimental structure lives in the leading axes (the
fig1c_pca_2x2 scatter is therefore a faithful 2D summary). Companion
scatter: fig1c_pca_2x2. Claim tier: PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_mapping_qc_viz.R` | `geom_col (project_theme + save_figure)` | `figures.width_column=3.5; figures.base_size_column=9; figures.base_size=16` | `03_results/02_eda/tables/fig1c_pca_varexp.csv` |

## tables/fig1c_pca_data.csv

Plot-ready tidy table of per-library PCA coordinates used to render `figures/fig1c_pca_2x2.png`. One row per library (n = 20); columns: `library_id` (full 12630-RS-0XX identifier), `lib` (short 3-digit suffix), `PC1`, `PC2` (principal component scores), `temp` (inferred temperature: "37" or "39"), `genotype` (inferred genotype: "WT" or "cGASKO"), `group` (combined condition: WT_37, cGASKO_37, WT_39, cGASKO_39). Temperature and genotype are marker-derived (label-blind), not deposited GEO labels. Emitted by `01_mapping_qc.R` (the compute half); consumed read-only by `01_mapping_qc_viz.R`.

**How to read:** Each row is one library. PC1/PC2 are the scores used for the x/y axes of the scatter plot. The `temp` and `genotype` columns drive color and shape encodings respectively. Sign of PC scores is arbitrary; only relative positions matter.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_mapping_qc.R` (compute) / `02_analysis/scripts/01_mapping_qc_viz.R` (consume) | `read.csv` (viz-only; no recomputation in viz script) | `—` | `03_results/objects/01_eda.rds` |

---

## tables/fig1c_pca_varexp.csv

Plot-ready tidy table of per-PC variance explained used to render `figures/fig1c_pca_scree.png`. One row per principal component (n = 20 rows); columns: `PC` (label, e.g. "PC1"), `pct_var` (percent variance explained), `n_top` (number of variable genes used for PCA; constant = 2000 across all rows). PC1 = 78.9%, PC2 = 3.2%, PC3 = 2.7%; PC20 ≈ 0 (rank deficiency). The `n_top` column also supplies the figure title text used in the viz script. Emitted by `01_mapping_qc.R`; consumed read-only by `01_mapping_qc_viz.R`.

**How to read:** Rows in descending variance order (PC1 first). The `pct_var` column is the y-axis of the scree plot. `n_top = 2000` for all rows reflects the variable-gene filter applied before PCA.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_mapping_qc.R` (compute) / `02_analysis/scripts/01_mapping_qc_viz.R` (consume) | `read.csv` (viz-only; no recomputation in viz script) | `—` | `03_results/objects/01_eda.rds` |

