# 02_eda — Exploratory Analysis

Label-blind PCA confirms that both experimental factors (temperature and genotype) are recoverable from expression data alone, supporting the inferred 2x2 design (WT/cGASKO x 37C/39C). This phase produces no statistics — it is a visual confirmation of the mapping QC verdict from phase 01_qc.

---

## figures/fig1c_pca_2x2.png

Temperature is the dominant axis of variation among the 20 libraries:
PC1 (57.3% var) separates 37C from 39C, and genotype tracks the much
smaller PC2 (4.4% var). Both experimental factors are recoverable
from expression alone, with the labels inferred rather than
deposited. Claim tier: PROVISIONAL (pending collaborator sample
sheet).

**How to read:** x = PC1 (% var in axis title); y = PC2. Point COLOR = inferred
temperature (blue = 37C, red = 39C); point SHAPE = inferred genotype
(circle = WT, triangle = cGASKO). Grey text = library id. Read the
2x2 as temperature separating along PC1 (the big axis) and genotype
along the much smaller PC2. GENE UNIVERSE: all 19,657 delivered
symbols carrying variance across the 20 libraries (22 constant genes
dropped), which is the same universe the DE stage models, so the %
var here is comparable to everything downstream. A PCA restricted to
the 2,000 most variable genes finds the SAME leading axis (|r|
0.99995 on PC1) but reports PC1 at 78.9%, because selecting on
variance concentrates variance into PC1; read the % var as a property
of the gene universe, not of the design. Labels are inferred
(label-blind), NOT deposited. Companion scree panel: fig1c_pca_scree.
Claim tier: PROVISIONAL.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/01_mapping_qc_viz.R` | `geom_point + geom_text_repel (project_theme + save_figure)` | `figures.width_column=3.5; figures.base_size_column=9; figures.base_size=16` | `03_results/02_eda/tables/fig1c_pca_data.csv` |

## figures/fig1c_pca_scree.png

One axis carries the experiment: PC1 captures 57.3% of the variance
and tracks temperature almost exactly, while PC2 (4.4%) sits within
about a point of PC3 (3.5%) and the rest of the tail. Temperature is
recoverable from a single component; genotype is not separated by the
leading axes at all. Claim tier: PROVISIONAL (pending collaborator
sample sheet).

**How to read:** x = principal component (descending variance); y = % variance
explained. Orange bar = PC1, the temperature axis; grey = PC2+. Only
PC1 stands clear of the tail, so read one dominant axis rather than
two design axes. GENE UNIVERSE: all 19,657 delivered symbols carrying
variance across the 20 libraries, matching the DE stage. On the 2,000
most variable genes alone PC1 would read 78.9% instead of 57.3%,
because variance selection concentrates variance into the leading
axis. Companion scatter: fig1c_pca_2x2. Claim tier: PROVISIONAL.

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

