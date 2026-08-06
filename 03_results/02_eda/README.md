# 02_eda — How much of the variance each design factor carries

Two panels place the twenty libraries in the space of their own variance. The question is a
sizing question: the 2×2 design has two factors, and this stage says how much of the
transcriptome each one moves.

The answer is asymmetric. PC1 carries 57.3% of the variance and tracks temperature almost
exactly (R² 0.98 against temperature, 0.00 against genotype). Genotype sits on PC2 at 4.4%,
within about a point of PC3. Both factors are recoverable from expression, and they are
recoverable to very different degrees.

The stage computes no statistics. It renders the decomposition behind the mapping check in
[`../01_qc/`](../01_qc/) and gives the variance budget every later contrast is read against.

---

## Figures

### `figures/fig1c_pca_2x2.png`

**Temperature is the leading axis of variation among the twenty libraries.**
Principal components of the delivered log2 CPM matrix. x, PC1 with its percent variance in the
axis title; y, PC2. Point colour gives temperature (blue, 37 °C; red, 39 °C) and point shape
gives genotype (circle, wild-type; triangle, cGAS-knockout); grey text labels each library.
PC1 separates the two temperatures cleanly and genotype separates along the much smaller PC2,
so the 2×2 reads on the plane with one factor dominant.

The gene universe is 19,657 symbols, every delivered symbol carrying variance across the twenty
libraries. The differential-expression figures report 19,679 because a principal-component
decomposition drops the 22 genes that are constant across all libraries; the two universes are
otherwise identical. Restricting to the 2,000 most variable genes gives the same leading axis
(|r| 0.99995) at 78.9%, so the percentage is a property of the gene universe.
*Source* `tables/fig1c_pca_data.csv` · `02_analysis/scripts/01_mapping_qc_viz.R`.

### `figures/fig1c_pca_scree.png`

**One component carries the experiment.**
Variance explained per principal component. x, component in descending variance order; y,
percent variance explained. Orange marks PC1, the temperature axis; grey marks PC2 onward.
PC1 at 57.3% stands clear of the tail while PC2 (4.4%) sits beside PC3 (3.5%), so the panel
reads as one dominant axis with a shoulder. Companion scatter: `fig1c_pca_2x2`.
*Source* `tables/fig1c_pca_varexp.csv` · `02_analysis/scripts/01_mapping_qc_viz.R`.

---

## Tables

### `tables/fig1c_pca_data.csv`

Per-library coordinates, twenty rows. `library_id` and its short `lib` suffix identify the
library; `PC1` and `PC2` are the scores the scatter draws; `temp`, `genotype` and `group` carry
the marker-derived condition and drive colour and shape. Principal-component signs are
arbitrary, so relative position is the readable quantity. Emitted by `01_mapping_qc.R` and
consumed read-only by the viz script.

### `tables/fig1c_pca_varexp.csv`

Variance explained per component, twenty rows in descending order. `PC` labels the component,
`pct_var` is the scree height, and `n_top` records the variable-gene count the decomposition
ran on. Written by `01_mapping_qc.R`; the viz script reads it and computes nothing.
