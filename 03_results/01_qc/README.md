# 01_qc — Recovering the design from expression

The deposited CPM table arrives with twenty libraries and a proposed condition assignment. This
stage tests that assignment against the data, using two markers whose behaviour is known in
advance: the heat-shock genes report temperature, and `Cgas` reports genotype. The call is made
label-blind, from expression alone, and then compared with the assignment under test.

All twenty libraries agree. The heat-shock thermometer separates the two temperature halves by
+1.12 log2CPM, `Cgas` sits above the knockout in wild-type within both halves, and the
resulting per-library call matches the owner's sample sheet on 20 of 20 libraries
(`tables/mapping_verdict.csv`).

The stage exists because a naive join would get it wrong. The deposited CPM columns run
temperature-major while the GEO accessions run genotype-major, so joining by position mislabels
ten libraries. `figures/fig1d_scramble.png` is that demonstration, and it is why every
downstream label in this compartment is marker-derived.

**Compute** writes the four plot-ready tables and the machine-checkable verdict.
**Visualisation** draws three panels: one per marker family, plus the scramble exhibit.

---

## Figures

### `figures/fig1a_thermometer.png`

**The heat-shock thermometer recovers temperature from expression alone.**
Per-library log2(CPM+0.5) for four heat-shock genes, one facet per gene. x, library
(12630-RS-021 through 040, in deposited column order); y, log2(CPM+0.5), free per facet. Fill
gives the marker-derived temperature call: blue, 37 °C; red, 39 °C. Libraries 031–040 stand
above 021–030 in every facet, which places the temperature split at the midpoint of the
deposited column order. Mean log2CPM across the four markers is 4.91 in the hot half against
3.79 in the cool half.
*Source* `tables/fig1a_thermometer_data.csv` · `02_analysis/scripts/01_mapping_qc_viz.R`.

### `figures/fig1b_cgas.png`

**`Cgas` expression recovers genotype within each temperature half.**
Per-library `Cgas` log2(CPM+0.5), faceted into the 37 °C and 39 °C halves. x, library; y,
`Cgas` log2(CPM+0.5). Fill gives the marker-derived genotype call: green, wild-type; purple,
cGAS-knockout. Read within a facet: wild-type sits above knockout in both, by 1.02 log2CPM at
37 °C and 1.24 at 39 °C. Faceting by temperature keeps the heat effect out of the genotype
call.
*Source* `tables/fig1b_cgas_data.csv` · `02_analysis/scripts/01_mapping_qc_viz.R`.

### `figures/fig1d_scramble.png`

**A positional accession join mislabels half the libraries.**
Two rows of tiles per library, comparing the marker-derived assignment (top) against a naive
positional GSM-to-column join (bottom). x, library; y, the two label sources. Tile fill gives
the assigned condition on a four-level scale. A black outline with an × marks a library where
the two sources disagree; ten do, and they run consecutively from 026 to 035, exactly where the
temperature-major and genotype-major orderings cross. This compartment uses the marker-derived
mapping, and this panel is the reason.
*Source* `tables/fig1d_scramble_data.csv` · `02_analysis/scripts/01_mapping_qc_viz.R`.

---

## Tables

| File | What it holds | How to read it |
|---|---|---|
| `tables/mapping_verdict.csv` | One row per library: the assignment under test, the marker-derived temperature call, the marker-derived genotype call, and whether the two agree. | `concordant` reads TRUE on all 20 rows. This is the machine-checkable form of the verdict; every other artifact in this stage is its illustration. |
| `tables/fig1a_thermometer_data.csv` | Long form, one row per library × heat-shock marker, with `log2cpm` and the temperature call. | `temp` is derived from these values, so the column is the panel's conclusion and its colour channel at once. |
| `tables/fig1b_cgas_data.csv` | One row per library: `Cgas` log2(CPM+0.5), the temperature half it sits in, and the genotype call. | Compare within a `temp` group. The wild-type minus knockout difference is the quantity the panel reports. |
| `tables/fig1d_scramble_data.csv` | One row per library × label source, with the condition each source assigns and a `discordant` flag. | Filter on `discordant` to list the ten libraries a positional join gets wrong. |

Both marker families are ordinary quantities from the delivered CPM matrix, so this stage
computes no model and carries no test statistic. Its output is an identity check.
