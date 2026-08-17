# The mouse 39 °C iTreg anchor — GSE329522

Bulk RNA-seq of induced regulatory T cells cultured at 37 °C and at 39 °C, in wild-type and
cGAS-knockout genotypes. This repository holds the analysis of that experiment and the gene
signatures it produces.

The experiment is the anchor of a larger project asking whether human inflammatory and
autoinflammatory disease states contain transcriptional programs consistent with a mouse
temperature-stress axis in T cells and with cGAS–STING-related biology. This repository is the
mouse half. It derives the signatures; the human compartments consume them.

## The experiment

Induced regulatory T cells were differentiated from primary murine splenic CD4⁺ T cells and
cultured in a 2×2 design — genotype (wild-type, cGAS-knockout) × temperature (37 °C, 39 °C) —
with five biological replicates per group, twenty libraries in total. Deposited as
[GSE329522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329522).

## What the analysis found

Temperature dominates the design. PC1 carries 57.3 % of the variance and tracks temperature
almost exactly, and warming changes 8,723 genes in wild-type and 8,901 in cGAS-knockout cells at
FDR 0.05. The two genotypes warm in step: across 10,418 heat-responsive genes the knockout
response is 0.99× the wild-type response at r = 0.95.

The cGAS-dependence test is the heat-by-genotype interaction, and it is the least-powered contrast
in the design at one degree of freedom and n = 5. Twenty-three genes clear FDR 0.05 there. A gene
failing that test has **no detectable cGAS-dependence at n = 5**, which is the phrasing every
artifact in this repository uses — it is not evidence of cGAS-independence.

The interferon and hypoxia arms behave differently under that test, and that asymmetry is the
result. The two Hallmark interferon-response sets carry the interaction, at +3.14 and +2.83 with
adjusted p 2.7e-24 and 5.0e-22. HALLMARK_HYPOXIA rises in both heat arms, +1.91 and +1.95, and
stays flat in the interaction at −1.27, adjusted p 0.17.

Each number above is reproduced in a table under `03_results/`, named in that stage's `README.md`.

## The deliverable

`03_results/human_projection/` is the contract the human compartments read: three frozen up arms
in human symbols, each shipping an up list, a down list, a signed ranked list, and a manifest
recording sizes and gates.

| Signature | Mouse symbols | Human symbols |
|---|---|---|
| `WT_heat_up` | 213 | 202 |
| `KO_heat_up` | 239 | 221 |
| `Interaction_up` | 9 | 7 |

Genes are gated on FDR and \|log2FC\| ≥ 1, and converted to human symbols with pinned offline
babelgene 22.9. Every gene's fate in that conversion is recorded, so a symbol missing downstream
can be told apart from a gene that was never there.

## Layout

| Path | Contents |
|---|---|
| `00_data/` | Read-only inputs: the normalized CPM matrix with its provenance record, and the symbol-alias map. Raw sequence data and the external reference sets are not redistributed here — the scripts fetch or cache them. |
| `02_analysis/scripts/` | 60 numbered scripts, 59 in R and one in Python. Every artifact under `03_results/` is reproducible from one of them; a `_viz` script never computes and a compute script never plots. |
| `02_analysis/config/` | `analysis_config.yaml` holds every path and threshold. Nothing is hardcoded in a script. |
| `02_analysis/helpers/` | Shared functions: figure styling, ortholog conversion, symbol-alias resolution, source hashing. |
| `03_results/` | One directory per analysis stage, each with `tables/`, `figures/` and a `README.md` captioning every file. `03_results/README.md` is the reading order. |
| `.devcontainer/` | The container the analysis runs in. |

`03_results/README.md` is the best entry point: it gives the reading order stage by stage, the
headline numbers with the table each comes from, and the provenance of every gene set, lens and
network the analysis reads.

## Two conventions that hold throughout

**Signs.** A positive statistic points to the numerator of its contrast: 39 °C for a heat
contrast, wild-type for a genotype contrast. For the interaction, positive means the heat response
is larger in wild-type.

**Naming.** A signature is named for how it was derived, never for the mechanism it is hoped to
represent. `WT_heat_up` is the set of genes higher at 39 °C than at 37 °C in wild-type cells, and
that description is checkable against the differential-expression tables. A gene set enriching
locates its gene content in a ranking; establishing that the program the set is named for is
actually present is a separate measurement, made in the decomposition stages.

## Reproducing

The scripts run in numeric order inside the provided container, reading paths and thresholds from
`02_analysis/config/analysis_config.yaml`. Two things are worth knowing before starting:

- **Figures are not tracked.** `.png` and `.pdf` are ignored repo-wide. Each stage's tables and its
  `README.md` are the tracked record, and every figure regenerates from the script its caption
  names.
- **`01_modules/` is empty in a source archive.** The two toolkits under it are git submodules, and
  GitHub's release tarballs do not include submodule contents. Clone with `--recursive` to get
  them.

## Data availability

The count matrix analysed here is deposited at
[GSE329522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329522). Every external gene
set, lens and network is listed with its version and retrieval method in
`03_results/README.md`, so each is traceable to its source rather than to this repository.

## License

Two licenses, split by what the file is:

- **Code** — everything under `02_analysis/`, `.devcontainer/` and the scripts anywhere in the
  tree: [MIT](LICENSE).
- **Results and prose** — the tables, figures and README text under `03_results/`, and the
  documentation: [CC BY 4.0](LICENSE-CC-BY-4.0.txt).

Reuse of either requires attribution. The underlying sequence data carries the terms of its GEO
deposition, not these.
