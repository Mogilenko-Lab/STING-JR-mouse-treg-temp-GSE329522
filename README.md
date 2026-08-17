# The mouse 39 °C iTreg anchor — GSE329522

Bulk RNA-seq of induced regulatory T cells cultured at 37 °C and at 39 °C, in wild-type and
cGAS-knockout genotypes.

This is the mouse half of a project asking whether human inflammatory and autoinflammatory disease
states carry transcriptional programs consistent with a temperature-stress axis in T cells and with
cGAS–STING-related biology. It derives the signatures; the human compartments score them.

## Design

Induced regulatory T cells differentiated from primary murine splenic CD4⁺ T cells, in a 2×2 design
of genotype (wild-type, cGAS-knockout) × temperature (37 °C, 39 °C), five biological replicates per
group, twenty libraries. Deposited as
[GSE329522](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE329522).

Seven contrasts follow from that design. The analysis fits them, reads the rankings against thirteen
gene-set databases and two regulatory networks, freezes the genes that rise at 39 °C into named
signatures, converts those to human symbols, and decomposes them against curated lenses.

## Findings

- **Temperature dominates.** PC1 carries 57.3 % of the variance and tracks temperature. Warming
  changes 8,723 genes in wild-type and 8,901 in cGAS-knockout cells at FDR 0.05.
- **Both genotypes warm in step.** Across 10,418 heat-responsive genes the knockout response is
  0.99× the wild-type response at r = 0.95.
- **The cGAS-dependence test is the weakest contrast** — the heat-by-genotype interaction, at one
  degree of freedom and n = 5, where 23 genes clear FDR 0.05. Every artifact here phrases a failure
  there as *no detectable cGAS-dependence at n = 5*.
- **Interferon and hypoxia diverge under that test.** The two Hallmark interferon-response sets carry
  the interaction at +3.14 and +2.83, adjusted p 2.7e-24 and 5.0e-22. HALLMARK_HYPOXIA rises in both
  heat arms at +1.91 and +1.95 and stays flat in the interaction at −1.27, adjusted p 0.17.

Each number is reproduced in a table under `03_results/`, named in that stage's `README.md`.

## The deliverable

`03_results/human_projection/` is the contract the human compartments read — three frozen arms in
human symbols, each with an up list, a down list, a signed ranked list, and a manifest of sizes and
gates. Genes are gated on FDR and \|log2FC\| ≥ 1 and converted with pinned offline babelgene 22.9,
with a ledger recording every gene's fate through that step.

| Signature | Mouse symbols | Human symbols |
|---|---|---|
| `WT_heat_up` | 213 | 202 |
| `KO_heat_up` | 239 | 221 |
| `Interaction_up` | 9 | 7 |

## Conventions

**Signs.** A positive statistic points to the numerator of its contrast: 39 °C for a heat contrast,
wild-type for a genotype contrast. For the interaction, positive means the heat response is larger in
wild-type.

**Naming.** A signature is named for how it was derived, so `WT_heat_up` is checkable against the
differential-expression tables. A gene set enriching locates its gene content in a ranking; whether
the program it is named for is present is measured separately, in the decomposition stages.

## Layout

| Path | Contents |
|---|---|
| `00_data/` | The normalized CPM matrix with its provenance record, and the symbol-alias map. Scripts fetch or cache the external reference sets; sequence data lives in GEO. |
| `02_analysis/scripts/` | 60 numbered scripts, 59 in R and one in Python. Compute scripts write tables; `_viz` scripts draw from them. |
| `02_analysis/config/` | `analysis_config.yaml` — every path and threshold the scripts read. |
| `02_analysis/helpers/` | Figure styling, ortholog conversion, symbol-alias resolution, source hashing. |
| `03_results/` | One directory per stage, each with `tables/`, `figures/` and a `README.md` captioning every file. |

`03_results/README.md` is the entry point: reading order stage by stage, each finding with the table
it comes from, and the provenance of every gene set, lens and network the analysis reads.

## Reproducing

Scripts run in numeric order inside the container under `.devcontainer/`, reading paths and
thresholds from `02_analysis/config/analysis_config.yaml`. Tables and stage READMEs are the tracked
record — figures regenerate from the script each caption names, since `.png` and `.pdf` are ignored
repo-wide. Clone with `--recursive` to populate `01_modules/`, which GitHub source tarballs carry as
empty directories.

## Environment

The analysis runs in `scdock-r-dev:v0.5.10`, pinned on the `dev-core` service in
`.devcontainer/docker-compose.yml`. That image is defined by
[scbio-docker](https://github.com/tony-zhelonkin/scbio-docker) at commit
[`5885cd3`](https://github.com/tony-zhelonkin/scbio-docker/commit/5885cd306ea908cb1949e7238b9186074b938953).

[RNAseq-toolkit](https://github.com/tony-zhelonkin/RNAseq-toolkit) supplies the GSEA plotters and the
prebuilt mouse gene-set objects — TransportDB, MitoCarta 3.0 and MitoXplorer 3.0 — read from its
reference tree. This release records it at commit
[`752481f`](https://github.com/tony-zhelonkin/RNAseq-toolkit/commit/752481fd13542ccb81d2b9b92ba57305cf13d6fc)
(`v0.2.0-9-g752481f`, on `dev`) under `01_modules/`.

## License

**MIT** for code — `02_analysis/`, `.devcontainer/`, and scripts anywhere in the tree
([LICENSE](LICENSE)). **CC BY 4.0** for results and prose — the tables, figures and README text under
`03_results/`, and the documentation ([LICENSE-CC-BY-4.0.txt](LICENSE-CC-BY-4.0.txt)). Both require
attribution. Sequence data carries the terms of its GEO deposition.
