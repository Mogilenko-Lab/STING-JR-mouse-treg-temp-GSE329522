# Results — the mouse 39 °C iTreg anchor

This is the experimental anchor of the project. Induced regulatory T cells were cultured at
37 °C and at 39 °C, in wild-type and cGAS-knockout genotypes, five libraries per cell of the
2×2 design (GSE329522). Everything downstream of this compartment consumes what these twenty
libraries produce.

The arc runs in one direction. Confirm the sample-to-condition mapping from expression alone,
fit the seven contrasts the design supports, read those rankings against thirteen gene-set
databases and two regulatory networks, then freeze the genes that rise at 39 °C into named
signatures and convert them to human symbols. The last third of the tree turns back on those
signatures and asks what they are made of.

## The reading order

| Stage | Question it answers |
|---|---|
| [`01_qc`](01_qc/) | Which library is which condition, read from marker expression alone. |
| [`02_eda`](02_eda/) | How much of the variance each design factor carries. |
| [`03_de`](03_de/) | Which genes move, in each of the seven contrasts. |
| [`04_tf`](04_tf/) | Which transcription factors the rankings nominate, and how far that nomination holds. |
| [`05_progeny`](05_progeny/) | The same question asked of fourteen signalling footprints. |
| [`06_gsea`](06_gsea/) | The full gene-set sweep: 13 databases × 7 contrasts. |
| [`07_synthesis`](07_synthesis/) | The four evidence arms joined into one cross-method table. |
| [`08_coresh`](08_coresh/) | Where else in public mouse data these genes co-vary. |
| [`09_gatom`](09_gatom/) | The metabolic subnetwork each contrast recruits. |
| [`10_signature`](10_signature/) | The candidate signatures, sized at two gates. |
| [`11_projection`](11_projection/) | Those signatures in human symbols, with every gene's fate recorded. |
| [`12_hsr_decomp`](12_hsr_decomp/) | What the frozen 39 °C set contains, against two curated lenses. |
| [`13_semantic_decomp`](13_semantic_decomp/) | What the same set looks like to the Gene Ontology graph. |
| [`14_signature_expression`](14_signature_expression/) | Each signature gene's behaviour in the four design cells. |
| [`15_go_decomposition`](15_go_decomposition/) | Which GO terms the arms organise into, against a depth-matched null. |
| [`16_arm_composition`](16_arm_composition/) | The arms partitioned against five frozen MSigDB collections. |
| [`human_projection`](human_projection/) | The frozen signature contract the human compartments read. |
| [`master`](master/) | Accumulator tables: every gene, set and factor in one place. |
| [`objects`](objects/) | Recomputable checkpoints, plus the gene-set vocabulary ledger. |
| [`interactive`](interactive/) | Browser dashboards over the same tables. |

## The headline result, and where to check it

Temperature dominates this design. PC1 carries 57.3% of the variance and tracks temperature
almost exactly ([`02_eda/tables/fig1c_pca_varexp.csv`](02_eda/tables/)), and warming changes
8,723 genes in wild-type and 8,901 in cGAS-knockout cells at FDR 0.05
([`03_de/tables/_overview/de_counts_summary.csv`](03_de/tables/_overview/)). The two genotypes
warm in step: over 10,418 heat-responsive genes the knockout response is 0.99× the wild-type
response at r = 0.95 ([`03_de/tables/_overview/cgas_dependence_stats.csv`](03_de/tables/_overview/)).

The cGAS-dependence test is the heat-by-genotype interaction, and it is the least-powered
contrast in the design at one degree of freedom and n = 5. Twenty-three genes clear FDR 0.05
there, all of them in one direction. A gene failing that test has **no detectable
cGAS-dependence at n = 5**, which is the phrasing every artifact in this tree uses.

The interferon and hypoxia arms behave differently under it, and that asymmetry is the mouse
result. The two Hallmark interferon-response sets carry the interaction (+3.14 and +2.83,
adjusted p 2.7e-24 and 5.0e-22), while HALLMARK_HYPOXIA rises in both heat arms and stays flat
in the interaction (+1.91 and +1.95 against −1.27 at adjusted p 0.17) — all four read from
[`master/master_gsea_table.csv`](master/).

## The frozen signatures

[`human_projection/`](human_projection/) is the contract. `WT_heat_up` carries 213 mouse
symbols at the FDR + |log2FC| ≥ 1 gate and arrives in human as 202
([`11_projection/tables/_overview/conversion_ledger.csv`](11_projection/tables/_overview/)).
`KO_heat_up` carries 239 and arrives as 221. `Interaction_up` carries 9 and arrives as 7.
Each ships as an up list, a down list and a signed ranked list, with a manifest recording sizes
and gates.

## Two conventions that hold everywhere in this tree

**Signs.** A positive statistic points to the numerator of its contrast: 39 °C for a heat
contrast, wild-type for a genotype contrast. For the interaction, positive means the heat
response is larger in wild-type. Orange draws positive and blue draws negative throughout.

**Naming.** A signature is named for how it was derived. `WT_heat_up` is the set of genes
higher at 39 °C than at 37 °C in wild-type cells, and that description is checkable against
[`03_de/tables/by_contrast/WT_heat/volcano.csv`](03_de/tables/by_contrast/). A gene set
enriching locates its gene content in a ranking; establishing that the program the set is named
for is present is a separate measurement, and [`15_go_decomposition/`](15_go_decomposition/)
and [`16_arm_composition/`](16_arm_composition/) are where it is made.

## Layout

Every stage holds `figures/` and `tables/`, with `_overview/` for cross-contrast artifacts and
`by_contrast/<contrast>/` where a stage runs per contrast. A figure's source table is its
same-stem neighbour under `tables/`, so `figures/_overview/foo.png` is drawn from
`tables/_overview/foo.csv` and nothing is computed at draw time.

`.png` and `.pdf` files are ignored repo-wide. Each stage's tables and its README are the
tracked record, and every figure regenerates from the script its caption names.
