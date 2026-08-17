# 09_gatom — The metabolic subnetwork each contrast recruits

GATOM searches an atom-transition graph of KEGG metabolism for the maximally-regulated connected
subnetwork: the set of enzymatic reactions whose genes carry the strongest coherent signal in a
contrast, joined into one graph. It runs on four contrasts — `WT_heat`, `KO_heat`, `Interaction`
and `Temp_main` — with gene scores from the limma differential expression and no metabolomics
input, so node size is uniform and every score on the graph is transcriptional.

Module sizes are small and they separate by contrast. `Temp_main` returns 6 nodes over 5 edges
at weight 48.32, `WT_heat` 7 over 6 at 47.01, and `KO_heat` 18 over 17 at 44.83. The interaction
returns an empty module. Only the KEGG reaction network was available for this run
(`met.combined.db.rds` is absent from `00_data/references/gatom/`), so the combined KEGG+Rhea
network contributes no bar to any panel.

A module is a statistically optimal connected subgraph of an atom-transition network. It records
where differential-expression signal concentrates in metabolic space. Metabolic flux is
untested.

---

## Figures

### `figures/by_contrast/<contrast>/module_graph_kegg.png` — four panels

**The recruited subnetwork drawn as a graph.**
Nodes are metabolites and atoms, at uniform size because no metabolomics data entered the run.
Edges are enzymatic reactions, labelled with the gene symbol encoding the enzyme. Edge colour
gives the gene's log2 fold change — orange where the enzyme rises in the contrast numerator,
blue where it falls — and edge width, where available, gives −log10 of the raw p-value from
GATOM's own scoring. The twenty most connected nodes carry labels. The subtitle reports the
module's node and edge count, its solution weight, and the mean signed log2 fold change across
its edges.

Module sizes: `WT_heat` V=7 E=6 w=47.01, `KO_heat` V=18 E=17 w=44.83, `Temp_main` V=6 E=5
w=48.32, `Interaction` empty.
*Source* `../objects/10_gatom_<contrast>.rds` · `02_analysis/scripts/14_gatom_viz.R`.

### `figures/_overview/module_sizes.png`

**How large a subnetwork each contrast recruits.**
Horizontal bars, one per contrast × network. x, number of reaction edges in the module. Blue
marks the KEGG reaction network, the only network in this run. An absent contrast is one whose
module came back empty.
*Source* `tables/_overview/module_sizes.csv` · `02_analysis/scripts/14_gatom_viz.R`.

### `figures/_overview/module_weights.png`

**The optimisation score behind each module.**
Bars, one per contrast × network. y, the solver's objective value, which combines module size
with the magnitude of its node and edge scores. A higher weight means more or larger enzyme
signal concentrated inside one connected subnetwork.
*Source* `tables/_overview/module_weights.csv` · `02_analysis/scripts/14_gatom_viz.R`.

### `figures/_overview/module_summary.png`

**Size and direction of every module on one canvas.**
Two stacked panels sharing the contrast axis. Top, bar height gives the module's reaction-edge
count. Bottom, dot position gives the mean log2 fold change across the module's enzyme edges, so
a dot above the dashed zero line marks a net up-regulated subnetwork in the contrast numerator.
Blue marks the KEGG network.
*Source* `tables/_overview/module_summary.csv` · `02_analysis/scripts/14_gatom_viz.R`.

---

## Tables

| File | What it holds |
|---|---|
| `tables/by_contrast/<contrast>/gatom_modules.csv` | The module for that contrast, edge by edge: reaction, gene symbol, log2 fold change, p-value and the node pair it joins. |
| `tables/_overview/module_sizes.csv` | One row per contrast × network with the node and edge counts. |
| `tables/_overview/module_weights.csv` | The same keys with the solver objective value. |
| `tables/_overview/module_summary.csv` | Size joined to the mean signed log2 fold change per module — the source of the two-panel summary. |

The three-row roll-up of these modules is carried in
[`../master/master_gatom_modules.csv`](../master/) under `database = GATOM_KEGG`, in the shared
schema, so they join the enrichment and activity tables.

---

## Signature provenance

| Resource | Provenance |
|---|---|
| The atom-transition networks | Frozen under `00_data/references/gatom/`. This run used the KEGG reaction network (`network.kegg.rds`, `met.kegg.db.rds`) with the mouse annotation `org.Mm.eg.gatom.anno.rds`. The combined KEGG+Rhea network needs `met.combined.db.rds`, which is absent from that directory, so it contributes no bar to any panel. |
| The gene scores | Per-gene log2 fold changes and p-values from the limma-trend fit on **GSE329522** — induced regulatory T cells from primary murine splenic CD4⁺ T cells, genotype × temperature, five libraries per cell. Made in [`../03_de/`](../03_de/). |
| Metabolite scores | None. The run carries no metabolomics input, so node size is uniform and every score on the graph is transcriptional. |
