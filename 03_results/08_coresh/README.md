# 08_coresh — Where else these gene programs co-move

The gene-set sweep asks whether curated programs enrich in this contrast. This stage asks a
different question of the same signatures: **across the public mouse transcriptome, in which
datasets do these genes co-vary?** CoReSh sweeps each query across a compendium of roughly 85
public mouse GEO datasets and scores how tightly the query's genes move together within each one
(`pctVar`, a principal-component-inspired co-regulation score). The top-ranked datasets are then
projected back into gene-level loadings to give data-driven modules, which run through the same
enrichment machinery every curated collection does.

Four queries seed the search, taken from the signature-derivation step: `WT_heat_up` and
`Interaction_up`, each at both stringency gates (`fdr_logfc`, `fdr_only`). The fourteen modules
they produce are named `CORESH_<query>_<GSE>` and are unique to this analysis. All fourteen ids
are disjoint from every curated collection id in the master table, which is why they get a
directory of their own here.

**What the recovered datasets are.** The compendium stores variance structure and accessions,
with no titles. Each recovered dataset's identity was therefore researched separately and frozen
in `tables/_overview/coresh_dataset_annotation.{csv,json}`: title, organism, tissue, the
perturbation compared, and neutral context flags for interferon, viral, cGAS-STING and thermal
biology, each with a source URL and a PubMed id.

Figure labels carry that identity, so a point reads `GSE89069 · viral infection · embryonic
brain`. The `WT_heat_up` query surfaced a broad stress, metabolic, developmental and tumour mix.
The `Interaction_up` query surfaced viral-infection and nucleic-acid-sensing datasets. The
annotation is descriptive and enters no statistic.

A derived module is a co-regulation neighbourhood mined from one public dataset's variance
structure. It is named for how it was made.

**Access.** The search step reads the ~20 GB mmu CoReSh compendium (syn66227307), mounted
read-only from the shared reference cache.

---

## Figures

### `figures/_overview/coresh_pctvar_overview.png`

**Which public datasets carry each query's co-regulation.**
Horizontal bars, one per public GEO dataset, faceted by query. x, `pctVar`, the share of that
dataset's variance the query genes jointly explain — a co-regulation score, unsigned and always
positive, higher meaning tighter co-movement. The row label gives the accession, the mapped
query size and the within-query rank. `pctVar` is normalised by how many query genes the
platform measures, so it compares within a query.
*Source* `tables/coresh_ranked.csv` · `02_analysis/scripts/15_coresh_viz.R`.

### `figures/_overview/coresh_nes_dotplot.png`

**The derived modules scored back across all seven contrasts.**
One circle per module × contrast. Rows are the fourteen derived sets, labelled with the
recovered dataset's researched identity (accession · context · tissue) and grouped into
left-side facet bands by the query that seeded them; within a band, rows order by median NES.
Columns are the seven contrasts. Fill gives NES, orange positive toward the numerator and blue
negative, clamped to ±3.5; size gives −log10(adjusted p); a black outline marks FDR < 0.05.
*Source* `../master/master_gsea_table.csv` rows where `database = CoReSh_derived` ·
`02_analysis/scripts/15_coresh_viz.R`.

### `figures/by_contrast/<contrast>/CoReSh_derived/` — the four-panel battery

Every contrast gets the same four panels the curated databases get in
[`../06_gsea/`](../06_gsea/), so the derived sets read at parity with them. Each panel is a
`.pdf` and `.png` pair built from a `gseaResult` reconstructed viz-side: NES and adjusted p come
verbatim from `master_gsea_table.csv`, the ranked list is the contrast's moderated-t vector, and
the running score is recomputed deterministically, so a curve and its tabulated NES agree.

- **`dotplot.png`** — x, GeneRatio (leading-edge fraction); size, −log10(adjusted p); fill, NES;
  black outline, FDR < 0.05.
- **`facet.png`** — the same dotplot split into an NES > 0 block and an NES < 0 block.
- **`barplot.png`** — NES bars for the FDR-significant sets, top twenty by |NES|.
- **`running_sum.png`** — three stacked panels for the top five by |NES|: the running enrichment
  score with y pinned to [−1, 1], one tick row per set at each member's rank, and the ranked
  moderated t beneath.

Row and legend labels carry the recovered dataset's identity. NES > 0 places the module's genes
toward the contrast numerator. On the interaction, positive means stronger induction in
wild-type.
*Source* `tables/by_contrast/<contrast>/coresh_gsea.csv` ·
`02_analysis/scripts/15_coresh_viz.R`.

---

## Tables

| File | What it holds |
|---|---|
| `tables/coresh_ranked.csv` | The full compendium sweep: every query × every GEO dataset with `pctVar`, matched query size and within-query rank. |
| `tables/coresh_provenance.csv` | The top-ranked datasets that became modules — query, accession, rank, `pctVar`. |
| `tables/_overview/coresh_dataset_annotation.csv` | What each seeded dataset is: title, organism, tissue, perturbation, context flags, source URL and PubMed id. The frozen research evidence sits beside it in `…annotation.json`; both are built by `08b_coresh_annotate.R`. |
| `tables/_overview/coresh_gsea_all.csv` | Every `CoReSh_derived` enrichment row across contrasts, with NES, adjusted p and leading edge. |
| `tables/_overview/coresh_gsea_summary.csv` | Per contrast, how many derived sets reach significance and which one leads. |
| `tables/_overview/coresh_pctvar_overview.csv` | The plotted bars of the ranking panel. |
| `tables/_overview/coresh_nes_dotplot.csv` | The plotted points of the module dot plot, with the dataset annotation joined on. |
| `tables/by_contrast/<contrast>/coresh_gsea.csv` | Per-contrast derived-set enrichment, ordered by adjusted p — the source of that contrast's battery. |

---

## Signature provenance

| Resource | Provenance |
|---|---|
| The compendium | The public mouse GEO compendium distributed on Synapse as **syn66227307** (mmu, roughly 85 chunks, ~20 GB), mounted read-only from the shared reference cache. It stores variance structure and accessions. |
| The four queries | `WT_heat_up` and `Interaction_up`, each at the `fdr_logfc` and `fdr_only` gates, derived here from **GSE329522** — induced regulatory T cells from primary murine splenic CD4⁺ T cells, genotype × temperature, five libraries per cell. [`../10_signature/`](../10_signature/) records the gates. |
| The fourteen `CORESH_<query>_<GSE>` modules | Gene-level loadings projected back from the top-ranked datasets of that sweep. Each module is named for the query that seeded it and the GEO accession it was mined from, and each is unique to this analysis. |
| The dataset annotation | Researched per accession and frozen in `tables/_overview/coresh_dataset_annotation.{csv,json}`, each row carrying its source URL and PubMed id. Built by `08b_coresh_annotate.R`. It is descriptive and enters no statistic. |
| The rankings scored against | Moderated-t vectors from the limma-trend fit on the same twenty libraries, made in [`../03_de/`](../03_de/). |
