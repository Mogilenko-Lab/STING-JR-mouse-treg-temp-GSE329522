# CoReSh compendium ranking + derived-set GSEA

This folder asks a specific question about the GSE329522 iTreg signatures: **where else,
across the public mouse transcriptome, do these same gene programs co-move?** CoReSh
(co-regulation search) sweeps each project signature across a compendium of ~85 mouse GEO
datasets and scores how strongly the signature's genes co-vary within each one (a
PCA-inspired `pctVar`). The top-ranked public datasets are then turned into data-driven
gene sets and run through the same GSEA machinery as every curated collection, so a reader
can place them side by side with GO / Hallmark / KEGG / Reactome.

Sample mapping for the mouse data is owner-confirmed (GSE329522, 2×2 genotype × temperature
iTreg), so nothing here carries a provisional-label caveat.

## What was queried

Four project signatures seed the search, exported from the signature-derivation step:

- **`WT_heat_up`** — genes up in wild-type iTregs at 39 °C vs 37 °C (a generic thermal
  program), at two stringencies (`fdr_logfc`, `fdr_only`).
- **`Interaction_up`** — genes whose 39-vs-37 °C induction is greater in wild-type than in
  cGAS-KO (the cGAS-dependent, genotype×temperature term), at the same two stringencies.

## Why the results live here and not with the GSEA databases

The sets scored in this folder (`database = CoReSh_derived`) are **data-driven
co-regulation modules** mined from public datasets, named `CORESH_<query>_<GSE>`. They exist
in no curated collection, so they are unique to this analysis — verified programmatically:
all 14 derived-set ids are disjoint from every GSEA-collection id in the master table. That
is why they populate this folder rather than the standard GSEA sweep.

## What the recovered datasets actually are

Because the compendium stores only variance structure (accession + platform, no titles),
each recovered dataset's identity was researched separately and frozen in
`tables/_overview/coresh_dataset_annotation.{csv,json}` — title, organism, tissue, the
perturbation compared, and neutral context flags (interferon / viral / cGAS-STING /
thermal), each with a source URL and PubMed id. Figure labels carry this identity, so a
point reads as `GSE89069 · viral infection · embryonic brain`, not a bare accession. The
`WT_heat_up` query surfaced a broad stress / metabolic / developmental / tumor mix; the
`Interaction_up` query surfaced viral-infection and nucleic-acid-sensing datasets. The
annotation is descriptive only and never feeds the enrichment statistics.

---

## figures/_overview/coresh_pctvar_overview.png

Top public mouse GEO datasets co-regulating each project signature
query, ranked by CoReSh pctVar.

**How to read:** Each bar = one public GEO dataset. Bar length = pctVar (% of that
dataset's variance explained by the query), a PCA-inspired
co-regulation score (higher = stronger co-regulation); label = GSE,
mapped query size k, rank. Facets separate the project signature
queries. Claim tier: L3-DE (compendium co-regulation score). pctVar
>= 0 (variance fraction, unsigned).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_pctvar_overview` | `coresh.top_n_hits=5; figures.top_pathways=20; thresholds.gsea_fdr=0.05` | `03_results/objects/coresh_ranked.rds` |

## figures/_overview/coresh_nes_dotplot.png

Cross-contrast enrichment of CoReSh-derived co-regulation sets (NES
fill, -log10(padj) size, FDR outline); rows grouped by seeding query
origin; labels carry each recovered dataset's identity.

**How to read:** Each circle = one CoReSh-derived gene set (row) in one contrast
(column). Row label = the recovered public GEO dataset (GSE · context
· tissue, from the dataset annotation). Fill = NES: orange positive
(up in numerator), blue negative; clamped at ±3.5. Size =
-log10(padj). Black outline = FDR < 0.05. Rows are grouped into
left-side facet bands by the project signature query that seeded each
set; within a band, sets are ordered by median NES (highest at top).
Claim tier: L3-DE (fgsea, BH-FDR).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_nes_dotplot` | `figures.top_pathways=20; figures.nes_cap=3.5; thresholds.gsea_fdr=0.05` | `03_results/master/master_gsea_table.csv rows database=CoReSh_derived` |

## figures/by_contrast/&lt;contrast&gt;/CoReSh_derived/ — the four-panel battery

For each of the seven contrasts, the `CoReSh_derived` sets get the same four panels every
curated database gets in the standard GSEA sweep, so CoReSh reads at parity with them. Each
panel is a vector PDF + raster PNG built from a `gseaResult` reconstructed viz-side: NES and
padj are taken verbatim from `master_gsea_table.csv`, the ranked list is the contrast's
limma-trend t-statistic vector, and the running score is computed deterministically (no
permutation re-run), so the curves never disagree with the tabulated NES.

- **`dotplot.png`** — every derived set, x = GeneRatio (leading-edge fraction), point size =
  −log10(q), fill = NES (orange up / blue down), black outline = FDR < 0.05.
- **`facet.png`** — the same dotplot split into Up (NES > 0) and Down (NES < 0) panels.
- **`barplot.png`** — NES bars for the FDR-significant sets only, top 20 by |NES|.
- **`running_sum.png`** — a real three-panel GSEA enrichment plot for the top 5 sets by
  |NES|: running enrichment score (top), gene-hit ticks at each set member's rank (middle,
  one colour per set), and the ranked t-statistic (bottom); ES y-range pinned to [−1, 1] so
  curves stay comparable across contrasts.

**How to read:** Row / legend label = the recovered public GEO dataset
(`GSE · context · tissue`). NES > 0 = the set is enriched among genes up in the contrast's
numerator; NES < 0 = enriched among the down side. Sign convention follows the contrast
(e.g. for `Interaction`, positive = stronger induction in wild-type than in cGAS-KO). Claim
tier: L3-DE (fgsea multilevel, BH-FDR). Per-contrast source table:
`tables/by_contrast/<contrast>/coresh_gsea.csv`.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `emit_coresh_cell` → `gsea_dotplot` / `gsea_dotplot_facet` / `gsea_barplot` / `gsea_running_sum_plot` | `figures.top_pathways=20; figures.nes_cap=3.5; figures.running_sum_top=5; thresholds.gsea_fdr=0.05` | `master_gsea_table.csv` (database=CoReSh_derived) + `02_de_results.rds` (ranks) + `coresh_derived_sets.rds` (sets) |

## Tables

| Table | What it is |
|---|---|
| `tables/coresh_ranked.csv` | Full compendium sweep: every query × every GEO dataset with pctVar, size, rank. |
| `tables/coresh_provenance.csv` | The top-ranked datasets that were turned into derived gene sets (query → GSE → rank/pctVar). |
| `tables/_overview/coresh_dataset_annotation.csv` | What each seeded dataset is (identity + context flags + PubMed / URL). Frozen research evidence: `…annotation.json`. Built by `08b_coresh_annotate.R`. |
| `tables/_overview/coresh_gsea_all.csv` | All `CoReSh_derived` GSEA rows across contrasts (NES, padj, leading edge). |
| `tables/_overview/coresh_gsea_summary.csv` | Per-contrast counts of significant derived sets and the top set. |
| `tables/by_contrast/<contrast>/coresh_gsea.csv` | Per-contrast derived-set GSEA (padj-ordered) — the source table for that contrast's battery. |

**Access note.** This analysis requires the ~20 GB mmu CoReSh compendium (syn66227307),
mounted read-only from the shared reference cache, before the search step can run.

