# 06_gsea — artifact captions

## Stage overview

06_gsea GSEA stage: MSigDB collections (Hallmark, KEGG, Reactome, WikiPathways,
GO_BP, GO_MF, GO_CC, TF_Targets) plus 4 custom databases (TransportDB, MitoPathways,
MitoXplorer, Lombardi2022_HIF) — 12 databases total — across 7 contrasts. Per-contrast
figures in `figures/by_contrast/<contrast>/<DB>/` are built with RNAseq-toolkit plotters
(`gsea_dotplot`, `gsea_dotplot_facet`, `gsea_barplot`, `gsea_running_sum_plot`) on a
`gseaResult` S4 object reconstructed viz-side from the master table + cached DE ranks
+ gene-set lists; the `running_sum` is a real `enrichplot::gseaplot2` three-panel ES
curve, never a lollipop proxy. Cross-contrast overview panels in `figures/_overview/`.
Produced by `12_gsea_viz.R` (VIZ-ONLY; GSEA computed by `05_gsea_msigdb_run.R` and
`06_gsea_custom_run.R`).

**Key biology framing:**
- IFN/ISG arm: expected cGAS-dependent — positive NES in WT_heat and positive Interaction NES.
- HIF/glycolysis arm: no detectable cGAS-dependence at n=5 — non-significant Interaction NES.
- NEVER write "cGAS-independent" or "cGAS-independent/parallel" — always "no detectable
  cGAS-dependence at n=5" (the 1-df Interaction term is underpowered at n=5/group).
- NEVER crown HIF1alpha/HIF2alpha as the driver; NES is an enrichment statistic (L3 tier).
- Lombardi-2022 48-gene consensus (doi:10.1016/j.celrep.2022.111652) is the published
  HIF orthogonal check; the bespoke 16-gene list is ~92% heat-shock/glycolytic and its
  hypoxia-diagnostic core (Pdk1/Bnip3/Bnip3l/Car9) is repressed, not induced.
- Sample labels: PROVISIONAL (groups inferred from marker genes Hspa1b + Cgas, pending
  collaborator sample sheet).

**NES sign convention (applies everywhere in this stage):**
NES > 0 = gene set enriched in the numerator of the contrast (hot = 39 °C, WT genotype).
NES < 0 = gene set enriched in the denominator (cold = 37 °C, cGAS-KO genotype).
Interaction NES > 0 = the enrichment is stronger in WT than cGAS-KO (cGAS-dependent induction).

**Glyph convention (applies everywhere in this stage):**
Orange = up / positive NES. Blue = down / negative NES. * or filled point = padj < 0.05.
Black ring/outline on dot = padj < 0.05 (significant). Open shape or no outline = not significant.

**Claim tier:** L3 (DE/enrichment statistics). Never L7 (mechanism) in a figure title.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `save_overview / save_figure` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv` |

---

## Per-contrast × per-database figure families

Each cell in the 7-contrast × 12-database grid holds exactly four panel files, each
saved as a `.pdf` (vector) and `.png` (raster) pair. The wildcard `<contrast>` ranges over:

| Contrast key | Display label | Design role |
|---|---|---|
| `WT_heat` | WT: heat (39 vs 37 °C) | Per-genotype heat effect, WT arm |
| `KO_heat` | cGAS-KO: heat (39 vs 37 °C) | Per-genotype heat effect, cGAS-KO arm |
| `Interaction` | Heat × genotype interaction (cGAS-dependence of the heat response) | 1-df interaction term; the cGAS-dependence payoff; lowest power (n=5) |
| `Geno_at_39` | Genotype (WT vs cGAS-KO) @ 39 °C | Genotype simple effect at heat |
| `Geno_at_37` | Genotype (WT vs cGAS-KO) @ 37 °C | Genotype simple effect at baseline |
| `Temp_main` | Heat main effect (pooled genotypes) | Pooled 39 °C-versus-37 °C response; HIF/glycolysis arm clearest here |
| `Geno_main` | Genotype main effect (pooled temps) | Average cGAS effect across temperatures |

The wildcard `<DB>` ranges over:

| Database | Source | Content |
|---|---|---|
| `Hallmark` | MSigDB | 50 curated hallmark gene sets (MSigDB H collection) |
| `KEGG` | MSigDB | KEGG pathway gene sets (MSigDB C2:CP:KEGG) |
| `Reactome` | MSigDB | Reactome pathway gene sets (MSigDB C2:CP:REACTOME) |
| `WikiPathways` | MSigDB | WikiPathways gene sets (MSigDB C2:CP:WIKIPATHWAYS) |
| `GO_BP` | MSigDB | Gene Ontology Biological Process (MSigDB C5:GO:BP) |
| `GO_MF` | MSigDB | Gene Ontology Molecular Function (MSigDB C5:GO:MF) |
| `GO_CC` | MSigDB | Gene Ontology Cellular Component (MSigDB C5:GO:CC) |
| `TF_Targets` | MSigDB | Transcription factor target gene sets (MSigDB C3:TFT) |
| `TransportDB` | Custom | Membrane transporter gene sets; metabolic/transport focus; no IFN/ISG sets |
| `MitoPathways` | Custom | Mitochondrial pathway gene sets; metabolic focus; no IFN/ISG sets |
| `MitoXplorer` | Custom | MitoXplorer metabolic atlas gene sets; metabolic focus; no IFN/ISG sets |
| `Lombardi2022_HIF` | Custom | Lombardi-2022 48-gene published HIF consensus (doi:10.1016/j.celrep.2022.111652) |

All 84 (7 × 12) cells were attempted; cells with 0 master rows or no gene-set list on disk
were skipped. Gene-set size filter: 15–500 genes (GSEA_MIN_SIZE / GSEA_MAX_SIZE from config).
Statistical threshold: FDR (Benjamini-Hochberg padj) < 0.05. Ranking metric: t-statistic
from limma-trend DE (02_de_results.rds). NES and padj values are taken verbatim from
`master_gsea_table.csv` (fgsea compute engine); running-sum ES curves are computed
deterministically by `enrichplot::gseaScores(geneList, geneSet, exponent=1)` on the
reconstructed gseaResult (NES-agreement check: max|dNES| < 0.04, r > 0.9999, 100% sign
agreement on WT_heat/Hallmark benchmark cell).

The four panel types are described once below. Each applies identically across all 84 cells.

---

## figures/by_contrast/<contrast>/<DB>/running_sum.png

Real three-panel GSEA enrichment curve for the top 5 pathways by |NES| in the given
(contrast × database) cell, rendered via `enrichplot::gseaplot2()` through the
RNAseq-toolkit `gsea_running_sum_plot()` wrapper. Each panel set shows:
panel 1 (top) = running enrichment score (ES) line with the leading-edge peak annotated;
panel 2 (middle) = gene-hit tick marks at each set member's position in the ranked list;
panel 3 (bottom) = the ranked t-statistic (ranking metric) across all genes. Only pathways
with a cached gene-set membership list (geneset_\<src\>_\<DB\>.rds) appear; pathways
without a list are excluded from the running-sum (they still appear in dotplot/barplot).

**How to read:** Each coloured trace is one pathway (up to top 5 by |NES|, labelled in the
legend with the formatted pathway name). The ES y-axis is pinned to [−1.0, 1.0] across all
databases so curves are directly comparable (config `figures.running_sum_ylim`). The
leading-edge peak (maximum |ES|) determines NES sign: peak to the left = positive NES
(genes concentrated at the top/numerator end of the ranking); peak to the right = negative
NES (genes concentrated at the bottom/denominator end). The gene-hit tick panel (middle)
shows where each set member falls in the ranked gene list — dense ticks near rank 0 confirm
a strongly up-enriched set. The ranked metric panel (bottom) shows the t-statistic landscape:
positive = up in numerator, negative = up in denominator. NES>0 = enriched in hot/WT
numerator side. Claim tier: L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_running_sum_plot` (wraps `enrichplot::gseaplot2`) | `figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5` | `03_results/master/master_gsea_table.csv` + `03_results/objects/02_de_results.rds` + `03_results/objects/geneset_<src>_<DB>.rds` |

---

## figures/by_contrast/<contrast>/<DB>/dotplot.png

Dot plot showing the top 20 pathways by padj for the given (contrast × database) cell,
rendered via the RNAseq-toolkit `gsea_dotplot()`. Each dot is one pathway. All pathways
from the master table are shown (significant and not), ordered by FDR. When no pathways
are FDR-significant, a "No significant pathways (FDR < 0.05)" label annotation is overlaid
in the upper-left corner.

**How to read:** x-axis = GeneRatio (number of leading-edge genes / gene-set size; a larger
GeneRatio means the leading edge is a larger fraction of the set). Point size = −log10(FDR
padj); larger dots = more significant. Point fill = NES on a diverging scale: orange
(positive NES = enriched in numerator/hot/WT) through white (NES ≈ 0) to blue (negative
NES = enriched in denominator/cold/KO); NES capped at ±3.5. Black ring outline = padj <
0.05 (FDR-significant). No outline = not significant. y-axis = pathway names (formatted by
`format_pathway_name()`; database prefix stripped). Pathways sorted by padj ascending
(most significant at top). Top 20 pathways shown (`figures.top_pathways=20`). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv` + `03_results/objects/02_de_results.rds` + `03_results/objects/geneset_<src>_<DB>.rds` |

---

## figures/by_contrast/<contrast>/<DB>/facet.png

Faceted dot plot splitting the top 20 pathways into an Up (NES > 0) panel and a Down
(NES < 0) panel, rendered via the RNAseq-toolkit `gsea_dotplot_facet()`. Glyph encoding
is identical to the `dotplot` sibling; the facet layout makes the directional structure
immediately visible when sets cluster in both quadrants.

**How to read:** Two-panel layout: left/top facet = pathways enriched in numerator (NES > 0,
orange gradient); right/bottom facet = pathways enriched in denominator (NES < 0, blue
gradient). x-axis = GeneRatio; point size = −log10(padj); fill = NES (diverging
orange–white–blue, capped at ±3.5); black ring = padj < 0.05. y-axis = formatted pathway
names. Top 20 pathways per direction shown. This view is most useful when both up- and
down-enriched sets are present simultaneously (e.g., WT_heat Hallmark where IFN sets are
up and metabolic-suppression sets are down). Claim tier: L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot_facet` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv` + `03_results/objects/02_de_results.rds` + `03_results/objects/geneset_<src>_<DB>.rds` |

---

## figures/by_contrast/<contrast>/<DB>/barplot.png

NES bar chart for FDR-significant pathways only (padj < 0.05), top 20 by |NES|, rendered
via the RNAseq-toolkit `gsea_barplot()`. Bars extend left (negative NES, blue) or right
(positive NES, orange) from zero. If no pathways are FDR-significant, the panel is blank
or carries a "no significant pathways" annotation.

**How to read:** x-axis = NES (capped at ±3.5); bars extend from 0 toward positive
(orange = enriched in numerator/hot/WT) or negative (blue = enriched in
denominator/cold/KO). y-axis = formatted pathway names sorted by NES descending (most
positive at top). Only padj < 0.05 sets are shown (at most top 20 by |NES|). Bar length
encodes effect size (NES magnitude); all bars are significant, so no outline-ring encoding
is needed here. This is the most compact high-signal view: if a cell's barplot is empty,
no pathways survived FDR correction in that (contrast × database) cell. Claim tier: L3.
PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_barplot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv` + `03_results/objects/02_de_results.rds` + `03_results/objects/geneset_<src>_<DB>.rds` |

---

## Cross-contrast overview figures

The two surviving overview panels below summarise results across all contrasts and/or
databases. Each is saved as a `.pdf` (vector) and `.png` (raster) pair.

## figures/_overview/gsea_asymmetry_panel.png

cGAS-dependence asymmetry in Hallmark GSEA (focal contrasts: WT_heat,
KO_heat, Interaction, Temp_main). IFN/immune sets (top block) are
expected to show positive NES in WT_heat and positive Interaction NES
(cGAS-dependent heat-induced IFN arm). HIF/glycolysis sets (bottom
block) are expected positive in WT_heat/KO_heat/Temp_main but NOT
significantly enriched in Interaction — interpreted as no detectable
cGAS-dependence at n=5, NOT as proven cGAS-independence (the 1-df
Interaction term is underpowered).

**How to read:** Two-row heatmap: top block = IFN/immune Hallmark sets; bottom block =
HIF/glycolysis sets. Color = NES (orange = NES>0 / enriched in hot/WT
numerator; blue = NES<0 / enriched in cold/KO denominator). * = padj
< 0.05. The Interaction contrast is the cGAS-dependence read-out:
positive/significant Interaction NES for IFN sets = cGAS-dependent;
non-significant Interaction NES for HIF sets = no detectable
cGAS-dependence at n=5 (NEVER 'cGAS-independent'). Claim tier: L3.
PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `geom_tile` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv` |

## figures/_overview/gsea_lombardi_vs_bespoke_hif.png

Comparison of HIF-signature GSEA NES across contrasts: published
Lombardi-2022 48-gene consensus (orange) vs Hallmark HYPOXIA
(purple). Both are expected to track the heat response in
WT_heat/KO_heat/Temp_main but NOT show significant Interaction NES —
consistent with no detectable cGAS-dependence of the HIF/glycolysis
arm at n=5. The Lombardi-2022 set provides an orthogonal published
benchmark absent from the original hand-curated 16-gene list (~92%
heat-shock/glycolytic; the hypoxia-diagnostic core repressed, not
induced).

**How to read:** Line-point traces show NES across contrasts for each HIF signature.
Filled points = padj < 0.05; open points = not significant. Positive
NES = enriched in numerator (hot/WT); negative = denominator. The
Interaction contrast (if present) is the cGAS-dependence read-out:
non-significant NES for the HIF arm = no detectable cGAS-dependence
at n=5 (NEVER 'cGAS-independent'). House rule: never crown
HIF1α/HIF2α as driver; NES is an enrichment statistic (L3 tier).
PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `geom_point + geom_line` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv` |

## figures/by_contrast/<contrast>/Hallmark/*.png

Hallmark GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_Hallmark.rds).
IFN/ISG and immune sets enrich positively in WT-time contrasts and
reverse sign in the STING/IFNAR1 interaction contrasts.

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_Hallmark.rds}` |

## figures/by_contrast/<contrast>/KEGG/*.png

KEGG GSEA per contrast (dotplot/facet/barplot/running_sum) via the
RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_KEGG.rds).
IFN/ISG and immune sets enrich positively in WT-time contrasts and
reverse sign in the STING/IFNAR1 interaction contrasts.

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_KEGG.rds}` |

## figures/by_contrast/<contrast>/Reactome/*.png

Reactome GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_Reactome.rds).
IFN/ISG and immune sets enrich positively in WT-time contrasts and
reverse sign in the STING/IFNAR1 interaction contrasts.

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_Reactome.rds}` |

## figures/by_contrast/<contrast>/WikiPathways/*.png

WikiPathways GSEA per contrast (dotplot/facet/barplot/running_sum)
via the RNAseq-toolkit plotters on a viz-side-reconstructed
gseaResult (NES/padj taken verbatim from master_gsea_table.csv;
ranked vector from 02_de_results.rds; gene sets from
geneset_msigdb_WikiPathways.rds). IFN/ISG and immune sets enrich
positively in WT-time contrasts and reverse sign in the STING/IFNAR1
interaction contrasts.

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_WikiPathways.rds}` |

## figures/by_contrast/<contrast>/GO_BP/*.png

GO_BP GSEA per contrast (dotplot/facet/barplot/running_sum) via the
RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_GO_BP.rds).
IFN/ISG and immune sets enrich positively in WT-time contrasts and
reverse sign in the STING/IFNAR1 interaction contrasts.

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_GO_BP.rds}` |

## figures/by_contrast/<contrast>/GO_MF/*.png

GO_MF GSEA per contrast (dotplot/facet/barplot/running_sum) via the
RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_GO_MF.rds).
IFN/ISG and immune sets enrich positively in WT-time contrasts and
reverse sign in the STING/IFNAR1 interaction contrasts.

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_GO_MF.rds}` |

## figures/by_contrast/<contrast>/GO_CC/*.png

GO_CC GSEA per contrast (dotplot/facet/barplot/running_sum) via the
RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_GO_CC.rds).
IFN/ISG and immune sets enrich positively in WT-time contrasts and
reverse sign in the STING/IFNAR1 interaction contrasts.

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_GO_CC.rds}` |

## figures/by_contrast/<contrast>/TF_Targets/*.png

TF_Targets GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from
geneset_msigdb_TF_Targets.rds). IFN/ISG and immune sets enrich
positively in WT-time contrasts and reverse sign in the STING/IFNAR1
interaction contrasts.

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_TF_Targets.rds}` |

## figures/by_contrast/<contrast>/TransportDB/*.png

TransportDB GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from
geneset_custom_TransportDB.rds). per-contrast set enrichment for this
metabolic/transport database — read each set's direction on its own,
NOT as an IFN-axis signature (these DBs carry no IFN/ISG sets).

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_TransportDB.rds}` |

## figures/by_contrast/<contrast>/MitoPathways/*.png

MitoPathways GSEA per contrast (dotplot/facet/barplot/running_sum)
via the RNAseq-toolkit plotters on a viz-side-reconstructed
gseaResult (NES/padj taken verbatim from master_gsea_table.csv;
ranked vector from 02_de_results.rds; gene sets from
geneset_custom_MitoPathways.rds). per-contrast set enrichment for
this metabolic/transport database — read each set's direction on its
own, NOT as an IFN-axis signature (these DBs carry no IFN/ISG sets).

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_MitoPathways.rds}` |

## figures/by_contrast/<contrast>/MitoXplorer/*.png

MitoXplorer GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from
geneset_custom_MitoXplorer.rds). per-contrast set enrichment for this
metabolic/transport database — read each set's direction on its own,
NOT as an IFN-axis signature (these DBs carry no IFN/ISG sets).

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_MitoXplorer.rds}` |

## figures/by_contrast/<contrast>/Lombardi2022_HIF/*.png

Lombardi2022_HIF GSEA per contrast
(dotplot/facet/barplot/running_sum) via the RNAseq-toolkit plotters
on a viz-side-reconstructed gseaResult (NES/padj taken verbatim from
master_gsea_table.csv; ranked vector from 02_de_results.rds; gene
sets from geneset_custom_Lombardi2022_HIF.rds). IFN/ISG and immune
sets enrich positively in WT-time contrasts and reverse sign in the
STING/IFNAR1 interaction contrasts.

**How to read:** dotplot: x = GeneRatio (leading-edge / set size), point size =
-log10(FDR), fill = NES (orange #B35806 up / blue #2166AC down),
black outline = padj < 0.05. facet: same dotplot split into Up
(NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets
only, top 20 by |NES|. running_sum: a REAL three-panel GSEA
enrichment curve (top = running enrichment score with the
leading-edge peak; middle = gene-hit ticks at each set member's rank;
bottom = the ranked t-statistic) for the top 5 sets by |NES|; ES
y-range pinned to [-1.0,1.0] so curves are comparable across
databases. NES>0 = enriched in numerator (hot/WT side). Claim tier:
L3. PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_Lombardi2022_HIF.rds}` |

## README overview

06_gsea GSEA stage: MSigDB 8 collections + 4 custom databases
(TransportDB, MitoPathways, MitoXplorer, Lombardi2022_HIF) across 7
contrasts. Per-contrast figures in by_contrast/<c>/<DB>/ are built
with the RNAseq-toolkit plotters
(gsea_dotplot/facet/barplot/running_sum) on a gseaResult
reconstructed viz-side from the master table + cached DE ranks +
gene-set lists; the running_sum is a real enrichplot::gseaplot2
three-panel ES curve. Cross-contrast overviews in _overview/.
Produced by 12_gsea_viz.R (VIZ-ONLY; GSEA computed by
05/06_gsea_*_run.R). Key biology: IFN/ISG arm is expected
cGAS-dependent (positive WT_heat + Interaction NES); HIF/glycolysis
arm shows no detectable cGAS-dependence at n=5 (non-significant
Interaction NES). NEVER write 'cGAS-independent' — always 'no
detectable cGAS-dependence at n=5'. Lombardi-2022 48-gene consensus
(doi:10.1016/j.celrep.2022.111652) is the published HIF orthogonal
check. Claim tier: L3. Sample labels: PROVISIONAL.

**How to read:** NES > 0 = enriched in numerator side of the contrast (hot=39 degrees
C, WT genotype). NES < 0 = enriched in denominator side (cold=37
degrees C, cGASKO). Glyph convention: orange = up/positive, blue =
down/negative (consistent with colors.diverging in config). * or
filled point = padj < 0.05. Interaction contrast: (WT_heat) -
(KO_heat); positive Interaction NES = cGAS-dependent induction.
PROVISIONAL sample labels: group identity inferred from marker genes
(Hspa1b thermometer + Cgas), pending collaborator sample sheet. Claim
tier: L3 (DE/enrichment statistics). Never L7 (mechanism) in a figure
title.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `save_overview / save_figure` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.5; colors.diverging` | `03_results/master/master_gsea_table.csv` |
