# 06_gsea: artifact captions

## Stage overview

06_gsea GSEA stage: 8 MSigDB collections (Hallmark, KEGG, Reactome, WikiPathways,
GO_BP, GO_MF, GO_CC, TF_Targets) plus 5 custom databases (TransportDB, MitoPathways,
MitoXplorer, HSR_lens, TCR_activation), 13 databases in total, across
7 contrasts. Per-contrast
figures in `figures/by_contrast/<contrast>/<DB>/` are built with RNAseq-toolkit plotters
(`gsea_dotplot`, `gsea_dotplot_facet`, `gsea_barplot`, `gsea_running_sum_plot`) on a
`gseaResult` S4 object reconstructed viz-side from the master table plus cached DE ranks
plus gene-set lists; the `running_sum` is a real `enrichplot::gseaplot2` three-panel ES
curve. Cross-contrast overview panels in `figures/_overview/`.
Produced by `12_gsea_viz.R` (VIZ-ONLY; GSEA computed by `05_gsea_msigdb_run.R` and
`06_gsea_custom_run.R`). The `gsea_pooled_overview*` panels in `_overview/` come from
`12b_gsea_overview_pooled_viz.R`, not from `12_gsea_viz.R`.

**What the Hallmark interferon and hypoxia/metabolic sets actually do here**
(read off `master_gsea_table.csv`, so these are checkable):
- The two interferon-response sets are NEGATIVE in both per-genotype heat contrasts
  (INTERFERON_ALPHA_RESPONSE: WT_heat NES −1.90 padj 1.7e-04, KO_heat −2.52 padj 8.0e-12;
  INTERFERON_GAMMA_RESPONSE: −1.77 and −2.20) and strongly POSITIVE in the interaction term
  (+3.14 padj 2.7e-24 and +2.83 padj 5.0e-22). The interaction is positive because the
  decrease is deeper in cGAS-KO, not because heat raises them in WT.
- HYPOXIA, GLYCOLYSIS and MTORC1_SIGNALING are positive and FDR-significant in both heat
  contrasts (HYPOXIA WT_heat +1.91 padj 4.2e-06, KO_heat +1.95 padj 2.5e-06) with a
  negative, non-significant interaction NES (HYPOXIA −1.27, padj 0.17).
  OXIDATIVE_PHOSPHORYLATION is negative and never significant.
- A non-significant interaction NES means **no detectable cGAS-dependence at n=5**. It is
  never evidence of cGAS-independence: the 1-df interaction term is the least-powered
  contrast in this design.
- A set enriching is not evidence that the program its name invokes is present. Composition
  has to be established separately (see `03_results/15_go_decomposition/` and
  `03_results/16_arm_composition/`).
- NEVER write "cGAS-independent" or "cGAS-independent/parallel". Always write "no detectable
  cGAS-dependence at n=5" (the 1-df Interaction term is underpowered at n=5/group).
- NEVER crown HIF1alpha/HIF2alpha as the driver; NES is an enrichment statistic (L3 tier).
  This extends to naming: a gene set or a block of gene sets is named by how it was curated,
  never after a regulator it is hoped to represent.
- The `Lombardi2022_HIF` database was removed from this stage. A pan-cancer consensus derived
  under hypoxia in cancer cells is the wrong reference for a 39 °C in-vitro contrast in
  iTregs, and it meets `WT_heat_up` at 7 of that arm's 213 genes, 3.3% (see
  `mouse_anchor/03_results/12_hsr_decomp/README.md`). Hypoxia in this stage is read off the
  curated, versioned collections instead; `figures/_overview/hypoxia_routes_by_contrast.png`
  is the panel that does it. The curation script and its outputs under
  `00_data/references/gene_sets/` stay, because the curated comparison lens in
  `mouse_anchor/03_results/12_hsr_decomp/` and the semantic decomposition read the set from
  there and never from this stage.
- **A gene-set size is read at runtime and never written as a literal.** An earlier revision
  printed "48-gene" on a figure face for a set that held 100, because the count had been
  hardcoded from a single worksheet of a supplement. A hardcoded size is how a wrong label
  survives a re-curation.
- The bespoke 16-gene HIF list is ~92% heat-shock/glycolytic and its hypoxia-diagnostic core
  (Pdk1/Bnip3/Bnip3l/Car9) is repressed. Hypoxia in this stage is read off the curated,
  versioned collections instead of off that list.
- Sample labels: confirmed against the owner's sample sheet (2026-07-22), 20 of 20
  libraries concordant with the label-blind marker call from Hspa1b and Cgas.

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
| `02_analysis/scripts/12_gsea_viz.R` | `save_overview / save_figure` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv` |

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
| `HSR_lens` | Custom | Curated heat-shock-response lens, 2 sets (`HSR_core`, `HSR_sensitivity`); built by `00d`/`00e` from msigdbr v2026.1.Hs |
| `TCR_activation` | Custom | Curated TCR/immediate-early T-cell-activation lens, 1 set; built by `00f` |

All 91 (7 × 13) cells were attempted and all 91 emitted their four panels, so nothing was
skipped for want of master rows or a cached gene-set list: 728 PDF/PNG files = 91 cells × 4
panels × 2 formats. `figures/_overview/` holds 12 more files (6 panels × 2 formats):
`gsea_asymmetry_panel` written by `12_gsea_viz.R`, the three `gsea_pooled_overview*` panels
written by `12b_gsea_overview_pooled_viz.R`, and the two hypoxia panels written by
`28_hypoxia_focus_viz.R`. Gene-set size filter:
15–500 genes (GSEA_MIN_SIZE / GSEA_MAX_SIZE from config).
Statistical threshold: FDR (Benjamini-Hochberg padj) < 0.05. Ranking metric: t-statistic
from limma-trend DE (02_de_results.rds). NES and padj values are taken verbatim from
`master_gsea_table.csv` (fgsea compute engine); running-sum ES curves are computed
deterministically by `enrichplot::gseaScores(geneList, geneSet, exponent=1)` on the
reconstructed gseaResult (NES-agreement check: max|dNES| < 0.04, r > 0.9999, 100% sign
agreement on WT_heat/Hallmark benchmark cell).

The four panel types are described once below. Each applies identically across all 91 cells.

**The four panels use two different ranking metrics, and that is what governs every absence
you will notice.** `dotplot` and `facet` select the top 20 by **adjusted p**; `barplot` and
`running_sum` select by **|NES|**. Worked example, reproducible from
`03_results/master/master_gsea_table.csv`: in the `WT_heat` × `Hallmark` cell,
`HALLMARK_HYPOXIA` is **6th by |NES| (+1.908)** and **4th by adjusted p (4.16e-06)**, so it
falls outside the running-sum's top 5 by |NES| while sitting comfortably inside the dotplot's
top 20 by adjusted p. `HALLMARK_ANGIOGENESIS` runs the other way: 2nd by NES (+2.119) but
10th by adjusted p (1.84e-04, behind the two interferon sets that tie at 1.67e-04). An
absence from one panel is a statement about that panel's ranking metric and nothing else.

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
shows where each set member falls in the ranked gene list; dense ticks near rank 0 confirm
a strongly up-enriched set. The ranked metric panel (bottom) shows the t-statistic landscape:
positive = up in numerator, negative = up in denominator. NES>0 = enriched in hot/WT
numerator side. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_running_sum_plot` (wraps `enrichplot::gseaplot2`) | `figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2` | `03_results/master/master_gsea_table.csv` + `03_results/objects/02_de_results.rds` + `03_results/objects/geneset_<src>_<DB>.rds` |

---

## figures/by_contrast/<contrast>/<DB>/dotplot.png

Dot plot of the top 20 pathways by adjusted p for the given (contrast × database) cell,
rendered via the RNAseq-toolkit `gsea_dotplot()`. Each dot is one pathway. Both significant
and non-significant sets appear. When no pathway is FDR-significant, a "No significant
pathways (FDR < 0.05)" label annotation is overlaid in the upper-left corner.

**How to read:** x-axis = GeneRatio (number of leading-edge genes / gene-set size; a larger
GeneRatio means the leading edge is a larger fraction of the set). Point size = −log10(FDR
padj); larger dots = more significant. Point fill = NES on a diverging scale: orange
(positive NES = enriched in numerator/hot/WT) through white (NES ≈ 0) to blue (negative
NES = enriched in denominator/cold/KO); NES capped at ±3.2. Black ring outline = padj <
0.05 (FDR-significant). No outline = not significant. y-axis = pathway names (formatted by
`format_pathway_name()`; database prefix stripped). **Selection and display order are two
different things here:** the 20 sets drawn are the top 20 by adjusted p
(`filterBy = "p.adjust"`, `figures.top_pathways=20`), but the y-axis is then ordered by
**GeneRatio descending** (`sortBy = "GeneRatio"`, the toolkit default), so vertical position
is not a significance ranking. Read significance off dot size and the black ring, never off
row order. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv` + `03_results/objects/02_de_results.rds` + `03_results/objects/geneset_<src>_<DB>.rds` |

---

## figures/by_contrast/<contrast>/<DB>/facet.png

Faceted dot plot splitting the cell's sets into an Up (NES > 0) panel and a Down (NES < 0)
panel, rendered via the RNAseq-toolkit `gsea_dotplot_facet()`. Glyph encoding is identical to
the `dotplot` sibling; the facet layout makes the directional structure immediately visible
when sets cluster in both quadrants.

**Capped at 10 sets per direction, drawn on the 13 × 6.5 in landscape geometry.** Reusing
`figures.top_pathways` (20) here put up to 40 wrapped set names on one canvas and the labels
collided into illegibility, worst in `KO_heat` × `WikiPathways`. Rows are dropped rather than
labels clipped or the canvas stretched vertically: a tall canvas reads at the back of a room
but fails the journal-column half of the same legibility standard. Capping alone was not
enough, because 20 rows only fit 6.5 in of height if each set name occupies **one** line, and
at 8.5 in wide the long Reactome, MitoPathways and WikiPathways names still broke onto two.
The panel is therefore drawn wide with `wrap_width = 85`, which buys the label column room to
keep every name on one line without truncating any of it. Verified against measured label
lengths on the three worst cells: `KO_heat` × `Reactome` (longest label 81 characters),
`Geno_at_37` × `MitoPathways` (71), `KO_heat` × `WikiPathways` (62).

**The panel prints its own cap and the full per-direction count** in the subtitle (for example
"Top 10 per direction of 271 up and 146 down"), so the reduction is visible on the figure
rather than buried here, and it names the table that holds every set.

**How to read:** Two-panel layout: left/top facet = pathways enriched in numerator (NES > 0,
orange gradient); right/bottom facet = pathways enriched in denominator (NES < 0, blue
gradient). x-axis = GeneRatio; point size = −log10(padj); fill = NES (diverging
orange–white–blue, capped at ±3.2); black ring = padj < 0.05. y-axis = formatted pathway
names. Top 10 per direction shown, at most 20 rows. **Nothing is lost by the cap:** every set
in the cell, capped-out or not, is in that contrast's own table,
`tables/by_contrast/<contrast>/gsea_msigdb.csv` for the eight MSigDB collections or
`gsea_custom.csv` for the six custom databases (verified: `WT_heat` × `Hallmark` has 50 rows in
the table and 50 in the master; `WT_heat` × `HSR_lens` has 2 and 2). The subtitle names that
path on every panel. This view is most useful when up- and
down-enriched sets are present simultaneously. In the `WT_heat` × `Hallmark` cell that is the
case, with the directions the other way round from what one might guess: the hypoxia,
glycolysis and mTORC1 sets sit in the Up facet (HYPOXIA +1.91, GLYCOLYSIS +1.61,
MTORC1_SIGNALING +1.39, all padj < 0.05) while both interferon-response sets sit in the Down
facet (INTERFERON_ALPHA_RESPONSE −1.90, INTERFERON_GAMMA_RESPONSE −1.77, padj 1.7e-04 for
both). Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot_facet` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv` + `03_results/objects/02_de_results.rds` + `03_results/objects/geneset_<src>_<DB>.rds` |

---

## figures/by_contrast/<contrast>/<DB>/barplot.png

NES bar chart for FDR-significant pathways only (padj < 0.05), top 20 by |NES|, rendered
via the RNAseq-toolkit `gsea_barplot()`. Bars extend left (negative NES, blue) or right
(positive NES, orange) from zero.

**The 16 cells with no FDR-significant set carry an explicit empty-state panel**, not a blank
page: a boxed statement naming the panel's selection rule and pointing at `dotplot`/`facet`,
which show the full ranking regardless of significance. Those cells are `WT_heat`/MitoXplorer,
`KO_heat`/MitoXplorer, `Interaction`/{TransportDB, MitoPathways, MitoXplorer,
HSR_lens, TCR_activation}, `Geno_at_39`/{TransportDB, HSR_lens, TCR_activation},
`Geno_at_37`/{TF_Targets, TCR_activation}, `Temp_main`/{MitoPathways,
MitoXplorer}, and `Geno_main`/{TF_Targets, TCR_activation}. An empty barplot is a result about
that cell's FDR, never a rendering failure.

**How to read:** x-axis = NES (capped at ±3.2); bars extend from 0 toward positive
(orange = enriched in numerator/hot/WT) or negative (blue = enriched in
denominator/cold/KO). y-axis = formatted pathway names sorted by NES descending (most
positive at top). Only padj < 0.05 sets are shown (at most top 20 by |NES|). Bar length
encodes effect size (NES magnitude); all bars are significant, so no outline-ring encoding
is needed here. This is the most compact high-signal view. A cell showing the boxed
empty-state statement instead of bars had no set survive FDR correction in that
(contrast × database) cell. Claim tier: L3. Sample-to-condition mapping confirmed against the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_barplot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv` + `03_results/objects/02_de_results.rds` + `03_results/objects/geneset_<src>_<DB>.rds` |

---

## Cross-contrast overview figures

`figures/_overview/` holds six panels, each saved as a `.pdf` (vector) and `.png` (raster)
pair. The one captioned immediately below is written by `12_gsea_viz.R`. The three
`gsea_pooled_overview*` panels are written by `12b_gsea_overview_pooled_viz.R` and the two
`hypoxia_*` panels by `28_hypoxia_focus_viz.R`; all five are captioned at the end of this
file.

## figures/_overview/gsea_asymmetry_panel.png

Hallmark GSEA NES for every set whose MSigDB name matches one of two
patterns, drawn across the contrasts WT_heat, KO_heat, Interaction,
Temp_main. Top block (pattern
INTERFERON|INNATE_IMMUNE|INFLAMMATORY|TNF|IL6|IL2|ALLOGRAFT): the two
interferon-response sets carry negative NES in both per-genotype heat
contrasts (WT_heat -1.90, KO_heat -2.52 for the alpha set; -1.77 and
-2.20 for the gamma set) and positive NES in the interaction term
(+3.14, padj 2.7e-24; +2.83, padj 5.0e-22), while the TNF, IL2/STAT5,
inflammatory and allograft sets are positive in both heat contrasts
with a non-significant interaction. Bottom block (pattern
HYPOXIA|GLYCOLYSIS|MTORC1|OXIDATIVE_PHOSPHORYLATION): HYPOXIA,
GLYCOLYSIS and MTORC1 are positive and FDR-significant in both heat
contrasts (HYPOXIA WT_heat +1.91 padj 4.2e-06, KO_heat +1.95 padj
2.5e-06) with a negative, non-significant interaction NES (HYPOXIA
-1.27, padj 0.17), and OXIDATIVE_PHOSPHORYLATION is negative and
non-significant throughout. A non-significant interaction NES reads
as no detectable cGAS-dependence at n=5; the 1-df interaction term is
the least-powered contrast in this design. An enrichment locates a
set's gene content in a ranking, and establishing that the program a
set's name invokes is present takes a separate composition test.

**How to read:** Two-block heatmap, one tile per (gene set x contrast). Rows are the
Hallmark sets whose MSigDB name matches the block's pattern: top
block =
'INTERFERON|INNATE_IMMUNE|INFLAMMATORY|TNF|IL6|IL2|ALLOGRAFT'; bottom
block = 'HYPOXIA|GLYCOLYSIS|MTORC1|OXIDATIVE_PHOSPHORYLATION'. Within
each block the top 12 sets by max |NES| across the plotted contrasts
are drawn, and that rule governs every absence here: a set can be
missing on a small |NES| while its adjusted p is among the smallest.
The cap does not bite here (Interferon and inflammatory sets = 7,
Hypoxia and metabolic sets = 4), so every Hallmark set matching
either pattern is on the panel and no absence needs explaining.

The two blocks group set names by pattern. Which regulator drives
either block is a separate question. Fill = NES: orange = NES > 0,
enriched in the contrast numerator, 39 °C or WT; blue = NES < 0,
enriched in the denominator, 37 °C or cGAS-KO; fill is squished at
+/-3.2. A star marks padj < 0.05. Row order within a block is by NES.

The Interaction column is the cGAS-dependence read-out. A positive
significant interaction NES is consistent with cGAS-dependent
induction. A non-significant one reads as no detectable
cGAS-dependence at n = 5, the 1-df interaction term being the
least-powered contrast in this design. An enrichment locates a set's
gene content in a ranking; establishing that the program the set's
name invokes is present takes a separate composition test. Claim
tier: L3. Sample-to-condition mapping confirmed against the owner's
sample sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `geom_tile` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv` |

## figures/by_contrast/<contrast>/Hallmark/*.png

Hallmark gene-set enrichment for each contrast, as four panels per
cell: dotplot, facet, barplot and running sum. Every NES and adjusted
p is read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from geneset_msigdb_Hallmark.rds,
so the panels carry the sweep's own numbers. Each set's direction is
read per contrast on its own. These panels rank sets within a cell;
which regulator drives a set is a separate question.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

The two rankings diverge, and the WT_heat Hallmark cell shows how.
HALLMARK_HYPOXIA sits 6th by |NES| (+1.91) and 4th by adjusted p
(4.2e-06), so it falls inside the dotplot's top 20 by adjusted p and
outside the running sum's top 5 by |NES|.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_Hallmark.rds}` |

## figures/by_contrast/<contrast>/KEGG/*.png

KEGG gene-set enrichment for each contrast, as four panels per cell:
dotplot, facet, barplot and running sum. Every NES and adjusted p is
read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from geneset_msigdb_KEGG.rds, so
the panels carry the sweep's own numbers. Each set's direction is
read per contrast on its own. These panels rank sets within a cell;
which regulator drives a set is a separate question.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_KEGG.rds}` |

## figures/by_contrast/<contrast>/Reactome/*.png

Reactome gene-set enrichment for each contrast, as four panels per
cell: dotplot, facet, barplot and running sum. Every NES and adjusted
p is read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from geneset_msigdb_Reactome.rds,
so the panels carry the sweep's own numbers. Each set's direction is
read per contrast on its own. These panels rank sets within a cell;
which regulator drives a set is a separate question.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_Reactome.rds}` |

## figures/by_contrast/<contrast>/WikiPathways/*.png

WikiPathways gene-set enrichment for each contrast, as four panels
per cell: dotplot, facet, barplot and running sum. Every NES and
adjusted p is read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from
geneset_msigdb_WikiPathways.rds, so the panels carry the sweep's own
numbers. Each set's direction is read per contrast on its own. These
panels rank sets within a cell; which regulator drives a set is a
separate question.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_WikiPathways.rds}` |

## figures/by_contrast/<contrast>/GO_BP/*.png

GO_BP gene-set enrichment for each contrast, as four panels per cell:
dotplot, facet, barplot and running sum. Every NES and adjusted p is
read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from geneset_msigdb_GO_BP.rds, so
the panels carry the sweep's own numbers. Each set's direction is
read per contrast on its own. These panels rank sets within a cell;
which regulator drives a set is a separate question.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_GO_BP.rds}` |

## figures/by_contrast/<contrast>/GO_MF/*.png

GO_MF gene-set enrichment for each contrast, as four panels per cell:
dotplot, facet, barplot and running sum. Every NES and adjusted p is
read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from geneset_msigdb_GO_MF.rds, so
the panels carry the sweep's own numbers. Each set's direction is
read per contrast on its own. These panels rank sets within a cell;
which regulator drives a set is a separate question.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_GO_MF.rds}` |

## figures/by_contrast/<contrast>/GO_CC/*.png

GO_CC gene-set enrichment for each contrast, as four panels per cell:
dotplot, facet, barplot and running sum. Every NES and adjusted p is
read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from geneset_msigdb_GO_CC.rds, so
the panels carry the sweep's own numbers. Each set's direction is
read per contrast on its own. These panels rank sets within a cell;
which regulator drives a set is a separate question.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_GO_CC.rds}` |

## figures/by_contrast/<contrast>/TF_Targets/*.png

TF_Targets gene-set enrichment for each contrast, as four panels per
cell: dotplot, facet, barplot and running sum. Every NES and adjusted
p is read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from
geneset_msigdb_TF_Targets.rds, so the panels carry the sweep's own
numbers. Each set's direction is read per contrast on its own. These
panels rank sets within a cell; which regulator drives a set is a
separate question.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_TF_Targets.rds}` |

## figures/by_contrast/<contrast>/TransportDB/*.png

TransportDB gene-set enrichment for each contrast, as four panels per
cell: dotplot, facet, barplot and running sum. Every NES and adjusted
p is read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from
geneset_custom_TransportDB.rds, so the panels carry the sweep's own
numbers. Membrane-transporter sets. This database holds transporters
alone, so its rows bear on transport and carry no interferon or
heat-shock term.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_TransportDB.rds}` |

## figures/by_contrast/<contrast>/MitoPathways/*.png

MitoPathways gene-set enrichment for each contrast, as four panels
per cell: dotplot, facet, barplot and running sum. Every NES and
adjusted p is read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from
geneset_custom_MitoPathways.rds, so the panels carry the sweep's own
numbers. MitoCarta MitoPathways sets: mitochondrial function and
metabolism. The collection holds no interferon or heat-shock term.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_MitoPathways.rds}` |

## figures/by_contrast/<contrast>/MitoXplorer/*.png

MitoXplorer gene-set enrichment for each contrast, as four panels per
cell: dotplot, facet, barplot and running sum. Every NES and adjusted
p is read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from
geneset_custom_MitoXplorer.rds, so the panels carry the sweep's own
numbers. MitoXplorer sets, the mitochondrial metabolic atlas. The
collection holds no interferon or heat-shock term.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_MitoXplorer.rds}` |

## figures/by_contrast/<contrast>/HSR_lens/*.png

HSR_lens gene-set enrichment for each contrast, as four panels per
cell: dotplot, facet, barplot and running sum. Every NES and adjusted
p is read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from geneset_custom_HSR_lens.rds,
so the panels carry the sweep's own numbers. Two curated
heat-shock-response terms, HSR_core and HSR_sensitivity, built by the
00d and 00e curation steps and derived independently of the empirical
heat arms. Their scope is proteotoxic stress in general; the 37/39 °C
contrast is the one that bears on temperature.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_HSR_lens.rds}` |

## figures/by_contrast/<contrast>/TCR_activation/*.png

TCR_activation gene-set enrichment for each contrast, as four panels
per cell: dotplot, facet, barplot and running sum. Every NES and
adjusted p is read from master_gsea_table.csv, the ranked vector from
02_de_results.rds and the gene sets from
geneset_custom_TCR_activation.rds, so the panels carry the sweep's
own numbers. One curated TCR and immediate-early T-cell-activation
term, built by the 00f curation step with the HSR core's genes
excluded. A heat arm overlapping this term reads as activation.

**How to read:** Each panel ranks the sets it draws, and the rules differ, so an
absence is a statement about that panel's own ranking metric. Dotplot
and facet take the top 20 by adjusted p. Barplot takes the
FDR-significant sets, top 20 by |NES|. The running sum takes the top
5 by |NES|. The dotplot selects by adjusted p and orders its y axis
by GeneRatio descending, so vertical position there is a gene-ratio
order.

Dotplot: x = GeneRatio (leading-edge genes / set size), point size =
-log10(adjusted p), fill = NES (orange #B35806 positive, blue #2166AC
negative), black outline = adjusted p < 0.05. Facet splits that
dotplot into an NES > 0 and an NES < 0 panel. Barplot: NES bars from
zero, y ordered by NES descending.

Running sum: three stacked panels on one x axis. Top, the running
enrichment score, whose peak is the enrichment score and whose set
members left of the peak are the leading edge; its y range is pinned
to [-1.0, 1.0] so a curve here compares with every running sum in the
project. Middle, one named tick row per set in that set's own colour,
each tick a member at its rank. Bottom, the ranked moderated t, the
statistic the ranking was built on. X is a gene's position as a
fraction of the list's length, because ranked universes differ in
length between compartments, and the axis title carries this one's
size.

NES > 0 = enriched in the contrast numerator, 39 °C or WT. NES < 0 =
enriched in the denominator, 37 °C or cGAS-KO. An enrichment locates
a set's gene content in a ranking; establishing that the program the
set's name invokes is present takes a separate composition test.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_TCR_activation.rds}` |

## README overview

Gene-set enrichment of the iTreg 39/37 °C x cGAS design across 7
contrasts, over 8 MSigDB collections (Hallmark, KEGG, Reactome,
WikiPathways, GO_BP, GO_MF, GO_CC, TF_Targets) and 5 curated
databases (TransportDB, MitoPathways, MitoXplorer, HSR_lens,
TCR_activation).

Each (contrast x database) cell has its own dotplot, facet, barplot
and running sum under by_contrast/<contrast>/<database>/;
cross-contrast panels sit in _overview/. The enrichment itself was
computed by 05_gsea_msigdb_run.R and 06_gsea_custom_run.R; this stage
draws it. Every gene-set size is read at run time from
master_gsea_table.csv, so a size on a figure is the size the test
used. Claim tier: L3. Sample-to-condition mapping confirmed against
the owner's sample sheet (2026-07-22): 20 of 20 libraries concordant
with the label-blind marker call.

**How to read:** NES > 0 = enriched in the numerator side of the contrast, 39 °C or
WT. NES < 0 = enriched in the denominator side, 37 °C or cGAS-KO.
Orange = positive, blue = negative (colors.diverging in the config).
A star or a filled point marks adjusted p < 0.05.

The Interaction contrast is (WT_heat) - (KO_heat), and it is the
cGAS-dependence read-out. A positive significant Interaction NES is
consistent with cGAS-dependent induction. A non-significant one reads
as no detectable cGAS-dependence at n = 5; the 1-df interaction term
is the least-powered contrast in this design.

Every panel here ranks the sets it draws, on adjusted p for the
dotplot and facet and on |NES| for the barplot and running sum. Read
an absence against the ranking rule in that panel's own caption.
Claim tier: L3, DE and enrichment statistics.

Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call. Group identity was recovered independently
from marker genes as well, the Hspa1b thermometer plus Cgas, and the
two agree.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `save_overview / save_figure` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; running_sum_x=rank/n_ranked; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv` |

## figures/_overview/gsea_pooled_overview_WT_heat.png

Cross-database GSEA overview for the WT_heat contrast, labelled WT:
heat (39 vs 37 °C), under the pooled Benjamini-Hochberg family of
4,792 tests across 13 databases. Per-database FDR < 0.05 gives 860
significant sets; the pooled FDR gives 865 (46 lost on pooling, 51
gained). The panel draws the top 5 sets per database for 13 of the 13
pooled databases. MitoXplorer carries no pooled-significant set at
all and keeps its facet so the null stays visible.

**How to read:** Lollipop panel, one facet per database, x = normalised enrichment
score (NES). Orange fill = positive NES (enriched in the contrast
numerator), blue = negative, squished at the configured NES cap. The
glyph says what becomes of a set when the per-database
Benjamini-Hochberg correction is replaced by one pooled correction
over every database: filled circle = significant under both, diamond
= per-database only (lost on pooling), triangle = pooled only (gained
on pooling), open circle = neither. The legend keys carry a neutral
grey fill because the panel spends fill on NES.

Each facet draws the top N sets in its database by padj_pooled, ties
broken on the raw p-value and then on |NES|, with N from
gsea_pooled_overview.top_n_per_db. That rank is the whole of the
selection. An earlier revision also pinned the Hallmark HYPOXIA and
INTERFERON sets and the three project-curated lenses at any rank,
which let the panel argue for a conclusion; the pin has been removed.

A facet header reads (n/N): n of the N sets in that database pass the
pooled FDR, counted over the whole database. A header reading (0/23)
marks a real null, and those facets keep their place on the canvas;
dropping a database for coming out empty would leave a panel selected
for positives. The panel omits whatever
gsea_pooled_overview.exclude_databases names, currently nothing, so
the facet count and the pooled-family database count agree. An
omitted database stays inside the pooled correction and in every
06_gsea table.

Read padj_pooled as a comparability device across databases. It is a
calibrated error rate only where the sets are independent, and GO
terms and pathway sets share genes. Claim tier: L3, enrichment
statistics. Sample mapping owner-confirmed.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12b_gsea_overview_pooled_viz.R` | `ggplot2 / geom_point + facet_wrap` | `thresholds.gsea_fdr=0.05; gsea_pooled_overview.top_n_per_db=5; gsea_pooled_overview.exclude_databases=[]; gsea_pooled_overview.facet_ncol=2; figures.nes_cap=3.2` | `03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv` |

## figures/_overview/gsea_pooled_overview_Interaction.png

Cross-database GSEA overview for the Interaction contrast, labelled
Heat × genotype interaction (cGAS-dependence of the heat response),
under the pooled Benjamini-Hochberg family of 4,792 tests across 13
databases. Per-database FDR < 0.05 gives 192 significant sets; the
pooled FDR gives 164 (32 lost on pooling, 4 gained). The panel draws
the top 5 sets per database for 13 of the 13 pooled databases. 5 of
them carry no pooled-significant set at all (HSR_lens, MitoPathways,
MitoXplorer, TCR_activation, TransportDB) and keep their facets so
the nulls stay visible.

**How to read:** Lollipop panel, one facet per database, x = normalised enrichment
score (NES). Orange fill = positive NES (enriched in the contrast
numerator), blue = negative, squished at the configured NES cap. The
glyph says what becomes of a set when the per-database
Benjamini-Hochberg correction is replaced by one pooled correction
over every database: filled circle = significant under both, diamond
= per-database only (lost on pooling), triangle = pooled only (gained
on pooling), open circle = neither. The legend keys carry a neutral
grey fill because the panel spends fill on NES.

Each facet draws the top N sets in its database by padj_pooled, ties
broken on the raw p-value and then on |NES|, with N from
gsea_pooled_overview.top_n_per_db. That rank is the whole of the
selection. An earlier revision also pinned the Hallmark HYPOXIA and
INTERFERON sets and the three project-curated lenses at any rank,
which let the panel argue for a conclusion; the pin has been removed.

A facet header reads (n/N): n of the N sets in that database pass the
pooled FDR, counted over the whole database. A header reading (0/23)
marks a real null, and those facets keep their place on the canvas;
dropping a database for coming out empty would leave a panel selected
for positives. The panel omits whatever
gsea_pooled_overview.exclude_databases names, currently nothing, so
the facet count and the pooled-family database count agree. An
omitted database stays inside the pooled correction and in every
06_gsea table.

Read padj_pooled as a comparability device across databases. It is a
calibrated error rate only where the sets are independent, and GO
terms and pathway sets share genes. Claim tier: L3, enrichment
statistics. Sample mapping owner-confirmed.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12b_gsea_overview_pooled_viz.R` | `ggplot2 / geom_point + facet_wrap` | `thresholds.gsea_fdr=0.05; gsea_pooled_overview.top_n_per_db=5; gsea_pooled_overview.exclude_databases=[]; gsea_pooled_overview.facet_ncol=2; figures.nes_cap=3.2` | `03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv` |

## figures/_overview/gsea_pooled_overview.png

The two per-contrast pooled-FDR overviews abreast on one canvas, WT
heat on the left and Interaction on the right, so the 13 displayed
databases line up block for block between the contrasts. Under the
shared pooled family of 4,792 tests across 13 databases, WT heat
carries 860 significant sets before pooling and 865 after,
Interaction 192 and 164. The per-contrast panels hold the full
numbers.

**How to read:** Two copies of the per-contrast panel abreast: WT heat on the left,
Interaction on the right, sharing the single legend at the foot. The
facet order is identical on both sides, so a database sits at the
same height in each; each panel's own subtitle carries its
pooled-family counts, which differ between contrasts. Lollipop panel,
one facet per database, x = normalised enrichment score (NES). Orange
fill = positive NES (enriched in the contrast numerator), blue =
negative, squished at the configured NES cap. The glyph says what
becomes of a set when the per-database Benjamini-Hochberg correction
is replaced by one pooled correction over every database: filled
circle = significant under both, diamond = per-database only (lost on
pooling), triangle = pooled only (gained on pooling), open circle =
neither. The legend keys carry a neutral grey fill because the panel
spends fill on NES.

Each facet draws the top N sets in its database by padj_pooled, ties
broken on the raw p-value and then on |NES|, with N from
gsea_pooled_overview.top_n_per_db. That rank is the whole of the
selection. An earlier revision also pinned the Hallmark HYPOXIA and
INTERFERON sets and the three project-curated lenses at any rank,
which let the panel argue for a conclusion; the pin has been removed.

A facet header reads (n/N): n of the N sets in that database pass the
pooled FDR, counted over the whole database. A header reading (0/23)
marks a real null, and those facets keep their place on the canvas;
dropping a database for coming out empty would leave a panel selected
for positives. The panel omits whatever
gsea_pooled_overview.exclude_databases names, currently nothing, so
the facet count and the pooled-family database count agree. An
omitted database stays inside the pooled correction and in every
06_gsea table.

Read padj_pooled as a comparability device across databases. It is a
calibrated error rate only where the sets are independent, and GO
terms and pathway sets share genes. Claim tier: L3, enrichment
statistics. Sample mapping owner-confirmed. This canvas is double
width by construction (two full facet grids at N sets per database);
for reading one contrast, prefer gsea_pooled_overview_WT_heat.png or
gsea_pooled_overview_Interaction.png.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12b_gsea_overview_pooled_viz.R` | `patchwork / p_wt beside p_int` | `thresholds.gsea_fdr=0.05; gsea_pooled_overview.top_n_per_db=5; gsea_pooled_overview.exclude_databases=[]` | `03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv` |

## figures/_overview/hypoxia_routes_by_contrast.png

Running-enrichment curves for four hypoxia-named gene sets, one per
database, across the three focal contrasts. WT: heat (39 vs 37 °C):
Hallmark NES +1.92 (padj 5.32e-06), GO MF NES +1.75 (padj 0.0374), GO
BP NES +1.73 (padj 5.96e-04), Reactome NES -0.91 (padj 0.723).
cGAS-KO: heat (39 vs 37 °C): Hallmark NES +1.96 (padj 3.01e-06), GO
MF NES +1.77 (padj 0.04), GO BP NES +1.56 (padj 0.00991), Reactome
NES -0.92 (padj 0.708). Heat × genotype interaction (cGAS-dependence
of the heat response): Hallmark NES -1.27 (padj 0.171), GO MF NES
+0.67 (padj 0.987), GO BP NES +1.34 (padj 0.306), Reactome NES -0.78
(padj 0.937). Three of the four sets carry a positive NES and reach
FDR < 0.05 in both per-genotype heat contrasts, and the fourth,
REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA, carries a negative NES and no
significance in any of the three; 6 of the 12 (set x contrast) cells
reach FDR < 0.05. Every Interaction NES is non-significant, which
reads as no detectable cGAS-dependence at n=5; the 1-df interaction
term is the least-powered contrast in this design. All four sets were
chosen by name, so the panel reports these four curves and leaves
ranking to the general panels.

**How to read:** One facet per contrast, and inside each facet four overlaid
running-enrichment curves, one per gene set, keyed by colour in the
shared legend below the panels. X is each gene's position in that
facet's own ranked list as a fraction of its length, most
numerator-shifted at 0 and most denominator-shifted at 1, so the
facets share an axis while each keeps its own ordering. The fraction
carries across ranked universes of different length, and the axis
title carries these rankings' size.

A curve steps up at each gene belonging to its set and decays between
them, so an early high peak places the members at the numerator end
and a curve near zero spreads them through the list. The y range is
pinned to [-1, 1] so every curve here compares with every other
running sum in this GSEA sweep. Inside each facet the NES and
adjusted p for that contrast are printed in each set's own colour, in
the legend's top-to-bottom order and keyed by database name, each of
the four sets coming from a different database. Four curves per facet
stay legible with the member ticks and the ranked-metric panel left
off; the companion three-panel figure draws them for one set. NES > 0
= enriched in the contrast numerator (39 °C or WT). NES < 0 =
enriched in the denominator (37 °C or cGAS-KO).

The Interaction facet is the cGAS-dependence read-out. A positive
significant Interaction NES is consistent with cGAS-dependent
induction; a non-significant one reads as no detectable
cGAS-dependence at n=5. These 4 sets were chosen by name, one per
database. Their adjusted p and |NES| played no part in the choice,
and the sweep's general panels select on those statistics alone. This
is a named-set read-out; where these four rank among the thousands of
sets they came from is a separate question.

Read the fourth curve as closely as the other three. Three of these
hypoxia-named sets carry a positive NES in both heat contrasts and
the Reactome set carries a negative one, so shared wording in two set
names leaves shared behaviour an open question. An enrichment locates
a set's gene content in a ranking. Establishing that the program the
set's name invokes is present takes a separate composition test. The
four names overlap in wording while their gene content differs, which
is why each has its own curve. Claim tier: L3 (DE and enrichment
statistics). Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/28_hypoxia_focus_viz.R` | `geom_line / facet_wrap` | `thresholds.gsea_fdr=0.05; figures.running_sum_ylim=[-1.0,1.0]; figures.running_sum_heights; running_sum_x=rank/n_ranked; figures.running_sum_top=5; figures.top_pathways=20; colors.okabe_ito` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_{Hallmark,GO_MF,GO_BP,Reactome}.rds}` |

## figures/_overview/runsum_HALLMARK_HYPOXIA.png

Running enrichment of the single Hallmark set HALLMARK_HYPOXIA
against each of the three focal contrasts' ranked t-statistics, on
one axis. In WT heat: NES +1.9240, p 4.256e-07, adjusted p 5.32e-06,
172 of the set's genes present in the 19,679-gene ranked universe, 57
of them in the leading edge, and the curve peaks at an enrichment
score of +0.458 in the top 12% of the ranking, so the set's members
are concentrated toward the 39 °C end. Across the three contrasts the
NES values are WT heat (39 vs 37 °C) +1.92 (FDR 5.32e-06), cGAS-KO
heat (39 vs 37 °C) +1.96 (FDR 3.01e-06), Interaction -1.27 (FDR
0.171), so the set separates warmed from unwarmed iTregs in both
genotypes while the interaction term carries none. The statistics are
taken verbatim from master_gsea_table.csv, and every curve is
recomputed deterministically from the same ranked vector the sweep
used, at exponent 1, with no permutation re-run. This set sits 6th by
|NES| and 4th by adjusted p in the WT_heat Hallmark cell, so the
|NES| quota of this sweep's general running-sum panels leaves it out.
Its enrichment locates where these 172 genes sit in a ranking, and
the gene content of the set is a separate question from the name the
set carries.

**How to read:** Three stacked panels sharing one x axis, which is each gene's
position in that contrast's ranked list as a FRACTION of its length,
most 39 °C-shifted at 0, most 37 °C-shifted at 1, because ranked
universes differ in length between compartments; the axis title
carries this one's size. Colour keys the CONTRAST, and the same three
colours run through all three panels. Top panel: the running
enrichment score, which steps up at each gene belonging to the set
and decays between them. Its peak is the enrichment score, and the
set members left of the peak are the leading edge. The y range is
pinned to [-1, 1] so the curves stay comparable to every other
running-sum figure in this GSEA sweep. Middle panel: one tick per set
member at that member's rank, in its own labelled row per contrast,
so where a contrast places the same genes can be read row against
row. Bottom panel: the ranked t-statistic each curve above was
computed on, which shows how much signal each rank carries; it is
drawn at every 25th rank plus both endpoints, which redraws a sorted
vector exactly at this size. The legend carries each contrast's NES,
adjusted p and genes present in the ranked universe. NES > 0 =
enriched in the contrast numerator (39 °C or WT). NES < 0 = enriched
in the denominator (37 °C or cGAS-KO). This set is 6th by |NES| and
4th by adjusted p in its WT_heat cell, so the general running-sum
panels of this GSEA sweep, which draw the top 5 per cell by |NES|,
leave it out. This set was chosen by name, and its adjusted p and
|NES| played no part in the choice. The Hallmark collection
contributed 50 sets to this sweep, and the two ranks quoted above are
the whole of what this panel says about where this one stands among
them. An enrichment locates a set's gene content in a ranking.
Establishing that the program the set's name invokes is present takes
a separate composition test. Claim tier: L3 (DE and enrichment
statistics). Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/28_hypoxia_focus_viz.R` | `top-level (f1_p_es / f1_p_rug / f1_p_met)` | `thresholds.gsea_fdr=0.05; figures.running_sum_ylim=[-1.0,1.0]; figures.running_sum_heights; running_sum_x=rank/n_ranked; figures.running_sum_top=5; figures.top_pathways=20; colors.okabe_ito` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_Hallmark.rds}` |

