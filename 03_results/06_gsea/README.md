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
  (Pdk1/Bnip3/Bnip3l/Car9) is repressed. That list is not plotted in this stage; see
  `03_results/04_tf/` for the comparison against it.
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
numerator side. Claim tier: L3. PROVISIONAL sample labels.

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
row order. Claim tier: L3. PROVISIONAL sample labels.

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
both). Claim tier: L3. PROVISIONAL sample labels.

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
(contrast × database) cell. Claim tier: L3. PROVISIONAL sample labels.

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
as no detectable cGAS-dependence at n=5 and never as proven
cGAS-independence: the 1-df interaction term is the least-powered
contrast in this design. Set composition is not established by this
panel, so a set's enrichment is not evidence that the program its
name invokes is present.

**How to read:** Two-block heatmap, one tile per (gene set x contrast). Rows are the
Hallmark sets whose MSigDB name matches the block's pattern: top
block =
'INTERFERON|INNATE_IMMUNE|INFLAMMATORY|TNF|IL6|IL2|ALLOGRAFT'; bottom
block = 'HYPOXIA|GLYCOLYSIS|MTORC1|OXIDATIVE_PHOSPHORYLATION'.
SELECTION RULE, which governs every absence you see: within each
block only the top 12 sets by max |NES| across the plotted contrasts
are drawn, so a set can be missing because its |NES| is small even
when its adjusted p is among the smallest. The cap does not bite here
(Interferon and inflammatory sets = 7, Hypoxia and metabolic sets =
4), so every Hallmark set matching either pattern is on the panel and
no absence needs explaining. The two blocks are NAME-PATTERN
groupings over set names and carry no claim about which regulator
drives either block. Fill = NES: orange = NES > 0 (enriched in the
contrast numerator, 39 °C or WT); blue = NES < 0 (enriched in the
denominator, 37 °C or cGAS-KO); fill is squished at +/-3.2. * = padj
< 0.05. Row order within a block is by NES. The Interaction column is
the cGAS-dependence read-out: a positive significant interaction NES
is consistent with cGAS-dependent induction; a non-significant
interaction NES means no detectable cGAS-dependence at n=5 and never
'cGAS-independent'. A set enriching is not evidence that the program
its name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `geom_tile` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv` |

## figures/by_contrast/<contrast>/Hallmark/*.png

Hallmark GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_Hallmark.rds).
Read each set's direction per contrast on its own; the panels rank
sets within a cell and make no statement about which regulator drives
any of them.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). Worked example of the two rankings diverging: in the
WT_heat Hallmark cell HALLMARK_HYPOXIA is 6th by |NES| (+1.91) and
4th by adjusted p (4.2e-06), so it is outside the running-sum's top 5
by |NES| while sitting well inside the dotplot's top 20 by adjusted
p. An absence from one panel is a statement about that panel's
ranking metric only. A set's enrichment is not evidence that the
program its name invokes is present; composition would have to be
established separately. Claim tier: L3. Sample-to-condition mapping
confirmed against the owner's sample sheet (2026-07-22): 20 of 20
libraries concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_Hallmark.rds}` |

## figures/by_contrast/<contrast>/KEGG/*.png

KEGG GSEA per contrast (dotplot/facet/barplot/running_sum) via the
RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_KEGG.rds). Read
each set's direction per contrast on its own; the panels rank sets
within a cell and make no statement about which regulator drives any
of them.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_KEGG.rds}` |

## figures/by_contrast/<contrast>/Reactome/*.png

Reactome GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_Reactome.rds).
Read each set's direction per contrast on its own; the panels rank
sets within a cell and make no statement about which regulator drives
any of them.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_Reactome.rds}` |

## figures/by_contrast/<contrast>/WikiPathways/*.png

WikiPathways GSEA per contrast (dotplot/facet/barplot/running_sum)
via the RNAseq-toolkit plotters on a viz-side-reconstructed
gseaResult (NES/padj taken verbatim from master_gsea_table.csv;
ranked vector from 02_de_results.rds; gene sets from
geneset_msigdb_WikiPathways.rds). Read each set's direction per
contrast on its own; the panels rank sets within a cell and make no
statement about which regulator drives any of them.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_WikiPathways.rds}` |

## figures/by_contrast/<contrast>/GO_BP/*.png

GO_BP GSEA per contrast (dotplot/facet/barplot/running_sum) via the
RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_GO_BP.rds).
Read each set's direction per contrast on its own; the panels rank
sets within a cell and make no statement about which regulator drives
any of them.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_GO_BP.rds}` |

## figures/by_contrast/<contrast>/GO_MF/*.png

GO_MF GSEA per contrast (dotplot/facet/barplot/running_sum) via the
RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_GO_MF.rds).
Read each set's direction per contrast on its own; the panels rank
sets within a cell and make no statement about which regulator drives
any of them.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_GO_MF.rds}` |

## figures/by_contrast/<contrast>/GO_CC/*.png

GO_CC GSEA per contrast (dotplot/facet/barplot/running_sum) via the
RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_msigdb_GO_CC.rds).
Read each set's direction per contrast on its own; the panels rank
sets within a cell and make no statement about which regulator drives
any of them.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_GO_CC.rds}` |

## figures/by_contrast/<contrast>/TF_Targets/*.png

TF_Targets GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from
geneset_msigdb_TF_Targets.rds). Read each set's direction per
contrast on its own; the panels rank sets within a cell and make no
statement about which regulator drives any of them.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_TF_Targets.rds}` |

## figures/by_contrast/<contrast>/TransportDB/*.png

TransportDB GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from
geneset_custom_TransportDB.rds). Membrane-transporter sets; this
database carries no interferon or heat-shock sets, so nothing here
bears on those axes.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_TransportDB.rds}` |

## figures/by_contrast/<contrast>/MitoPathways/*.png

MitoPathways GSEA per contrast (dotplot/facet/barplot/running_sum)
via the RNAseq-toolkit plotters on a viz-side-reconstructed
gseaResult (NES/padj taken verbatim from master_gsea_table.csv;
ranked vector from 02_de_results.rds; gene sets from
geneset_custom_MitoPathways.rds). MitoCarta MitoPathways sets
(mitochondrial function and metabolism); this database carries no
interferon or heat-shock sets.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_MitoPathways.rds}` |

## figures/by_contrast/<contrast>/MitoXplorer/*.png

MitoXplorer GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from
geneset_custom_MitoXplorer.rds). MitoXplorer sets (mitochondrial
metabolic atlas); this database carries no interferon or heat-shock
sets.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_MitoXplorer.rds}` |

## figures/by_contrast/<contrast>/HSR_lens/*.png

HSR_lens GSEA per contrast (dotplot/facet/barplot/running_sum) via
the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult
(NES/padj taken verbatim from master_gsea_table.csv; ranked vector
from 02_de_results.rds; gene sets from geneset_custom_HSR_lens.rds).
Two curated heat-shock-response terms (HSR_core, HSR_sensitivity),
built by 00d/00e and held anchor-independent of the empirical heat
arms. Proteotoxic-stress-general in scope; only the 37/39 °C contrast
bears on temperature.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_HSR_lens.rds}` |

## figures/by_contrast/<contrast>/TCR_activation/*.png

TCR_activation GSEA per contrast (dotplot/facet/barplot/running_sum)
via the RNAseq-toolkit plotters on a viz-side-reconstructed
gseaResult (NES/padj taken verbatim from master_gsea_table.csv;
ranked vector from 02_de_results.rds; gene sets from
geneset_custom_TCR_activation.rds). One curated TCR and
immediate-early T-cell-activation term, built by 00f and disjoint
from the HSR core by construction. Overlap of a heat arm with this
term reads as activation.

**How to read:** SELECTION RULES, one per panel, which govern which sets appear:
dotplot and facet draw the top 20 sets by ADJUSTED P; barplot draws
FDR-significant sets only, top 20 by |NES|; running_sum draws the top
5 sets by |NES|. Note that the dotplot SELECTS by adjusted p but
ORDERS its y-axis by GeneRatio descending, so vertical position on
the dotplot is not a significance ranking. dotplot: x = GeneRatio
(leading-edge genes / set size), point size = -log10(padj), fill =
NES (orange #B35806 = positive / blue #2166AC = negative), black
outline = padj < 0.05. facet: the same dotplot split into an NES>0
and an NES<0 panel. barplot: NES bars from zero, y-axis ordered by
NES descending. running_sum: a three-panel enrichment curve per set
(top = running enrichment score with the leading-edge peak; middle =
gene-hit ticks at each member's rank; bottom = the ranked
t-statistic), ES y-range pinned to [-1.0,1.0] so curves stay
comparable across databases. NES > 0 = enriched in the contrast
numerator (39 °C or WT); NES < 0 = enriched in the denominator (37 °C
or cGAS-KO). A set's enrichment is not evidence that the program its
name invokes is present; composition would have to be established
separately. Claim tier: L3. Sample-to-condition mapping confirmed
against the owner's sample sheet (2026-07-22): 20 of 20 libraries
concordant with the label-blind marker call.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_custom_TCR_activation.rds}` |

## README overview

06_gsea GSEA stage: 8 MSigDB collections (Hallmark, KEGG, Reactome,
WikiPathways, GO_BP, GO_MF, GO_CC, TF_Targets) plus 5 custom
databases (TransportDB, MitoPathways, MitoXplorer, HSR_lens,
TCR_activation) across 7 contrasts. Per-contrast figures in
by_contrast/<c>/<DB>/ are built with the RNAseq-toolkit plotters
(gsea_dotplot/facet/barplot/running_sum) on a gseaResult
reconstructed viz-side from the master table plus cached DE ranks
plus gene-set lists; the running_sum is a real enrichplot::gseaplot2
three-panel ES curve. Cross-contrast overviews in _overview/.
Produced by 12_gsea_viz.R (VIZ-ONLY; GSEA computed by
05/06_gsea_*_run.R). Standing constraints: write 'no detectable
cGAS-dependence at n=5' and never 'cGAS-independent'; do not name a
gene set or a block of sets after a regulator it is hoped to
represent; read every gene-set size at runtime from the master table.
Claim tier: L3. Sample-to-condition mapping confirmed against the
owner's sample sheet (2026-07-22): 20 of 20 libraries concordant with
the label-blind marker call.

**How to read:** NES > 0 = enriched in the numerator side of the contrast (39 °C or WT
genotype). NES < 0 = enriched in the denominator side (37 °C or
cGAS-KO). Glyph convention: orange = positive, blue = negative
(colors.diverging in config). * or filled point = padj < 0.05.
Interaction contrast: (WT_heat) - (KO_heat); a positive Interaction
NES is consistent with cGAS-dependent induction, and a
non-significant one means no detectable cGAS-dependence at n=5. Every
panel in this stage ranks the sets it draws, and the ranking metric
differs between panels (adjusted p for dotplot/facet, |NES| for
barplot/running_sum), so read an absence against the ranking rule
named in that panel's own caption before reading it as a null.
Sample-to-condition mapping confirmed against the owner's sample
sheet (2026-07-22): 20 of 20 libraries concordant with the
label-blind marker call. Group identity was also recovered
independently from marker genes (the Hspa1b thermometer plus Cgas),
and the two agree. Claim tier: L3 (DE/enrichment statistics). Never
L7 (mechanism) in a figure title.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12_gsea_viz.R` | `save_overview / save_figure` | `thresholds.gsea_fdr=0.05; figures.top_pathways=20; figures.running_sum_top=5; figures.running_sum_ylim=[-1.0,1.0]; figures.nes_cap=3.2; colors.diverging` | `03_results/master/master_gsea_table.csv` |

## figures/_overview/gsea_pooled_overview_WT_heat.png

Cross-database GSEA overview for the WT_heat contrast, labelled WT:
heat (39 vs 37 °C), under the pooled Benjamini-Hochberg family of
4,786 tests across 13 databases. Per-database FDR < 0.05 gives 854
significant sets; the pooled FDR gives 853 (44 lost on pooling, 43
gained). The panel draws the top 5 sets per database for 13 of the 13
pooled databases. MitoXplorer carries no pooled-significant set at
all and keeps its facet so the null stays visible.

**How to read:** Lollipop panel, one facet per database, x = Normalized Enrichment
Score (NES). Orange fill = positive NES (enriched in the contrast
numerator), blue = negative, squished at the configured NES cap.
GLYPH = what becomes of a set when the per-database
Benjamini-Hochberg correction is replaced by one pooled correction
over every database: filled circle = significant under both, diamond
= per-database only (lost on pooling), triangle = pooled only (gained
on pooling), open circle = neither. The legend keys carry a neutral
grey fill because the panel spends fill on NES. SELECTION RULE: top N
per database by padj_pooled, ties broken by the raw p-value then by
|NES|; N is gsea_pooled_overview.top_n_per_db. Nothing is pinned on
top of that rank. An earlier revision force-included the Hallmark
HYPOXIA and INTERFERON sets and the three project-curated lenses
whatever their rank, which let the panel argue for a conclusion; that
pin is gone. FACET HEADER (n/N): n of the N sets in that database
pass the pooled FDR, counted over the whole database. A header
reading (0/23) marks a real null, and those facets stay on the
canvas, since dropping a database because it came out empty would
leave a panel selected for positives. EXCLUDED FROM THE PANEL:
whatever gsea_pooled_overview.exclude_databases names, which is
currently empty, so the facet count and the pooled-family database
count agree. An excluded database stays inside the pooled correction
and in every 06_gsea table. CAVEAT: padj_pooled is a comparability
device across databases rather than a calibrated error rate, since GO
terms and pathway sets share genes. Claim tier: L3 (enrichment
statistics). Sample mapping owner-confirmed.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12b_gsea_overview_pooled_viz.R` | `ggplot2 / geom_point + facet_wrap` | `thresholds.gsea_fdr=0.05; gsea_pooled_overview.top_n_per_db=5; gsea_pooled_overview.exclude_databases=[]; gsea_pooled_overview.facet_ncol=2; figures.nes_cap=3.2` | `03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv` |

## figures/_overview/gsea_pooled_overview_Interaction.png

Cross-database GSEA overview for the Interaction contrast, labelled
Heat × genotype interaction (cGAS-dependence of the heat response),
under the pooled Benjamini-Hochberg family of 4,786 tests across 13
databases. Per-database FDR < 0.05 gives 186 significant sets; the
pooled FDR gives 154 (36 lost on pooling, 4 gained). The panel draws
the top 5 sets per database for 13 of the 13 pooled databases. 5 of
them carry no pooled-significant set at all (HSR_lens, MitoPathways,
MitoXplorer, TCR_activation, TransportDB) and keep their facets so
the nulls stay visible.

**How to read:** Lollipop panel, one facet per database, x = Normalized Enrichment
Score (NES). Orange fill = positive NES (enriched in the contrast
numerator), blue = negative, squished at the configured NES cap.
GLYPH = what becomes of a set when the per-database
Benjamini-Hochberg correction is replaced by one pooled correction
over every database: filled circle = significant under both, diamond
= per-database only (lost on pooling), triangle = pooled only (gained
on pooling), open circle = neither. The legend keys carry a neutral
grey fill because the panel spends fill on NES. SELECTION RULE: top N
per database by padj_pooled, ties broken by the raw p-value then by
|NES|; N is gsea_pooled_overview.top_n_per_db. Nothing is pinned on
top of that rank. An earlier revision force-included the Hallmark
HYPOXIA and INTERFERON sets and the three project-curated lenses
whatever their rank, which let the panel argue for a conclusion; that
pin is gone. FACET HEADER (n/N): n of the N sets in that database
pass the pooled FDR, counted over the whole database. A header
reading (0/23) marks a real null, and those facets stay on the
canvas, since dropping a database because it came out empty would
leave a panel selected for positives. EXCLUDED FROM THE PANEL:
whatever gsea_pooled_overview.exclude_databases names, which is
currently empty, so the facet count and the pooled-family database
count agree. An excluded database stays inside the pooled correction
and in every 06_gsea table. CAVEAT: padj_pooled is a comparability
device across databases rather than a calibrated error rate, since GO
terms and pathway sets share genes. Claim tier: L3 (enrichment
statistics). Sample mapping owner-confirmed.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12b_gsea_overview_pooled_viz.R` | `ggplot2 / geom_point + facet_wrap` | `thresholds.gsea_fdr=0.05; gsea_pooled_overview.top_n_per_db=5; gsea_pooled_overview.exclude_databases=[]; gsea_pooled_overview.facet_ncol=2; figures.nes_cap=3.2` | `03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv` |

## figures/_overview/gsea_pooled_overview.png

The two per-contrast pooled-FDR overviews abreast on one canvas, WT
heat on the left and Interaction on the right, so the 13 displayed
databases line up block for block between the contrasts. Under the
shared pooled family of 4,786 tests across 13 databases, WT heat
carries 854 significant sets before pooling and 853 after,
Interaction 186 and 154. The per-contrast panels hold the full
numbers.

**How to read:** Two copies of the per-contrast panel abreast: WT heat on the left,
Interaction on the right, sharing the single legend at the foot. The
facet order is identical on both sides, so a database sits at the
same height in each; each panel's own subtitle carries its
pooled-family counts, which differ between contrasts. Lollipop panel,
one facet per database, x = Normalized Enrichment Score (NES). Orange
fill = positive NES (enriched in the contrast numerator), blue =
negative, squished at the configured NES cap. GLYPH = what becomes of
a set when the per-database Benjamini-Hochberg correction is replaced
by one pooled correction over every database: filled circle =
significant under both, diamond = per-database only (lost on
pooling), triangle = pooled only (gained on pooling), open circle =
neither. The legend keys carry a neutral grey fill because the panel
spends fill on NES. SELECTION RULE: top N per database by
padj_pooled, ties broken by the raw p-value then by |NES|; N is
gsea_pooled_overview.top_n_per_db. Nothing is pinned on top of that
rank. An earlier revision force-included the Hallmark HYPOXIA and
INTERFERON sets and the three project-curated lenses whatever their
rank, which let the panel argue for a conclusion; that pin is gone.
FACET HEADER (n/N): n of the N sets in that database pass the pooled
FDR, counted over the whole database. A header reading (0/23) marks a
real null, and those facets stay on the canvas, since dropping a
database because it came out empty would leave a panel selected for
positives. EXCLUDED FROM THE PANEL: whatever
gsea_pooled_overview.exclude_databases names, which is currently
empty, so the facet count and the pooled-family database count agree.
An excluded database stays inside the pooled correction and in every
06_gsea table. CAVEAT: padj_pooled is a comparability device across
databases rather than a calibrated error rate, since GO terms and
pathway sets share genes. Claim tier: L3 (enrichment statistics).
Sample mapping owner-confirmed. This canvas is double width by
construction (two full facet grids at N sets per database); for
reading one contrast, prefer gsea_pooled_overview_WT_heat.png or
gsea_pooled_overview_Interaction.png.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/12b_gsea_overview_pooled_viz.R` | `patchwork / p_wt beside p_int` | `thresholds.gsea_fdr=0.05; gsea_pooled_overview.top_n_per_db=5; gsea_pooled_overview.exclude_databases=[]` | `03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv` |


## figures/_overview/hypoxia_running_sum_wt_heat.png

Running enrichment of the single Hallmark set HALLMARK_HYPOXIA
against the WT 39-vs-37 °C ranked t-statistic. NES +1.9078, p
3.331e-07, adjusted p 4.164e-06, 171 of the set's genes present in
the 19,679-gene ranked universe, 56 of them in the leading edge. The
curve peaks at an enrichment score of +0.455 at rank 2,297, which is
the top 12% of the ranking, so the set's members are concentrated
toward the 39 °C end of the WT ordering. The statistics on the figure
are taken verbatim from master_gsea_table.csv, and the curve geometry
is recomputed deterministically from the same ranked vector the sweep
used, at exponent 1, with no permutation re-run. This set sits 6th by
|NES| and 4th by adjusted p in the WT_heat Hallmark cell, so it falls
outside the general running-sum panels of this GSEA sweep, which draw
the top 5 sets per cell by |NES|. Its enrichment locates where these
171 genes sit in one ranking, and the gene content of the set is a
separate question from the name the set carries.

**How to read:** Three stacked panels sharing one x axis: the rank of every gene in
the WT_heat ranked list, most 39 °C-shifted on the left, most 37
°C-shifted on the right. Top panel: the running enrichment score,
which steps up at each gene belonging to the set and decays between
them. Its peak is the enrichment score, and the set members left of
the peak are the leading edge. The y range is pinned to [-1, 1] so
the curve stays comparable to every other running-sum figure in this
GSEA sweep. Middle panel: one tick per set member at that member's
rank, over a band showing where the ranking crosses zero. Bottom
panel: the ranked t-statistic, which shows how much signal each rank
carries. The legend carries the set name, its genes present in the
ranked universe, its NES and its adjusted p. NES > 0 = enriched in
the contrast numerator (39 °C or WT). NES < 0 = enriched in the
denominator (37 °C or cGAS-KO). This set is 6th by |NES| and 4th by
adjusted p in its cell, so the general running-sum panels of this
GSEA sweep, which draw the top 5 per cell by |NES|, leave it out by
construction. SELECTION RULE: this set was chosen BY NAME, and its
adjusted p and |NES| played no part in that choice. The general
per-database and pooled panels of this GSEA sweep pin no set by name,
and no pin from this figure is carried back into them. This panel
draws one set out of the 50 that the Hallmark collection contributed
to this contrast, so read the two ranks quoted above as the whole of
what it says about where this set stands among them. A set's
enrichment is not evidence that the program its name invokes is
present; composition would have to be established separately. What
this curve locates is where the set's member genes sit in one
ranking. Claim tier: L3 (DE and enrichment statistics). PROVISIONAL
sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/28_hypoxia_focus_viz.R` | `gsea_running_sum_plot` | `thresholds.gsea_fdr=0.05; figures.running_sum_ylim=[-1.0,1.0]; figures.running_sum_heights; figures.running_sum_top=5; figures.top_pathways=20; colors.okabe_ito` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_Hallmark.rds}` |

## figures/_overview/hypoxia_routes_by_contrast.png

Running-enrichment curves for four hypoxia-named gene sets, one per
database, across the three focal contrasts. WT: heat (39 vs 37 °C):
Hallmark NES +1.91 (padj 4.16e-06), GO MF NES +1.75 (padj 0.0408), GO
BP NES +1.71 (padj 0.00113), Reactome NES -0.91 (padj 0.723).
cGAS-KO: heat (39 vs 37 °C): Hallmark NES +1.95 (padj 2.48e-06), GO
MF NES +1.77 (padj 0.0431), GO BP NES +1.54 (padj 0.0162), Reactome
NES -0.92 (padj 0.716). Heat × genotype interaction (cGAS-dependence
of the heat response): Hallmark NES -1.27 (padj 0.169), GO MF NES
+0.67 (padj 0.987), GO BP NES +1.34 (padj 0.31), Reactome NES -0.78
(padj 0.942). Three of the four sets carry a positive NES and reach
FDR < 0.05 in both per-genotype heat contrasts, and the fourth,
REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA, carries a negative NES and no
significance in any of the three; 6 of the 12 (set x contrast) cells
reach FDR < 0.05. Every Interaction NES is non-significant, which
reads as no detectable cGAS-dependence at n=5 and never as proven
cGAS-independence: the 1-df interaction term is the least-powered
contrast in this design. All four sets were chosen by name, so the
panel reports these four curves and ranks nothing.

**How to read:** One facet per contrast, and within each facet four overlaid
running-enrichment curves, one per gene set, keyed by colour in the
shared legend below the panels. The x axis is the rank of every gene
in THAT contrast's ranked list, most numerator-shifted on the left,
most denominator-shifted on the right, so the facets share an axis
length and each carries its own ordering. A curve steps up at each
gene belonging to its set and decays between them: an early high peak
means the members are packed at the numerator end, and a curve near
zero means they are spread through the list. The y range is pinned to
[-1, 1] so every curve stays comparable to every other running-sum
figure in this GSEA sweep. Inside each facet the NES and adjusted p
for that contrast are printed in each set's own colour, in the
legend's own top-to-bottom order, keyed by database name because each
of the four sets comes from a different database. The member ticks
and the ranked-metric panel are left off to keep four curves per
facet legible; the companion figure hypoxia_running_sum_wt_heat.png
draws all three panels for one set. NES > 0 = enriched in the
contrast numerator (39 °C or WT). NES < 0 = enriched in the
denominator (37 °C or cGAS-KO). The Interaction facet is the
cGAS-dependence read-out: a positive significant Interaction NES is
consistent with cGAS-dependent induction, and a non-significant one
means no detectable cGAS-dependence at n=5. SELECTION RULE: these 4
sets were chosen BY NAME, one per database, with no input from their
adjusted p or |NES|. The general per-database and pooled panels of
this GSEA sweep pin no set by name, and no pin from this figure is
carried back into them. So this is a named-set read-out and says
nothing about how these four rank among the thousands of sets they
came from. Read the fourth curve as carefully as the other three:
three of these hypoxia-named sets carry a positive NES in both heat
contrasts and the Reactome set carries a negative one, so shared
wording in two set names leaves shared behaviour an open question. A
set's enrichment is not evidence that the program its name invokes is
present; composition would have to be established separately. The
four names overlap in wording while their gene content differs, which
is why each is drawn on its own curve. Claim tier: L3 (DE and
enrichment statistics). PROVISIONAL sample labels.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/28_hypoxia_focus_viz.R` | `geom_line / facet_wrap` | `thresholds.gsea_fdr=0.05; figures.running_sum_ylim=[-1.0,1.0]; figures.running_sum_heights; figures.running_sum_top=5; figures.top_pathways=20; colors.okabe_ito` | `03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_{Hallmark,GO_MF,GO_BP,Reactome}.rds}` |

