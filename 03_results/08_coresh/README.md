# 08_coresh — artifact captions

## figures/_overview/coresh_pctvar_overview.png

Top public mouse GEO datasets co-regulating the project's WT_heat_up
signature (Q_sig_WT_heat_up), ranked by CoReSh pctVar.

**How to read:** Each bar = one public GEO dataset. Bar length = pctVar (% of variance
in that dataset explained by the query), a PCA-inspired co-regulation
score (higher = stronger co-regulation); label = GSE accession,
mapped query size k, and rank. The query is the mouse-native
WT_heat_up signature exported from 17_signature_sets.rds. Claim tier:
L3-DE (data-driven, compendium coregulation score); sample labels
PROVISIONAL. Sign convention: pctVar >= 0 (variance fraction, not
signed).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_pctvar_overview` | `coresh.top_n_hits=5; figures.top_n=20; thresholds.gsea_fdr=0.05; colors.diverging; colors.okabe_ito` | `03_results/objects/coresh_ranked.rds` |

## figures/_overview/coresh_nes_dotplot.png

Cross-contrast enrichment of CoReSh-derived co-regulation sets (NES
fill, -log10(padj) size, FDR<0.05 outline); left strip = seeding
query origin. All expected sets are derived from the WT_heat_up
signature query.

**How to read:** Each circle = one CoReSh-derived gene set (row) in one contrast
(column). Fill color = NES (Normalized Enrichment Score): orange =
positive enrichment (up in numerator), blue = negative enrichment
(down); scale clamped at ±3.5 for visual comparability. Circle size =
-log10(padj); larger = more significant. Black outline = FDR < 0.05
(significant); no outline = not significant. Left color strip =
biological origin of the query that seeded the derived set
(WT_heat_up signature for the current CoReSh arm). Sets ordered by
median NES across contrasts (highest at top). Claim tier: L3-DE
(fgsea, BH-FDR). Sample labels PROVISIONAL (inferred from heat-shock
thermometer + Cgas expression).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_coresh_nes_dotplot` | `figures.top_n=20; figures.nes_cap=3.5; thresholds.gsea_fdr=0.05; colors.diverging; colors.okabe_ito` | `03_results/objects/gsea_coresh_<contrast>.rds (fallback: 03_results/master/master_gsea_table.csv rows database=CoReSh_derived)` |

## figures/by_contrast/WT_heat/coresh_gsea_lollipop.png

WT_heat: top CoReSh-derived WT_heat_up co-regulation sets by |NES|;
positive NES (orange ↑ up) = set enriched among up-regulated genes
(numerator-high); negative NES (blue ↓ down) = enriched among
down-regulated. Filled dots = FDR < 0.05; open = n.s.

**How to read:** Each row = one CoReSh-derived gene set (named CORESH_<query>_<GSE>).
Horizontal bar length = NES magnitude; direction = sign: rightward
(orange) = up-enriched, leftward (blue) = down-enriched in this
contrast. x-axis clamped at ±3.5. Filled dot = FDR < 0.05
(significant). Set name format: 'Signature: WT heat up <GSE>' — GSE
is the public dataset whose co-regulation pattern was used to derive
the gene set. Claim tier: L3-DE (fgsea multilevel, BH-FDR). Sign
convention: NES > 0 = ↑ up; NES < 0 = ↓ down. Sample labels
PROVISIONAL (inferred thermometer + Cgas). No detectable
cGAS-dependence at n=5 (not cGAS-independent).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_coresh_lollipop` | `figures.top_n=20; figures.nes_cap=3.5; thresholds.gsea_fdr=0.05; colors.diverging` | `03_results/objects/gsea_coresh_WT_heat.rds` |

## figures/by_contrast/KO_heat/coresh_gsea_lollipop.png

KO_heat: top CoReSh-derived WT_heat_up co-regulation sets by |NES|;
positive NES (orange ↑ up) = set enriched among up-regulated genes
(numerator-high); negative NES (blue ↓ down) = enriched among
down-regulated. Filled dots = FDR < 0.05; open = n.s.

**How to read:** Each row = one CoReSh-derived gene set (named CORESH_<query>_<GSE>).
Horizontal bar length = NES magnitude; direction = sign: rightward
(orange) = up-enriched, leftward (blue) = down-enriched in this
contrast. x-axis clamped at ±3.5. Filled dot = FDR < 0.05
(significant). Set name format: 'Signature: WT heat up <GSE>' — GSE
is the public dataset whose co-regulation pattern was used to derive
the gene set. Claim tier: L3-DE (fgsea multilevel, BH-FDR). Sign
convention: NES > 0 = ↑ up; NES < 0 = ↓ down. Sample labels
PROVISIONAL (inferred thermometer + Cgas). No detectable
cGAS-dependence at n=5 (not cGAS-independent).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_coresh_lollipop` | `figures.top_n=20; figures.nes_cap=3.5; thresholds.gsea_fdr=0.05; colors.diverging` | `03_results/objects/gsea_coresh_KO_heat.rds` |

## figures/by_contrast/Interaction/coresh_gsea_lollipop.png

Interaction: top CoReSh-derived WT_heat_up co-regulation sets by
|NES|; positive NES (orange ↑ up) = set enriched among up-regulated
genes (numerator-high); negative NES (blue ↓ down) = enriched among
down-regulated. Filled dots = FDR < 0.05; open = n.s.

**How to read:** Each row = one CoReSh-derived gene set (named CORESH_<query>_<GSE>).
Horizontal bar length = NES magnitude; direction = sign: rightward
(orange) = up-enriched, leftward (blue) = down-enriched in this
contrast. x-axis clamped at ±3.5. Filled dot = FDR < 0.05
(significant). Set name format: 'Signature: WT heat up <GSE>' — GSE
is the public dataset whose co-regulation pattern was used to derive
the gene set. Claim tier: L3-DE (fgsea multilevel, BH-FDR). Sign
convention: NES > 0 = ↑ up; NES < 0 = ↓ down. Sample labels
PROVISIONAL (inferred thermometer + Cgas). No detectable
cGAS-dependence at n=5 (not cGAS-independent).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_coresh_lollipop` | `figures.top_n=20; figures.nes_cap=3.5; thresholds.gsea_fdr=0.05; colors.diverging` | `03_results/objects/gsea_coresh_Interaction.rds` |

## figures/by_contrast/Geno_at_39/coresh_gsea_lollipop.png

Geno_at_39: top CoReSh-derived WT_heat_up co-regulation sets by
|NES|; positive NES (orange ↑ up) = set enriched among up-regulated
genes (numerator-high); negative NES (blue ↓ down) = enriched among
down-regulated. Filled dots = FDR < 0.05; open = n.s.

**How to read:** Each row = one CoReSh-derived gene set (named CORESH_<query>_<GSE>).
Horizontal bar length = NES magnitude; direction = sign: rightward
(orange) = up-enriched, leftward (blue) = down-enriched in this
contrast. x-axis clamped at ±3.5. Filled dot = FDR < 0.05
(significant). Set name format: 'Signature: WT heat up <GSE>' — GSE
is the public dataset whose co-regulation pattern was used to derive
the gene set. Claim tier: L3-DE (fgsea multilevel, BH-FDR). Sign
convention: NES > 0 = ↑ up; NES < 0 = ↓ down. Sample labels
PROVISIONAL (inferred thermometer + Cgas). No detectable
cGAS-dependence at n=5 (not cGAS-independent).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_coresh_lollipop` | `figures.top_n=20; figures.nes_cap=3.5; thresholds.gsea_fdr=0.05; colors.diverging` | `03_results/objects/gsea_coresh_Geno_at_39.rds` |

## figures/by_contrast/Geno_at_37/coresh_gsea_lollipop.png

Geno_at_37: top CoReSh-derived WT_heat_up co-regulation sets by
|NES|; positive NES (orange ↑ up) = set enriched among up-regulated
genes (numerator-high); negative NES (blue ↓ down) = enriched among
down-regulated. Filled dots = FDR < 0.05; open = n.s.

**How to read:** Each row = one CoReSh-derived gene set (named CORESH_<query>_<GSE>).
Horizontal bar length = NES magnitude; direction = sign: rightward
(orange) = up-enriched, leftward (blue) = down-enriched in this
contrast. x-axis clamped at ±3.5. Filled dot = FDR < 0.05
(significant). Set name format: 'Signature: WT heat up <GSE>' — GSE
is the public dataset whose co-regulation pattern was used to derive
the gene set. Claim tier: L3-DE (fgsea multilevel, BH-FDR). Sign
convention: NES > 0 = ↑ up; NES < 0 = ↓ down. Sample labels
PROVISIONAL (inferred thermometer + Cgas). No detectable
cGAS-dependence at n=5 (not cGAS-independent).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_coresh_lollipop` | `figures.top_n=20; figures.nes_cap=3.5; thresholds.gsea_fdr=0.05; colors.diverging` | `03_results/objects/gsea_coresh_Geno_at_37.rds` |

## figures/by_contrast/Temp_main/coresh_gsea_lollipop.png

Temp_main: top CoReSh-derived WT_heat_up co-regulation sets by |NES|;
positive NES (orange ↑ up) = set enriched among up-regulated genes
(numerator-high); negative NES (blue ↓ down) = enriched among
down-regulated. Filled dots = FDR < 0.05; open = n.s.

**How to read:** Each row = one CoReSh-derived gene set (named CORESH_<query>_<GSE>).
Horizontal bar length = NES magnitude; direction = sign: rightward
(orange) = up-enriched, leftward (blue) = down-enriched in this
contrast. x-axis clamped at ±3.5. Filled dot = FDR < 0.05
(significant). Set name format: 'Signature: WT heat up <GSE>' — GSE
is the public dataset whose co-regulation pattern was used to derive
the gene set. Claim tier: L3-DE (fgsea multilevel, BH-FDR). Sign
convention: NES > 0 = ↑ up; NES < 0 = ↓ down. Sample labels
PROVISIONAL (inferred thermometer + Cgas). No detectable
cGAS-dependence at n=5 (not cGAS-independent).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_coresh_lollipop` | `figures.top_n=20; figures.nes_cap=3.5; thresholds.gsea_fdr=0.05; colors.diverging` | `03_results/objects/gsea_coresh_Temp_main.rds` |

## figures/by_contrast/Geno_main/coresh_gsea_lollipop.png

Geno_main: top CoReSh-derived WT_heat_up co-regulation sets by |NES|;
positive NES (orange ↑ up) = set enriched among up-regulated genes
(numerator-high); negative NES (blue ↓ down) = enriched among
down-regulated. Filled dots = FDR < 0.05; open = n.s.

**How to read:** Each row = one CoReSh-derived gene set (named CORESH_<query>_<GSE>).
Horizontal bar length = NES magnitude; direction = sign: rightward
(orange) = up-enriched, leftward (blue) = down-enriched in this
contrast. x-axis clamped at ±3.5. Filled dot = FDR < 0.05
(significant). Set name format: 'Signature: WT heat up <GSE>' — GSE
is the public dataset whose co-regulation pattern was used to derive
the gene set. Claim tier: L3-DE (fgsea multilevel, BH-FDR). Sign
convention: NES > 0 = ↑ up; NES < 0 = ↓ down. Sample labels
PROVISIONAL (inferred thermometer + Cgas). No detectable
cGAS-dependence at n=5 (not cGAS-independent).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/15_coresh_viz.R` | `create_coresh_lollipop` | `figures.top_n=20; figures.nes_cap=3.5; thresholds.gsea_fdr=0.05; colors.diverging` | `03_results/objects/gsea_coresh_Geno_main.rds` |

