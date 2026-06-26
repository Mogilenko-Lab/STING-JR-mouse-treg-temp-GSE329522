# 09_gatom — artifact captions

## figures/by_contrast/WT_heat/module_graph_kegg.png

WT: heat (39 vs 37 °C) (kegg): GATOM SGMWCS metabolic module — nodes
are metabolites/atoms; edges are enzymatic reactions colored by gene
log2FC (orange = up-regulated at 39°C, blue = down-regulated; no
metabolomics data supplied so node size is uniform). Module V=7 E=6
w=47.01. This is an expression-only GATOM run (k.met=NULL). Claim
tier L3 (enrichment statistics; mechanism is interpretive text only).
PROVISIONAL - inferred sample mapping (Hspa1b/Hsph1 thermometer +
Cgas), pending collaborator sample sheet

**How to read:** Nodes = metabolites/atoms (uniform size — no metabolomics). Edges =
enzymatic reactions; gene symbol labels on edges. Edge color: orange
= up-regulated enzyme (log2FC > 0), blue = down-regulated. Edge width
(when available): -log10(raw p-value) from GATOM BUM scoring. Sign
convention: positive log2FC = higher in numerator of the contrast.
The direction cue (subtitle) reflects the mean signed log2FC of
module edges. Top-20 most connected nodes are labeled to reduce
overplotting. Claim tier L3: the module is a statistically optimal
connected subgraph, NOT a direct measurement of metabolic flux.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_gatom_viz.R` | `build_module_ggraph` | `colors.diverging; figures.top_n=20; thresholds.gatom_k_gene=50` | `03_results/objects/10_gatom_WT_heat.rds` |

## figures/by_contrast/KO_heat/module_graph_kegg.png

cGAS-KO: heat (39 vs 37 °C) (kegg): GATOM SGMWCS metabolic module —
nodes are metabolites/atoms; edges are enzymatic reactions colored by
gene log2FC (orange = up-regulated at 39°C, blue = down-regulated; no
metabolomics data supplied so node size is uniform). Module V=18 E=17
w=44.83. This is an expression-only GATOM run (k.met=NULL). Claim
tier L3 (enrichment statistics; mechanism is interpretive text only).
PROVISIONAL - inferred sample mapping (Hspa1b/Hsph1 thermometer +
Cgas), pending collaborator sample sheet

**How to read:** Nodes = metabolites/atoms (uniform size — no metabolomics). Edges =
enzymatic reactions; gene symbol labels on edges. Edge color: orange
= up-regulated enzyme (log2FC > 0), blue = down-regulated. Edge width
(when available): -log10(raw p-value) from GATOM BUM scoring. Sign
convention: positive log2FC = higher in numerator of the contrast.
The direction cue (subtitle) reflects the mean signed log2FC of
module edges. Top-20 most connected nodes are labeled to reduce
overplotting. Claim tier L3: the module is a statistically optimal
connected subgraph, NOT a direct measurement of metabolic flux.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_gatom_viz.R` | `build_module_ggraph` | `colors.diverging; figures.top_n=20; thresholds.gatom_k_gene=50` | `03_results/objects/10_gatom_KO_heat.rds` |

## figures/by_contrast/Interaction/module_graph_kegg.png

Heat × genotype interaction (cGAS-dependence of the heat response)
(kegg): GATOM SGMWCS metabolic module — nodes are metabolites/atoms;
edges are enzymatic reactions colored by gene log2FC (orange =
up-regulated at 39°C, blue = down-regulated; no metabolomics data
supplied so node size is uniform). Module V=n/a E=n/a w=n/a. This is
an expression-only GATOM run (k.met=NULL). Claim tier L3 (enrichment
statistics; mechanism is interpretive text only). PROVISIONAL -
inferred sample mapping (Hspa1b/Hsph1 thermometer + Cgas), pending
collaborator sample sheet

**How to read:** Nodes = metabolites/atoms (uniform size — no metabolomics). Edges =
enzymatic reactions; gene symbol labels on edges. Edge color: orange
= up-regulated enzyme (log2FC > 0), blue = down-regulated. Edge width
(when available): -log10(raw p-value) from GATOM BUM scoring. Sign
convention: positive log2FC = higher in numerator of the contrast.
The direction cue (subtitle) reflects the mean signed log2FC of
module edges. Top-20 most connected nodes are labeled to reduce
overplotting. Claim tier L3: the module is a statistically optimal
connected subgraph, NOT a direct measurement of metabolic flux.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_gatom_viz.R` | `build_module_ggraph` | `colors.diverging; figures.top_n=20; thresholds.gatom_k_gene=50` | `03_results/objects/10_gatom_Interaction.rds` |

## figures/by_contrast/Temp_main/module_graph_kegg.png

Heat main effect (pooled genotypes) (kegg): GATOM SGMWCS metabolic
module — nodes are metabolites/atoms; edges are enzymatic reactions
colored by gene log2FC (orange = up-regulated at 39°C, blue =
down-regulated; no metabolomics data supplied so node size is
uniform). Module V=6 E=5 w=48.32. This is an expression-only GATOM
run (k.met=NULL). Claim tier L3 (enrichment statistics; mechanism is
interpretive text only). PROVISIONAL - inferred sample mapping
(Hspa1b/Hsph1 thermometer + Cgas), pending collaborator sample sheet

**How to read:** Nodes = metabolites/atoms (uniform size — no metabolomics). Edges =
enzymatic reactions; gene symbol labels on edges. Edge color: orange
= up-regulated enzyme (log2FC > 0), blue = down-regulated. Edge width
(when available): -log10(raw p-value) from GATOM BUM scoring. Sign
convention: positive log2FC = higher in numerator of the contrast.
The direction cue (subtitle) reflects the mean signed log2FC of
module edges. Top-20 most connected nodes are labeled to reduce
overplotting. Claim tier L3: the module is a statistically optimal
connected subgraph, NOT a direct measurement of metabolic flux.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_gatom_viz.R` | `build_module_ggraph` | `colors.diverging; figures.top_n=20; thresholds.gatom_k_gene=50` | `03_results/objects/10_gatom_Temp_main.rds` |

## figures/_overview/module_sizes.png

GATOM metabolic module size (reaction edges) per contrast × network:
contrasts with the most connected metabolic sub-networks have the
highest edge count; WT_heat and Temp_main are expected to show the
largest modules (heat program recruits glycolytic + mitochondrial
enzymes). Claim tier L3.

**How to read:** Horizontal bar = one (contrast × network) combination; longer bar =
bigger module. Length = number of reaction edges (enzyme-encoded
enzymatic steps) in the SGMWCS module. Blue bars = KEGG reaction
network (the only network in this run). The Combined KEGG+Rhea
network is unavailable (met.combined.db.rds absent in
00_data/references/gatom/), so 10_gatom_modules.R degraded to
KEGG-only — there is no orange (combined) bar by design. An absent
contrast means its module was empty (GATOM returned 0 edges). Claim
tier L3: module size reflects DE signal density in the atom-graph,
not metabolic flux.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_gatom_viz.R` | `module_sizes_bar` | `colors.diverging; thresholds.gatom_k_gene=50` | `03_results/objects/10_gatom_WT_heat.rds; 03_results/objects/10_gatom_KO_heat.rds; 03_results/objects/10_gatom_Interaction.rds; 03_results/objects/10_gatom_Temp_main.rds` |

## figures/_overview/module_weights.png

GATOM MWCS solution weight per contrast × network: the objective
score of the SGMWCS problem integrates both module size and edge/node
score magnitudes; contrasts with strong uniform DE in metabolic
enzymes produce the highest weights. Claim tier L3.

**How to read:** Each bar = one (contrast × network) combination. Height = SGMWCS
objective value (dimensionless; from mwcsr::solve_mwcsp). A higher
weight indicates more/larger-magnitude enzyme DE signals concentrated
in a connected metabolic sub-network. Claim tier L3: does not imply
metabolic flux or enzyme activity directly.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_gatom_viz.R` | `module_weights_bar` | `colors.diverging; thresholds.gatom_k_gene=50` | `03_results/objects/10_gatom_WT_heat.rds; 03_results/objects/10_gatom_KO_heat.rds; 03_results/objects/10_gatom_Interaction.rds; 03_results/objects/10_gatom_Temp_main.rds` |

## figures/_overview/module_summary.png

Cross-contrast GATOM module summary: top panel shows module size
(reaction edges) per contrast × network; bottom panel shows mean
module enzyme log2FC (pseudo-NES), indicating whether the recruited
metabolic sub-network is net up- or down-regulated. Corroborates the
MitoCarta-anchored Complex-I / metabolic-pseudohypoxia mechanism at
the enrichment-statistics tier (L3). PROVISIONAL - inferred sample
mapping (Hspa1b/Hsph1 thermometer + Cgas), pending collaborator
sample sheet

**How to read:** Top panel: bar height = number of enzymatic reaction edges in the
SGMWCS module. Bottom panel: dot position = mean log2FC across all
enzyme-edges in the module. Blue = KEGG-only network (the only
network in this run; Combined KEGG+Rhea was unavailable —
met.combined.db.rds absent). Positive pseudo-NES (above dashed line)
= net up-regulation of module enzymes in the contrast numerator. Sign
convention: positive log2FC = higher in numerator (e.g. 39°C arm for
heat contrasts). Claim tier L3: module statistics, not metabolic flux
measurement. PROVISIONAL: sample-group labels are marker-inferred,
pending collaborator sample sheet.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/14_gatom_viz.R` | `module_summary_panel` | `colors.diverging; figures.top_n=20; thresholds.gatom_k_gene=50` | `03_results/objects/10_gatom_WT_heat.rds; 03_results/objects/10_gatom_KO_heat.rds; 03_results/objects/10_gatom_Interaction.rds; 03_results/objects/10_gatom_Temp_main.rds` |

