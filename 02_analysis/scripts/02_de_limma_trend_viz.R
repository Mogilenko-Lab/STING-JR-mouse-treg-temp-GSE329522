#!/usr/bin/env Rscript
# =============================================================================
# 02_de_limma_trend_viz.R  --  Phase 2 VIZ: cGAS-dependence figure + volcanoes
# =============================================================================
# Phase:    2 (stage 03_de)
# Role:     VISUALIZE half of the "normalize-then-visualize" split. Reads the
#           plot-ready tidy table + the DE checkpoint emitted by
#           02_de_limma_trend.R and renders figures. Performs NO statistics
#           (no lmFit/eBayes/topTable/makeContrasts/p.adjust); it only plots
#           already-computed columns. Runs STANDALONE after the compute script.
#
# Inputs:   03_results/03_de/tables/fig2_marker_means.csv   (per-group means +
#                                           Interaction logFC/adjP, tidy)
#           03_results/objects/02_de_results.rds             (7 topTables, for
#                                           volcanoes -- existing columns only)
# Outputs:  03_results/03_de/figures/fig2_cgas_dependence_markers.pdf
#           03_results/03_de/figures/volcano_{WT_heat,KO_heat,Interaction}.pdf
#
# Dependencies: ggplot2, dplyr
# =============================================================================

source("02_analysis/config/config.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# -----------------------------------------------------------------------------
# 1. Read plot-ready inputs (NO recomputation)
# -----------------------------------------------------------------------------
fig2_means_path <- file.path(stage_dir("03_de", "tables"), "fig2_marker_means.csv")
if (!file.exists(fig2_means_path)) {
  stop(sprintf("Missing %s -- run 02_de_limma_trend.R (compute) first.", fig2_means_path))
}
gm_long <- read.csv(fig2_means_path, stringsAsFactors = FALSE, check.names = FALSE)

de_path <- file.path(DIR_OBJECTS, "02_de_results.rds")
if (!file.exists(de_path)) {
  stop(sprintf("Missing %s -- run 02_de_limma_trend.R (compute) first.", de_path))
}
de_results <- readRDS(de_path)

# -----------------------------------------------------------------------------
# 2. FIG 2 -- cGAS-dependence marker readout (the money panel)
# -----------------------------------------------------------------------------
# Per-group mean log2CPM, points + lines connecting 37->39 within genotype,
# faceted into IFN/ISG arm vs HIF/glycolysis arm; HIF-specific vs shared-glyco
# distinguished by linetype. Each gene annotated with Interaction adj.P.
# All numeric inputs come straight from fig2_marker_means.csv -- no statistics.

# Reconstruct the plotting subclass (linetype key) from hif_class: ISG-arm genes
# carry hif_class = NA in the tidy table and map to the "ISG" linetype class.
gm_long <- gm_long %>%
  mutate(
    subclass = ifelse(is.na(hif_class) | hif_class == "", "ISG", hif_class),
    temp     = factor(temp, levels = c("37C", "39C")),
    genotype = factor(genotype, levels = c("WT", "cGASKO")),
    arm      = factor(arm, levels = c("IFN / ISG arm", "HIF / glycolysis arm"))
  )

# gene ordering preserved from the table (ISG arm first, then HIF/glyco arm)
fig_genes <- unique(gm_long$gene)

# per-gene Interaction adj.P annotation (placed at top of each gene's panel).
# adjP is read from the tidy table (constant per gene) -- not recomputed.
ann <- gm_long %>%
  group_by(gene, arm) %>%
  summarise(adjP = inter_adjP[1],
            subclass = subclass[1],
            ymax = max(mean_log2cpm),
            ymin = min(mean_log2cpm),
            .groups = "drop") %>%
  mutate(lab = sprintf("int.adjP=%.2g%s", adjP, ifelse(adjP < DE_FDR, " *", " ns")))
ann$arm <- factor(ann$arm, levels = levels(gm_long$arm))
# y position for label = top of each gene's range
ann$ylab <- ann$ymax + 0.12 * (max(gm_long$mean_log2cpm) - min(gm_long$mean_log2cpm)) / 6 + 0.3

gcols <- c(WT = "#B35806", cGASKO = "#2166AC")  # diverging-ish: WT warm, KO cool

p <- ggplot(gm_long, aes(x = temp, y = mean_log2cpm,
                         color = genotype, group = genotype, linetype = subclass)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.9) +
  facet_wrap(~ arm + gene, scales = "free_y", ncol = 4,
             labeller = labeller(.multi_line = FALSE)) +
  geom_text(data = ann, aes(x = 1.5, y = ylab, label = lab),
            inherit.aes = FALSE, size = 2.5, vjust = 0, color = "grey20") +
  scale_color_manual(values = gcols, name = "genotype") +
  scale_linetype_manual(values = c("ISG" = "solid",
                                   "HIF-specific" = "solid",
                                   "shared-glycolytic" = "22"),
                        name = "gene class") +
  labs(
    title = "cGAS-dependence of the heat response: ISG arm vs HIF/glycolysis arm",
    subtitle = "Per-group mean log2(CPM+0.5); lines connect 37C->39C within genotype. Interaction adj.P annotated per gene (* = adj.P < 0.05).",
    x = "temperature", y = "mean log2(CPM+0.5)",
    caption = paste0(provisional_caption(),
                     "\nDirection labels rest on the marker prior (KO lowers Cgas; heat raises HSPs), pending sample sheet.",
                     "\nn=5/group; the Interaction is the lowest-powered term -- a non-significant HIF interaction means NO DETECTABLE cGAS-dependence at n=5, NOT proof of independence.")
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        plot.caption = element_text(size = 6.5, hjust = 0),
        strip.text = element_text(size = 7))

fig2_path <- file.path(stage_dir("03_de", "figures"), "fig2_cgas_dependence_markers.pdf")
ggsave(fig2_path, p, width = 12, height = 7)
cat("[SAVE]", fig2_path, "\n")

# -----------------------------------------------------------------------------
# 3. Per-contrast volcanoes (toolkit create_standard_volcano)
# -----------------------------------------------------------------------------
# Plots existing logFC / P.Value / adj.P.Val columns from 02_de_results.rds.
# The toolkit helper does its own thresholding for COLORING only -- no DE stats
# are (re)computed here.
volc_ok <- tryCatch({
  source("01_modules/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R")
  for (cn in c("WT_heat", "KO_heat", "Interaction")) {
    df <- de_results[[cn]]
    rownames(df) <- make.unique(df$gene_symbol)
    g <- create_standard_volcano(
      de_results = df[, c("logFC", "P.Value", "adj.P.Val")],
      decision_by = "fdr", p_cutoff = DE_FDR, fc_cutoff = 1, top_n = 8,
      title = sprintf("Volcano: %s (limma-trend)", cn),
      subtitle = provisional_caption()
    )
    vp <- file.path(stage_dir("03_de", "figures"), sprintf("volcano_%s.pdf", cn))
    ggsave(vp, g, width = 8, height = 6)
    cat("[SAVE]", vp, "\n")
  }
  TRUE
}, error = function(e) { cat("[WARN] volcano step skipped:", conditionMessage(e), "\n"); FALSE })

cat("\n=== Phase 2 VIZ complete | volcano_ok =", volc_ok, "===\n")
