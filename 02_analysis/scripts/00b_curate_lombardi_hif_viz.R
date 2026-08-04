#!/usr/bin/env Rscript
# =============================================================================
# 00b_curate_lombardi_hif_viz.R  --  VIZ: Lombardi-2022 HIF signature derivation
# =============================================================================
# Project: GSE329522 STING/cGAS Hyperthermia iTreg
# Role:    VISUALIZE half of the "normalize-then-visualize" split. Reads ONLY the
#          plot-ready tidy tables emitted by the curation COMPUTE script
#          (00b_curate_lombardi_hif.R) and renders ONE single-claim figure.
#          Performs NO statistics (no recurrence counting, no Poisson-binomial
#          null, no p.adjust); it only plots already-computed columns. Runs
#          STANDALONE after the compute script.
#
# Inputs (00_data/references/gene_sets/tables/):
#   lombardi_recurrence_data.csv  (per-gene: recurrence_count, null_padj, kept, ...)
#   lombardi_recurrence_null.csv  (per recurrence-level 0..32: observed count,
#                                  null_expected_at_count, clears_consensus_alpha)
#
# Output (00_data/references/gene_sets/figures/):
#   lombardi_recurrence.pdf
#
# CLAIM: documents exactly how the published Lombardi HIF signature was
#        reconstructed -- observed cross-TCGA recurrence vs a Poisson-binomial
#        null; the consensus genes are the conserved right tail (recurrence >= 6
#        of 32 cancers), sitting far above the null expectation.
# Dependencies: config.R; ggplot2, dplyr
# =============================================================================

source("02_analysis/config/config.R")
source("02_analysis/helpers/figure_style.R")       # contract: project_theme + save_overview; FIG_CFG
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

STAGE   <- "../00_data/references/gene_sets"
SCRIPT  <- "02_analysis/scripts/00b_curate_lombardi_hif_viz.R"

GENESET_DIR <- file.path(PROJECT_ROOT, "00_data/references/gene_sets")
TBL_DIR     <- file.path(GENESET_DIR, "tables")
FIG_DIR     <- file.path(GENESET_DIR, "figures")

DOWN <- DIVERGING_COLORS$negative   # blue
UP   <- DIVERGING_COLORS$positive   # orange
cap  <- sample_mapping_stamp()


rd <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)

null_df <- rd("lombardi_recurrence_null.csv")
gene_df <- rd("lombardi_recurrence_data.csv")

# --- derive PLOT LABELS only (no statistics) -------------------------------
null_df$clears_consensus_alpha <- as.logical(null_df$clears_consensus_alpha)
gene_df$kept <- as.logical(gene_df$kept)

# Chosen threshold = the lowest recurrence level that clears the consensus alpha
# (read straight from the precomputed clears_consensus_alpha flag).
threshold <- min(null_df$recurrence_level[null_df$clears_consensus_alpha])
n_datasets <- max(null_df$recurrence_level)
n_kept <- sum(gene_df$kept)

# Long form: observed vs null-expected counts per recurrence level (same y-axis).
obs <- data.frame(recurrence_level = null_df$recurrence_level,
                  count = null_df$n_genes_observed, series = "observed",
                  stringsAsFactors = FALSE)
exp <- data.frame(recurrence_level = null_df$recurrence_level,
                  count = null_df$null_expected_at_count, series = "null-expected",
                  stringsAsFactors = FALSE)
plot_df <- rbind(obs, exp)
plot_df$series <- factor(plot_df$series, levels = c("observed", "null-expected"))
plot_df$above  <- plot_df$recurrence_level >= threshold

series_cols <- c("observed" = UP, "null-expected" = DOWN)

# log1p y-axis so the conserved right tail (tiny counts) and the huge background
# (50k genes at recurrence 0) are both legible on the SAME axis.
fig <- ggplot(plot_df, aes(x = recurrence_level, y = count, fill = series)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_vline(xintercept = threshold - 0.5, color = "grey25", linetype = "22", linewidth = 0.7) +
  annotate("rect", xmin = threshold - 0.5, xmax = n_datasets + 0.5,
           ymin = 0, ymax = Inf, fill = UP, alpha = 0.06) +
  annotate("label", x = threshold - 0.5, y = max(plot_df$count), hjust = 0.0, vjust = 1,
           label = sprintf("consensus threshold: recurrence >= %d / %d\n%d genes retained (conserved right tail)",
                           threshold, n_datasets, n_kept),
           size = 3, fill = "grey96", color = "grey15") +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0, 1, 5, 10, 50, 100, 500, 5000, 50000)) +
  scale_fill_manual(values = series_cols, name = NULL) +
  scale_x_continuous(breaks = seq(0, n_datasets, by = 4)) +
  labs(
    title = "Lombardi-2022 HIF signature: observed cross-TCGA recurrence vs Poisson-binomial null",
    subtitle = sprintf("Genes ranked by how many of %d TCGA cancers call them HIF-correlated. Observed (orange) detaches from the\nnull expectation (blue) in the right tail; the consensus set is everything at recurrence >= %d / %d (shaded).",
                       n_datasets, threshold, n_datasets),
    x = sprintf("recurrence (number of %d TCGA cancers)", n_datasets),
    y = "genes at this recurrence level (log1p axis)",
    caption = cap
  ) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.minor = element_blank(),
        legend.position = "top")

# Remove stale flat-dir PDF if it exists
.stale_pdf <- file.path(FIG_DIR, "lombardi_recurrence.pdf")
if (file.exists(.stale_pdf)) file.remove(.stale_pdf)

# Clean up stale deck asset PNG if it exists in the wrong stage directory
.stale_deck_png <- file.path(DIR_RESULTS, "04_tf", "deck_assets", "lombardi_recurrence.png")
if (file.exists(.stale_deck_png)) file.remove(.stale_deck_png)

save_overview(
  fig, STAGE, "lombardi_recurrence",
  table = gene_df[, c("human_symbol", "recurrence_count", "null_padj", "kept")],
  finding = paste0(
    "Consensus reconstruction of the published Lombardi-2022 HIF signature: observed ",
    "cross-TCGA recurrence (orange) detaches from the Poisson-binomial null expectation ",
    "(blue) in the right tail; consensus genes are retained at recurrence >= 6 of 32 cancers. ",
    cap),
  script = SCRIPT, fn = "save_overview",
  config_kv = "figures.base_size=16; figures.base_size_column=9; consensus_threshold=6",
  input = "00_data/references/gene_sets/tables/lombardi_recurrence_{data,null}.csv",
  how_to_read = paste0(
    "Observed (orange) vs Poisson-binomial null-expected (blue) gene counts at each ",
    "recurrence level (log1p y-axis). Consensus threshold (recurrence >= 6) is shaded. ",
    "Claim tier: curation methodology documentation (illustrative null test)."),
  width = 9, height = 6, config = FIG_CFG)

cat("[DONE] 00b_curate_lombardi_hif_viz.R -- recurrence derivation figure rendered.\n")
