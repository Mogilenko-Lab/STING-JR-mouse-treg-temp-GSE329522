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
load_packages()

GENESET_DIR <- file.path(PROJECT_ROOT, "00_data/references/gene_sets")
TBL_DIR     <- file.path(GENESET_DIR, "tables")
FIG_DIR     <- file.path(GENESET_DIR, "figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

DOWN <- DIVERGING_COLORS$negative   # blue
UP   <- DIVERGING_COLORS$positive   # orange
cap  <- provisional_caption()

# DECK EXPORT (TASK 4): canonical stamped PDF always; clean 300-dpi PNG (no
# embargo caption) into 03_results/04_tf/deck_assets/ when DECK_EXPORT is set.
DECK_EXPORT <- nzchar(Sys.getenv("DECK_EXPORT"))
DECK_DIR    <- file.path(DIR_RESULTS, "04_tf", "deck_assets")
if (DECK_EXPORT) dir.create(DECK_DIR, recursive = TRUE, showWarnings = FALSE)

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
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey25"),
        plot.caption  = element_text(size = 7, color = "grey45", hjust = 0, face = "italic"),
        legend.position = "top")

ggsave(file.path(FIG_DIR, "lombardi_recurrence.pdf"), fig, width = 9, height = 6)
cat("  wrote", file.path(FIG_DIR, "lombardi_recurrence.pdf"), "\n")
if (DECK_EXPORT) {
  ggsave(file.path(DECK_DIR, "lombardi_recurrence.png"), fig + labs(caption = NULL),
         width = 10, height = 5.625, dpi = 300)
  cat("  wrote deck_assets/lombardi_recurrence.png (clean, no embargo stamp)\n")
}
cat("[DONE] 00b_curate_lombardi_hif_viz.R -- recurrence derivation figure rendered.\n")
