#!/usr/bin/env Rscript
# =============================================================================
# 03b_decoupler_method_comparison_viz.R  --  VIZ: decoupleR method comparison
# =============================================================================
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2: genotype x temperature)
# Phase:   3b (stage 04_tf)
#
# Role:    VISUALIZE half of the "normalize-then-visualize" split.  Reads ONLY
#          the plot-ready tidy tables emitted by 03b_decoupler_method_comparison.R
#          and renders ONE single-claim figure per PDF.  Performs NO statistics
#          (no run_ulm / run_mlm / run_consensus / p.adjust / cor / prcomp /
#          decouple).  Runs STANDALONE after the compute script.
#
# Inputs  (03_results/04_tf/tables/):
#   fig3j_allmethods_topTF_data.csv
#   fig3k_method_rank_divergence_data.csv
#   fig3k_method_rank_spearman.csv
#
# Outputs (03_results/04_tf/figures/_overview/), dual .print.pdf + .screen.png
# per stem via the figure-style contract (save_overview):
#   fig3j_topTF_allmethods_WT_heat, fig3j_topTF_allmethods_Interaction,
#   fig3k_method_rank_divergence  (+ same-stem source tables under tables/_overview/).
#
# Figure discipline (AGENTS.md): ONE claim = ONE dedicated figure.
# Dependencies: figure_style.R (contract); config.R (palette tokens); ggplot2, dplyr
# =============================================================================

source("02_analysis/config/config.R")          # palette tokens: AXIS_COLORS/AXIS_LABELS/DIVERGING_COLORS
source("02_analysis/helpers/figure_style.R")  # contract: FIG_CFG, project_theme, save_overview, ...
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
cap <- sample_mapping_stamp()

STAGE   <- "04_tf"
SCRIPT  <- "02_analysis/scripts/03b_decoupler_method_comparison_viz.R"
TBL_DIR <- stage_dir("04_tf", "tables")
TOP_N   <- as.integer(FIG_CFG$figures$top_n %||% 20L)

DOWN <- DIVERGING_COLORS$negative   # "#2166AC"  blue  (down)
NEUT <- DIVERGING_COLORS$neutral    # "#F7F7F7"  white
UP   <- DIVERGING_COLORS$positive   # "#B35806"  orange (up)

# Small read helper -- read pre-computed tidy tables only; zero statistics here.
rd <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE,
                            stringsAsFactors = FALSE)

# One-time relocation cleanup: this script previously wrote flat-dir PDFs into
# 04_tf/figures/. Now every owned stem lands in figures/_overview/. Remove the
# stale flat-dir artifacts for the owned stems (save_overview purges _overview/).
.FLAT_FIG_DIR <- stage_dir("04_tf", "figures")
.OWNED_STEMS  <- c("fig3j_topTF_allmethods_WT_heat",
                   "fig3j_topTF_allmethods_Interaction",
                   "fig3k_method_rank_divergence")
for (.stale in list.files(.FLAT_FIG_DIR, full.names = TRUE)) {
  if (any(startsWith(basename(.stale), .OWNED_STEMS)) &&
      grepl("\\.(pdf|png)$", .stale)) {
    file.remove(.stale)
    cat(sprintf("  removed stale flat-dir artifact %s\n", basename(.stale)))
  }
}

# Axis color map -- sourced from config.R (AXIS_COLORS) as the single source of truth.
# HIF = orange, IFN/NFkB = blue, other = grey55.
axis_cols <- AXIS_COLORS

# Canonical statistic order for facet display
STAT_ORDER <- c("ulm", "mlm", "wsum", "norm_wsum", "corr_wsum", "consensus")

# Human-readable labels for statistic facets
STAT_LABELS <- c(
  ulm       = "ULM",
  mlm       = "MLM",
  wsum      = "wsum",
  norm_wsum = "norm_wsum",
  corr_wsum = "corr_wsum",
  consensus = "consensus"
)

# =============================================================================
# HELPER: render a faceted top-TF lollipop figure for ONE contrast
# =============================================================================

make_fig3j <- function(fig3j_df, chosen_contrast, fig_label, out_fname) {
  d <- fig3j_df[fig3j_df$contrast == chosen_contrast, ]

  d$statistic <- factor(d$statistic, levels = STAT_ORDER, labels = STAT_LABELS)
  d$axis      <- factor(d$axis, levels = c("HIF", "IFN", "other"))
  d$key_tf    <- as.logical(d$key_tf)

  # Within each facet, order TFs by score so lollipops read cleanly.
  # We use a per-facet ordering key: score within (statistic, direction).
  # facet_wrap with scales="free_y" requires row_key as factor per subset.
  # We build a global row_key by concatenating statistic + sorted order.

  d <- d[order(d$statistic, d$score), ]
  d$row_key <- paste0(d$statistic, "_",
                      formatC(seq_len(nrow(d)), width = 4, flag = "0"))
  d$row_key <- factor(d$row_key, levels = d$row_key, labels = d$source)

  # Label only key TFs
  d_key <- d[d$key_tf, ]

  fig <- ggplot(d, aes(x = score, y = row_key)) +
    geom_vline(xintercept = 0, color = "grey65", linewidth = 0.3) +
    geom_segment(aes(xend = 0, yend = row_key, color = axis), linewidth = 0.55) +
    geom_point(aes(color = axis, size = key_tf)) +
    geom_text(data = d_key,
              aes(label = source, x = score),
              nudge_x = ifelse(d_key$score > 0, 0.12, -0.12),
              hjust   = ifelse(d_key$score > 0, 0,     1),
              size    = 2.2, color = "black") +
    facet_wrap(~ statistic, scales = "free", ncol = 3) +
    scale_color_manual(values = axis_cols, name = "TF axis",
                       labels = AXIS_LABELS) +
    scale_size_manual(values = c("TRUE" = 2.6, "FALSE" = 1.3), guide = "none") +
    labs(
      title    = sprintf("%s. Top activated / suppressed TFs under all 6 decoupleR statistics (%s)",
                         fig_label, chosen_contrast),
      subtitle = paste0(
        "Each facet shows the top-", 12, " up + top-", 12, " down TFs under one statistic.  ",
        "Axes differ across statistics (scales='free') -- the point is RANK/IDENTITY, not absolute score.\n",
        "HIF axis = orange, IFN/NFkB axis = blue.  Key TFs labelled.  ",
        "Shared network: CollecTRI (CACHED; get_collectri() disabled)."
      ),
      x = "TF activity score (scale not comparable across facets)",
      y = NULL
    ) +
    project_theme(config = FIG_CFG) +
    theme(axis.text.y = element_text(size = 5.5),
          panel.spacing = unit(0.9, "lines"))

  save_overview(
    fig, STAGE, sub("\\.pdf$", "", out_fname),
    table = d[, c("source", "statistic", "score", "axis", "key_tf")],
    finding = paste0(
      "Method comparison for contrast ", chosen_contrast, ": lollipops show the top ",
      "activated / suppressed TFs under all 6 decoupleR statistics. ", cap),
    script = SCRIPT, fn = "make_fig3j",
    config_kv = "figures.base_size=16; figures.base_size_column=9; statistics",
    input = "03_results/04_tf/tables/fig3j_allmethods_topTF_data.csv",
    how_to_read = paste0(
      "Lollipops = per-TF activity, faceted by statistic (ULM, MLM, wsum, norm_wsum, ",
      "corr_wsum, consensus), colored by axis (HIF orange; IFN/NFkB blue; other grey). ",
      "Key TFs are labelled. Note that y-axes differ (scales='free') to compare rank. ",
      "Claim tier: methodological comparison (ranking identity)."),
    width = 14, height = 11, config = FIG_CFG)
}

# =============================================================================
# FIG 3j-A -- Faceted top-TF lollipop for WT_heat (6 statistics)
# CLAIM: shows where HIF1a and IFN-axis TFs sit under every method for the
#        primary heat-response contrast; rank, not score, is the signal.
# =============================================================================

fig3j_df <- rd("fig3j_allmethods_topTF_data.csv")

make_fig3j(fig3j_df,
           chosen_contrast = "WT_heat",
           fig_label       = "Fig 3j-A",
           out_fname       = "fig3j_topTF_allmethods_WT_heat.pdf")

# =============================================================================
# FIG 3j-B -- Faceted top-TF lollipop for Interaction (6 statistics)
# CLAIM: shows whether the cGAS-dependence signal (Interaction) is consistent
#        across all methods; IFN-axis TFs should top every facet.
# =============================================================================

make_fig3j(fig3j_df,
           chosen_contrast = "Interaction",
           fig_label       = "Fig 3j-B",
           out_fname       = "fig3j_topTF_allmethods_Interaction.pdf")

# =============================================================================
# FIG 3k -- MLM is the structural outlier across BOTH axes (rank heatmap)
# CLAIM (audited against fig3k_method_rank_spearman.csv + divergence data):
#   MLM is the ONLY multivariate / de-confounding estimator; its full-ranking
#   Spearman vs ULM is 0.62, while the whole univariate family sits at 0.97-1.00.
#   That single structural difference reshuffles BOTH axes at once: among the
#   univariate top set, Rela #3->#350, Nfkb1 #6->#86, Hif1a #12->#142 under MLM.
#   So HIF1a's "collapse" is one INSTANCE of a general MLM reshuffle, NOT a
#   HIF-specific quirk -- and the IFN axis is demoted by MLM too, not stable.
# =============================================================================

cat("[3k] Building fig3k_method_rank_divergence.pdf ...\n")

div_df  <- rd("fig3k_method_rank_divergence_data.csv")
spe_df  <- rd("fig3k_method_rank_spearman.csv")

# -- Spearman annotation string per statistic (for axis labels) --
# Flag MLM as the lone multivariate estimator so the column reads as "different".
spe_rho <- setNames(spe_df$spearman_vs_ulm, spe_df$statistic)
spe_label <- setNames(
  vapply(spe_df$statistic, function(s) {
    tag <- if (s == "mlm") "  [multivariate]" else ""
    sprintf("%s%s\nrho=%.2f", STAT_LABELS[s], tag, spe_rho[s])
  }, character(1)),
  spe_df$statistic
)

# Focus on WT_heat for the main divergence panel
wt_div <- div_df[div_df$contrast == "WT_heat", ]
wt_div$statistic <- factor(wt_div$statistic, levels = STAT_ORDER,
                            labels = spe_label[STAT_ORDER])

# Order TFs by their mean rank across stats (most activated at top)
tf_mean_rank <- tapply(wt_div$rank, wt_div$source, mean)
tf_order     <- names(sort(tf_mean_rank))
wt_div$source <- factor(wt_div$source, levels = tf_order)

# Annotate HIF axis TFs
wt_div$tf_axis <- ifelse(wt_div$source %in% c("Hif1a", "Epas1"), "HIF",
                  ifelse(wt_div$source %in% c("Irf1","Irf3","Irf7",
                                               "Stat1","Stat2","Nfkb1","Rela"), "IFN",
                         "other"))

# x position of the MLM column (it is the 2nd factor level in STAT_ORDER) so we
# can box it -- the visual cue that makes the eye land on the outlier.
mlm_x <- which(STAT_ORDER == "mlm")

# -------- heatmap tile: rank as fill color --------------------------------
# PALETTE FIX (TASK 3): the old gradient2(midpoint=max_rank/2) put almost every
# focused TF below the midpoint, washing the whole grid to one hue so MLM's shift
# was invisible. Switch to a BINNED viridis scale on rank quantiles: top ranks
# pop bright, demoted ranks go dark -- so an MLM cell jumping from bright (top)
# to dark (demoted) is immediately obvious. (No statistics: just display bins.)
fig3k <- ggplot(wt_div, aes(x = statistic, y = source, fill = rank)) +
  geom_tile(color = "white", linewidth = 0.35) +
  # Boxed MLM column = the structural outlier the eye should land on.
  annotate("rect",
           xmin = mlm_x - 0.5, xmax = mlm_x + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = NA, color = "grey15", linewidth = 1.1) +
  geom_text(aes(label = rank,
                color = ifelse(rank <= 50, "white", "grey15")),
            size = 2.4) +
  scale_fill_stepsn(
    colours = c("#4B0C6B", "#781C6D", "#BB3754", "#ED6925", "#FCB519", "#FCFFA4"),
    breaks  = c(20, 50, 100, 200, 400),
    name    = "rank\n(lower = more\nactivated)"
  ) +
  scale_color_identity() +
  # Facet to add a side-strip for axis membership
  facet_grid(tf_axis ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(tf_axis = c(
               "HIF"   = "HIF axis",
               "IFN"   = "IFN/NFkB axis",
               "other" = "other TFs"
             ))) +
  labs(
    title    = "Fig 3k. MLM (the only multivariate estimator) reshuffles BOTH axes (WT_heat)",
    subtitle = paste0(
      "Rank of each TF under the 6 decoupleR statistics; boxed MLM column is the outlier.  ",
      "Spearman vs ULM (full 658-TF ranking) is 0.62 for MLM\n",
      "but 0.97-1.00 for every univariate statistic.  Under MLM, the univariate top set is demoted across ",
      "BOTH axes: Rela 3->350, Nfkb1 6->86, Hif1a 12->142.\n",
      "HIF1a's collapse is ONE instance of this general MLM de-confounding reshuffle, not a HIF-specific quirk."
    ),
    x = "decoupleR statistic (rho = Spearman rank-corr vs ULM, full WT_heat ranking)",
    y = NULL
  ) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x  = element_text(size = 8, angle = 10, hjust = 0.5, vjust = 1,
                                    lineheight = 0.88),
        axis.text.y  = element_text(size = 7.5),
        strip.text.y = element_text(angle = 0, face = "bold", size = 8),
        legend.key.height = unit(1.0, "cm"),
        panel.spacing.y   = unit(0.4, "lines"))

save_overview(
  fig3k, STAGE, "fig3k_method_rank_divergence",
  table = wt_div[, c("source", "statistic", "rank", "tf_axis")],
  finding = paste0(
    "MLM (the only multivariate estimator) reshuffles both axes on WT_heat: rank ",
    "heatmap shows MLM is a structural outlier with a Spearman rho of 0.62 vs ULM, ",
    "while univariate statistics sit at 0.97-1.00. Top univariate TFs are demoted ",
    "under MLM: Rela 3->350, Nfkb1 6->86, Hif1a 12->142. ", cap),
  script = SCRIPT, fn = "save_overview",
  config_kv = "figures.base_size=16; figures.base_size_column=9; statistics",
  input = "03_results/04_tf/tables/fig3k_method_rank_divergence_data.csv; fig3k_method_rank_spearman.csv",
  how_to_read = paste0(
    "Rank heatmap: rows = key TFs, columns = statistics. Fill = binned rank quantile ",
    "(bright yellow = highly activated, dark purple = demoted/inactive). Cell label = ",
    "absolute rank. The boxed column marks MLM (multivariate outlier). ",
    "Claim tier: methodological rank divergence."),
  width = 11, height = 9, config = FIG_CFG)

cat("[DONE] 03b VIZ complete -- 3 single-claim Phase-3b figures rendered from tidy tables.\n")
