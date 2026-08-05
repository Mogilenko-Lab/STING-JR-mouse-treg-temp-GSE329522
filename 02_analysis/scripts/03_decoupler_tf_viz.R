#!/usr/bin/env Rscript
# =============================================================================
# 03_decoupler_tf_viz.R: TF activity and the Hif1a rank-cascade forensics
# =============================================================================
# Eight single-claim figures for the GSE329522 iTreg design (genotype x temperature),
# each drawn from the tidy tables the compute script wrote. This script runs no
# statistics; it orders factors, facets, and writes label text.
#
# Reads 03_results/04_tf/tables/:
#   fig3a_tf_interaction_axes_data.csv    fig3e_mlm_collinearity_data.csv
#   fig3b_top_tf_by_contrast_data.csv     fig3e_mlm_collinearity_summary.csv
#   fig3c_hif1a_rank_cascade_data.csv     fig3e_score_collapse.csv
#   fig3d_regulon_swap_data.csv           fig3f_consensus_zmix_data.csv
#   fig3d_regulon_swap_summary.csv        fig3g_target_decomposition_data.csv
#   fig3i_interaction_primer_data.csv     fig3g_target_decomposition_summary.csv
#
# Writes, one save_overview() call per stem, under 03_results/04_tf/:
#   figures/_overview/<stem>.{pdf,png}   tables/_overview/<stem>.csv   README.md caption
#
# One claim per figure. fig3b is the single multi-panel, and its four facets share the
# ULM-score x axis.
#
# Run from the compartment root:
#   Rscript 02_analysis/scripts/03_decoupler_tf_viz.R
# =============================================================================

library(ggrepel)
source("02_analysis/helpers/figure_style.R")  # contract: FIG_CFG, project_theme, save_overview, purge_figures, ...
source("02_analysis/config/config.R")          # palette tokens: AXIS_COLORS/MODULE_COLORS/DIVERGING_COLORS
load_packages()  # brings in ggplot2 and dplyr

STAGE   <- "04_tf"
SCRIPT  <- "02_analysis/scripts/03_decoupler_tf_viz.R"
TBL_DIR <- stage_dir("04_tf", "tables")
TOP_N   <- as.integer(FIG_CFG$figures$top_n %||% 20L)

DOWN <- DIVERGING_COLORS$negative   # blue   (down / suppressed)
NEUT <- DIVERGING_COLORS$neutral    # white  (neutral)
UP   <- DIVERGING_COLORS$positive   # orange (up / activated)

rd <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)

# save_overview purges figures/_overview/ per stem, so the run owns its whole namespace
# once the flat-directory artifacts of an earlier layout are cleared here.
.FLAT_FIG_DIR <- stage_dir("04_tf", "figures")
.OWNED_STEMS  <- c("fig3a_tf_interaction_axes", "fig3b_top_tf_by_contrast",
                   "fig3c_hif1a_rank_cascade", "fig3d_regulon_swap",
                   "fig3e_mlm_collinearity", "fig3f_consensus_zmix",
                   "fig3g_target_decomposition",
                   "fig3i_interaction_primer", "fig3b_hif1a_robustness")
for (.stale in list.files(.FLAT_FIG_DIR, full.names = TRUE)) {
  if (any(startsWith(basename(.stale), .OWNED_STEMS)) &&
      grepl("\\.(pdf|png)$", .stale)) {
    file.remove(.stale)
    cat(sprintf("  removed stale flat-dir artifact %s\n", basename(.stale)))
  }
}

# Every stem here is a cross-method forensic overview, so all of them route through
# save_overview() to figures/_overview/. The screen .png is the clean raster for deck
# packaging, with the caption in the README.

# Fixed so the fig3d jitter lands in the same place on every run. Without it the panel
# redraws differently each time and a re-render reads as a changed result.
set.seed(42)

# =============================================================================
# FIG 3a: TF activity on the cGAS x heat Interaction contrast (CollecTRI ULM)
# =============================================================================
fig3a_df <- rd("fig3a_tf_interaction_axes_data.csv")
fig3a_df$sig  <- as.logical(fig3a_df$sig)
fig3a_df$axis <- factor(fig3a_df$axis, levels = c("HIF", "IFN", "other"))
fig3a_df <- fig3a_df[order(fig3a_df$score), ]
fig3a_df$source <- factor(fig3a_df$source, levels = fig3a_df$source)

axis_shapes <- c("HIF" = 17, "IFN" = 15, "other" = 16)
axis_sizes  <- c("HIF" = 3.4, "IFN" = 3.4, "other" = 1.9)

# Segments and points share one family colour, from AXIS_COLORS, so the HIF cluster near
# zero and the positive IFN/NFkB cluster separate at a glance.
# Drawn text (2026-08-05): the previous subtitle stated the finding and used
# "->" as prose. It now says only what the glyphs are; the flat/NS HIF axis and
# the interaction-positive IFN members are in the README caption below.
FIG3A_W     <- 11   # must match the width passed to save_overview()
FIG3A_TITLE <- "TF activity on the cGAS x heat Interaction contrast (CollecTRI ULM)"
FIG3A_SUBTITLE <- paste0(
  "One lollipop per TF; x = its ULM activity score on the Interaction contrast, zero rule marked.\n",
  "Colour and shape mark the TF axis: orange triangle HIF, blue square IFN/NFkB, grey other. * = BH padj < 0.05 (n=5).")
fits_canvas(FIG3A_TITLE,    FIG_CFG$figures$title_size    %||% 16, "bold",  FIG3A_W, "fig3a title")
fits_canvas(FIG3A_SUBTITLE, FIG_CFG$figures$subtitle_size %||% 11, "plain", FIG3A_W, "fig3a subtitle")

fig3a <- ggplot(fig3a_df, aes(x = score, y = source)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_segment(aes(xend = 0, yend = source, color = axis), linewidth = 0.7) +
  geom_point(aes(color = axis, shape = axis, size = axis)) +
  geom_text(data = subset(fig3a_df, sig & axis != "other"),
            aes(label = "*"), nudge_x = 0.28, size = 5, color = "black") +
  scale_color_manual(values = AXIS_COLORS, name = "TF axis",
                     labels = AXIS_LABELS) +
  scale_shape_manual(values = axis_shapes, name = "TF axis",
                     labels = AXIS_LABELS) +
  scale_size_manual(values = axis_sizes, guide = "none") +
  guides(color = guide_legend(override.aes = list(size = 3, shape = c(17, 15, 16))),
         shape = "none") +
  labs(
    title = FIG3A_TITLE,
    subtitle = FIG3A_SUBTITLE,
    x = "TF activity (ULM score; Interaction = WT_heat - KO_heat)", y = NULL
  ) +
  project_theme(config = FIG_CFG)

save_overview(
  fig3a, STAGE, "fig3a_tf_interaction_axes",
  table     = fig3a_df,
  finding   = paste0("On the cGAS x heat Interaction contrast, the HIF axis (Hif1a/Epas1) is ",
                     "flat/NS -- no detectable cGAS-dependence -- while IFN members ",
                     "Irf3/Stat2/Stat1 are interaction-positive."),
  script    = SCRIPT, fn = "save_overview",
  config_kv = "colors.diverging; design.axis_colors",
  input     = "03_results/04_tf/tables/fig3a_tf_interaction_axes_data.csv",
  how_to_read = paste0(
    "Lollipop x = TF activity (ULM score; Interaction = WT_heat - KO_heat); ",
    "color/shape = TF axis (orange triangle = HIF; blue square = IFN/NFkB; grey = other). ",
    "Glyph: * = BH padj < 0.05. Asymmetry, not proven independence; claim tier L3 (n=5)."),
  config    = FIG_CFG, width = 11, height = 9)

# =============================================================================
# FIG 3b: top activated and suppressed TFs per contrast (CollecTRI ULM)
# =============================================================================
# Four comparable lollipop facets on one shared ULM-score x axis, which is what makes a
# facet_wrap the right frame for what each contrast captures.
fig3b_df <- rd("fig3b_top_tf_by_contrast_data.csv")
# Cap to figures.top_n per contrast and direction so the facets hold the font floor with
# every label intact.
fig3b_df <- do.call(rbind, lapply(
  split(fig3b_df, list(fig3b_df$contrast, fig3b_df$direction), drop = TRUE),
  function(g) g[order(-abs(g$score)), ][seq_len(min(nrow(g), TOP_N)), ]))
fig3b_df$contrast <- factor(fig3b_df$contrast,
                            levels = c("WT_heat", "KO_heat", "Temp_main", "Interaction"))
fig3b_df$axis     <- factor(fig3b_df$axis, levels = c("HIF", "IFN", "other"))
fig3b_df$key_tf   <- as.logical(fig3b_df$key_tf)
# Order within each facet by score, carried on a row key rather than through tidytext.
fig3b_df <- fig3b_df[order(fig3b_df$contrast, fig3b_df$score), ]
fig3b_df$row_key <- factor(seq_len(nrow(fig3b_df)),
                           levels = seq_len(nrow(fig3b_df)),
                           labels = fig3b_df$source)

fig3b <- ggplot(fig3b_df, aes(x = score, y = row_key)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_segment(aes(xend = 0, yend = row_key, color = axis), linewidth = 0.6) +
  geom_point(aes(color = axis, size = key_tf)) +
  geom_text(data = subset(fig3b_df, key_tf),
            aes(label = source), nudge_x = 0.15, hjust = 0, size = 2.4, color = "black") +
  facet_wrap(~ contrast, scales = "free_y", ncol = 2,
             labeller = labeller(contrast = function(x) contrast_label(x))) +
  scale_color_manual(values = AXIS_COLORS, name = "TF axis", labels = AXIS_LABELS) +
  scale_size_manual(values = c("TRUE" = 2.8, "FALSE" = 1.5), guide = "none") +
  labs(
    title = "Top activated and suppressed TFs per contrast (CollecTRI ULM)",
    subtitle = "Shared x-axis (ULM score) across all four facets. Key TFs (HIF orange, IFN/NFkB blue) labelled.\nReads out what each contrast captures: heat arms vs the temperature main effect vs the cGAS x heat Interaction.",
    x = "TF activity (ULM score)", y = NULL
  ) +
  project_theme(config = FIG_CFG)

save_overview(
  fig3b, STAGE, "fig3b_top_tf_by_contrast",
  table     = fig3b_df,
  finding   = paste0("Top activated/suppressed TFs (CollecTRI ULM) per contrast on a shared ",
                     "ULM-score x-axis: reads out what each heat arm, the temperature main ",
                     "effect, and the cGAS x heat Interaction capture."),
  script    = SCRIPT, fn = "save_overview",
  config_kv = paste0("figures.top_n=", TOP_N, "; design.axis_colors"),
  input     = "03_results/04_tf/tables/fig3b_top_tf_by_contrast_data.csv",
  how_to_read = paste0(
    "Four facets share the x-axis (TF activity, ULM score). Lollipop color = TF axis ",
    "(orange = HIF; blue = IFN/NFkB; grey = other); HIF/IFN watchlist TFs labelled. ",
    "Capped to top-", TOP_N, " per direction per contrast by |score|. Claim tier L3 (n=5)."),
  config    = FIG_CFG, width = 16, height = 14, wide = TRUE)

# =============================================================================
# FIG 3c: Hif1a's rank across four method and network configurations
# =============================================================================
# The cascade this figure and the three that follow decompose: #1 DoRothEA -> #12 on the
# CollecTRI swap -> #142 under MLM -> #8 under the consensus.
fig3c_df <- rd("fig3c_hif1a_rank_cascade_data.csv")
fig3c_df <- fig3c_df[order(fig3c_df$method_order), ]
fig3c_df$method <- factor(fig3c_df$method, levels = fig3c_df$method)
fig3c_df$lab    <- sprintf("#%d / %d\nscore=%.2f", fig3c_df$rank, fig3c_df$n_tfs, fig3c_df$score)

fig3c <- ggplot(fig3c_df, aes(x = method, y = rank, color = score, group = 1)) +
  geom_line(color = "grey55", linewidth = 0.6) +
  geom_point(size = 5) +
  geom_text(aes(label = lab), vjust = -0.6, size = 3, color = "grey15", lineheight = 0.9) +
  scale_color_gradient2(low = DOWN, mid = NEUT, high = UP, midpoint = 0, name = "ULM/consensus\nscore") +
  scale_y_reverse(expand = expansion(mult = c(0.10, 0.22))) +
  labs(
    title = sprintf("Hif1a rank across inference method and network, %s",
                    contrast_label("WT_heat")),
    subtitle = "Rank 1 is the most activated. Each point is labelled with the rank and the score behind it.",
    x = NULL, y = "Hif1a rank among activated TFs (lower = higher activity)"
  ) +
  project_theme(config = FIG_CFG)

save_overview(
  fig3c, STAGE, "fig3c_hif1a_rank_cascade",
  table     = fig3c_df,
  finding   = paste0("Hif1a's #1 DoRothEA ranking is method/network-fragile: it collapses ",
                     "under the CollecTRI swap (#12) and multivariate MLM de-confounding ",
                     "(#142), then is partially re-inflated by the consensus (#8)."),
  script    = SCRIPT, fn = "save_overview",
  config_kv = "colors.diverging",
  input     = "03_results/04_tf/tables/fig3c_hif1a_rank_cascade_data.csv",
  how_to_read = paste0(
    "x = inference method/network config (left -> right reading order); y = Hif1a rank ",
    "among activated TFs (reversed; rank #1 at top = most activated). Point color = ULM/",
    "consensus score (orange up / blue down). Labels give rank/total and score. Tier L3 (n=5)."),
  config    = FIG_CFG, width = 11, height = 9)

# =============================================================================
# FIG 3d: the regulon swap behind #1 -> #12 (DoRothEA -> CollecTRI)
# =============================================================================
# One axis, the per-target t_wt distribution by membership; the regulon sizes and overlap
# counts ride as an annotation.
fig3d_df  <- rd("fig3d_regulon_swap_data.csv")
fig3d_sum <- rd("fig3d_regulon_swap_summary.csv")
fig3d_df$membership <- factor(fig3d_df$membership,
                              levels = c("dorothea_only", "both", "collectri_only"),
                              labels = c("DoRothEA-only", "shared (both)", "CollecTRI-only"))

# Membership palette from MODULE_COLORS tokens: teal DoRothEA-only, grey shared, and the
# stress purple for the CollecTRI-only targets the swap adds.
mem_cols <- c(
  "DoRothEA-only"  = unname(MODULE_COLORS["hif1a_hypoxic_core"]),
  "shared (both)"  = unname(MODULE_COLORS["other_unclassified"]),
  "CollecTRI-only" = unname(MODULE_COLORS["heatshock_stress"]))

dor <- fig3d_sum[fig3d_sum$network == "DoRothEA", ]
ct  <- fig3d_sum[fig3d_sum$network == "CollecTRI", ]
swap_note <- sprintf(
  "DoRothEA regulon n=%d (mean|t|=%.2f) -> CollecTRI regulon n=%d (mean|t|=%.2f)\n%d shared; %d DoRothEA-only; %d CollecTRI-only (high |t|, mixed sign) targets added",
  dor$n, dor$mean_abs_t, ct$n, ct$mean_abs_t,
  dor$n_intersection[1], dor$n_dorothea_only[1], dor$n_collectri_only[1])

# Drawn text (2026-08-05): the previous subtitle carried the whole dilution
# mechanism -- mixed-sign CollecTRI-only targets cancelling and dragging the
# mor-weighted mean down. That is the explanation for the #1 -> #12 step and it
# is already stated in the README caption below, so the subtitle keeps only the
# encoding.
FIG3D_TITLE <- "Per-target heat t-statistic by Hif1a regulon membership"
FIG3D_SUBTITLE <- paste0(
  "Box and jitter of the per-target WT_heat t-statistic, split by which network lists the target in\n",
  "Hif1a's regulon. The label gives the three regulon sizes behind the DoRothEA-to-CollecTRI swap.")
fits_canvas(FIG3D_TITLE,    FIG_CFG$figures$title_size    %||% 16, "bold",  11, "fig3d title")
fits_canvas(FIG3D_SUBTITLE, FIG_CFG$figures$subtitle_size %||% 11, "plain", 11, "fig3d subtitle")

fig3d <- ggplot(fig3d_df, aes(x = membership, y = t_wt, fill = membership)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.16, height = 0, size = 1, alpha = 0.5, color = "grey20") +
  scale_fill_manual(values = mem_cols, guide = "none") +
  annotate("label", x = 0.6, y = max(fig3d_df$t_wt) * 0.96, hjust = 0, vjust = 1,
           label = swap_note, size = 2.7, fill = "grey96", color = "grey15") +
  labs(
    title = FIG3D_TITLE,
    subtitle = FIG3D_SUBTITLE,
    x = "Hif1a target membership across networks", y = "WT_heat target t-statistic"
  ) +
  project_theme(config = FIG_CFG)

save_overview(
  fig3d, STAGE, "fig3d_regulon_swap",
  table     = fig3d_df,
  finding   = paste0("The DoRothEA -> CollecTRI regulon swap dilutes Hif1a's ULM signal: the ",
                     "CollecTRI-only targets carry high |t| but mixed sign, so they cancel and ",
                     "drag the regulon mean down (explains the #1 -> #12 drop)."),
  script    = SCRIPT, fn = "save_overview",
  config_kv = "design.module_colors",
  input     = "03_results/04_tf/tables/fig3d_regulon_swap_data.csv",
  how_to_read = paste0(
    "Box+jitter of per-target WT_heat t-statistic by regulon membership ",
    "(DoRothEA-only / shared / CollecTRI-only). The CollecTRI-only box straddles 0 ",
    "(mixed sign), diluting the mor-weighted mean. Annotation gives regulon sizes. Tier L3 (n=5)."),
  config    = FIG_CFG, width = 11, height = 9)

# =============================================================================
# FIG 3e: target collinearity behind #12 -> #142 (ULM -> MLM)
# =============================================================================
# One axis, the count of other TFs sharing each Hif1a target; the ULM-to-MLM score
# collapse rides as an annotation.
fig3e_df  <- rd("fig3e_mlm_collinearity_data.csv")
fig3e_sc  <- rd("fig3e_score_collapse.csv")
fig3e_sum <- rd("fig3e_mlm_collinearity_summary.csv")   # aggregates from the compute stage
fig3e_df$co_regulated <- as.logical(fig3e_df$co_regulated)

pct_co  <- 100 * fig3e_sum$frac_co_regulated
mean_ot <- fig3e_sum$mean_other_TFs
ulm <- fig3e_sc[fig3e_sc$method == "ULM", ]
mlm <- fig3e_sc[fig3e_sc$method == "MLM", ]
collapse_note <- sprintf(
  "Hif1a ULM score %.2f (rank #%d)\n  -> MLM score %.2f (rank #%d)\n%.0f%% of targets co-regulated; mean %.0f other TFs/target",
  ulm$score, ulm$rank, mlm$score, mlm$rank, pct_co, mean_ot)

# Drawn text (2026-08-05): the previous subtitle carried the load-bearing
# sentence of the whole rank cascade -- nearly every Hif1a target is
# co-regulated, so the multivariate MLM spreads its signal across the regulons
# that share it. That sentence is preserved in the README caption below, which
# is where the contract asks for it; the subtitle keeps only the encoding.
FIG3E_TITLE <- "Hif1a targets by how many other CollecTRI regulons claim them"
FIG3E_SUBTITLE <- paste0(
  "Histogram of Hif1a's CollecTRI targets by the number of other CollecTRI regulons that also list\n",
  "each target. The dashed line marks the mean; the label gives Hif1a's ULM and MLM scores and ranks.")
fits_canvas(FIG3E_TITLE,    FIG_CFG$figures$title_size    %||% 16, "bold",  11, "fig3e title")
fits_canvas(FIG3E_SUBTITLE, FIG_CFG$figures$subtitle_size %||% 11, "plain", 11, "fig3e subtitle")

fig3e <- ggplot(fig3e_df, aes(x = n_other_TFs)) +
  geom_histogram(binwidth = 10, fill = DOWN, color = "white", boundary = 0) +
  geom_vline(xintercept = mean_ot, color = UP, linewidth = 0.8, linetype = "22") +
  annotate("label", x = max(fig3e_df$n_other_TFs) * 0.98, y = Inf, hjust = 1, vjust = 1.1,
           label = collapse_note, size = 2.9, fill = "grey96", color = "grey15") +
  labs(
    title = FIG3E_TITLE,
    subtitle = FIG3E_SUBTITLE,
    x = "number of other CollecTRI TFs sharing each Hif1a target", y = "Hif1a targets (count)"
  ) +
  project_theme(config = FIG_CFG)

save_overview(
  fig3e, STAGE, "fig3e_mlm_collinearity",
  table     = fig3e_df,
  finding   = paste0("Nearly every Hif1a target is co-regulated by many other CollecTRI TFs, so ",
                     "the multivariate MLM re-attributes Hif1a's signal away from it -- collapsing ",
                     "its score (explains the #12 -> #142 drop)."),
  script    = SCRIPT, fn = "save_overview",
  config_kv = "colors.diverging",
  input     = "03_results/04_tf/tables/fig3e_mlm_collinearity_data.csv",
  how_to_read = paste0(
    "Histogram: for each Hif1a target, how many OTHER CollecTRI TFs also regulate it. ",
    "Dashed orange line = mean. Annotation gives the ULM -> MLM score/rank collapse. ",
    "High co-regulation is why MLM de-confounds Hif1a's signal away. Tier L3 (n=5)."),
  config    = FIG_CFG, width = 11, height = 8)

# =============================================================================
# FIG 3f: the consensus z-mix behind #142 -> #8 (MLM -> consensus)
# =============================================================================
fig3f_df <- rd("fig3f_consensus_zmix_data.csv")
stat_order <- c("ulm", "mlm", "wsum", "norm_wsum", "corr_wsum", "consensus")
fig3f_df$statistic <- factor(fig3f_df$statistic, levels = stat_order)
fig3f_df$family    <- factor(fig3f_df$family,
                             levels = c("univariate", "multivariate", "consensus"))
fig3f_df <- fig3f_df[order(fig3f_df$statistic), ]

consensus_z <- fig3f_df$z_score[fig3f_df$statistic == "consensus"]
fam_cols <- c("univariate" = UP, "multivariate" = DOWN, "consensus" = "grey25")

fig3f <- ggplot(fig3f_df, aes(x = statistic, y = z_score, fill = family)) +
  geom_hline(yintercept = consensus_z, color = "grey40", linetype = "22", linewidth = 0.6) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.2f", z_score)), vjust = -0.4, size = 3, color = "grey20") +
  annotate("text", x = 0.6, y = consensus_z, hjust = 0, vjust = -0.4,
           label = sprintf("consensus mean z = %.2f", consensus_z), size = 2.9, color = "grey35") +
  scale_fill_manual(values = fam_cols, name = "statistic family") +
  labs(
    title = "Hif1a folded z per decoupleR statistic, against the consensus mean",
    subtitle = "The step from rank 142 back to rank 8. The four univariate statistics stay high and mlm alone is low,\nand the consensus is their mean (dashed line), which carries Hif1a back up.",
    x = "decoupleR statistic", y = "Hif1a folded z-score"
  ) +
  project_theme(config = FIG_CFG)

save_overview(
  fig3f, STAGE, "fig3f_consensus_zmix",
  table     = fig3f_df,
  finding   = paste0("Only the multivariate MLM gives Hif1a a low folded-z; the four univariate ",
                     "statistics stay high and the consensus (their mean) re-inflates Hif1a back ",
                     "to #8 (explains the #142 -> #8 recovery)."),
  script    = SCRIPT, fn = "save_overview",
  config_kv = "colors.diverging",
  input     = "03_results/04_tf/tables/fig3f_consensus_zmix_data.csv",
  how_to_read = paste0(
    "Bars = Hif1a folded z per decoupleR statistic; fill = statistic family ",
    "(orange univariate / blue multivariate / grey consensus). Dashed line = consensus ",
    "mean z. The lone low MLM bar is outvoted by the univariate family. Tier L3 (n=5)."),
  config    = FIG_CFG, width = 11, height = 8)

# =============================================================================
# FIG 3g: where Hif1a's WT_heat score comes from, over all 353 regulon members
# =============================================================================
# What this regulon reads out is a heat-induced glycolytic and stress program that
# partially overlaps HIF targets, so every label here names the class rather than "the HIF
# program".
fig3g_df  <- rd("fig3g_target_decomposition_data.csv")
fig3g_sum <- rd("fig3g_target_decomposition_summary.csv")

# The stress genes that carry the right tail of the 'other' class.
STRESS_GENES <- c("Timp1", "Sdc1", "Cdkn1a", "Serpine1", "Eno2", "Hspa1a", "Spp1")

# Labelled on the face: the top stress drivers, and all seven HIF-specific members.
LABEL_GENES <- c(
  "Timp1", "Sdc1", "Cdkn1a", "Serpine1", "Eno2", "Hspa1a",
  "Vegfa", "Egln3", "Slc2a1", "Car9", "Bnip3", "Bnip3l", "Pdk1"
)

# Ascending contribution gives the left-to-right waterfall order.
fig3g_df <- fig3g_df[order(fig3g_df$contrib), ]
fig3g_df$rank <- seq_len(nrow(fig3g_df))

# Classes come from the compute table and share fig3l's vocabulary. A stress gene overrides
# to the heat-shock hue; the rest of 'other' stays grey.
fig3g_df$display_color <- ifelse(
  fig3g_df$class == "hypoxic HIF core", MODULE_COLORS["hif1a_hypoxic_core"],   # repressed-core hue
  ifelse(
    fig3g_df$class == "shared/glycolytic", MODULE_COLORS["shared_angio_glucose"],  # light orange (UP)
    ifelse(
      fig3g_df$target %in% STRESS_GENES, MODULE_COLORS["heatshock_stress"],        # purple
      MODULE_COLORS["other_unclassified"]                                            # grey80
    )
  )
)

# The four display categories, keyed for scale_color_manual.
FIG3G_COLS <- c(
  "hypoxic HIF core (repressed)" = unname(MODULE_COLORS["hif1a_hypoxic_core"]),
  "shared/glycolytic"            = unname(MODULE_COLORS["shared_angio_glucose"]),
  "heat-shock/stress"            = unname(MODULE_COLORS["heatshock_stress"]),
  "other (unclassified)"         = unname(MODULE_COLORS["other_unclassified"])
)

fig3g_df$legend_class <- ifelse(
  fig3g_df$class == "hypoxic HIF core",    "hypoxic HIF core (repressed)",
  ifelse(
    fig3g_df$class == "shared/glycolytic", "shared/glycolytic",
    ifelse(
      fig3g_df$target %in% STRESS_GENES, "heat-shock/stress",
      "other (unclassified)"
    )
  )
)

# The labelled genes, which save_overview writes below as this figure's same-stem table.
fig3g_labels_df <- fig3g_df[fig3g_df$target %in% LABEL_GENES,
                             c("source","target","mor","t_wt","contrib","class","rank","legend_class")]

fig3g_lab_sub <- fig3g_df[fig3g_df$target %in% LABEL_GENES, ]

# Annotation text, read from the summary table.
get_pct_g  <- function(cl) fig3g_sum$pct_of_total[fig3g_sum$class == cl]
get_sum_g  <- function(cl) fig3g_sum$sum_contrib[fig3g_sum$class == cl]
annot_text <- sprintf(
  "Score provenance (%%  of aggregate signal):\n  other (unclassified):    %.1f%%\n  shared/glycolytic:        %.1f%%\n  hypoxic HIF core:         %.1f%% (below zero)",
  get_pct_g("other"), get_pct_g("shared/glycolytic"), get_pct_g("hypoxic HIF core"))

# A coarse three-class view of the regulon fig3l decomposes finer, which is why the two
# figures report different percentages for the same genes.
sub_g_expanded <- paste0(
  "All 353 Hif1a regulon members, ranked by signed contribution (sign(mor) x t_wt); colour marks the coarse regulon class, and\n",
  sprintf("the four members of the hypoxic HIF core (Pdk1/Bnip3/Bnip3l/Car9, teal) sum to %+.2f.\n",
          get_sum_g("hypoxic HIF core")),
  sprintf("Class shares of the aggregate signal: %.1f%% other | %.1f%% shared/glycolytic | %.1f%% hypoxic HIF core (fig3l splits 'other' finer).",
          get_pct_g("other"), get_pct_g("shared/glycolytic"), get_pct_g("hypoxic HIF core"))
)
# fig3g is saved at width = 16 (wide = TRUE); guard the drawn lines against it.
fits_canvas(sub_g_expanded, FIG_CFG$figures$subtitle_size %||% 11, "plain", 16, "fig3g subtitle")

fig3g <- ggplot(fig3g_df, aes(x = rank, y = contrib, color = legend_class)) +
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  geom_segment(aes(xend = rank, yend = 0), linewidth = 0.35, alpha = 0.7) +
  # The three named classes get the larger head so they read against 340 grey stems.
  geom_point(aes(size = legend_class %in% c("hypoxic HIF core (repressed)", "heat-shock/stress", "shared/glycolytic")),
             alpha = 0.85) +
  geom_text_repel(data = fig3g_lab_sub,
            aes(x = rank, y = contrib, label = target),
            size = 2.6, fontface = "bold", color = "grey10",
            max.overlaps  = Inf,
            min.segment.length = 0,
            box.padding   = 0.45,
            point.padding = 0.25,
            segment.size  = 0.3,
            segment.color = "grey45",
            direction     = "both",
            seed          = 42) +
  annotate("label",
           x = 1, y = max(fig3g_df$contrib) * 0.95,
           label = annot_text,
           hjust = 0, vjust = 1,
           size = 2.5, fill = "grey97", color = "grey20") +
  scale_color_manual(
    values = FIG3G_COLS,
    name   = "regulon class",
    guide  = guide_legend(override.aes = list(size = 3, alpha = 1))
  ) +
  scale_size_manual(values = c("TRUE" = 2.2, "FALSE" = 0.8), guide = "none") +
  scale_x_continuous(breaks = NULL) +
  labs(
    title    = sprintf("Hif1a target contribution across the ranked regulon, %s",
                       contrast_label("WT_heat")),
    subtitle = sub_g_expanded,
    x        = "regulon members ranked by contribution (left = most negative; right = most positive)",
    y        = "contribution to WT_heat ULM score (sign(mor) x t_wt)"
  ) +
  project_theme(config = FIG_CFG)

save_overview(
  fig3g, STAGE, "fig3g_target_decomposition",
  table     = fig3g_labels_df,
  finding   = paste0("Hif1a's positive WT_heat ULM score is built almost entirely from its ",
                     "'other' regulon members -- generic stress/ECM genes tower over the HIF-",
                     "specific core, which is repressed. This is a heat-induced glycolytic/stress ",
                     "program partially overlapping HIF targets, not the canonical HIF program."),
  script    = SCRIPT, fn = "save_overview",
  config_kv = "design.module_colors",
  input     = "03_results/04_tf/tables/fig3g_target_decomposition_data.csv",
  how_to_read = paste0(
    "Lollipop landscape: all 353 Hif1a regulon members ranked by signed contribution ",
    "(sign(mor) x t_wt; + = pushes the ULM score up). Color = regulon class (purple heat-",
    "shock/stress; orange shared/glycolytic; teal repressed hypoxic core; grey other). ",
    "Labelled = top stress drivers + the 7 HIF-specific members (the source table). Tier L3 (n=5)."),
  config    = FIG_CFG, width = 16, height = 9, wide = TRUE)

# =============================================================================
# FIG 3i: reading the cGAS x heat Interaction as a difference of heat slopes
# =============================================================================
# Interaction = (WT_39 - WT_37) - (cGASKO_39 - cGASKO_37). Two facets on free y, since the
# two genes sit at different log2 CPM levels while sharing the metric.
fig3i_df <- rd("fig3i_interaction_primer_data.csv")
fig3i_df$gene     <- factor(fig3i_df$gene, levels = c("Ifit1", "Vegfa"))
fig3i_df$temp     <- factor(fig3i_df$temp, levels = c("37", "39"))
fig3i_df$genotype <- factor(fig3i_df$genotype, levels = c("WT", "cGASKO"))

# The contrast values ride in the facet strip through a labeller, which keeps them clear of
# the plotted slopes. Every number is read from the pre-computed contrast columns.
strip_labels <- vapply(c("Ifit1", "Vegfa"), function(g) {
  r <- fig3i_df[fig3i_df$gene == g, ][1, ]
  arm_tag <- if (g == "Ifit1") "IFN arm" else "shared HIF/glycolytic target"
  sprintf(
    "%s  [%s]\nWT heat %+.2f  /  KO heat %+.2f\nInteraction %+.2f (adjP %.3f)",
    g, arm_tag, r$wt_heat, r$ko_heat, r$interaction, r$inter_adjP
  )
}, character(1))

# Genotype palette from config tokens: WT teal, cGASKO the house up-orange.
GENO_COLS <- c("WT" = unname(MODULE_COLORS["hif1a_hypoxic_core"]), "cGASKO" = UP)

fig3i <- ggplot(fig3i_df, aes(x = temp, y = mean_log2cpm,
                               color = genotype, group = genotype,
                               linetype = genotype)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 3.2) +
  facet_wrap(~ gene, scales = "free_y", ncol = 2,
             labeller = labeller(gene = strip_labels)) +
  scale_color_manual(values = GENO_COLS, name = "Genotype") +
  scale_linetype_manual(values = c("WT" = "solid", "cGASKO" = "dashed"), name = "Genotype") +
  scale_x_discrete(labels = c("37" = "37 C", "39" = "39 C")) +
  # Headroom so slopes never crowd the strip band above.
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(
    title    = "The cGAS x heat Interaction as the difference of the two genotype heat-slopes",
    subtitle = paste0(
      "Interaction = (WT_39 - WT_37) - (KO_39 - KO_37).  WT solid, cGAS-KO dashed.\n",
      "Diverging slopes = cGAS-dependent (Ifit1, IFN arm).  ",
      "Parallel slopes = no detectable cGAS-dependence (Vegfa, shared HIF/glycolytic target)."
    ),
    x = "Temperature", y = "mean log2 CPM (group mean)"
  ) +
  project_theme(config = FIG_CFG)

save_overview(
  fig3i, STAGE, "fig3i_interaction_primer",
  table     = fig3i_df,
  finding   = paste0("Reading the cGAS x heat Interaction as the difference of genotype heat-",
                     "slopes: Ifit1 (IFN arm) slopes diverge -> large positive Interaction ",
                     "(cGAS-dependent); Vegfa (shared HIF/glycolytic target) slopes stay parallel ",
                     "-> Interaction ~0 (no detectable cGAS-dependence)."),
  script    = SCRIPT, fn = "save_overview",
  config_kv = "design.module_colors; colors.diverging",
  input     = "03_results/04_tf/tables/fig3i_interaction_primer_data.csv",
  how_to_read = paste0(
    "Two facets (free y; same metric log2 CPM): x = temperature, lines = genotype ",
    "(WT teal solid, cGAS-KO orange dashed). Diverging slopes = cGAS-dependent; ",
    "parallel slopes = no detectable cGAS-dependence. Strip text carries the contrast ",
    "values (WT heat / KO heat / Interaction + adjP). Tier L3 (n=5)."),
  config    = FIG_CFG, width = 12, height = 8)

cat("[DONE] 03_decoupler_tf_viz.R -- 8 single-claim Phase-3 figures rendered from tidy tables.\n")
