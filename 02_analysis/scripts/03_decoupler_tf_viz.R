#!/usr/bin/env Rscript
# =============================================================================
# 03_decoupler_tf_viz.R  --  PHASE 3 VIZ: TF activity + Hif1a rank-cascade forensics
# =============================================================================
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2: genotype x temperature)
# Phase:   3 (stage 04_tf)
# Role:    VISUALIZE half of the "normalize-then-visualize" split. Reads ONLY the
#          plot-ready tidy tables emitted by the Phase-3 COMPUTE script (Agent 2)
#          and renders ONE single-claim figure per file. Performs NO statistics
#          (no run_ulm/run_mlm/run_consensus, no p.adjust/prcomp/cor/decouple); it
#          only plots already-computed columns, doing cosmetic reshaping (factor
#          ordering, faceting, label text). Runs STANDALONE after compute.
#
# Inputs (03_results/04_tf/tables/, written by the compute script):
#   fig3a_tf_interaction_axes_data.csv        fig3e_mlm_collinearity_data.csv
#   fig3b_top_tf_by_contrast_data.csv         fig3e_score_collapse.csv
#   fig3c_hif1a_rank_cascade_data.csv         fig3f_consensus_zmix_data.csv
#   fig3d_regulon_swap_data.csv               fig3g_target_decomposition_data.csv
#   fig3d_regulon_swap_summary.csv            fig3g_target_decomposition_summary.csv
#                                             fig3h_lombardi_vs_phylo_data.csv
#
# Outputs (03_results/04_tf/figures/), one ggsave each:
#   fig3a_tf_interaction_axes.pdf    fig3e_mlm_collinearity.pdf
#   fig3b_top_tf_by_contrast.pdf     fig3f_consensus_zmix.pdf
#   fig3c_hif1a_rank_cascade.pdf     fig3g_target_decomposition.pdf
#   fig3d_regulon_swap.pdf           fig3h_lombardi_vs_phylo.pdf
#
# Figure discipline (AGENTS.md): ONE claim = ONE dedicated figure. The ONLY
# multi-panel here is fig3b, and ONLY because all 4 facets share the SAME x-axis
# (ULM score). Every figure carries the provisional_caption() stamp.
# Dependencies: config.R; ggplot2, dplyr
# =============================================================================

library(ggrepel)
source("02_analysis/config/config.R")
load_packages()  # limma/edgeR pulled in but unused; ggplot2 + dplyr are what we need

TBL_DIR <- stage_dir("04_tf", "tables")
FIG_DIR <- stage_dir("04_tf", "figures")

DOWN <- DIVERGING_COLORS$negative   # blue   (down / suppressed)
NEUT <- DIVERGING_COLORS$neutral    # white  (neutral)
UP   <- DIVERGING_COLORS$positive   # orange (up / activated)

cap <- provisional_caption()

rd <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# DECK EXPORT (TASK 4): when DECK_EXPORT is set, ALSO write a clean 300-dpi PNG
# (no embargo caption line) to 03_results/04_tf/deck_assets/ for pptx packaging.
# Canonical stamped PDFs are ALWAYS written, unchanged. The clean PNG keeps all
# real scientific subtitles/annotations -- only labs(caption=NULL) is dropped.
# -----------------------------------------------------------------------------
DECK_EXPORT <- nzchar(Sys.getenv("DECK_EXPORT"))
DECK_DIR    <- file.path(DIR_RESULTS, "04_tf", "deck_assets")
if (DECK_EXPORT) dir.create(DECK_DIR, recursive = TRUE, showWarnings = FALSE)

#' Write the canonical stamped PDF, and (when DECK_EXPORT) a clean 16:9 PNG.
#' @param plot   ggplot object carrying the provisional caption.
#' @param fname  PDF filename, e.g. "fig3a_tf_interaction_axes.pdf".
#' @param width,height  PDF dimensions (canonical, unchanged).
#' @param deck_h  PNG height in inches (PNG width fixed at 10in, ~16:9-friendly).
save_fig <- function(plot, fname, width, height, deck_h = 5.625) {
  ggsave(file.path(FIG_DIR, fname), plot, width = width, height = height)
  cat(sprintf("  wrote %s\n", fname))
  if (DECK_EXPORT) {
    png_name <- sub("\\.pdf$", ".png", fname)
    ggsave(file.path(DECK_DIR, png_name), plot + labs(caption = NULL),
           width = 10, height = deck_h, dpi = 300)
    cat(sprintf("  wrote deck_assets/%s (clean, no embargo stamp)\n", png_name))
  }
}

base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey25"),
        plot.caption  = element_text(size = 7, color = "grey45", hjust = 0, face = "italic"))

# Retire the illegible composite that this script replaces.
old_composite <- file.path(FIG_DIR, "fig3b_hif1a_robustness.pdf")
if (file.exists(old_composite)) {
  file.remove(old_composite)
  cat(sprintf("  removed retired composite %s\n", old_composite))
}

# =============================================================================
# FIG 3a -- TF activity on the cGAS x heat Interaction contrast
# CLAIM: HIF axis (Hif1a/Epas1) flat/NS; IFN/NFkB axis (Irf3/Stat1/Stat2)
#        interaction-positive.
# =============================================================================
fig3a_df <- rd("fig3a_tf_interaction_axes_data.csv")
fig3a_df$sig  <- as.logical(fig3a_df$sig)
fig3a_df$axis <- factor(fig3a_df$axis, levels = c("HIF", "IFN", "other"))
fig3a_df <- fig3a_df[order(fig3a_df$score), ]
fig3a_df$source <- factor(fig3a_df$source, levels = fig3a_df$source)

axis_shapes <- c("HIF" = 17, "IFN" = 15, "other" = 16)
axis_sizes  <- c("HIF" = 3.4, "IFN" = 3.4, "other" = 1.9)

# Color by TF axis family (single source of truth: AXIS_COLORS from config.R).
# Segments and points share the same family color so the HIF cluster (orange,
# near zero) vs IFN/NFkB cluster (blue, positive) is visible at a glance.
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
    title = "Fig 3a. TF activity on the cGAS x heat Interaction contrast (CollecTRI ULM)",
    subtitle = "HIF axis (orange triangles) is flat/NS -> no detectable cGAS-dependence; the IFN members Irf3/Stat2/Stat1 are interaction-positive\n(* = BH padj < 0.05). NFkB members (Nfkb1/Rela) lean positive but are NS at n=5.",
    x = "TF activity (ULM score; Interaction = WT_heat - KO_heat)", y = NULL,
    caption = cap
  ) +
  base_theme

save_fig(fig3a, "fig3a_tf_interaction_axes.pdf", width = 8.5, height = 7.5, deck_h = 7.0)

# =============================================================================
# FIG 3b -- Top up/down TF per contrast (THE ONLY multi-panel)
# Justification: all 4 facets share the SAME x-axis (ULM score) and are directly
# comparable lollipops -> a same-axis facet_wrap, NOT a heterogeneous composite.
# CLAIM: illustrates what each contrast (WT_heat/KO_heat/Temp_main/Interaction)
#        captures.
# =============================================================================
fig3b_df <- rd("fig3b_top_tf_by_contrast_data.csv")
fig3b_df$contrast <- factor(fig3b_df$contrast,
                            levels = c("WT_heat", "KO_heat", "Temp_main", "Interaction"))
fig3b_df$axis     <- factor(fig3b_df$axis, levels = c("HIF", "IFN", "other"))
fig3b_df$key_tf   <- as.logical(fig3b_df$key_tf)
# Within-facet ordering by score; tidytext-free ordering via interaction key.
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
  facet_wrap(~ contrast, scales = "free_y", ncol = 2) +
  scale_color_manual(values = AXIS_COLORS, name = "TF axis", labels = AXIS_LABELS) +
  scale_size_manual(values = c("TRUE" = 2.8, "FALSE" = 1.5), guide = "none") +
  labs(
    title = "Fig 3b. Top activated / suppressed TFs per contrast (CollecTRI ULM)",
    subtitle = "Shared x-axis (ULM score) across all four facets. Key TFs (HIF orange, IFN/NFkB blue) labelled.\nReads out what each contrast captures: heat arms vs the temperature main effect vs the cGAS x heat Interaction.",
    x = "TF activity (ULM score)", y = NULL,
    caption = cap
  ) +
  base_theme +
  theme(axis.text.y = element_text(size = 6.5),
        strip.text  = element_text(face = "bold", size = 9))

save_fig(fig3b, "fig3b_top_tf_by_contrast.pdf", width = 11, height = 9, deck_h = 8.0)

# =============================================================================
# FIG 3c -- Hif1a rank cascade across 4 method/network configs
# CLAIM: the #1 DoRothEA ranking is method/network-fragile (#1->#12->#142->#8).
# =============================================================================
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
    title = "Fig 3c. Hif1a rank cascade across inference method / network (WT_heat)",
    subtitle = "rank #1 (top) = most activated. The DoRothEA #1 collapses under the CollecTRI swap (#12), the\nmultivariate MLM de-confounding (#142), and is partially re-inflated by the consensus (#8): method/network-fragile.",
    x = NULL, y = "Hif1a rank among activated TFs (lower = higher activity)",
    caption = cap
  ) +
  base_theme +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))

save_fig(fig3c, "fig3c_hif1a_rank_cascade.pdf", width = 8.5, height = 6.5, deck_h = 6.0)

# =============================================================================
# FIG 3d -- Regulon swap: explains #1 -> #12 (DoRothEA -> CollecTRI)
# DESIGN: single claim-panel = per-target t_wt distribution by membership
#         (boxplot + jittered points). The regulon-size / overlap counts (131 vs
#         353, 101 shared) are an annotation, NOT a co-equal second axis.
# CLAIM: CollecTRI dilutes the Hif1a regulon with many low-|t| collectri_only
#        targets, dragging the ULM score down.
# =============================================================================
fig3d_df  <- rd("fig3d_regulon_swap_data.csv")
fig3d_sum <- rd("fig3d_regulon_swap_summary.csv")
fig3d_df$membership <- factor(fig3d_df$membership,
                              levels = c("dorothea_only", "both", "collectri_only"),
                              labels = c("DoRothEA-only", "shared (both)", "CollecTRI-only"))

mem_cols <- c("DoRothEA-only" = "#1B7837", "shared (both)" = "grey45", "CollecTRI-only" = "#762A83")

dor <- fig3d_sum[fig3d_sum$network == "DoRothEA", ]
ct  <- fig3d_sum[fig3d_sum$network == "CollecTRI", ]
swap_note <- sprintf(
  "DoRothEA regulon n=%d (mean|t|=%.2f) -> CollecTRI regulon n=%d (mean|t|=%.2f)\n%d shared; %d DoRothEA-only; %d CollecTRI-only (high |t|, mixed sign) targets added",
  dor$n, dor$mean_abs_t, ct$n, ct$mean_abs_t,
  dor$n_intersection[1], dor$n_dorothea_only[1], dor$n_collectri_only[1])

fig3d <- ggplot(fig3d_df, aes(x = membership, y = t_wt, fill = membership)) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55) +
  geom_jitter(width = 0.16, height = 0, size = 1, alpha = 0.5, color = "grey20") +
  scale_fill_manual(values = mem_cols, guide = "none") +
  annotate("label", x = 0.6, y = max(fig3d_df$t_wt) * 0.96, hjust = 0, vjust = 1,
           label = swap_note, size = 2.7, fill = "grey96", color = "grey15") +
  labs(
    title = "Fig 3d. Regulon swap dilutes Hif1a's ULM signal (explains #1 -> #12)",
    subtitle = "Per-target WT_heat t-statistic by membership. The CollecTRI-only targets the swap ADDS carry high |t| but\nMIXED sign (their box straddles 0), so they cancel and add little net directional signal -> diluting the\nmor-weighted regulon mean and dragging the Hif1a ULM score below the tighter DoRothEA value.",
    x = "Hif1a target membership across networks", y = "WT_heat target t-statistic",
    caption = cap
  ) +
  base_theme

save_fig(fig3d, "fig3d_regulon_swap.pdf", width = 8.5, height = 6.5, deck_h = 6.0)

# =============================================================================
# FIG 3e -- MLM collinearity: explains #12 -> #142 (ULM -> MLM)
# DESIGN: single claim-panel = histogram of n_other_TFs across Hif1a targets.
#         The ULM->MLM score collapse (5.11 -> 1.13) is an annotation, NOT a
#         second mixed axis.
# CLAIM: Hif1a's targets are massively co-regulated, so the multivariate MLM
#        attributes their signal away from Hif1a, collapsing its score.
# =============================================================================
fig3e_df  <- rd("fig3e_mlm_collinearity_data.csv")
fig3e_sc  <- rd("fig3e_score_collapse.csv")
fig3e_sum <- rd("fig3e_mlm_collinearity_summary.csv")   # pre-computed aggregates (no stats in viz)
fig3e_df$co_regulated <- as.logical(fig3e_df$co_regulated)

pct_co  <- 100 * fig3e_sum$frac_co_regulated
mean_ot <- fig3e_sum$mean_other_TFs
ulm <- fig3e_sc[fig3e_sc$method == "ULM", ]
mlm <- fig3e_sc[fig3e_sc$method == "MLM", ]
collapse_note <- sprintf(
  "Hif1a ULM score %.2f (rank #%d)\n  -> MLM score %.2f (rank #%d)\n%.0f%% of targets co-regulated; mean %.0f other TFs/target",
  ulm$score, ulm$rank, mlm$score, mlm$rank, pct_co, mean_ot)

fig3e <- ggplot(fig3e_df, aes(x = n_other_TFs)) +
  geom_histogram(binwidth = 10, fill = "#2166AC", color = "white", boundary = 0) +
  geom_vline(xintercept = mean_ot, color = UP, linewidth = 0.8, linetype = "22") +
  annotate("label", x = max(fig3e_df$n_other_TFs) * 0.98, y = Inf, hjust = 1, vjust = 1.1,
           label = collapse_note, size = 2.9, fill = "grey96", color = "grey15") +
  labs(
    title = "Fig 3e. Target collinearity collapses Hif1a under MLM (explains #12 -> #142)",
    subtitle = "For each Hif1a target: how many OTHER CollecTRI TFs also regulate it. The dashed line marks the mean.\nBecause nearly every Hif1a target is co-regulated, the multivariate MLM re-attributes its signal away from Hif1a.",
    x = "number of OTHER CollecTRI TFs sharing each Hif1a target", y = "Hif1a targets (count)",
    caption = cap
  ) +
  base_theme

save_fig(fig3e, "fig3e_mlm_collinearity.pdf", width = 8.5, height = 6, deck_h = 5.625)

# =============================================================================
# FIG 3f -- Consensus z-mix: explains #142 -> #8 (MLM -> consensus)
# CLAIM: only mlm is low; the univariate-family statistics outvote it and the
#        consensus mean (marked) re-inflates Hif1a back to #8.
# =============================================================================
fig3f_df <- rd("fig3f_consensus_zmix_data.csv")
stat_order <- c("ulm", "mlm", "wsum", "norm_wsum", "corr_wsum", "consensus")
fig3f_df$statistic <- factor(fig3f_df$statistic, levels = stat_order)
fig3f_df$family    <- factor(fig3f_df$family,
                             levels = c("univariate", "multivariate", "consensus"))
fig3f_df <- fig3f_df[order(fig3f_df$statistic), ]

consensus_z <- fig3f_df$z_score[fig3f_df$statistic == "consensus"]
fam_cols <- c("univariate" = UP, "multivariate" = "#2166AC", "consensus" = "grey25")

fig3f <- ggplot(fig3f_df, aes(x = statistic, y = z_score, fill = family)) +
  geom_hline(yintercept = consensus_z, color = "grey40", linetype = "22", linewidth = 0.6) +
  geom_col(width = 0.66) +
  geom_text(aes(label = sprintf("%.2f", z_score)), vjust = -0.4, size = 3, color = "grey20") +
  annotate("text", x = 0.6, y = consensus_z, hjust = 0, vjust = -0.4,
           label = sprintf("consensus mean z = %.2f", consensus_z), size = 2.9, color = "grey35") +
  scale_fill_manual(values = fam_cols, name = "statistic family") +
  labs(
    title = "Fig 3f. The consensus outvotes MLM, re-inflating Hif1a (explains #142 -> #8)",
    subtitle = "Hif1a's per-statistic folded-z. Only mlm (multivariate) is low; the four univariate statistics stay high\nand the consensus is their mean (dashed line) -> Hif1a is pulled back up to #8 despite the MLM verdict.",
    x = "decoupleR statistic", y = "Hif1a folded z-score",
    caption = cap
  ) +
  base_theme

save_fig(fig3f, "fig3f_consensus_zmix.pdf", width = 8.5, height = 6, deck_h = 5.625)

# =============================================================================
# FIG 3g -- EXPANDED: full ranked contribution landscape of all 353 Hif1a
#           regulon members (WT_heat).
# OBJECT: "a heat-induced glycolytic/stress program partially overlapping HIF targets"
#         (NEVER "the HIF program").
# CLAIM: Hif1a's positive WT_heat ULM score is built almost entirely from its
#        340 'other' regulon members -- generic stress/ECM/proliferation genes
#        tower over the HIF-specific core. The right tail of 'other' is heat-
#        shock/stress (Timp1/Sdc1/Cdkn1a/Serpine1/Eno2/Hspa1a), visually
#        connecting to fig3l's module bucketing.
# VIZ ONLY -- no statistics. All data from fig3g_target_decomposition_data.csv
#             (unchanged) and fig3g_target_decomposition_summary.csv.
# =============================================================================
fig3g_df  <- rd("fig3g_target_decomposition_data.csv")
fig3g_sum <- rd("fig3g_target_decomposition_summary.csv")

# --- stress-contaminant highlight genes (spec §2 / fig3m note) ---------------
STRESS_GENES <- c("Timp1", "Sdc1", "Cdkn1a", "Serpine1", "Eno2", "Hspa1a", "Spp1")

# --- genes to text-label (spec §2 fig3g): top stress drivers + HIF-specific ---
LABEL_GENES <- c(
  # stress drivers: the 6 named in the spec
  "Timp1", "Sdc1", "Cdkn1a", "Serpine1", "Eno2", "Hspa1a",
  # HIF-specific members (all 7, since there are only 7)
  "Vegfa", "Egln3", "Slc2a1", "Car9", "Bnip3", "Bnip3l", "Pdk1"
)

# --- rank rows by contrib (ascending = left-to-right waterfall order) ---------
fig3g_df <- fig3g_df[order(fig3g_df$contrib), ]
fig3g_df$rank <- seq_len(nrow(fig3g_df))

# --- assign display color:  HIF-specific / shared-glycolytic keep their class
#     colors; stress contaminants override to MODULE_COLORS["heatshock_stress"];
#     remaining 'other' = grey80. ------------------------------------------
# RECONCILED VOCAB (item 4): classes come from the compute table and now match
# fig3l's finer modules -- "hypoxic HIF core" (REPRESSED diagnostic core) is
# colored with the distinct repressed-core hue; "shared/glycolytic" (UP HIF
# targets + glycolytic) is the light-orange up-set; heat-shock = purple; other = grey.
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

# Build a named color vector for scale_color_manual ---------------------------
# Four display categories used on the plot.
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

# --- viz-local convenience table of labelled genes (re-rank, NOT a statistic) -
fig3g_labels_df <- fig3g_df[fig3g_df$target %in% LABEL_GENES,
                             c("source","target","mor","t_wt","contrib","class","rank","legend_class")]
write.csv(fig3g_labels_df,
          file.path(TBL_DIR, "fig3g_target_landscape_labels.csv"),
          row.names = FALSE)
cat("  wrote fig3g_target_landscape_labels.csv\n")

# --- text-label data: label position nudged above stem-head ------------------
fig3g_lab_sub <- fig3g_df[fig3g_df$target %in% LABEL_GENES, ]
# nudge direction: right of zero point upward, left of zero point downward
fig3g_lab_sub$nudge_x <- ifelse(fig3g_lab_sub$contrib >= 0,  0.6, -0.6)
fig3g_lab_sub$hjust   <- ifelse(fig3g_lab_sub$contrib >= 0, 0, 1)

# --- annotation text from summary table --------------------------------------
get_pct_g  <- function(cl) fig3g_sum$pct_of_total[fig3g_sum$class == cl]
get_sum_g  <- function(cl) fig3g_sum$sum_contrib[fig3g_sum$class == cl]
annot_text <- sprintf(
  "Score provenance (%%  of aggregate signal):\n  other (unclassified):    %.1f%%\n  shared/glycolytic:        %.1f%%\n  hypoxic HIF core:         %.1f%% (REPRESSED)",
  get_pct_g("other"), get_pct_g("shared/glycolytic"), get_pct_g("hypoxic HIF core"))

# Coarse 3-class view of the SAME regulon fig3l decomposes finer (fig3l carves the
# 7 heat-shock genes out of this "other" lump) -- the %s differ by construction, not error.
sub_g_expanded <- paste0(
  "All 353 Hif1a regulon members, ranked by signed contribution (sign(mor) x t_wt). Stress/ECM genes (purple) tower over the\n",
  sprintf("regulon; the HIF-diagnostic hypoxic core (Pdk1/Bnip3/Bnip3l/Car9, teal) is REPRESSED (%+.2f total).\n",
          get_sum_g("hypoxic HIF core")),
  sprintf("Coarse 3-class provenance: %.1f%% other | %.1f%% shared/glycolytic | %.1f%% hypoxic HIF core (fig3l decomposes 'other' finer).",
          get_pct_g("other"), get_pct_g("shared/glycolytic"), get_pct_g("hypoxic HIF core"))
)

fig3g <- ggplot(fig3g_df, aes(x = rank, y = contrib, color = legend_class)) +
  # zero-line
  geom_hline(yintercept = 0, color = "grey55", linewidth = 0.35) +
  # vertical stems (lollipop tails)
  geom_segment(aes(xend = rank, yend = 0), linewidth = 0.35, alpha = 0.7) +
  # lollipop heads
  geom_point(aes(size = legend_class %in% c("hypoxic HIF core (repressed)", "heat-shock/stress", "shared/glycolytic")),
             alpha = 0.85) +
  # gene labels (message-carrying only) -- ggrepel for collision-free placement
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
  # annotation box (provenance breakdown, from summary table)
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
    title    = "Fig 3g. WHERE Hif1a's score comes from: full ranked regulon contribution landscape (WT_heat)",
    subtitle = sub_g_expanded,
    x        = "regulon members ranked by contribution (left = most negative; right = most positive)",
    y        = "contribution to WT_heat ULM score (sign(mor) x t_wt)",
    caption  = cap
  ) +
  base_theme +
  theme(
    legend.position  = "right",
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 8, face = "bold"),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    panel.grid.major.x = element_blank()
  )

save_fig(fig3g, "fig3g_target_decomposition.pdf", width = 11, height = 7, deck_h = 6.5)

# =============================================================================
# FIG 3h -- Lombardi vs Phylo HIF signatures x 3 contrasts
# CLAIM: HIF arm is cGAS-independent by a curation-independent measure -- both
#        signatures up in BOTH heat arms, Interaction NS for both.
# (Title/legend text taken from the DYNAMIC sig_label column; never hardcoded.)
# =============================================================================
fig3h_df <- rd("fig3h_lombardi_vs_phylo_data.csv")
fig3h_df$condition <- factor(fig3h_df$condition,
                             levels = c("WT_heat", "KO_heat", "Interaction"))
fig3h_df$sig <- as.logical(fig3h_df$sig)
# Build legend labels from the dynamic sig_label per signature.
sig_lab_map <- tapply(fig3h_df$sig_label, fig3h_df$signature, function(x) x[1])

fig3h <- ggplot(fig3h_df, aes(x = condition, y = score, fill = signature)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.3) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(aes(label = ifelse(sig, "*", "ns")),
            position = position_dodge(width = 0.75), vjust = -0.35, size = 3.4) +
  scale_fill_manual(values = setNames(c(UP, "#7B3294"), names(sig_lab_map)),
                    labels = sig_lab_map, name = "HIF signature") +
  labs(
    title = "Fig 3h. Provenance / hypoxia control: a field-consensus and a hand-made HIF list both rise with heat",
    subtitle = sprintf("Field-consensus %s and hand-made %s both go UP with heat in BOTH genotypes (Interaction NS for both) -- a hypoxia/provenance\ncheck that the hand-made Phylo list isn't idiosyncratic, NOT thermal-HIF proof. CAVEAT: each metagene COLLAPSES HIF1/HIF2 and is\nDIRECTION-BLIND -- it cannot see the repressed hypoxic core (Pdk1/Bnip3/Bnip3l/Car9; fig3l). * = BH padj < 0.05; ns otherwise.",
                       sig_lab_map[["Lombardi2022_HIF"]], sig_lab_map[["Phylo16_HIF"]]),
    x = NULL, y = "signature activity (ULM score)",
    caption = cap
  ) +
  base_theme

save_fig(fig3h, "fig3h_lombardi_vs_phylo.pdf", width = 8.5, height = 6, deck_h = 5.625)

# =============================================================================
# FIG 3i -- Interaction primer: how the cGAS x heat Interaction is read
# CLAIM: The Interaction = (WT_39-WT_37) - (cGASKO_39-cGASKO_37).
#        Ifit1 (IFN arm, cGAS-dependent): slopes DIVERGE -> large positive Interaction.
#        Vegfa (HIF arm, cGAS-independent): slopes PARALLEL -> Interaction ~0.
# Design: two facets (free_y; same metric log2CPM, different absolute levels),
#         WT solid / cGAS-KO dashed, annotated with contrast values from table.
# ZERO statistics in this block: all numbers come from fig3i_interaction_primer_data.csv.
# =============================================================================
fig3i_df <- rd("fig3i_interaction_primer_data.csv")
fig3i_df$gene     <- factor(fig3i_df$gene, levels = c("Ifit1", "Vegfa"))
fig3i_df$temp     <- factor(fig3i_df$temp, levels = c("37", "39"))
fig3i_df$genotype <- factor(fig3i_df$genotype, levels = c("WT", "cGASKO"))

# LAYOUT FIX (TASK 2): the contrast numbers used to be geom_label'd INTO the
# panel (x=1.5, y=-Inf), overlapping the lines. Move them OUT of the data region
# by folding them into the per-facet STRIP text via a labeller -> they now sit in
# the strip band ABOVE each panel, never over the plotted slopes. (No stats here:
# every number is read straight from the pre-computed contrast columns.)
strip_labels <- vapply(c("Ifit1", "Vegfa"), function(g) {
  r <- fig3i_df[fig3i_df$gene == g, ][1, ]
  arm_tag <- if (g == "Ifit1") "IFN arm" else "shared HIF/glycolytic target"
  sprintf(
    "%s  [%s]\nWT heat %+.2f  /  KO heat %+.2f\nInteraction %+.2f (adjP %.3f)",
    g, arm_tag, r$wt_heat, r$ko_heat, r$interaction, r$inter_adjP
  )
}, character(1))

# Qualitative genotype palette: WT = teal/dark, cGASKO = orange (house UP color).
GENO_COLS <- c("WT" = "#1B7837", "cGASKO" = "#D6604D")

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
    title    = "Fig 3i. Reading the cGAS x heat Interaction: difference of the two genotype heat-slopes",
    subtitle = paste0(
      "Interaction = (WT_39 - WT_37) - (KO_39 - KO_37).  WT solid, cGAS-KO dashed.\n",
      "Diverging slopes = cGAS-dependent (Ifit1, IFN arm).  ",
      "Parallel slopes = no detectable cGAS-dependence (Vegfa, shared HIF/glycolytic target)."
    ),
    x = "Temperature", y = "mean log2 CPM (group mean)",
    caption = cap
  ) +
  base_theme +
  theme(
    plot.title    = element_text(size = 11, face = "bold"),
    plot.subtitle = element_text(size = 8.5, lineheight = 1.25),
    strip.text    = element_text(face = "bold", size = 9, lineheight = 1.1),
    strip.background = element_rect(fill = "grey94", color = NA),
    panel.spacing = unit(1.4, "lines"),
    plot.margin   = margin(t = 6, r = 8, b = 4, l = 8, unit = "mm")
  )

save_fig(fig3i, "fig3i_interaction_primer.pdf", width = 9, height = 6, deck_h = 5.625)

cat("[DONE] 03_decoupler_tf_viz.R -- 9 single-claim Phase-3 figures rendered from tidy tables.\n")
