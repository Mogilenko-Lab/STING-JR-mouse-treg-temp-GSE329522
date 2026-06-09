#!/usr/bin/env Rscript
# =============================================================================
# 03d_ulm_mechanic_viz.R - PHASE 3 (stage 04_tf) VIZ: fig3m
#   "How decoupleR-ULM nominates a TF" -- the SCORING MECHANIC, pedagogical.
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# Role:    VISUALIZE half of the "normalize-then-visualize" split. Reads ONLY the
#          plot-ready tidy tables emitted by 03d_ulm_mechanic.R and renders the
#          figure. Performs NO statistics (no run_ulm/run_mlm/p.adjust/prcomp/
#          cor); it only plots already-computed columns with cosmetic reshaping.
#          Runs STANDALONE after the compute script.
#
# Design (design-spec sec 2 fig3m): TWO panels SHARING the x-axis = signed/
#   aligned target contribution (aligned_contrib). The shared axis is what
#   licenses the side-by-side (one-axis rule).
#     Panel A = Hif1a (PROMISCUOUS) -- broad regulon; generic stress genes push
#               the pile RIGHT. Contaminants highlighted in MODULE_COLORS
#               ["heatshock_stress"] + labelled; HIF-specific members shown in
#               MODULE_COLORS["hif1a_hypoxic_core"] for continuity with fig3g/3l.
#     Panel B = Stat2 (SPECIFIC) -- small IFN-restricted regulon; without
#               coherent movement it does NOT pile up (the foil).
#   Each TF's aggregate ULM score is drawn as a labelled diamond + vertical rule
#   at the regulon's WEIGHTED CENTER (mean aligned contribution) -- visually
#   "the score is the pile-up". One literal annotation states the ULM model:
#   t_gene ~ mor_TF -> score.
#
# Inputs (03_results/04_tf/tables/, written by 03d_ulm_mechanic.R):
#   fig3m_ulm_mechanic_data.csv      fig3m_ulm_mechanic_summary.csv
#
# Outputs:
#   03_results/04_tf/figures/fig3m_ulm_mechanic.pdf         (canonical, stamped)
#   03_results/04_tf/deck_assets/fig3m_ulm_mechanic.png     (DECK_EXPORT=1; clean)
#
# Object name (HARD GUARDRAIL, sec 0): the Hif1a panel describes "a heat-induced
#   glycolytic/stress program partially overlapping HIF targets" -- NEVER "the
#   HIF program" / "HIF1a the TF". The comparator NEVER implies HIF2a anything.
#   This figure is about the SCORING MECHANIC, not a biology claim.
# Dependencies: config.R; ggplot2, dplyr (no statistics)
# =============================================================================

source("02_analysis/config/config.R")
load_packages()  # ggplot2 + dplyr; limma/edgeR pulled in but unused

TBL_DIR <- stage_dir("04_tf", "tables")
FIG_DIR <- stage_dir("04_tf", "figures")

cap <- provisional_caption()
rd  <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)

# Cross-figure palette continuity (single source of truth: config.R MODULE_COLORS).
COL_STRESS <- unname(MODULE_COLORS["heatshock_stress"])    # purple -- contamination
COL_HIF    <- unname(MODULE_COLORS["hif1a_hypoxic_core"])  # full HIF-orange -- specific core
COL_OTHER  <- "grey70"                                     # faint context members
COL_SCORE  <- "grey15"                                     # the aggregate-score diamond/rule

# -----------------------------------------------------------------------------
# DECK EXPORT: when DECK_EXPORT is set, ALSO write a clean 300-dpi PNG (no
# embargo caption line) to 03_results/04_tf/deck_assets/. Canonical stamped PDF
# is ALWAYS written, unchanged. Only labs(caption=NULL) is dropped for the PNG.
# -----------------------------------------------------------------------------
DECK_EXPORT <- nzchar(Sys.getenv("DECK_EXPORT"))
DECK_DIR    <- file.path(DIR_RESULTS, "04_tf", "deck_assets")
if (DECK_EXPORT) dir.create(DECK_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(plot, fname, width, height, deck_h = 5.625) {
  ggsave(file.path(FIG_DIR, fname), plot, width = width, height = height)
  cat(sprintf("  wrote %s\n", fname))
  if (DECK_EXPORT) {
    png_name <- sub("\\.pdf$", ".png", fname)
    ggsave(file.path(DECK_DIR, png_name), plot + labs(caption = NULL),
           width = width, height = deck_h, dpi = 300)
    cat(sprintf("  wrote deck_assets/%s (clean, no embargo stamp)\n", png_name))
  }
}

base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey25"),
        plot.caption  = element_text(size = 7, color = "grey45", hjust = 0, face = "italic"),
        strip.text    = element_text(face = "bold", size = 10),
        legend.position = "bottom")

# =============================================================================
# READ tidy tables
# =============================================================================
dat  <- rd("fig3m_ulm_mechanic_data.csv")
summ <- rd("fig3m_ulm_mechanic_summary.csv")
dat$is_stress_contaminant <- as.logical(dat$is_stress_contaminant)

# Panel labels carry the regulon class + member count (object-name disciplined).
n_by_tf  <- setNames(summ$n_targets, summ$tf)
panel_lab <- function(tf) {
  cls <- summ$regulon_class[summ$tf == tf]
  sprintf("%s regulon (%s, n=%d members)", tf, cls, n_by_tf[[tf]])
}
TF_A <- summ$tf[summ$regulon_class == "promiscuous"][1]   # Hif1a
TF_B <- summ$tf[summ$regulon_class == "specific"][1]      # Stat2
dat$panel <- factor(ifelse(dat$tf == TF_A, panel_lab(TF_A), panel_lab(TF_B)),
                    levels = c(panel_lab(TF_A), panel_lab(TF_B)))

# Per-member point class drives color (shared legend across panels):
#   stress contaminant (Hif1a) -> purple;  HIF-specific core (Hif1a) -> orange;
#   everything else -> faint grey.
dat$member_class <- "regulon member"
dat$member_class[dat$is_stress_contaminant] <- "heat-shock/stress contaminant"
# HIF-specific members: the curated hypoxic core + shared/feedback subset that
# fig3g/fig3l label as the HIF-diagnostic fraction (orange for cross-fig continuity).
HIF_SPECIFIC <- c("Vegfa", "Egln3", "Slc2a1", "Car9", "Bnip3", "Bnip3l", "Pdk1")
dat$member_class[dat$tf == TF_A & dat$target %in% HIF_SPECIFIC] <- "HIF-specific member"
dat$member_class <- factor(dat$member_class,
  levels = c("heat-shock/stress contaminant", "HIF-specific member", "regulon member"))

member_cols <- c("heat-shock/stress contaminant" = COL_STRESS,
                 "HIF-specific member"            = COL_HIF,
                 "regulon member"                 = COL_OTHER)

# Aggregate-score diamonds (the weighted center = mean aligned contribution).
score_df <- data.frame(
  panel = factor(c(panel_lab(TF_A), panel_lab(TF_B)),
                 levels = c(panel_lab(TF_A), panel_lab(TF_B))),
  x     = summ$mean_aligned_contrib[match(c(TF_A, TF_B), summ$tf)],
  tf    = c(TF_A, TF_B),
  sign  = summ$ulm_score_sign[match(c(TF_A, TF_B), summ$tf)]
)
score_df$score_label <- sprintf("ULM score (pile-up)\n%s = %+.2f", score_df$tf, score_df$x)
# DECLUTTER (item 6): the score label used to sit on the point cloud and collide
# with the formula + dashed rule at panel center. Park each panel's label in the
# emptier top corner (Hif1a's pile is far-right, so its label goes top-LEFT; Stat2
# is near 0 so its label goes top-RIGHT), connected to the diamond by a thin leader.
xmin <- min(dat$aligned_contrib); xmax <- max(dat$aligned_contrib)
score_df$label_x <- ifelse(score_df$tf == TF_A, xmin * 0.98, xmax * 0.98)
score_df$label_hjust <- ifelse(score_df$tf == TF_A, 0, 1)
score_df$label_y <- 1.82

# Genes to text-label in Panel A: the 7 stress contaminants pushing the pile right.
label_df <- dat[dat$tf == TF_A & dat$is_stress_contaminant, ]

# Reproducible jitter for the beeswarm-style strip.
set.seed(42)

# Formula annotation moved OFF the swarm to a bottom FOOTER band (low y, left edge)
# so it never collides with the center diamond/rule/label. Panel A only.
formula_df <- data.frame(
  panel = factor(panel_lab(TF_A), levels = levels(dat$panel)),
  x = xmin, y = 0.42,
  lab = "ULM model:  t_gene ~ mor_TF -> score   (score = regulon-weighted mean of aligned contributions)"
)

# In-panel verdict cue for the Stat2 (specific) foil panel: a dispersed cloud with
# no coherent pile-up reads as a LOW score. Pinned to the top band, away from points.
stat2_cue_df <- data.frame(
  panel = factor(panel_lab(TF_B), levels = levels(dat$panel)),
  x = xmin, y = 1.82,
  lab = "no pile-up of aligned targets -> low score\n(targets straddle 0; no coherent heat-MAIN signal)"
)

# =============================================================================
# FIG 3m -- the scoring mechanic (two panels, SHARED x = aligned contribution)
# =============================================================================
fig3m <- ggplot(dat, aes(x = aligned_contrib, y = 1)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  # regulon members as a jittered strip/beeswarm on the SHARED contribution axis
  geom_jitter(aes(color = member_class, size = member_class, alpha = member_class),
              height = 0.32, width = 0) +
  # the aggregate ULM score = the pile-up: a vertical rule + labelled diamond.
  # The diamond/rule stay at the weighted center; the LABEL is parked at the right
  # margin (item-6 declutter) with a thin leader so it no longer sits on the swarm.
  geom_vline(data = score_df, aes(xintercept = x),
             color = COL_SCORE, linetype = "22", linewidth = 0.5) +
  geom_point(data = score_df, aes(x = x, y = 1),
             shape = 23, size = 5, fill = COL_SCORE, color = "white", stroke = 0.6) +
  geom_segment(data = score_df, aes(x = x, xend = label_x, y = 1.30, yend = label_y - 0.06),
               color = COL_SCORE, linewidth = 0.3, linetype = "12") +
  geom_text(data = score_df, aes(x = label_x, y = label_y, label = score_label, hjust = label_hjust),
            size = 2.9, color = COL_SCORE, lineheight = 0.9, fontface = "bold") +
  # Stat2 (specific) foil: in-panel "no pile-up -> low score" verdict cue.
  geom_text(data = stat2_cue_df, aes(x = x, y = y, label = lab),
            hjust = 0, size = 2.8, color = "grey25", lineheight = 0.95, fontface = "italic") +
  # label the stress contaminants dragging the Hif1a pile to the right
  ggrepel::geom_text_repel(
    data = label_df, aes(x = aligned_contrib, y = 1, label = target),
    color = COL_STRESS, size = 2.7, fontface = "bold",
    direction = "both",
    ylim = c(NA, 0.78),          # push labels below the point strip (away from diamond/formula)
    xlim = c(NA, NA),
    box.padding        = 0.40,
    point.padding      = 0.30,
    min.segment.length = 0,
    segment.size       = 0.30,
    segment.color      = COL_STRESS,
    max.overlaps       = Inf,
    seed               = 42) +
  # the literal ULM formula, stated once (small + plain)
  geom_text(data = formula_df, aes(x = x, y = y, label = lab),
            hjust = 0, size = 2.7, color = "grey25", fontface = "italic") +
  scale_color_manual(values = member_cols, name = NULL, drop = FALSE) +
  scale_size_manual(values = c("heat-shock/stress contaminant" = 2.6,
                               "HIF-specific member" = 2.6,
                               "regulon member" = 1.5), guide = "none") +
  scale_alpha_manual(values = c("heat-shock/stress contaminant" = 0.95,
                                "HIF-specific member" = 0.95,
                                "regulon member" = 0.45), guide = "none") +
  scale_y_continuous(limits = c(0.3, 1.90), breaks = NULL) +
  facet_wrap(~ panel, ncol = 1, scales = "free_y") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Fig 3m. How decoupleR-ULM nominates a TF: the score is the regulon's pile-up",
    subtitle = paste0(
      "A broad (promiscuous) regulon scores high off any high-|t| members it happens to contain -- here generic heat-shock/stress genes\n",
      "(purple) push the ", TF_A, " pile right, scoring a heat-induced glycolytic/stress program that partially overlaps HIF targets.\n",
      "A small, specific regulon (", TF_B, ") only piles up when its OWN targets coherently move -- on this contrast they do not."),
    x = "Aligned (signed) target contribution  =  sign(mor) × t_gene  (shared axis)",
    y = NULL,
    caption = cap
  ) +
  base_theme

save_fig(fig3m, "fig3m_ulm_mechanic.pdf", width = 9.5, height = 7.0, deck_h = 6.2)

cat("[DONE] 03d_ulm_mechanic_viz.R -- fig3m rendered (no statistics).\n")
