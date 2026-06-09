#!/usr/bin/env Rscript
# =============================================================================
# 03c_hif_program_attribution_viz.R - PHASE 3 (stage 04_tf) VIZ: fig3l
#   the attribution CENTREPIECE -- module-bucketed, sign-aware lollipop
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# Role:    VISUALIZE half of the "normalize-then-visualize" split. Reads ONLY
#          the tidy tables emitted by 03c_hif_program_attribution.R and renders
#          ONE single-claim figure. Performs NO statistics (no run_ulm/run_mlm,
#          no p.adjust/prcomp/cor); contrib + module + direction are read as-is.
#          Runs STANDALONE after the compute script.
#
# Inputs (03_results/04_tf/tables/):
#   - fig3l_hif_attribution_data.csv   (per-member module/contrib/direction)
#   - fig3l_module_summary.csv         (per-module aggregates; for the callout)
#
# Outputs:
#   - 03_results/04_tf/figures/fig3l_hif_attribution.pdf  (canonical, STAMPED)
#   - 03_results/04_tf/deck_assets/fig3l_hif_attribution.png  (DECK_EXPORT=1;
#       clean 300-dpi, caption-stripped only)
#
# Object name (HARD GUARDRAIL, design-spec sec 0): "a heat-induced glycolytic/
#   stress program partially overlapping HIF targets" -- NEVER "the HIF
#   program" / "HIF1a the TF". HIF2a/Epas1 never crowned.
#
# Encoding: horizontal sign-aware lollipop, x = signed contrib, faceted by
#   module, facets ordered heatshock_stress -> shared_angio_glucose ->
#   autoreg_feedback -> hif1a_hypoxic_core (punchline at bottom: hypoxic core
#   left of zero). MODULE_COLORS from config.R (single source of truth).
# Dependencies: config.R; ggplot2, dplyr
# =============================================================================

source("02_analysis/config/config.R")
load_packages()   # ggplot2 + dplyr

TBL_DIR <- stage_dir("04_tf", "tables")
FIG_DIR <- stage_dir("04_tf", "figures")

cap <- provisional_caption()
rd  <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# DECK EXPORT: when DECK_EXPORT is set, ALSO write a clean 300-dpi PNG (no
# embargo caption) to 03_results/04_tf/deck_assets/. The stamped PDF is always
# written, unchanged; only labs(caption=NULL) is dropped for the PNG.
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
           width = 10, height = deck_h, dpi = 300)
    cat(sprintf("  wrote deck_assets/%s (clean, no embargo stamp)\n", png_name))
  }
}

base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey25"),
        plot.caption  = element_text(size = 7, color = "grey45", hjust = 0, face = "italic"))

# =============================================================================
# FIG 3l -- module-bucketed, sign-aware attribution of Hif1a's WT_heat "activity"
# CLAIM: re-bucketed by biological module, the positive Hif1a "activity" is
#        carried by a heat-shock/stress fraction (all UP) while the HIF1a-
#        selective hypoxic-survival core (Pdk1/Bnip3/Bnip3l/Car9) is REPRESSED
#        (left of zero); only the shared/glycolytic/feedback subset rises. This
#        is a heat-induced glycolytic/stress program partially overlapping HIF
#        targets -- NOT a canonical hypoxic-HIF output.
# =============================================================================
dat  <- rd("fig3l_hif_attribution_data.csv")
summ <- rd("fig3l_module_summary.csv")

# Facet order = the reading order (stress UP ... shared UP ... feedback UP ...
# hypoxic core DOWN). other_unclassified is dropped to a single summarized row
# so the 4 curated buckets stay legible (it is the faint context bucket).
curated <- c("heatshock_stress", "shared_angio_glucose",
             "autoreg_feedback", "hif1a_hypoxic_core")
module_labels <- c(
  heatshock_stress     = "heat-shock / stress\n(non-HIF; the contamination)",
  shared_angio_glucose = "shared angio / glucose\n(shared HIF1/HIF2)",
  autoreg_feedback     = "autoreg feedback\n(Egln3/PHD3 brake)",
  hif1a_hypoxic_core   = "HIF1a-selective\nhypoxic core (DIAGNOSTIC)"
)

plot_df <- dat[dat$module %in% curated, ]
plot_df$module <- factor(plot_df$module, levels = curated)
# Within facet: order by contrib so the lollipops read cleanly.
plot_df <- plot_df[order(plot_df$module, plot_df$contrib), ]
plot_df$target <- factor(plot_df$target, levels = unique(plot_df$target))

# Callout numbers (read from the summary table -- no recomputation).
stress_sum <- summ$sum_contrib[summ$module == "heatshock_stress"]
core_sum   <- summ$sum_contrib[summ$module == "hif1a_hypoxic_core"]
# "net HIF-specific core" = shared + feedback + hypoxic-core (the +14.27 the
# synthesis quotes: Vegfa/Egln3/Slc2a1 up minus Pdk1/Bnip3l/Bnip3/Car9 down).
net_hif_core <- sum(summ$sum_contrib[summ$module %in%
                      c("shared_angio_glucose", "autoreg_feedback", "hif1a_hypoxic_core")])

# N1 reconcile: fig3l carves the 7 heat-shock genes OUT of fig3g's 340-member
# "other" lump (91.8%) into a named bucket, leaving a 339-member residual "other"
# (90.7%) -- the two adjacent figures differ by construction, not by error.
callout <- sprintf(
  "heat-shock/stress bucket sum = %+.2f\n>> net HIF-specific core (%+.2f)\nHIF1a-selective hypoxic core = %+.2f (REPRESSED)\n(heat-shock genes carved OUT of fig3g's 91.8%% 'other' -> 90.7%% residual here)",
  stress_sum, net_hif_core, core_sum)

xr <- max(abs(plot_df$contrib)) * 1.18

# Callout placed ONCE, in the top (heat-shock) facet only, in the empty left
# (negative) region so it never overlaps the right-leaning stress lollipops.
callout_df <- data.frame(
  module = factor("heatshock_stress", levels = curated),
  x = -xr, y = Inf, label = callout, stringsAsFactors = FALSE)

fig3l <- ggplot(plot_df, aes(x = contrib, y = target, color = module)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_segment(aes(xend = 0, yend = target), linewidth = 0.7) +
  geom_point(size = 3.2) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(module = module_labels)) +
  scale_color_manual(values = MODULE_COLORS, guide = "none") +
  scale_x_continuous(limits = c(-xr, xr),
                     expand = expansion(mult = c(0.02, 0.02))) +
  geom_label(data = callout_df, inherit.aes = FALSE,
             aes(x = x, y = y, label = label), hjust = 0, vjust = 1.1,
             size = 2.8, fill = "grey96", color = "grey15") +
  labs(
    title = "Fig 3l. A heat-induced glycolytic/stress program partially overlapping HIF targets",
    subtitle = paste0(
      "Hif1a's positive WT_heat 'activity', re-bucketed by curated biological module (signed contribution = sign(mor) x t_wt,\n",
      "copied from fig3g -- no recomputation). Heat-shock/stress UP, shared/glycolytic UP, feedback UP, but the HIF1a-selective\n",
      "hypoxic-survival core (Pdk1/Bnip3/Bnip3l/Car9) is REPRESSED (left of zero). NOT a canonical hypoxic-HIF output."),
    x = "contribution to WT_heat ULM signal (signed t * mor)", y = NULL,
    caption = cap
  ) +
  base_theme +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8, lineheight = 0.95),
        panel.spacing = unit(0.6, "lines"))

save_fig(fig3l, "fig3l_hif_attribution.pdf", width = 9, height = 7.5, deck_h = 6.5)

cat("[DONE] 03c_hif_program_attribution_viz.R -- fig3l rendered from tidy tables (no statistics).\n")
