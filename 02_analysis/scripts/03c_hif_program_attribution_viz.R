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
source("02_analysis/helpers/figure_style.R")       # contract: project_theme + save_overview; FIG_CFG
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

STAGE   <- "04_tf"
SCRIPT  <- "02_analysis/scripts/03c_hif_program_attribution_viz.R"
TBL_DIR <- stage_dir("04_tf", "tables")

cap <- provisional_caption()
rd  <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)


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
# The punchline numbers ride in the subtitle (kept in the deck PNG); the carve-out
# reconciliation note rides in the caption -- NO in-plot text box (it overlaid the
# panel). "net HIF-specific core" = shared + feedback + hypoxic-core.
quant_line <- sprintf(
  "Heat-shock/stress bucket %+.2f  >>  net HIF-specific core %+.2f;  HIF1a-selective hypoxic core %+.2f (REPRESSED).",
  stress_sum, net_hif_core, core_sum)
carve_note <- "Heat-shock genes carved out of fig3g's 91.8% 'other' -> 90.7% residual here."

xr <- max(abs(plot_df$contrib)) * 1.06

fig3l <- ggplot(plot_df, aes(x = contrib, y = target, color = module)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_segment(aes(xend = 0, yend = target), linewidth = 0.7) +
  geom_point(size = 3.2) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(module = module_labels)) +
  scale_color_manual(values = MODULE_COLORS, guide = "none") +
  scale_x_continuous(limits = c(-xr, xr),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "Fig 3l. A heat-induced glycolytic/stress program partially overlapping HIF targets",
    subtitle = paste0(
      "Hif1a's positive WT_heat 'activity', re-bucketed by curated biological module (signed contribution = sign(mor) x t_wt,\n",
      "copied from fig3g -- no recomputation). Heat-shock/stress UP, shared/glycolytic UP, feedback UP, but the HIF1a-selective\n",
      "hypoxic-survival core (Pdk1/Bnip3/Bnip3l/Car9) is REPRESSED (left of zero). NOT a canonical hypoxic-HIF output.\n",
      quant_line),
    x = "contribution to WT_heat ULM signal (signed t * mor)", y = NULL
  ) +
  project_theme(config = FIG_CFG) +
  theme(strip.text.y = element_text(angle = 0, face = "bold", size = 8, lineheight = 0.95),
        panel.spacing = unit(0.6, "lines"))

# One-time relocation cleanup: this script previously wrote flat-dir PDFs into 04_tf/figures/.
.FLAT_FIG_DIR <- stage_dir("04_tf", "figures")
.stale_pdf <- file.path(.FLAT_FIG_DIR, "fig3l_hif_attribution.pdf")
if (file.exists(.stale_pdf)) file.remove(.stale_pdf)

save_overview(
  fig3l, STAGE, "fig3l_hif_attribution",
  table = plot_df[, c("target", "module", "contrib", "direction")],
  finding = paste0(
    "Re-bucketed by curated biological module, Hif1a's positive WT_heat activity is ",
    "carried by a heat-shock/stress fraction (all UP) while the HIF1a-selective ",
    "hypoxic-survival core is repressed; only the shared/glycolytic/feedback subset ",
    "rises. This is a heat-induced glycolytic/stress program partially overlapping ",
    "HIF targets -- NOT a canonical hypoxic-HIF output. ", cap),
  script = SCRIPT, fn = "save_overview",
  config_kv = "figures.base_size=16; figures.base_size_column=9; modules",
  input = "03_results/04_tf/tables/fig3l_hif_attribution_data.csv",
  how_to_read = paste0(
    "Horizontal lollipop plot: x = signed contribution = sign(mor) x t_wt of Hif1a ",
    "targets, grouped and colored by module. Facets are ordered from stress/contamination ",
    "(top) to hypoxic-survival core (bottom; diagnostic of true hypoxia, all repressed). ",
    "Claim tier: descriptive target attribution. ",
    "Note: ", carve_note),
  width = 9, height = 7.5, config = FIG_CFG)

cat("[DONE] 03c_hif_program_attribution_viz.R -- fig3l rendered (no statistics).\n")
