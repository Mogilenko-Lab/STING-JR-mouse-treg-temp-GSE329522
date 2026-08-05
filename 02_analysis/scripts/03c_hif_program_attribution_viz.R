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
#   autoreg_feedback -> hif1a_hypoxic_core (the negative bucket reads last).
#   MODULE_COLORS + MODULE_LABELS from config.R (single source of truth for
#   both the bucket colours and the bucket names, so the human counterpart of
#   this panel reuses the same definitions instead of retyping them).
#
# Drawn text (2026-08-05 rebuild): the title NAMES THE PANEL and states no
#   finding; the subtitle is two short descriptive lines; the strip labels name
#   each bucket by its curated membership and carry that bucket's own n and sum.
#   Every claim about what the pattern MEANS lives in the README caption, which
#   save_overview() writes from `finding=`. No internal figure tags, no
#   capitals-for-emphasis, no "->" as prose on the canvas. fits_canvas() below
#   (from figure_style.R) makes a line that would overflow the canvas a hard
#   error rather than a silently clipped PNG -- that clipping is the defect this
#   rebuild fixes.
# Dependencies: config.R; figure_style.R (fits_canvas); ggplot2, dplyr
# =============================================================================

source("02_analysis/config/config.R")
source("02_analysis/helpers/figure_style.R")       # contract: project_theme + save_overview; FIG_CFG
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

STAGE   <- "04_tf"
SCRIPT  <- "02_analysis/scripts/03c_hif_program_attribution_viz.R"
TBL_DIR <- stage_dir("04_tf", "tables")

cap <- sample_mapping_stamp()
rd  <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)


# =============================================================================
# FIG 3l -- Hif1a's WT_heat target contributions, split by curated module
#
# What the panel shows: the SAME signed contributions already decomposed for the
# Hif1a regulon, regrouped into four curated modules and drawn on one signed
# axis. It recomputes nothing.
#
# What the panel is NOT allowed to say: the reading of that pattern -- which
# module dominates, and what a negative HIF1a-selective core implies -- is a
# finding, and findings live in the README caption (`finding=` below), not in
# drawn text. The bucket names on the canvas describe curated MEMBERSHIP only.
# =============================================================================
dat  <- rd("fig3l_hif_attribution_data.csv")
summ <- rd("fig3l_module_summary.csv")

# Facet order = the reading order, ending on the one module whose members sit
# left of zero. other_unclassified is left out so the 4 curated buckets stay
# legible (it is the faint context bucket; its members are in the table).
curated <- c("heatshock_stress", "shared_angio_glucose",
             "autoreg_feedback", "hif1a_hypoxic_core")

plot_df <- dat[dat$module %in% curated, ]
plot_df$module <- factor(plot_df$module, levels = curated)
# Within facet: order by contrib so the lollipops read cleanly.
plot_df <- plot_df[order(plot_df$module, plot_df$contrib), ]
plot_df$target <- factor(plot_df$target, levels = unique(plot_df$target))

# Per-module n and sum, read from the summary table -- no recomputation. Each
# bucket carries its own two numbers in its own strip, next to its own points,
# so no reader has to map a list of sums back onto a facet order.
mod_n   <- setNames(summ$n_targets,   summ$module)
mod_sum <- setNames(summ$sum_contrib, summ$module)
# Bucket NAMES come from config.R::MODULE_LABELS (shared with the human
# counterpart panel); this script only appends the per-bucket n and sum.
module_labels <- vapply(curated, function(m) {
  sprintf("%s\n%d %s, sum %+.2f", MODULE_LABELS[[m]], mod_n[[m]],
          if (mod_n[[m]] == 1L) "gene" else "genes", mod_sum[[m]])
}, character(1))
names(module_labels) <- curated

# The one number that spans facets rather than sitting inside one: the three
# HIF-annotated modules (shared + feedback + hypoxic core) summed together.
net_hif_core <- sum(mod_sum[c("shared_angio_glucose", "autoreg_feedback",
                              "hif1a_hypoxic_core")])
stress_sum <- mod_sum[["heatshock_stress"]]
core_sum   <- mod_sum[["hif1a_hypoxic_core"]]

# Reconciliation note for the caption (not the canvas): this panel carves the 7
# heat-shock genes out of the coarser fig3g "other" lump (91.8%), leaving a
# 339-member residual "other" (90.7%). The two adjacent figures differ by
# construction, not by error.
carve_note <- paste0(
  "The seven heat-shock genes here are carved out of the coarser 91.8% 'other' ",
  "lump in fig3g, leaving a 339-member 90.7% residual: the two figures differ ",
  "by construction, not by error.")

xr <- max(abs(plot_df$contrib)) * 1.06

# --- drawn text -------------------------------------------------------------
# Title NAMES the panel. Subtitle says what a point is and what x means, in two
# short lines. Neither states a finding, carries a figure tag, or shouts.
FIG_W <- 11; FIG_H <- 7.5
fig_title <- "Hif1a target contributions to the WT_heat ULM signal, by module"
fig_subtitle <- paste0(
  "Each point is one Hif1a CollecTRI target; x is its signed contribution, sign(mor) x t_wt.\n",
  sprintf("Each module carries its own n and sum; the three HIF-annotated modules total %+.2f.",
          net_hif_core))
fig_xlab <- "signed contribution to the WT_heat ULM signal, sign(mor) x t_wt"

.f <- FIG_CFG$figures %||% list()
fits_canvas(fig_title,    .f$title_size    %||% 16, "bold",  FIG_W, "title")
fits_canvas(fig_subtitle, .f$subtitle_size %||% 11, "plain", FIG_W, "subtitle")

fig3l <- ggplot(plot_df, aes(x = contrib, y = target, color = module)) +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.4) +
  geom_segment(aes(xend = 0, yend = target), linewidth = 0.7) +
  geom_point(size = 3.2) +
  facet_grid(module ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(module = module_labels)) +
  scale_color_manual(values = MODULE_COLORS, guide = "none") +
  scale_x_continuous(limits = c(-xr, xr),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(title = fig_title, subtitle = fig_subtitle, x = fig_xlab, y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(strip.text.y  = element_text(angle = 0, face = "bold", size = 11,
                                     lineheight = 1.05),
        panel.spacing = unit(0.6, "lines"))

# One-time relocation cleanup: this script previously wrote flat-dir PDFs into 04_tf/figures/.
.FLAT_FIG_DIR <- stage_dir("04_tf", "figures")
.stale_pdf <- file.path(.FLAT_FIG_DIR, "fig3l_hif_attribution.pdf")
if (file.exists(.stale_pdf)) file.remove(.stale_pdf)

save_overview(
  fig3l, STAGE, "fig3l_hif_attribution",
  table = plot_df[, c("target", "module", "contrib", "direction")],
  # The finding the panel no longer states on its canvas lands here, in full and
  # with its numbers, which is where the results contract asks for it.
  finding = sprintf(paste0(
    "Regrouping Hif1a's WT_heat target contributions by curated module splits the ",
    "positive signal into a heat-shock/stress fraction summing to %+.2f, all seven ",
    "members up, against a HIF-annotated remainder of %+.2f (shared angio/glucose ",
    "%+.2f, autoregulatory feedback %+.2f, HIF1a-selective hypoxic core %+.2f). ",
    "Pdk1, Bnip3, Bnip3l and Car9 are canonical HIF1a-induced targets, yet all ",
    "four sit left of zero on the contrast where Hif1a scores positive. The ",
    "positive score is thus consistent with a heat-induced stress/glycolytic ",
    "program overlapping the HIF1a regulon, not with canonical hypoxic HIF1a ",
    "output. Correlative: attribution of an inferred activity score, not evidence ",
    "about HIF1a protein. %s"),
    stress_sum, net_hif_core,
    mod_sum[["shared_angio_glucose"]], mod_sum[["autoreg_feedback"]], core_sum,
    cap),
  script = SCRIPT, fn = "save_overview",
  config_kv = "figures.title_size=16; figures.subtitle_size=11; MODULE_COLORS; MODULE_LABELS",
  input = "03_results/04_tf/tables/fig3l_hif_attribution_data.csv",
  how_to_read = paste0(
    "Horizontal lollipop, one row per gene. x = the signed contribution ",
    "sign(mor) x t_wt this gene makes to Hif1a's WT_heat ULM score; the grey rule ",
    "marks zero, so points right of it are genes the contrast moves up and points ",
    "left of it genes it moves down, on a symmetric scale. Facets are the four ",
    "curated modules, ordered heat-shock/stress, shared angio/glucose, ",
    "autoregulatory feedback, HIF1a-selective hypoxic core; each strip carries ",
    "that module's gene count and summed contribution. Colour encodes module only ",
    "(config.R::MODULE_COLORS). Per-gene isoform attributions and the 339 ",
    "unclassified regulon members sit in fig3l_hif_attribution_data.csv. Claim ",
    "tier: descriptive attribution; no statistics are computed here. ", carve_note),
  width = FIG_W, height = FIG_H, config = FIG_CFG)

cat("[DONE] 03c_hif_program_attribution_viz.R -- fig3l rendered (no statistics).\n")
