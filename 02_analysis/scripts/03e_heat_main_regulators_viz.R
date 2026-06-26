#!/usr/bin/env Rscript
# =============================================================================
# 03e_heat_main_regulators_viz.R - PHASE 3 (stage 04_tf) VIZ: fig3n, the HSF1 gap
# =============================================================================
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2: genotype x temperature)
# Role:    VISUALIZE half of the split. Reads ONLY the tidy table emitted by
#          03e_heat_main_regulators.R and renders the single-claim fig3n
#          lollipop. Performs NO statistics (no run_ulm/p.adjust/...); only
#          cosmetic factor ordering + plotting of already-computed columns.
#          Runs STANDALONE after compute.
#
# Input  (03_results/04_tf/tables/):  fig3n_heat_main_regulators_data.csv
# Outputs (figure-style contract, dual variants):
#   03_results/04_tf/figures/_overview/fig3n_heat_main_regulators.{print.pdf,screen.png}
#   03_results/04_tf/tables/_overview/fig3n_heat_main_regulators.csv
#   03_results/04_tf/README.md caption (via save_overview)
#
# CLAIM (fig3n): On the heat-MAIN contrasts, Hsf1 is significantly CO-ELEVATED
#   (Temp_main 3.20, padj 0.015) alongside the HIF axis (Hif1a +5.14 / Epas1
#   +4.17). Hsf1 moves WITH the HIF axis on every heat contrast yet was not
#   foregrounded earlier in this analysis. Hsf1 is colored heat-shock purple
#   (HEAT_AXIS_COLORS) to connect visually to the regulon contaminant in the
#   fig3g/fig3l arc.
#   Guardrail: co-elevation + prior neglect ONLY -- Hsf1 is NOT claimed to cause
#   or outrank HIF, and the equal-weight lollipop crowns NO single master TF.
#   Axis IS the decoupleR-ULM score, the SAME estimator as fig3a/fig3c -- it is
#   cross-quotable with them; do NOT add any GSEA-NES / different-estimator caveat.
#
# Dependencies: config.R (palettes/labels); figure_style.R (theme+save); ggplot2, dplyr
# =============================================================================

source("02_analysis/config/config.R")              # HEAT_AXIS_COLORS, AXIS_LABELS
source("02_analysis/helpers/figure_style.R")       # contract: project_theme + save_overview; FIG_CFG
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })

STAGE  <- "04_tf"
SCRIPT <- "02_analysis/scripts/03e_heat_main_regulators_viz.R"
TBL_DIR <- stage_dir("04_tf", "tables")

PROV <- provisional_caption()

rd <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)

# =============================================================================
# FIG 3n -- heat-MAIN regulators on the ULM axis, grouped by axis, faceted by contrast
# =============================================================================
fig3n_df <- rd("fig3n_heat_main_regulators_data.csv")
fig3n_df$sig      <- fig3n_df$padj < 0.05
fig3n_df$axis     <- factor(fig3n_df$axis, levels = c("heatshock", "HIF", "IFN", "other"))
fig3n_df$contrast <- factor(fig3n_df$contrast, levels = c("Temp_main", "WT_heat", "KO_heat"))

# Order TFs by axis then score so the heat-shock Hsf1 sits with the HIF block at
# the top; consistent y-order across facets (free_y not needed -- same 8 TFs).
tf_order <- fig3n_df %>%
  filter(contrast == "Temp_main") %>%
  arrange(axis, score) %>%
  pull(tf)
fig3n_df$tf <- factor(fig3n_df$tf, levels = tf_order)

# Pretty facet (contrast) labels.
contrast_labs <- c(
  Temp_main = "Temp_main (heat, both genotypes)",
  WT_heat   = "WT_heat (39 vs 37, WT)",
  KO_heat   = "KO_heat (39 vs 37, cGAS-KO)"
)

# Axis legend labels -- extend the shared AXIS_LABELS with the heat-shock entry.
heat_axis_labels <- c("heatshock" = "heat-shock axis (Hsf1)", AXIS_LABELS)

fig3n <- ggplot(fig3n_df, aes(x = score, y = tf)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.3) +
  geom_segment(aes(xend = 0, yend = tf, color = axis), linewidth = 0.7) +
  geom_point(aes(color = axis), size = 3.1) +
  geom_text(data = subset(fig3n_df, sig),
            aes(label = "*"), nudge_x = 0.18, size = 6, color = "black") +
  facet_wrap(~ contrast, ncol = 3, labeller = labeller(contrast = contrast_labs)) +
  scale_color_manual(values = HEAT_AXIS_COLORS, name = "TF axis",
                     labels = heat_axis_labels, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.14))) +
  labs(
    title = "Fig 3n. The heat-shock regulator Hsf1 is co-elevated with HIF on heat-MAIN (CollecTRI ULM)",
    subtitle = paste0(
      "Hsf1 (purple) is significantly UP on all three heat-MAIN contrasts (Temp_main 3.20, padj 0.015), co-elevated with the HIF\n",
      "axis (Hif1a +5.14 / Epas1 +4.17, orange) -- not foregrounded earlier in this analysis. Same decoupleR-ULM estimator as\n",
      "fig3a/fig3c (cross-quotable). * = BH padj < 0.05. Co-elevation only: NOT a claim that Hsf1 causes or outranks HIF."
    ),
    x = "TF activity (CollecTRI ULM score)", y = NULL
  ) +
  project_theme(config = FIG_CFG) +
  theme(panel.spacing = unit(1.1, "lines"))

save_overview(
  fig3n, STAGE, "fig3n_heat_main_regulators",
  table = fig3n_df[, c("tf", "axis", "contrast", "score", "padj", "sig")],
  finding = paste0(
    "On the heat-MAIN contrasts the heat-shock regulator Hsf1 is significantly ",
    "co-elevated (Temp_main 3.20, padj 0.015) alongside the HIF axis (Hif1a ",
    "+5.14 / Epas1 +4.17) -- co-elevation only, NOT a claim that Hsf1 causes or ",
    "outranks HIF, and no single master TF is crowned. The axis is the ",
    "decoupleR-ULM score, the same estimator as fig3a/fig3c and cross-quotable ",
    "with them. ", PROV),
  script = SCRIPT, fn = "save_overview",
  config_kv = "figures.base_size=16; figures.base_size_column=9",
  input = "03_results/04_tf/tables/fig3n_heat_main_regulators_data.csv",
  how_to_read = paste0(
    "Lollipops = per-TF activity, faceted by heat-MAIN contrast, colored by axis ",
    "(heat-shock purple = Hsf1; HIF orange; IFN blue; other grey). x = decoupleR-ULM ",
    "score (CollecTRI), the SAME estimator as fig3a/3c -- cross-quotable, no ",
    "GSEA-NES caveat. * = BH padj < 0.05. Claim tier: descriptive co-elevation; ",
    "the equal visual weight crowns no master TF."),
  config = FIG_CFG)

cat("[DONE] 03e_heat_main_regulators_viz.R -- fig3n rendered from tidy table.\n")
