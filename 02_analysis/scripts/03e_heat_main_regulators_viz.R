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
# Output (03_results/04_tf/figures/): fig3n_heat_main_regulators.pdf
#          (+ deck_assets/fig3n_heat_main_regulators.png when DECK_EXPORT set)
#
# CLAIM (fig3n): On the heat-MAIN contrasts, Hsf1 is significantly CO-ELEVATED
#   (Temp_main 3.20, padj 0.015) alongside the HIF axis (Hif1a +5.14 / Epas1
#   +4.17). Hsf1 moves WITH the HIF axis on every heat contrast yet was not
#   foregrounded earlier in this analysis. Hsf1 is colored heat-shock purple
#   (HEAT_AXIS_COLORS) to connect visually to the regulon contaminant in the
#   fig3g/fig3l arc.
#   Guardrail: co-elevation + prior neglect ONLY -- Hsf1 is NOT claimed to cause
#   or outrank HIF. Axis is decoupleR-ULM, the SAME estimator as fig3a/fig3c.
#
# Dependencies: config.R; ggplot2, dplyr
# =============================================================================

source("02_analysis/config/config.R")
load_packages()  # ggplot2 + dplyr are what we need

TBL_DIR <- stage_dir("04_tf", "tables")
FIG_DIR <- stage_dir("04_tf", "figures")

cap <- provisional_caption()

rd <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# DECK EXPORT: when DECK_EXPORT is set, ALSO write a clean 300-dpi PNG (no embargo
# caption line) to 03_results/04_tf/deck_assets/. Canonical stamped PDF is ALWAYS
# written, unchanged. Matches 03_decoupler_tf_viz.R's save_fig() contract.
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

# House theme -- identical to 03_decoupler_tf_viz.R.
base_theme <- theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title    = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9, color = "grey25"),
        plot.caption  = element_text(size = 7, color = "grey45", hjust = 0, face = "italic"))

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
            aes(label = "*"), nudge_x = 0.18, size = 5, color = "black") +
  facet_wrap(~ contrast, ncol = 3, labeller = labeller(contrast = contrast_labs)) +
  scale_color_manual(values = HEAT_AXIS_COLORS, name = "TF axis",
                     labels = heat_axis_labels, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.14))) +
  labs(
    title = "Fig 3n. The heat-shock regulator Hsf1 is co-elevated with HIF on heat-MAIN (CollecTRI ULM)",
    subtitle = paste0(
      "Hsf1 (purple) is significantly UP on all three heat-MAIN contrasts (Temp_main 3.20, padj 0.015), co-elevated with the HIF\n",
      "axis (Hif1a +5.14 / Epas1 +4.17, orange) -- not foregrounded earlier in this analysis. Same decoupleR-ULM estimator as\n",
      "fig3a/fig3c. * = BH padj < 0.05. Co-elevation only: NOT a claim that Hsf1 causes or outranks HIF."
    ),
    x = "TF activity (CollecTRI ULM score)", y = NULL,
    caption = cap
  ) +
  base_theme +
  theme(strip.text  = element_text(face = "bold", size = 8.5),
        axis.text.y = element_text(size = 8.5),
        panel.spacing = unit(1.1, "lines"))

save_fig(fig3n, "fig3n_heat_main_regulators.pdf", width = 11, height = 5.5, deck_h = 5.0)

cat("[DONE] 03e_heat_main_regulators_viz.R -- fig3n rendered from tidy table.\n")
