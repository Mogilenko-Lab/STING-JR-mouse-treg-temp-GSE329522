#!/usr/bin/env Rscript
# 19_hsr_decomposition_viz.R -- VIZ
# =============================================================================
# Mouse three-lens heat-response decomposition figures (stage 12_hsr_decomp).
#
# VIZ ONLY. Reads only the tables written by 19_hsr_decomposition.R; computes no
# enrichment, overlap, rank, or attribution statistics. Figures are saved through
# the project figure-style contract (save_figure), and source tables are written
# as flat same-stem CSV neighbors under 03_results/12_hsr_decomp/tables/.
#
# Honest ceiling: HSR_core is proteotoxic-stress-general, not fever-specific.
# The 37/39 contrast is confirmatory for this experimental perturbation but does
# not make WT_heat_up causal for fever. Use correlative language only.
#
# Run from project root:
#   Rscript 02_analysis/scripts/19_hsr_decomposition_viz.R
# =============================================================================

source("02_analysis/helpers/figure_style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})
options(stringsAsFactors = FALSE)

STAGE <- "12_hsr_decomp"
SCRIPT <- "02_analysis/scripts/19_hsr_decomposition_viz.R"
TBL_DIR <- file.path("03_results", STAGE, "tables")
FIG_DIR <- file.path("03_results", STAGE, "figures")
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ATTR_LEVELS <- c("thermal_HSR", "activation", "shared_both", "neither")
ATTR_LABELS <- c(
  thermal_HSR = "Thermal HSR",
  activation = "Activation",
  shared_both = "Shared",
  neither = "Neither"
)
CONTRAST_LEVELS <- c("WT_heat", "Temp_main", "KO_heat")
TERM_LEVELS <- c("HSR_core", "TCR_activation")
HONEST_CEILING <- paste(
  "HSR core is proteotoxic-stress-general, not fever-specific;",
  "read fever/thermal-cause language as correlative."
)

need <- file.path(TBL_DIR, c("hsr_decomp_summary.csv", "hsr_decomp_lens_nes.csv"))
for (p in need)
  if (!file.exists(p)) stop("[19_viz] Missing table: ", p, " -- run 19_hsr_decomposition.R first.")

summary_df <- readr::read_csv(file.path(TBL_DIR, "hsr_decomp_summary.csv"),
                              show_col_types = FALSE, progress = FALSE)
lens_nes <- readr::read_csv(file.path(TBL_DIR, "hsr_decomp_lens_nes.csv"),
                            show_col_types = FALSE, progress = FALSE)

attr_src <- summary_df %>%
  dplyr::filter(lens == "HSR_core") %>%
  dplyr::mutate(
    attribution = factor(attribution, levels = ATTR_LEVELS),
    attribution_label = ATTR_LABELS[as.character(attribution)],
    label = sprintf("%d (%.1f%%)", n, 100 * fraction)
  ) %>%
  dplyr::arrange(attribution)
readr::write_csv(round_numeric_cols(attr_src, sig = 9),
                 file.path(TBL_DIR, "wtheatup_attribution.csv"))

fill_vals <- c(
  thermal_HSR = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  activation = FIG_CFG$colors$okabe_ito$orange %||% "#E69F00",
  shared_both = FIG_CFG$colors$okabe_ito$reddish_purple %||% "#CC79A7",
  neither = "grey70"
)

fig_attr <- ggplot(attr_src, aes(x = "WT_heat_up", y = fraction, fill = attribution)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.4) +
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5),
            size = (FIG_CFG$figures$label_size %||% 4) * 0.9,
            colour = "grey10") +
  scale_fill_manual(values = fill_vals, labels = ATTR_LABELS, name = "Attribution") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "WT_heat_up decomposes into HSR, activation, shared, and neither buckets",
    subtitle = "HSR_core lens versus disjoint TCR/IEG activation lens; empirical WT_heat_up is activation-heavy by construction.",
    x = NULL,
    y = "fraction of WT_heat_up genes",
    caption = paste0("Attribution: thermal_HSR = HSR_core only; activation = TCR_activation only; shared = both; neither = neither. ",
                     "NES/t > 0 = higher at 39 C. Claim tier: confirmatory for WT 39-vs-37 contrast; correlative for fever causality. ",
                     HONEST_CEILING)
  ) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right")

purge_figures(STAGE, "wtheatup_attribution", config = FIG_CFG)
save_figure(fig_attr, STAGE, "wtheatup_attribution",
            width = 8.5, height = 6.5, config = FIG_CFG)

nes_src <- lens_nes %>%
  dplyr::filter(term %in% TERM_LEVELS, contrast %in% CONTRAST_LEVELS) %>%
  dplyr::mutate(
    term = factor(term, levels = TERM_LEVELS),
    contrast = factor(contrast, levels = CONTRAST_LEVELS),
    direction = dplyr::case_when(nes > 0 ~ "up at 39 C", nes < 0 ~ "down at 39 C", TRUE ~ "zero"),
    label = sprintf("%.2f", nes)
  ) %>%
  dplyr::arrange(contrast, term)
readr::write_csv(round_numeric_cols(nes_src, sig = 9),
                 file.path(TBL_DIR, "lens_nes_by_contrast.csv"))

term_fills <- c(
  HSR_core = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  TCR_activation = FIG_CFG$colors$okabe_ito$orange %||% "#E69F00"
)

fig_nes <- ggplot(nes_src, aes(x = contrast, y = nes, fill = term)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey35") +
  geom_col(position = position_dodge(width = 0.78), width = 0.68) +
  geom_text(aes(label = label),
            position = position_dodge(width = 0.78),
            vjust = ifelse(nes_src$nes >= 0, -0.35, 1.2),
            size = (FIG_CFG$figures$label_size %||% 4) * 0.85,
            colour = "grey15") +
  scale_fill_manual(values = term_fills,
                    labels = c(HSR_core = "HSR core", TCR_activation = "TCR/IEG activation"),
                    name = "Lens") +
  scale_x_discrete(labels = function(x) contrast_label(x, short = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.12))) +
  labs(
    title = "Thermal HSR and activation lenses across heat-response contrasts",
    subtitle = "Grouped bars show signed GSEA NES for the curated HSR core and TCR/IEG activation lens.",
    x = NULL,
    y = "GSEA NES (signed t ranking)",
    caption = paste0("NES > 0 means the lens is enriched toward genes higher in the numerator / 39 C direction for heat contrasts. ",
                     "Claim tier: confirmatory for the temperature contrast; correlative for fever/thermal-cause language. ",
                     HONEST_CEILING)
  ) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 0),
        legend.position = "right")

purge_figures(STAGE, "lens_nes_by_contrast", config = FIG_CFG)
save_figure(fig_nes, STAGE, "lens_nes_by_contrast",
            width = 9.5, height = 6.5, config = FIG_CFG)

expected <- file.path(FIG_DIR, c(
  "wtheatup_attribution.pdf", "wtheatup_attribution.png",
  "lens_nes_by_contrast.pdf", "lens_nes_by_contrast.png"
))
missing <- expected[!file.exists(expected)]
if (length(missing) > 0L)
  stop("[19_viz] Missing expected figure output(s): ", paste(missing, collapse = ", "))

message("[19_viz] COMPLETE: wrote flat figures and same-stem source tables for stage ", STAGE)
