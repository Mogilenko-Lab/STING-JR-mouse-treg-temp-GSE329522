#!/usr/bin/env Rscript
# =============================================================================
# 03f_experiment_tiers_viz.R - PHASE 3 (stage 04_tf) VIZ: fig3o
#   The conceptual tier-ladder infographic -- "what to run next".
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# Role:    VISUALIZE half of the split for fig3o. Reads ONLY the hand-authored
#          config table emitted by 03f_experiment_tiers.R and renders a clean
#          conceptual LADDER (one tier per ascending band, boxed arms, gate
#          arrows). This is NOT a data plot -- it carries NO statistics
#          (no run_ulm/run_mlm/p.adjust/prcomp/cor); it draws fixed boxes/arrows
#          from the lookup table. Stamped "CONCEPTUAL - not data" prominently.
#          Runs STANDALONE after the compute script.
#
# Input  (03_results/04_tf/tables/):
#   - fig3o_experiment_tiers_data.csv
# Outputs:
#   - 03_results/04_tf/figures/fig3o_experiment_tiers.pdf      (canonical, stamped)
#   - 03_results/04_tf/deck_assets/fig3o_experiment_tiers.png  (DECK_EXPORT=1; clean)
#
# Framing constraints (bind the figure text): necessity = "some activity the drug
#   perturbs is required at 39C", NEVER "HIF1a is required"; belzutifan = a HIF2a-
#   selective cross-check tool the drug can't exclude, never a crowning of HIF2a;
#   Tier 2 genetics GATED, never led-with.
# Dependencies: config.R; ggplot2, dplyr (no statistics)
# =============================================================================

source("02_analysis/config/config.R")
load_packages()   # ggplot2 + dplyr; limma/edgeR pulled in but unused

TBL_DIR <- stage_dir("04_tf", "tables")
FIG_DIR <- stage_dir("04_tf", "figures")

cap <- provisional_caption()
rd  <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# DECK EXPORT: when DECK_EXPORT is set, ALSO write a clean 300-dpi PNG (no
# embargo caption) to deck_assets/. The stamped PDF is always written.
# -----------------------------------------------------------------------------
DECK_EXPORT <- nzchar(Sys.getenv("DECK_EXPORT"))
DECK_DIR    <- file.path(DIR_RESULTS, "04_tf", "deck_assets")
if (DECK_EXPORT) dir.create(DECK_DIR, recursive = TRUE, showWarnings = FALSE)

save_fig <- function(plot, fname, width, height, deck_h = 6.5) {
  ggsave(file.path(FIG_DIR, fname), plot, width = width, height = height)
  cat(sprintf("  wrote %s\n", fname))
  if (DECK_EXPORT) {
    png_name <- sub("\\.pdf$", ".png", fname)
    ggsave(file.path(DECK_DIR, png_name), plot + labs(caption = NULL),
           width = width, height = deck_h, dpi = 300)
    cat(sprintf("  wrote deck_assets/%s (clean, no embargo stamp)\n", png_name))
  }
}

# =============================================================================
# READ the config table
# =============================================================================
dat <- rd("fig3o_experiment_tiers_data.csv")
dat$tier <- as.integer(dat$tier)

# Tier-band colors: ascending cost = ascending saturation. Reuse the house
# diverging/family hues so the ladder reads "cheap -> expensive" by warmth.
tier_band_cols <- c("0" = "#E6F2EC", "1" = "#FCEBD2", "2" = "#F3DCE8")  # mint -> cream -> mauve
tier_edge_cols <- c("0" = "#1B7837", "1" = "#E08214", "2" = "#762A83")  # green -> orange -> purple
cost_label     <- c("0" = "cost: cheap", "1" = "cost: medium", "2" = "cost: expensive (gated)")

# -----------------------------------------------------------------------------
# LADDER GEOMETRY (fixed coordinates; no computation on data values).
# Three stacked tier bands; Tier 0 at the BOTTOM (start here) ascending to Tier 2
# at the TOP (last). Within each band, arms are boxes laid left-to-right.
# -----------------------------------------------------------------------------
n_arm_max <- max(table(dat$tier))
band_h    <- 2.20         # vertical height of each tier band (room for header/boxes/gate)
band_gap  <- 0.55         # gap between bands (room for the ascent arrows)
box_w     <- 0.92         # arm-box width (within a unit cell)
box_h     <- 1.28
# Within a band: header sits at the TOP edge, boxes in the MIDDLE, gate at the
# BOTTOM edge -- so none of the three text layers collide with the arm boxes.
hdr_off   <- band_h/2 - 0.10   # header y offset above band center
gate_off  <- band_h/2 - 0.06   # gate y offset below band center
box_off   <- -0.16             # box center y offset (slightly below band center)

# y center of each tier band: tier 0 lowest, tier 2 highest (ascending ladder).
band_y <- function(t) t * (band_h + band_gap)

# Per-arm box frame (boxes centered in the middle band of each tier).
dat$bx_cx <- dat$row_in_tier - 1                  # 0-indexed column
dat$bx_cy <- band_y(dat$tier) + box_off
dat$xmin  <- dat$bx_cx - box_w/2
dat$xmax  <- dat$bx_cx + box_w/2
dat$ymin  <- dat$bx_cy - box_h/2
dat$ymax  <- dat$bx_cy + box_h/2
dat$tier_chr <- as.character(dat$tier)

# Wrap helper for the arm text (no stats; pure string formatting).
wrap2 <- function(x, width) vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"), character(1))
dat$arm_text <- sprintf("%s\n\n%s",
                        wrap2(dat$experiment, 26),
                        wrap2(dat$what_it_excludes_or_tests, 30))

# -----------------------------------------------------------------------------
# Tier band rectangles (one per tier) spanning all arm columns + a label strip.
# -----------------------------------------------------------------------------
x_lo <- -0.5 - 0.10
x_hi <- (n_arm_max - 1) + 0.5 + 0.10
band_df <- do.call(rbind, lapply(sort(unique(dat$tier)), function(t) {
  data.frame(
    tier = t, tier_chr = as.character(t),
    xmin = x_lo - 0.05, xmax = x_hi + 0.05,
    ymin = band_y(t) - band_h/2, ymax = band_y(t) + band_h/2,
    tier_label = dat$tier_label[dat$tier == t][1],
    cost = cost_label[as.character(t)],
    gate = dat$gate[dat$tier == t][1],
    stringsAsFactors = FALSE
  )
}))

# Ascent arrows between bands (Tier 0 -> 1 -> 2) on the left margin, with the gate.
arrow_df <- data.frame(
  x    = x_lo - 0.35,
  ylo  = c(band_y(0) + band_h/2, band_y(1) + band_h/2),
  yhi  = c(band_y(1) - band_h/2, band_y(2) - band_h/2),
  stringsAsFactors = FALSE
)

# -----------------------------------------------------------------------------
# CONCEPTUAL stamp -- unmistakable diagonal-style banner, top-right.
# -----------------------------------------------------------------------------
stamp_df <- data.frame(
  x = x_hi + 0.05, y = band_y(2) + band_h/2 + 0.30,
  lab = "CONCEPTUAL - not data", stringsAsFactors = FALSE)

# =============================================================================
# FIG 3o -- the tier ladder
# =============================================================================
fig3o <- ggplot() +
  # tier bands
  geom_rect(data = band_df,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = tier_chr),
            color = NA, alpha = 0.65) +
  # band header: tier label + cost (left-aligned, ABOVE the arm boxes)
  geom_text(data = band_df,
            aes(x = x_lo, y = ymin + band_h - 0.08,
                label = sprintf("%s   [%s]", tier_label, cost),
                color = tier_chr),
            hjust = 0, vjust = 1, fontface = "bold", size = 3.8) +
  # gate line per band (small italic, BELOW the arm boxes, spanning width)
  geom_text(data = band_df,
            aes(x = x_lo, y = ymin + 0.08, label = wrap2(gate, 95)),
            hjust = 0, vjust = 0, size = 2.6, color = "grey30", fontface = "italic") +
  # arm boxes
  geom_rect(data = dat,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = tier_chr),
            fill = "white", linewidth = 0.7) +
  geom_text(data = dat,
            aes(x = bx_cx, y = bx_cy, label = arm_text),
            size = 2.5, lineheight = 0.95, color = "grey10") +
  # ascent arrows (gated progression)
  geom_segment(data = arrow_df,
               aes(x = x, xend = x, y = ylo, yend = yhi),
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"),
               linewidth = 0.8, color = "grey35") +
  # CONCEPTUAL stamp
  geom_label(data = stamp_df, aes(x = x, y = y, label = lab),
             hjust = 1, vjust = 0, fill = "#FFF3B0", color = "#B35806",
             fontface = "bold", size = 4.2, label.size = 0.6) +
  scale_fill_manual(values = tier_band_cols, guide = "none") +
  scale_color_manual(values = tier_edge_cols, guide = "none") +
  coord_cartesian(xlim = c(x_lo - 0.55, x_hi + 0.10),
                  ylim = c(band_y(0) - band_h/2 - 0.10, band_y(2) + band_h/2 + 0.70),
                  clip = "off") +
  labs(
    title = "Fig 3o. What to run next: a cost-tiered plan that excludes the confound before naming the activity",
    subtitle = paste0(
      "The drug shows that SOME activity it perturbs is required at 39C; the deck shows the molecular label is not identifiable.\n",
      "Climb the ladder (bottom -> top): Tier 0 exclude the heat x drug-cytotoxicity confound -> Tier 1 'is it HIF at all?'\n",
      "-> Tier 2 isoform genetics. Tier 2 is GATED -- do NOT lead with floxed-Hif1a vs floxed-Epas1. (belzutifan = a\n",
      "HIF2a-selective cross-check the drug cannot exclude, not a crowning of HIF2a.)"),
    x = NULL, y = NULL,
    caption = cap
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title    = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 8.5, color = "grey25", lineheight = 1.2, hjust = 0),
    plot.caption  = element_text(size = 7, color = "grey45", hjust = 0, face = "italic"),
    plot.margin   = margin(t = 6, r = 10, b = 6, l = 8, unit = "mm")
  )

save_fig(fig3o, "fig3o_experiment_tiers.pdf", width = 12, height = 8, deck_h = 6.8)

cat("[DONE] 03f_experiment_tiers_viz.R -- fig3o conceptual tier ladder rendered (no statistics).\n")
