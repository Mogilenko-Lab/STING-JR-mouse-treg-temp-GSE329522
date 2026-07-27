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
# The one number this script does produce is the eulerr fit residual (stress and
# diagError). That is a property of the DRAWING, not of the data: it measures how
# far the fitted circle areas fall from the counts the compute script already
# established. It is generated where the circles are generated, reported on the
# figure, and stored in that figure's same-stem table.
#
# Set-membership panels. The same seven region counts are drawn three ways so no
# single geometry has to carry the whole claim:
#   hsr_lens_membership_euler - area-proportional (eulerr), residual printed
#   hsr_lens_membership_venn  - conventional fixed layout, counts carry quantity
#   hsr_lens_membership_upset - intersection bars, empty intersections shown as 0
# All three read hsr_lens_membership.csv. Region geometry is FITTED or FIXED and
# label anchors are SOLVED, never hand-placed: a hand-laid diagram can assert a
# containment the counts contradict.
#
# Honest ceiling: HSR_core is proteotoxic-stress-general, not fever-specific.
# The 37/39 contrast is confirmatory for this experimental perturbation but does
# not make WT_heat_up causal for fever. Use correlative language only.
#
# Run from project root:
#   Rscript 02_analysis/scripts/19_hsr_decomposition_viz.R
# Missing renderers: Rscript 02_analysis/helpers/install_deps.R
# =============================================================================

source("02_analysis/helpers/figure_style.R")

suppressPackageStartupMessages({
  library(ComplexUpset)
  library(dplyr)
  library(eulerr)
  library(ggforce)
  library(ggplot2)
  library(ggridges)
  library(ggvenn)
  library(readr)
})
options(stringsAsFactors = FALSE)

STAGE <- "12_hsr_decomp"
SCRIPT <- "02_analysis/scripts/19_hsr_decomposition_viz.R"
TBL_DIR <- file.path("03_results", STAGE, "tables")
FIG_DIR <- file.path("03_results", STAGE, "figures")
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

ATTR_LEVELS <- c("HSR_core_only", "TCR_activation_only", "shared_both", "neither")
ATTR_LABELS <- c(
  HSR_core_only = "HSR_core only",
  TCR_activation_only = "TCR_activation only",
  shared_both = "Shared",
  neither = "Neither"
)
CONTRAST_LEVELS <- c("WT_heat", "Temp_main", "KO_heat")
TERM_LEVELS <- c("HSR_core", "TCR_activation")
SET_ORDER <- c("WT_heat_up", "HSR_core", "TCR_activation")
# eulerr's optimiser is iterative; pin the seed so the fitted layout is stable.
EULER_SEED <- 20260727L
# Label-anchor solver resolution (points per axis over the diagram bounding box).
POI_GRID <- 600L
HONEST_CEILING <- paste(
  "HSR_core is a curated stress-response reference;",
  "read human fever-cause language as unsupported."
)

need <- file.path(TBL_DIR, c("hsr_decomp_summary.csv", "hsr_decomp_lens_nes.csv",
                             "hsr_rank_position_panel.csv", "hsr_lens_membership.csv"))
for (p in need)
  if (!file.exists(p)) stop("[19_viz] Missing table: ", p, " -- run 19_hsr_decomposition.R first.")

summary_df <- readr::read_csv(file.path(TBL_DIR, "hsr_decomp_summary.csv"),
                              show_col_types = FALSE, progress = FALSE)
lens_nes <- readr::read_csv(file.path(TBL_DIR, "hsr_decomp_lens_nes.csv"),
                            show_col_types = FALSE, progress = FALSE)
rank_position <- readr::read_csv(file.path(TBL_DIR, "hsr_rank_position_panel.csv"),
                                 show_col_types = FALSE, progress = FALSE)
membership <- readr::read_csv(file.path(TBL_DIR, "hsr_lens_membership.csv"),
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
  HSR_core_only = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  TCR_activation_only = FIG_CFG$colors$okabe_ito$orange %||% "#E69F00",
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
    caption = paste0("Attribution: HSR_core_only = HSR_core only; TCR_activation_only = TCR_activation only; shared = both; neither = neither. ",
                     "NES/t > 0 = higher at 39 °C. Claim tier: confirmatory for WT 39-vs-37 contrast; correlative for fever causality. ",
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
    direction = dplyr::case_when(nes > 0 ~ "up at 39 °C", nes < 0 ~ "down at 39 °C", TRUE ~ "zero"),
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
    title = "HSR_core and TCR_activation lenses across 39 °C response contrasts",
    subtitle = "Grouped bars show signed GSEA NES for the curated HSR core and TCR/IEG activation lens.",
    x = NULL,
    y = "GSEA NES (signed t ranking)",
    caption = paste0("NES > 0 means the lens is enriched toward genes higher in the numerator / 39 °C direction for heat contrasts. ",
                     "Claim tier: confirmatory for the temperature contrast; correlative for fever/thermal-cause language. ",
                     HONEST_CEILING)
  ) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 0),
        legend.position = "right")

purge_figures(STAGE, "lens_nes_by_contrast", config = FIG_CFG)
save_figure(fig_nes, STAGE, "lens_nes_by_contrast",
            width = 9.5, height = 6.5, config = FIG_CFG)

rank_order <- c("TCR_activation", "HSR_core", "WT_heat_up")
rank_labels <- rank_position %>%
  dplyr::distinct(set, set_label) %>%
  dplyr::mutate(set = factor(set, levels = rank_order)) %>%
  dplyr::arrange(set)
rank_position <- rank_position %>%
  dplyr::mutate(set = factor(set, levels = rank_order),
                set_label = factor(set_label, levels = rank_labels$set_label))
rank_summary <- rank_position %>%
  dplyr::distinct(set, set_label, n_in_ranking, median_rank_pct, median_t,
                  gate_min_rank_pct, gate_max_rank_pct, gate_pass_n,
                  gate_inner_total_n, gate_deepest_gene, gate_deepest_t,
                  gate_deepest_logFC, ranking_total_n)
gate_meta <- rank_summary %>% dplyr::slice(1)

set_cols <- c(
  WT_heat_up = FIG_CFG$colors$diverging$up %||% "#B35806",
  HSR_core = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  TCR_activation = FIG_CFG$colors$okabe_ito$orange %||% "#E69F00"
)

fig_rank <- ggplot(rank_position, aes(x = rank_pct, y = set_label)) +
  annotate("rect", xmin = gate_meta$gate_min_rank_pct[1], xmax = gate_meta$gate_max_rank_pct[1],
           ymin = -Inf, ymax = Inf, alpha = 0.12,
           fill = FIG_CFG$colors$okabe_ito$yellow) +
  geom_vline(xintercept = gate_meta$gate_max_rank_pct[1], linetype = "dashed",
             linewidth = 0.45, colour = "grey30") +
  geom_density_ridges(aes(fill = set), alpha = 0.35, scale = 0.85,
                      rel_min_height = 0.01, colour = "grey25", linewidth = 0.35) +
  geom_point(aes(colour = set), position = position_jitter(height = 0.055, width = 0),
             size = 1.35, alpha = 0.75) +
  geom_text(data = rank_summary,
            aes(x = 99, y = set_label,
                label = sprintf("median %.1f%%\nt = %.2f", median_rank_pct, median_t)),
            inherit.aes = FALSE, hjust = 1, size = (FIG_CFG$figures$label_size %||% 4) * 0.72,
            colour = "grey15") +
  annotate("text", x = 1.1, y = 3.45,
           label = sprintf("gate span\n%d of %d passed",
                           gate_meta$gate_pass_n[1], gate_meta$gate_inner_total_n[1]),
           hjust = 0, size = (FIG_CFG$figures$label_size %||% 4) * 0.72, colour = "grey20") +
  annotate("text", x = gate_meta$gate_max_rank_pct[1] + 1.1, y = 0.72,
           label = sprintf("%s\nrank %d; t %.2f; logFC %.2f",
                           gate_meta$gate_deepest_gene[1], gate_meta$gate_inner_total_n[1],
                           gate_meta$gate_deepest_t[1], gate_meta$gate_deepest_logFC[1]),
           hjust = 0, size = (FIG_CFG$figures$label_size %||% 4) * 0.66, colour = "grey20") +
  scale_colour_manual(values = set_cols, guide = "none") +
  scale_fill_manual(values = set_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title = "WT_heat rank positions for the gate output and two curated lenses",
    subtitle = sprintf("The WT_heat signed-t ranking contains %s genes; the shaded band is the empirical span of the export gate.",
                       format(gate_meta$ranking_total_n[1], big.mark = ",")),
    x = "rank percentile in WT_heat signed-t ranking",
    y = NULL,
    caption = paste0("Each point is one gene; ridgelines show density along the same rank-percentile axis.\n",
                     "Shaded band = empirical gate span; only 213 of the first 2,010 ranked genes passed because logFC >= 1 was also required.")
  ) +
  project_theme(config = FIG_CFG)

purge_figures(STAGE, "hsr_rank_position_panel", config = FIG_CFG)
save_figure(fig_rank, STAGE, "hsr_rank_position_panel",
            width = 11.5, height = 6.5, config = FIG_CFG)

# =============================================================================
# Three-set membership: WT_heat_up x HSR_core x TCR_activation, drawn three ways.
# =============================================================================

# ---- 1. Resolve each region of hsr_lens_membership.csv to a set-membership triple.
# The compute script names regions in prose ("WT_heat_up only", "WT_heat_up n
# HSR_core", "shared by all three"), so read membership off the name, and take the
# all-three row from its glyph rather than from a string that names no set.
regions <- data.frame(
  component = as.character(membership$component),
  n = as.integer(membership$n_mouse),
  glyph = as.character(membership$glyph),
  stringsAsFactors = FALSE
)
is_triple <- regions$glyph == "empty_triple"
for (s in SET_ORDER)
  regions[[s]] <- is_triple | grepl(s, regions$component, fixed = TRUE)
member_mat <- as.matrix(regions[, SET_ORDER, drop = FALSE])
regions$degree <- as.integer(rowSums(member_mat))
regions$is_empty <- regions$n == 0L
regions$combo <- apply(member_mat, 1, function(r) paste(SET_ORDER[r], collapse = "&"))

# The seven regions must add back up to the three marginal set sizes the compute
# script recorded. This is the check the previous hand-laid panel did not have:
# it is what fails loudly if a region is ever mis-assigned to a set.
marginals <- c(WT_heat_up = membership$wt_heat_up_n[1],
               HSR_core = membership$hsr_core_n[1],
               TCR_activation = membership$tcr_activation_n[1])
for (s in SET_ORDER) {
  got <- sum(regions$n[regions[[s]]])
  if (!isTRUE(all.equal(as.numeric(got), as.numeric(marginals[[s]]))))
    stop(sprintf("[19_viz] Region counts for %s sum to %d but the table records %s.",
                 s, got, format(marginals[[s]])))
}
if (nrow(regions) != 7L || sum(regions$degree == 0L) != 0L)
  stop("[19_viz] hsr_lens_membership.csv must carry exactly the seven non-empty-degree regions.")

set_cols_membership <- c(
  WT_heat_up = FIG_CFG$colors$diverging$up %||% "#B35806",
  HSR_core = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  TCR_activation = FIG_CFG$colors$okabe_ito$orange %||% "#E69F00"
)
pretty_region <- function(component) gsub(" ∩ ", " with ", component, fixed = TRUE)

# Emitted widths for the three membership panels, declared once so save_figure()
# and the caption wrapper below cannot drift apart.
FIG_W_EULER <- 11.5
FIG_W_VENN <- 11.0
FIG_W_UPSET <- 11.0
# project_theme() sets the caption tier's size and lineheight but cannot re-flow
# text. coord_fixed() letterboxes the Euler and Venn panels, and the caption is
# laid out against the PANEL, so the usable caption width is narrower than the
# emitted figure width and differs per panel. CAPTION_WRAP is set to fit inside
# the narrowest of the three rather than tracking each one separately.
CAPTION_WRAP <- 100L
wrap_caption <- function(...) paste(strwrap(paste0(...), width = CAPTION_WRAP), collapse = "\n")

# Human counterpart, carried only when the external compartment table was readable.
human_ok <- is.finite(membership$human_wt_hsr_intersect[1]) &&
  is.finite(membership$human_wt_heatup_n[1])
human_line <- if (human_ok) {
  sprintf("Human counterpart after ortholog projection: %d of %d WT_heat_up genes fall in the curated HSR core.",
          as.integer(membership$human_wt_hsr_intersect[1]),
          as.integer(membership$human_wt_heatup_n[1]))
} else {
  "Human counterpart unavailable: the JIA sorted-Treg overlap table was absent or malformed."
}
cond_line <- sprintf(
  "Conditioning either curated lens on the other moves NES by about %+.3f.",
  mean(c(membership$conditional_delta_nes_HSR_core[1],
         membership$conditional_delta_nes_TCR_activation[1])))
membership_ceiling <- paste(
  "Set membership only: a lens enriching is a separate measurement from that lens being",
  "contained in WT_heat_up.", HONEST_CEILING)

# ---- 2. Solve a label anchor for every region, empty ones included.
# For a point p and circle i, margin_i(p) = r_i - |p - c_i| when the region requires
# being INSIDE circle i, and |p - c_i| - r_i when it requires being OUTSIDE. The
# anchor is the p maximising min_i margin_i. Where the region exists that is its
# pole of inaccessibility and the margin is the inscribed radius; where the region
# is empty the margin goes negative and the anchor is the point that comes closest
# to satisfying the region's definition -- i.e. exactly where the count would sit
# if it were not zero. One rule, seven regions, no typed coordinates.
solve_anchors <- function(circles, regions_df, n_grid = POI_GRID) {
  pad <- 0.08 * max(circles$r)
  gx <- seq(min(circles$x - circles$r) - pad, max(circles$x + circles$r) + pad, length.out = n_grid)
  gy <- seq(min(circles$y - circles$r) - pad, max(circles$y + circles$r) + pad, length.out = n_grid)
  g <- expand.grid(x = gx, y = gy, KEEP.OUT.ATTRS = FALSE)
  inside_margin <- vapply(seq_len(nrow(circles)), function(i)
    circles$r[i] - sqrt((g$x - circles$x[i])^2 + (g$y - circles$y[i])^2),
    numeric(nrow(g)))
  colnames(inside_margin) <- circles$set
  out <- lapply(seq_len(nrow(regions_df)), function(k) {
    m <- vapply(SET_ORDER, function(s)
      if (isTRUE(regions_df[[s]][k])) inside_margin[, s] else -inside_margin[, s],
      numeric(nrow(g)))
    worst <- do.call(pmin, as.data.frame(m))
    j <- which.max(worst)
    data.frame(anchor_x = g$x[j], anchor_y = g$y[j], anchor_margin = worst[j])
  })
  dplyr::bind_cols(regions_df, dplyr::bind_rows(out))
}

# ---- 3. PANEL A -- area-proportional Euler, with the fit residual reported.
set.seed(EULER_SEED)
euler_fit <- eulerr::euler(stats::setNames(regions$n, regions$combo), shape = "circle")
euler_stress <- as.numeric(euler_fit$stress)
euler_diag_error <- as.numeric(euler_fit$diagError)

circles <- data.frame(
  set = rownames(euler_fit$ellipses),
  x = euler_fit$ellipses$h,
  y = euler_fit$ellipses$k,
  r = euler_fit$ellipses$a,
  stringsAsFactors = FALSE
)
if (!setequal(circles$set, SET_ORDER))
  stop("[19_viz] eulerr returned circles for an unexpected set of names.")
circles <- circles[match(SET_ORDER, circles$set), , drop = FALSE]

euler_src <- solve_anchors(circles, regions) %>%
  dplyr::left_join(
    data.frame(
      combo = names(euler_fit$original.values),
      original_value = as.numeric(euler_fit$original.values),
      fitted_value = as.numeric(euler_fit$fitted.values),
      stringsAsFactors = FALSE),
    by = "combo") %>%
  dplyr::mutate(
    residual = .data$fitted_value - .data$original_value,
    region_exists = .data$anchor_margin > 0,
    euler_stress = euler_stress,
    euler_diag_error = euler_diag_error,
    label = pretty_region(.data$component))

# Degree-1 regions are large enough to label in place; the overlaps carry the
# load-bearing numbers and get named callouts outside the diagram, so the reader
# never has to infer which sliver a bare digit belongs to. Radial callout placement
# collides here -- the two empty regions sit at almost the same bearing -- so the
# callouts are stacked in a right-hand column, ordered by anchor height, which
# separates them and keeps the leader lines from crossing.
bbox_x <- range(c(circles$x - circles$r, circles$x + circles$r))
bbox_y <- range(c(circles$y - circles$r, circles$y + circles$r))
span_x <- diff(bbox_x)
span_y <- diff(bbox_y)

euler_src$placement <- ifelse(euler_src$degree == 1L, "inside", "callout")
euler_src$callout_x <- NA_real_
euler_src$callout_y <- NA_real_
call_idx <- which(euler_src$placement == "callout")
call_idx <- call_idx[order(euler_src$anchor_y[call_idx], decreasing = TRUE)]
euler_src$callout_x[call_idx] <- bbox_x[2] + 0.10 * span_x
euler_src$callout_y[call_idx] <- seq(bbox_y[2] - 0.04 * span_y,
                                     bbox_y[1] + 0.04 * span_y,
                                     length.out = length(call_idx))

euler_inside <- euler_src[euler_src$placement == "inside", , drop = FALSE]
euler_callout <- euler_src[euler_src$placement == "callout", , drop = FALSE]

readr::write_csv(round_numeric_cols(
  euler_src %>% dplyr::select("component", "n", "degree", "is_empty", "combo",
                              dplyr::all_of(SET_ORDER), "original_value", "fitted_value",
                              "residual", "anchor_x", "anchor_y", "anchor_margin",
                              "region_exists", "euler_stress", "euler_diag_error"),
  sig = 9), file.path(TBL_DIR, "hsr_lens_membership_euler.csv"))

fig_euler <- ggplot() +
  ggforce::geom_circle(data = circles, aes(x0 = x, y0 = y, r = r, fill = set, colour = set),
                       alpha = 0.18, linewidth = 0.8) +
  geom_segment(data = euler_callout,
               aes(x = anchor_x, y = anchor_y, xend = callout_x, yend = callout_y),
               colour = "grey45", linewidth = 0.35) +
  geom_point(data = euler_callout, aes(x = anchor_x, y = anchor_y, shape = is_empty),
             size = 2.6, stroke = 1.0, colour = "grey20", fill = "white") +
  geom_label(data = euler_callout,
             aes(x = callout_x, y = callout_y,
                 label = sprintf("%s\n%d %s", label, n,
                                 ifelse(is_empty, "genes - EMPTY", "genes"))),
             hjust = 0, size = (FIG_CFG$figures$label_size %||% 4) * 0.78, lineheight = 0.95,
             colour = "grey10", fill = "white", linewidth = 0.25,
             label.padding = unit(0.18, "lines")) +
  geom_text(data = euler_inside, aes(x = anchor_x, y = anchor_y, label = n),
            size = (FIG_CFG$figures$label_size %||% 4) * 1.05, fontface = "bold",
            colour = "grey10") +
  scale_fill_manual(values = set_cols_membership, breaks = SET_ORDER, name = NULL) +
  scale_colour_manual(values = set_cols_membership, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") +
  coord_fixed(xlim = c(bbox_x[1] - 0.04 * span_x, bbox_x[2] + 0.62 * span_x),
              ylim = c(bbox_y[1] - 0.06 * span_y, bbox_y[2] + 0.06 * span_y),
              expand = FALSE, clip = "off") +
  labs(
    title = "Area-proportional membership of the 39 °C-derived up arm and two curated lenses",
    subtitle = sprintf(paste("Circle areas are fitted to the counts; eulerr residual stress = %.2g,",
                             "diagError = %.2g (0 = an exact fit)."),
                       euler_stress, euler_diag_error),
    x = NULL, y = NULL,
    caption = wrap_caption(
      "Bold numbers sit inside the region they count; each callout names an overlap and its leader line ",
      "points at it. Open circles mark the two overlaps that exist; open diamonds mark the two EMPTY ",
      "intersections, anchored where the count would sit if it were not zero. ",
      human_line, " ", cond_line, " ", membership_ceiling)
  ) +
  project_theme(config = FIG_CFG) +
  # coord_fixed() letterboxes the panel; anchor the caption to the PLOT edge so it
  # gets the full emitted width the wrap above was computed against.
  theme(legend.position = "bottom", plot.caption.position = "plot")

# ---- 4. PANEL B -- conventional fixed-layout Venn; only the printed counts scale.
# Expand the region counts back to one row per gene so the renderer re-derives the
# seven region sizes itself; the stopifnot is that round trip closing.
venn_df <- do.call(rbind, lapply(seq_len(nrow(regions)), function(k) {
  if (regions$n[k] == 0L) return(NULL)
  matrix(rep(unlist(regions[k, SET_ORDER]), regions$n[k]),
         ncol = length(SET_ORDER), byrow = TRUE)
}))
venn_df <- as.data.frame(venn_df)
names(venn_df) <- SET_ORDER
stopifnot(nrow(venn_df) == sum(regions$n),
          all(colSums(venn_df) == marginals[SET_ORDER]))

venn_src <- regions %>%
  dplyr::select("component", "n", "degree", "is_empty", "combo", dplyr::all_of(SET_ORDER)) %>%
  dplyr::mutate(layout = "conventional fixed three-circle Venn",
                area_is_quantitative = FALSE,
                drawn_region_count = nrow(regions))
readr::write_csv(round_numeric_cols(venn_src, sig = 9),
                 file.path(TBL_DIR, "hsr_lens_membership_venn.csv"))

fig_venn <- ggvenn::ggvenn(
  venn_df, columns = SET_ORDER, show_percentage = FALSE,
  fill_color = unname(set_cols_membership[SET_ORDER]), fill_alpha = 0.28,
  stroke_color = "grey25", stroke_size = 0.6,
  text_size = (FIG_CFG$figures$label_size %||% 4) * 1.15,
  set_name_size = (FIG_CFG$figures$label_size %||% 4) * 1.0) +
  coord_fixed(clip = "off") +
  labs(
    title = "The same seven counts, drawn where geometry carries no quantity",
    subtitle = paste("Fixed-layout Venn: every region is drawn, so a printed 0 marks a region that",
                     "exists and is empty."),
    x = NULL, y = NULL,
    caption = wrap_caption(
      "All three circles are equal and equally overlapping by convention and only the printed counts ",
      "vary, so this layout cannot assert a containment the counts contradict. HSR_core with ",
      "TCR_activation and the three-way region both read 0. ",
      human_line, " ", cond_line, " ", membership_ceiling)
  ) +
  project_theme(config = FIG_CFG) +
  theme(legend.position = "none", plot.caption.position = "plot")

# ---- 5. PANEL C -- UpSet over the same three sets, empty intersections included.
# ComplexUpset drops a zero-size intersection from its label layer, so the two
# empty combinations are annotated explicitly at their own axis positions; their
# order is read from `upset_order`, not typed in.
upset_order <- regions %>% dplyr::arrange(dplyr::desc(.data$n), .data$degree, .data$component)
upset_intersections <- lapply(seq_len(nrow(upset_order)), function(k)
  SET_ORDER[unlist(upset_order[k, SET_ORDER])])
zero_pos <- which(upset_order$n == 0L)

upset_src <- upset_order %>%
  dplyr::select("component", "n", "degree", "is_empty", "combo", dplyr::all_of(SET_ORDER)) %>%
  dplyr::mutate(axis_position = dplyr::row_number(),
                set_size_WT_heat_up = marginals[["WT_heat_up"]],
                set_size_HSR_core = marginals[["HSR_core"]],
                set_size_TCR_activation = marginals[["TCR_activation"]])
readr::write_csv(round_numeric_cols(upset_src, sig = 9),
                 file.path(TBL_DIR, "hsr_lens_membership_upset.csv"))

upset_label_size <- (FIG_CFG$figures$label_size %||% 4) * 0.95
fig_upset <- ComplexUpset::upset(
  venn_df, intersect = SET_ORDER,
  intersections = upset_intersections,
  sort_intersections = FALSE, sort_sets = FALSE,
  keep_empty_groups = TRUE, min_size = 0,
  name = "intersection of the three sets",
  height_ratio = 0.42, width_ratio = 0.26,
  matrix = ComplexUpset::intersection_matrix(
    geom = geom_point(size = (FIG_CFG$figures$point_size %||% 2.4) * 1.15)),
  base_annotations = list(
    "genes in intersection" = ComplexUpset::intersection_size(
      text = list(vjust = -0.4, size = upset_label_size),
      text_colors = c(on_background = "grey10", on_bar = "white")) +
      annotate("text", x = zero_pos, y = 0, label = "0", vjust = -0.6,
               size = upset_label_size, colour = "grey10") +
      annotate("text", x = zero_pos, y = 0, label = "empty", vjust = -2.1,
               size = upset_label_size * 0.8, colour = "grey40")),
  themes = ComplexUpset::upset_default_themes(
    text = element_text(size = FIG_CFG$figures$base_size %||% 14),
    axis.text = element_text(size = FIG_CFG$figures$axis_text_size %||% 11),
    axis.title = element_text(size = FIG_CFG$figures$axis_title_size %||% 13))) +
  # upset() returns a patchwork: a bare labs() would land on the last sub-panel,
  # so the figure-level text goes through plot_annotation with the project theme.
  patchwork::plot_annotation(
    title = "UpSet view: both empty intersections are drawn as zero, not left out",
    subtitle = paste("Bars are intersection sizes over the same three sets; the two zero-height bars",
                     "are labelled rather than omitted."),
    caption = wrap_caption(
      "Dots mark which sets each bar belongs to; the left bars are the three set sizes. Every ",
      "combination of the three sets gets a column, so an intersection with no genes appears as a ",
      "zero-height bar labelled 'empty' rather than as a missing column. ",
      human_line, " ", cond_line, " ", membership_ceiling),
    theme = project_theme(config = FIG_CFG) + theme(plot.caption.position = "plot"))

purge_figures(STAGE, "hsr_lens_membership", config = FIG_CFG)
save_figure(fig_euler, STAGE, "hsr_lens_membership_euler",
            width = FIG_W_EULER, height = 7.5, void = TRUE, config = FIG_CFG)
save_figure(fig_venn, STAGE, "hsr_lens_membership_venn",
            width = FIG_W_VENN, height = 7.5, void = TRUE, config = FIG_CFG)
save_figure(fig_upset, STAGE, "hsr_lens_membership_upset",
            width = FIG_W_UPSET, height = 7.0, config = FIG_CFG)

message(sprintf("[19_viz] eulerr fit: stress = %.3g, diagError = %.3g, max |region residual| = %.3g",
                euler_stress, euler_diag_error, max(abs(euler_src$residual))))

expected <- file.path(FIG_DIR, c(
  "wtheatup_attribution.pdf", "wtheatup_attribution.png",
  "lens_nes_by_contrast.pdf", "lens_nes_by_contrast.png",
  "hsr_rank_position_panel.pdf", "hsr_rank_position_panel.png",
  "hsr_lens_membership_euler.pdf", "hsr_lens_membership_euler.png",
  "hsr_lens_membership_venn.pdf", "hsr_lens_membership_venn.png",
  "hsr_lens_membership_upset.pdf", "hsr_lens_membership_upset.png"
))
missing <- expected[!file.exists(expected)]
if (length(missing) > 0L)
  stop("[19_viz] Missing expected figure output(s): ", paste(missing, collapse = ", "))

message("[19_viz] COMPLETE: wrote flat figures and same-stem source tables for stage ", STAGE)
