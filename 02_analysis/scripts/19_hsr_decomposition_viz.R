#!/usr/bin/env Rscript
# 19_hsr_decomposition_viz.R -- VIZ
# =============================================================================
# Mouse set-construction figures (stage 12_hsr_decomp): what the thresholded
# WT_heat_up set is made of, and how it hands off to the human compartments.
#
# VIZ ONLY. Reads only the tables written by 19_hsr_decomposition.R; computes no
# enrichment, overlap, rank, attribution or census statistic. Figures go to
# 03_results/12_hsr_decomp/figures/_overview/ through the project figure-style
# contract (save_figure), and each figure's same-stem source table sits beside it
# under 03_results/12_hsr_decomp/tables/_overview/.
#
# Every panel says on its own face whether it is asking a question or answering
# one, and a panel that only corroborates an answer already given says that too.
# The four questions this stage walks, in reader order:
#
#   wtheatup_attribution    ANSWERS  what is the thresholded set made of?
#                                    -> 198 of 213 genes in neither curated lens.
#   lens_nes_by_contrast    ANSWERS  so is the heat-shock response simply absent
#                                    from the 39 °C response? -> no, it is induced
#                                    in both genotypes.
#   hsr_rank_position_panel ANSWERS  why is an induced response not in the set?
#                                    -> it is a moderate mover and the gate kept
#                                    only extreme ones.
#   gate_projection_bridge  ANSWERS  how do the 213-gene mouse gate and the
#                                    199-gene human set relate, and how much of
#                                    the human set survives downstream?
#   hsr_lens_membership_*   CORROBORATE the first answer, three ways.
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
# Legibility contract. project_theme() sets type sizes but cannot re-flow text, so
# every title, subtitle and caption is wrapped here against the width the figure is
# actually emitted at. Dense caveats live in the stage README; the face carries at
# most one warning line.
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
TBL_OVW <- file.path(TBL_DIR, "_overview")
FIG_DIR <- file.path("03_results", STAGE, "figures")
FIG_OVW <- file.path(FIG_DIR, "_overview")
dir.create(TBL_OVW, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_OVW, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Text fitting. Characters that fit across an emitted figure, per type tier, at
# the sizes project_theme() sets. Derived once from the rendered panels rather
# than guessed, and deliberately conservative: an under-filled line costs nothing,
# a clipped one costs the sentence.
# ---------------------------------------------------------------------------
CHARS_PER_INCH <- c(title = 8.5, subtitle = 11, caption = 13)
fit_text <- function(..., width_in, tier = "caption") {
  w <- max(20L, as.integer(floor(width_in * CHARS_PER_INCH[[tier]])))
  paste(strwrap(paste0(...), width = w), collapse = "\n")
}
# One shared way of saying which job a panel is doing, so the reader can walk the
# stage and never wonder why a panel is there.
asks <- function(q) paste0("Question — ", q)
answers <- function(q, a) paste0("Question — ", q, "\nAnswer — ", a)
corroborates <- function(a) paste0("Corroborates — ", a)

# The stage used to emit its figures flat under figures/ and their same-stem source
# tables flat under tables/. Clear both, so a reader cannot open a panel or a table
# that no run maintains any more. Only the stems this script owns are removed; the
# compute script's own tables/ files are untouched.
OWNED_STEMS <- c("wtheatup_attribution", "lens_nes_by_contrast", "gate_projection_bridge",
                 "hsr_lens_membership_euler", "hsr_lens_membership_venn",
                 "hsr_lens_membership_upset")
for (stale in list.files(FIG_DIR, pattern = "\\.(png|pdf)$", full.names = TRUE))
  file.remove(stale)
for (stem in OWNED_STEMS) {
  stale <- file.path(TBL_DIR, paste0(stem, ".csv"))
  if (file.exists(stale)) file.remove(stale)
}

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

need <- c(file.path(TBL_DIR, c("hsr_decomp_summary.csv", "hsr_decomp_lens_nes.csv",
                               "hsr_lens_membership.csv")),
          file.path(TBL_OVW, c("hsr_rank_position_panel.csv", "gate_projection_bridge.csv")))
for (p in need)
  if (!file.exists(p)) stop("[19_viz] Missing table: ", p, " -- run 19_hsr_decomposition.R first.")

summary_df <- readr::read_csv(file.path(TBL_DIR, "hsr_decomp_summary.csv"),
                              show_col_types = FALSE, progress = FALSE)
lens_nes <- readr::read_csv(file.path(TBL_DIR, "hsr_decomp_lens_nes.csv"),
                            show_col_types = FALSE, progress = FALSE)
rank_position <- readr::read_csv(file.path(TBL_OVW, "hsr_rank_position_panel.csv"),
                                 show_col_types = FALSE, progress = FALSE)
membership <- readr::read_csv(file.path(TBL_DIR, "hsr_lens_membership.csv"),
                              show_col_types = FALSE, progress = FALSE)
bridge <- readr::read_csv(file.path(TBL_OVW, "gate_projection_bridge.csv"),
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
                 file.path(TBL_OVW, "wtheatup_attribution.csv"))

attr_n <- function(a) attr_src$n[attr_src$attribution == a][1]
attr_pct <- function(a) 100 * attr_src$fraction[attr_src$attribution == a][1]
attr_total <- attr_src$denominator[1]

fill_vals <- c(
  HSR_core_only = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  TCR_activation_only = FIG_CFG$colors$okabe_ito$orange %||% "#E69F00",
  shared_both = FIG_CFG$colors$okabe_ito$reddish_purple %||% "#CC79A7",
  neither = "grey70"
)

# The stacked bar IS the answer, so the title states the count rather than naming
# the operation. What the set is MADE OF and what the 39 °C response ENRICHES FOR
# are two different measurements: this panel is the membership one, and its
# subtitle says so, because a rank-enrichment property read onto a membership
# partition is exactly the inversion the composition-before-enrichment rule exists
# to prevent.
W_ATTR <- 11.0
fig_attr <- ggplot(attr_src, aes(x = "WT_heat_up", y = fraction, fill = attribution)) +
  geom_col(width = 0.55, colour = "white", linewidth = 0.4) +
  geom_text(data = dplyr::filter(attr_src, .data$fraction >= 0.05),
            aes(label = label), position = position_stack(vjust = 0.5),
            size = (FIG_CFG$figures$label_size %||% 4) * 1.15,
            fontface = "bold", colour = "grey10") +
  # The three small segments cannot hold text at any legible size, so they are
  # called out to the right of the bar at even spacing instead of being crowded
  # on top of one another.
  geom_text(data = dplyr::filter(attr_src, .data$fraction < 0.05) %>%
              dplyr::mutate(y_call = seq(0.985, 0.945, length.out = dplyr::n())),
            aes(y = y_call, label = sprintf("%s  %s", ATTR_LABELS[as.character(attribution)], label)),
            x = 1.34, hjust = 0, size = (FIG_CFG$figures$label_size %||% 4) * 0.95,
            colour = "grey15", inherit.aes = FALSE) +
  scale_fill_manual(values = fill_vals, labels = ATTR_LABELS, name = "Attribution") +
  scale_x_discrete(expand = expansion(add = c(0.55, 1.55))) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = fit_text(sprintf("%d of the %d genes in WT_heat_up fall in neither curated lens",
                             attr_n("neither"), attr_total),
                     width_in = W_ATTR, tier = "title"),
    subtitle = fit_text(
      answers("is the thresholded WT_heat_up set a restatement of either curated lens?",
              sprintf(paste("no. Whole-set membership is %d of %d in TCR_activation (%.1f%%),",
                            "%d of %d in HSR_core (%.1f%%), none in both, and %.0f%% unassigned.",
                            "This partitions MEMBERSHIP, not rank enrichment."),
                      attr_n("TCR_activation_only"), attr_total, attr_pct("TCR_activation_only"),
                      attr_n("HSR_core_only"), attr_total, attr_pct("HSR_core_only"),
                      attr_pct("neither"))),
      width_in = W_ATTR, tier = "subtitle"),
    x = NULL,
    y = "fraction of WT_heat_up genes",
    caption = fit_text(
      "Membership is counted over the whole arm, never over a leading edge. ",
      "The unassigned remainder is reported as a remainder: it is not named, and it is ",
      "evidence for no mechanism. Claim tier: a direct count over frozen sets. ",
      HONEST_CEILING,
      width_in = W_ATTR, tier = "caption")
  ) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        plot.caption.position = "plot")

save_figure(fig_attr, STAGE, "wtheatup_attribution", overview = TRUE,
            width = W_ATTR, height = 7.0, config = FIG_CFG)

nes_src <- lens_nes %>%
  dplyr::filter(term %in% TERM_LEVELS, contrast %in% CONTRAST_LEVELS) %>%
  dplyr::mutate(
    term = factor(term, levels = TERM_LEVELS),
    contrast = factor(contrast, levels = CONTRAST_LEVELS),
    direction = dplyr::case_when(nes > 0 ~ "up at 39 °C", nes < 0 ~ "down at 39 °C", TRUE ~ "zero"),
    label = sprintf("%.2f", nes),
    # Effective set size travels with every NES, on the face. Here the two lenses
    # are wholly present in the ranked list, so nominal and effective coincide --
    # which is itself worth printing, because it is what makes these three
    # contrasts comparable to each other and is not true downstream.
    n_nominal = set_size,
    n_effective = set_size,
    size_fdr_label = sprintf("n %d of %d  ·  FDR %s",
                             n_effective, n_nominal,
                             ifelse(padj < 1e-3, formatC(padj, format = "e", digits = 1),
                                    formatC(padj, format = "f", digits = 3)))
  ) %>%
  dplyr::arrange(contrast, term)
readr::write_csv(round_numeric_cols(nes_src, sig = 9),
                 file.path(TBL_OVW, "lens_nes_by_contrast.csv"))

term_fills <- c(
  HSR_core = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  TCR_activation = FIG_CFG$colors$okabe_ito$orange %||% "#E69F00"
)

# Horizontal, so the effective set size and FDR can sit in a clean right-hand
# column the way the JIA purge panel prints them, instead of being squeezed above
# a vertical bar.
W_NES <- 12.0
NES_MAX <- max(nes_src$nes)
X_ANNOT <- NES_MAX * 1.22
X_LIM <- NES_MAX * 1.95
hsr_wt <- dplyr::filter(nes_src, term == "HSR_core", contrast == "WT_heat")
hsr_ko <- dplyr::filter(nes_src, term == "HSR_core", contrast == "KO_heat")

fig_nes <- ggplot(nes_src, aes(y = contrast, x = nes, fill = term)) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey35") +
  geom_col(position = position_dodge(width = 0.74), width = 0.64) +
  geom_text(aes(label = label), position = position_dodge(width = 0.74),
            hjust = -0.25, size = (FIG_CFG$figures$label_size %||% 4) * 1.0,
            colour = "grey15") +
  geom_text(aes(x = X_ANNOT, label = size_fdr_label),
            position = position_dodge(width = 0.74), hjust = 0,
            size = (FIG_CFG$figures$label_size %||% 4) * 0.9, colour = "grey25") +
  scale_fill_manual(values = term_fills,
                    labels = c(HSR_core = "HSR core (curated)",
                               TCR_activation = "TCR/IEG activation (curated)"),
                    name = "Lens") +
  scale_y_discrete(limits = rev(CONTRAST_LEVELS),
                   labels = function(x) contrast_label(x, short = TRUE)) +
  scale_x_continuous(limits = c(0, X_LIM), breaks = scales::pretty_breaks(4),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title = fit_text("Warming to 39 °C induces the curated heat-shock core, with or without cGAS",
                     width_in = W_NES, tier = "title"),
    subtitle = fit_text(
      answers("if almost no HSR_core gene is in the set, is the heat-shock response simply absent from the 39 °C response?",
              sprintf(paste("no. HSR_core is the strongest of the two lenses in all three contrasts",
                            "(WT NES %+.2f, FDR %.1e; cGAS-KO NES %+.2f, FDR %.1e)."),
                      hsr_wt$nes[1], hsr_wt$padj[1], hsr_ko$nes[1], hsr_ko$padj[1])),
      width_in = W_NES, tier = "subtitle"),
    x = "GSEA NES (signed t ranking)",
    y = NULL,
    caption = fit_text(
      "NES > 0 means the lens is enriched toward genes higher at 39 °C. Effective set size is ",
      "printed as n present of n nominal; both lenses are wholly recovered in these rankings. ",
      "Enrichment of a lens is a separate measurement from that lens being contained in the ",
      "thresholded set. Claim tier: confirmatory for this temperature contrast; correlative for ",
      "any fever reading. ", HONEST_CEILING,
      width_in = W_NES, tier = "caption")
  ) +
  project_theme(config = FIG_CFG) +
  theme(legend.position = "bottom",
        panel.grid.major.y = element_blank(),
        plot.caption.position = "plot")

save_figure(fig_nes, STAGE, "lens_nes_by_contrast", overview = TRUE,
            width = W_NES, height = 6.5, config = FIG_CFG)

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
W_RANK <- 12.0

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
            inherit.aes = FALSE, hjust = 1, vjust = -0.15,
            size = (FIG_CFG$figures$label_size %||% 4) * 0.95,
            fontface = "bold", colour = "grey15") +
  annotate("text", x = 13.0, y = 3.42,
           label = sprintf("export gate: %d of the first %s ranked genes passed",
                           gate_meta$gate_pass_n[1],
                           format(gate_meta$gate_inner_total_n[1], big.mark = ",")),
           hjust = 0, size = (FIG_CFG$figures$label_size %||% 4) * 0.95, colour = "grey20") +
  annotate("text", x = gate_meta$gate_max_rank_pct[1] + 1.4, y = 0.68,
           label = sprintf("deepest gene admitted: %s\nrank %s; t %.2f; logFC %.2f",
                           gate_meta$gate_deepest_gene[1],
                           format(gate_meta$gate_inner_total_n[1], big.mark = ","),
                           gate_meta$gate_deepest_t[1], gate_meta$gate_deepest_logFC[1]),
           hjust = 0, size = (FIG_CFG$figures$label_size %||% 4) * 0.9, colour = "grey20") +
  scale_colour_manual(values = set_cols, guide = "none") +
  scale_fill_manual(values = set_cols, guide = "none") +
  scale_x_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 75, 100),
                     labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0))) +
  labs(
    title = fit_text("The heat-shock core moved moderately; the export gate kept only extreme movers",
                     width_in = W_RANK, tier = "title"),
    subtitle = fit_text(
      answers("the heat-shock response is induced, so why is almost none of it in the set?",
              sprintf(paste("because the gate is a threshold on effect size, not on pathway.",
                            "%s genes sit at median rank percentile %.1f%% (median t %.2f) while %s sits at",
                            "%.1f%% (median t %.2f), and the deepest gene the gate admitted sits at %.2f%%."),
                      rank_summary$set[rank_summary$set == "WT_heat_up"][1],
                      rank_summary$median_rank_pct[rank_summary$set == "WT_heat_up"][1],
                      rank_summary$median_t[rank_summary$set == "WT_heat_up"][1],
                      rank_summary$set[rank_summary$set == "HSR_core"][1],
                      rank_summary$median_rank_pct[rank_summary$set == "HSR_core"][1],
                      rank_summary$median_t[rank_summary$set == "HSR_core"][1],
                      gate_meta$gate_max_rank_pct[1])),
      width_in = W_RANK, tier = "subtitle"),
    x = "rank percentile in WT_heat signed-t ranking",
    y = NULL,
    caption = fit_text(
      sprintf(paste("Each point is one gene; ridgelines show density along the same axis. The shaded band",
                    "is the empirical gate span, %.4f%% to %.2f%% (rank 1 to %s of %s), and only %d genes",
                    "inside it passed, because the gate also required logFC >= 1. The WT_heat_up row is",
                    "left-shifted by construction and is drawn as the gate's output, not as a result."),
              gate_meta$gate_min_rank_pct[1], gate_meta$gate_max_rank_pct[1],
              format(gate_meta$gate_inner_total_n[1], big.mark = ","),
              format(gate_meta$ranking_total_n[1], big.mark = ","),
              gate_meta$gate_pass_n[1]),
      width_in = W_RANK, tier = "caption")
  ) +
  project_theme(config = FIG_CFG) +
  theme(plot.caption.position = "plot")

save_figure(fig_rank, STAGE, "hsr_rank_position_panel", overview = TRUE,
            width = W_RANK, height = 6.8, config = FIG_CFG)

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
# coord_fixed() letterboxes the Euler and Venn panels and the caption is laid out
# against the PANEL, so the usable caption width is narrower than the emitted
# figure width. plot.caption.position = "plot" recovers the full width; the wrap
# is still set to fit inside the narrowest of the three rather than tracking each.
CAPTION_WRAP <- 100L
wrap_caption <- function(...) paste(strwrap(paste0(...), width = CAPTION_WRAP), collapse = "\n")

# These three panels re-draw one fact that the attribution bar has already
# answered, so each says on its face that it corroborates rather than answers.
# The dense material -- the human counterpart, the conditional-NES check, the
# anchor-solver contract -- moved to the stage README, which is where a caveat
# survives being shrunk into a journal column. Exactly one warning stays on the
# face, and it is the one a reader can misread the geometry without.
membership_ceiling <- paste(
  "Set membership only: a lens enriching is a separate measurement from that lens being",
  "contained in WT_heat_up.")
corroboration_line <- sprintf(
  "Corroborates the membership answer: %d of %d WT_heat_up genes are in neither curated lens.",
  as.integer(membership$n_mouse[membership$component == "WT_heat_up only"][1]),
  as.integer(membership$wt_heat_up_n[1]))

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
  sig = 9), file.path(TBL_OVW, "hsr_lens_membership_euler.csv"))

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
    # coord_fixed + void letterbox the panel, so the subtitle's usable width is
    # narrower than the emitted figure; keep both lines short and pre-broken.
    subtitle = paste(corroboration_line,
                     sprintf(paste("Circle areas are fitted to those counts (eulerr stress %.1g),",
                                   "so geometry is evidence here."),
                             euler_stress),
                     sep = "\n"),
    x = NULL, y = NULL,
    caption = wrap_caption(
      "Bold numbers sit inside the region they count; each callout names an overlap and its leader line ",
      "points at it. Open diamonds mark the two intersections that exist and are EMPTY. ",
      membership_ceiling)
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
                 file.path(TBL_OVW, "hsr_lens_membership_venn.csv"))

fig_venn <- ggvenn::ggvenn(
  venn_df, columns = SET_ORDER, show_percentage = FALSE,
  fill_color = unname(set_cols_membership[SET_ORDER]), fill_alpha = 0.28,
  stroke_color = "grey25", stroke_size = 0.6,
  text_size = (FIG_CFG$figures$label_size %||% 4) * 1.15,
  set_name_size = (FIG_CFG$figures$label_size %||% 4) * 1.0) +
  coord_fixed(clip = "off") +
  labs(
    title = "The same seven counts, drawn where geometry carries no quantity",
    # coord_fixed letterboxes this panel too; keep both subtitle lines short.
    subtitle = paste(corroboration_line,
                     "Fixed layout: a printed 0 marks a region that exists and is empty.",
                     sep = "\n"),
    x = NULL, y = NULL,
    caption = wrap_caption(
      "All three circles are equal and equally overlapping by convention and only the printed counts ",
      "vary, so this layout cannot assert a containment the counts contradict. ",
      membership_ceiling)
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
                 file.path(TBL_OVW, "hsr_lens_membership_upset.csv"))

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
    subtitle = paste(corroboration_line,
                     paste("Bars are intersection sizes over the same three sets; the two zero-height",
                           "bars are labelled rather than omitted."),
                     sep = "\n"),
    caption = wrap_caption(
      "Dots mark which sets each bar belongs to; the left bars are the three set sizes. Every ",
      "combination of the three sets gets a column, so an intersection with no genes appears as a ",
      "zero-height bar labelled 'empty' rather than as a missing column. ",
      membership_ceiling),
    theme = project_theme(config = FIG_CFG) + theme(plot.caption.position = "plot"))

save_figure(fig_euler, STAGE, "hsr_lens_membership_euler", overview = TRUE,
            width = FIG_W_EULER, height = 7.5, void = TRUE, config = FIG_CFG)
save_figure(fig_venn, STAGE, "hsr_lens_membership_venn", overview = TRUE,
            width = FIG_W_VENN, height = 7.5, void = TRUE, config = FIG_CFG)
save_figure(fig_upset, STAGE, "hsr_lens_membership_upset", overview = TRUE,
            width = FIG_W_UPSET, height = 7.0, config = FIG_CFG)

message(sprintf("[19_viz] eulerr fit: stress = %.3g, diagError = %.3g, max |region residual| = %.3g",
                euler_stress, euler_diag_error, max(abs(euler_src$residual))))

# =============================================================================
# The handoff panel: mouse gate 213 -> human projection 199 -> what is testable.
#
# Two numbers travel through this project side by side. This is the last place a
# reader meets the mouse 213, so it is where the two are reconciled, and where the
# reader is told that the 199 is not the same 199 in every human compartment.
# =============================================================================
W_BRIDGE <- 13.0
COMPARTMENT_LABELS <- c(
  human_treg_arthritis = "JIA sorted T cells",
  human_pbmc_febrile = "Kawasaki PBMC",
  human_ra_synovium = "RA synovium",
  sting_positive_control = "SAVI PBMC"
)
side_cols <- c(mouse = FIG_CFG$colors$diverging$up %||% "#B35806",
               human = FIG_CFG$colors$okabe_ito$blue %||% "#0072B2")

funnel <- bridge %>%
  dplyr::filter(.data$block == "funnel") %>%
  dplyr::arrange(.data$step)
lens_rows <- bridge %>% dplyr::filter(.data$block == "lens")
lens_note <- lens_rows %>%
  dplyr::group_by(.data$step) %>%
  dplyr::summarise(txt = paste(sprintf("%s %d",
                                       sub("^of WT_heat_up, in curated ([A-Za-z_]+).*$", "\\1", .data$label),
                                       .data$n_genes), collapse = "  ·  "),
                   .groups = "drop")
funnel <- funnel %>%
  dplyr::left_join(lens_note, by = "step") %>%
  dplyr::mutate(
    label = factor(.data$label, levels = rev(.data$label)),
    point_label = ifelse(is.na(.data$txt),
                         format(.data$n_genes, big.mark = ","),
                         sprintf("%s\nin a curated lens:  %s",
                                 format(.data$n_genes, big.mark = ","), .data$txt)))

X_FLOOR <- 8
fig_funnel <- ggplot(funnel, aes(y = .data$label, x = .data$n_genes, colour = .data$side)) +
  geom_segment(aes(x = X_FLOOR, xend = .data$n_genes, yend = .data$label), linewidth = 1.6) +
  geom_point(size = (FIG_CFG$figures$point_size %||% 2.4) * 1.9) +
  geom_text(aes(label = .data$point_label), hjust = -0.14, vjust = 0.42,
            lineheight = 0.95, size = (FIG_CFG$figures$label_size %||% 4) * 0.95,
            colour = "grey15") +
  scale_colour_manual(values = side_cols,
                      labels = c(mouse = "mouse gene symbols", human = "human gene symbols"),
                      name = NULL) +
  scale_x_log10(limits = c(X_FLOOR, 1.1e6), breaks = c(10, 100, 1000, 10000),
                labels = function(x) format(x, big.mark = ",", scientific = FALSE),
                expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(labels = function(x) sub(" \\(", "\n(", x)) +
  labs(title = "What the mouse gate hands over",
       x = "genes (log scale)", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(legend.position = "bottom", panel.grid.major.y = element_blank())

down_min <- bridge %>% dplyr::filter(.data$block == "downstream_min")
down_max <- bridge %>% dplyr::filter(.data$block == "downstream_max")
n_human <- funnel$denominator[funnel$step == 5][1]
n_everywhere <- funnel$n_genes[funnel$step == 5][1]
n_lists_total <- funnel$n_lists[funnel$step == 5][1]

if (nrow(down_min) > 0L) {
  down <- down_min %>%
    dplyr::transmute(compartment = .data$label, n_lists = .data$n_lists,
                     list_len_min = .data$list_len_min, list_len_max = .data$list_len_max,
                     n_eff_min = .data$n_genes) %>%
    dplyr::left_join(down_max %>% dplyr::transmute(compartment = .data$label,
                                                   n_eff_max = .data$n_genes),
                     by = "compartment") %>%
    dplyr::mutate(
      display = unname(COMPARTMENT_LABELS[.data$compartment]),
      display = ifelse(is.na(.data$display), .data$compartment, .data$display),
      annot = sprintf("%d of %d to %d of %d   ·   %d ranked lists, %s to %s genes long",
                      .data$n_eff_min, n_human, .data$n_eff_max, n_human, .data$n_lists,
                      format(.data$list_len_min, big.mark = ","),
                      format(.data$list_len_max, big.mark = ","))) %>%
    dplyr::arrange(.data$n_eff_max) %>%
    dplyr::mutate(display = factor(.data$display, levels = .data$display))

  fig_down <- ggplot(down, aes(y = .data$display)) +
    geom_vline(xintercept = n_everywhere, linetype = "dashed",
               linewidth = 0.5, colour = "grey35") +
    geom_segment(aes(x = .data$n_eff_min, xend = .data$n_eff_max, yend = .data$display),
                 linewidth = 1.6, colour = side_cols[["human"]]) +
    geom_point(aes(x = .data$n_eff_min), size = (FIG_CFG$figures$point_size %||% 2.4) * 1.5,
               colour = side_cols[["human"]]) +
    geom_point(aes(x = .data$n_eff_max), size = (FIG_CFG$figures$point_size %||% 2.4) * 1.5,
               colour = side_cols[["human"]]) +
    geom_text(aes(x = 2, label = .data$annot), hjust = 0, vjust = -1.15,
              size = (FIG_CFG$figures$label_size %||% 4) * 0.9, colour = "grey20") +
    annotate("text", x = n_everywhere + 3, y = 0.4,
             label = sprintf("dashed rule: %d genes are present in all %d lists",
                             n_everywhere, n_lists_total),
             hjust = 0, size = (FIG_CFG$figures$label_size %||% 4) * 0.9, colour = "grey25") +
    scale_x_continuous(limits = c(0, n_human), breaks = scales::pretty_breaks(5),
                       expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_discrete(expand = expansion(add = c(1.0, 0.75))) +
    labs(title = sprintf("How much of the %d each compartment can test", n_human),
         x = sprintf("genes of the %d present in a ranked list of that compartment", n_human),
         y = NULL) +
    project_theme(config = FIG_CFG) +
    theme(panel.grid.major.y = element_blank())
} else {
  fig_down <- ggplot() +
    annotate("text", x = 0, y = 0, hjust = 0.5, size = (FIG_CFG$figures$label_size %||% 4) * 1.1,
             colour = "grey35",
             label = paste("No human ranked list was reachable from this checkout,",
                           "so downstream recovery is\nnot drawn. It is a real question",
                           "with no answer available here, not an absence of one.")) +
    labs(title = "How much of the projected set each compartment can test", x = NULL, y = NULL) +
    project_theme(config = FIG_CFG) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank())
}

bridge_headline <- sprintf(
  paste("the mouse gate output is %d genes and the human set is %d, and they are the same object",
        "seen either side of ortholog projection. Only %d of the %d are present in all %d human",
        "ranked lists, so the %d is not the same %d in every compartment."),
  funnel$n_genes[funnel$step == 3][1], n_human, n_everywhere, n_human,
  n_lists_total, n_human, n_human)

# Stacked, not side by side: each half carries long right-hand annotations, and
# two half-width panels clip them. One column gives both the full emitted width.
fig_bridge <- patchwork::wrap_plots(fig_funnel, fig_down, ncol = 1, heights = c(1, 0.82)) +
  patchwork::plot_annotation(
    title = fit_text("Mouse gate 213, human projection 199, and 17 genes every human list can test",
                     width_in = W_BRIDGE, tier = "title"),
    subtitle = fit_text(
      answers("what does the human side actually receive from the mouse gate?", bridge_headline),
      width_in = W_BRIDGE, tier = "subtitle"),
    caption = fit_text(
      "Top: each step is a count, drawn on a log axis because the first and last differ by three ",
      "orders of magnitude. The curated-lens counts are membership against the same two lenses in ",
      "each species, and the two lenses behave differently under ortholog projection: ",
      "TCR_activation is a hand-curated panel with a strictly 1:1 human-to-mouse map, so it is 66 ",
      "genes on both sides and its 12 are literally the same 12 genes; HSR_core is assembled ",
      "paralog-complete from Reactome and GO with a 1:many map, so it is 47 mouse against 56 human ",
      "and its mouse Hspa1a/Hspa1b pair collapses to one human symbol in the projected set. Read ",
      "the mouse and human lens counts as counterparts, not as one measurement. Bottom: the bar spans the ",
      "fewest to the most of the projected genes present in any one ranked list of that ",
      "compartment, and the dashed rule marks the genes present in every list. Recovery is a ",
      "property of the analysis that produced each ranked list, not of the biology, so a magnitude ",
      "comparison across compartments is a different measurement on a smaller set. Claim tier: ",
      "direct counts over frozen sets and published ranked lists.",
      width_in = W_BRIDGE, tier = "caption"),
    theme = project_theme(config = FIG_CFG) + theme(plot.caption.position = "plot"))

save_figure(fig_bridge, STAGE, "gate_projection_bridge", overview = TRUE,
            width = W_BRIDGE, height = 10.5, config = FIG_CFG)

expected <- file.path(FIG_OVW, c(
  "wtheatup_attribution.pdf", "wtheatup_attribution.png",
  "lens_nes_by_contrast.pdf", "lens_nes_by_contrast.png",
  "hsr_rank_position_panel.pdf", "hsr_rank_position_panel.png",
  "hsr_lens_membership_euler.pdf", "hsr_lens_membership_euler.png",
  "hsr_lens_membership_venn.pdf", "hsr_lens_membership_venn.png",
  "hsr_lens_membership_upset.pdf", "hsr_lens_membership_upset.png",
  "gate_projection_bridge.pdf", "gate_projection_bridge.png"
))
missing <- expected[!file.exists(expected)]
if (length(missing) > 0L)
  stop("[19_viz] Missing expected figure output(s): ", paste(missing, collapse = ", "))
stray <- list.files(FIG_DIR, pattern = "\\.(png|pdf)$", full.names = TRUE)
if (length(stray) > 0L)
  stop("[19_viz] Figures must live under figures/_overview/; found: ", paste(stray, collapse = ", "))

message("[19_viz] COMPLETE: wrote ", length(expected) / 2,
        " _overview figures and their same-stem source tables for stage ", STAGE)
