#!/usr/bin/env Rscript
# 21_semantic_decomposition_viz.R -- VIZ ONLY (no computing)
# =============================================================================
# Draws the overview panel for stage 13_semantic_decomp from the frozen object
# written by 20_semantic_decomposition.R. Every number on the face is read from
# that object.
#
#   wtheatup_lens_proximity -- graded per-gene proximity to each curated lens,
#     against a background matched on annotation depth
#
# Run from project root:
#   Rscript 02_analysis/scripts/21_semantic_decomposition_viz.R
# =============================================================================

source("02_analysis/config/config.R")
source("02_analysis/helpers/figure_style.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
options(stringsAsFactors = FALSE)

STAGE  <- "13_semantic_decomp"
SCRIPT <- "02_analysis/scripts/21_semantic_decomposition_viz.R"
set_paper_style(config = FIG_CFG)

obj <- readRDS(file.path(DIR_OBJECTS, "20_semantic_decomp.rds"))
coverage <- obj$coverage
prox     <- obj$proximity
lens_sum <- obj$lens_summary
prov     <- setNames(obj$provenance$value, obj$provenance$key)

# Lenses the frozen object scores but this panel does not draw. Lombardi2022_HIF is
# scored against 72 of the set's genes where the other two are scored against 153,
# because its own depth band saturates for the rest; two lenses answering over one
# denominator and a third over a different one is not a comparison a reader can make
# on one axis. It stays in the stage's per-gene and coverage-loss tables, where the
# denominator travels with it.
LENS_NOT_DRAWN <- "Lombardi2022_HIF"
prox     <- prox[!prox$lens %in% LENS_NOT_DRAWN, , drop = FALSE]
lens_sum <- lens_sum[!lens_sum$lens %in% LENS_NOT_DRAWN, , drop = FALSE]

OI <- FIG_CFG$colors$okabe_ito
LAB <- (FIG_CFG$figures$label_size %||% 4)
PT  <- (FIG_CFG$figures$point_size %||% 2.4)

n_annot   <- coverage$n_genes[coverage$step == 2L]

## ---------------------------------------------------------------------------
## Figure 2 -- graded proximity to each lens
## ---------------------------------------------------------------------------
lens_levels <- lens_sum$lens[order(-lens_sum$n_members_in_query)]
prox2 <- prox %>%
  dplyr::filter(!is.na(.data$pct_depth_matched_background)) %>%
  dplyr::mutate(lens = factor(.data$lens, levels = lens_levels),
                status = dplyr::case_when(
                  !.data$counted        ~ "too sparsely annotated to score",
                  .data$is_lens_member  ~ "in the lens",
                  TRUE                  ~ "not in the lens"))
# The counts go in the axis labels, over in-panel text: three blocks centred
# on three categories collide as soon as any of them is wider than the category slot.
axis_lab <- setNames(
  sprintf("%s\n%d of the set are in it\n%d of %d scored sit above the 95th\npercentile, against %.1f by chance\n%d too sparsely annotated to score",
          lens_sum$lens, lens_sum$n_members_in_query, lens_sum$n_above_p95,
          lens_sum$n_scored, lens_sum$n_expected_by_chance,
          lens_sum$n_excluded_saturated_band),
  lens_sum$lens)

pD <- ggplot(prox2, aes(x = .data$lens, y = .data$pct_depth_matched_background)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", linewidth = 0.55,
             colour = "grey35") +
  geom_point(aes(colour = .data$status, size = .data$status,
                 alpha = .data$status, shape = .data$status),
             position = position_jitter(width = 0.22, height = 0, seed = 7L)) +
  scale_x_discrete(labels = axis_lab) +
  scale_colour_manual(values = c("in the lens" = OI$vermillion,
                                 "not in the lens" = "grey55",
                                 "too sparsely annotated to score" = OI$sky_blue),
                      name = NULL) +
  scale_shape_manual(values = c("in the lens" = 16, "not in the lens" = 16,
                                "too sparsely annotated to score" = 4), name = NULL) +
  scale_size_manual(values = c("in the lens" = PT * 1.9,
                               "not in the lens" = PT * 0.75,
                               "too sparsely annotated to score" = PT * 0.75),
                    guide = "none") +
  scale_alpha_manual(values = c("in the lens" = 1, "not in the lens" = 0.55,
                                "too sparsely annotated to score" = 0.5),
                     guide = "none") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 0.95),
                     expand = expansion(add = c(0.04, 0.04))) +
  labs(title = "For the genes no lens contains, how close does each one sit?",
       subtitle = paste0("One dot per gene, at its annotation's similarity to the nearest other member of that lens. Red dots are the lens's own members and\n",
                         "act as the control; crosses are genes too sparsely annotated to rank."),
       x = NULL, y = "percentile against a\ndepth-matched background") +
  project_theme(config = FIG_CFG) +
  theme(legend.position = "bottom", panel.grid.major.x = element_blank())

conf <- prox %>% dplyr::mutate(lens = factor(.data$lens, levels = lens_levels))
rho_txt <- lens_sum %>%
  dplyr::mutate(lens = factor(.data$lens, levels = lens_levels),
                txt = sprintf("Spearman %+.2f", .data$rho_max_vs_n_terms))
pE <- ggplot(conf, aes(x = .data$n_terms, y = .data$max_to_lens)) +
  geom_point(colour = "grey55", size = PT * 0.7, alpha = 0.6) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 0.9,
              colour = OI$blue) +
  geom_text(data = rho_txt, aes(x = Inf, y = -Inf, label = .data$txt),
            hjust = 1.08, vjust = -0.7, size = LAB * 0.95, colour = OI$blue) +
  facet_wrap(~ lens, nrow = 1) +
  scale_x_log10() +
  labs(title = "Why the percentiles above are taken within a depth band",
       subtitle = paste0("A gene annotated to more terms sits closer to every reference set. That is how much has been recorded about the gene, not\n",
                         "biology, and pooling the background over all depths would read it as signal."),
       x = sprintf("GO %s terms annotated to the gene (log scale)", prov[["ontology"]]),
       y = "similarity to the\nnearest lens member") +
  project_theme(config = FIG_CFG)

fig2 <- (pD / pE) +
  patchwork::plot_layout(heights = c(1.35, 1)) +
  patchwork::plot_annotation(
    title = "Membership is binary, and the genes it leaves out are not all equally far away",
    subtitle = sprintf(paste0("The two lenses contain %d of the %d annotated genes. This asks the graded question of the rest: how close does each gene's\n",
                              "GO %s annotation sit to each lens? Scoring high means a gene resembles the lens in what has been recorded about it, which is\n",
                              "not evidence it participates in the process the lens is named for."),
                      sum(lens_sum$n_members_in_query), n_annot, prov[["ontology"]]),
    theme = project_theme(config = FIG_CFG))

save_overview(
  fig2, STAGE, "wtheatup_lens_proximity",
  table = lens_sum,
  finding = paste0(sprintf("Of the %d annotated genes in the mouse 39 °C-derived up arm, the two curated lenses drawn here contain %d. ",
                           n_annot, sum(lens_sum$n_members_in_query)),
                   paste(sprintf("%s holds %d, and of the %d genes scorable against a depth-matched background %d sit above its 95th percentile against a chance count of %.1f",
                                 lens_sum$lens, lens_sum$n_members_in_query,
                                 lens_sum$n_scored, lens_sum$n_above_p95,
                                 lens_sum$n_expected_by_chance),
                         collapse = "; "), "."),
  script = SCRIPT, fn = "top-level (pD/pE)",
  config_kv = sprintf("semantic.depth_bands = %s, semantic.bg_sample_size = %s, semantic.measure = %s, lens_not_drawn = %s",
                      prov[["depth_bands"]], prov[["bg_sample"]], prov[["measure"]],
                      LENS_NOT_DRAWN),
  input = "03_results/objects/20_semantic_decomp.rds",
  how_to_read = paste0(
    "Upper panel: one dot per gene per lens, at the percentile of its similarity to the nearest other lens member. ",
    "No gene is scored against itself, so the red dots control that the measure works. The dashed rule is the 95th ",
    "percentile; each axis label gives the observed count above it beside the count expected by chance, and a lens ",
    "whose two counts match has no graded signal beyond its listed members. Crosses are genes whose depth band ",
    "already saturates: drawn, excluded from the counts, tallied on the axis. Lower panel: why the banding is ",
    "needed. A third scored lens, ", LENS_NOT_DRAWN, ", is left off over a different denominator and stays in the ",
    "stage's per-gene and coverage-loss tables."),
  config = FIG_CFG, wide = TRUE, height = 12.5)

message("21_semantic_decomposition_viz: VIZ DONE.")
