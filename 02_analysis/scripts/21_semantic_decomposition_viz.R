#!/usr/bin/env Rscript
# 21_semantic_decomposition_viz.R -- VIZ ONLY (no computing)
# =============================================================================
# Draws the two overview panels for stage 13_semantic_decomp from the frozen
# object written by 20_semantic_decomposition.R. Every number on every face is
# read from that object.
#
#   wtheatup_semantic_coherence -- can the ontology see the set, is the set one
#     process, and can it be split into modules
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
coh      <- obj$coherence
prox     <- obj$proximity
lens_sum <- obj$lens_summary
scale_v  <- setNames(obj$scale$value, obj$scale$quantity)
lens_lost <- obj$lens_lost
prov     <- setNames(obj$provenance$value, obj$provenance$key)

OI <- FIG_CFG$colors$okabe_ito
LAB <- (FIG_CFG$figures$label_size %||% 4)
PT  <- (FIG_CFG$figures$point_size %||% 2.4)
MAX_K_VIZ <- as.integer(YAML_CONFIG$semantic$max_k_partition %||% 15L)

FAM_COL <- c(query              = OI$vermillion,
             single_go_term     = OI$bluish_green,
             null_depth_matched = OI$blue,
             null_uniform       = OI$sky_blue,
             curated_lens       = OI$reddish_purple,
             gate_sibling       = OI$orange)
FAM_LAB <- c(query              = "WT_heat_up (the set under test)",
             single_go_term     = "one GO term and its offspring",
             null_depth_matched = "null matched on annotation depth",
             null_uniform       = "null drawn uniformly",
             curated_lens       = "curated lens",
             gate_sibling       = "sibling gate set")

n_nominal <- coverage$n_genes[coverage$step == 1L]
n_annot   <- coverage$n_genes[coverage$step == 2L]

## ---------------------------------------------------------------------------
## Figure 1, panel A -- what the ontology can and cannot see
## ---------------------------------------------------------------------------
short <- coverage %>%
  dplyr::filter(.data$step >= 3L, .data$step <= 5L, .data$n_genes > 0L) %>%
  dplyr::mutate(short_label = sub("^no BP: ", "", .data$label))
recov <- coverage[coverage$step == 6L, ]
lad <- data.frame(
  label = c(sprintf("in the set\n%d genes", n_nominal),
            sprintf("GO %s can place it\n%d genes", prov[["ontology"]], n_annot)),
  n = c(n_nominal, n_annot),
  kind = c("nominal", "placed"))
lad$label <- factor(lad$label, levels = rev(lad$label))

lost_note <- lens_lost %>%
  dplyr::filter(nzchar(.data$lost_genes)) %>%
  dplyr::mutate(txt = sprintf("%s (%s of %d in the set)", .data$lost_genes,
                              .data$lens, .data$n_overlap_nominal))
shortfall_txt <- paste(sprintf("%d %s", short$n_genes, short$short_label),
                       collapse = "\n")

pA <- ggplot(lad, aes(x = .data$n, y = .data$label, fill = .data$kind)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = format(.data$n, big.mark = ",")), hjust = -0.25,
            size = LAB * 1.05, fontface = "bold", colour = "grey15") +
  annotate("text", x = 4, y = 0.62, hjust = 0, vjust = 1, colour = "grey25",
           size = LAB * 0.9,
           label = sprintf("the %d it cannot place:\n%s",
                           n_nominal - n_annot, shortfall_txt)) +
  scale_fill_manual(values = c(nominal = "grey72", placed = OI$bluish_green),
                    guide = "none") +
  scale_x_continuous(limits = c(0, n_nominal * 1.24),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(add = c(1.9, 0.6))) +
  labs(title = sprintf("Can GO %s see the set?", prov[["ontology"]]),
       subtitle = paste0(
         sprintf("Curated evidence only, with %s annotation set aside%s. ",
                 prov[["drop_evidence"]],
                 if (nrow(recov) && recov$n_genes[1] > 0)
                   sprintf(" (keeping it would recover %d)", recov$n_genes[1]) else ""),
         if (nrow(lost_note))
           sprintf("The gap is not neutral:\n%s carries no %s term in mouse or in human, so a lens overlap drawn\nelsewhere contains a gene this panel has to drop.",
                   lost_note$lost_genes[1], prov[["ontology"]])
         else "Every curated-lens member in the set carries an annotation."),
       x = "genes", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank())

## ---------------------------------------------------------------------------
## Figure 1, panel B -- is the set one process?
## ---------------------------------------------------------------------------
theme_group <- function(lbl) sub("^(GO:[0-9]+ [^\\[]+)\\[.*$", "\\1", lbl)
cohB <- coh %>%
  dplyr::mutate(grp = dplyr::case_when(
    .data$family == "single_go_term" ~ trimws(theme_group(.data$label)),
    .data$family == "curated_lens"   ~ .data$label,
    .data$family == "gate_sibling"   ~ .data$label,
    .data$family == "query"          ~ "WT_heat_up",
    .data$family == "null_uniform"   ~ "null, drawn uniformly",
    .data$family == "null_depth_matched" ~ "null, matched on annotation depth"))
ord <- cohB %>% dplyr::group_by(.data$grp, .data$family) %>%
  dplyr::summarise(m = mean(.data$mean_sim), n = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(.data$m)
cohB$grp <- factor(cohB$grp, levels = ord$grp)
q_mean <- scale_v[["query_mean_sim"]]

pB <- ggplot(cohB, aes(x = .data$mean_sim, y = .data$grp, colour = .data$family)) +
  geom_vline(xintercept = q_mean, linetype = "dashed", linewidth = 0.55,
             colour = OI$vermillion) +
  geom_point(size = PT * 1.25, alpha = 0.85,
             position = position_jitter(height = 0.13, width = 0, seed = 1L)) +
  # orientation = "y": summarise x within each y group, not the reverse
  stat_summary(fun = mean, geom = "point", shape = 124, size = PT * 3.2,
               colour = "grey20", orientation = "y") +
  # One key serves both panels and it lives under the last one. patchwork cannot merge
  # these two guides because the panels carry different layer sets, so collecting them
  # emits the legend twice instead of once.
  scale_colour_manual(values = FAM_COL, labels = FAM_LAB, name = NULL,
                      breaks = names(FAM_LAB)[names(FAM_LAB) %in% cohB$family],
                      guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.06))) +
  labs(title = "Is it one process, or many?",
       subtitle = sprintf(paste0("Mean pairwise Wang similarity, every set matched to the query's %d genes. ",
                                 "Against a null\nmatched on annotation depth the set sits %+.2f SD (p %s). ",
                                 "Against a uniform null it sits %+.2f SD,\nand that difference is the confound, not the biology: the gate admits genes carrying a median\n%d %s terms where the expressed background carries %d."),
                          n_annot, scale_v[["z_vs_matched"]],
                          format.pval(scale_v[["p_vs_matched"]], digits = 2, eps = 1e-3),
                          scale_v[["z_vs_uniform"]],
                          round(scale_v[["query_med_terms_per_gene"]]),
                          prov[["ontology"]],
                          round(scale_v[["uniform_null_med_terms_per_gene"]])),
       x = "mean pairwise semantic similarity within the set", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(legend.position = "none", panel.grid.major.y = element_blank())

## ---------------------------------------------------------------------------
## Figure 1, panel C -- can it be split into modules?
## ---------------------------------------------------------------------------
qrow <- coh %>% dplyr::filter(.data$family == "query")
# Read the family means off the data rather than describing the panel from memory:
# the depth-matched null is as unsplittable as the query, and only the uniform null
# separates. Saying otherwise would contradict the points drawn right beside the text.
fam_big <- tapply(coh$largest_cluster_frac, coh$family, mean)
pct <- function(f) sprintf("%.0f%%", 100 * fam_big[[f]])
pC <- ggplot(coh, aes(x = .data$best_sil, y = .data$largest_cluster_frac,
                      colour = .data$family)) +
  geom_hline(yintercept = 0.6, linetype = "dotted", linewidth = 0.5,
             colour = "grey45") +
  geom_point(size = PT * 1.2, alpha = 0.8) +
  geom_point(data = qrow, size = PT * 2.6, shape = 21, stroke = 1.5,
             fill = OI$vermillion, colour = "grey15") +
  # Repelled rather than nudged by a fixed offset. The query's coordinates move with
  # the measure, and a hand-tuned hjust/vjust that clears the neighbouring points at
  # one value drops the label on top of them at the next.
  ggrepel::geom_text_repel(
    data = qrow, aes(label = "WT_heat_up"), size = LAB * 1.0, fontface = "bold",
    colour = "grey15", show.legend = FALSE, seed = 3L,
    min.segment.length = 0, segment.colour = "grey35", segment.size = 0.4,
    box.padding = 0.7, point.padding = 0.5, direction = "both",
    nudge_x = 0.006, nudge_y = -0.09) +
  scale_colour_manual(values = FAM_COL, labels = FAM_LAB, name = NULL,
                      breaks = names(FAM_LAB)[names(FAM_LAB) %in% coh$family]) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.02), expand = expansion(mult = c(0.02, 0.04))) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.14))) +
  labs(title = "Can it be split into modules?",
       subtitle = sprintf(paste0("Best of three linkages plus k-medoids, k = 2 to %d, by average silhouette. No gene-level partition is reported for this set:\n",
                                 "its best split leaves %.0f%% of the genes in one group. A null matched on annotation depth is no better (%s), so being\n",
                                 "unsplittable is a property of well-annotated gene sets and not a finding about this one. The uniformly drawn null does\n",
                                 "separate (%s), because genes carrying one or two terms are semantic islands rather than because those sets have structure.\n",
                                 "Single-GO-term references, coherent by construction, also sit high (%s). Forcing k on any of these would draw an artefact."),
                          MAX_K_VIZ, 100 * qrow$largest_cluster_frac[1],
                          pct("null_depth_matched"), pct("null_uniform"),
                          pct("single_go_term")),
       x = "best average silhouette achieved (higher = more separable)",
       y = "genes in the\nlargest cluster") +
  project_theme(config = FIG_CFG) +
  theme(legend.position = "bottom")

fig1 <- (pA / pB / pC) +
  patchwork::plot_layout(heights = c(0.72, 1.15, 1.0)) +
  patchwork::plot_annotation(
    title = "What the mouse 39 °C-derived up arm is made of, asked of the ontology",
    subtitle = sprintf(paste0("GO %s semantic similarity over the %d gene set, read instead of the expression (%s measure, %s combiner, GOSemSim %s, GO.db %s).\n",
                              "Nothing here reads expression, so no grouping shown can have been derived from an enrichment result. Semantic proximity means two\n",
                              "genes carry similar annotation. It is not evidence of co-regulation, and an enriching set is not evidence that the process it is named for is present."),
                       prov[["ontology"]], n_nominal, prov[["measure"]],
                       prov[["combine"]], prov[["GOSemSim"]], prov[["GO.db"]]),
    theme = project_theme(config = FIG_CFG))

save_overview(
  fig1, STAGE, "wtheatup_semantic_coherence",
  table = coh,
  finding = sprintf(paste0("GO %s places %d of the %d genes in the mouse 39 °C-derived up arm. ",
                           "Its mean pairwise semantic similarity is %.3f, which is %+.2f SD from a null matched on annotation depth (p %s) ",
                           "and so not distinguishable from a random set of equally well-studied genes, ",
                           "while single-GO-term references reach %.3f and the uniform null sits at %.3f. ",
                           "No gene-level partition is reported: the best available split leaves %.0f%% of the genes in one group, and the depth-matched null behaves the same way (%s)."),
                    prov[["ontology"]], n_annot, n_nominal,
                    scale_v[["query_mean_sim"]], scale_v[["z_vs_matched"]],
                    format.pval(scale_v[["p_vs_matched"]], digits = 2, eps = 1e-3),
                    mean(coh$mean_sim[coh$family == "single_go_term"]),
                    scale_v[["uniform_null_mean"]],
                    100 * qrow$largest_cluster_frac[1], pct("null_depth_matched")),
  script = SCRIPT, fn = "top-level (pA/pB/pC)",
  config_kv = sprintf("semantic.measure = %s, semantic.combine = %s, semantic.ontology = %s, semantic.n_null_matched = %s, semantic.seed = %s",
                      prov[["measure"]], prov[["combine"]], prov[["ontology"]],
                      prov[["n_null_matched"]], prov[["seed"]]),
  input = "03_results/objects/20_semantic_decomp.rds",
  how_to_read = paste0(
    "Top, coverage: grey is the set as the gate produced it, green how much carries a GO Biological Process ",
    "annotation, and the inset names what it cannot place. Middle, coherence: each dot is one set's ",
    "mean pairwise similarity, the tick a family mean, the dashed rule the set under test. Both nulls are drawn on ",
    "purpose, and only the depth-matched one carries the comparison, because the uniform one is confounded by the ",
    "gate admitting well-studied genes. Bottom, separability: best silhouette over three linkages, k-medoids and k ",
    "to 15, against the share of genes that split puts in one cluster. A set that decomposes lands right and low. ",
    "This one lands top-left, so no partition is drawn."),
  config = FIG_CFG, wide = TRUE, height = 16.5)

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
# The counts go in the axis labels rather than as in-panel text: three blocks centred
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
       subtitle = paste0("Each grey dot is one gene of the set, placed by how its annotation compares with the nearest OTHER member of that lens, scored against\n",
                         "background genes carrying a similar number of terms. Red dots are genes the lens already contains, and they act as the control: a\n",
                         "working measure puts them high without being shown the answer. Crosses are genes too sparsely annotated to rank within their band."),
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
       subtitle = paste0("A gene annotated to more terms sits closer to every reference set, whatever the set is. That is a property of how much has\n",
                         "been recorded about the gene, not of the biology, and pooling the background over all depths would read it as signal."),
       x = sprintf("GO %s terms annotated to the gene (log scale)", prov[["ontology"]]),
       y = "similarity to the\nnearest lens member") +
  project_theme(config = FIG_CFG)

fig2 <- (pD / pE) +
  patchwork::plot_layout(heights = c(1.35, 1)) +
  patchwork::plot_annotation(
    title = "Membership is binary, and the genes it leaves out are not all equally far away",
    subtitle = sprintf(paste0("The curated lenses contain %d of the %d genes between them. ",
                              "This asks the graded question of the rest: how close does each gene's GO %s\nannotation sit to each lens, ",
                              "measured against background genes annotated to a similar number of terms. ",
                              "A gene scoring high here resembles\nthe lens in what has been recorded about it. ",
                              "That is not evidence it participates in the process the lens is named for."),
                      sum(lens_sum$n_members_in_query), n_annot, prov[["ontology"]]),
    theme = project_theme(config = FIG_CFG))

save_overview(
  fig2, STAGE, "wtheatup_lens_proximity",
  table = lens_sum,
  finding = paste0(sprintf("Of the %d annotated genes in the mouse 39 °C-derived up arm, the curated lenses contain %d. ",
                           n_annot, sum(lens_sum$n_members_in_query)),
                   paste(sprintf("%s holds %d, and of the %d genes scorable against a depth-matched background %d sit above its 95th percentile against a chance count of %.1f",
                                 lens_sum$lens, lens_sum$n_members_in_query,
                                 lens_sum$n_scored, lens_sum$n_above_p95,
                                 lens_sum$n_expected_by_chance),
                         collapse = "; "), "."),
  script = SCRIPT, fn = "top-level (pD/pE)",
  config_kv = sprintf("semantic.depth_bands = %s, semantic.bg_sample_size = %s, semantic.measure = %s",
                      prov[["depth_bands"]], prov[["bg_sample"]], prov[["measure"]]),
  input = "03_results/objects/20_semantic_decomp.rds",
  how_to_read = paste0(
    "Upper panel: each dot is one gene against one lens, at the percentile of its similarity to the nearest OTHER ",
    "lens member, scored against background genes carrying a similar number of terms. No gene is scored against ",
    "itself, so the red dots act as a control that the measure works. The dashed rule is the 95th percentile, and ",
    "each axis label gives the observed count above it beside the count expected by chance. A lens whose two counts ",
    "match has no graded signal beyond its listed members. Crosses are genes whose depth band already reaches the ",
    "ceiling: drawn, excluded from the counts, tallied on the axis. Lower panel: why the banding is needed."),
  config = FIG_CFG, wide = TRUE, height = 12.5)

message("21_semantic_decomposition_viz: VIZ DONE.")
