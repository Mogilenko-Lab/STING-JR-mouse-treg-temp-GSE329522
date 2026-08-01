#!/usr/bin/env Rscript
# 25_go_decomposition_viz.R -- VIZ ONLY (no computing)
# =============================================================================
# Draws the five overview panels for the ontology-wide over-representation stage
# (15_go_decomposition) from the frozen object 24_go_decomposition.R writes. Every
# number on every face is read from that object; nothing here recomputes a test,
# a null, or a cluster.
#
#   wtheatup_term_blocks        the 35 Wang-similarity blocks the primary arm's
#                               enriched terms collapse to, sized by term count
#   arms_observed_vs_null       both headline counts against the depth-matched
#                               permutation null, per arm. Read this one first.
#   wtheatup_null_recurrence    how often a depth-matched random draw reaches
#                               significance for each term, against the term's rank
#   wtheatup_proteostasis_probe every configured proteostasis probe with its size,
#                               its gene hits and its hypergeometric result
#   arms_coverage_ladder        the three-way per-gene coverage split, per arm,
#                               with the coverage null beside it
#
# Run from project root (after 24_go_decomposition.R):
#   Rscript 02_analysis/scripts/25_go_decomposition_viz.R
# =============================================================================

source("02_analysis/config/config.R")
source("02_analysis/helpers/figure_style.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(ggrepel)
  library(scales)
})
options(stringsAsFactors = FALSE)

STAGE  <- "15_go_decomposition"
SCRIPT <- "02_analysis/scripts/25_go_decomposition_viz.R"
set_paper_style(config = FIG_CFG)

OBJ_PATH <- file.path(DIR_OBJECTS, "24_go_decomposition.rds")
stopifnot(file.exists(OBJ_PATH))
obj <- readRDS(OBJ_PATH)

ARM  <- obj$primary_arm
ONT  <- obj$ontology_primary
IEA  <- obj$iea_primary
ARMS <- obj$arm_order
PCUT <- as.numeric(obj$cutoffs$p)
prov <- stats::setNames(obj$provenance$value, obj$provenance$key)

OI  <- FIG_CFG$colors$okabe_ito
LAB <- (FIG_CFG$figures$label_size %||% 4)
PT  <- (FIG_CFG$figures$point_size %||% 2.4)

CFG_KV <- sprintf(
  "go_decomposition.primary_ontology = %s, go_decomposition.primary_iea = %s, go_decomposition.p_cutoff = %s, go_decomposition.min_gs_size = %s, go_decomposition.max_gs_size = %s, go_decomposition.n_null = %s, go_decomposition.seed = %s",
  ONT, IEA, prov[["p_cutoff"]], prov[["min_gs_size"]], prov[["max_gs_size"]],
  prov[["n_null"]], prov[["seed"]])

split_genes <- function(x) unique(unlist(strsplit(x, "/", fixed = TRUE)))
pfmt <- function(p) format.pval(p, digits = 2, eps = 1 / (as.numeric(prov[["n_null"]]) + 1))

message("25_go_decomposition_viz: primary arm ", ARM, " (", ONT, ", ", IEA, ")")

## ===========================================================================
## Figure 1 -- the term blocks of the primary arm
## ===========================================================================

enr_p <- obj$enriched %>%
  dplyr::filter(.data$arm == ARM, .data$ontology == ONT, .data$iea_variant == IEA)

blocks <- enr_p %>%
  dplyr::group_by(.data$cluster_id, .data$cluster_representative) %>%
  dplyr::summarise(
    n_terms = dplyr::n(),
    n_genes = length(split_genes(.data$geneID)),
    best_p_adjust = min(.data$p.adjust),
    n_terms_p_matched = sum(.data$p_matched < PCUT, na.rm = TRUE),
    median_recurrence = stats::median(.data$frac_matched_reaching_q, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(frac_terms_p_matched = .data$n_terms_p_matched / .data$n_terms) %>%
  dplyr::arrange(dplyr::desc(.data$n_terms), .data$best_p_adjust) %>%
  as.data.frame()

# Two blocks can carry the same representative label; the block id disambiguates them
# so no two rows of the axis read identically.
blocks$label <- blocks$cluster_representative
dupl <- blocks$label %in% blocks$label[duplicated(blocks$label)]
blocks$label[dupl] <- sprintf("%s (block %d)", blocks$label[dupl], blocks$cluster_id[dupl])
blocks$label <- factor(blocks$label, levels = rev(blocks$label))

n_blocks <- nrow(blocks)
arm_n <- obj$arms[[ARM]]$n_nominal

p1a <- ggplot(blocks, aes(x = .data$n_terms, y = .data$label,
                          fill = -log10(.data$best_p_adjust))) +
  geom_col(width = 0.74) +
  geom_text(aes(label = sprintf("%d genes", .data$n_genes)), hjust = -0.16,
            size = LAB * 0.8, colour = "grey20") +
  scale_fill_gradient(low = "grey82", high = OI$blue,
                      name = "best adjusted p\nin the block, -log10",
                      limits = c(0, max(-log10(blocks$best_p_adjust)))) +
  scale_x_continuous(limits = c(0, max(blocks$n_terms) * 1.26),
                     expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(labels = scales::label_wrap(34)) +
  labs(title = sprintf("Enriched terms per block, %d blocks", n_blocks),
       x = sprintf("GO %s terms in the block", ONT), y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
        legend.key.width = unit(1.6, "lines"))

p1b <- ggplot(blocks, aes(x = .data$frac_terms_p_matched, y = .data$label)) +
  geom_segment(aes(x = 0, xend = .data$frac_terms_p_matched, yend = .data$label),
               colour = "grey78", linewidth = 1.0) +
  geom_point(size = PT * 1.5, colour = OI$vermillion) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1),
                     labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.08))) +
  # Short on purpose: patchwork left-aligns a panel title over the panel's own region and
  # a longer one runs off the right edge of the canvas.
  labs(title = "Null verdict",
       x = sprintf("share with p_matched below %.2f", PCUT), y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig1 <- (p1a | p1b) +
  patchwork::plot_layout(widths = c(2.7, 1)) +
  patchwork::plot_annotation(
    title = sprintf("GO %s term blocks for %s, by term count and best adjusted p",
                    ONT, ARM),
    subtitle = sprintf(paste0(
      "The %d enriched terms collapse to %d blocks at Wang similarity %s, average linkage cut at %s (GOSemSim %s). ",
      "Bar length is the number of terms in\nthe block, the bar-end count the union of genes those terms hold, out of the arm's %d. ",
      "The right panel gives the share of the block's terms whose\nhit count a depth-matched random draw of the same size matches in under %.0f%% of %s replicates."),
      nrow(enr_p), n_blocks, prov[["simplify_cutoff"]], prov[["cluster_height"]],
      prov[["GOSemSim_version"]], arm_n, 100 * PCUT, prov[["n_null"]]),
    theme = project_theme(config = FIG_CFG))

save_overview(
  fig1, STAGE, "wtheatup_term_blocks",
  table = blocks[, c("cluster_id", "cluster_representative", "n_terms", "n_genes",
                     "best_p_adjust", "n_terms_p_matched", "frac_terms_p_matched",
                     "median_recurrence")],
  finding = sprintf(paste0(
    "The %d GO %s terms enriched for %s collapse to %d blocks. The largest is %s with %d terms over %d of the arm's %d genes, ",
    "and %.0f%% of its terms have a depth-matched p_matched under %.2f. Across all blocks the union of genes the enriched terms hold is %d."),
    nrow(enr_p), ONT, ARM, n_blocks, blocks$cluster_representative[1],
    blocks$n_terms[1], blocks$n_genes[1], arm_n,
    100 * blocks$frac_terms_p_matched[1], PCUT,
    length(split_genes(enr_p$geneID))),
  script = SCRIPT, fn = "top-level (p1a | p1b)",
  config_kv = CFG_KV, input = "03_results/objects/24_go_decomposition.rds",
  how_to_read = paste0(
    "One row per block. Bar length is the number of enriched terms in the block, fill the smallest adjusted p ",
    "among them, the bar-end text the arm genes those terms cover between them. The right panel is the same ",
    "block's depth-matched verdict: the share of its terms whose observed hit count a matched random draw ",
    "equals or beats in under 5% of 2000 replicates. A long bar with a low dot is a large block a random draw ",
    "reproduces. Blocks are a similarity cut over term annotation, so a gene can sit in several and the gene ",
    "counts do not sum to the arm. Claim tier: descriptive, hypothesis-generating."),
  config = FIG_CFG, wide = TRUE, height = 17)

## ===========================================================================
## Figure 2 -- both headline counts against the depth-matched null
## ===========================================================================

Q_TERMS <- "GO BP terms reaching significance"
Q_GENES <- "arm genes in a term reaching significance"

draws <- obj$null_draws %>%
  dplyr::filter(.data$ontology == ONT, .data$iea_variant == IEA,
                .data$null == "depth_matched")
ts <- obj$term_summary %>%
  dplyr::filter(.data$ontology == ONT, .data$iea_variant == IEA)

obs_base <- draws %>%
  dplyr::distinct(.data$arm, .data$n_signif_obs, .data$n_covered_obs,
                  .data$n_query_annotated)

marks <- dplyr::bind_rows(
  dplyr::left_join(obs_base, ts[, c("arm", "p_n_signif_matched", "null_median_n_signif",
                                "null_q95_n_signif", "uniform_median_n_signif")],
                   by = "arm") %>%
    dplyr::transmute(arm = .data$arm, quantity = Q_TERMS,
                     observed = .data$n_signif_obs, p_value = .data$p_n_signif_matched,
                     null_median = .data$null_median_n_signif,
                     null_q95 = .data$null_q95_n_signif,
                     uniform_median = .data$uniform_median_n_signif,
                     n_drawn = .data$n_query_annotated),
  dplyr::left_join(obs_base, ts[, c("arm", "p_covered_matched", "null_median_covered",
                                "null_q95_covered", "uniform_median_covered")],
                   by = "arm") %>%
    dplyr::transmute(arm = .data$arm, quantity = Q_GENES,
                     observed = .data$n_covered_obs, p_value = .data$p_covered_matched,
                     null_median = .data$null_median_covered,
                     null_q95 = .data$null_q95_covered,
                     uniform_median = .data$uniform_median_covered,
                     n_drawn = .data$n_query_annotated))

long <- dplyr::bind_rows(
  draws %>% dplyr::transmute(arm = .data$arm, quantity = Q_TERMS,
                             value = .data$n_signif),
  draws %>% dplyr::transmute(arm = .data$arm, quantity = Q_GENES,
                             value = .data$covered))
long$arm <- factor(long$arm, levels = ARMS)
long$quantity <- factor(long$quantity, levels = c(Q_TERMS, Q_GENES))
marks$arm <- factor(marks$arm, levels = ARMS)
marks$quantity <- factor(marks$quantity, levels = c(Q_TERMS, Q_GENES))
marks$note <- sprintf("observed %d, permutation p %s\nnull median %.0f, 95th percentile %.0f",
                      marks$observed, pfmt(marks$p_value),
                      marks$null_median, marks$null_q95)

arm_lab <- stats::setNames(
  sprintf("%s\n%d genes drawn", obs_base$arm, obs_base$n_query_annotated), obs_base$arm)

p2 <- ggplot(long, aes(x = .data$value)) +
  geom_histogram(bins = 34, fill = "grey72", colour = NA) +
  geom_vline(data = marks, aes(xintercept = .data$observed, colour = "observed"),
             linewidth = 1.1) +
  geom_point(data = marks, aes(x = .data$uniform_median, y = 0,
                               colour = "uniform null median"),
             shape = 17, size = PT * 1.5) +
  # A label rather than plain text: the observed rule spans the whole panel height, so in
  # the arms where it lands early it runs straight through the annotation.
  geom_label(data = marks, aes(x = -Inf, y = Inf, label = .data$note),
             hjust = -0.02, vjust = 1.05, size = LAB * 0.82, colour = "grey20",
             fill = "white", linewidth = 0, label.padding = unit(0.18, "lines"),
             lineheight = 1.05) +
  facet_grid(arm ~ quantity, scales = "free", switch = "y",
             labeller = labeller(arm = arm_lab, quantity = label_wrap_gen(34))) +
  scale_colour_manual(values = c("observed" = OI$vermillion,
                                 "uniform null median" = OI$sky_blue), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.42))) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.06))) +
  labs(title = "Observed counts against the depth-matched permutation null, per arm",
       subtitle = sprintf(paste0(
         "Each histogram is %s draws of a gene set the size of that arm's annotated part, matched to it band by band on how many GO %s terms each gene\n",
         "carries, scored through the same hypergeometric and the same BH adjustment as the arm. The left column counts terms reaching q below %.2f, the\n",
         "right counts drawn genes falling in at least one of them. The triangle marks the median of a uniformly drawn null of the same size, which is the\n",
         "comparison the annotation-depth matching replaces. Seed %s."),
         prov[["n_null"]], ONT, PCUT, prov[["seed"]]),
       x = "count in one replicate", y = "replicates") +
  project_theme(config = FIG_CFG) +
  theme(legend.position = "bottom", panel.spacing.x = unit(1.1, "lines"),
        panel.spacing.y = unit(0.9, "lines"),
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, hjust = 1))

save_overview(
  p2, STAGE, "arms_observed_vs_null",
  table = marks[, c("arm", "quantity", "n_drawn", "observed", "null_median",
                    "null_q95", "uniform_median", "p_value")],
  finding = sprintf(paste0(
    "For %s the %d enriched GO %s terms sit against a depth-matched null median of %.0f and 95th percentile of %.0f, permutation p %s over %s replicates. ",
    "The %d arm genes in at least one enriched term sit against a null median of %.0f and 95th percentile of %.0f, permutation p %s. ",
    "A uniformly drawn null of the same size reaches a median of %.0f terms, so annotation depth accounts for most of the term count before any biology is read."),
    ARM,
    marks$observed[marks$arm == ARM & marks$quantity == Q_TERMS],
    ONT,
    marks$null_median[marks$arm == ARM & marks$quantity == Q_TERMS],
    marks$null_q95[marks$arm == ARM & marks$quantity == Q_TERMS],
    pfmt(marks$p_value[marks$arm == ARM & marks$quantity == Q_TERMS]),
    prov[["n_null"]],
    marks$observed[marks$arm == ARM & marks$quantity == Q_GENES],
    marks$null_median[marks$arm == ARM & marks$quantity == Q_GENES],
    marks$null_q95[marks$arm == ARM & marks$quantity == Q_GENES],
    pfmt(marks$p_value[marks$arm == ARM & marks$quantity == Q_GENES]),
    marks$uniform_median[marks$arm == ARM & marks$quantity == Q_TERMS]),
  script = SCRIPT, fn = "top-level (p2)",
  config_kv = CFG_KV, input = "03_results/objects/24_go_decomposition.rds",
  how_to_read = paste0(
    "Grey is the null over 2000 depth-matched draws, the red rule the arm's own value, the blue triangle the ",
    "median of a uniformly drawn null of the same size. The permutation p is the share of draws reaching the ",
    "observed value or more, floored at 1/2001. Rows are arms and carry the genes drawn per replicate, which is ",
    "the arm's annotated size. Left column: terms clearing the adjusted-p gate. Right column: drawn genes ",
    "landing in one of those terms, the null for the coverage number the rest of this directory rests on. ",
    "Claim tier: this panel bounds what the rest of the stage can support."),
  config = FIG_CFG, wide = TRUE, height = 12.5)

## ===========================================================================
## Figure 3 -- null recurrence against rank
## ===========================================================================

dn <- obj$depth_null %>%
  dplyr::filter(.data$arm == ARM, .data$ontology == ONT, .data$iea_variant == IEA) %>%
  dplyr::arrange(.data$p_adjust_hypergeometric) %>%
  dplyr::mutate(rank = dplyr::row_number(),
                survives = ifelse(.data$p_matched < PCUT,
                                  sprintf("p_matched below %.2f", PCUT),
                                  sprintf("p_matched at or above %.2f", PCUT))) %>%
  as.data.frame()

med_rec <- stats::median(dn$frac_matched_reaching_q)
n_20 <- sum(dn$frac_matched_reaching_q > 0.20)
n_50 <- sum(dn$frac_matched_reaching_q > 0.50)
head_lab <- dn[seq_len(min(8L, nrow(dn))), ]
y_top <- max(dn$frac_matched_reaching_q) * 1.08

REC_COL <- c(OI$blue, "grey62")
names(REC_COL) <- c(sprintf("p_matched below %.2f", PCUT),
                    sprintf("p_matched at or above %.2f", PCUT))

p3a <- ggplot(dn, aes(x = .data$rank, y = .data$frac_matched_reaching_q)) +
  geom_hline(yintercept = med_rec, linetype = "dashed", linewidth = 0.6,
             colour = "grey35") +
  geom_hline(yintercept = c(0.20, 0.50), linetype = "dotted", linewidth = 0.5,
             colour = "grey55") +
  geom_point(aes(colour = .data$survives), size = PT * 0.8, alpha = 0.75) +
  geom_smooth(method = "loess", formula = y ~ x, se = FALSE, linewidth = 1.1,
              colour = OI$vermillion) +
  ggrepel::geom_text_repel(
    data = head_lab, aes(label = .data$Description), size = LAB * 0.85,
    colour = "grey15", seed = 11L, min.segment.length = 0,
    segment.colour = "grey45", segment.size = 0.4, box.padding = 0.55,
    point.padding = 0.4, direction = "both", xlim = c(1, nrow(dn) * 0.34),
    max.overlaps = Inf) +
  scale_colour_manual(values = REC_COL, name = NULL) +
  # coord_cartesian rather than scale limits: a scale limit drops the histogram bin whose
  # lower edge stat_bin puts below zero, and ggplot reports that as removed rows.
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.03))) +
  coord_cartesian(ylim = c(0, y_top)) +
  labs(title = sprintf("Recurrence against rank, %d terms", nrow(dn)),
       x = sprintf("rank in the %s term list by adjusted p", ARM),
       y = "share of depth-matched\ndraws reaching q below 0.05") +
  project_theme(config = FIG_CFG) +
  theme(legend.position = "bottom")

p3b <- ggplot(dn, aes(y = .data$frac_matched_reaching_q)) +
  geom_histogram(bins = 30, fill = "grey72", colour = NA) +
  geom_hline(yintercept = med_rec, linetype = "dashed", linewidth = 0.6,
             colour = "grey35") +
  geom_hline(yintercept = c(0.20, 0.50), linetype = "dotted", linewidth = 0.5,
             colour = "grey55") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.02))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0, y_top)) +
  labs(title = "All terms", x = "terms", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

fig3 <- (p3a | p3b) +
  patchwork::plot_layout(widths = c(3.2, 1)) +
  patchwork::plot_annotation(
    title = sprintf("Depth-matched null recurrence against rank in the %s term list", ARM),
    subtitle = sprintf(paste0(
      "For each of the %d enriched terms, the share of %s depth-matched random draws that also take it to q below %.2f. Median %.3f, ",
      "%d terms above 20%%\nand %d above 50%%. The dashed rule is the median and the dotted rules are 20%% and 50%%. ",
      "The eight most significant terms are named."),
      nrow(dn), prov[["n_null"]], PCUT, med_rec, n_20, n_50),
    theme = project_theme(config = FIG_CFG))

save_overview(
  fig3, STAGE, "wtheatup_null_recurrence",
  table = dn[, c("rank", "ID", "Description", "term_size", "k_obs",
                 "p_adjust_hypergeometric", "p_matched", "frac_matched_reaching_q",
                 "k_null_matched_mean", "simplify_kept")],
  finding = sprintf(paste0(
    "Over the %d GO %s terms enriched for %s, a depth-matched random draw of the same size reaches q below %.2f for the median term in %.1f%% of %s replicates, ",
    "for %d terms in more than 20%% and for %d terms in more than 50%%. Recurrence rises towards the head of the list: %s."),
    nrow(dn), ONT, ARM, PCUT, 100 * med_rec, prov[["n_null"]], n_20, n_50,
    paste(sprintf("%s %.3f", head_lab$Description[seq_len(min(5L, nrow(head_lab)))],
                  head_lab$frac_matched_reaching_q[seq_len(min(5L, nrow(head_lab)))]),
          collapse = ", ")),
  script = SCRIPT, fn = "top-level (p3a | p3b)",
  config_kv = CFG_KV, input = "03_results/objects/24_go_decomposition.rds",
  how_to_read = paste0(
    "One dot per enriched term, at its rank by adjusted p against how often a gene set matched to the arm on ",
    "annotation depth also takes it to significance over 2000 draws. Blue dots are terms whose observed hit ",
    "count a matched draw equals or beats in under 5% of replicates. The red curve is a loess fit, and the ",
    "right panel is the same y axis as a histogram over all terms. A term high on this axis is one the ",
    "ontology hands to any well-annotated gene set of this size. Claim tier: a property of the term and of ",
    "the background, bounding how the block panel reads."),
  config = FIG_CFG, wide = TRUE, height = 9)

## ===========================================================================
## Figure 4 -- the proteostasis probe as explicit nulls
## ===========================================================================

pr <- obj$proteostasis %>% dplyr::filter(.data$arm == ARM) %>% as.data.frame()
dup_ids <- unique(pr$go_id[duplicated(pr$go_id)])
pr$found_by <- ifelse(pr$go_id %in% dup_ids, "configured id and name sweep",
                      ifelse(pr$probe_source == "name_sweep", "name sweep",
                             "configured id"))
pr <- pr[!(pr$go_id %in% dup_ids & pr$probe_source == "name_sweep"), ]

pr$verdict <- ifelse(!pr$tested, "never entered the test",
                     ifelse(!is.na(pr$p_adjust) & pr$p_adjust < PCUT &
                              pr$enriched %in% TRUE,
                            sprintf("tested, adjusted p below %.2f", PCUT),
                            "tested, adjusted p above the cutoff"))
pr$x <- ifelse(pr$tested, -log10(pr$pvalue), 0)
pr$hits <- ifelse(nzchar(pr$genes_hit), gsub("/", ", ", pr$genes_hit), "no arm gene")
pr$right <- ifelse(pr$tested,
                   sprintf("%d of %d genes: %s", pr$k_arm, pr$term_size_in_universe,
                           pr$hits),
                   # The stage's own status token with the underscores opened out, so the
                   # face and the status column of the table read as the same value.
                   sprintf("%s (%s genes in the universe)",
                           gsub("_", " ", pr$status),
                           ifelse(is.na(pr$term_size_in_universe), "no",
                                  as.character(pr$term_size_in_universe))))
pr$right <- vapply(pr$right, function(s)
  paste(strwrap(s, width = 38), collapse = "\n"), character(1), USE.NAMES = FALSE)
pr$row_label <- sprintf("%s  %s", pr$go_id, pr$go_term)
pr <- pr[order(pr$ontology, pr$tested, -pr$x), ]
pr$row_label <- factor(pr$row_label, levels = rev(unique(pr$row_label)))

VERD_COL <- c("grey45", OI$vermillion, OI$sky_blue)
names(VERD_COL) <- c("tested, adjusted p above the cutoff",
                     sprintf("tested, adjusted p below %.2f", PCUT),
                     "never entered the test")
VERD_SHP <- c(16, 16, 4)
names(VERD_SHP) <- names(VERD_COL)
present <- names(VERD_COL)[names(VERD_COL) %in% pr$verdict]

x_top <- max(c(pr$x, -log10(PCUT))) * 1.14

p4a <- ggplot(pr, aes(x = .data$x, y = .data$row_label)) +
  geom_segment(aes(x = 0, xend = .data$x, yend = .data$row_label),
               colour = "grey85", linewidth = 0.8) +
  geom_vline(xintercept = -log10(PCUT), linetype = "dashed", linewidth = 0.6,
             colour = "grey35") +
  geom_point(aes(colour = .data$verdict, shape = .data$verdict), size = PT * 1.6,
             stroke = 1.1) +
  scale_colour_manual(values = VERD_COL, breaks = present, name = NULL) +
  scale_shape_manual(values = VERD_SHP, breaks = present, name = NULL) +
  scale_x_continuous(limits = c(0, x_top), expand = expansion(mult = c(0.01, 0.02))) +
  scale_y_discrete(labels = scales::label_wrap(40)) +
  facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
  labs(title = "Hypergeometric p for each probe",
       x = "raw hypergeometric p, -log10", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

p4b <- ggplot(pr, aes(x = 0, y = .data$row_label, label = .data$right)) +
  geom_text(hjust = 0, size = LAB * 0.8, colour = "grey20", lineheight = 1.02) +
  facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0.02, 0)),
                     breaks = NULL) +
  labs(title = "Arm genes in the term", x = NULL, y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid = element_blank(), axis.line = element_blank(),
        axis.text.y = element_blank(), axis.ticks = element_blank(),
        strip.text = element_blank())

n_tested <- sum(pr$tested)
sig <- pr[pr$verdict == sprintf("tested, adjusted p below %.2f", PCUT), ]
sig_txt <- if (nrow(sig)) paste(sprintf("%s (%d genes: %s, p_matched %s)", sig$go_term,
                                        sig$k_arm, gsub("/", ", ", sig$genes_hit),
                                        format(round(sig$p_matched, 3))),
                                collapse = " and ") else "none"
chap <- sort(intersect(c("HSPA1A", "HSPH1", "HSPA1B", "DNAJB1"),
                       split_genes(paste(pr$genes_hit, collapse = "/"))))

fig4 <- (p4a | p4b) +
  patchwork::plot_layout(widths = c(2.1, 1)) +
  patchwork::plot_annotation(
    title = sprintf("Hypergeometric result for every proteostasis probe term in %s", ARM),
    subtitle = sprintf(paste0(
      "%d of the %d distinct probe terms entered the hypergeometric. The %d that did not are marked with their reason, which is a property of the term's size in\n",
      "the %s-symbol background or of its gene content, and carries no information about the arm. Two terms clear the adjusted-p cutoff, and the right\n",
      "panel names the arm genes carrying each of them. The dashed rule is a raw p of %.2f. Probe list and name sweep are both configured."),
      n_tested, nrow(pr), nrow(pr) - n_tested, prov[["background_n"]], PCUT),
    theme = project_theme(config = FIG_CFG))

save_overview(
  fig4, STAGE, "wtheatup_proteostasis_probe",
  table = pr[, c("go_id", "go_term", "ontology", "found_by", "status", "tested",
                 "term_size_in_universe", "k_arm", "expected_k", "pvalue", "p_adjust",
                 "p_matched", "frac_matched_reaching_q", "enriched", "genes_hit")],
  finding = sprintf(paste0(
    "Of %d distinct proteostasis probe terms looked up against %s, %d entered the hypergeometric and %d cleared the adjusted-p cutoff. ",
    "Those %d are %s, and the genes carrying them are inflammatory mediators. ",
    "Every chaperone, folding and unfolded-protein term that entered returned an adjusted p above the cutoff, and the chaperone genes present in the arm are %s."),
    nrow(pr), ARM, n_tested, nrow(sig), nrow(sig), sig_txt,
    paste(chap, collapse = " and ")),
  script = SCRIPT, fn = "top-level (p4a | p4b)",
  config_kv = CFG_KV, input = "03_results/objects/24_go_decomposition.rds",
  how_to_read = paste0(
    "One row per probe term, split by ontology. The dot is the raw hypergeometric p on a -log10 axis and the ",
    "dashed rule is p = 0.05. A cross at zero is a term the test never saw, and the right panel gives the ",
    "reason: above the 500-gene cap, below the 10-gene floor, or holding no arm gene at all. For every tested ",
    "term the right panel names the arm genes it holds, so the two terms that clear the cutoff can be read by ",
    "their gene content. Claim tier: these are searched-and-absent rows, reported so an absence carries a ",
    "number."),
  config = FIG_CFG, wide = TRUE, height = 11)

## ===========================================================================
## Figure 5 -- the coverage ladder, three-way per arm
## ===========================================================================

cs <- obj$coverage_summary
CLASSES <- c("in a term reaching significance", "annotated, in no such term",
             sprintf("no GO %s annotation", ONT))
ladder <- dplyr::bind_rows(
  data.frame(arm = cs$arm, class = CLASSES[1], n = cs$n_in_enriched_term),
  data.frame(arm = cs$arm, class = CLASSES[2], n = cs$n_annotated_no_enriched_term),
  data.frame(arm = cs$arm, class = CLASSES[3], n = cs$n_no_annotation))
ladder$arm <- factor(ladder$arm, levels = rev(ARMS))
ladder$class <- factor(ladder$class, levels = CLASSES)
ladder <- ladder %>%
  dplyr::group_by(.data$arm) %>%
  dplyr::arrange(.data$class, .by_group = TRUE) %>%
  dplyr::mutate(mid = cumsum(.data$n) - .data$n / 2) %>%
  dplyr::ungroup() %>% as.data.frame()

tot <- data.frame(arm = factor(cs$arm, levels = rev(ARMS)), n = cs$n_genes,
                  covered = cs$n_in_enriched_term,
                  null_median = cs$null_median_covered,
                  null_q95 = cs$null_q95_covered,
                  p_covered = cs$p_covered_matched)
tot$note <- sprintf("%d genes, %d in a term, null median %.0f, p %s",
                    tot$n, tot$covered, tot$null_median, pfmt(tot$p_covered))
LABEL_FLOOR <- 20L   # below this a count printed inside its segment overflows it

CLS_COL <- c(OI$bluish_green, OI$orange, "grey72")
names(CLS_COL) <- CLASSES

p5 <- ggplot(ladder, aes(x = .data$n, y = .data$arm, fill = .data$class)) +
  # reverse = TRUE stacks in FACTOR order, which is the order `mid` was accumulated in.
  # Left at the default the drawn stack runs the other way and every in-bar count lands
  # over the wrong segment, which reads as a mislabelled figure rather than as an error.
  geom_col(width = 0.62, position = position_stack(reverse = TRUE)) +
  geom_text(data = ladder[ladder$n >= LABEL_FLOOR, ],
            aes(x = .data$mid, label = .data$n), size = LAB * 0.9,
            colour = "grey10", fontface = "bold") +
  geom_point(data = tot, aes(x = .data$null_median, y = .data$arm,
                             shape = "depth-matched null median coverage"),
             inherit.aes = FALSE, size = PT * 3.4, colour = "grey10") +
  geom_text(data = tot, aes(x = .data$n, y = .data$arm, label = .data$note),
            inherit.aes = FALSE, hjust = -0.06, size = LAB * 0.82, colour = "grey25") +
  scale_fill_manual(values = CLS_COL, drop = FALSE, name = NULL) +
  scale_shape_manual(values = c("depth-matched null median coverage" = 124),
                     name = NULL) +
  scale_x_continuous(limits = c(0, max(tot$n) * 1.62),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = sprintf("Where each arm's genes sit relative to the enriched GO %s terms",
                       ONT),
       subtitle = sprintf(paste0(
         "Every gene of every arm lands in exactly one of three classes against the %s-symbol projection background. The tick inside each bar is the median\n",
         "coverage of %s gene sets matched to that arm on annotation depth, and the p beside the bar is the share of those draws covering as many genes as\n",
         "the arm does. Coverage is a count of genes the enriched terms touch, and a gene is counted once however many terms hold it."),
         prov[["background_n"]], prov[["n_null"]]),
       x = "genes", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
        legend.box = "vertical", legend.spacing.y = unit(1, "pt"))

save_overview(
  p5, STAGE, "arms_coverage_ladder",
  table = cs[, c("arm", "n_genes", "n_no_annotation", "n_annotated_no_enriched_term",
                 "n_in_enriched_term", "n_no_bp_but_other_ontology", "n_arm_annotated",
                 "n_covered_obs", "null_median_covered", "null_q95_covered",
                 "uniform_median_covered", "p_covered_matched")],
  finding = sprintf(paste0(
    "%d of %s's %d genes sit in at least one enriched GO %s term, %d carry an annotation and sit in none, and %d carry no %s annotation at all (%d of those are visible in %s). ",
    "A depth-matched random draw of %d genes covers a median of %.0f, so the permutation p on the coverage number is %s."),
    cs$n_in_enriched_term[cs$arm == ARM], ARM, cs$n_genes[cs$arm == ARM], ONT,
    cs$n_annotated_no_enriched_term[cs$arm == ARM], cs$n_no_annotation[cs$arm == ARM],
    ONT, cs$n_no_bp_but_other_ontology[cs$arm == ARM],
    setdiff(unlist(strsplit(prov[["ontologies"]], ",")), ONT)[1],
    cs$n_arm_annotated[cs$arm == ARM], cs$null_median_covered[cs$arm == ARM],
    pfmt(cs$p_covered_matched[cs$arm == ARM])),
  script = SCRIPT, fn = "top-level (p5)",
  config_kv = CFG_KV, input = "03_results/objects/24_go_decomposition.rds",
  how_to_read = paste0(
    "One stacked bar per arm, in gene counts. Green is genes held by at least one term that reached ",
    "significance, orange is genes the ontology can place that no such term holds, grey is genes with no ",
    "annotation in this ontology. Segment counts are printed inside the bar when the segment is wide enough. ",
    "The black tick is the median coverage over 2000 depth-matched draws of the same size, and the text at the ",
    "bar end carries the arm total, its coverage, that median and the permutation p. Read the tick before the ",
    "green segment: coverage this high is what the ontology gives any set of equally well-annotated genes. ",
    "Claim tier: descriptive."),
  config = FIG_CFG, wide = TRUE, height = 7.5)

message("25_go_decomposition_viz: VIZ DONE.")
