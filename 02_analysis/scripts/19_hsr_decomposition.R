#!/usr/bin/env Rscript
# 19_hsr_decomposition.R -- COMPUTE
# =============================================================================
# Mouse three-lens heat-response decomposition (stage 12_hsr_decomp).
#
# Purpose
#   Ask what the thresholded WT_heat response (WT activated iTregs, 39 °C vs 37 °C)
#   is made of, by measuring it against two curated lenses:
#     1. WT_heat_up, the empirical up signature produced by the export gate;
#     2. HSR_core / HSR_sensitivity, curated heat-shock-response lenses;
#     3. TCR_activation, the activation / IEG pole.
#   Membership is counted over the WHOLE arm, never over a leading edge: classifying
#   only leading-edge genes and re-testing those subsets would test genes selected
#   because they enriched. What the set is made of and what the response enriches for
#   are two different measurements and this script keeps them apart.
#
# Inputs (frozen; read only)
#   03_results/objects/02_de_results.rds
#   03_results/objects/17_signature_sets.rds
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_mouse_lens.rds
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_human_lens.rds
#   00_data/references/gene_sets/tcr_activation_lens/tcr_activation_mouse.rds
#   00_data/references/gene_sets/tcr_activation_lens/tcr_activation_human.rds
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_gene_taxonomy.csv
#   03_results/human_projection/signatures/WT_heat/WT_heat_up.txt
#   ../<compartment>/03_results/*/tables/ranked_*.tsv   (optional; census only)
#
# Outputs
#   03_results/12_hsr_decomp/tables/hsr_decomp_lens_nes.csv
#   03_results/12_hsr_decomp/tables/hsr_decomp_wtheatup_attribution.csv
#   03_results/12_hsr_decomp/tables/hsr_decomp_summary.csv
#   03_results/12_hsr_decomp/tables/hsr_decomp_overlap.csv
#   03_results/12_hsr_decomp/tables/hsr_decomp_rank_concordance.csv
#   03_results/12_hsr_decomp/tables/hsr_decomp_conditional.csv
#   03_results/12_hsr_decomp/tables/hsr_lens_membership.csv
#   03_results/12_hsr_decomp/tables/_overview/hsr_rank_position_panel.csv
#   03_results/12_hsr_decomp/tables/_overview/gate_projection_bridge.csv
#   03_results/objects/19_hsr_decomp_gsea.rds
#
# Method
#   Uses frozen limma topTables only; no DE model is re-fit. Ranked GSEA uses the
#   pipeline-standard signed t statistic (NES/t > 0 = up at 39 °C) with
#   clusterProfiler::GSEA(by = "fgsea") on symbol-space TERM2GENE and a low
#   minGSSize floor so the 47-gene HSR_core is scored. RNG is seeded with GSEA_SEED.
#   Gene-set attribution and overlap are direct set operations.
#
# Honest ceiling
#   Even the clean HSR core is proteotoxic-stress-general, not fever-specific. The
#   37/39 experimental contrast makes the HSR lens confirmatory for this temperature
#   perturbation, but it does NOT make WT_heat_up causal for fever. Use correlative
#   language only.
#
# Run from project root:
#   Rscript 02_analysis/scripts/19_hsr_decomposition.R
# =============================================================================

source("02_analysis/config/config.R")
source("02_analysis/helpers/de_gsea_helpers.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(tibble)
})
options(stringsAsFactors = FALSE)

STAGE <- "12_hsr_decomp"
SCRIPT <- "02_analysis/scripts/19_hsr_decomposition.R"
TBL_DIR <- stage_dir(STAGE, "tables")
# Same-stem neighbours of the _overview figures live beside them, under tables/_overview/.
TBL_OVW <- file.path(TBL_DIR, "_overview")
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_OVW, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("WT_heat", "Temp_main", "KO_heat")
LENS_TERMS <- c("HSR_core", "HSR_sensitivity", "TCR_activation")
MIN_SIZE <- min(10L, GSEA_MIN_SIZE)
MAX_SIZE <- GSEA_MAX_SIZE
SEED <- GSEA_SEED
RANK_COL <- RANK_METRIC
HONEST_CEILING <- paste(
  "The HSR lens is a curated stress-response reference;",
  "the 37/39 contrast supports the experimental response interpretation, not a human fever claim."
)
ATTR_LEVELS <- c("HSR_core_only", "TCR_activation_only", "shared_both", "neither")

message("=================================================================")
message("19_hsr_decomposition: three-lens decomposition (12_hsr_decomp)")
message("  rank metric: ", RANK_COL, " (NES/t > 0 = up at 39 °C)")
message("  GSEA seed: ", SEED, "  minGSSize: ", MIN_SIZE, "  maxGSSize: ", MAX_SIZE)
message("=================================================================")

de <- load_de_results()
missing_contrasts <- setdiff(CONTRASTS, names(de))
if (length(missing_contrasts) > 0L)
  stop("[19] Missing DE contrasts: ", paste(missing_contrasts, collapse = ", "))

ranked <- lapply(CONTRASTS, function(co) build_ranked_vector(de[[co]], metric = RANK_COL))
names(ranked) <- CONTRASTS
wt_ranked <- ranked[["WT_heat"]]

sig_path <- file.path(DIR_OBJECTS, "17_signature_sets.rds")
if (!file.exists(sig_path)) stop("[19] Missing frozen signature object: ", sig_path)
sig <- readRDS(sig_path)
wt_heat_up <- unique(as.character(sig$sets$WT_heat$up$fdr_logfc))
wt_heat_up <- wt_heat_up[!is.na(wt_heat_up) & wt_heat_up != ""]
if (length(wt_heat_up) == 0L) stop("[19] WT_heat_up fdr_logfc set is empty.")

hsr_path <- "00_data/references/gene_sets/temp_hsr_lens/temp_hsr_mouse_lens.rds"
tcr_path <- "00_data/references/gene_sets/tcr_activation_lens/tcr_activation_mouse.rds"
tax_path <- "00_data/references/gene_sets/temp_hsr_lens/temp_hsr_gene_taxonomy.csv"
for (p in c(hsr_path, tcr_path, tax_path))
  if (!file.exists(p)) stop("[19] Missing frozen input: ", p)

hsr <- readRDS(hsr_path)
tcr <- readRDS(tcr_path)
sets <- list(
  HSR_core = unique(as.character(hsr$HSR_core)),
  HSR_sensitivity = unique(as.character(hsr$HSR_sensitivity)),
  TCR_activation = unique(as.character(tcr$TCR_activation))
)
sets <- lapply(sets, function(x) x[!is.na(x) & x != ""])

if (length(intersect(sets$HSR_core, sets$TCR_activation)) != 0L) {
  warning("[19] HSR_core and TCR_activation are not disjoint; overlap will be shown in hsr_decomp_overlap.csv.")
}

build_term2gene <- function(gene_sets) {
  data.frame(
    term = rep(names(gene_sets), lengths(gene_sets)),
    gene = unlist(gene_sets, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

run_lens_gsea <- function(gene_list, gene_sets, min_size = MIN_SIZE) {
  t2g <- build_term2gene(gene_sets)
  t2g <- dplyr::filter(t2g, gene %in% names(gene_list))
  if (nrow(t2g) == 0L) stop("[19] No gene-set genes overlap the ranked vector.")

  set.seed(SEED)
  res <- clusterProfiler::GSEA(
    geneList = gene_list,
    TERM2GENE = t2g,
    minGSSize = min_size,
    maxGSSize = MAX_SIZE,
    pvalueCutoff = 1,
    pAdjustMethod = "fdr",
    eps = 0,
    by = "fgsea",
    nPermSimple = 1e5,
    seed = TRUE,
    verbose = FALSE
  )
  res@geneSets <- split(t2g$gene, t2g$term)
  res
}

leading_edge_vec <- function(x) {
  if (is.null(x) || is.na(x) || !nzchar(x)) character(0) else strsplit(x, "/", fixed = TRUE)[[1]]
}

tidy_gsea <- function(gsea_obj, contrast) {
  rr <- as.data.frame(gsea_obj@result, stringsAsFactors = FALSE)
  if (nrow(rr) == 0L) {
    return(data.frame(
      term = character(), contrast = character(), set_size = integer(),
      nes = numeric(), pvalue = numeric(), padj = numeric(),
      leading_edge_n = integer(), leading_edge_genes = character(),
      stringsAsFactors = FALSE
    ))
  }
  rr$core_enrichment <- rr$core_enrichment %||% ""
  le <- lapply(rr$core_enrichment, leading_edge_vec)
  data.frame(
    term = rr$ID,
    contrast = contrast,
    set_size = rr$setSize,
    nes = rr$NES,
    pvalue = rr$pvalue,
    padj = rr$p.adjust,
    leading_edge_n = lengths(le),
    leading_edge_genes = vapply(le, paste, character(1), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

gsea_by_contrast <- lapply(CONTRASTS, function(co) {
  message("[19] GSEA: ", co)
  run_lens_gsea(ranked[[co]], sets[LENS_TERMS])
})
names(gsea_by_contrast) <- CONTRASTS
saveRDS(gsea_by_contrast, file.path(DIR_OBJECTS, "19_hsr_decomp_gsea.rds"))

lens_nes <- dplyr::bind_rows(lapply(names(gsea_by_contrast), function(co)
  tidy_gsea(gsea_by_contrast[[co]], co))) %>%
  dplyr::mutate(term = factor(term, levels = LENS_TERMS),
                contrast = factor(contrast, levels = CONTRASTS)) %>%
  dplyr::arrange(contrast, term) %>%
  dplyr::mutate(term = as.character(term), contrast = as.character(contrast))
readr::write_csv(round_numeric_cols(lens_nes, sig = 9),
                 file.path(TBL_DIR, "hsr_decomp_lens_nes.csv"))

taxonomy <- readr::read_csv(tax_path, show_col_types = FALSE, progress = FALSE)
tax_long <- taxonomy %>%
  dplyr::select(category, mouse_symbols) %>%
  dplyr::filter(!is.na(mouse_symbols), mouse_symbols != "") %>%
  tidyr::separate_rows(mouse_symbols, sep = "[,;/|[:space:]]+") %>%
  dplyr::filter(mouse_symbols != "") %>%
  dplyr::distinct(mouse_symbols, category) %>%
  dplyr::group_by(mouse_symbols) %>%
  dplyr::summarise(hsr_category = paste(sort(unique(category)), collapse = ";"),
                   .groups = "drop")
tax_map <- stats::setNames(tax_long$hsr_category, tax_long$mouse_symbols)

attrib_for <- function(genes, thermal_set, activation_set) {
  in_hsr <- genes %in% thermal_set
  in_act <- genes %in% activation_set
  dplyr::case_when(
    in_hsr & !in_act ~ "HSR_core_only",
    !in_hsr & in_act ~ "TCR_activation_only",
    in_hsr & in_act ~ "shared_both",
    TRUE ~ "neither"
  )
}

attrib <- data.frame(
  gene = wt_heat_up,
  in_HSR_core = wt_heat_up %in% sets$HSR_core,
  in_TCR_activation = wt_heat_up %in% sets$TCR_activation,
  in_HSR_sensitivity = wt_heat_up %in% sets$HSR_sensitivity,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    attribution = attrib_for(gene, sets$HSR_core, sets$TCR_activation),
    hsr_category = unname(tax_map[gene]),
    hsr_category = ifelse(is.na(hsr_category), "", hsr_category),
    wt_heat_t = unname(wt_ranked[gene]),
    attribution = factor(attribution, levels = ATTR_LEVELS)
  ) %>%
  dplyr::arrange(attribution, dplyr::desc(wt_heat_t), gene) %>%
  dplyr::mutate(attribution = as.character(attribution))
readr::write_csv(round_numeric_cols(attrib, sig = 9),
                 file.path(TBL_DIR, "hsr_decomp_wtheatup_attribution.csv"))

summarise_attrib <- function(thermal_name, thermal_set) {
  att <- attrib_for(wt_heat_up, thermal_set, sets$TCR_activation)
  total <- length(wt_heat_up)
  data.frame(lens = thermal_name, attribution = ATTR_LEVELS, stringsAsFactors = FALSE) %>%
    dplyr::left_join(as.data.frame(table(factor(att, levels = ATTR_LEVELS)),
                                   stringsAsFactors = FALSE) %>%
                       dplyr::rename(attribution = Var1, n = Freq),
                     by = "attribution") %>%
    dplyr::mutate(n = as.integer(n %||% 0L),
                  fraction = n / total,
                  denominator = total,
                  honest_ceiling = HONEST_CEILING)
}

summary_df <- dplyr::bind_rows(
  summarise_attrib("HSR_core", sets$HSR_core),
  summarise_attrib("HSR_sensitivity", sets$HSR_sensitivity)
)
readr::write_csv(round_numeric_cols(summary_df, sig = 9),
                 file.path(TBL_DIR, "hsr_decomp_summary.csv"))

overlap_sets <- list(
  WT_heat_up = wt_heat_up,
  HSR_core = sets$HSR_core,
  HSR_sensitivity = sets$HSR_sensitivity,
  TCR_activation = sets$TCR_activation
)
pairs <- utils::combn(names(overlap_sets), 2, simplify = FALSE)
overlap_df <- dplyr::bind_rows(lapply(pairs, function(p) {
  a <- unique(overlap_sets[[p[1]]])
  b <- unique(overlap_sets[[p[2]]])
  n_i <- length(intersect(a, b))
  data.frame(
    set_a = p[1], set_b = p[2],
    n_a = length(a), n_b = length(b),
    n_intersect = n_i,
    jaccard = n_i / length(union(a, b)),
    stringsAsFactors = FALSE
  )
}))
readr::write_csv(round_numeric_cols(overlap_df, sig = 9),
                 file.path(TBL_DIR, "hsr_decomp_overlap.csv"))

rank_df <- data.frame(
  gene = names(wt_ranked),
  rank = seq_along(wt_ranked),
  t = unname(wt_ranked),
  stringsAsFactors = FALSE
)
rank_df$rank_pct <- if (nrow(rank_df) > 1L) (rank_df$rank - 1) / (nrow(rank_df) - 1) else 0
rank_summary_one <- function(set_name, genes) {
  sub <- rank_df[rank_df$gene %in% genes, , drop = FALSE]
  data.frame(
    set = set_name,
    n_in_ranking = nrow(sub),
    median_rank = stats::median(sub$rank),
    median_rank_pct = stats::median(sub$rank_pct),
    median_t = stats::median(sub$t),
    stringsAsFactors = FALSE
  )
}
rank_conc <- dplyr::bind_rows(
  rank_summary_one("HSR_core", sets$HSR_core),
  rank_summary_one("TCR_activation", sets$TCR_activation),
  rank_summary_one("WT_heat_up", wt_heat_up)
)
hsr_ranks <- rank_df$rank[rank_df$gene %in% sets$HSR_core]
tcr_ranks <- rank_df$rank[rank_df$gene %in% sets$TCR_activation]
wilcox_p <- stats::wilcox.test(hsr_ranks, tcr_ranks, alternative = "two.sided", exact = FALSE)$p.value
higher_set <- if (stats::median(hsr_ranks) < stats::median(tcr_ranks)) {
  "HSR_core"
} else if (stats::median(tcr_ranks) < stats::median(hsr_ranks)) {
  "TCR_activation"
} else {
  "tie"
}
rank_conc <- rank_conc %>%
  dplyr::mutate(wilcox_p = wilcox_p, higher_set = higher_set)
readr::write_csv(round_numeric_cols(rank_conc, sig = 9),
                 file.path(TBL_DIR, "hsr_decomp_rank_concordance.csv"))

# F2 source table: gene-level positions in the WT_heat signed-t ranking plus
# row-summary and empirical gate-span metadata. The guard proves the plotted gate
# and frozen WT_heat_up set are the same object.
md_path <- file.path("03_results", "03_de", "tables", "by_contrast", "WT_heat", "md.csv")
if (!file.exists(md_path)) stop("[19] Missing WT_heat DE table for rank-position panel: ", md_path)
md <- readr::read_csv(md_path, show_col_types = FALSE, progress = FALSE) %>%
  dplyr::arrange(dplyr::desc(.data$t)) %>%
  dplyr::mutate(rank = dplyr::row_number(),
                n_ranking = dplyr::n(),
                rank_pct = 100 * .data$rank / .data$n_ranking)
if (nrow(md) != 19679L)
  stop("[19] WT_heat ranking has ", nrow(md), " genes; expected 19679 for the F2 axis.")

gate_md <- md %>%
  dplyr::filter(.data$adj.P.Val < sig$thresholds$de_fdr,
                .data$logFC >= sig$thresholds$de_logfc)
if (nrow(gate_md) != 213L)
  stop("[19] WT_heat_up gate guard failed: md.csv gives ", nrow(gate_md), " genes, expected 213.")
if (!setequal(gate_md$gene_symbol, wt_heat_up))
  stop("[19] WT_heat_up gate guard failed: md.csv gate genes do not match 17_signature_sets.rds.")

gate_deepest <- gate_md %>%
  dplyr::arrange(dplyr::desc(.data$rank)) %>%
  dplyr::slice(1)
if (!identical(gate_deepest$gene_symbol[1], "Cpne6") || gate_deepest$rank[1] != 2010L)
  stop("[19] WT_heat_up gate span changed: deepest gate gene is ",
       gate_deepest$gene_symbol[1], " at rank ", gate_deepest$rank[1],
       "; expected Cpne6 at rank 2010.")

rank_sets <- list(
  WT_heat_up = wt_heat_up,
  HSR_core = sets$HSR_core,
  TCR_activation = sets$TCR_activation
)
rank_detail <- dplyr::bind_rows(lapply(names(rank_sets), function(nm) {
  md %>%
    dplyr::filter(.data$gene_symbol %in% rank_sets[[nm]]) %>%
    dplyr::transmute(set = nm, gene = .data$gene_symbol, rank = .data$rank,
                     rank_pct = .data$rank_pct, t = .data$t,
                     logFC = .data$logFC, adj.P.Val = .data$adj.P.Val)
}))
rank_position_panel <- rank_detail %>%
  dplyr::left_join(rank_conc %>%
                     dplyr::transmute(set, n_in_ranking, median_rank,
                                      median_rank_pct = 100 * .data$median_rank_pct,
                                      median_t),
                   by = "set") %>%
  dplyr::mutate(
    set_label = dplyr::case_when(
      .data$set == "WT_heat_up" ~ sprintf("WT_heat_up (gate output; n=%d)", .data$n_in_ranking),
      .data$set == "HSR_core" ~ sprintf("HSR_core (n=%d)", .data$n_in_ranking),
      .data$set == "TCR_activation" ~ sprintf("TCR_activation (n=%d)", .data$n_in_ranking),
      TRUE ~ .data$set
    ),
    gate_pass_n = nrow(gate_md),
    gate_span_n = gate_deepest$rank[1],
    gate_min_rank = min(gate_md$rank),
    gate_max_rank = max(gate_md$rank),
    gate_min_rank_pct = min(gate_md$rank_pct),
    gate_max_rank_pct = max(gate_md$rank_pct),
    gate_deepest_gene = gate_deepest$gene_symbol[1],
    gate_deepest_t = gate_deepest$t[1],
    gate_deepest_logFC = gate_deepest$logFC[1],
    gate_inner_pass_n = nrow(gate_md),
    gate_inner_total_n = gate_deepest$rank[1],
    ranking_total_n = nrow(md)
  ) %>%
  dplyr::arrange(factor(.data$set, levels = c("WT_heat_up", "HSR_core", "TCR_activation")),
                 .data$rank, .data$gene)
readr::write_csv(round_numeric_cols(rank_position_panel, sig = 9),
                 file.path(TBL_OVW, "hsr_rank_position_panel.csv"))
# The stage previously wrote this table flat; remove the stale copy so a reader
# cannot pick up a file no run maintains any more.
stale_rank <- file.path(TBL_DIR, "hsr_rank_position_panel.csv")
if (file.exists(stale_rank)) file.remove(stale_rank)

get_nes <- function(gsea_obj, term) {
  rr <- as.data.frame(gsea_obj@result, stringsAsFactors = FALSE)
  rr <- rr[rr$ID == term, , drop = FALSE]
  if (nrow(rr) == 0L) return(list(nes = NA_real_, padj = NA_real_))
  list(nes = rr$NES[1], padj = rr$p.adjust[1])
}

conditional_rows <- list()
for (target in c("HSR_core", "TCR_activation")) {
  conditioned_on <- if (identical(target, "HSR_core")) "TCR_activation" else "HSR_core"
  uncond <- get_nes(gsea_by_contrast$WT_heat, target)
  cond_rank <- wt_ranked[!(names(wt_ranked) %in% sets[[conditioned_on]])]
  cond_gsea <- run_lens_gsea(cond_rank, sets[target])
  cond <- get_nes(cond_gsea, target)
  conditional_rows[[length(conditional_rows) + 1L]] <- data.frame(
    term = target,
    conditioned_on = conditioned_on,
    nes_uncond = uncond$nes,
    nes_cond = cond$nes,
    delta_nes = cond$nes - uncond$nes,
    padj_cond = cond$padj,
    stringsAsFactors = FALSE
  )
}
conditional <- dplyr::bind_rows(conditional_rows)
readr::write_csv(round_numeric_cols(conditional, sig = 9),
                 file.path(TBL_DIR, "hsr_decomp_conditional.csv"))

# F3 source table: three-set membership counts for WT_heat_up / HSR_core /
# TCR_activation, with the human WT_heat_up-HSR_core overlap carried only when the
# external compartment table is present and well-formed.
overlap_n <- function(a, b) {
  rr <- overlap_df %>%
    dplyr::filter((.data$set_a == a & .data$set_b == b) |
                    (.data$set_a == b & .data$set_b == a))
  if (nrow(rr) != 1L) stop("[19] Missing overlap row for ", a, " and ", b)
  rr$n_intersect[1]
}
wt_hsr <- overlap_n("WT_heat_up", "HSR_core")
wt_tcr <- overlap_n("WT_heat_up", "TCR_activation")
hsr_tcr <- overlap_n("HSR_core", "TCR_activation")
triple <- length(Reduce(intersect, list(wt_heat_up, sets$HSR_core, sets$TCR_activation)))

human_overlap_path <- "/workspaces/STING-JR/human_treg_arthritis/03_results/10_hsr_lens/tables/hsr_wtheatup_overlap.csv"
human_wt_heatup_n <- NA_real_
human_hsr_n <- NA_real_
human_wt_hsr_intersect <- NA_real_
human_overlap_note <- "human overlap table absent or malformed"
if (file.exists(human_overlap_path)) {
  human_overlap <- tryCatch(readr::read_csv(human_overlap_path, show_col_types = FALSE, progress = FALSE),
                            error = function(e) NULL)
  if (!is.null(human_overlap) &&
      all(c("set_a", "set_b", "n_a", "n_b", "n_intersect") %in% names(human_overlap))) {
    hr <- human_overlap %>%
      dplyr::filter((.data$set_a == "WT_heat_up" & .data$set_b == "HSR_core") |
                      (.data$set_a == "HSR_core" & .data$set_b == "WT_heat_up")) %>%
      dplyr::slice(1)
    if (nrow(hr) == 1L) {
      human_wt_heatup_n <- ifelse(hr$set_a[1] == "WT_heat_up", hr$n_a[1], hr$n_b[1])
      human_hsr_n <- ifelse(hr$set_a[1] == "HSR_core", hr$n_a[1], hr$n_b[1])
      human_wt_hsr_intersect <- hr$n_intersect[1]
      human_overlap_note <- human_overlap_path
    }
  }
}

conditional_wide <- conditional %>%
  dplyr::select("term", "delta_nes", "nes_uncond", "nes_cond") %>%
  tidyr::pivot_wider(names_from = "term",
                     values_from = c("delta_nes", "nes_uncond", "nes_cond"))

membership <- data.frame(
  component = c("WT_heat_up only", "WT_heat_up ∩ HSR_core", "WT_heat_up ∩ TCR_activation",
                "HSR_core only", "TCR_activation only", "HSR_core ∩ TCR_activation",
                "shared by all three"),
  n_mouse = c(length(wt_heat_up) - wt_hsr - wt_tcr + triple,
              wt_hsr - triple,
              wt_tcr - triple,
              length(sets$HSR_core) - wt_hsr - hsr_tcr + triple,
              length(sets$TCR_activation) - wt_tcr - hsr_tcr + triple,
              hsr_tcr - triple,
              triple),
  glyph = c("region", "overlap", "overlap", "region", "region", "empty_lens_overlap", "empty_triple"),
  wt_heat_up_n = length(wt_heat_up),
  hsr_core_n = length(sets$HSR_core),
  tcr_activation_n = length(sets$TCR_activation),
  human_wt_heatup_n = human_wt_heatup_n,
  human_hsr_n = human_hsr_n,
  human_wt_hsr_intersect = human_wt_hsr_intersect,
  human_overlap_source = human_overlap_note,
  stringsAsFactors = FALSE
) %>%
  dplyr::mutate(
    hsr_core_tcr_activation_intersect = hsr_tcr,
    conditional_delta_nes_HSR_core = conditional_wide$delta_nes_HSR_core[1],
    conditional_delta_nes_TCR_activation = conditional_wide$delta_nes_TCR_activation[1],
    conditional_nes_uncond_HSR_core = conditional_wide$nes_uncond_HSR_core[1],
    conditional_nes_cond_HSR_core = conditional_wide$nes_cond_HSR_core[1],
    conditional_nes_uncond_TCR_activation = conditional_wide$nes_uncond_TCR_activation[1],
    conditional_nes_cond_TCR_activation = conditional_wide$nes_cond_TCR_activation[1]
  )
readr::write_csv(round_numeric_cols(membership, sig = 9),
                 file.path(TBL_DIR, "hsr_lens_membership.csv"))

# =============================================================================
# The handoff ledger: mouse gate -> human projection -> what is testable downstream.
#
# Two numbers travel through this project side by side and are easy to conflate:
# the 213-gene mouse gate output and the 199-gene human set it becomes. This table
# reconciles them, carries each side's curated-lens membership against the SAME two
# lenses in each species, and then measures how much of the human set is actually
# present in the human ranked lists that consume it.
#
# The downstream census is a read-only sweep of the sibling compartments and is
# OPTIONAL: a standalone clone of this repository has no siblings, so the census
# degrades to zero lists and the panel says so rather than failing. Whatever it
# finds is recorded with its file count and length range, so the number is
# checkable rather than remembered.
# =============================================================================
human_set_path <- file.path("03_results", "human_projection", "signatures",
                            "WT_heat", "WT_heat_up.txt")
if (!file.exists(human_set_path))
  stop("[19] Missing frozen human projection set: ", human_set_path)
human_wt_up <- unique(trimws(readLines(human_set_path, warn = FALSE)))
human_wt_up <- human_wt_up[nzchar(human_wt_up)]
if (length(human_wt_up) == 0L) stop("[19] Human WT_heat_up set is empty: ", human_set_path)

hsr_human_path <- "00_data/references/gene_sets/temp_hsr_lens/temp_hsr_human_lens.rds"
tcr_human_path <- "00_data/references/gene_sets/tcr_activation_lens/tcr_activation_human.rds"
for (p in c(hsr_human_path, tcr_human_path))
  if (!file.exists(p)) stop("[19] Missing frozen human lens: ", p)
clean_set <- function(x) {
  x <- unique(as.character(x))
  x[!is.na(x) & nzchar(x)]
}
human_lens <- list(
  HSR_core = clean_set(readRDS(hsr_human_path)$HSR_core),
  TCR_activation = clean_set(readRDS(tcr_human_path)$TCR_activation)
)

# Read-only census of every ranked list the human compartments have published.
# Column 1 is the gene symbol; the files carry no header.
census_roots <- c(
  human_treg_arthritis = file.path("..", "human_treg_arthritis", "03_results",
                                   "03_pseudobulk", "tables"),
  human_pbmc_febrile   = file.path("..", "human_pbmc_febrile", "03_results",
                                   "03_pseudobulk", "tables"),
  human_ra_synovium    = file.path("..", "human_ra_synovium", "03_results",
                                   "04_pseudobulk_de", "tables"),
  sting_positive_control = file.path("..", "sting_positive_control", "03_results",
                                     "03_pseudobulk", "tables")
)
# A ranked list shorter than this is a partial write, not a short list; skip it
# rather than let a half-written file move a published count.
MIN_RANKED_LEN <- 1000L
read_ranked <- function(path) {
  con <- tryCatch(readr::read_tsv(path, col_names = FALSE, col_types = readr::cols(.default = "c"),
                                  progress = FALSE),
                  error = function(e) NULL)
  if (is.null(con) || nrow(con) < MIN_RANKED_LEN) return(NULL)
  clean_set(trimws(con[[1]]))
}
census_lists <- list()
for (comp in names(census_roots)) {
  files <- sort(list.files(census_roots[[comp]], pattern = "^ranked_.*\\.tsv$",
                           full.names = TRUE))
  for (f in files) {
    genes <- read_ranked(f)
    if (is.null(genes)) next
    census_lists[[f]] <- list(compartment = comp, genes = genes)
  }
}
n_census_lists <- length(census_lists)
if (n_census_lists == 0L)
  message("[19] Downstream census: no ranked list found under ../<compartment>/03_results/. ",
          "The bridge ledger will record the mouse->human half only.")

census_rows <- if (n_census_lists > 0L) {
  per_list <- dplyr::bind_rows(lapply(names(census_lists), function(f) {
    g <- census_lists[[f]]$genes
    data.frame(compartment = census_lists[[f]]$compartment,
               list_file = basename(f),
               list_len = length(g),
               n_eff = length(intersect(human_wt_up, g)),
               stringsAsFactors = FALSE)
  }))
  per_list %>%
    dplyr::group_by(.data$compartment) %>%
    dplyr::summarise(n_lists = dplyr::n(),
                     list_len_min = min(.data$list_len),
                     list_len_max = max(.data$list_len),
                     n_eff_min = min(.data$n_eff),
                     n_eff_max = max(.data$n_eff),
                     .groups = "drop")
} else {
  data.frame(compartment = character(), n_lists = integer(),
             list_len_min = integer(), list_len_max = integer(),
             n_eff_min = integer(), n_eff_max = integer(),
             stringsAsFactors = FALSE)
}
common_universe <- if (n_census_lists > 0L)
  Reduce(intersect, lapply(census_lists, `[[`, "genes")) else character(0)
n_testable_everywhere <- length(intersect(human_wt_up, common_universe))

bridge_row <- function(block, step, side, label, n_genes, denominator, note,
                       n_lists = NA_integer_, list_len_min = NA_integer_,
                       list_len_max = NA_integer_) {
  data.frame(
    block = block, step = step, side = side, label = label,
    n_genes = as.integer(n_genes),
    denominator = as.integer(denominator),
    pct_of_denominator = if (is.na(denominator) || denominator == 0L) NA_real_
                         else 100 * n_genes / denominator,
    n_lists = as.integer(n_lists),
    list_len_min = as.integer(list_len_min),
    list_len_max = as.integer(list_len_max),
    note = note,
    stringsAsFactors = FALSE
  )
}

bridge <- dplyr::bind_rows(
  bridge_row("funnel", 1L, "mouse", "WT_heat ranking (every measured gene)",
             nrow(md), nrow(md),
             "signed-t ranking of the WT 39-vs-37 °C contrast"),
  bridge_row("funnel", 2L, "mouse", "inside the gate's rank span (rank 1 to 2,010)",
             gate_deepest$rank[1], nrow(md),
             sprintf("deepest admitted gene %s at rank %d", gate_deepest$gene_symbol[1],
                     gate_deepest$rank[1])),
  bridge_row("funnel", 3L, "mouse", "passed the gate: WT_heat_up",
             nrow(gate_md), gate_deepest$rank[1],
             sprintf("gate is adj.P < %g and logFC >= %g",
                     sig$thresholds$de_fdr, sig$thresholds$de_logfc)),
  bridge_row("funnel", 4L, "human", "after mouse-to-human ortholog projection",
             length(human_wt_up), nrow(gate_md),
             "frozen at 03_results/human_projection/signatures/WT_heat/WT_heat_up.txt"),
  bridge_row("funnel", 5L, "human",
             sprintf("present in every one of the %d human ranked lists", n_census_lists),
             n_testable_everywhere, length(human_wt_up),
             if (n_census_lists > 0L)
               sprintf("common universe across the %d lists is %d symbols",
                       n_census_lists, length(common_universe))
             else "no ranked list was reachable from this checkout",
             n_lists = n_census_lists),
  bridge_row("lens", 3L, "mouse",
             sprintf("of WT_heat_up, in curated HSR_core (%d genes)", length(sets$HSR_core)),
             wt_hsr, nrow(gate_md), "mouse lens, direct set membership"),
  bridge_row("lens", 3L, "mouse",
             sprintf("of WT_heat_up, in curated TCR_activation (%d genes)",
                     length(sets$TCR_activation)),
             wt_tcr, nrow(gate_md), "mouse lens, direct set membership"),
  bridge_row("lens", 4L, "human",
             sprintf("of WT_heat_up, in curated HSR_core (%d genes)",
                     length(human_lens$HSR_core)),
             length(intersect(human_wt_up, human_lens$HSR_core)), length(human_wt_up),
             "human ortholog lens, direct set membership"),
  bridge_row("lens", 4L, "human",
             sprintf("of WT_heat_up, in curated TCR_activation (%d genes)",
                     length(human_lens$TCR_activation)),
             length(intersect(human_wt_up, human_lens$TCR_activation)), length(human_wt_up),
             "human ortholog lens, direct set membership"),
  if (nrow(census_rows) > 0L) dplyr::bind_rows(lapply(seq_len(nrow(census_rows)), function(i) {
    r <- census_rows[i, ]
    dplyr::bind_rows(
      bridge_row("downstream_min", 5L, "human", r$compartment, r$n_eff_min,
                 length(human_wt_up), "fewest of the 199 present in any one of this compartment's lists",
                 n_lists = r$n_lists, list_len_min = r$list_len_min, list_len_max = r$list_len_max),
      bridge_row("downstream_max", 5L, "human", r$compartment, r$n_eff_max,
                 length(human_wt_up), "most of the 199 present in any one of this compartment's lists",
                 n_lists = r$n_lists, list_len_min = r$list_len_min, list_len_max = r$list_len_max))
  })) else NULL
)
readr::write_csv(round_numeric_cols(bridge, sig = 9),
                 file.path(TBL_OVW, "gate_projection_bridge.csv"))

expected <- c(file.path(TBL_DIR, c(
  "hsr_decomp_lens_nes.csv",
  "hsr_decomp_wtheatup_attribution.csv",
  "hsr_decomp_summary.csv",
  "hsr_decomp_overlap.csv",
  "hsr_decomp_rank_concordance.csv",
  "hsr_decomp_conditional.csv",
  "hsr_lens_membership.csv")),
  file.path(TBL_OVW, c("hsr_rank_position_panel.csv", "gate_projection_bridge.csv")))
missing <- expected[!file.exists(expected)]
if (length(missing) > 0L) stop("[19] Missing expected output(s): ", paste(missing, collapse = ", "))

message(sprintf("[19] Bridge ledger: mouse gate %d -> human %d; HSR_core %d mouse / %d human; ",
                nrow(gate_md), length(human_wt_up), wt_hsr,
                length(intersect(human_wt_up, human_lens$HSR_core))),
        sprintf("TCR_activation %d mouse / %d human; %d of %d testable across %d ranked lists.",
                wt_tcr, length(intersect(human_wt_up, human_lens$TCR_activation)),
                n_testable_everywhere, length(human_wt_up), n_census_lists))
message("[19] COMPLETE: wrote ", length(expected), " tables and 03_results/objects/19_hsr_decomp_gsea.rds")
print(summary_df)
print(dplyr::filter(lens_nes, contrast == "WT_heat", term %in% c("HSR_core", "TCR_activation")))
