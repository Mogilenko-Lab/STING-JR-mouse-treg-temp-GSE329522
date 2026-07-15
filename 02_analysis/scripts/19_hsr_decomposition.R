#!/usr/bin/env Rscript
# 19_hsr_decomposition.R -- COMPUTE
# =============================================================================
# Mouse three-lens heat-response decomposition (stage 12_hsr_decomp).
#
# Purpose
#   Decompose the WT_heat response (WT activated iTregs, 39 C vs 37 C) into:
#     1. WT_heat_up, the empirical up signature and activation-heavy by construction;
#     2. HSR_core / HSR_sensitivity, curated thermal proteostasis lenses;
#     3. TCR_activation, the activation / IEG pole.
#
# Inputs (frozen; read only)
#   03_results/objects/02_de_results.rds
#   03_results/objects/17_signature_sets.rds
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_mouse_lens.rds
#   00_data/references/gene_sets/tcr_activation_lens/tcr_activation_mouse.rds
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_gene_taxonomy.csv
#
# Outputs
#   03_results/12_hsr_decomp/tables/hsr_decomp_lens_nes.csv
#   03_results/12_hsr_decomp/tables/hsr_decomp_wtheatup_attribution.csv
#   03_results/12_hsr_decomp/tables/hsr_decomp_summary.csv
#   03_results/12_hsr_decomp/tables/hsr_decomp_overlap.csv
#   03_results/12_hsr_decomp/tables/hsr_decomp_rank_concordance.csv
#   03_results/12_hsr_decomp/tables/hsr_decomp_conditional.csv
#   03_results/objects/19_hsr_decomp_gsea.rds
#
# Method
#   Uses frozen limma topTables only; no DE model is re-fit. Ranked GSEA uses the
#   pipeline-standard signed t statistic (NES/t > 0 = up at 39 C) with
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
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)

CONTRASTS <- c("WT_heat", "Temp_main", "KO_heat")
LENS_TERMS <- c("HSR_core", "HSR_sensitivity", "TCR_activation")
MIN_SIZE <- min(10L, GSEA_MIN_SIZE)
MAX_SIZE <- GSEA_MAX_SIZE
SEED <- GSEA_SEED
RANK_COL <- RANK_METRIC
HONEST_CEILING <- paste(
  "The HSR lens is proteotoxic-stress-general, not fever-specific;",
  "37/39 contrast supports correlative temperature-response interpretation, not causal fever claims."
)
ATTR_LEVELS <- c("thermal_HSR", "activation", "shared_both", "neither")

message("=================================================================")
message("19_hsr_decomposition: three-lens decomposition (12_hsr_decomp)")
message("  rank metric: ", RANK_COL, " (NES/t > 0 = up at 39 C)")
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
    in_hsr & !in_act ~ "thermal_HSR",
    !in_hsr & in_act ~ "activation",
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

expected <- file.path(TBL_DIR, c(
  "hsr_decomp_lens_nes.csv",
  "hsr_decomp_wtheatup_attribution.csv",
  "hsr_decomp_summary.csv",
  "hsr_decomp_overlap.csv",
  "hsr_decomp_rank_concordance.csv",
  "hsr_decomp_conditional.csv"
))
missing <- expected[!file.exists(expected)]
if (length(missing) > 0L) stop("[19] Missing expected output(s): ", paste(missing, collapse = ", "))

message("[19] COMPLETE: wrote ", length(expected), " tables and 03_results/objects/19_hsr_decomp_gsea.rds")
print(summary_df)
print(dplyr::filter(lens_nes, contrast == "WT_heat", term %in% c("HSR_core", "TCR_activation")))
