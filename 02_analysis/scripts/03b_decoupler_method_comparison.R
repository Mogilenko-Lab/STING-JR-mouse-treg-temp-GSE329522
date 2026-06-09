#!/usr/bin/env Rscript
# =============================================================================
# 03b_decoupler_method_comparison.R  --  COMPUTE: decoupleR method comparison
# =============================================================================
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   04_tf
#
# Role:    COMPUTE half of the "normalize-then-visualize" split.  Runs ALL six
#          decoupleR statistics (ulm, mlm, wsum, norm_wsum, corr_wsum,
#          consensus) on the four key contrasts (WT_heat, KO_heat, Temp_main,
#          Interaction) using the CACHED CollecTRI network.  Emits three tidy
#          CSVs and prints a Spearman rank-correlation summary.  Contains NO
#          ggplot() / ggsave() and NO plotting.
#
# Inputs:
#   - 03_results/objects/02_de_results.rds        (t-stat per gene, all contrasts)
#   - 03_results/objects/net_collectri_mouse.rds  (source / target / mor)
#
# Outputs (03_results/04_tf/tables/):
#   - fig3j_allmethods_topTF_data.csv
#   - fig3k_method_rank_divergence_data.csv
#   - fig3k_method_rank_spearman.csv
#
# NETWORK NOTE: decoupleR::get_collectri() is BROKEN here (OmnipathR errors).
#   We load the CACHED net from 00c_prepare_networks.R.
#
# Dependencies: decoupleR, dplyr, tidyr, readr, tibble, purrr
# =============================================================================

source("02_analysis/config/config.R")
load_packages(extra = c("decoupleR", "tidyr", "readr", "tibble", "purrr"))

set.seed(GSEA_SEED %||% 123)

TBL_DIR <- stage_dir("04_tf", "tables")

# Key contrasts for this sub-analysis
KEY_CONTRASTS <- c("WT_heat", "KO_heat", "Temp_main", "Interaction")

# Axis membership (from config)
HIF_AXIS <- c("Hif1a", "Epas1")
IFN_AXIS <- c("Irf1", "Irf3", "Irf7", "Stat1", "Stat2", "Nfkb1", "Rela")

# Six statistics emitted by decouple() + consensus_score = TRUE
ALL_STATS <- c("ulm", "mlm", "wsum", "norm_wsum", "corr_wsum", "consensus")

# How many top TFs to surface per direction (up/down) per (contrast, statistic)
TOPN <- 12

cat("=================================================================\n")
cat("03b: decoupleR method comparison (all 6 statistics, 4 contrasts)\n")
cat("=================================================================\n\n")

# =============================================================================
# 1. LOAD INPUTS + BUILD t-STAT MATRIX
# =============================================================================

cat("[1] Loading DE results + CollecTRI network ...\n")

de <- readRDS(file.path(DIR_OBJECTS, "02_de_results.rds"))
stopifnot(all(KEY_CONTRASTS %in% names(de)))

genes <- de[[1]]$gene_symbol
stopifnot(!any(duplicated(genes)))

mat <- sapply(de, function(d) d$t)
rownames(mat) <- genes
mat[is.na(mat)] <- 0

cat(sprintf("  t-stat matrix: %d genes x %d contrasts\n", nrow(mat), ncol(mat)))

net_ct <- readRDS(file.path(DIR_OBJECTS, "net_collectri_mouse.rds"))
cat(sprintf("  CollecTRI: %d edges, %d TFs (pre-minsize)\n",
            nrow(net_ct), length(unique(net_ct$source))))

# =============================================================================
# 2. RUN decouple() ON EACH KEY CONTRAST (all 6 statistics + consensus)
# =============================================================================

cat("\n[2] Running decouple() with ulm + mlm + wsum (+ consensus) on 4 contrasts ...\n")

# Helper: run decouple on a single-column matrix and return ranked tidy frame
run_all_stats <- function(contrast_name) {
  cat(sprintf("  contrast: %s ...\n", contrast_name))
  m <- mat[, contrast_name, drop = FALSE]

  dec <- decouple(
    mat          = m,
    network      = net_ct,
    .source      = "source",
    .target      = "target",
    statistics   = c("ulm", "mlm", "wsum"),
    args         = list(
      ulm  = list(.mor = "mor", minsize = 5),
      mlm  = list(.mor = "mor", minsize = 5),
      wsum = list(.mor = "mor", minsize = 5)
    ),
    consensus_score = TRUE,
    minsize         = 5,
    show_toy_call   = FALSE
  )

  # Add rank within each statistic (descending score = rank 1 = most activated)
  dec %>%
    mutate(contrast = contrast_name) %>%
    group_by(statistic) %>%
    mutate(
      rank_up   = rank(-score, ties.method = "first"),   # 1 = top activated
      rank_down = rank( score, ties.method = "first")    # 1 = top suppressed
    ) %>%
    ungroup()
}

all_dec <- bind_rows(lapply(KEY_CONTRASTS, run_all_stats))

# Verify we have the expected statistics
got_stats <- sort(unique(all_dec$statistic))
exp_stats  <- sort(ALL_STATS)
stopifnot(all(exp_stats %in% got_stats))
cat(sprintf("  Statistics recovered: %s\n\n", paste(got_stats, collapse = ", ")))

# =============================================================================
# 3. AXIS AND KEY-TF ANNOTATION
# =============================================================================

axis_of <- function(tf) {
  ifelse(tf %in% HIF_AXIS, "HIF",
  ifelse(tf %in% IFN_AXIS, "IFN", "other"))
}

all_dec <- all_dec %>%
  mutate(
    direction = ifelse(score > 0, "Up", "Down"),
    axis      = axis_of(source),
    key_tf    = source %in% KEY_TFS   # KEY_TFS from config.R
  )

# =============================================================================
# 4. fig3j TABLE -- top-12 up + top-12 down per (contrast, statistic)
# =============================================================================

cat("[4] Building fig3j_allmethods_topTF_data.csv ...\n")

fig3j <- bind_rows(lapply(KEY_CONTRASTS, function(cc) {
  bind_rows(lapply(ALL_STATS, function(st) {
    d <- all_dec %>% filter(contrast == cc, statistic == st)
    top_up   <- d %>% arrange(rank_up)   %>% head(TOPN) %>%
                  mutate(direction = "Up",   rank = rank_up)
    top_down <- d %>% arrange(rank_down) %>% head(TOPN) %>%
                  mutate(direction = "Down", rank = rank_down)
    bind_rows(top_up, top_down)
  }))
})) %>%
  select(contrast, statistic, source, score, rank, direction, axis, key_tf) %>%
  distinct(contrast, statistic, source, direction, .keep_all = TRUE) %>%
  arrange(contrast, statistic, direction, rank)

cat(sprintf("  %d rows: top/bottom %d per (contrast x statistic)\n", nrow(fig3j), TOPN))

# =============================================================================
# 5. fig3k TABLE -- focused TF set rank under all statistics on WT_heat + Interaction
# =============================================================================

cat("[5] Building fig3k_method_rank_divergence_data.csv ...\n")

# Focused TF set: union of each statistic's top-10 on WT_heat + HIF/IFN axis always
wt_top10_per_stat <- all_dec %>%
  filter(contrast == "WT_heat") %>%
  group_by(statistic) %>%
  slice_min(rank_up, n = 10) %>%
  pull(source) %>% unique()

focused_tfs <- unique(c(wt_top10_per_stat, HIF_AXIS, IFN_AXIS))
cat(sprintf("  Focused TF set: %d TFs\n", length(focused_tfs)))

fig3k_div <- all_dec %>%
  filter(contrast %in% c("WT_heat", "Interaction"),
         source %in% focused_tfs) %>%
  mutate(rank = rank_up) %>%     # rank ascending by score (1 = most activated)
  select(contrast, source, statistic, rank, score) %>%
  arrange(contrast, source, statistic)

cat(sprintf("  fig3k_method_rank_divergence_data.csv: %d rows\n", nrow(fig3k_div)))

# =============================================================================
# 6. fig3k SPEARMAN -- rank correlation of each statistic vs ULM on WT_heat
# =============================================================================

cat("\n[6] Spearman rank-correlation vs ULM on WT_heat ...\n")

wt_dec <- all_dec %>% filter(contrast == "WT_heat")

# Build a common TF set (those ranked by ULM and each statistic)
ulm_ranks <- wt_dec %>%
  filter(statistic == "ulm") %>%
  select(source, rank_ulm = rank_up)

spearman_rows <- lapply(ALL_STATS, function(st) {
  st_ranks <- wt_dec %>%
    filter(statistic == st) %>%
    select(source, rank_st = rank_up)

  shared <- inner_join(ulm_ranks, st_ranks, by = "source")
  rho <- cor(shared$rank_ulm, shared$rank_st, method = "spearman")
  data.frame(statistic = st, spearman_vs_ulm = round(rho, 4),
             n_tfs = nrow(shared), stringsAsFactors = FALSE)
})

spearman_tab <- do.call(rbind, spearman_rows)
spearman_tab <- spearman_tab[order(-spearman_tab$spearman_vs_ulm), ]

cat("\n  Spearman rank-correlation vs ULM (WT_heat, full TF ranking):\n")
print(spearman_tab, row.names = FALSE)
cat("\n")

# Sanity check: MLM should be the clear outlier (lowest correlation vs ULM)
mlm_rho <- spearman_tab$spearman_vs_ulm[spearman_tab$statistic == "mlm"]
other_rhos <- spearman_tab$spearman_vs_ulm[spearman_tab$statistic != "mlm"]
cat(sprintf("  MLM Spearman=%.4f; all other stats range [%.4f, %.4f] -> MLM outlier: %s\n",
            mlm_rho, min(other_rhos), max(other_rhos),
            ifelse(mlm_rho < min(other_rhos), "YES", "NO")))

# =============================================================================
# 7. Hif1a rank under all 6 statistics on WT_heat + Interaction
# =============================================================================

cat("\n[7] Hif1a ranks across all statistics + contrasts:\n")
hif_summary <- all_dec %>%
  filter(source == "Hif1a",
         contrast %in% c("WT_heat", "Interaction")) %>%
  select(contrast, statistic, rank = rank_up, score) %>%
  arrange(contrast, match(statistic, ALL_STATS))
print(as.data.frame(hif_summary), row.names = FALSE)
cat("\n")

# =============================================================================
# 8. EMIT TIDY TABLES
# =============================================================================

cat("[8] Writing tidy tables ...\n")
emit <- function(df, fname) {
  p <- file.path(TBL_DIR, fname)
  write.csv(df, p, row.names = FALSE)
  cat(sprintf("  [SAVE] %-50s %4d rows x %d cols\n", fname, nrow(df), ncol(df)))
}

emit(fig3j,       "fig3j_allmethods_topTF_data.csv")
emit(fig3k_div,   "fig3k_method_rank_divergence_data.csv")
emit(spearman_tab,"fig3k_method_rank_spearman.csv")

cat("\n[DONE] 03b COMPUTE complete. Run 03b_decoupler_method_comparison_viz.R for figures.\n")
