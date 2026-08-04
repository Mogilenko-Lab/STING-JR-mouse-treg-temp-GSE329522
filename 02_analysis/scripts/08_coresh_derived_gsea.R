# 08_coresh_derived_gsea.R — COMPUTE
## GSEA on CoReSh-derived gene sets (stage 08_coresh) against the per-contrast DE ranking.
##
## Inputs:
##   03_results/objects/02_de_results.rds          — per-contrast limma-trend topTables
##                                                   (Symbol rownames, `t` rank metric; hub)
##   03_results/objects/coresh_derived_sets.rds    — named list of fgsea-format gene sets
##                                                   produced by 07_coresh_search.R
##                                                   SHAPE: named list<character vector of mouse symbols>
##                                                     list("CORESH_<query>_<GSE>" = c("Ifit1","Isg15",...), ...)
##                                                   Names follow build_coresh_gmt() convention:
##                                                   "CORESH_<query_name>_<GSE_id>" (e.g. CORESH_Q_curated_isg_ifn_GSE12345).
##                                                   Any non-empty named list of symbol vectors is valid.
##
## Outputs (compute-only; no plots):
##   03_results/08_coresh/tables/by_contrast/<contrast>/coresh_gsea.csv  — per-contrast GSEA (padj-ordered)
##   03_results/08_coresh/tables/_overview/coresh_gsea_all.csv           — all contrasts combined
##   03_results/master/master_gsea_table.csv                             — idempotent append (database="CoReSh_derived")
##   03_results/objects/gsea_coresh_<contrast>.rds                       — per-contrast fgsea cache
##
## Run from project root:
##   Rscript 02_analysis/scripts/08_coresh_derived_gsea.R
##
## GATED: the CoReSh arm requires the ~20 GB mmu Synapse compendium (consumed read-only
##   from the shared reference cache); if it has not been run this script stops loudly.

# ============================================================================
# 0. Environment setup (config.R FIRST, then de_gsea_helpers.R — per contract)
# ============================================================================

source("02_analysis/config/config.R")           # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS,
                                                  # DIR_MASTER, SPECIES, stage_dir(),
                                                  # load_or_compute (config.R variant),
                                                  # GSEA_MIN_SIZE/MAX_SIZE/SEED/NPERM/FDR,
                                                  # RANK_METRIC = "t", %||%
source("02_analysis/helpers/de_gsea_helpers.R") # path-keyed load_or_compute, load_de_results,
                                                  # build_ranked_vector, run_fgsea,
                                                  # append_master_table, round_numeric_cols,
                                                  # direction_from_nes

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})
options(stringsAsFactors = FALSE)

STAGE   <- "08_coresh"
tbl_dir <- stage_dir(STAGE, "tables")            # 03_results/08_coresh/tables/ (asserted in config)

# Layout subdirectory names from figures: config block (with fallback).
BY_CONTRAST_DIR <- YAML_CONFIG$figures$by_contrast_dir %||% "by_contrast"
OVERVIEW_DIR    <- YAML_CONFIG$figures$overview_dir    %||% "_overview"

# ============================================================================
# 1. GUARD — the CoReSh-derived sets must exist (NEVER fabricate)
# ============================================================================
## This script is DOWNSTREAM of 07_coresh_search.R, which requires the ~20 GB mmu
## Synapse compendium (syn66227307), consumed read-only from the shared reference cache.
## Until 07 has run, coresh_derived_sets.rds will not exist and we stop loudly.

coresh_sets_fp <- file.path(DIR_OBJECTS, "coresh_derived_sets.rds")

if (!file.exists(coresh_sets_fp)) {
  stop(
    "CoReSh-derived sets not found — run 07_coresh_search.R first ",
    "(requires the mmu compendium mounted read-only from the shared reference cache)."
  )
}

# ============================================================================
# 2. Load the CoReSh-derived gene sets
# ============================================================================
## EXPECTED SHAPE: a named list<character vector> (the fgsea `pathways` format).
##   * Names — build_coresh_gmt() convention: "CORESH_<query_name>_<GSE_id>"
##     (e.g. "CORESH_Q_curated_isg_ifn_GSE12345"). Any non-empty named list of
##     mouse-symbol vectors is valid; the name is used as pathway_id downstream.
##   * Values — character vectors of mouse gene symbols (already size-filtered to
##     [gsea_min_size, gsea_max_size] and Jaccard-deduped in 07_coresh_search.R).
##   * The list may be empty (length 0) if all derived sets failed the size/Jaccard
##     filter — this is a soft-exit case handled below.

coresh_sets <- readRDS(coresh_sets_fp)

## Contract: named list. Error loudly if the shape is wrong (prevents silent 0-hit GSEA).
if (!is.list(coresh_sets) || is.null(names(coresh_sets)) || any(!nzchar(names(coresh_sets)))) {
  stop(
    "coresh_derived_sets.rds has unexpected shape (expected a NAMED list of symbol vectors). ",
    "Re-run 07_coresh_search.R to rebuild."
  )
}

## Soft-exit: 07 may have produced an empty list if no sets passed size/Jaccard filter.
if (length(coresh_sets) == 0L) {
  warning(
    "coresh_derived_sets.rds is empty — 07_coresh_search.R produced 0 derived sets ",
    "(all failed size/Jaccard filter). Writing empty summary; nothing to enrich."
  )
  empty_summary <- tibble::tibble(
    contrast      = character(),
    n_sets_tested = integer(),
    n_sig_fdr     = integer(),
    n_up_sig      = integer(),
    n_down_sig    = integer(),
    top_set       = character(),
    top_nes       = double()
  )
  readr::write_csv(empty_summary, file.path(tbl_dir, OVERVIEW_DIR, "coresh_gsea_all.csv"))
  message("08_coresh_derived_gsea: 0 derived sets — done (empty summary written).")
  quit(save = "no", status = 0)
}

n_sets <- length(coresh_sets)
n_genes <- length(unique(unlist(coresh_sets, use.names = FALSE)))
message(sprintf("[08_coresh_derived_gsea] CoReSh-derived sets: %d sets, %d unique symbols.",
                n_sets, n_genes))

## Drop any sets with 0 symbols (fgsea errors on empty pathways).
nonempty <- vapply(coresh_sets, function(g) length(g) > 0L, logical(1))
if (any(!nonempty)) {
  warning(sprintf("Dropping %d empty set(s) from coresh_derived_sets.", sum(!nonempty)))
  coresh_sets <- coresh_sets[nonempty]
}

# ============================================================================
# 3. Load the DE hub (NEVER re-fit DE)
# ============================================================================
## 02_de_results.rds: named list<topTable data.frame>.
##   Each element has gene-SYMBOL rownames + columns: gene_symbol, logFC, AveExpr,
##   t, P.Value, adj.P.Val, B, contrast.
##   Verified by load_de_results() (Symbol-rowname contract + `t` column asserts).

de_results <- load_de_results()
stopifnot(length(de_results) >= 1L)
message(sprintf("[08_coresh_derived_gsea] DE hub: %d contrasts (%s).",
                length(de_results), paste(names(de_results), collapse = ", ")))

## Symbol-rownames guard (CoReSh sets are symbol-keyed; Ensembl rownames → 0 hits).
## load_de_results() already asserts this but we echo it explicitly for CoReSh.
ex_rn <- rownames(de_results[[1L]])
if (mean(grepl("^ENSMUSG", ex_rn)) > 0.5) {
  stop(
    "DE topTable rownames look like Ensembl IDs (e.g. ", ex_rn[1L], "). ",
    "CoReSh-derived sets use mouse symbols — re-key the topTables to gene symbols and re-run."
  )
}

# ============================================================================
# 4. GSEA parameters (from config via helpers; never hardcoded)
# ============================================================================
## .gsea_*() helpers defined in de_gsea_helpers.R read config.R constants first,
## then fall back to YAML_CONFIG$thresholds — single source of truth.

## Use the config.R constant directly; GSEA_FDR_CUTOFF is set there from YAML_CONFIG$thresholds$gsea_fdr.
## Re-bind under a local shorter alias for readability; both names refer to the same value.
GSEA_FDR <- GSEA_FDR_CUTOFF   # = YAML_CONFIG$thresholds$gsea_fdr %||% 0.05  (display threshold, not run cutoff)

## pick contrasts: all YAML-declared contrasts that also appear in the DE hub
focal <- if (!is.null(YAML_CONFIG$design$contrasts)) {
  yaml_names <- vapply(YAML_CONFIG$design$contrasts, function(x) x$name %||% "", character(1))
  intersect(yaml_names[nzchar(yaml_names)], names(de_results))
} else {
  names(de_results)
}
message(sprintf("[08_coresh_derived_gsea] contrasts to run: %s", paste(focal, collapse = ", ")))

# ============================================================================
# 5. Per-contrast fgsea (cached) -> per-contrast CSV
# ============================================================================
## For each contrast:
##   (a) build a ranked vector (gene symbol -> t-statistic, decreasing).
##   (b) run_fgsea() from de_gsea_helpers.R using the CoReSh-derived sets.
##       run_fgsea() calls fgsea::fgseaMultilevel internally (falls back to fgsea::fgsea
##       on older fgsea installs), seeds from gsea_seed, and returns a tidy data.frame
##       in the master_gsea_table schema (lowercase `nes`, direction).
##   (c) cache the result as gsea_coresh_<contrast>.rds (load_or_compute).
##   (d) write 03_results/08_coresh/tables/by_contrast/<contrast>/coresh_gsea.csv.

set.seed(GSEA_SEED)   # GSEA_SEED from config.R; run_fgsea also sets the seed internally
all_results <- list()

for (co in focal) {

  # (a) ranked vector: gene symbol -> t-statistic
  ranked <- build_ranked_vector(de_results[[co]])
  if (is.null(ranked) || length(ranked) == 0L) {
    message(sprintf("  [%s] skip: build_ranked_vector returned empty vector.", co))
    next
  }

  # (b)+(c) run fgsea — cached per contrast
  cache_fp <- file.path(DIR_OBJECTS, sprintf("gsea_coresh_%s.rds", co))

  # load_or_compute keys on filename only — it cannot see that the CoReSh-derived set
  # inventory changed (07 re-derived them from a different query set). Force a recompute
  # when the cache does not cover every current derived set (mirrors 07's force_sweep).
  force_gsea <- file.exists(cache_fp) &&
    !all(names(coresh_sets) %in% unique(as.character(readRDS(cache_fp)$pathway_id)))

  gsea_df <- load_or_compute(cache_fp, function() {
    run_fgsea(
      ranked    = ranked,
      gene_sets = coresh_sets,
      database  = "CoReSh_derived",
      contrast  = co
      # minSize/maxSize/eps/nperm/seed resolved from config inside run_fgsea()
    )
  }, force = force_gsea)

  if (is.null(gsea_df) || nrow(gsea_df) == 0L) {
    message(sprintf("  [%s] 0 CoReSh-derived sets enriched.", co))
    next
  }

  n_sig <- sum(gsea_df$padj < GSEA_FDR, na.rm = TRUE)
  message(sprintf("  [%s] %d sets, %d sig (padj < %.2f).", co, nrow(gsea_df), n_sig, GSEA_FDR))
  all_results[[co]] <- gsea_df

  # (d) per-contrast CSV
  cdir <- file.path(tbl_dir, BY_CONTRAST_DIR, co)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(
    round_numeric_cols(dplyr::arrange(gsea_df, padj)),
    file.path(cdir, "coresh_gsea.csv")
  )
}

# ============================================================================
# 6. Soft-exit if no contrasts produced results
# ============================================================================
if (length(all_results) == 0L) {
  warning(
    "No CoReSh-derived sets enriched in any contrast ",
    "(all fgsea calls returned 0 rows). Writing empty overview; master unchanged."
  )
  empty_summary <- tibble::tibble(
    contrast      = character(),
    n_sets_tested = integer(),
    n_sig_fdr     = integer(),
    n_up_sig      = integer(),
    n_down_sig    = integer(),
    top_set       = character(),
    top_nes       = double()
  )
  overview_dir <- file.path(tbl_dir, OVERVIEW_DIR)
  dir.create(overview_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(empty_summary, file.path(overview_dir, "coresh_gsea_all.csv"))
  message("08_coresh_derived_gsea: no sets enriched in any contrast — done.")
  quit(save = "no", status = 0)
}

# ============================================================================
# 7. Combine across contrasts -> _overview/coresh_gsea_all.csv
# ============================================================================

all_df <- dplyr::bind_rows(all_results)
all_df <- round_numeric_cols(all_df)

overview_dir <- file.path(tbl_dir, OVERVIEW_DIR)
dir.create(overview_dir, recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  dplyr::arrange(all_df, contrast, padj),
  file.path(overview_dir, "coresh_gsea_all.csv")
)
message(sprintf("[08_coresh_derived_gsea] overview -> %s (%d rows across %d contrasts).",
                file.path(overview_dir, "coresh_gsea_all.csv"),
                nrow(all_df), length(all_results)))

# ============================================================================
# 8. Append to master_gsea_table.csv (idempotent by database = "CoReSh_derived")
# ============================================================================
## append_master_table() (de_gsea_helpers.R):
##   1. Validates df against schemas.master_gsea_table.required_columns (config).
##   2. Drops any existing rows where database == "CoReSh_derived" (idempotent replace).
##   3. Binds the new rows, rounds to 9 sig figs (byte-stable), writes back.
## master_gsea_table required schema (from config schemas.master_gsea_table):
##   pathway_id, pathway_name, database, nes, pvalue, padj, set_size,
##   core_enrichment, contrast, direction   (all present in run_fgsea() output).

append_master_table(all_df, "master_gsea_table.csv", key_col = "database")
message(sprintf("[08_coresh_derived_gsea] master_gsea_table.csv updated (%d CoReSh_derived rows).",
                nrow(all_df)))

# ============================================================================
# 9. Summary table (across contrasts)
# ============================================================================

summary_df <- all_df %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(
    n_sets_tested = dplyr::n(),
    n_sig_fdr     = sum(padj < GSEA_FDR, na.rm = TRUE),
    n_up_sig      = sum(direction == "Up"   & padj < GSEA_FDR, na.rm = TRUE),
    n_down_sig    = sum(direction == "Down" & padj < GSEA_FDR, na.rm = TRUE),
    top_set       = pathway_id[which.max(abs(nes))],
    top_nes       = nes[which.max(abs(nes))],
    .groups       = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_sig_fdr))

readr::write_csv(summary_df, file.path(overview_dir, "coresh_gsea_summary.csv"))
message(sprintf("[08_coresh_derived_gsea] summary -> %s.",
                file.path(overview_dir, "coresh_gsea_summary.csv")))

# ============================================================================
# 10. Final structural asserts (fail loudly rather than silently skip)
# ============================================================================
## Every contracted output file must exist after this script runs.

stopifnot(
  # per-contrast CSVs (one per contrast that produced results)
  all(vapply(names(all_results), function(co)
    file.exists(file.path(tbl_dir, BY_CONTRAST_DIR, co, "coresh_gsea.csv")),
    logical(1))),
  # overview table
  file.exists(file.path(overview_dir, "coresh_gsea_all.csv")),
  file.exists(file.path(overview_dir, "coresh_gsea_summary.csv")),
  # master table
  file.exists(file.path(DIR_MASTER, "master_gsea_table.csv")),
  # per-contrast caches
  all(vapply(names(all_results), function(co)
    file.exists(file.path(DIR_OBJECTS, sprintf("gsea_coresh_%s.rds", co))),
    logical(1)))
)

message(sprintf(
  "[08_coresh_derived_gsea] COMPLETE — %d contrasts, %d CoReSh_derived rows appended to master. ",
  length(all_results), nrow(all_df)),
  sprintf("(captions floor at L3-DE)"))
