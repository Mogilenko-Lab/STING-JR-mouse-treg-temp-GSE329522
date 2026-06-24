# 05_gsea_msigdb_run.R — COMPUTE
## fgsea over the 8 MSigDB collections × all 7 contrasts, ranked by the limma-trend t.
## Stage: 06_gsea.  Run from project root AFTER 04_gsea_set_prep.R and 11_emit_universe.R:
##   Rscript 02_analysis/scripts/05_gsea_msigdb_run.R
##
## ASSUMPTIONS
##   * 04_gsea_set_prep.R has already cached geneset_msigdb_<name>.rds under DIR_OBJECTS.
##     Any collection whose cache file is absent at run time is WARN-skipped (not an error).
##   * 11_emit_universe.R has produced 03_results/objects/gene_universe.txt (background
##     universe, ~19 679 MGI symbols). Presence is checked at start; absence stops with a
##     clear message.
##   * 02_de_limma_trend.R has produced 03_results/objects/02_de_results.rds (named list
##     of 7 limma topTables, one per design.contrast, Symbol rownames, column `t` present).
##   * Inputs are READ-ONLY. Only 03_results/{06_gsea,objects,master}/ are written.
##
## OUTPUTS
##   Objects (one per contrast, full named-list of raw fgsea results):
##     03_results/objects/gsea_msigdb_<contrast>.rds
##   Per-contrast tidy tables (all 8 collections merged):
##     03_results/06_gsea/tables/by_contrast/<contrast>/gsea_msigdb.csv
##   Overview table (all contrasts × all collections):
##     03_results/06_gsea/tables/_overview/gsea_msigdb_all.csv
##   Master table (idempotent append, keyed on `database`):
##     03_results/master/master_gsea_table.csv
##
## IDEMPOTENCY
##   load_or_compute() returns from cache on re-runs (no recompute).
##   append_master_table() drops existing rows keyed on `database` before re-appending
##   (so MSigDB rows replace only their own slice, leaving custom/coresh rows untouched).
##   All numeric columns are rounded to 9 sig figs for byte-stable CSV writes.
##
## COMPUTE ONLY — no ggplot / ggsave / plotting of any kind.

# ============================================================================
# 0. Environment setup
# ============================================================================

source("02_analysis/config/config.R")       # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS,
                                             # DIR_MASTER, SPECIES, GSEA_*, RANK_METRIC, %||%
source("02_analysis/helpers/de_gsea_helpers.R")  # load_or_compute (path-keyed), load_de_results,
                                                  # build_ranked_vector, run_fgsea,
                                                  # append_master_table, round_numeric_cols,
                                                  # .empty_gsea_df

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})
options(stringsAsFactors = FALSE)

STAGE <- "06_gsea"

# ============================================================================
# 1. Pre-flight: validate required upstream inputs
# ============================================================================

## (a) gene_universe.txt — contrast-invariant enrichment background
UNIVERSE_TXT <- file.path(DIR_OBJECTS, "gene_universe.txt")
if (!file.exists(UNIVERSE_TXT)) {
  stop("05_gsea_msigdb_run: gene_universe.txt not found at ", UNIVERSE_TXT,
       "\n  Run 11_emit_universe.R first.")
}
gene_universe <- readLines(UNIVERSE_TXT)
message(sprintf("[05] Gene universe loaded: %d symbols from %s", length(gene_universe), UNIVERSE_TXT))

## (b) 02_de_results.rds — limma-trend topTables
## load_de_results() validates Symbol rownames + presence of `t` column; stops on failure.
de <- load_de_results()
contrasts <- names(de)
message(sprintf("[05] DE results loaded: %d contrasts (%s).",
                length(contrasts), paste(contrasts, collapse = ", ")))

## (c) geneset_manifest.csv — read to discover which MSigDB collections are actually cached.
##     We log this but do NOT stop if manifest is absent: we'll warn-skip individually below.
manifest_path <- file.path(DIR_OBJECTS, "geneset_manifest.csv")
if (file.exists(manifest_path)) {
  manifest_df <- readr::read_csv(manifest_path, show_col_types = FALSE, progress = FALSE)
  message(sprintf("[05] geneset_manifest.csv loaded: %d collection rows.", nrow(manifest_df)))
} else {
  manifest_df <- NULL
  warning("[05] geneset_manifest.csv not found at ", manifest_path,
          " — continuing without manifest (will still attempt to load each cached RDS).")
}

## (d) GSEA parameters (all from config.R constants, never hardcoded here)
MINSZ <- GSEA_MIN_SIZE   # 15    (thresholds.gsea_min_size)
MAXSZ <- GSEA_MAX_SIZE   # 500   (thresholds.gsea_max_size)
SEED  <- GSEA_SEED       # 123   (thresholds.gsea_seed)
NPERM <- GSEA_NPERM      # 1e5   (thresholds.gsea_nperm)  — used only for legacy fgsea fallback
FDR   <- GSEA_FDR_CUTOFF # 0.05  (thresholds.gsea_fdr)    — display only (not used to filter)
RM    <- RANK_METRIC      # "t"   (config.R constant; NEVER logFC)

## (e) Confirm the rank metric column is present in at least one topTable (loud guard)
ex_tt   <- de[[contrasts[[1]]]]
ex_rn   <- rownames(ex_tt)
if (mean(grepl("^ENSMUSG", ex_rn)) > 0.5)
  stop("[05] DE rownames look like Ensembl IDs (", ex_rn[1],
       ") — GSEA needs gene Symbols. Re-run 02_de_limma_trend.R with Symbol row-keying.")
if (!RM %in% colnames(ex_tt))
  stop("[05] Rank metric column '", RM, "' absent from '", contrasts[[1]], "' topTable. ",
       "Check 02_de_results.rds. Columns: ", paste(colnames(ex_tt), collapse = ", "))
message("[05] Symbol rownames confirmed (", ex_rn[1], "); rank metric '", RM, "' present.")

# ============================================================================
# 2. MSigDB collection names (from config; single source of truth)
# ============================================================================

## YAML_CONFIG$databases$msigdb is a list of {category, subcategory, name} records.
## We iterate in config order; names are the cache-file stems.
msigdb_cfg  <- YAML_CONFIG$databases$msigdb
msigdb_names <- vapply(msigdb_cfg, function(m) m$name, character(1))

message(sprintf("[05] MSigDB collections configured: %s", paste(msigdb_names, collapse = ", ")))

# ============================================================================
# 3. Helper: load one cached gene-set RDS or warn+return NULL
# ============================================================================

## ASSUMPTION: 04_gsea_set_prep.R cached each collection as a named list of
## mouse gene-symbol vectors (fgsea `pathways` format), already size-filtered.
## We re-apply minSize/maxSize here as a defensive guard in case the script is
## run against a cache built with different threshold settings.

load_geneset_cache <- function(name) {
  fp <- file.path(DIR_OBJECTS, sprintf("geneset_msigdb_%s.rds", name))
  if (!file.exists(fp)) {
    warning(sprintf("[05] geneset_msigdb_%s.rds not found at %s — skipping collection '%s'.",
                    name, fp, name),
            call. = FALSE)
    return(NULL)
  }
  gs <- readRDS(fp)
  if (!is.list(gs) || length(gs) == 0) {
    warning(sprintf("[05] geneset_msigdb_%s.rds loaded but is empty — skipping collection '%s'.",
                    name, name),
            call. = FALSE)
    return(NULL)
  }
  ## Defensive size-refilter (idempotent when thresholds match 04's run)
  sizes <- vapply(gs, length, integer(1))
  gs    <- gs[sizes >= MINSZ & sizes <= MAXSZ]
  if (length(gs) == 0) {
    warning(sprintf("[05] '%s' has 0 sets after size-filter [%d,%d] — skipping.", name, MINSZ, MAXSZ),
            call. = FALSE)
    return(NULL)
  }
  gs
}

# ============================================================================
# 4. Per-contrast loop: compute (cached) → collect tidy rows
# ============================================================================

## Accumulator: all tidy result rows, keyed for later splitting
all_rows     <- list()   # named paste(contrast,db) -> tidy data.frame
summary_rows <- list()   # named paste(contrast,db) -> one-row summary tibble

## Stage tables root (06_gsea/tables) — stage_dir() creates it
tbl_root <- stage_dir(STAGE, "tables")

## Helper: ensure a subdirectory exists and return its path
ensure_dir <- function(...) {
  p <- file.path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

## Helper: strip genes with NA/empty names from a topTable before building the ranked vector.
## fgsea rejects ranked vectors with NA or empty-string names.
clean_de_table <- function(tt) {
  rn   <- rownames(tt)
  keep <- !is.na(rn) & nzchar(rn)
  tt[keep, , drop = FALSE]
}

for (co in contrasts) {

  message(sprintf("\n== contrast: %s ==", co))

  de_clean <- clean_de_table(de[[co]])

  ## ---- 4a. Build ranked vector once per contrast ----
  ## build_ranked_vector() deduplicates (keeps largest |t|), drops NA/Inf/empty names,
  ## sorts decreasing. Fallback ladder: t -> stat -> sign(logFC)*-log10(P.Value).
  ranked <- build_ranked_vector(de_clean, RM)
  message(sprintf("  ranked vector: %d genes, range [%.2f, %.2f]",
                  length(ranked), min(ranked), max(ranked)))

  ## ---- 4b. Compute cached raw fgsea results (all 8 collections, one list per contrast) ----
  ## Cache key: "gsea_msigdb_<contrast>.rds" — a named list (db_name -> tidy data.frame).
  ## Storing the tidy data.frame (not the raw fgsea data.table) keeps the cache format
  ## stable across fgsea versions and avoids list-column (leadingEdge) serialisation issues.
  ##
  ## NOTE: the cache stores TIDY rows (output of run_fgsea()), not raw fgsea objects,
  ## because run_fgsea() coerces leadingEdge to slash-joined strings (the master contract)
  ## and applies the schema renaming in one place. Re-running always reads from cache unless
  ## force=TRUE is passed to load_or_compute().

  cache_file <- sprintf("gsea_msigdb_%s.rds", co)

  gsea_list <- load_or_compute(
    cache_file,
    function() {
      ## Compute one run_fgsea() call per collection; warn+NULL on missing cache / empty sets.
      ## setNames(lapply(...), msigdb_names) ensures every slot is present (even if NULL),
      ## so a future load_or_compute() cache-hit can check expected_keys cleanly.
      setNames(
        lapply(msigdb_names, function(nm) {
          gs <- load_geneset_cache(nm)
          if (is.null(gs)) return(NULL)     # warn already emitted by load_geneset_cache()
          message(sprintf("  [fgsea] %s (%d sets) ...", nm, length(gs)))
          run_fgsea(
            ranked    = ranked,
            gene_sets = gs,
            database  = nm,
            contrast  = co,
            minSize   = MINSZ,
            maxSize   = MAXSZ,
            seed      = SEED,
            nperm     = NPERM
          )
        }),
        msigdb_names
      )
    }
  )  ## end load_or_compute

  ## ---- 4c. Accumulate tidy rows + summary stats ----
  co_rows <- list()   # this contrast's rows across all collections (for by_contrast CSV)

  for (nm in msigdb_names) {
    rows <- gsea_list[[nm]]

    ## rows is either a well-typed data.frame (possibly 0-row) or NULL (collection skipped)
    is_null  <- is.null(rows)
    is_empty <- !is_null && nrow(rows) == 0L

    if (is_null || is_empty) {
      status <- if (is_null) "skipped" else "empty"
      summary_rows[[paste(co, nm)]] <- tibble::tibble(
        contrast          = co,
        database          = nm,
        n_sets_tested     = 0L,
        n_sig_up          = 0L,
        n_sig_down        = 0L,
        top_pathway_up    = NA_character_,
        top_nes_up        = NA_real_,
        top_pathway_down  = NA_character_,
        top_nes_down      = NA_real_,
        status            = status
      )
      next
    }

    ## Accumulate for master + overview
    key <- paste(co, nm)
    all_rows[[key]] <- rows
    co_rows[[nm]]   <- rows

    ## Summary (display FDR only; does not filter master rows)
    sig  <- dplyr::filter(rows, padj < FDR)
    up   <- dplyr::filter(sig, nes > 0)
    dn   <- dplyr::filter(sig, nes < 0)

    top_up_pw  <- if (nrow(up) > 0) up$pathway_name[which.max(up$nes)]  else NA_character_
    top_up_nes <- if (nrow(up) > 0) max(up$nes)                         else NA_real_
    top_dn_pw  <- if (nrow(dn) > 0) dn$pathway_name[which.min(dn$nes)]  else NA_character_
    top_dn_nes <- if (nrow(dn) > 0) min(dn$nes)                         else NA_real_

    summary_rows[[paste(co, nm)]] <- tibble::tibble(
      contrast          = co,
      database          = nm,
      n_sets_tested     = nrow(rows),
      n_sig_up          = nrow(up),
      n_sig_down        = nrow(dn),
      top_pathway_up    = top_up_pw,
      top_nes_up        = top_up_nes,
      top_pathway_down  = top_dn_pw,
      top_nes_down      = top_dn_nes,
      status            = "ok"
    )

    message(sprintf("  %s: %d sets, %d up / %d down at padj<%.2g",
                    nm, nrow(rows), nrow(up), nrow(dn), FDR))
  }  ## end collection loop

  ## ---- 4d. Write per-contrast table ----
  ## Layout: 03_results/06_gsea/tables/by_contrast/<contrast>/gsea_msigdb.csv
  ## All 8 (or fewer if some skipped) collections for this contrast, sorted by padj.

  by_co_dir <- ensure_dir(tbl_root, "by_contrast", co)

  if (length(co_rows) > 0) {
    co_tbl <- dplyr::bind_rows(co_rows) %>%
      dplyr::arrange(padj, dplyr::desc(abs(nes)))
  } else {
    ## All collections were skipped/empty — write a well-typed 0-row table (not a crash)
    co_tbl <- .empty_gsea_df()
  }

  readr::write_csv(round_numeric_cols(co_tbl), file.path(by_co_dir, "gsea_msigdb.csv"))
  message(sprintf("  -> written: by_contrast/%s/gsea_msigdb.csv  (%d rows)", co, nrow(co_tbl)))

}  ## end contrast loop

# ============================================================================
# 5. Overview table: all contrasts × all collections
# ============================================================================

overview_dir <- ensure_dir(tbl_root, "_overview")

all_tbl <- dplyr::bind_rows(all_rows) %>%
  dplyr::arrange(contrast, padj, dplyr::desc(abs(nes)))

readr::write_csv(round_numeric_cols(all_tbl),
                 file.path(overview_dir, "gsea_msigdb_all.csv"))
message(sprintf("\n[05] Overview: %d total rows -> _overview/gsea_msigdb_all.csv", nrow(all_tbl)))

# ============================================================================
# 6. Summary TSV (per contrast × collection)
# ============================================================================

summary_df <- dplyr::bind_rows(summary_rows) %>%
  dplyr::arrange(contrast, database)

readr::write_csv(round_numeric_cols(summary_df),
                 file.path(tbl_root, "gsea_msigdb_summary.csv"))
message(sprintf("[05] Summary: %d rows -> %s/gsea_msigdb_summary.csv", nrow(summary_df), STAGE))

# ============================================================================
# 7. Append to master_gsea_table.csv (idempotent, keyed on `database`)
# ============================================================================
## append_master_table() CONTRACT (from de_gsea_helpers.R):
##   Key: `database` (default key_col). On each call it drops every row in the master
##   file whose `database` value matches the incoming df$database, then re-appends.
##   This means each database's rows are atomically replaced on re-run, while rows
##   from other databases (e.g. custom DBs added by 06_gsea_custom_run.R) are
##   preserved. MSigDB rows never clobber custom/CoReSh rows and vice versa.
##
## APPEND-KEY STRATEGY (why `database` not `contrast`):
##   Using `database` as the key groups all contrasts for one collection into a single
##   replaceable slice. A re-run of this script for new contrasts or new collections
##   therefore:
##     - replaces ALL rows for each MSigDB database it writes (Hallmark, KEGG, ...)
##     - does NOT touch rows from databases NOT written here (TransportDB, CoReSh, etc.)
##   This is correct because this script owns the full cross-contrast sweep for each
##   MSigDB collection. If only a subset of contrasts were re-run, the old contrast rows
##   for that database would be wiped — acceptable for a compute-only idempotent script.

if (nrow(all_tbl) > 0) {
  ## Write per-database to give append_master_table() the correct granularity
  for (db_nm in msigdb_names) {
    db_rows <- dplyr::filter(all_tbl, database == db_nm)
    if (nrow(db_rows) == 0L) next
    append_master_table(db_rows, "master_gsea_table.csv", key_col = "database")
    message(sprintf("[05] master_gsea_table.csv <- %s: %d rows appended", db_nm, nrow(db_rows)))
  }
} else {
  warning("[05] No MSigDB GSEA rows produced — master_gsea_table.csv NOT updated. ",
          "Check Symbol rownames in 02_de_results.rds and msigdbr install.")
}

# ============================================================================
# 8. Final structural asserts (stop loudly if contract violated)
# ============================================================================

## (a) Per-contrast tables present
for (co in contrasts) {
  p <- file.path(tbl_root, "by_contrast", co, "gsea_msigdb.csv")
  stopifnot("by_contrast gsea_msigdb.csv missing" = file.exists(p))
}

## (b) Overview and summary present
stopifnot(file.exists(file.path(overview_dir, "gsea_msigdb_all.csv")))
stopifnot(file.exists(file.path(tbl_root, "gsea_msigdb_summary.csv")))

## (c) Each contrast's raw-result cache written under DIR_OBJECTS
for (co in contrasts) {
  p <- file.path(DIR_OBJECTS, sprintf("gsea_msigdb_%s.rds", co))
  stopifnot("gsea_msigdb cache RDS missing" = file.exists(p))
}

## (d) Master table schema validation — only if we wrote rows
master_p <- file.path(DIR_MASTER, "master_gsea_table.csv")
if (file.exists(master_p)) {
  req_cols <- YAML_CONFIG$schemas$master_gsea_table$required_columns
  if (!is.null(req_cols)) {
    gt <- readr::read_csv(master_p, show_col_types = FALSE, progress = FALSE)
    missing_cols <- setdiff(req_cols, colnames(gt))
    if (length(missing_cols) > 0)
      stop("[05] master_gsea_table.csv is missing required schema columns: ",
           paste(missing_cols, collapse = ", "))
    message(sprintf("[05] master_gsea_table.csv schema OK (%d rows, %d cols).",
                    nrow(gt), ncol(gt)))
  }
}

# ============================================================================
# 9. Biology sanity checks (warn, never stop)
# ============================================================================
## These are find-not-crash checks. A failure here means "check the data",
## not "the pipeline is broken". For the cGAS-STING / hyperthermia design:
##   WT_heat: genes DE in WT at 39°C vs 37°C — the ISG/IFN arm should appear here.
##   Interaction: (WT_heat) - (KO_heat) — IFN arm should be cGAS-dependent (positive NES);
##                HIF/glycolysis arm should NOT (no significant enrichment in Interaction).
##
## NOTE: contrast names come from config (design.contrasts[*].name); "WT_heat" and
## "Interaction" must match exactly. If they are absent the check silently skips.

if (nrow(all_tbl) > 0 && "WT_heat" %in% contrasts) {
  wt_h_hall <- dplyr::filter(all_tbl,
    contrast == "WT_heat",
    database == "Hallmark",
    grepl("INTERFERON", pathway_id, ignore.case = TRUE)
  )
  if (nrow(wt_h_hall) > 0) {
    ifna <- dplyr::filter(wt_h_hall, grepl("ALPHA", pathway_id, ignore.case = TRUE))
    if (nrow(ifna) > 0) {
      if (ifna$nes[1] > 0 && ifna$padj[1] < FDR)
        message(sprintf("[SANITY OK] WT_heat: HALLMARK_INTERFERON_ALPHA_RESPONSE NES=%.2f padj=%.2e",
                        ifna$nes[1], ifna$padj[1]))
      else
        message(sprintf("[SANITY WARN] WT_heat: HALLMARK_INTERFERON_ALPHA_RESPONSE NES=%.2f padj=%.2e — expected positive & padj<%.2g; check ranking sign / symbol rownames.",
                        ifna$nes[1], ifna$padj[1], FDR))
    } else {
      message("[SANITY WARN] WT_heat: HALLMARK_INTERFERON_ALPHA_RESPONSE not found in Hallmark results.")
    }
  } else {
    message("[SANITY WARN] WT_heat: No Hallmark IFN sets found at all — check msigdbr install or gene-set cache.")
  }
}

if (nrow(all_tbl) > 0 && "Interaction" %in% contrasts) {
  int_hall <- dplyr::filter(all_tbl,
    contrast == "Interaction",
    database == "Hallmark",
    grepl("INTERFERON_ALPHA", pathway_id, ignore.case = TRUE)
  )
  if (nrow(int_hall) > 0) {
    if (int_hall$nes[1] > 0 && int_hall$padj[1] < FDR)
      message(sprintf("[SANITY OK] Interaction: HALLMARK_INTERFERON_ALPHA_RESPONSE NES=%.2f (positive = cGAS-dependent IFN arm, expected)",
                      int_hall$nes[1]))
    else
      message(sprintf("[SANITY WARN] Interaction: HALLMARK_INTERFERON_ALPHA_RESPONSE NES=%.2f padj=%.2e — expected positive at padj<%.2g if IFN arm is cGAS-dependent.",
                      int_hall$nes[1], int_hall$padj[1], FDR))
  }
}

# ============================================================================
# 10. Done
# ============================================================================

message(sprintf(
  "\n05_gsea_msigdb_run complete.\n  Contrasts: %d  Collections: %d  Total rows: %d\n  Objects: %s/gsea_msigdb_<contrast>.rds\n  Tables:  %s/by_contrast/<contrast>/gsea_msigdb.csv\n           %s/_overview/gsea_msigdb_all.csv\n  Master:  %s/master_gsea_table.csv",
  length(contrasts),
  length(msigdb_names),
  nrow(all_tbl),
  DIR_OBJECTS,
  tbl_root,
  tbl_root,
  DIR_MASTER
))
