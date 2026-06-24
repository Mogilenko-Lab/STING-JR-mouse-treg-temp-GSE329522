# 06_gsea_custom_run.R — COMPUTE
## clusterProfiler::GSEA (by="fgsea") over the custom databases (TransportDB,
## MitoPathways, MitoXplorer, Lombardi2022_HIF) × all 7 contrasts. Caches gseaResult
## S4 objects (consumable by the RNAseq-toolkit GSEA plotters); derives the tidy master
## rows via normalize_gsea_results().
##
## EACH DATABASE IS RUN SEPARATELY — its own GSEA call, its own `database` tag,
## its own append key in master_gsea_table.csv — so custom rows coexist with
## MSigDB rows from 05_gsea_msigdb_run.R without clobbering them.
##
## Run from project root (AFTER 04_gsea_set_prep.R, 05_gsea_msigdb_run.R):
##   Rscript 02_analysis/scripts/06_gsea_custom_run.R
##
## Outputs
## -------
##   03_results/06_gsea/tables/by_contrast/<contrast>/gsea_custom.csv
##       One file per contrast containing results from ALL custom DBs present,
##       combined and sorted by padj. The `database` column identifies the DB.
##   03_results/06_gsea/tables/_overview/gsea_custom_all.csv
##       Same rows as all per-contrast files combined, plus the contrast column.
##   03_results/master/master_gsea_table.csv
##       Idempotent APPEND, keyed per-database (e.g. "TransportDB"), so MSigDB
##       rows written by 05 are never touched and each custom DB replaces only
##       its own rows on re-run.
##   03_results/objects/gsea_custom_<contrast>.rds
##       Cached result for the contrast (named list: db_name -> clusterProfiler gseaResult S4,
##       with @geneSets patched in so the running-sum viz can walk custom-set membership).
##
## Idempotency:  load_or_compute() returns cached RDS on re-run.
## Byte-stability: round_numeric_cols() (9 sig figs) before every CSV write.
## Compute-only:  no ggplot/ggsave anywhere in this file.
##
## ASSUMPTIONS (verify in-container; noted where consequential):
##   [A1] 02_de_results.rds is a named list of limma-trend topTables, each with
##        gene-SYMBOL ROWNAMES and a `t` column. Produced by 02_de_limma_trend.R.
##   [A2] geneset_custom_<Name>.rds files are named-list gene-set objects produced
##        by 04_gsea_set_prep.R (shape (b) of load_custom_geneset: a named list of
##        character vectors, already size-filtered to gsea_min/max_size). The RDS
##        name key uses the `name` field from databases.custom in the YAML config.
##   [A3] master_gsea_table.csv may or may not pre-exist (05 writes it first; this
##        script appends). If absent, append_master_table() creates it.
##   [A4] gene_universe.txt (from 11_emit_universe.R) is advisory background; it is
##        read and reported but fgsea ranks on the full DE t-vector (no gene-set
##        subsetting by universe here — fgsea naturally intersects with each pathway).
##   [A5] Lombardi2022_HIF is an ORTHOGONAL published-HIF check (Ratcliffe/Mole,
##        doi:10.1016/j.celrep.2022.111652; 48-gene consensus, human→mouse,
##        protein-coding, provenance-stamped by 00b_curate_lombardi_hif.R). Its NES
##        on the Interaction contrast is the direct contrast to the hand-made 16-gene
##        Biomni list. Per house rules: results floored at L3; never crown HIF1α/2α;
##        "no detectable cGAS-dependence at n=5" — NEVER "cGAS-independent".

## ============================================================================
## 0. Setup
## ============================================================================

source("02_analysis/config/config.R")       # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS,
                                             # DIR_MASTER, GSEA_*, RANK_METRIC, %||%,
                                             # load_or_compute, stage_dir
source("02_analysis/helpers/de_gsea_helpers.R")  # load_de_results, build_ranked_vector,
                                                  # run_fgsea (now unused), append_master_table,
                                                  # load_custom_geneset, round_numeric_cols,
                                                  # direction_from_nes, .empty_gsea_df

## Toolkit GSEA processing helpers (gseaResult -> tidy normalizer + pathway-name prettifier).
RTK <- file.path(PROJECT_ROOT, YAML_CONFIG$paths$rnaseq_toolkit %||% "01_modules/RNAseq-toolkit")
source(file.path(RTK, "scripts", "GSEA", "GSEA_plotting", "format_pathway_names.R"))   # format_pathway_name()
source(file.path(RTK, "scripts", "GSEA", "GSEA_processing", "normalize_gsea.R"))       # normalize_gsea_results(), empty_gsea_tibble()

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr)
  library(clusterProfiler)   # GSEA() -> gseaResult S4 (by="fgsea")
  library(org.Mm.eg.db)      # mouse annotation backing clusterProfiler
})
options(stringsAsFactors = FALSE)

## ============================================================================
## 1. Config constants (single source of truth — never hardcoded)
## ============================================================================

STAGE     <- "06_gsea"
FDR_DISP  <- GSEA_FDR_CUTOFF      # display threshold (0.05); run always keeps all (pvalue_cutoff=1)
MINSZ     <- GSEA_MIN_SIZE         # 15  (from config.R / YAML thresholds)
MAXSZ     <- GSEA_MAX_SIZE         # 500
SEED      <- GSEA_SEED             # 123
NPERM     <- GSEA_NPERM            # 100000  (nPermSimple for clusterProfiler::GSEA)
RANK_COL  <- (YAML_CONFIG$gsea$rank_metric %||% RANK_METRIC)   # "t"
## clusterProfiler::GSEA run-time engine settings (gsea: block; §5 of the blast-radius doc)
PCUT      <- YAML_CONFIG$gsea$pvalue_cutoff_run %||% 1       # keep ALL pathways at run time
PADJM     <- YAML_CONFIG$gsea$padj_method       %||% "fdr"   # pAdjustMethod
EPS       <- YAML_CONFIG$gsea$eps                %||% 0       # exact fgsea p-values

## The cached object CLASS changed (tidy data.frame -> clusterProfiler gseaResult S4).
## load_or_compute() is path-keyed and would return a stale tidy cache, so delete any
## pre-existing gsea_custom_<contrast>.rds up front to force a clean recompute.
.stale_custom <- list.files(DIR_OBJECTS, pattern = "^gsea_custom_.*\\.rds$", full.names = TRUE)
if (length(.stale_custom) > 0) {
  file.remove(.stale_custom)
  message(sprintf("[06] Removed %d stale gsea_custom_*.rds cache(s) (object class changed to gseaResult).",
                  length(.stale_custom)))
}

## Stage table directories
tbl_dir     <- stage_dir(STAGE, "tables")          # 03_results/06_gsea/tables/ (created)
overview_dir <- file.path(tbl_dir, "_overview")
dir.create(overview_dir, recursive = TRUE, showWarnings = FALSE)

## ============================================================================
## 2. Custom-DB registry — built from databases.custom in the YAML config.
##    Each entry: list(name = "<DisplayName>", rds = "<absolute cache path>").
##    04_gsea_set_prep.R stores the cache as geneset_custom_<name>.rds where
##    <name> is the `name` field from the YAML (case-preserved, e.g. "Lombardi2022_HIF").
## ============================================================================

custom_cfg <- YAML_CONFIG$databases$custom
if (is.null(custom_cfg) || length(custom_cfg) == 0)
  stop("databases.custom is empty in analysis_config.yaml. Register the 4 custom DBs.")

## Build the registry, skipping any whose 04-produced cache is absent.
db_registry <- list()
for (entry in custom_cfg) {
  nm       <- entry$name                                    # e.g. "TransportDB"
  rds_path <- file.path(DIR_OBJECTS, sprintf("geneset_custom_%s.rds", nm))
  if (!file.exists(rds_path)) {
    warning(sprintf(
      "[06_gsea_custom_run] SKIP '%s': cached gene set not found at\n  %s\n  Run 04_gsea_set_prep.R first.",
      nm, rds_path
    ))
    next
  }
  db_registry[[nm]] <- rds_path
}
DB_NAMES <- names(db_registry)   # e.g. c("TransportDB","MitoPathways","MitoXplorer","Lombardi2022_HIF")

if (length(DB_NAMES) == 0)
  stop("No custom gene-set caches found under ", DIR_OBJECTS,
       " — run 04_gsea_set_prep.R before this script.")

message(sprintf("[06] Custom DBs available: %s", paste(DB_NAMES, collapse = ", ")))

## Load each gene-set collection once (named list of symbol-vector sets).
## 04 already applied size filters; we pass minSize=1/maxSize=Inf to run_fgsea
## so those filters are NOT applied a second time (belt-and-suspenders comment
## below in the loop explains why we still pass the config values there).
collections <- lapply(db_registry, function(rds_p) {
  obj <- readRDS(rds_p)
  ## Normalise: both T2G/T2N shape and direct named-list shape are supported
  ## by load_custom_geneset(); here the 04 cache already converted to the
  ## named-list shape (b).  But call it defensively so this script is robust
  ## if the owner ever swaps the cache format.
  if (is.list(obj) && !is.null(names(obj)) &&
      !any(c("T2G", "T2N") %in% names(obj))) {
    ## shape (b): already a named list of character vectors
    lapply(obj, function(g) unique(as.character(g)))
  } else {
    ## shape (a): T2G/T2N toolkit list — normalise to named list
    load_custom_geneset(rds_p)
  }
})
names(collections) <- DB_NAMES

## Report set counts
for (nm in DB_NAMES)
  message(sprintf("[06]   %-25s  %d gene sets loaded", nm, length(collections[[nm]])))

## ----------------------------------------------------------------------------
## Per-(contrast x custom-DB) GSEA — direct clusterProfiler, NOT run_fgsea().
## Mirrors the 14839 run_one_custom_db() idiom: build TERM2GENE (term,gene) +
## TERM2NAME (term,name) from the cached named-list collection, run GSEA(by="fgsea",
## eps=0), then APPLY THE REQUIRED PATCH res@geneSets <- split(T2G$gene, T2G$term)
## so the viz running-sum/gseaplot2 can walk membership for custom sets (clusterProfiler
## does not always populate @geneSets for a supplied TERM2GENE).
## Returns a gseaResult (or NULL on no overlap / error / no result).
## ----------------------------------------------------------------------------
build_term2gene <- function(gs) {
  data.frame(
    term = rep(names(gs), lengths(gs)),
    gene = unlist(gs, use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

run_one_custom_db <- function(ranked, gs, dbn) {
  T2G <- build_term2gene(gs)
  ## TERM2NAME: the set name IS the display name for custom DBs (no separate description
  ## table in the cached named-list shape); pass term==name so clusterProfiler keeps the ID.
  T2N <- data.frame(term = names(gs), name = names(gs), stringsAsFactors = FALSE)
  T2G_f <- dplyr::filter(T2G, gene %in% names(ranked))
  if (nrow(T2G_f) == 0) { message(sprintf("    [%s] no overlapping genes", dbn)); return(NULL) }

  res <- tryCatch({
    set.seed(SEED)
    clusterProfiler::GSEA(
      geneList      = ranked,
      TERM2GENE     = T2G_f,
      TERM2NAME     = T2N,
      minGSSize     = MINSZ, maxGSSize = MAXSZ,
      pvalueCutoff  = PCUT,          # 1 — keep all
      pAdjustMethod = PADJM,         # "fdr"
      eps           = EPS,           # 0 — exact p-values
      by            = "fgsea",
      nPermSimple   = NPERM,
      seed          = TRUE,
      verbose       = FALSE)
  }, error = function(e) { message(sprintf("    [%s] GSEA error: %s", dbn, conditionMessage(e))); NULL })
  if (is.null(res) || nrow(res@result) == 0) return(res)
  ## REQUIRED patch — the running-sum viz (12_gsea_viz.R via enrichplot::gseaplot2) reads
  ## @geneSets; custom TERM2GENE GSEA does not always populate it.
  res@geneSets <- split(T2G_f$gene, T2G_f$term)
  res
}

## ============================================================================
## 3. DE results (the ranking hub)
## ============================================================================

de <- load_de_results()   # validates Symbol rownames + t column; path = DIR_OBJECTS/02_de_results.rds
CONTRASTS <- names(de)    # 7 contrasts: WT_heat, KO_heat, Interaction, ...
message(sprintf("[06] DE results: %d contrasts (%s)", length(CONTRASTS),
                paste(CONTRASTS, collapse = ", ")))

## Advisory: gene universe from 11_emit_universe.R (read for message, not subsetting)
universe_txt <- file.path(DIR_OBJECTS, "gene_universe.txt")
if (file.exists(universe_txt)) {
  universe_genes <- readLines(universe_txt)
  message(sprintf("[06] Gene universe: %d symbols (from gene_universe.txt; advisory only — fgsea intersects internally).",
                  length(universe_genes)))
} else {
  warning("[06] gene_universe.txt not found. Run 11_emit_universe.R for the background universe. Proceeding without it.")
  universe_genes <- character(0)
}

## ============================================================================
## 4. Per-contrast GSEA loop
##    For each contrast: build ranked vector → clusterProfiler::GSEA per DB → collect.
##    Each contrast's full result (named list db -> gseaResult S4) is cached as
##    gsea_custom_<contrast>.rds (load_or_compute) so the viz can drive the toolkit
##    plotters. The tidy master rows are DERIVED from each gseaResult below via
##    normalize_gsea_results(), then projected to the strict 10-col master schema.
## ============================================================================

all_results <- list()   # accumulates TIDY data.frames for every (contrast × DB) pair

MASTER_COLS <- YAML_CONFIG$schemas$master_gsea_table$required_columns

## gseaResult -> tidy master-schema rows (toolkit normalizer, projected to 10 cols).
tidy_from_gsea <- function(g, dbn, co) {
  if (is.null(g) || nrow(g@result) == 0) return(.empty_gsea_df())
  nt <- normalize_gsea_results(g, database = dbn, contrast = co)   # toolkit tibble (NES uppercase)
  if ("NES" %in% colnames(nt) && !"nes" %in% colnames(nt))
    nt <- dplyr::rename(nt, nes = NES)
  ## drop toolkit extras (leading_edge_size, gene_ratio, genes_full_set, neg_log_padj)
  as.data.frame(nt[, MASTER_COLS, drop = FALSE], stringsAsFactors = FALSE)
}

for (co in CONTRASTS) {
  message(sprintf("[06] == contrast: %s ==", co))

  ## Ranked vector by t-statistic (the config rank metric; Symbol -> t, sorted decreasing)
  ranked <- build_ranked_vector(de[[co]], metric = RANK_COL)

  ## Cache key = contrast name. Cached value = named list db -> gseaResult S4 (with @geneSets
  ## patched). setNames(lapply(...)) keeps every DB as a named slot even when a DB returns NULL.
  gsea_list <- load_or_compute(
    sprintf("gsea_custom_%s.rds", co),
    compute_fn = function() {
      setNames(
        lapply(DB_NAMES, function(dbn) {
          gsets <- collections[[dbn]]
          message(sprintf("  [%s] GSEA on %s (%d sets)...", co, dbn, length(gsets)))
          g <- run_one_custom_db(ranked = ranked, gs = gsets, dbn = dbn)
          if (!is.null(g) && nrow(g@result) > 0) {
            n_sig <- sum(g@result$p.adjust < FDR_DISP, na.rm = TRUE)
            message(sprintf("    -> %d pathways tested, %d at padj<%.2g (NES>0: %d, NES<0: %d)",
                            nrow(g@result), n_sig, FDR_DISP,
                            sum(g@result$p.adjust < FDR_DISP & g@result$NES > 0, na.rm = TRUE),
                            sum(g@result$p.adjust < FDR_DISP & g@result$NES < 0, na.rm = TRUE)))
          } else {
            message(sprintf("    -> 0 pathways (no overlap / no result for %s sets)", dbn))
          }
          g
        }),
        DB_NAMES
      )
    },
    force = FALSE,
    desc  = sprintf("custom GSEA for contrast '%s'", co)
  )

  ## Derive tidy rows per DB from the cached gseaResult objects
  tidy_list <- setNames(lapply(DB_NAMES, function(dbn) tidy_from_gsea(gsea_list[[dbn]], dbn, co)), DB_NAMES)

  ## Accumulate for CSV writes (bind all DB results for this contrast)
  contrast_rows <- dplyr::bind_rows(tidy_list)
  if (nrow(contrast_rows) > 0) {
    ## Per-contrast output: one combined CSV (all custom DBs, sorted by padj).
    ## Layout: 03_results/06_gsea/tables/by_contrast/<contrast>/gsea_custom.csv
    ## The `database` column in each row identifies which DB it came from.
    ## Convention (§5 global): contrast is the DIRECTORY, not a filename suffix.
    cdir <- file.path(tbl_dir, "by_contrast", co)
    dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
    out_contrast <- dplyr::arrange(round_numeric_cols(contrast_rows), padj)
    readr::write_csv(out_contrast, file.path(cdir, "gsea_custom.csv"))
    message(sprintf("  -> wrote by_contrast/%s/gsea_custom.csv  (%d rows)", co, nrow(out_contrast)))
  } else {
    message(sprintf("  -> no rows for contrast %s (all DBs returned empty or no overlap)", co))
  }

  all_results[[co]] <- tidy_list   # store per-DB TIDY rows for master-table writes + asserts
}

## ============================================================================
## 5. Overview CSV: _overview/gsea_custom_all.csv — ALL custom rows, ALL contrasts.
##    Same schema; contrast is already a column. Sorted by contrast, then padj.
##    This is the flat "full sweep" file for downstream viz (12_gsea_viz.R).
## ============================================================================

all_rows <- dplyr::bind_rows(
  lapply(all_results, function(res_for_contrast) dplyr::bind_rows(res_for_contrast))
)

if (nrow(all_rows) > 0) {
  overview_rows <- dplyr::arrange(round_numeric_cols(all_rows), contrast, padj)
  readr::write_csv(overview_rows, file.path(overview_dir, "gsea_custom_all.csv"))
  message(sprintf("[06] overview: _overview/gsea_custom_all.csv  (%d rows total)", nrow(overview_rows)))
} else {
  warning("[06] No custom GSEA rows produced at all — check DE results, gene symbols, and collection overlap.")
}

## ============================================================================
## 6. Append to master_gsea_table.csv — PER-DATABASE key.
##
##    append_master_table(df, "master_gsea_table.csv", key_col = "database")
##    drops and replaces every row whose `database` column value matches any
##    value in df$database, then writes the full combined table.
##
##    By calling it ONCE PER DATABASE (not once for all custom rows combined),
##    each DB key ("TransportDB", "MitoPathways", "MitoXplorer", "Lombardi2022_HIF")
##    replaces exactly its own rows and leaves all other rows — including the
##    MSigDB rows written by 05_gsea_msigdb_run.R ("Hallmark", "KEGG", etc.) —
##    completely untouched. This is the idempotency guarantee.
## ============================================================================

for (dbn in DB_NAMES) {
  ## Collect all contrasts' rows for this DB
  db_rows <- dplyr::bind_rows(
    lapply(CONTRASTS, function(co) {
      r <- all_results[[co]][[dbn]]
      if (is.null(r) || nrow(r) == 0L) return(NULL)
      r
    })
  )
  if (is.null(db_rows) || nrow(db_rows) == 0L) {
    message(sprintf("[06] master: %s — no rows to append (skipping)", dbn))
    next
  }
  ## Validate that the `database` column matches the expected tag (sanity check)
  stopifnot(all(db_rows$database == dbn))
  ## append_master_table() validates schema, enforces lowercase `nes`, does the
  ## drop-and-replace for key_col="database", rounds to 9 sig figs, writes CSV.
  append_master_table(db_rows, "master_gsea_table.csv", key_col = "database")
  message(sprintf("[06] master: appended %d rows for '%s' (key_col=database; MSigDB rows preserved).",
                  nrow(db_rows), dbn))
}

## ============================================================================
## 7. Summary statistics (printed to log; no file written — keep it light)
## ============================================================================

message("\n[06] ===== SUMMARY =====")
message(sprintf("  Contrasts   : %d (%s)", length(CONTRASTS), paste(CONTRASTS, collapse = ", ")))
message(sprintf("  Custom DBs  : %d (%s)", length(DB_NAMES),  paste(DB_NAMES,  collapse = ", ")))
if (nrow(all_rows) > 0) {
  message(sprintf("  Total rows  : %d (all contrasts x all custom DBs)", nrow(all_rows)))
  for (dbn in DB_NAMES) {
    db_sub <- dplyr::filter(all_rows, database == dbn)
    n_sig  <- sum(db_sub$padj < FDR_DISP, na.rm = TRUE)
    message(sprintf("    %-25s  %d rows total, %d at padj<%.2g",
                    dbn, nrow(db_sub), n_sig, FDR_DISP))
  }
  message(sprintf("  Sign note: NES > 0 = enriched in the NUMERATOR side of the contrast."))
  message(sprintf("  For WT_heat/KO_heat: NES>0 = enriched at 39C; NES<0 = enriched at 37C."))
  message(sprintf("  For Interaction:     NES>0 = heat response LARGER in WT than cGASKO"))
  message(sprintf("     (the cGAS-dependence payoff; 1 df, lowest power — phrase as"))
  message(sprintf("      'no detectable cGAS-dependence at n=5', NEVER 'cGAS-independent')."))
  message(sprintf("  Lombardi2022_HIF NES on Interaction = orthogonal HIF check vs Biomni's"))
  message(sprintf("    16-gene list; never crown HIF1a/HIF2a; floor claims at L3."))
}
message("[06_gsea_custom_run] COMPLETE.")

## ============================================================================
## 8. Final structural asserts (fail loud if outputs are missing)
## ============================================================================

## Each contrast that produced any rows must have a gsea_custom.csv
for (co in CONTRASTS) {
  co_rows <- dplyr::bind_rows(all_results[[co]])
  cdir    <- file.path(tbl_dir, "by_contrast", co)
  if (nrow(co_rows) > 0) {
    stopifnot(
      "gsea_custom.csv missing for contrast with rows" =
        file.exists(file.path(cdir, "gsea_custom.csv"))
    )
  }
}

## Overview CSV exists if any rows were produced
if (nrow(all_rows) > 0)
  stopifnot(file.exists(file.path(overview_dir, "gsea_custom_all.csv")))

## Each per-contrast RDS cache exists
for (co in CONTRASTS)
  stopifnot(file.exists(file.path(DIR_OBJECTS, sprintf("gsea_custom_%s.rds", co))))

## master_gsea_table.csv has the required schema columns
master_path <- file.path(DIR_MASTER, "master_gsea_table.csv")
if (file.exists(master_path)) {
  master <- readr::read_csv(master_path, show_col_types = FALSE, progress = FALSE)
  req    <- YAML_CONFIG$schemas$master_gsea_table$required_columns
  if (!is.null(req)) {
    missing_cols <- setdiff(req, colnames(master))
    if (length(missing_cols) > 0)
      stop("master_gsea_table.csv missing required schema columns: ",
           paste(missing_cols, collapse = ", "))
  }
  ## Confirm custom DB rows are present (for any DB that produced rows)
  dbs_with_rows <- DB_NAMES[vapply(DB_NAMES, function(dbn) {
    r <- dplyr::bind_rows(lapply(CONTRASTS, function(co) all_results[[co]][[dbn]]))
    nrow(r) > 0
  }, logical(1))]
  if (length(dbs_with_rows) > 0) {
    in_master <- intersect(dbs_with_rows, unique(master$database))
    if (length(in_master) < length(dbs_with_rows)) {
      missing_dbs <- setdiff(dbs_with_rows, in_master)
      warning("Custom DBs with rows not found in master_gsea_table.csv: ",
              paste(missing_dbs, collapse = ", "))
    }
  }
}

message("[06_gsea_custom_run] Structural asserts passed.")
