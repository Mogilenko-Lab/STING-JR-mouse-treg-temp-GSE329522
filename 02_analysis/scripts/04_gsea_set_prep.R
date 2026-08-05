# 04_gsea_set_prep.R — COMPUTE
## Load + cache ALL gene-set collections for the GSEA arm (MSigDB 8 collections +
## custom DBs configured in databases.custom), ready for 05_gsea_msigdb_run.R and
## 06_gsea_custom_run.R. Produces named-list .rds objects (gene symbol -> character
## vector) in 03_results/objects/ plus a manifest CSV enumerating every available
## collection.
##
## Run from project root:
##   Rscript 02_analysis/scripts/04_gsea_set_prep.R
##
## SYMBOL VOCABULARY. Reference collections ship CURRENT MGI symbols; this project's
## matrix is frozen to GENCODE vM25's vintage, where 2,341 of 19,679 symbols are no longer
## current. Every collection is therefore resolved through the committed alias map
## (00_data/references/symbol_alias/symbol_alias_map.csv, built by 00_symbol_alias_map.R)
## BEFORE the size filter — before, because resolving after it would trim
## MITOPATHWAYS_OXPHOS.Complex_V.CV_subunits for the wrong reason. A resolved set carries
## BOTH spellings: the accepted reference symbol is absent from this vocabulary and so
## matches nothing, and the published set_size is fgsea's post-intersection `size` rather
## than the nominal length, so the extra spelling is inert where it counts. Substituting
## instead — which reads cleaner — shrinks the 12 sets that carry both vintages of one gene
## below gsea_min_size and DROPS them, which the fix must not do (see resolve_sets).
##
## The reporting contract is a ledger with a bucket per cause, never a pass/fail recovery
## floor: geneset_symbol_ledger.csv carries matched / matched-via-alias / flagged /
## rejected / absent per SET, and geneset_manifest.csv rolls the same columns up per
## collection. `n_overlap_before_alias` is kept beside `n_matched` so the pre-fix number
## stays legible next to the new one.
##
## Outputs
## -------
##   03_results/objects/geneset_msigdb_<name>.rds   — one per databases.msigdb entry
##   03_results/objects/geneset_custom_<name>.rds   — one per databases.custom entry (if path exists)
##   03_results/objects/geneset_manifest.csv        — database, n_sets, type, source, ledger counts
##   03_results/objects/geneset_symbol_ledger.csv   — the same counts per gene SET
##   03_results/objects/geneset_alias_applied.csv   — every (set, reference->matrix) pair applied
##
## Assumptions
## -----------
##   * config.R sets PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS, SPECIES, GSEA_MIN_SIZE,
##     GSEA_MAX_SIZE, RANK_METRIC, %||%, load_or_compute().
##   * de_gsea_helpers.R (sourced after config.R) provides the path-keyed variant of
##     load_or_compute(), load_msigdb_collection(), load_custom_geneset().
##   * symbol_alias.R provides resolve_sets(), symbol_ledger(), assert_ledger_closes().
##   * databases.msigdb is a list of {category, subcategory, name} maps.
##   * databases.custom is a list of {name, path} maps; paths are relative to
##     PROJECT_ROOT. Project RDS (paths starting 00_data/) are built by the 00-series
##     curation scripts; toolkit RDS (paths under 01_modules/) are built by
##     the RNAseq-toolkit reference-processing pipeline. If any path is absent the
##     script warns and skips that DB — it does NOT error.
##   * msigdbr is a lazy dependency (not sourced at the top). If missing, MSigDB
##     collections are skipped with a clear install message.
##   * Idempotent: load_or_compute() returns from cache on re-runs unless force = TRUE.
##     The alias resolution happens INSIDE the cached closure, so a change to the alias map
##     requires deleting the caches (or FORCE_GENESET_REBUILD=1) to take effect — the
##     script says so on every run.

# ============================================================================
# 0. Environment setup
# ============================================================================

source("02_analysis/config/config.R")        # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS, SPECIES,
                                              # GSEA_MIN_SIZE/MAX_SIZE, load_or_compute, %||%
source("02_analysis/helpers/de_gsea_helpers.R")  # path-keyed load_or_compute, load_msigdb_collection,
                                                  # load_custom_geneset
source("02_analysis/helpers/symbol_alias.R")     # resolve_sets, symbol_ledger, assert_ledger_closes

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})
options(stringsAsFactors = FALSE)

# Validate guard: require key config blocks populated
if (is.null(YAML_CONFIG$databases$msigdb) || length(YAML_CONFIG$databases$msigdb) == 0) {
  stop("databases.msigdb is empty in analysis_config.yaml. Check config before running 04.")
}

min_size <- GSEA_MIN_SIZE  # from config.R (15)
max_size <- GSEA_MAX_SIZE  # from config.R (500)
species  <- SPECIES        # "Mus musculus"
FORCE    <- nzchar(Sys.getenv("FORCE_GENESET_REBUILD"))

message(sprintf("[04_gsea_set_prep] species=%s  min_size=%d  max_size=%d  force=%s",
                species, min_size, max_size, FORCE))

# ============================================================================
# 0b. The alias map + the two vocabulary layers
# ============================================================================
## MATRIX_VOCABULARY is gene_universe.txt, the modelled DE background. In this project it
## is ALSO exactly the set of symbols the delivered CPM table carries — the duplicate-symbol
## collapse at 01_mapping_qc.R drops Ensembl ids, never symbols, so 19,679 unique gene_name
## values survive as 19,679 modelled symbols. There is therefore exactly ONE vocabulary
## layer here and n_expression_filtered is structurally 0; the column is kept so the ledger
## columns match the human compartment's, and the two layers are compared below so a future
## divergence surfaces rather than hides.
##
## What this project CANNOT say: n_absent_from_reference. That needs the GENCODE vM25
## feature list, and the collaborators delivered a CPM table rather than a quantification
## against a tracked GTF, so `reference_vocabulary_available` is FALSE and an unmatched
## reference gene lands in n_below_detection — "not in the delivered quantification" — which
## is the strongest honest statement available. It is NOT "absent from the reference".

SA <- YAML_CONFIG$symbol_alias
if (is.null(SA))
  stop("[04] analysis_config.yaml has no `symbol_alias:` block — see 00_symbol_alias_map.R.")
VOCAB_TXT <- file.path(PROJECT_ROOT, SA$matrix_vocabulary)
MAP_CSV   <- file.path(PROJECT_ROOT, SA$map_path)
if (!file.exists(VOCAB_TXT))
  stop("[04] matrix vocabulary not found at ", VOCAB_TXT, " — run 11_emit_universe.R first.")
if (!file.exists(MAP_CSV))
  stop("[04] symbol alias map not found at ", MAP_CSV,
       " — run 02_analysis/scripts/00_symbol_alias_map.R first.")

MATRIX_VOCABULARY <- unique(trimws(readLines(VOCAB_TXT, warn = FALSE)))
MATRIX_VOCABULARY <- MATRIX_VOCABULARY[nzchar(MATRIX_VOCABULARY)]
ALIAS_MAP <- readr::read_csv(MAP_CSV, show_col_types = FALSE, progress = FALSE)
PAIRS <- accepted_pairs(ALIAS_MAP)

QUANT_CSV <- file.path(PROJECT_ROOT, SA$quantified_vocabulary %||% "")
QUANTIFIED_VOCABULARY <- MATRIX_VOCABULARY
if (nzchar(SA$quantified_vocabulary %||% "") && file.exists(QUANT_CSV)) {
  qcol <- SA$quantified_vocabulary_column %||% "gene_name"
  q <- readr::read_csv(QUANT_CSV, col_select = dplyr::all_of(qcol),
                       show_col_types = FALSE, progress = FALSE)
  QUANTIFIED_VOCABULARY <- unique(as.character(q[[qcol]]))
  QUANTIFIED_VOCABULARY <- QUANTIFIED_VOCABULARY[!is.na(QUANTIFIED_VOCABULARY) &
                                                   nzchar(QUANTIFIED_VOCABULARY)]
} else {
  message("[04] delivered CPM table absent; the quantified vocabulary falls back to the ",
          "modelled universe, so n_expression_filtered cannot be non-zero by construction.")
}
message(sprintf(
  "[04] vocabulary: %d modelled symbols; %d quantified in the delivered table (%d quantified but not modelled). Alias map: %d accepted pairs, %d flagged, %d rejected.",
  length(MATRIX_VOCABULARY), length(QUANTIFIED_VOCABULARY),
  length(setdiff(QUANTIFIED_VOCABULARY, MATRIX_VOCABULARY)), length(PAIRS),
  sum(ALIAS_MAP$resolution == "flagged_for_review"),
  sum(!ALIAS_MAP$resolution %in% c("accepted", "flagged_for_review"))))
if (!FORCE)
  message("[04] cached collections are reused as-is. After editing the alias map, re-run ",
          "with FORCE_GENESET_REBUILD=1 or the resolution will not reach the caches.")

# ============================================================================
# 1. Helpers
# ============================================================================

## Size-filter a named-list gene-set collection (fgsea pathways format).
## Removes sets whose symbol count is outside [min, max]. Returns the filtered list.
filter_by_size <- function(gsets, min_sz, max_sz) {
  sizes <- vapply(gsets, length, integer(1))
  gsets[sizes >= min_sz & sizes <= max_sz]
}

## Resolve one raw collection into this matrix's vocabulary, then size-filter. Emits the
## per-set ledger and the applied-pair rows as a side effect into LEDGERS/APPLIED, so the
## published ledger covers every set the reference ships — including the ones the size
## filter then drops, which is where a vocabulary loss would otherwise become invisible.
LEDGERS <- list(); APPLIED <- list(); COLLAPSES <- list(); MANY_TO_ONE <- list()

resolve_and_filter <- function(raw, database) {
  res <- resolve_sets(raw, MATRIX_VOCABULARY, ALIAS_MAP)
  led <- symbol_ledger(raw, ALIAS_MAP,
                       ranked_vocabulary = MATRIX_VOCABULARY,
                       matrix_vocabulary = QUANTIFIED_VOCABULARY,
                       reference_vocabulary = NULL)
  assert_ledger_closes(led, sprintf("geneset ledger for %s", database))
  # Invariant 1: exact-string matching may not move. The alias fix may only ADD.
  n_overlap_before <- vapply(raw, function(g)
    length(intersect(unique(g), MATRIX_VOCABULARY)), integer(1))
  stopifnot("alias resolution changed an exact-match count" =
              all(led$n_matched == n_overlap_before[led$gene_set]))
  # Invariant 3: a set may not gain more matching genes than pairs applied, may not shed a
  # nominal member, and may not acquire a duplicate the reference did not already ship.
  n_after <- vapply(res$sets, function(g)
    length(intersect(unique(g), MATRIX_VOCABULARY)), integer(1))
  stopifnot(
    "a set gained more matching genes than pairs applied" =
      all(n_after[led$gene_set] - led$n_matched <= led$n_matched_via_alias),
    "resolution shortened a set's nominal length" =
      all(lengths(res$sets) >= lengths(raw)[names(res$sets)]),
    "resolution introduced a duplicate the reference did not ship" =
      all(lengths(res$sets) - vapply(res$sets, function(s) length(unique(s)), integer(1)) <=
            lengths(raw)[names(res$sets)] -
              vapply(raw, function(s) length(unique(s)), integer(1))[names(res$sets)]))

  led <- led %>% mutate(database = database,
                        n_overlap_before_alias = as.integer(n_overlap_before[.data$gene_set]),
                        .before = 1)
  LEDGERS[[database]] <<- led
  if (nrow(res$applied))
    APPLIED[[database]] <<- res$applied %>% mutate(database = database, .before = 1)
  if (nrow(res$collapsed))
    COLLAPSES[[database]] <<- res$collapsed %>% mutate(database = database, .before = 1)
  if (nrow(res$many_to_one))
    MANY_TO_ONE[[database]] <<- res$many_to_one %>% mutate(database = database, .before = 1)

  # Nominal-filter membership is the OLD membership plus whatever resolution newly admits,
  # by construction, so the fix cannot cost a set. It has to be spelled out rather than left
  # to filter_by_size(): some msigdbr sets ship a symbol twice (one gs_name, two Ensembl
  # ids), filter_by_size() counts the raw length while a resolved set is de-duplicated, and
  # REACTOME_ENDOSOMAL_VACUOLAR_PATHWAY at raw length 15 / 14 distinct symbols would
  # otherwise fall out of a 15-floor for a reason that has nothing to do with aliases. The
  # de-duplication is reported below rather than absorbed.
  before  <- filter_by_size(raw, min_size, max_size)
  admits  <- filter_by_size(res$sets, min_size, max_size)
  keep    <- union(names(before), names(admits))
  crossed <- setdiff(names(admits), names(before))
  dedup   <- names(before)[!names(before) %in% names(admits)]
  filt    <- res$sets[keep]
  stopifnot("the nominal size filter lost a set" = all(names(before) %in% names(filt)))
  message(sprintf("  [alias] %s: %d sets, %d pairs applied over %d sets; %d newly admitted at the nominal size filter, %d retained only because de-duplication may not cost a set",
                  database, length(raw), nrow(res$applied),
                  if (nrow(res$applied)) dplyr::n_distinct(res$applied$gene_set) else 0L,
                  length(crossed), length(dedup)))
  if (length(crossed))
    message("    [alias] newly admitted: ", paste(crossed, collapse = ", "))
  if (length(dedup))
    message("    [alias] duplicate-symbol sets retained: ", paste(dedup, collapse = ", "))
  filt
}

## Build a one-row manifest tibble for a single collection. The six ledger counts are the
## per-collection roll-up of geneset_symbol_ledger.csv, summed over the collection's UNIQUE
## reference symbols rather than over sets, so a gene in twelve sets is counted once.
manifest_row <- function(database, type, n_sets, source, cache_path, raw = NULL) {
  base <- tibble::tibble(
    database     = as.character(database),
    type         = as.character(type),
    n_sets       = as.integer(n_sets),
    source       = as.character(source),
    cache_path   = as.character(cache_path)
  )
  if (is.null(raw)) return(base)
  u <- list(all = unique(unlist(raw, use.names = FALSE)))
  led <- symbol_ledger(u, ALIAS_MAP,
                       ranked_vocabulary = MATRIX_VOCABULARY,
                       matrix_vocabulary = QUANTIFIED_VOCABULARY,
                       reference_vocabulary = NULL)
  assert_ledger_closes(led, sprintf("manifest roll-up for %s", database))
  dplyr::bind_cols(base, led %>% dplyr::select(
    n_unique_set_genes, n_matched, n_matched_via_alias, n_alias_flagged_for_review,
    n_alias_rejected_ambiguous, n_expression_filtered, n_below_detection,
    n_absent_from_reference, reference_vocabulary_available))
}

# ============================================================================
# 2. MSigDB collections (8 configured)
# ============================================================================
## Each collection is cached as a NAMED LIST of mouse gene-symbol vectors, the
## fgsea `pathways` format expected by run_fgsea() in 05 and 06.
## Size-filtering (min_size/max_size) is applied AFTER fetching so that:
##   (a) the cached object only contains sets actually usable at the configured thresholds;
##   (b) downstream scripts pass the cached list directly to fgsea with minSize=1/maxSize=Inf
##       (the thresholds have already been applied here).

msigdb_manifest <- list()
RAW <- list()   # the unresolved collections, kept for the manifest roll-up

for (m in YAML_CONFIG$databases$msigdb) {
  nm    <- m$name
  cat_  <- m$category
  sub_  <- m$subcategory %||% ""

  cache_file <- sprintf("geneset_msigdb_%s.rds", nm)
  cache_path <- file.path(DIR_OBJECTS, cache_file)

  fetch_raw <- function() {
    message(sprintf("  [msigdb] fetching %s (category=%s subcategory='%s')", nm, cat_, sub_))
    load_msigdb_collection(category = cat_, subcategory = sub_, species = species)
  }
  gsets <- load_or_compute(cache_path, function() {
    raw <- fetch_raw()
    if (length(raw) == 0) {
      warning(sprintf("  [msigdb] %s returned 0 gene sets — stored as empty list.", nm))
      return(list())
    }
    RAW[[nm]] <<- raw
    filt <- resolve_and_filter(raw, nm)
    message(sprintf("  [msigdb] %s: %d sets raw -> %d kept (size %d-%d)",
                    nm, length(raw), length(filt), min_size, max_size))
    filt
  }, force = FORCE)
  # The ledger must be published even when the cache short-circuited the closure, or a
  # cached run would silently drop the accounting it exists to provide.
  if (is.null(RAW[[nm]])) {
    raw <- fetch_raw()
    if (length(raw)) { RAW[[nm]] <- raw; invisible(resolve_and_filter(raw, nm)) }
  }

  n <- length(gsets)
  src <- if (nzchar(sub_)) sprintf("MSigDB %s/%s", cat_, sub_) else sprintf("MSigDB %s", cat_)
  msigdb_manifest[[nm]] <- manifest_row(nm, "msigdb", n, src, cache_path, raw = RAW[[nm]])
  message(sprintf("[04] MSigDB %-15s -> %d sets  (%s)", nm, n, cache_file))
}

# ============================================================================
# 3. Custom gene-set databases
# ============================================================================
## Each custom DB is cached as a NAMED LIST of mouse gene-symbol vectors (same
## fgsea pathways format). load_custom_geneset() handles both toolkit T2G/T2N
## shape and direct named-list shape (including the provenance-stamped project rds).
## Size-filtering is applied at the same thresholds as MSigDB.
##
## MISSING PATH HANDLING:
##   If the configured path does not exist (relative to PROJECT_ROOT), the script
##   emits a WARNING message (including the expected absolute path) and skips that
##   DB entirely. The manifest records the entry with n_sets = NA and status "MISSING".
##   The script continues with the remaining DBs and does NOT stop.

custom_manifest <- list()

for (db in YAML_CONFIG$databases$custom) {
  nm   <- db$name
  rel  <- db$path   # relative to PROJECT_ROOT as declared in config

  # Resolve to absolute path (config paths are project-root-relative)
  abs_path <- file.path(PROJECT_ROOT, rel)

  cache_file <- sprintf("geneset_custom_%s.rds", nm)
  cache_path <- file.path(DIR_OBJECTS, cache_file)

  if (!file.exists(abs_path)) {
    ## Warn clearly; record as MISSING in manifest; continue.
    warning(sprintf(
      "[04_gsea_set_prep] MISSING custom DB '%s': file not found at\n    %s\n  (config path: %s)\n  Skipping this DB. Run its builder script first (e.g. 00b or RNAseq-toolkit reference pipeline).",
      nm, abs_path, rel
    ))
    custom_manifest[[nm]] <- manifest_row(
      nm, "custom_MISSING", NA_integer_, rel, NA_character_
    )
    next  # skip to next DB — do NOT error
  }

  fetch_raw <- function() {
    message(sprintf("  [custom] loading %s from %s", nm, abs_path))
    load_custom_geneset(abs_path)
  }
  gsets <- load_or_compute(cache_path, function() {
    raw <- fetch_raw()
    if (length(raw) == 0) {
      warning(sprintf("  [custom] %s returned 0 gene sets — stored as empty list.", nm))
      return(list())
    }
    RAW[[nm]] <<- raw
    filt <- resolve_and_filter(raw, nm)
    message(sprintf("  [custom] %s: %d sets raw -> %d kept (size %d-%d)",
                    nm, length(raw), length(filt), min_size, max_size))
    filt
  }, force = FORCE)
  if (is.null(RAW[[nm]])) {
    raw <- fetch_raw()
    if (length(raw)) { RAW[[nm]] <- raw; invisible(resolve_and_filter(raw, nm)) }
  }

  n <- length(gsets)
  custom_manifest[[nm]] <- manifest_row(nm, "custom", n, rel, cache_path, raw = RAW[[nm]])
  message(sprintf("[04] custom  %-15s -> %d sets  (%s)", nm, n, cache_file))
}

# ============================================================================
# 4. Write manifest CSV
# ============================================================================
## geneset_manifest.csv is the authoritative index of all available gene-set
## collections for 05_gsea_msigdb_run.R, 06_gsea_custom_run.R, and the viz arm.
## Columns: database, type, n_sets, source, cache_path + the vocabulary ledger roll-up.
## type "msigdb"         -> load the cache_path rds directly (already size-filtered).
## type "custom"         -> load the cache_path rds directly (already size-filtered).
## type "custom_MISSING" -> n_sets is NA; the build for that DB has not run yet.
##
## How to read the ledger columns. They count the collection's UNIQUE reference symbols and
## close on
##   n_unique_set_genes == n_matched + n_matched_via_alias + n_alias_flagged_for_review +
##                         n_alias_rejected_ambiguous + n_expression_filtered +
##                         n_below_detection + n_absent_from_reference
## n_matched is the exact string match — the pre-fix number, and it must never move.
## n_matched_via_alias is what the fix recovered. n_alias_rejected_ambiguous is the
## ownership guard firing (the candidate is another gene's official symbol; Gck->Gk,
## Tacr1->Spr, Slc6a4->Htt are the kind of trap it catches). n_expression_filtered is
## structurally 0 in this project — the delivered CPM table and the modelled universe are
## the same 19,679 symbols. n_below_detection means "not in the delivered quantification",
## and n_absent_from_reference stays 0 with reference_vocabulary_available = FALSE because
## the vM25 feature list is not persisted here; that split is the one thing this compartment
## cannot make, and it is reported as unavailable rather than guessed.

manifest_all <- dplyr::bind_rows(c(msigdb_manifest, custom_manifest))
manifest_path <- file.path(DIR_OBJECTS, "geneset_manifest.csv")
readr::write_csv(manifest_all, manifest_path)
message(sprintf("[04] manifest written -> %s  (%d collections)", manifest_path, nrow(manifest_all)))

# ---- the per-SET ledger, and the applied pairs ------------------------------------------
## Set granularity is mandatory rather than nice: every claim in this family is per-set
## (MITOPATHWAYS_OXPHOS.Complex_V, HSR_core, HALLMARK_UNFOLDED_PROTEIN_RESPONSE), and a
## per-collection roll-up cannot carry them. Every set the reference ships gets a row,
## including the ones the size filter drops — that is exactly where a vocabulary loss
## would otherwise become invisible.

ledger_all <- dplyr::bind_rows(LEDGERS)
if (nrow(ledger_all)) {
  assert_ledger_closes(ledger_all, "geneset_symbol_ledger")
  ledger_path <- file.path(DIR_OBJECTS, "geneset_symbol_ledger.csv")
  readr::write_csv(ledger_all %>% arrange(.data$database, .data$gene_set), ledger_path)
  message(sprintf("[04] set-level ledger written -> %s  (%d sets over %d collections)",
                  ledger_path, nrow(ledger_all), dplyr::n_distinct(ledger_all$database)))

  applied_all <- dplyr::bind_rows(APPLIED)
  readr::write_csv(applied_all, file.path(DIR_OBJECTS, "geneset_alias_applied.csv"))
  gained <- ledger_all %>% filter(.data$n_matched_via_alias > 0) %>%
    mutate(pct = 100 * .data$n_matched_via_alias / pmax(.data$n_matched, 1L)) %>%
    arrange(desc(.data$pct))
  message(sprintf("[04] %d sets gained >=1 gene; %d gained >=10%% of their matched size.",
                  nrow(gained), sum(gained$pct >= 10)))
  print(as.data.frame(gained %>% head(12) %>%
                        transmute(.data$database, .data$gene_set,
                                  matched = .data$n_matched,
                                  via_alias = .data$n_matched_via_alias,
                                  resolved = .data$set_size_resolved,
                                  curated = .data$n_unique_set_genes)), row.names = FALSE)

  # Sets that cross the min_size floor because of the fix are the headline: before it they
  # had no row in master_gsea_table.csv at all, in any contrast.
  crossed <- ledger_all %>%
    filter(.data$n_matched < min_size, .data$set_size_resolved >= min_size)
  message(sprintf("[04] %d sets cross the gsea_min_size=%d floor and become testable:",
                  nrow(crossed), min_size))
  if (nrow(crossed))
    print(as.data.frame(crossed %>% transmute(.data$database, .data$gene_set,
                                              matched = .data$n_matched,
                                              resolved = .data$set_size_resolved)),
          row.names = FALSE)

  collapse_all <- dplyr::bind_rows(COLLAPSES)
  if (nrow(collapse_all)) {
    message("[04] sets where a substitution collapsed onto a symbol already present:")
    print(as.data.frame(collapse_all), row.names = FALSE)
  }
  m2o <- dplyr::bind_rows(MANY_TO_ONE)
  if (nrow(m2o)) {
    message("[04] sets where two reference symbols resolved onto ONE matrix symbol ",
            "(two members become one measurement — reported, never silently merged):")
    print(as.data.frame(m2o), row.names = FALSE)
  }
}

# ============================================================================
# 5. Summary
# ============================================================================

n_msigdb  <- sum(manifest_all$type == "msigdb")
n_custom  <- sum(manifest_all$type == "custom")
n_missing <- sum(manifest_all$type == "custom_MISSING")

message(sprintf(
  "[04_gsea_set_prep] COMPLETE — %d MSigDB collections, %d custom DBs cached, %d custom MISSING.",
  n_msigdb, n_custom, n_missing
))
if (n_missing > 0) {
  missing_names <- manifest_all$database[manifest_all$type == "custom_MISSING"]
  message(sprintf("  Missing DBs (check path + builder scripts): %s",
                  paste(missing_names, collapse = ", ")))
}
