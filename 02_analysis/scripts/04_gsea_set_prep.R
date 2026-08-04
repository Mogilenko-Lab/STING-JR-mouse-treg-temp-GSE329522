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
## Outputs
## -------
##   03_results/objects/geneset_msigdb_<name>.rds   — one per databases.msigdb entry
##   03_results/objects/geneset_custom_<name>.rds   — one per databases.custom entry (if path exists)
##   03_results/objects/geneset_manifest.csv        — database, n_sets, type, source
##
## Assumptions
## -----------
##   * config.R sets PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS, SPECIES, GSEA_MIN_SIZE,
##     GSEA_MAX_SIZE, RANK_METRIC, %||%, load_or_compute().
##   * de_gsea_helpers.R (sourced after config.R) provides the path-keyed variant of
##     load_or_compute(), load_msigdb_collection(), load_custom_geneset().
##   * databases.msigdb is a list of {category, subcategory, name} maps.
##   * databases.custom is a list of {name, path} maps; paths are relative to
##     PROJECT_ROOT. Project RDS (paths starting 00_data/) are built by the 00-series
##     curation scripts; toolkit RDS (paths under 01_modules/) are built by
##     the RNAseq-toolkit reference-processing pipeline. If any path is absent the
##     script warns and skips that DB — it does NOT error.
##   * msigdbr is a lazy dependency (not sourced at the top). If missing, MSigDB
##     collections are skipped with a clear install message.
##   * Idempotent: load_or_compute() returns from cache on re-runs unless force = TRUE.

# ============================================================================
# 0. Environment setup
# ============================================================================

source("02_analysis/config/config.R")        # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS, SPECIES,
                                              # GSEA_MIN_SIZE/MAX_SIZE, load_or_compute, %||%
source("02_analysis/helpers/de_gsea_helpers.R")  # path-keyed load_or_compute, load_msigdb_collection,
                                                  # load_custom_geneset

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

message(sprintf("[04_gsea_set_prep] species=%s  min_size=%d  max_size=%d",
                species, min_size, max_size))

# ============================================================================
# 1. Helpers
# ============================================================================

## Size-filter a named-list gene-set collection (fgsea pathways format).
## Removes sets whose symbol count is outside [min, max]. Returns the filtered list.
filter_by_size <- function(gsets, min_sz, max_sz) {
  sizes <- vapply(gsets, length, integer(1))
  gsets[sizes >= min_sz & sizes <= max_sz]
}

## Build a one-row manifest tibble for a single collection.
manifest_row <- function(database, type, n_sets, source, cache_path) {
  tibble::tibble(
    database     = as.character(database),
    type         = as.character(type),
    n_sets       = as.integer(n_sets),
    source       = as.character(source),
    cache_path   = as.character(cache_path)
  )
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

for (m in YAML_CONFIG$databases$msigdb) {
  nm    <- m$name
  cat_  <- m$category
  sub_  <- m$subcategory %||% ""

  cache_file <- sprintf("geneset_msigdb_%s.rds", nm)
  cache_path <- file.path(DIR_OBJECTS, cache_file)

  gsets <- load_or_compute(cache_path, function() {
    message(sprintf("  [msigdb] fetching %s (category=%s subcategory='%s')", nm, cat_, sub_))
    raw <- load_msigdb_collection(category = cat_, subcategory = sub_, species = species)
    if (length(raw) == 0) {
      warning(sprintf("  [msigdb] %s returned 0 gene sets — stored as empty list.", nm))
      return(list())
    }
    filt <- filter_by_size(raw, min_size, max_size)
    message(sprintf("  [msigdb] %s: %d sets raw -> %d kept (size %d-%d)",
                    nm, length(raw), length(filt), min_size, max_size))
    filt
  })

  n <- length(gsets)
  src <- if (nzchar(sub_)) sprintf("MSigDB %s/%s", cat_, sub_) else sprintf("MSigDB %s", cat_)
  msigdb_manifest[[nm]] <- manifest_row(nm, "msigdb", n, src, cache_path)
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

  gsets <- load_or_compute(cache_path, function() {
    message(sprintf("  [custom] loading %s from %s", nm, abs_path))
    raw <- load_custom_geneset(abs_path)
    if (length(raw) == 0) {
      warning(sprintf("  [custom] %s returned 0 gene sets — stored as empty list.", nm))
      return(list())
    }
    filt <- filter_by_size(raw, min_size, max_size)
    message(sprintf("  [custom] %s: %d sets raw -> %d kept (size %d-%d)",
                    nm, length(raw), length(filt), min_size, max_size))
    filt
  })

  n <- length(gsets)
  custom_manifest[[nm]] <- manifest_row(nm, "custom", n, rel, cache_path)
  message(sprintf("[04] custom  %-15s -> %d sets  (%s)", nm, n, cache_file))
}

# ============================================================================
# 4. Write manifest CSV
# ============================================================================
## geneset_manifest.csv is the authoritative index of all available gene-set
## collections for 05_gsea_msigdb_run.R, 06_gsea_custom_run.R, and the viz arm.
## Columns: database, type, n_sets, source, cache_path.
## type "msigdb"         -> load the cache_path rds directly (already size-filtered).
## type "custom"         -> load the cache_path rds directly (already size-filtered).
## type "custom_MISSING" -> n_sets is NA; the build for that DB has not run yet.

manifest_all <- dplyr::bind_rows(c(msigdb_manifest, custom_manifest))
manifest_path <- file.path(DIR_OBJECTS, "geneset_manifest.csv")
readr::write_csv(manifest_all, manifest_path)
message(sprintf("[04] manifest written -> %s  (%d collections)", manifest_path, nrow(manifest_all)))

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
