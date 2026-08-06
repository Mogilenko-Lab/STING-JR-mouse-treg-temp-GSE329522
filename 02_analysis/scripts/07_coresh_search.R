# 07_coresh_search.R — COMPUTE
## CoReSh (COregulation REgulation SearcH) — rank the public MOUSE (mmu) GEO compendium by
## how strongly the configured matrix of project signatures is CO-REGULATED in each public
## dataset (PCA-inspired pctVar score), then derive coregulation gene sets from the top hits
## for downstream GSEA (08_coresh_derived_gsea.R).
##
## Per the `coresh-signature-search` skill
## (01_modules/SciAgent-toolkit/skills/coresh-signature-search/SKILL.md):
##   * search engine  : coresh_batch()        (thin wrapper around fgsea::geseca; coreshMatch())
##   * derived sets   : build_coresh_gmt()     (CoReSh-to-GSEA bridge; gene loadings on top hits)
##   * ID conversion  : sym2ent()/ent2sym(species = "mouse")   (chunk rownames are integer mouse Entrez)
##   * data dependency: preprocessed_chunks/mmu/*_full_objects.qs2   (~20 GB Synapse compendium, syn66227307)
## Skill scripts are sourced from the canonical skill dir; 01_modules/.ref carries the 14839 layout.
##
## Run from project root:
##   Rscript 02_analysis/scripts/07_coresh_search.R
##
## ============================================================================
## !!! DATA DEPENDENCY — READ BEFORE RUNNING !!!
## ----------------------------------------------------------------------------
## (A) THE ~20 GB mmu COMPENDIUM.
##     The mouse compendium (Synapse syn66227307, ~85 *_full_objects.qs2 chunks) is consumed
##     read-only from the shared reference cache, resolved from CORESH_CHUNKS or the chunk-dir
##     convention below. Acquiring it needs a personal Synapse access token and happens
##     out-of-band. Absent chunks stop the script loudly.
##
## (B) QUERY CONTRACT.
##     Queries are coresh.query_signatures in analysis_config.yaml. Each sets entry resolves
##     sets[[contrast]][[direction]][[gate]] from the shared
##     03_results/objects/17_signature_sets.rds object.
##
## (C) CHUNK PATH ASSUMED.
##     Default chunk dir = paths.coresh_chunks where set, else the 14839 convention
##     00_data/references/coresh/current/preprocessed_chunks/mmu, overridable via the
##     CORESH_CHUNKS env var per the skill's path convention. Owner: confirm in-container.
##
## ============================================================================
##
## Outputs (compute-only; lazy heavy deps; idempotent via load_or_compute)
## ----------------------------------------------------------------------------
##   03_results/objects/coresh_ranked.rds              — full per-query ranking (data.table; cache)
##   03_results/objects/coresh_derived_sets.rds        — CoReSh-derived gene sets (named list, fgsea
##                                                        pathways format) for 08_coresh_derived_gsea.R
##   03_results/objects/coresh_query_entrez.rds        — the exact queries searched with (provenance)
##   03_results/08_coresh/tables/coresh_ranked.csv     — the ranked-compendium table
##   03_results/08_coresh/tables/coresh_provenance.csv — derived-set -> source GSE / pctVar trace

# ============================================================================
# 0. Environment setup  (config.R FIRST, then de_gsea_helpers.R — per contract)
# ============================================================================

source("02_analysis/config/config.R")            # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS, SPECIES,
                                                  # stage_dir(), load_or_compute (config.R variant), %||%,
                                                  # project config constants
source("02_analysis/helpers/de_gsea_helpers.R")  # path-keyed load_or_compute, round_numeric_cols
                                                  # (build_ranked_vector / gene_universe are GSEA-arm
                                                  #  helpers; CoReSh queries are integer-Entrez sets,
                                                  #  not ranked vectors, so they are not used here)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})
options(stringsAsFactors = FALSE)

STAGE   <- "08_coresh"                            # declared in analysis_config.yaml:stages
tbl_dir <- stage_dir(STAGE, "tables")            # 03_results/08_coresh/tables/ (created)

# ----------------------------------------------------------------------------
# 0a. CoReSh config (DEFENSIVE — coresh: block may be absent; see assumption (B))
# ----------------------------------------------------------------------------
coresh_cfg <- YAML_CONFIG$coresh %||% list()
CORESH_SPECIES <- coresh_cfg$species        %||% "mouse"   # sym2ent/ent2sym species; mmu chunks => mouse
TOP_N          <- as.integer(coresh_cfg$top_n_hits      %||% 5L)   # top datasets/query -> derived sets
MIN_Q          <- as.integer(coresh_cfg$min_query_size  %||% 3L)   # coresh_batch asserts length(q) >= 3
N_CORES        <- as.integer(coresh_cfg$n_cores         %||% 4L)
USE_PVALUES    <- isTRUE(coresh_cfg$pvalues)                       # variance-only (FALSE) by default
MIN_SET        <- as.integer(GSEA_MIN_SIZE)                        # derived-set size floor (15)
MAX_SET        <- as.integer(GSEA_MAX_SIZE)                        # derived-set size ceiling (500)
N_DERIVE       <- as.integer(coresh_cfg$n_derive        %||% 50L)  # genes/derived set BEFORE size filter
JACCARD_THRESH <- as.numeric(coresh_cfg$jaccard         %||% 0.8)  # dedupe near-identical derived sets

if (!identical(CORESH_SPECIES, "mouse"))
  stop("07_coresh_search: CORESH_SPECIES must be 'mouse' (mmu chunks demand mouse Entrez); got '",
       CORESH_SPECIES, "'. Fix coresh.species in analysis_config.yaml.")

message(sprintf(
  "[07_coresh_search] species=%s top_n_hits=%d min_query_size=%d n_cores=%d pvalues=%s | derived: n=%d size=[%d,%d] jaccard=%.2f",
  CORESH_SPECIES, TOP_N, MIN_Q, N_CORES, USE_PVALUES, N_DERIVE, MIN_SET, MAX_SET, JACCARD_THRESH))

# ============================================================================
# 1. Source the canonical CoReSh skill scripts (skill API; sym2ent FIRST)
# ============================================================================
## Order matters: extract_gene_loadings.R stop()s unless sym2ent/ent2sym are already defined.

CORESH_LIB <- file.path(PROJECT_ROOT,
  "01_modules/SciAgent-toolkit/skills/coresh-signature-search/scripts")
if (!dir.exists(CORESH_LIB))
  stop("07_coresh_search: CoReSh skill scripts dir not found: ", CORESH_LIB,
       " — the coresh-signature-search skill must be present (01_modules/SciAgent-toolkit).")

source(file.path(CORESH_LIB, "symbols_to_entrez.R"))     # sym2ent(), ent2sym()          — FIRST
source(file.path(CORESH_LIB, "coresh_batch.R"))          # coresh_batch(), coreshMatch()
source(file.path(CORESH_LIB, "extract_gene_loadings.R")) # build_coresh_gmt()             — needs sym2ent

# ============================================================================
# 2. Pre-flight: GUARD the data dependency (STOP loudly; NEVER fabricate)
# ============================================================================
## Resolve the mmu chunk dir from (priority): CORESH_CHUNKS env var > paths.coresh_chunks
## > the 14839 convention. If the resolved dir is absent or contains no *_full_objects.qs2,
## STOP with the exact provisioning pointer. See assumption (A)/(C) above.

DEFAULT_CHUNKS_REL <- "00_data/references/coresh/current/preprocessed_chunks/mmu"
chunk_dir <- Sys.getenv("CORESH_CHUNKS", unset = "")
if (!nzchar(chunk_dir)) {
  rel <- YAML_CONFIG$paths$coresh_chunks %||% DEFAULT_CHUNKS_REL
  chunk_dir <- if (startsWith(rel, "/")) rel else file.path(PROJECT_ROOT, rel)
}

PROVISION_NOTE <- "docs/_internal/reasoning/2026-06-24_02_coresh-provisioning.md"
if (!dir.exists(chunk_dir)) {
  stop(
    "CoReSh mmu compendium NOT provisioned — chunk dir absent:\n    ", chunk_dir,
    "\n  The ~20 GB mouse compendium (Synapse syn66227307) must be downloaded IN-CONTAINER by",
    "\n  the owner (needs a personal Synapse access token; Claude cannot do this).",
    "\n  See: ", PROVISION_NOTE,
    "\n  and: 01_modules/SciAgent-toolkit/skills/coresh-signature-search/references/synapse-data-setup.md",
    "\n  Then set paths.coresh_chunks (or the CORESH_CHUNKS env var) to the mmu/ dir and re-run."
  )
}
chunk_files <- list.files(chunk_dir, pattern = "_full_objects\\.qs2$", full.names = TRUE)
if (length(chunk_files) == 0L &&
    dir.exists(file.path(chunk_dir, "mmu")) &&
    length(list.files(file.path(chunk_dir, "mmu"), pattern = "_full_objects\\.qs2$", full.names = TRUE)) > 0L) {
  chunk_dir <- file.path(chunk_dir, "mmu")
  chunk_files <- list.files(chunk_dir, pattern = "_full_objects\\.qs2$", full.names = TRUE)
}
if (length(chunk_files) == 0L) {
  stop(
    "CoReSh mmu compendium NOT provisioned — no *_full_objects.qs2 chunks in:\n    ", chunk_dir,
    "\n  Directory exists but is empty / mis-shaped. Expect ~85 mmu chunks.",
    "\n  See: ", PROVISION_NOTE,
    "\n  and: 01_modules/SciAgent-toolkit/skills/coresh-signature-search/references/synapse-data-setup.md"
  )
}
message(sprintf("[07_coresh_search] CoReSh mmu compendium: %d chunk(s) at %s",
                length(chunk_files), chunk_dir))

# ============================================================================
# 3. Build one signature query per configured matrix entry
# ============================================================================

qsigs <- coresh_cfg$query_signatures %||% list()
if (is.null(qsigs$object) || !nzchar(as.character(qsigs$object)))
  stop("07_coresh_search: coresh.query_signatures.object is required.")
if (is.null(qsigs$sets) || !is.list(qsigs$sets) || length(qsigs$sets) == 0L)
  stop("07_coresh_search: coresh.query_signatures.sets must be a non-empty list.")

sig_fp <- if (startsWith(qsigs$object, "/")) qsigs$object else file.path(PROJECT_ROOT, qsigs$object)
if (!file.exists(sig_fp))
  stop("07_coresh_search: query signature object not found: ", sig_fp)

signature_obj <- readRDS(sig_fp)
signature_sets <- signature_obj$sets %||% signature_obj

all_queries <- list()
for (i in seq_along(qsigs$sets)) {
  spec <- qsigs$sets[[i]]
  required_spec <- c("contrast", "direction", "gate")
  missing_spec <- required_spec[!nzchar(vapply(required_spec, function(k) spec[[k]] %||% "", character(1)))]
  if (length(missing_spec) > 0L)
    stop("07_coresh_search: coresh.query_signatures.sets[[", i,
         "]] missing required field(s): ", paste(missing_spec, collapse = ", "))

  contrast <- as.character(spec$contrast)
  direction <- as.character(spec$direction)
  gate <- as.character(spec$gate)
  slot_label <- paste0("sets[['", contrast, "']][['", direction, "']][['", gate, "']]")

  sig_symbols <- tryCatch(
    signature_sets[[contrast]][[direction]][[gate]],
    error = function(e) NULL
  )
  if (is.null(sig_symbols))
    stop("07_coresh_search: missing query signature slot ", slot_label, " in ", sig_fp)
  if (!is.character(sig_symbols) || length(sig_symbols) == 0L)
    stop("07_coresh_search: query signature slot ", slot_label,
         " must be a non-empty character vector; got ",
         paste(class(sig_symbols), collapse = "/"), " length ", length(sig_symbols))

  sig_symbols <- unique(sig_symbols[nzchar(sig_symbols)])
  if (length(sig_symbols) == 0L)
    stop("07_coresh_search: query signature slot ", slot_label,
         " contains no non-empty gene symbols.")

  query_name <- paste0("Q_sig_", contrast, "_", direction, "_", gate)
  if (query_name %in% names(all_queries))
    stop("07_coresh_search: duplicate query name from coresh.query_signatures: ", query_name)

  query_entrez <- sym2ent(sig_symbols, species = CORESH_SPECIES)
  all_queries[[query_name]] <- query_entrez
  message(sprintf("  %-40s %d symbols -> %d mouse Entrez",
                  query_name, length(sig_symbols), length(query_entrez)))
}

# ============================================================================
# 6. Combine, coerce to integer, drop sub-min queries, run the cached sweep
# ============================================================================
## coresh_batch() asserts is.integer(q) && length(q) >= 3 — coerce + drop singletons/pairs.

all_queries <- lapply(all_queries, as.integer)
too_small   <- vapply(all_queries, length, integer(1)) < MIN_Q
if (any(too_small)) {
  message("  dropping ", sum(too_small), " query/queries with k < ", MIN_Q, ": ",
          paste(names(all_queries)[too_small], collapse = ", "))
  all_queries <- all_queries[!too_small]
}
if (length(all_queries) == 0L)
  stop("07_coresh_search: no queries pass the size filter (k >= ", MIN_Q,
       "). Check coresh.query_signatures and mouse symbol-to-Entrez mapping.")
message(sprintf("[07_coresh_search] %d queries pass size filter: %s",
                length(all_queries), paste(names(all_queries), collapse = ", ")))

## Persist the exact queries searched with (provenance / reproducibility).
saveRDS(all_queries, file.path(DIR_OBJECTS, "coresh_query_entrez.rds"))

## ---- cached sweep over the mmu chunks (variance-only by default; rule 4 / >1 min) ----
## The path-keyed load_or_compute (from de_gsea_helpers.R) keys on filename only — it
## cannot see that the query SET changed. Force a recompute when the cache does not cover
## every current query.
sweep_fp <- file.path(DIR_OBJECTS, "coresh_ranked.rds")
force_sweep <- file.exists(sweep_fp) &&
  !all(names(all_queries) %in% unique(as.character(readRDS(sweep_fp)$query_name)))
if (force_sweep)
  message("  cached sweep missing queries (",
          paste(setdiff(names(all_queries), unique(as.character(readRDS(sweep_fp)$query_name))),
                collapse = ", "), ") — recomputing sweep + derived sets.")

results <- load_or_compute(sweep_fp, function() {
  coresh_batch(queries = all_queries, chunk_dir = chunk_dir,
               n_cores = N_CORES, pvalues = USE_PVALUES)
}, force = force_sweep)
stopifnot(all(c("query_name", "gse", "gpl", "pctVar", "pval", "size", "rank") %in% colnames(results)))
message(sprintf("[07_coresh_search] CoReSh sweep: %d rows across %d queries.",
                nrow(results), length(all_queries)))

res_df <- as.data.frame(results, stringsAsFactors = FALSE)

# ============================================================================
# 7. Derive CoReSh coregulation gene sets (CoReSh-to-GSEA bridge) -> .rds for stage 08
# ============================================================================
## Each top-ranking GSE implicitly defines a coregulation module (query + co-moving genes).
## build_coresh_gmt() projects every gene onto the query direction and keeps the top-|loading|
## genes, then converts Entrez->mouse symbols, size-filters [MIN_SET, MAX_SET], Jaccard-dedupes.
## Returns GMT lines ("name\t-\tgene1\tgene2..."). We parse them into the fgsea `pathways`
## named-list format and save as coresh_derived_sets.rds — the contract input for
## 08_coresh_derived_gsea.R (an .rds, matching the rest of this project's geneset_*.rds).

top_hits <- results[results$rank <= TOP_N, ]     # top-N GEO datasets per query (data.table)
message(sprintf("[07_coresh_search] extracting loadings from top-%d hits/query (%d GSEs).",
                TOP_N, nrow(top_hits)))

gmt_lines <- load_or_compute(file.path(DIR_OBJECTS, "coresh_gmt_lines.rds"), function() {
  build_coresh_gmt(
    top_hits          = top_hits,
    queries           = all_queries,        # named integer-Entrez list; names match top_hits$query_name
    chunk_dir         = chunk_dir,
    species           = CORESH_SPECIES,      # "mouse" -> ent2sym(..., "mouse") inside
    n_top             = N_DERIVE,
    min_size          = MIN_SET,
    max_size          = MAX_SET,
    jaccard_threshold = JACCARD_THRESH)
}, force = force_sweep)

## Parse GMT lines ("name\t-\tgenes...") into a named list of symbol vectors (fgsea pathways).
parse_gmt_lines <- function(lines) {
  if (length(lines) == 0L) return(list())
  out <- lapply(lines, function(ln) {
    toks <- strsplit(ln, "\t", fixed = TRUE)[[1]]
    if (length(toks) < 3L) return(NULL)
    genes <- toks[-(1:2)]                          # drop name + description ("-")
    unique(genes[nzchar(genes)])
  })
  names(out) <- vapply(lines, function(ln) strsplit(ln, "\t", fixed = TRUE)[[1]][1], character(1))
  Filter(function(g) !is.null(g) && length(g) > 0L, out)
}
derived_sets <- parse_gmt_lines(gmt_lines)

derived_fp <- file.path(DIR_OBJECTS, "coresh_derived_sets.rds")
saveRDS(derived_sets, derived_fp)                 # the contract object for 08_coresh_derived_gsea.R
if (length(derived_sets) == 0L)
  warning("[07_coresh_search] 0 CoReSh-derived sets passed size/Jaccard filters — ",
          "08_coresh_derived_gsea.R will have nothing to score. Check query coverage / compendium.")
message(sprintf("[07_coresh_search] derived sets: %d -> %s", length(derived_sets), derived_fp))

# ============================================================================
# 8. Ranked-compendium table -> 03_results/08_coresh/tables/coresh_ranked.csv
# ============================================================================
## The contract deliverable: the ranked public-mmu compendium for the study signature.
## Ordered by query then rank.

ranked_out <- res_df[order(res_df$query_name, res_df$rank), , drop = FALSE]
ranked_out <- round_numeric_cols(ranked_out)      # byte-stable doubles (de_gsea_helpers.R)
readr::write_csv(ranked_out, file.path(tbl_dir, "coresh_ranked.csv"))
message(sprintf("[07_coresh_search] ranked compendium -> %s (%d rows, %d datasets)",
                file.path(tbl_dir, "coresh_ranked.csv"),
                nrow(ranked_out), length(unique(ranked_out$gse))))

# ============================================================================
# 9. Provenance side-table (derived set -> source GSE / query / rank / pctVar)
# ============================================================================

set_names <- names(derived_sets)
prov <- if (length(set_names) == 0L) {
  tibble::tibble(set_name = character(), query_name = character(), gse = character(),
                 rank = integer(), pctVar = double())
} else dplyr::bind_rows(lapply(set_names, function(nm) {
  ## set names are "CORESH_<query_name>_<GSE>"; query_name itself may contain "_", so locate the GSE token
  toks    <- strsplit(sub("^CORESH_", "", nm), "_", fixed = TRUE)[[1]]
  gse_pos <- which(grepl("^GSE", toks))
  gse_idx <- if (length(gse_pos)) max(gse_pos) else NA_integer_
  gse     <- if (!is.na(gse_idx)) toks[gse_idx] else NA_character_
  qn      <- if (!is.na(gse_idx) && gse_idx > 1L) paste(toks[seq_len(gse_idx - 1L)], collapse = "_") else NA_character_
  hit     <- res_df[res_df$query_name == qn & res_df$gse == gse, , drop = FALSE]
  tibble::tibble(
    set_name   = nm,
    query_name = qn,
    gse        = gse,
    rank       = if (nrow(hit)) as.integer(hit$rank[1])   else NA_integer_,
    pctVar     = if (nrow(hit)) as.numeric(hit$pctVar[1]) else NA_real_)
}))
readr::write_csv(prov, file.path(tbl_dir, "coresh_provenance.csv"))
if (nrow(prov) > 0L && any(!prov$query_name %in% names(all_queries)))
  stop("07_coresh_search: derived-set provenance parser produced unknown query_name(s): ",
       paste(unique(prov$query_name[!prov$query_name %in% names(all_queries)]), collapse = ", "))
message(sprintf("[07_coresh_search] provenance -> %s (%d derived sets)",
                file.path(tbl_dir, "coresh_provenance.csv"), nrow(prov)))

# ============================================================================
# 10. Final structural asserts
# ============================================================================

stopifnot(
  file.exists(file.path(DIR_OBJECTS, "coresh_ranked.rds")),
  file.exists(file.path(DIR_OBJECTS, "coresh_derived_sets.rds")),
  file.exists(file.path(DIR_OBJECTS, "coresh_query_entrez.rds")),
  file.exists(file.path(tbl_dir, "coresh_ranked.csv")),
  file.exists(file.path(tbl_dir, "coresh_provenance.csv"))
)
message(sprintf(
  "[07_coresh_search] COMPLETE — %d queries, %d ranking rows, %d derived sets. (captions floor at L3-DE)",
  length(all_queries), nrow(results), length(derived_sets)))
