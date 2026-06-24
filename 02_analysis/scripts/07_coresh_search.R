# 07_coresh_search.R — COMPUTE
## CoReSh (COregulation REgulation SearcH) — rank the public MOUSE (mmu) GEO compendium
## by how strongly the study's heat / IFN-ISG / HIF signature is CO-REGULATED in each
## public dataset (PCA-inspired pctVar score), then derive coregulation gene sets from the
## top hits for downstream GSEA (08_coresh_derived_gsea.R).
##
## Per the `coresh-signature-search` skill
## (01_modules/SciAgent-toolkit/skills/coresh-signature-search/SKILL.md):
##   * search engine  : coresh_batch()        (thin wrapper around fgsea::geseca; coreshMatch())
##   * derived sets   : build_coresh_gmt()     (CoReSh-to-GSEA bridge; gene loadings on top hits)
##   * ID conversion  : sym2ent()/ent2sym(species = "mouse")   (chunk rownames are integer mouse Entrez)
##   * data dependency: preprocessed_chunks/mmu/*_full_objects.qs2   (~20 GB Synapse compendium, syn66227307)
## Skill scripts sourced from the canonical skill dir (NOT 01_modules/.ref — that is the 14839 layout).
##
## Run from project root:
##   Rscript 02_analysis/scripts/07_coresh_search.R
##
## ============================================================================
## !!! LOUD ASSUMPTIONS / PROVISIONING — READ BEFORE RUNNING (this is the most
##     uncertain sweep arm; the recon flagged the CoReSh compendium may NOT yet be
##     provisioned in STING) !!!
## ----------------------------------------------------------------------------
## (A) DATA DEPENDENCY NOT YET PROVISIONED IN STING.
##     There is NO coresh chunk dir under 00_data/references/ at write time
##     (only gatom/, gene_sets/, networks/ exist). The ~20 GB mmu compendium must be
##     downloaded from Synapse (syn66227307) by the OWNER in-container — Claude/this
##     environment CANNOT do the Synapse download (needs a personal access token).
##     The script GUARDS this with an explicit stop() and points at the provisioning note:
##         docs/_internal/reasoning/2026-06-24_02_coresh-provisioning.md
##     NO compendium is fabricated. If you see the stop(), provision the data first.
##
## (B) CONFIG KEYS NOT YET DECLARED.
##     analysis_config.yaml has NO `coresh:`, NO `key_genes:`, and `paths:` has NO
##     `coresh_chunks` at write time (14839 has all three). This script therefore reads
##     every such key DEFENSIVELY via `%||%` with documented in-script defaults, and
##     falls back to the config.R gene-vector constants (ISG_MARKERS / HIF_GLYCO_MARKERS /
##     THERMO_MARKERS / KEY_TFS) for the curated queries. To make this a first-class arm,
##     the owner should add a `coresh:` block + `paths.coresh_chunks` (see the note). Until
##     then the in-script defaults reproduce the 14839 behaviour (species=mouse, top_n=5,
##     min_query_size=3, n_cores=4).
##
## (C) CHUNK PATH ASSUMED.
##     Default chunk dir = paths.coresh_chunks (if set) else the 14839 convention
##     00_data/references/coresh/current/preprocessed_chunks/mmu (overridable via the
##     CORESH_CHUNKS env var, per the skill's path convention). Owner: confirm in-container.
##
## (D) DE HUB CONTRACT.
##     Consumes 03_results/objects/02_de_results.rds (named list of 7 topTables keyed by the
##     YAML contrast names, gene-SYMBOL rownames + `t`). We NEVER re-fit DE here.
## ============================================================================
##
## Outputs (compute-only; lazy heavy deps; idempotent via load_or_compute)
## ----------------------------------------------------------------------------
##   03_results/objects/coresh_ranked.rds              — full per-query ranking (data.table; cache)
##   03_results/objects/coresh_derived_sets.rds        — CoReSh-derived gene sets (named list, fgsea
##                                                        pathways format) for 08_coresh_derived_gsea.R
##   03_results/objects/coresh_query_entrez.rds        — the exact queries searched with (provenance)
##   03_results/08_coresh/tables/coresh_ranked.csv     — the ranked-compendium table (curated + DE-derived)
##   03_results/08_coresh/tables/coresh_provenance.csv — derived-set -> source GSE / pctVar trace
##   03_results/08_coresh/tables/by_contrast/<c>/coresh_ranked.csv — per-contrast DE-derived ranking

# ============================================================================
# 0. Environment setup  (config.R FIRST, then de_gsea_helpers.R — per contract)
# ============================================================================

source("02_analysis/config/config.R")            # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS, SPECIES,
                                                  # stage_dir(), load_or_compute (config.R variant), %||%,
                                                  # ISG_MARKERS, HIF_GLYCO_MARKERS, THERMO_MARKERS, KEY_TFS
source("02_analysis/helpers/de_gsea_helpers.R")  # path-keyed load_or_compute, load_de_results,
                                                  # load_custom_geneset, round_numeric_cols
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
by_contrast_dir_name <- YAML_CONFIG$figures$by_contrast_dir %||% "by_contrast"

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
DE_TOP_N       <- as.integer(coresh_cfg$de_top_n        %||% 20L)  # top up-genes/contrast for DE-derived queries

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
# 3. Load the DE hub (the GSEA/CoReSh ranking source — NEVER re-fit DE here)
# ============================================================================

de_results <- load_de_results()                  # asserts Symbol rownames + `t` per contrast
stopifnot(length(de_results) >= 1L)
message(sprintf("[07_coresh_search] DE hub: %d contrasts (%s)",
                length(de_results), paste(names(de_results), collapse = ", ")))

# ============================================================================
# 4. CURATED queries — the heat / IFN-ISG / HIF signatures of interest
# ============================================================================
## Mirrors 14839's curated CoReSh queries (Q_curated_<group>) but built from STING's
## biology (the science-design brief §5): the IFN/ISG arm (cGAS-dependent), the HIF /
## glycolysis arm (no detectable cGAS-dependence at n=5), the heat-shock thermometer,
## the cGAS-STING sensor TFs, and — the FLAGSHIP query — the PUBLISHED Lombardi 2022
## 48-gene consensus HIF set. The recon's headline rationale for running CoReSh at all:
## it EXPOSES Biomni's hand-made 16-gene HIF list by quantifying whether a HIF signature
## is a real co-regulated module (positive control = a HIF query ranks hypoxia/VHL/
## tumor-hypoxia GEO studies in the top 20). The Lombardi-48 set is the curation-
## independent benchmark for exactly that test.
##
## Source of the curated symbol vectors (config, not hardcoded):
##   * If a `key_genes:` block exists in config, use it (14839 idiom).
##   * Else fall back to the config.R analysis constants (the single source of truth here).

key_genes <- YAML_CONFIG$key_genes
if (is.null(key_genes)) {
  key_genes <- list(
    isg_ifn      = ISG_MARKERS,        # Ifit1, Isg15, Irf7, Oasl2, Mx1, Stat1, Cxcl10  (cGAS-dependent arm)
    hif_glyco    = HIF_GLYCO_MARKERS,  # Slc2a1, Vegfa, Egln3, Bnip3, Pgk1, Ldha, Aldoa, Hk2
    heat_shock   = THERMO_MARKERS,     # Hspa1b, Hsph1, Hspa1a, Dnajb1  (the true thermal program)
    sensor_tfs   = KEY_TFS             # Irf3/7/1, Stat1/2, Hif1a, Epas1, Nfkb1, Rela
  )
  message("[07_coresh_search] coresh.key_genes absent — using config.R constants ",
          "(ISG_MARKERS / HIF_GLYCO_MARKERS / THERMO_MARKERS / KEY_TFS).")
} else {
  message("[07_coresh_search] using config key_genes block (", length(key_genes), " groups).")
}
key_genes <- lapply(key_genes, function(v) unique(as.character(unlist(v))))

## FLAGSHIP curated query: the Lombardi-2022 48-gene consensus HIF set (the published,
## curation-independent HIF benchmark — built by 00b_curate_lombardi_hif.R). Loaded from
## databases.custom if declared; warned-and-skipped if absent (do not error the whole arm).
lombardi_genes <- tryCatch({
  custom <- YAML_CONFIG$databases$custom %||% list()
  hit <- Filter(function(d) identical(d$name, "Lombardi2022_HIF"), custom)
  if (length(hit) == 0L) {
    message("[07_coresh_search] Lombardi2022_HIF not in databases.custom — skipping the flagship HIF query.")
    character(0)
  } else {
    p <- file.path(PROJECT_ROOT, hit[[1]]$path)
    if (!file.exists(p)) {
      warning("[07_coresh_search] Lombardi2022_HIF rds missing (", p,
              ") — run 00b_curate_lombardi_hif.R first. Skipping the flagship HIF query.")
      character(0)
    } else {
      sets <- load_custom_geneset(p)            # list(Lombardi2022_HIF = <mouse symbols>)
      unique(as.character(unlist(sets, use.names = FALSE)))
    }
  }
}, error = function(e) {
  warning("[07_coresh_search] could not load Lombardi2022_HIF: ", conditionMessage(e)); character(0)
})
if (length(lombardi_genes) > 0L) key_genes$lombardi_hif <- lombardi_genes

## symbol -> integer mouse Entrez (sym2ent warns on unmapped; does not error)
queries_curated <- lapply(key_genes, sym2ent, species = CORESH_SPECIES)
names(queries_curated) <- paste0("Q_curated_", names(key_genes))
for (nm in names(queries_curated))
  message(sprintf("  %-28s %d symbols -> %d mouse Entrez", nm,
                  length(key_genes[[sub("^Q_curated_", "", nm)]]), length(queries_curated[[nm]])))

# ============================================================================
# 5. DE-DERIVED queries — top-N up genes per contrast (covariation seed)
# ============================================================================
## Mirror 14839 §5.5: seed the sweep from EVERY contrast (covariation is meaningful for
## all directions). Headline layout is WT_heat | KO_heat | Interaction (+ Temp_main), but
## we run all 7 so the per-contrast tables are complete. Filter MOUSE housekeeping
## (ribosomal Rps/Rpl, mito-encoded mt-, haemoglobin Hb*) before conversion — these
## co-regulate everywhere and would confound pctVar (skill Common-Pitfalls).

top_n_entrez <- function(de_df, label, n = DE_TOP_N, species = CORESH_SPECIES) {
  if (is.null(de_df) || nrow(de_df) == 0L) { message("  skip ", label, ": empty DE"); return(NULL) }
  stat_col <- intersect(c("t", "stat", "statistic"), colnames(de_df))[1]
  if (is.na(stat_col)) { message("  skip ", label, ": no t/stat column"); return(NULL) }
  ord  <- order(de_df[[stat_col]], decreasing = TRUE)             # top up-regulated (numerator-high)
  syms <- rownames(de_df)
  if (is.null(syms) && "gene_symbol" %in% colnames(de_df)) syms <- as.character(de_df$gene_symbol)
  syms <- syms[ord][seq_len(min(n, nrow(de_df)))]
  syms <- syms[!grepl("^(Rp[sl]\\d+|mt-|Hb[abdeg]\\d*)", syms)]   # MOUSE Title-case housekeeping
  if (length(syms) == 0L) { message("  skip ", label, ": all housekeeping"); return(NULL) }
  e <- sym2ent(syms, species = species)
  message(sprintf("  %-28s %d symbols -> %d mouse Entrez", label, length(syms), length(e)))
  e
}

queries_derived <- list()
for (co in names(de_results)) {
  q <- top_n_entrez(de_results[[co]], paste0("Q_de_", co))
  if (!is.null(q)) queries_derived[[paste0("Q_de_", co)]] <- q
}

# ============================================================================
# 6. Combine, coerce to integer, drop sub-min queries, run the cached sweep
# ============================================================================
## coresh_batch() asserts is.integer(q) && length(q) >= 3 — coerce + drop singletons/pairs.

all_queries <- c(queries_curated, queries_derived)
all_queries <- lapply(all_queries, as.integer)
too_small   <- vapply(all_queries, length, integer(1)) < MIN_Q
if (any(too_small)) {
  message("  dropping ", sum(too_small), " query/queries with k < ", MIN_Q, ": ",
          paste(names(all_queries)[too_small], collapse = ", "))
  all_queries <- all_queries[!too_small]
}
if (length(all_queries) == 0L)
  stop("07_coresh_search: no queries pass the size filter (k >= ", MIN_Q,
       "). Check key_genes / DE hub.")
message(sprintf("[07_coresh_search] %d queries pass size filter: %s",
                length(all_queries), paste(names(all_queries), collapse = ", ")))

## Persist the exact queries searched with (provenance / reproducibility).
saveRDS(all_queries, file.path(DIR_OBJECTS, "coresh_query_entrez.rds"))

## ---- cached sweep over the mmu chunks (variance-only by default; rule 4 / >1 min) ----
## The path-keyed load_or_compute (from de_gsea_helpers.R) keys on filename only — it
## cannot see that the query SET grew. Force a recompute when the cache does not cover
## every current query (e.g. newly added contrasts / a newly provisioned Lombardi set).
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
## 08_coresh_derived_gsea.R (NOT a .gmt, matching the rest of this project's geneset_*.rds).

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
## Ordered by query then rank so curated (Q_curated_*) and DE-derived (Q_de_*) blocks are
## contiguous and human-readable.

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
message(sprintf("[07_coresh_search] provenance -> %s (%d derived sets)",
                file.path(tbl_dir, "coresh_provenance.csv"), nrow(prov)))

# ============================================================================
# 10. Per-contrast DE-derived rankings -> tables/by_contrast/<contrast>/coresh_ranked.csv
# ============================================================================
## One ranking per seeded contrast (the headline three-column read WT_heat | KO_heat |
## Interaction plus Temp_main are all here). Curated cross-contrast rankings stay at the
## tables/ root (coresh_ranked.csv above). Mirrors 14839 §5.8b + the house Results convention.

for (co in names(de_results)) {
  qn <- paste0("Q_de_", co)
  if (!qn %in% res_df$query_name) next            # query dropped (sub-min / housekeeping) — no empty dir
  rk  <- res_df[res_df$query_name == qn, , drop = FALSE]
  rk  <- round_numeric_cols(rk[order(rk$rank), , drop = FALSE])
  cdir <- file.path(tbl_dir, by_contrast_dir_name, co)
  dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(rk, file.path(cdir, "coresh_ranked.csv"))
}
message("[07_coresh_search] per-contrast DE-derived rankings written under tables/",
        by_contrast_dir_name, "/")

# ============================================================================
# 11. Final structural asserts
# ============================================================================

stopifnot(
  file.exists(file.path(DIR_OBJECTS, "coresh_ranked.rds")),
  file.exists(file.path(DIR_OBJECTS, "coresh_derived_sets.rds")),
  file.exists(file.path(DIR_OBJECTS, "coresh_query_entrez.rds")),
  file.exists(file.path(tbl_dir, "coresh_ranked.csv")),
  file.exists(file.path(tbl_dir, "coresh_provenance.csv"))
)
message(sprintf(
  "[07_coresh_search] COMPLETE — %d queries, %d ranking rows, %d derived sets. (provisional sample labels; captions floor at L3-DE)",
  length(all_queries), nrow(results), length(derived_sets)))
