## ortholog_utils.R — mouse→human ortholog mapping helper (babelgene, pinned/offline)
## ===========================================================================
## PROJECT HELPER, NOT A SKILL. Ortholog mapping has no skill in the active
## base + pathway-signature catalog, so the mapping policy lives here as functions
## (per the phase-0 scaffold rule: encode the collision policy in a helper, never
## inline it in the numbered scripts). The REVERSE of 00b_curate_lombardi_hif.R,
## which maps human→mouse (babelgene::orthologs(..., human = TRUE)); this file maps
## mouse→human (human = FALSE) with the SAME pinned, offline babelgene package — no
## new dependency, no network. Downstream: 18_projection_export.R freezes the human
## contract with these functions; 17_signature_derive.R uses ortholog_coverage() for
## a display-only dry-run preview.
##
## COMPUTE ONLY — no plotting, no file writes (the babelgene call is its only effect).
## Sourced AFTER 02_analysis/config/config.R (needs `%||%`; defensively self-defines).
##
## babelgene::orthologs(genes, species = "mouse", human = FALSE) returns per-EDGE rows
##   with columns: symbol (mouse), ensembl (mouse), human_symbol, ... (see the raw
##   frame). One mouse gene may yield several rows (one→many human); several mouse
##   genes may share one human_symbol (many→one human). `min_support` is babelgene's
##   orthology-source vote floor (default 3); recorded in provenance.
##
## COLLISION POLICY (defaults; overridable via analysis_config.yaml::decisions.projection):
##   one mouse → many human : binary up/down sets take the UNION of human orthologs;
##                            the ranked list assigns the mouse metric to EACH human
##                            ortholog (map_ranked_list expands the edge).
##   many mouse → one human : the ranked list keeps the entry with MAX |metric| per
##                            human symbol (the more extreme rank, never averaged);
##                            binary sets dedupe by union.
##   no human ortholog      : dropped and logged (every dropped gene auditable, the
##                            00b audit-log precedent).
## Signs/directions are preserved end to end: sets are split by DE direction upstream;
## the ranked metric (signed limma t) is carried through unchanged onto each human gene.
## ===========================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
}

## ---------------------------------------------------------------------------
## (0) default_ortholog_policy() — the collision policy as a plain list, so the
##     numbered scripts and the decisions block share ONE source of truth. Any
##     value here may be overridden by analysis_config.yaml::decisions.projection.
## ---------------------------------------------------------------------------
default_ortholog_policy <- function() {
  list(
    one_mouse_to_many_human = "union",     # binary sets: union; ranked: metric to each human
    many_mouse_to_one_human = "max_abs_t", # ranked: keep max|metric| per human; binary: union
    no_human_ortholog       = "drop",      # dropped + logged (auditable)
    min_support             = 3L           # babelgene::orthologs(min_support=)
  )
}

## ---------------------------------------------------------------------------
## (1) babelgene_provenance() — version + bundled orthology-data date, for the
##     provenance stamp the frozen contract (SIGNATURES.md) must record. `min_support`
##     echoes the value actually used so the map is reproducible.
## ---------------------------------------------------------------------------
babelgene_provenance <- function(min_support = 3L) {
  if (!requireNamespace("babelgene", quietly = TRUE))
    stop("ortholog_utils: package 'babelgene' is required (same pinned/offline pkg as 00b).")
  data_date <- tryCatch(as.character(utils::packageDescription("babelgene")$Date),
                        error = function(e) NA_character_)
  list(
    package        = "babelgene",
    version        = as.character(utils::packageVersion("babelgene")),
    data_date      = data_date,
    mapping_dir    = "mouse->human (species='mouse', human=FALSE)",
    min_support    = as.integer(min_support)
  )
}

## ---------------------------------------------------------------------------
## (2) build_ortholog_map(mouse_symbols, min_support = 3L) — the applied edge table.
##     Returns a tidy data.frame, one row per (mouse_symbol, human_symbol) edge, with a
##     per-edge `mapping_type` and the babelgene version. Mutually-exclusive types:
##       one2many  : this mouse gene maps to >1 human ortholog       (n_human > 1)
##       many2one  : exactly 1 human, but >1 mouse map to it         (n_human == 1 & n_mouse > 1)
##       one2one   : exactly 1 human and no other mouse maps to it    (n_human == 1 & n_mouse == 1)
##     Unmapped mouse genes (no ortholog at this min_support) are simply absent from
##     the table — callers detect them by set difference against the input.
## ---------------------------------------------------------------------------
build_ortholog_map <- function(mouse_symbols, min_support = 3L) {
  if (!requireNamespace("babelgene", quietly = TRUE))
    stop("ortholog_utils: package 'babelgene' is required (same pinned/offline pkg as 00b).")
  ms <- unique(as.character(mouse_symbols))
  ms <- ms[!is.na(ms) & ms != ""]
  if (length(ms) == 0L) return(.empty_ortholog_map())

  raw <- tryCatch(
    babelgene::orthologs(genes = ms, species = "mouse", human = FALSE,
                         min_support = min_support),
    error = function(e) {
      warning("ortholog_utils::build_ortholog_map: babelgene failed: ",
              conditionMessage(e), call. = FALSE); NULL
    })
  if (is.null(raw) || nrow(raw) == 0L) return(.empty_ortholog_map())

  edges <- data.frame(
    mouse_symbol  = as.character(raw$symbol),
    mouse_ensembl = as.character(raw$ensembl),
    human_symbol  = as.character(raw$human_symbol),
    stringsAsFactors = FALSE)
  edges <- edges[!is.na(edges$human_symbol) & edges$human_symbol != "" &
                 !is.na(edges$mouse_symbol) & edges$mouse_symbol != "", , drop = FALSE]
  # collapse identical (mouse,human) edges that differ only in support metadata
  edges <- edges[!duplicated(edges[c("mouse_symbol", "human_symbol")]), , drop = FALSE]
  if (nrow(edges) == 0L) return(.empty_ortholog_map())

  n_human <- table(edges$mouse_symbol)   # humans per mouse gene
  n_mouse <- table(edges$human_symbol)   # mouse genes per human symbol
  nh <- as.integer(n_human[edges$mouse_symbol])
  nm <- as.integer(n_mouse[edges$human_symbol])
  edges$mapping_type <- ifelse(nh > 1L, "one2many",
                        ifelse(nm > 1L, "many2one", "one2one"))
  edges$babelgene_version <- as.character(utils::packageVersion("babelgene"))
  edges <- edges[order(edges$mouse_symbol, edges$human_symbol), , drop = FALSE]
  rownames(edges) <- NULL
  edges
}

.empty_ortholog_map <- function() {
  data.frame(mouse_symbol = character(), mouse_ensembl = character(),
             human_symbol = character(), mapping_type = character(),
             babelgene_version = character(), stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## (3) ortholog_coverage(mouse_symbols, omap = NULL, min_support = 3L) — the dry-run
##     mapping-loss preview, PER INPUT GENE (mouse side). Classes each input gene as:
##       mapped_1to1 : maps to exactly one human ortholog (the input gene's own view;
##                     a human-side many→one collision still counts here — the loss it
##                     causes is on the human side, surfaced by build_ortholog_map()).
##       one_to_many : maps to >1 human ortholog.
##       unmapped    : no human ortholog at this min_support (dropped downstream).
##     Pass a prebuilt `omap` (e.g. one map over the whole DE universe) to avoid a
##     second babelgene call. Returns a one-row data.frame of counts.
## ---------------------------------------------------------------------------
ortholog_coverage <- function(mouse_symbols, omap = NULL, min_support = 3L) {
  ms <- unique(as.character(mouse_symbols))
  ms <- ms[!is.na(ms) & ms != ""]
  if (is.null(omap)) omap <- build_ortholog_map(ms, min_support = min_support)
  n_human_per_mouse <- table(omap$mouse_symbol[omap$mouse_symbol %in% ms])
  mapped   <- ms %in% names(n_human_per_mouse)
  many_nm  <- names(n_human_per_mouse)[n_human_per_mouse > 1L]
  one2many <- ms %in% many_nm
  data.frame(
    n_input     = length(ms),
    mapped_1to1 = sum(mapped & !one2many),
    one_to_many = sum(one2many),
    unmapped    = sum(!mapped),
    stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## (4) map_binary_set(mouse_symbols, omap) — apply the map to a binary up/down set.
##     Collision policy: UNION of human orthologs (one→many contributes all its human
##     genes; many→one collapses naturally via unique()); unmapped genes dropped.
##     Returns a sorted, de-duplicated character vector of human symbols.
## ---------------------------------------------------------------------------
map_binary_set <- function(mouse_symbols, omap) {
  ms <- unique(as.character(mouse_symbols))
  hs <- omap$human_symbol[omap$mouse_symbol %in% ms]
  hs <- hs[!is.na(hs) & hs != ""]
  sort(unique(hs))
}

## ---------------------------------------------------------------------------
## (5) map_ranked_list(ranked_df, omap, mouse_col, stat_col) — apply the map to a
##     ranked list. Collision policy:
##       one mouse → many human : the mouse metric is assigned to EACH human ortholog
##                                (edge expansion via the join).
##       many mouse → one human : keep the row with MAX |metric| per human symbol.
##       no human ortholog      : dropped (inner join).
##     Returns a 2-col data.frame (human_symbol, t) sorted by signed metric DESCENDING —
##     the fgsea/decoupleR .rnk shape. Column is named `t` (the project rank metric).
## ---------------------------------------------------------------------------
map_ranked_list <- function(ranked_df, omap, mouse_col = "gene_symbol", stat_col = "t") {
  if (is.null(ranked_df) || nrow(ranked_df) == 0L)
    return(data.frame(human_symbol = character(), t = numeric(), stringsAsFactors = FALSE))
  m <- data.frame(
    mouse_symbol = as.character(ranked_df[[mouse_col]]),
    stat         = suppressWarnings(as.numeric(ranked_df[[stat_col]])),
    stringsAsFactors = FALSE)
  m <- m[!is.na(m$mouse_symbol) & m$mouse_symbol != "" & !is.na(m$stat), , drop = FALSE]
  # edge-expand: one mouse -> many human duplicates the metric onto each human gene.
  j <- merge(m, omap[c("mouse_symbol", "human_symbol")], by = "mouse_symbol")
  j <- j[!is.na(j$human_symbol) & j$human_symbol != "", , drop = FALSE]
  if (nrow(j) == 0L)
    return(data.frame(human_symbol = character(), t = numeric(), stringsAsFactors = FALSE))
  # many mouse -> one human: keep the most extreme |metric| per human symbol.
  j <- j[order(-abs(j$stat)), , drop = FALSE]
  j <- j[!duplicated(j$human_symbol), , drop = FALSE]
  out <- data.frame(human_symbol = j$human_symbol, t = j$stat, stringsAsFactors = FALSE)
  out[order(-out$t), , drop = FALSE]
}

## ---------------------------------------------------------------------------
## End of ortholog_utils.R.
## Provides: default_ortholog_policy, babelgene_provenance, build_ortholog_map,
##   ortholog_coverage, map_binary_set, map_ranked_list  (+ internal .empty_ortholog_map).
## ---------------------------------------------------------------------------
message("[ortholog_utils] loaded (mouse->human babelgene mapping; lazy dep: babelgene).")
