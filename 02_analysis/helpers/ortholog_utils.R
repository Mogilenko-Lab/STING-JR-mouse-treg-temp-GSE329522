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
##   stale query symbol     : normalised to the current MGI symbol via org.Mm.eg.db
##                            BEFORE babelgene sees it, then mapped back to the symbol
##                            the data carries (see normalise_mouse_query below).
## Signs/directions are preserved end to end: sets are split by DE direction upstream;
## the ranked metric (signed limma t) is carried through unchanged onto each human gene.
##
## WHY THE QUERY SIDE NEEDS NORMALISING AT ALL. babelgene 22.9 keys on CURRENT MGI
## symbols. This matrix was quantified against GRCm38 + GENCODE vM25, so 2,341 of its
## 19,679 symbols are no longer current. babelgene simply does not know those keys, they
## are absent from the returned edge table, and every caller detects them by set
## difference against the input — which is to say they were counted as
## `n_dropped_no_ortholog` and rendered in SIGNATURES.md as "no human ortholog: dropped".
## That label was wrong for 146 of the 6,510: they are one-to-one stale aliases whose
## current MGI symbol babelgene maps perfectly well, i.e. a VOCABULARY loss wearing an
## ORTHOLOGY label. Ddx58 is the one to remember — RIG-I, a cytosolic nucleic-acid sensor,
## significantly down at 39 °C in WT and at logFC -1.28 in KO_heat, and it was missing
## from the human projection background entirely because its current symbol is Rigi.
## Which stale symbols babelgene happens to know is idiosyncratic (Atp5a1, Atp5b, Gars,
## Kars, Nrd1 all map; Wars, Ddx58, Skp1a, Ero1l, Mpp5 do not), so it is measured here and
## never assumed.
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
    # The query-side vocabulary policy, DECLARED rather than implied: a matrix symbol
    # babelgene could not key AND that is no longer a current MGI symbol is re-asked under
    # its current symbol via org.Mm.eg.db, subject to the one-to-one + ownership guards in
    # normalise_mouse_query(). Set to "none" to reproduce the pre-fix behaviour exactly.
    stale_query_symbol      = "normalise_via_org_db",
    # Which mouse row supplies the ranked metric when a human symbol is claimed by both a
    # primary-query edge and a recovered one. The recovered rows are real mouse paralogs
    # (Gm9774 = Adrm1b, Gm4737 = Ahcyl, 2610002M06Rik = Chmp1b2, Lilr4b = Lilrb4b) and two
    # are pseudogene-named (Gm10282 = Hmgn2-ps, Gm14698 = Styx-ps), so babelgene maps them
    # onto a human symbol the canonical row already supplied. Left to max|t| alone the
    # recovery would REPLACE the canonical metric for 5 human symbols — MYEF2 -1.06 -> +3.82
    # and NSMCE3 +1.42 -> -1.70 are sign flips, sourced from a pseudogene row. That is a
    # change to the frozen contract with no coverage benefit, and the point of the recovery
    # is to restore genes that were lost, not to re-arbitrate genes that were not. So the
    # primary query wins, the recovered edge is still PUBLISHED in ortholog_map.tsv (the
    # many->one collision is real and was previously invisible), and the ranked lists are
    # strictly additive: 142 human symbols gained, 0 values changed.
    recovered_edge_precedence = "primary_query_wins",
    min_support             = 3L           # babelgene::orthologs(min_support=)
  )
}

## ---------------------------------------------------------------------------
## (0b) normalise_mouse_query(mouse_symbols, vocabulary, db, flagged_pairs) — the
##      query-side RECOVERY pass. Lifts a STALE matrix symbol UP to its current MGI symbol
##      so babelgene can key it, and returns the back-mapping so the published edge table
##      keeps the symbol the DATA carries.
##
##      This is the OPPOSITE traversal from helpers/symbol_alias.R::build_alias_map(),
##      which resolves a current reference symbol DOWN into this matrix's older
##      vocabulary. The one-to-one safety condition means something different in each
##      direction, so the two are deliberately separate functions. The collision policy
##      for orthology lives in this file, which is why the recovery pass does too.
##
##      RECOVERY, NOT REWRITE — and this is the thing that has to be got right. Rewriting
##      the query wholesale, i.e. handing babelgene the current symbol for every stale one,
##      is a REGRESSION: measured over this matrix's 19,679 symbols it recovers 144 genes
##      and LOSES 178, because babelgene 22.9's bundled orthology data predates
##      org.Mm.eg.db 3.22 and knows several legacy names it does not know under their
##      current ones (Atp5a1, Atp5b, Gars, Kars, Nrd1 all map as delivered). It also
##      silently CHANGED one already-mapped gene's human ortholog. So the primary query is
##      left exactly as it was, and only the symbols it could not key are re-asked. The
##      map may then only ever grow, which is asserted by the caller.
##
##      Guards, every one of which produces a counted row rather than a silent drop:
##        stale               a symbol that is already a current official MGI symbol needs
##                            no normalisation and is not a candidate. This is also the
##                            prefilter that makes the ownership question well-posed.
##        one alias -> one Entrez id
##                            a stale symbol resolving to >1 Entrez id names >1 gene here.
##                            rejected_query_symbol_ambiguous_in_org_db.
##        ownership           the current symbol must be the official symbol of THAT
##                            Entrez id and no other. rejected_symbol_belongs_to_another_gene.
##        no two queries -> one current symbol
##                            otherwise one babelgene edge would map back onto two matrix
##                            rows. rejected_multiple_query_symbols_for_current_symbol.
##        target not already a matrix row
##                            if the current symbol is ALSO its own row in this matrix, the
##                            matrix carries two features for one gene; normalising would
##                            silently merge them under the many-mouse->one-human rule and
##                            could change an existing gene's ranked metric.
##                            rejected_current_symbol_already_in_vocabulary.
##        flagged_for_review  pairs withheld by human decision, not by the guard.
##
##      @param mouse_symbols the CANDIDATES to normalise (the symbols the primary query
##                     failed to key), not the whole matrix.
##      @param vocabulary every symbol the matrix carries, for the occupancy guard. Defaults
##                     to `mouse_symbols`, which is only correct when they are the same set.
##      @return list(query = character  the symbols to re-ask babelgene (accepted pairs
##                     only, so an empty result means nothing to recover),
##                   ledger = data.frame(matrix_symbol, current_symbol, entrez_id,
##                     resolution)  one row per CANDIDATE, accepted or not)
## ---------------------------------------------------------------------------
normalise_mouse_query <- function(mouse_symbols, vocabulary = mouse_symbols, db = NULL,
                                  flagged_pairs = character()) {
  ms <- unique(as.character(mouse_symbols))
  ms <- ms[!is.na(ms) & ms != ""]
  vocab <- unique(as.character(vocabulary))
  vocab <- vocab[!is.na(vocab) & vocab != ""]
  empty <- data.frame(matrix_symbol = character(), current_symbol = character(),
                      entrez_id = character(), resolution = character(),
                      stringsAsFactors = FALSE)
  out <- function(ledger) list(query = character(0), ledger = ledger)

  if (is.null(db)) {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE)) {
      message("[ortholog_utils] org.Mm.eg.db unavailable — query symbols passed to ",
              "babelgene unnormalised; stale symbols will read as 'no ortholog'.")
      return(out(empty))
    }
    db <- org.Mm.eg.db::org.Mm.eg.db
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(out(empty))

  symbol_keys <- AnnotationDbi::keys(db, keytype = "SYMBOL")
  alias_keys  <- AnnotationDbi::keys(db, keytype = "ALIAS")
  # AnnotationDbi::select() ERRORS when none of its keys is valid for the keytype, so every
  # call is prefiltered against the live keyspace and returns early on zero length.
  sel <- function(keys, keytype, columns) {
    keys <- unique(keys[!is.na(keys) & nzchar(keys)])
    keys <- intersect(keys, if (identical(keytype, "SYMBOL")) symbol_keys else alias_keys)
    if (!length(keys)) return(NULL)
    suppressMessages(AnnotationDbi::select(db, keys = keys, keytype = keytype,
                                           columns = columns))
  }

  stale <- setdiff(ms, symbol_keys)
  if (!length(stale)) return(out(empty))
  sl <- sel(stale, "ALIAS", c("ENTREZID", "SYMBOL"))
  if (is.null(sl)) return(out(empty))
  sl <- sl[!is.na(sl$ENTREZID) & !is.na(sl$SYMBOL), , drop = FALSE]
  if (!nrow(sl)) return(out(empty))
  sl <- sl[!duplicated(sl[c("ALIAS", "ENTREZID", "SYMBOL")]), , drop = FALSE]

  cand <- data.frame(matrix_symbol = as.character(sl$ALIAS),
                     current_symbol = as.character(sl$SYMBOL),
                     entrez_id = as.character(sl$ENTREZID),
                     resolution = NA_character_, stringsAsFactors = FALSE)

  # (i) one alias -> one Entrez id
  n_eg <- tapply(cand$entrez_id, cand$matrix_symbol, function(x) length(unique(x)))
  cand$resolution[is.na(cand$resolution) & n_eg[cand$matrix_symbol] > 1L] <-
    "rejected_query_symbol_ambiguous_in_org_db"

  # (ii) ownership: the current symbol must belong to THAT Entrez id
  own <- sel(unique(cand$current_symbol), "SYMBOL", "ENTREZID")
  if (!is.null(own)) {
    own <- own[!is.na(own$ENTREZID), , drop = FALSE]
    owner <- tapply(own$ENTREZID, own$SYMBOL, function(x) paste(unique(x), collapse = "|"))
    bad <- !is.na(owner[cand$current_symbol]) &
      owner[cand$current_symbol] != cand$entrez_id
    bad[is.na(bad)] <- FALSE
    cand$resolution[is.na(cand$resolution) & bad] <-
      "rejected_symbol_belongs_to_another_gene"
  }

  # (iii) no two surviving queries may claim one current symbol
  live <- is.na(cand$resolution)
  n_q <- table(cand$current_symbol[live])
  cand$resolution[live & cand$current_symbol %in% names(n_q)[n_q > 1L]] <-
    "rejected_multiple_query_symbols_for_current_symbol"

  # (iv) the current symbol must not already be its own row in this matrix
  cand$resolution[is.na(cand$resolution) & cand$current_symbol %in% vocab] <-
    "rejected_current_symbol_already_in_vocabulary"

  # (v) the human-decision withholding, applied last so it is visible as a decision
  flagged <- paste0(cand$matrix_symbol, "->", cand$current_symbol) %in% flagged_pairs
  cand$resolution[is.na(cand$resolution) & flagged] <- "flagged_for_review"
  cand$resolution[is.na(cand$resolution)] <- "accepted"

  acc <- cand[cand$resolution == "accepted", , drop = FALSE]
  stopifnot("normalise_mouse_query: accepted pairs are not one-to-one" =
              !any(duplicated(acc$matrix_symbol)) && !any(duplicated(acc$current_symbol)))
  cand <- cand[order(match(cand$resolution, c("accepted", "flagged_for_review")),
                     cand$resolution, cand$matrix_symbol), , drop = FALSE]
  rownames(cand) <- NULL
  list(query = acc$current_symbol, ledger = cand)
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
##
##     QUERY NORMALISATION (policy stale_query_symbol = "normalise_via_org_db"). Two
##     passes, in this order, and the order is the safety property:
##       pass 1  every input symbol exactly as the data carries it — unchanged from the
##               pre-fix behaviour, so every edge it produced is still produced.
##       pass 2  only the symbols pass 1 could not key AND that are not current MGI
##               symbols are re-asked under their current symbol (normalise_mouse_query),
##               and the returned edges are mapped BACK onto the matrix symbol.
##     So `mouse_symbol` always stays the symbol the data actually carries, and the map can
##     only ever GROW — asserted below, because the naive one-pass rewrite loses 178 genes
##     (see normalise_mouse_query). Two columns record what happened per edge:
##       matrix_symbol_normalised_to  the current symbol pass 2 keyed on, or NA
##       normalisation_source         "org.Mm.eg.db" where recovered, "" otherwise
##     The per-candidate decision ledger is returned as the `query_normalisation`
##     attribute — 18_projection_export.R publishes it so a recovery is never inferred
##     from an absence. Pass normalise = FALSE to reproduce the pre-fix behaviour exactly.
## ---------------------------------------------------------------------------
build_ortholog_map <- function(mouse_symbols, min_support = 3L, normalise = TRUE,
                               flagged_pairs = character()) {
  if (!requireNamespace("babelgene", quietly = TRUE))
    stop("ortholog_utils: package 'babelgene' is required (same pinned/offline pkg as 00b).")
  ms <- unique(as.character(mouse_symbols))
  ms <- ms[!is.na(ms) & ms != ""]
  if (length(ms) == 0L) return(.empty_ortholog_map())

  ask <- function(genes) {
    if (!length(genes)) return(NULL)
    tryCatch(
      babelgene::orthologs(genes = genes, species = "mouse", human = FALSE,
                           min_support = min_support),
      error = function(e) {
        warning("ortholog_utils::build_ortholog_map: babelgene failed: ",
                conditionMessage(e), call. = FALSE); NULL
      })
  }
  as_edges <- function(raw, key_map = NULL) {
    if (is.null(raw) || nrow(raw) == 0L) return(NULL)
    e <- data.frame(queried_symbol = as.character(raw$symbol),
                    mouse_ensembl  = as.character(raw$ensembl),
                    human_symbol   = as.character(raw$human_symbol),
                    stringsAsFactors = FALSE)
    # Home again: an edge is re-keyed onto the matrix symbol. babelgene can return a symbol
    # we never asked for only if it case-normalises, so an unmatched key is dropped rather
    # than guessed at.
    e$mouse_symbol <- if (is.null(key_map)) e$queried_symbol
                      else unname(key_map[e$queried_symbol])
    e$matrix_symbol_normalised_to <- if (is.null(key_map)) NA_character_ else e$queried_symbol
    e <- e[!is.na(e$mouse_symbol) & e$mouse_symbol != "" &
             !is.na(e$human_symbol) & e$human_symbol != "", , drop = FALSE]
    if (!nrow(e)) return(NULL)
    e[, c("mouse_symbol", "mouse_ensembl", "human_symbol",
          "matrix_symbol_normalised_to"), drop = FALSE]
  }

  edges <- as_edges(ask(ms))

  # ---- pass 2: recover only what pass 1 could not key --------------------------------
  qledger <- .empty_query_ledger()
  if (isTRUE(normalise)) {
    unmapped <- setdiff(ms, if (is.null(edges)) character(0) else edges$mouse_symbol)
    norm <- normalise_mouse_query(unmapped, vocabulary = ms,
                                  flagged_pairs = flagged_pairs)
    qledger <- norm$ledger
    if (length(norm$query)) {
      acc <- qledger[qledger$resolution == "accepted", , drop = FALSE]
      back <- setNames(acc$matrix_symbol, acc$current_symbol)
      rec <- as_edges(ask(norm$query), key_map = back)
      if (!is.null(rec)) edges <- rbind(edges, rec)
    }
  }
  if (is.null(edges) || nrow(edges) == 0L)
    return(.empty_ortholog_map(query_normalisation = qledger))

  edges$normalisation_source <- ifelse(is.na(edges$matrix_symbol_normalised_to),
                                       "", "org.Mm.eg.db")
  # collapse identical (mouse,human) edges that differ only in support metadata. Pass-1
  # rows come first, so a duplicate can only ever be dropped from pass 2 — the recovery
  # never overwrites an edge the original query already produced.
  edges <- edges[!duplicated(edges[c("mouse_symbol", "human_symbol")]), , drop = FALSE]

  n_human <- table(edges$mouse_symbol)   # humans per mouse gene
  n_mouse <- table(edges$human_symbol)   # mouse genes per human symbol
  nh <- as.integer(n_human[edges$mouse_symbol])
  nm <- as.integer(n_mouse[edges$human_symbol])
  edges$mapping_type <- ifelse(nh > 1L, "one2many",
                        ifelse(nm > 1L, "many2one", "one2one"))
  edges$babelgene_version <- as.character(utils::packageVersion("babelgene"))
  edges <- edges[order(edges$mouse_symbol, edges$human_symbol), , drop = FALSE]
  edges <- edges[, c("mouse_symbol", "mouse_ensembl", "human_symbol", "mapping_type",
                     "matrix_symbol_normalised_to", "normalisation_source",
                     "babelgene_version"), drop = FALSE]
  rownames(edges) <- NULL
  structure(edges, query_normalisation = qledger)
}

.empty_query_ledger <- function() {
  data.frame(matrix_symbol = character(), current_symbol = character(),
             entrez_id = character(), resolution = character(),
             stringsAsFactors = FALSE)
}

.empty_ortholog_map <- function(query_normalisation = NULL) {
  structure(
    data.frame(mouse_symbol = character(), mouse_ensembl = character(),
               human_symbol = character(), mapping_type = character(),
               matrix_symbol_normalised_to = character(),
               normalisation_source = character(),
               babelgene_version = character(), stringsAsFactors = FALSE),
    query_normalisation = query_normalisation %||% .empty_query_ledger())
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
##       many mouse → one human : keep the row with MAX |metric| per human symbol —
##                                EXCEPT that a primary-query edge outranks a recovered
##                                one regardless of |metric|
##                                (policy recovered_edge_precedence = "primary_query_wins";
##                                the reasoning is in default_ortholog_policy). A human
##                                symbol the map already carried therefore keeps the
##                                metric it already had, and the recovery can only add.
##       no human ortholog      : dropped (inner join).
##     Returns a 2-col data.frame (human_symbol, t) sorted by signed metric DESCENDING —
##     the fgsea/decoupleR .rnk shape. Column is named `t` (the project rank metric).
## ---------------------------------------------------------------------------
map_ranked_list <- function(ranked_df, omap, mouse_col = "gene_symbol", stat_col = "t",
                            primary_query_wins = TRUE) {
  if (is.null(ranked_df) || nrow(ranked_df) == 0L)
    return(data.frame(human_symbol = character(), t = numeric(), stringsAsFactors = FALSE))
  m <- data.frame(
    mouse_symbol = as.character(ranked_df[[mouse_col]]),
    stat         = suppressWarnings(as.numeric(ranked_df[[stat_col]])),
    stringsAsFactors = FALSE)
  m <- m[!is.na(m$mouse_symbol) & m$mouse_symbol != "" & !is.na(m$stat), , drop = FALSE]
  edges <- omap[c("mouse_symbol", "human_symbol")]
  edges$recovered <- if (isTRUE(primary_query_wins) &&
                         "normalisation_source" %in% names(omap))
    as.integer(omap$normalisation_source == "org.Mm.eg.db") else 0L
  # edge-expand: one mouse -> many human duplicates the metric onto each human gene.
  j <- merge(m, edges, by = "mouse_symbol")
  j <- j[!is.na(j$human_symbol) & j$human_symbol != "", , drop = FALSE]
  if (nrow(j) == 0L)
    return(data.frame(human_symbol = character(), t = numeric(), stringsAsFactors = FALSE))
  # many mouse -> one human: primary edges first, then the most extreme |metric|.
  j <- j[order(j$recovered, -abs(j$stat)), , drop = FALSE]
  j <- j[!duplicated(j$human_symbol), , drop = FALSE]
  out <- data.frame(human_symbol = j$human_symbol, t = j$stat, stringsAsFactors = FALSE)
  out[order(-out$t), , drop = FALSE]
}

## ---------------------------------------------------------------------------
## End of ortholog_utils.R.
## Provides: default_ortholog_policy, normalise_mouse_query, babelgene_provenance,
##   build_ortholog_map, ortholog_coverage, map_binary_set, map_ranked_list
##   (+ internal .empty_ortholog_map, .empty_query_ledger).
## ---------------------------------------------------------------------------
message("[ortholog_utils] loaded (mouse->human babelgene mapping + org.Mm.eg.db query ",
        "normalisation; lazy deps: babelgene, org.Mm.eg.db).")
