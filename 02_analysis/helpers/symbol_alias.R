#!/usr/bin/env Rscript
# symbol_alias.R — resolve a reference set's MGI symbols into this matrix's vocabulary.
# =============================================================================
# PROJECT HELPER, NOT A SKILL. The compartments of this super-repo are separate git
# submodules with separate remotes and neither may import from the other, so this file is
# a deliberate sibling of human_treg_arthritis/02_analysis/helpers/symbol_alias.R rather
# than a shared module. Same interface, same guard, different annotation database.
#
# The defect this exists for. The collaborators' delivered CPM table was quantified
# against GRCm38 + GENCODE vM25, so this compartment's count matrix is frozen to that
# build's MGI vintage. 2,341 of the 19,679 modelled symbols are no longer current
# official MGI symbols, and they cluster: the whole ATP-synthase block arrives as
# Atp5a1 Atp5b Atp5c1 Atp5d Atp5e Atp5g1..g3 Atp5h Atp5j Atp5j2 Atp5k Atp5l Atp5o Atpif1
# while MitoCarta 3.0 ships Atp5f1a..Atp5po. Every match in this compartment is an exact
# string match, so those 15 rows leave MITOPATHWAYS_OXPHOS.Complex_V silently: the set
# matched 7 of its 22 genes, fell under gsea_min_size = 15, and was therefore absent from
# master_gsea_table.csv entirely — no contrast, no row, no NES. For a compartment whose
# stated mechanism is mitochondrial/metabolic at 39 °C, that is not a rounding error.
#
# The direction of resolution. The reference symbol is NEWER and the matrix symbol is
# OLDER, so a reference symbol is resolved DOWN into the vocabulary the data carries.
# The OPPOSITE traversal — a stale matrix symbol lifted UP to current, which is what
# babelgene needs on its query side — means something different under the one-to-one
# safety condition and lives in ortholog_utils.R::normalise_mouse_query(). Do not try to
# make one function do both.
#
# The hazard that shapes the guard. Retired symbols were often reassigned as the official
# symbol of a DIFFERENT gene, so accepting a pair blindly attaches one gene's expression
# to another gene's set membership. A candidate that is the official symbol of any other
# Entrez id is therefore rejected, counted, and published beside the accepted pairs.
#
# Alias resolution is a correctness fix and never a way to grow a set. Nothing is accepted
# that does not survive the ownership guard, and the reporting contract is a ledger with a
# bucket per cause — matched, matched-via-alias, expression-filtered, below-detection —
# never a pass/fail recovery floor.
#
# Provides:
#   build_alias_map(reference_symbols, matrix_vocabulary, db, flagged_pairs)
#   accepted_pairs(alias_map)   the reference_symbol -> matrix_symbol lookup, flagged excluded
#   resolve_sets(sets, matrix_vocabulary, alias_map, replace)
#   symbol_ledger(sets, alias_map, ranked_vocabulary, matrix_vocabulary, reference_vocabulary)
#   assert_ledger_closes(ledger, label)
#
# Source from the project root, AFTER config.R:
#   source("02_analysis/helpers/symbol_alias.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

# The resolutions a candidate pair can carry. Every value except `accepted` withholds the
# pair, and `flagged_for_review` is the one that is withheld by human decision rather than
# by the automated guard.
ALIAS_RESOLUTIONS <- c(
  "accepted",
  "flagged_for_review",
  "rejected_symbol_belongs_to_another_gene",
  "rejected_multiple_aliases_in_vocabulary",
  "rejected_reference_symbol_ambiguous_in_org_db")

.alias_map_empty <- function() {
  tibble(reference_symbol = character(), matrix_symbol = character(),
         entrez_id = character(), n_aliases_in_vocabulary = integer(),
         resolution = character())
}

#' Candidate pairs mapping a reference symbol onto the vocabulary the data carries.
#'
#' A pair is a candidate when the reference symbol is absent from `matrix_vocabulary`,
#' resolves to exactly one Entrez id, and exactly one alias of that same Entrez id is
#' present in the vocabulary. Candidates are then filtered by the ownership guard.
#'
#' @param reference_symbols character  symbols as the reference set ships them.
#' @param matrix_vocabulary character  the target vocabulary to resolve into.
#' @param db  an org.*.eg.db object; org.Mm.eg.db for this compartment.
#' @param flagged_pairs character  "REF->MATRIX" strings withheld for human review.
#' @return tibble(reference_symbol, matrix_symbol, entrez_id, n_aliases_in_vocabulary,
#'   resolution), one row per candidate pair, carrying attributes `summary` (the
#'   no-candidate counts) and `n_rejected_ambiguous`.
build_alias_map <- function(reference_symbols, matrix_vocabulary,
                            db = NULL, flagged_pairs = character()) {
  reference_symbols <- unique(reference_symbols[!is.na(reference_symbols) &
                                                  nzchar(reference_symbols)])
  matrix_vocabulary <- unique(matrix_vocabulary)
  missing_symbols <- sort(setdiff(reference_symbols, matrix_vocabulary))
  empty <- .alias_map_empty()
  summ <- list(n_reference_symbols = length(reference_symbols),
               n_absent_from_vocabulary = length(missing_symbols),
               n_not_in_org_db = NA_integer_,
               n_no_alias_in_vocabulary = NA_integer_,
               n_candidates = 0L, n_accepted = 0L, n_flagged = 0L, n_rejected = 0L)
  finish <- function(map) {
    summ$n_candidates <- nrow(map)
    summ$n_accepted   <- sum(map$resolution == "accepted")
    summ$n_flagged    <- sum(map$resolution == "flagged_for_review")
    summ$n_rejected   <- sum(!map$resolution %in% c("accepted", "flagged_for_review"))
    structure(map, summary = summ,
              n_rejected_ambiguous = sum(
                map$resolution == "rejected_symbol_belongs_to_another_gene"))
  }

  if (is.null(db)) {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE) ||
        !requireNamespace("AnnotationDbi", quietly = TRUE)) {
      message("  org.Mm.eg.db unavailable, alias recovery skipped and 0 symbols recovered.")
      return(finish(empty))
    }
    db <- org.Mm.eg.db::org.Mm.eg.db
  }
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(finish(empty))
  if (!length(missing_symbols)) return(finish(empty))

  # AnnotationDbi::select() throws outright when NONE of its keys is valid for the
  # keytype, so every call is prefiltered against the live keyspace and returns early when
  # nothing survives. With a whole collection's ~700 unmatched symbols this never fires;
  # on a single small gene set it does, and a set whose only unmatched member is already a
  # legacy name would take the ownership guard down rather than report a clean zero.
  symbol_keys <- AnnotationDbi::keys(db, keytype = "SYMBOL")
  sel <- function(keys, keytype, columns) {
    keys <- unique(keys[!is.na(keys) & nzchar(keys)])
    if (identical(keytype, "SYMBOL")) keys <- intersect(keys, symbol_keys)
    if (!length(keys)) return(NULL)
    suppressMessages(AnnotationDbi::select(db, keys = keys, keytype = keytype,
                                           columns = columns))
  }

  eg <- sel(missing_symbols, "SYMBOL", "ENTREZID")
  summ$n_not_in_org_db <- length(setdiff(missing_symbols, symbol_keys))
  if (is.null(eg)) return(finish(empty))
  eg <- eg[!is.na(eg$ENTREZID), , drop = FALSE]
  # A reference symbol carrying more than one Entrez id names more than one gene here, so
  # it is withheld — but only reported once a candidate for it actually exists.
  ambiguous <- unique(eg$SYMBOL[duplicated(eg$SYMBOL)])
  eg <- eg[!(eg$SYMBOL %in% ambiguous), , drop = FALSE]
  if (!nrow(eg)) return(finish(empty))

  al <- sel(unique(eg$ENTREZID), "ENTREZID", "ALIAS")
  if (is.null(al)) return(finish(empty))
  al <- al[al$ALIAS %in% matrix_vocabulary, , drop = FALSE]
  if (!nrow(al)) return(finish(empty))
  # Two aliases of one gene present in the vocabulary leaves no unique target, so the pair
  # is withheld with both names recorded rather than silently returned as NA.
  cand <- al %>% group_by(.data$ENTREZID) %>%
    summarise(matrix_symbol = paste(sort(unique(.data$ALIAS)), collapse = "/"),
              n_aliases_in_vocabulary = n_distinct(.data$ALIAS), .groups = "drop")
  hits <- eg %>% inner_join(cand, by = "ENTREZID") %>%
    transmute(reference_symbol = .data$SYMBOL, matrix_symbol = .data$matrix_symbol,
              entrez_id = .data$ENTREZID,
              n_aliases_in_vocabulary = as.integer(.data$n_aliases_in_vocabulary))
  summ$n_no_alias_in_vocabulary <- length(setdiff(eg$SYMBOL, hits$reference_symbol))
  if (!nrow(hits)) return(finish(empty))

  multi <- hits %>% filter(.data$n_aliases_in_vocabulary > 1L) %>%
    mutate(resolution = "rejected_multiple_aliases_in_vocabulary")
  hits <- hits %>% filter(.data$n_aliases_in_vocabulary == 1L)
  if (!nrow(hits)) return(finish(multi))

  # The ownership guard: is the candidate the official symbol of some other gene?
  own <- sel(unique(hits$matrix_symbol), "SYMBOL", "ENTREZID")
  taken <- rep(FALSE, nrow(hits))
  if (!is.null(own)) {
    own <- own[!is.na(own$ENTREZID), , drop = FALSE]
    owner <- setNames(own$ENTREZID, own$SYMBOL)
    taken <- !is.na(owner[hits$matrix_symbol]) &
      owner[hits$matrix_symbol] != hits$entrez_id
    taken[is.na(taken)] <- FALSE
  }
  hits$resolution <- ifelse(taken, "rejected_symbol_belongs_to_another_gene", "accepted")
  # Pairs a human has to decide on are withheld here rather than downstream, so the
  # exclusion travels with the map and is visible in every consumer's ledger.
  flagged <- paste0(hits$reference_symbol, "->", hits$matrix_symbol) %in% flagged_pairs
  hits$resolution[flagged & hits$resolution == "accepted"] <- "flagged_for_review"

  amb <- if (length(ambiguous))
    tibble(reference_symbol = intersect(ambiguous, missing_symbols),
           matrix_symbol = NA_character_, entrez_id = NA_character_,
           n_aliases_in_vocabulary = NA_integer_,
           resolution = "rejected_reference_symbol_ambiguous_in_org_db")
  else .alias_map_empty()

  finish(bind_rows(hits, multi, amb) %>% arrange(.data$reference_symbol))
}

#' The reference_symbol -> matrix_symbol lookup an alias map licenses.
#'
#' Only `accepted` rows are returned; `flagged_for_review` is withheld by construction so
#' no consumer can apply a pair a human has not signed off on.
accepted_pairs <- function(alias_map) {
  if (is.null(alias_map) || !nrow(alias_map)) return(setNames(character(), character()))
  a <- alias_map[alias_map$resolution == "accepted", , drop = FALSE]
  setNames(a$matrix_symbol, a$reference_symbol)
}

#' Resolve every set in an fgsea-shaped `pathways` list into one vocabulary.
#'
#' @param sets named list of character vectors.
#' @param matrix_vocabulary character  the vocabulary the sets are being matched against.
#' @param alias_map tibble from build_alias_map(), or a named reference->matrix vector.
#' @param replace logical  substitute the reference symbol with the matrix symbol, or keep
#'   both spellings (FALSE, the default).
#'
#'   Why BOTH spellings are kept, which took a measurement to settle. An accepted reference
#'   symbol is absent from the vocabulary by construction, so keeping it matches nothing and
#'   is inert for fgsea — the published `set_size` comes from fgsea's post-intersection
#'   `size`, never from the nominal length. The tempting alternative, substituting, keeps the
#'   nominal size equal to the curated size and reads cleaner. It is also a REGRESSION:
#'   04_gsea_set_prep.R's filter_by_size() runs on the nominal length, and 12 sets across
#'   Reactome/GO_BP/GO_MF/GO_CC carry BOTH vintages of one gene (the MHC class Ib,
#'   β2-microglobulin and TAP-binding families above all), so substitution shrinks them by
#'   one, they fall under gsea_min_size, and the fix DROPS 12 sets that were previously
#'   tested. Keeping both spellings makes the nominal size monotone — it can only grow — so
#'   the size filter can only ever admit a set, never lose one. Sets that gained a genuinely
#'   inflated nominal size are still filtered on the true matched size by fgsea itself, and
#'   `n_alias_collapsed` in the ledger names every set where the two vintages met.
#'
#'   For the same reason the input vector is APPENDED to rather than re-`unique()`d. Some
#'   msigdbr sets ship one symbol twice under one gs_name, and both filter_by_size() here and
#'   05_gsea_msigdb_run.R's defensive re-filter count the raw length — so de-duplicating
#'   quietly shrinks those 12 sets by one, pushes them under gsea_min_size, and costs the
#'   pipeline sets that have nothing to do with aliases. The multiplicity the reference shipped
#'   is left exactly as it was; fgsea computes `size` over unique members regardless, so
#'   nothing statistical rides on it.
#' @return list(sets = resolved, applied = per set x pair,
#'   collapsed = sets that gained fewer genes than pairs applied,
#'   many_to_one = two reference symbols resolving onto one matrix symbol)
resolve_sets <- function(sets, matrix_vocabulary, alias_map, replace = FALSE) {
  pairs <- if (is.data.frame(alias_map)) accepted_pairs(alias_map) else alias_map
  matrix_vocabulary <- unique(matrix_vocabulary)
  applied <- list(); collapsed <- list(); many_to_one <- list()

  out <- lapply(names(sets), function(nm) {
    raw <- sets[[nm]]
    g <- unique(raw)
    hit <- names(pairs)[names(pairs) %in% g]
    hit <- hit[!hit %in% matrix_vocabulary]           # a set may carry both vintages
    tgt <- unname(pairs[hit])
    hit <- hit[tgt %in% matrix_vocabulary]
    tgt <- tgt[tgt %in% matrix_vocabulary]
    if (length(hit)) {
      applied[[nm]] <<- tibble(gene_set = nm, reference_symbol = hit, matrix_symbol = tgt)
      # Two reference symbols landing on one matrix symbol merges two set members into one
      # measurement. It is reported, never silently merged.
      dup_tgt <- unique(tgt[duplicated(tgt)])
      if (length(dup_tgt))
        many_to_one[[nm]] <<- tibble(
          gene_set = nm, matrix_symbol = dup_tgt,
          reference_symbols = vapply(dup_tgt,
                                     function(s) paste(sort(hit[tgt == s]), collapse = "/"),
                                     character(1)))
    }
    res <- if (replace) unique(c(setdiff(g, hit), tgt))
           else c(raw, setdiff(tgt, g))   # append only; the reference's own multiplicity stands
    # A set that does not grow by the number of pairs applied already carried the matrix
    # symbol under its other vintage, or two pairs collapsed onto one target. Under
    # replace = TRUE the reference spelling is also removed, so the yardstick is the count
    # of MATCHING members, not the raw length.
    gained <- length(intersect(res, matrix_vocabulary)) -
      length(intersect(g, matrix_vocabulary))
    if (length(hit) && gained != length(hit))
      collapsed[[nm]] <<- tibble(gene_set = nm,
                                 n_matching_before = length(intersect(g, matrix_vocabulary)),
                                 n_matching_after = length(intersect(res, matrix_vocabulary)),
                                 n_pairs_applied = length(hit),
                                 n_collapsed = length(hit) - gained)
    res
  })
  names(out) <- names(sets)
  list(sets = out, applied = bind_rows(applied), collapsed = bind_rows(collapsed),
       many_to_one = bind_rows(many_to_one))
}

#' The per-set symbol ledger: one bucket per cause, and they are kept apart.
#'
#' Conflating vocabulary loss with never-detected and with expression-filtered is the whole
#' defect, so each unmatched gene lands in exactly one bucket and the buckets close against
#' the set's size. The vocabulary layers in this compartment, and what each one states:
#'   reference_vocabulary  the GENCODE vM25 feature list. NOT persisted in this
#'                         compartment — the collaborators delivered a CPM table, not a
#'                         quantification against a tracked GTF — so it is normally
#'                         absent, `reference_vocabulary_available` is FALSE, and
#'                         n_absent_from_reference stays 0 rather than being fabricated.
#'   matrix_vocabulary     the delivered CPM table's gene_name column, i.e. every symbol
#'                         the experiment actually quantified. Below it is a DETECTION
#'                         statement about this experiment, not a vocabulary one.
#'   ranked_vocabulary     the modelled DE universe (gene_universe.txt), post
#'                         duplicate-symbol collapse and constant-row drop. Below it is a
#'                         statement about what the model could test.
#'
#' @return tibble, one row per set, closing on
#'   n_unique_set_genes == n_matched + n_matched_via_alias + n_alias_flagged_for_review +
#'     n_alias_rejected_ambiguous + n_expression_filtered + n_below_detection +
#'     n_absent_from_reference
#'
#' The buckets count REFERENCE genes, so they close against the set as the reference ships
#' it. `set_size_resolved` counts MATRIX symbols and is therefore smaller than
#' n_matched + n_matched_via_alias wherever pairs collapsed; `n_alias_collapsed` carries
#' that difference so the two never have to be reconciled by hand.
symbol_ledger <- function(sets, alias_map, ranked_vocabulary, matrix_vocabulary,
                          reference_vocabulary = NULL) {
  pairs <- accepted_pairs(alias_map)
  flag  <- alias_map[alias_map$resolution == "flagged_for_review", , drop = FALSE]
  rej   <- alias_map[alias_map$resolution ==
                       "rejected_symbol_belongs_to_another_gene", , drop = FALSE]
  ranked_vocabulary <- unique(ranked_vocabulary)
  matrix_vocabulary <- unique(matrix_vocabulary)
  have_ref <- !is.null(reference_vocabulary) && length(reference_vocabulary) > 0
  ref_vocab <- if (have_ref) unique(reference_vocabulary) else character()

  bind_rows(lapply(names(sets), function(nm) {
    g <- unique(sets[[nm]])
    matched <- intersect(g, ranked_vocabulary)
    rest <- setdiff(g, matched)
    # Precedence: a gene recoverable through an accepted pair is a vocabulary result and is
    # counted as such before any statement about detection or power is made.
    via <- rest[rest %in% names(pairs) & unname(pairs[rest]) %in% ranked_vocabulary]
    rest <- setdiff(rest, via)
    flagged <- rest[rest %in% flag$reference_symbol]
    rest <- setdiff(rest, flagged)
    ambig <- rest[rest %in% rej$reference_symbol]
    rest <- setdiff(rest, ambig)
    # The remaining buckets are read off the symbol the DATA would carry, so a retired
    # reference name whose matrix twin was dropped before modelling is reported as
    # expression-filtered rather than as never detected. Bucketing on the reference name
    # there would put a power statement in a detection bucket.
    eff <- ifelse(rest %in% names(pairs), unname(pairs[rest]), rest)
    expr_filtered <- rest[eff %in% matrix_vocabulary]
    rest <- setdiff(rest, expr_filtered)
    eff <- ifelse(rest %in% names(pairs), unname(pairs[rest]), rest)
    # With no vM25 feature list to appeal to, everything left is "not in the delivered CPM
    # table", and the delivered table is itself the detection-filtered quantification —
    # so the honest terminal bucket is below-detection, and absent-from-reference is a
    # claim this compartment cannot make. When a feature list IS supplied the two split.
    absent_ref <- if (have_ref) setdiff(rest, rest[eff %in% ref_vocab]) else character()
    below <- setdiff(rest, absent_ref)
    tibble(
      gene_set = nm, n_unique_set_genes = length(g),
      n_matched = length(matched),
      n_matched_via_alias = length(via),
      n_alias_flagged_for_review = length(flagged),
      n_alias_rejected_ambiguous = length(ambig),
      n_expression_filtered = length(expr_filtered),
      n_below_detection = length(below),
      n_absent_from_reference = length(absent_ref),
      # The true resolved size, de-duplicated. Two things collapse here and both are
      # reported rather than absorbed: a set carrying both vintages of one gene, and
      # several reference paralogs whose current names all resolve onto one matrix row.
      # Either way the set grows by fewer genes than pairs applied, and
      # `n_alias_collapsed` is that shortfall.
      set_size_resolved = length(unique(c(matched, unname(pairs[via])))),
      n_alias_collapsed = length(via) -
        (length(unique(c(matched, unname(pairs[via])))) - length(matched)),
      reference_vocabulary_available = have_ref,
      alias_pairs_applied = paste(sort(paste0(via, "->", unname(pairs[via]))),
                                  collapse = "/"),
      alias_pairs_flagged = paste(sort(paste0(
        flagged, "->", flag$matrix_symbol[match(flagged, flag$reference_symbol)])),
        collapse = "/"),
      alias_pairs_rejected = paste(sort(paste0(
        ambig, "->", rej$matrix_symbol[match(ambig, rej$reference_symbol)])),
        collapse = "/"))
  }))
}

#' Hard closure check on a ledger, asserted in-script rather than trusted in review.
assert_ledger_closes <- function(ledger, label = "symbol ledger") {
  s <- with(ledger, n_matched + n_matched_via_alias + n_alias_flagged_for_review +
              n_alias_rejected_ambiguous + n_expression_filtered + n_below_detection +
              n_absent_from_reference)
  bad <- which(s != ledger$n_unique_set_genes)
  if (length(bad))
    stop(sprintf(paste0("[symbol_alias] %s does not close for %d set(s), e.g. %s ",
                        "(%d genes, buckets sum to %d). Every unmatched gene must land ",
                        "in exactly one cause bucket."),
                 label, length(bad), ledger$gene_set[bad[1]],
                 ledger$n_unique_set_genes[bad[1]], s[bad[1]]))
  invisible(TRUE)
}

message("[symbol_alias] loaded (reference MGI symbols -> this matrix's vocabulary; lazy dep: org.Mm.eg.db).")
