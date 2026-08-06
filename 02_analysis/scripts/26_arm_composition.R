#!/usr/bin/env Rscript
# 26_arm_composition.R -- COMPUTE
# =============================================================================
# What fraction of each projected arm each pathway accounts for
# (stage 16_arm_composition).
#
# Why this stage exists
#   Two decompositions of these arms already exist, and each answers a different question
#   from the one a composition figure asks.
#
#   The nine-lens membership decomposition is close to a partition — 100 memberships over
#   67 of WT_heat_up's 199 genes, at most 4 programs per gene — with a panel picked by
#   whoever was asking, and 132 genes falling outside all of it.
#
#   15_go_decomposition removes the hand-picking and stops short of a partition. At term
#   level its 386 enriched GO BP terms give 3,621 gene-term memberships over a union of 159
#   genes, with 9 genes in exactly one term. Collapsing to its 35 blocks gives 972
#   block-memberships over the same 159 genes, 21 genes in exactly one block, 26 in twelve
#   or more. And a random depth-matched gene set calls the five terms that would headline a
#   figure significant 28% to 63% of the time.
#
#   This stage keeps the partition property and drops the hand-picking. Sets come from five
#   frozen MSigDB collections. Selection is by depth-matched permutation. Redundant sets are
#   pruned on how much of the ARM they share. Attribution is written two ways, so a reader
#   can see whether the composition is determined at all.
#
# Its relation to 15_go_decomposition
#   That stage tests a propagated gene-to-term map built from org.Hs.eg.db with an explicit
#   IEA switch. MSigDB's C5 GO:BP is a differently curated and differently filtered product.
#   The two GO_BP results are different tests of related content, and section 9 reports
#   their concordance.
#
# Reads
#   03_results/human_projection/signatures/... (the four arms + the ranked-list background)
#   ../human_treg_arthritis/03_results/objects/14_genesets.rds  (five frozen collections)
#   03_results/15_go_decomposition/tables/_overview/go_enriched_terms.csv  (concordance only)
#
# Writes 03_results/16_arm_composition/tables/_overview/
#   composition_ora_terms.csv       every enriched set, all arms x collections, null columns
#   composition_selected.csv        which sets entered, and why each was dropped
#   composition_gene_assignment.csv per arm x variant x gene: category and weight
#   composition_shares.csv          the bar substrate: arm x variant x accounting x category
#   composition_remainder.csv       the residue, split by reason
#   composition_hypoxia_sources.csv the anchored variant's sets and their mutual overlap
#   composition_null_summary.csv    per arm x collection, including the coverage null and
#                                   a not_tested_reason for any pair with no test
#   composition_permutation_floor.csv  how many selected sets sit on the p_matched floor,
#                                   the resolution limit of the primary key
#   composition_go_concordance.csv  this stage's GO_BP against 15_go_decomposition's
#   composition_provenance.csv
#   plus 03_results/objects/26_arm_composition.rds for the viz script.
#
# Run from project root (~30 min):
#   Rscript 02_analysis/scripts/26_arm_composition.R
# =============================================================================

source("02_analysis/config/config.R")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(Matrix)
  library(dplyr)
})
options(stringsAsFactors = FALSE)

STAGE   <- "16_arm_composition"
SCRIPT  <- "02_analysis/scripts/26_arm_composition.R"
TBL_OVW <- file.path(stage_dir(STAGE, "tables"), "_overview")
dir.create(TBL_OVW, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OBJECTS, recursive = TRUE, showWarnings = FALSE)

AC <- YAML_CONFIG$arm_composition
GO <- YAML_CONFIG$go_decomposition
stopifnot(!is.null(AC), !is.null(GO))

COLLECTIONS  <- as.character(unlist(AC$collections))
P_CUT        <- as.numeric(AC$p_cutoff              %||% 0.05)
Q_CUT        <- as.numeric(AC$q_cutoff              %||% 0.2)
PADJ_METHOD  <- as.character(AC$p_adjust_method      %||% "BH")
MIN_GS       <- as.integer(AC$min_gs_size            %||% 10L)
MAX_GS       <- as.integer(AC$max_gs_size            %||% 500L)
PM_CUT       <- as.numeric(AC$p_matched_cutoff       %||% 0.05)
MAX_RECUR    <- as.numeric(AC$max_null_recurrence    %||% 0.25)
PRUNE_J      <- as.numeric(AC$prune_hit_jaccard      %||% 0.5)
TOP_N        <- as.integer(AC$top_n_per_collection   %||% 10L)
ACCOUNTINGS  <- as.character(unlist(AC$accountings))
N_NULL       <- as.integer(AC$n_null                 %||% 2000L)
BAND_Q       <- as.numeric(unlist(AC$depth_band_quantiles))
SEED         <- as.integer(AC$seed                   %||% 20260731L)
HYPOXIA      <- AC$hypoxia_anchor

RESIDUAL <- "not claimed by any selected set"

message("=================================================================")
message("26_arm_composition: pathway shares of each arm (", STAGE, ")")
message("  collections ", paste(COLLECTIONS, collapse = "/"),
        "   top ", TOP_N, " per collection")
message("  selection: p_matched < ", PM_CUT, ", recurrence <= ", MAX_RECUR,
        ", hit-Jaccard prune at ", PRUNE_J)
message("  nulls ", N_NULL, " per arm per collection   seed ", SEED)
message("=================================================================")

## ---------------------------------------------------------------------------
## 0. Helpers
## ---------------------------------------------------------------------------

read_symbols <- function(p) unique(trimws(readLines(p, warn = FALSE)))

# p.adjust()'s default n counts NAs, which would deflate every adjusted p in a null
# replicate because most sets have no query hit and are not in the tested family.
bh_col <- function(p) {
  ok <- !is.na(p)
  out <- rep(NA_real_, length(p))
  out[ok] <- stats::p.adjust(p[ok], method = PADJ_METHOD)
  out
}

collapse_genes <- function(g) if (!length(g)) "" else paste(sort(unique(g)), collapse = "/")

# The variant vocabulary is defined ONCE and derived everywhere. It was previously a literal in
# two places, the two drifted under a rename, and the published `composition_selected.csv`
# disagreed with every other table in the directory. A reader filtering on the wrong label gets
# zero rows and no error, so the check before writing is not decoration.
variant_label <- function(anchored) if (anchored) "hypoxia_anchored" else "unpinned"
VARIANT_LEVELS <- c(variant_label(FALSE), variant_label(TRUE))

jaccard <- function(a, b) {
  u <- length(union(a, b))
  if (!u) return(0)
  length(intersect(a, b)) / u
}

## ---------------------------------------------------------------------------
## 1. Inputs: arms, background, frozen collections
## ---------------------------------------------------------------------------

BG_PATH <- file.path(PROJECT_ROOT, GO$background)
stopifnot(file.exists(BG_PATH))
BACKGROUND <- unique(utils::read.delim(
  BG_PATH, header = FALSE, col.names = c("gene", "score"))$gene)
message("[1] background: ", length(BACKGROUND), " human symbols")

ARMS <- lapply(GO$arms, function(a) {
  p <- file.path(PROJECT_ROOT, a$file)
  stopifnot(file.exists(p))
  g <- read_symbols(p)
  list(name = as.character(a$name), role = as.character(a$role %||% "comparator"),
       genes = g, n_nominal = length(g),
       n_in_background = length(intersect(g, BACKGROUND)))
})
names(ARMS) <- vapply(ARMS, function(a) a$name, character(1))
ARM_ORDER <- names(ARMS)
for (a in ARMS)
  message("[1] arm ", a$name, ": ", a$n_nominal, " nominal, ",
          a$n_in_background, " in background")

GS_PATH <- normalizePath(file.path(PROJECT_ROOT, AC$genesets_object), mustWork = TRUE)
GS_MD5  <- unname(tools::md5sum(GS_PATH))
GS_RAW  <- readRDS(GS_PATH)
missing_coll <- setdiff(COLLECTIONS, names(GS_RAW))
if (length(missing_coll))
  stop("collections absent from the frozen object: ", paste(missing_coll, collapse = ", "))
SETS <- lapply(COLLECTIONS, function(cl) GS_RAW[[cl]]$sets)
names(SETS) <- COLLECTIONS
SET_SOURCE <- vapply(COLLECTIONS,
                     function(cl) as.character(GS_RAW[[cl]]$source %||% NA_character_),
                     character(1))
for (cl in COLLECTIONS)
  message("[1] collection ", cl, ": ", length(SETS[[cl]]), " sets  (", SET_SOURCE[[cl]], ")")

## ---------------------------------------------------------------------------
## 2. One (term, gene) table per collection, restricted to the background
##
## Restriction matters: a set's size in the hypergeometric must be its size IN the
## universe being tested, which differs from its nominal MSigDB size, or every p is computed against
## a universe the genes are not in.
## ---------------------------------------------------------------------------

T2G <- lapply(COLLECTIONS, function(cl) {
  s <- lapply(SETS[[cl]], function(g) intersect(unique(g), BACKGROUND))
  s <- s[lengths(s) > 0]
  data.frame(term = rep(names(s), lengths(s)), gene = unlist(s, use.names = FALSE))
})
names(T2G) <- COLLECTIONS

SET_SIZES <- lapply(COLLECTIONS, function(cl) {
  tb <- table(T2G[[cl]]$term)
  data.frame(collection = cl, term = names(tb), size_in_background = as.integer(tb))
})
SET_SIZES <- dplyr::bind_rows(SET_SIZES)
for (cl in COLLECTIONS) {
  z <- SET_SIZES[SET_SIZES$collection == cl, ]
  message("[2] ", cl, ": ", nrow(z), " sets in background, ",
          sum(z$size_in_background >= MIN_GS & z$size_in_background <= MAX_GS),
          " testable at ", MIN_GS, "-", MAX_GS)
}

## ---------------------------------------------------------------------------
## 3. Over-representation: every arm x collection
## ---------------------------------------------------------------------------

ORA <- list()
for (arm in ARM_ORDER) {
  for (cl in COLLECTIONS) {
    e <- suppressMessages(clusterProfiler::enricher(
      gene = ARMS[[arm]]$genes, universe = BACKGROUND,
      pvalueCutoff = P_CUT, pAdjustMethod = PADJ_METHOD, qvalueCutoff = Q_CUT,
      minGSSize = MIN_GS, maxGSSize = MAX_GS, TERM2GENE = T2G[[cl]]))
    df <- if (is.null(e)) data.frame() else as.data.frame(e)
    if (nrow(df)) {
      df$arm <- arm; df$collection <- cl
      df$n_query_in_background <- ARMS[[arm]]$n_in_background
      df$hit_genes <- df$geneID
      df$n_hits <- df$Count
    }
    ORA[[paste(arm, cl, sep = "|")]] <- df
    message("[3] ", arm, " x ", cl, ": ", nrow(df), " enriched")
  }
}

## ---------------------------------------------------------------------------
## 4. Annotation-depth-matched permutation null, per arm x collection
##
## The engine is 15_go_decomposition's, with one addition: per replicate it also records
## how many of the DRAWN genes fall in at least one set reaching significance. That is the
## null the coverage and share numbers need, and the earlier stage does not have it, so
## "the ontology covers 80% of the arm" was unfalsifiable there.
##
## Both stages must change together if this engine changes.
## ---------------------------------------------------------------------------

null_engine <- function(t2g, query, n_null, seed) {
  genes <- sort(unique(t2g$gene))
  terms <- sort(unique(t2g$term))
  M0 <- Matrix::sparseMatrix(i = match(t2g$term, terms), j = match(t2g$gene, genes),
                             x = 1L, dims = c(length(terms), length(genes)))
  depth <- stats::setNames(Matrix::colSums(M0), genes)
  sz0 <- Matrix::rowSums(M0)
  keep <- sz0 >= MIN_GS & sz0 <= MAX_GS
  M <- M0[keep, , drop = FALSE]; sz <- sz0[keep]; terms <- terms[keep]
  NG <- length(genes)
  qa <- intersect(query, genes)
  n <- length(qa)
  if (n < 2L || !length(terms)) return(NULL)

  brk <- unique(c(-1, stats::quantile(depth, BAND_Q, names = FALSE), Inf))
  qband <- cut(depth[qa], brk)
  pool <- setdiff(genes, query)
  pband <- cut(depth[pool], brk)
  need <- table(qband); lev <- levels(qband)
  pool_by_band <- split(pool, pband)

  draw_matched <- function() {
    used <- character(0); borrowed <- 0L
    for (bi in seq_along(lev)) {
      k <- as.integer(need[[bi]]); if (k == 0L) next
      p <- setdiff(pool_by_band[[lev[bi]]] %||% character(0), used)
      if (length(p) >= k) {
        used <- c(used, sample(p, k))
      } else {
        used <- c(used, p)
        short <- k - length(p); borrowed <- borrowed + short
        for (bj in order(abs(seq_along(lev) - bi))) {
          if (short == 0L) break
          if (bj == bi) next
          extra <- setdiff(pool_by_band[[lev[bj]]] %||% character(0), used)
          take <- min(short, length(extra))
          if (take > 0L) { used <- c(used, sample(extra, take)); short <- short - take }
        }
      }
    }
    stopifnot(length(used) == n, !anyDuplicated(used))
    list(genes = used, borrowed = borrowed)
  }

  set.seed(seed)
  matched <- lapply(seq_len(n_null), function(i) draw_matched())
  uniform <- lapply(seq_len(n_null),
                    function(i) list(genes = sample(pool, n), borrowed = 0L))

  score <- function(sets) {
    idx <- lapply(sets, function(s) match(s$genes, genes))
    B <- Matrix::sparseMatrix(i = unlist(idx), j = rep(seq_along(idx), lengths(idx)),
                              x = 1L, dims = c(NG, length(idx)))
    K <- as.matrix(M %*% B)
    P <- stats::phyper(K - 1L, sz, NG - sz, n, lower.tail = FALSE)
    P[K == 0L] <- NA_real_
    Q <- apply(P, 2, bh_col)
    # THE COVERAGE NULL: per replicate, how many drawn genes are in >=1 significant set.
    covered <- vapply(seq_along(idx), function(j) {
      sig <- which(!is.na(Q[, j]) & Q[, j] < P_CUT)
      if (!length(sig)) return(0L)
      as.integer(sum(Matrix::colSums(M[sig, idx[[j]], drop = FALSE]) > 0))
    }, integer(1))
    list(K = K, Q = Q, n_signif = colSums(Q < P_CUT, na.rm = TRUE),
         family = colSums(K > 0L), covered = covered,
         borrowed = vapply(sets, function(s) s$borrowed, integer(1)))
  }
  rm_ <- score(matched); ru <- score(uniform)

  kobs <- as.integer(M %*% Matrix::sparseMatrix(
    i = match(qa, genes), j = rep(1L, n), x = 1L, dims = c(NG, 1L)))
  pobs <- stats::phyper(kobs - 1L, sz, NG - sz, n, lower.tail = FALSE)
  pobs[kobs == 0L] <- NA_real_
  qobs <- bh_col(pobs)
  sig_obs <- which(!is.na(qobs) & qobs < P_CUT)
  covered_obs <- if (!length(sig_obs)) 0L else as.integer(
    sum(Matrix::colSums(M[sig_obs, match(qa, genes), drop = FALSE]) > 0))

  per_term <- data.frame(
    term = terms, size_tested = as.integer(sz), k_obs = kobs,
    p_obs_recomputed = pobs, q_obs_recomputed = qobs,
    k_null_matched_mean = rowMeans(rm_$K),
    p_matched = (1 + rowSums(rm_$K >= kobs)) / (n_null + 1),
    frac_matched_reaching_q = rowMeans(rm_$Q < P_CUT, na.rm = TRUE),
    p_uniform = (1 + rowSums(ru$K >= kobs)) / (n_null + 1))
  per_term$frac_matched_reaching_q[is.nan(per_term$frac_matched_reaching_q)] <- 0

  list(per_term = per_term, n_query_annotated = n, n_universe = NG,
       n_terms_testable = length(terms), covered_obs = covered_obs,
       summary = data.frame(
         null = c("depth_matched", "uniform"), n_replicates = n_null,
         mean_family_size = c(mean(rm_$family), mean(ru$family)),
         median_n_signif  = c(stats::median(rm_$n_signif), stats::median(ru$n_signif)),
         q95_n_signif     = c(stats::quantile(rm_$n_signif, .95, names = FALSE),
                              stats::quantile(ru$n_signif, .95, names = FALSE)),
         n_observed_signif = length(sig_obs),
         median_covered   = c(stats::median(rm_$covered), stats::median(ru$covered)),
         q95_covered      = c(stats::quantile(rm_$covered, .95, names = FALSE),
                              stats::quantile(ru$covered, .95, names = FALSE)),
         covered_obs      = covered_obs,
         p_covered        = c((1 + sum(rm_$covered >= covered_obs)) / (n_null + 1),
                              (1 + sum(ru$covered >= covered_obs)) / (n_null + 1)),
         max_genes_borrowed = c(max(rm_$borrowed), 0)))
}

NULLS <- list()
NOT_TESTED <- list()
for (arm in ARM_ORDER) {
  for (cl in COLLECTIONS) {
    key <- paste(arm, cl, sep = "|")
    message("[4] null: ", arm, " x ", cl, " x ", N_NULL)
    NULLS[[key]] <- null_engine(T2G[[cl]], ARMS[[arm]]$genes, N_NULL, SEED)
    if (is.null(NULLS[[key]])) {
      # A silent NULL leaves an arm x collection cell simply absent from the summary,
      # which reads as an oversight. Record which branch.
      gset <- unique(T2G[[cl]]$gene)
      tb <- table(T2G[[cl]]$term)
      n_ok <- sum(tb >= MIN_GS & tb <= MAX_GS)
      n_q <- length(intersect(ARMS[[arm]]$genes, gset))
      # `null` carries an explicit level, over NA. A consumer subscripting with
      # base-R `df[df$null == "depth_matched", ]` gets an all-NA ROW back for every NA in
      # the comparison, so a single NA here silently lengthens every downstream lookup.
      # That is not hypothetical: it broke a caption in 27_arm_composition_viz.R.
      NOT_TESTED[[key]] <- data.frame(
        null = "not_tested", n_replicates = NA_integer_,
        arm = arm, collection = cl, n_query_annotated = n_q,
        n_terms_testable = n_ok,
        not_tested_reason = if (n_q < 2L)
          sprintf("only %d arm gene(s) annotated in this collection, need 2", n_q)
        else sprintf("no set in the %d-%d size window", MIN_GS, MAX_GS))
      message("      not tested: ", NOT_TESTED[[key]]$not_tested_reason)
    }
  }
}

# Seam: the matrix recomputation must agree with enricher() set by set.
SEAM <- dplyr::bind_rows(lapply(names(NULLS), function(nm) {
  ne <- NULLS[[nm]]; ob <- ORA[[nm]]
  if (is.null(ne) || !nrow(ob))
    return(data.frame(key = nm, n_compared = 0L, max_abs_dk = NA_real_,
                      max_abs_dq = NA_real_))
  j <- dplyr::inner_join(ob[, c("ID", "Count", "p.adjust")], ne$per_term,
                         by = c("ID" = "term"))
  data.frame(key = nm, n_compared = nrow(j),
             max_abs_dk = max(abs(j$Count - j$k_obs)),
             max_abs_dq = max(abs(j$p.adjust - j$q_obs_recomputed)))
}))
print(SEAM)
stopifnot(all(is.na(SEAM$max_abs_dk) | SEAM$max_abs_dk == 0))

# Attach the null columns to the enriched tables.
for (nm in names(NULLS)) {
  ne <- NULLS[[nm]]
  if (is.null(ne) || !nrow(ORA[[nm]])) next
  ORA[[nm]] <- dplyr::left_join(
    ORA[[nm]],
    ne$per_term[, c("term", "size_tested", "k_null_matched_mean", "p_matched",
                    "frac_matched_reaching_q", "p_uniform")],
    by = c("ID" = "term"))
}
ENRICHED <- dplyr::bind_rows(ORA)

## ---------------------------------------------------------------------------
## 5. Selection, and a full audit trail of what was dropped and why
## ---------------------------------------------------------------------------

hits_of <- function(row) strsplit(row$hit_genes, "/", fixed = TRUE)[[1]]

# Every candidate set for an arm: the enriched ones, plus -- in the anchored variant --
# the pinned hypoxia sets whether or not they enriched.
#
# This is the whole point of the anchored variant. Pinning only among enriched sets would
# let hypoxia appear only when hypoxia had already won on its own, which answers nothing.
# For WT_heat_up exactly one of the four hypoxia sets enriches; the other three are tested
# and null, and their share is what a reader wants to see. A pinned set that never enriched
# carries NA significance and is labelled, so its bar cannot be read as a result.
candidates_for <- function(arm, anchored) {
  cand <- ENRICHED[ENRICHED$arm == arm, , drop = FALSE]
  base <- if (!nrow(cand)) data.frame() else
    data.frame(arm = arm, collection = cand$collection, ID = cand$ID,
               hit_genes = cand$hit_genes, Count = cand$Count,
               p.adjust = cand$p.adjust, p_matched = cand$p_matched,
               # The raw p and the fold enrichment travel with every candidate because the
               # winner-take-all tie-break needs a quantity comparable across collections;
               # p.adjust is not one. See the ordering in section 6.
               pvalue = cand$pvalue, fold_enrichment = cand$FoldEnrichment,
               frac_matched_reaching_q = cand$frac_matched_reaching_q,
               size_tested = cand$size_tested, enriched = TRUE)
  if (!anchored) return(base)
  arm_genes <- intersect(ARMS[[arm]]$genes, BACKGROUND)
  extra <- dplyr::bind_rows(lapply(names(HYPOXIA), function(cl) {
    ids <- as.character(unlist(HYPOXIA[[cl]] %||% character(0)))
    ids <- setdiff(ids, if (nrow(base)) base$ID else character(0))
    dplyr::bind_rows(lapply(ids, function(id) {
      sz <- SET_SIZES$size_in_background[SET_SIZES$collection == cl &
                                           SET_SIZES$term == id]
      # A pin outside the tested size window cannot carry a hypergeometric share and is
      # reported as an absence of opportunity by section 8, distinct from a zero here.
      if (!length(sz) || sz < MIN_GS || sz > MAX_GS) return(NULL)
      g <- intersect(T2G[[cl]]$gene[T2G[[cl]]$term == id], arm_genes)
      if (!length(g)) return(NULL)
      data.frame(arm = arm, collection = cl, ID = id, hit_genes = collapse_genes(g),
                 Count = length(g), p.adjust = NA_real_, p_matched = NA_real_,
                 pvalue = NA_real_, fold_enrichment = NA_real_,
                 frac_matched_reaching_q = NA_real_, size_tested = as.integer(sz),
                 enriched = FALSE)
    }))
  }))
  dplyr::bind_rows(base, extra)
}

select_sets <- function(arm, anchored) {
  cand <- candidates_for(arm, anchored)
  if (!nrow(cand)) return(list(kept = data.frame(), audit = data.frame()))
  pinned <- character(0)
  if (anchored) {
    for (cl in names(HYPOXIA)) {
      pinned <- c(pinned, as.character(unlist(HYPOXIA[[cl]])))
    }
  }
  audit <- cand[, c("arm", "collection", "ID", "Count", "p.adjust", "p_matched",
                    "frac_matched_reaching_q", "size_tested", "enriched")]
  audit$pinned <- audit$ID %in% pinned
  audit$fails_p_matched <- !audit$pinned & (is.na(audit$p_matched) | audit$p_matched >= PM_CUT)
  audit$fails_recurrence <- !audit$pinned &
    (is.na(audit$frac_matched_reaching_q) | audit$frac_matched_reaching_q > MAX_RECUR)
  live <- audit[!audit$fails_p_matched & !audit$fails_recurrence, , drop = FALSE]
  # Ordered by evidence, NOT by pin status: a pin must not gain priority, or it could
  # prune away a category that earned its place. NA p_matched (a pin that never enriched)
  # sorts last for the same reason.
  live <- live[order(live$p_matched, live$p.adjust, na.last = TRUE), , drop = FALSE]

  # Redundancy prune on ARM-hit overlap: two sets claiming nearly the same arm genes
  # cannot both carry a share. Pinned sets are exempt from being dropped, and are also
  # never used AS a pruner, so anchoring can add a category but never remove one.
  keep_ids <- character(0); pruners <- character(0)
  pruned_by <- character(0); pruned_j <- numeric(0)
  hitmap <- setNames(lapply(live$ID, function(id)
    hits_of(cand[cand$ID == id, ][1, ])), live$ID)
  for (i in seq_len(nrow(live))) {
    id <- live$ID[i]
    drop_because <- NA_character_; drop_j <- NA_real_
    if (!live$pinned[i]) {
      for (k in pruners) {
        j <- jaccard(hitmap[[id]], hitmap[[k]])
        if (j > PRUNE_J) { drop_because <- k; drop_j <- j; break }
      }
    }
    if (is.na(drop_because)) {
      keep_ids <- c(keep_ids, id)
      if (!live$pinned[i]) pruners <- c(pruners, id)
    }
    pruned_by <- c(pruned_by, drop_because); pruned_j <- c(pruned_j, drop_j)
  }
  live$pruned_by <- pruned_by; live$pruned_hit_jaccard <- pruned_j
  live$survived_prune <- is.na(live$pruned_by)

  # Top N per collection, pins exempt from the cap.
  surv <- live[live$survived_prune, , drop = FALSE]
  surv <- surv %>%
    dplyr::group_by(.data$collection) %>%
    dplyr::mutate(rank_in_collection = dplyr::row_number()) %>%
    dplyr::ungroup() %>% as.data.frame()
  surv$selected <- surv$pinned | surv$rank_in_collection <= TOP_N
  live <- dplyr::left_join(
    live, surv[, c("ID", "rank_in_collection", "selected")], by = "ID")
  live$selected[is.na(live$selected)] <- FALSE

  audit <- dplyr::left_join(
    audit,
    live[, c("ID", "pruned_by", "pruned_hit_jaccard", "rank_in_collection", "selected")],
    by = "ID")
  audit$selected[is.na(audit$selected)] <- FALSE
  audit$variant <- variant_label(anchored)
  list(kept = audit[audit$selected, , drop = FALSE], audit = audit)
}

SELECTED <- list(); AUDIT <- list(); CANDS <- list()
for (arm in ARM_ORDER) {
  for (anch in c(FALSE, TRUE)) {
    v <- variant_label(anch)
    key <- paste(arm, v, sep = "|")
    r <- select_sets(arm, anch)
    SELECTED[[key]] <- r$kept
    AUDIT[[key]] <- r$audit
    CANDS[[key]] <- candidates_for(arm, anch)
    n_pin <- if (nrow(r$kept)) sum(r$kept$pinned) else 0L
    n_pin_null <- if (nrow(r$kept)) sum(r$kept$pinned & !r$kept$enriched) else 0L
    message("[5] ", arm, " / ", v, ": ", nrow(r$kept), " sets selected of ",
            nrow(r$audit), " candidates", if (anch) paste0(
              "  (", n_pin, " pinned, ", n_pin_null, " of them never enriched)") else "")
  }
}

## ---------------------------------------------------------------------------
## 6. Attribution, two ways
##
## fractional      a gene in k selected sets gives 1/k to each; sums to the arm size
## winner_take_all each gene to its single best selected set; a true partition
##
## Both sum to the arm's in-background size, so a share is always a share of the whole
## arm and the two are directly comparable.
## ---------------------------------------------------------------------------

ASSIGN <- list()
for (key in names(SELECTED)) {
  parts <- strsplit(key, "|", fixed = TRUE)[[1]]
  arm <- parts[1]; variant <- parts[2]
  kept <- SELECTED[[key]]
  arm_genes <- intersect(ARMS[[arm]]$genes, BACKGROUND)
  cand <- CANDS[[key]]

  memb <- if (!nrow(kept)) data.frame(gene = character(0), ID = character(0),
                                      collection = character(0), p_matched = numeric(0),
                                      p.adjust = numeric(0), pvalue = numeric(0),
                                      fold_enrichment = numeric(0), enriched = logical(0))
  else dplyr::bind_rows(lapply(seq_len(nrow(kept)), function(i) {
    id <- kept$ID[i]
    g <- intersect(hits_of(cand[cand$ID == id, ][1, ]), arm_genes)
    if (!length(g)) return(NULL)
    r <- cand[cand$ID == id, ][1, ]
    data.frame(gene = g, ID = id, collection = kept$collection[i],
               p_matched = kept$p_matched[i], p.adjust = kept$p.adjust[i],
               pvalue = r$pvalue %||% NA_real_,
               fold_enrichment = r$fold_enrichment %||% NA_real_,
               enriched = kept$enriched[i])
  }))

  k_per_gene <- if (nrow(memb)) table(memb$gene) else integer(0)
  rows <- if (nrow(memb)) {
    memb$n_sets_for_gene <- as.integer(k_per_gene[memb$gene])
    memb$weight_fractional <- 1 / memb$n_sets_for_gene
    # Winner-take-all: a pinned set that never enriched must never take a gene from a set
    # that did, so `enriched` leads the ordering before any p-value.
    #
    # The tie-break AFTER p_matched must be the RAW hypergeometric p, upstream of the adjusted
    # one. enricher() is called once per collection, so p.adjust is Benjamini-Hochberg
    # within five separate families of very different size (Hallmark 50 testable sets,
    # KEGG 177, GO MF 928, Reactome 1295, GO BP 4449). BH scales with family size, so the
    # same raw p becomes roughly 89-fold smaller in Hallmark than in GO BP, and ordering on
    # it hands shared genes to the smaller collections for no reason to do with the data.
    #
    # This is not a corner case. p_matched floors at 1/(n_null+1) and 30 of the 40 sets
    # selected for WT_heat_up sit exactly on that floor, so the primary key is tied for
    # most comparisons: 75 of the 83 genes held by more than one selected set had their
    # winner decided by the key that follows it. The raw p is the same test in every
    # collection and is comparable across them; FoldEnrichment is the final fallback.
    best <- memb[order(memb$gene, !memb$enriched, memb$p_matched, memb$pvalue,
                       -memb$fold_enrichment, na.last = TRUE), ]
    best <- best[!duplicated(best$gene), c("gene", "ID")]
    names(best)[2] <- "winner_id"
    memb <- dplyr::left_join(memb, best, by = "gene")
    memb$is_winner <- memb$ID == memb$winner_id
    memb
  } else memb

  unclaimed <- setdiff(arm_genes, if (nrow(rows)) unique(rows$gene) else character(0))
  res <- if (length(unclaimed)) data.frame(
    gene = unclaimed, ID = RESIDUAL, collection = RESIDUAL,
    p_matched = NA_real_, p.adjust = NA_real_, pvalue = NA_real_,
    fold_enrichment = NA_real_, n_sets_for_gene = 0L,
    weight_fractional = 1, winner_id = RESIDUAL, is_winner = TRUE) else NULL
  out <- dplyr::bind_rows(rows, res)
  out$arm <- arm; out$variant <- variant
  ASSIGN[[key]] <- out

  # Invariant: both accountings must sum to the arm's in-background size.
  stopifnot(abs(sum(out$weight_fractional) - length(arm_genes)) < 1e-8,
            sum(out$is_winner) == length(arm_genes))
}
ASSIGNMENT <- dplyr::bind_rows(ASSIGN)

SHARES <- dplyr::bind_rows(lapply(names(ASSIGN), function(key) {
  a <- ASSIGN[[key]]
  n <- length(intersect(ARMS[[a$arm[1]]]$genes, BACKGROUND))
  frac <- a %>% dplyr::group_by(.data$ID, .data$collection) %>%
    dplyr::summarise(n_genes = dplyr::n_distinct(.data$gene),
                     weight = sum(.data$weight_fractional), .groups = "drop") %>%
    dplyr::mutate(accounting = "fractional", share = .data$weight / n)
  wta <- a[a$is_winner, ] %>% dplyr::group_by(.data$ID, .data$collection) %>%
    dplyr::summarise(n_genes = dplyr::n_distinct(.data$gene),
                     weight = dplyr::n(), .groups = "drop") %>%
    dplyr::mutate(accounting = "winner_take_all", share = .data$weight / n)
  out <- dplyr::bind_rows(frac, wta)
  out$arm <- a$arm[1]; out$variant <- a$variant[1]; out$n_arm_in_background <- n
  out
}))

## ---------------------------------------------------------------------------
## 7. The residue, split by reason
## ---------------------------------------------------------------------------

annotated_anywhere <- unique(unlist(lapply(T2G, function(t) t$gene)))
testable_sets <- SET_SIZES[SET_SIZES$size_in_background >= MIN_GS &
                             SET_SIZES$size_in_background <= MAX_GS, ]
annotated_testable <- unique(unlist(lapply(COLLECTIONS, function(cl) {
  t <- T2G[[cl]]
  t$gene[t$term %in% testable_sets$term[testable_sets$collection == cl]]
})))

REMAINDER <- dplyr::bind_rows(lapply(names(ASSIGN), function(key) {
  a <- ASSIGN[[key]]
  arm <- a$arm[1]
  unc <- a$gene[a$ID == RESIDUAL]
  enr_any <- unique(ENRICHED$hit_genes[ENRICHED$arm == arm])
  enr_genes <- unique(unlist(strsplit(enr_any, "/", fixed = TRUE)))
  data.frame(
    arm = arm, variant = a$variant[1],
    n_arm_in_background = length(intersect(ARMS[[arm]]$genes, BACKGROUND)),
    n_unclaimed = length(unc),
    n_unclaimed_in_an_enriched_set = length(intersect(unc, enr_genes)),
    n_unclaimed_annotated_testable = length(intersect(unc, annotated_testable)),
    n_unclaimed_annotated_only_untestable_size =
      length(setdiff(intersect(unc, annotated_anywhere), annotated_testable)),
    n_unclaimed_no_annotation_at_all = length(setdiff(unc, annotated_anywhere)),
    genes_no_annotation_at_all = collapse_genes(setdiff(unc, annotated_anywhere)))
}))

## ---------------------------------------------------------------------------
## 8. The anchored variant's hypoxia sets, and their mutual overlap
##
## "Hypoxia" is not one category across collections. Reporting a single hypoxia bar would
## double count: in WT_heat_up the Hallmark and GO hypoxia hits share very few genes.
## Two of the five collections carry no hypoxia set at all, which is an absence of
## opportunity and not a zero share.
## ---------------------------------------------------------------------------

HYP <- dplyr::bind_rows(lapply(ARM_ORDER, function(arm) {
  arm_genes <- intersect(ARMS[[arm]]$genes, BACKGROUND)
  dplyr::bind_rows(lapply(COLLECTIONS, function(cl) {
    ids <- as.character(unlist(HYPOXIA[[cl]] %||% character(0)))
    if (!length(ids))
      return(data.frame(arm = arm, collection = cl, term = NA_character_,
                        status = "collection carries no hypoxia set",
                        size_in_background = NA_integer_, n_arm_hits = NA_integer_,
                        hits = NA_character_, p.adjust = NA_real_, p_matched = NA_real_))
    dplyr::bind_rows(lapply(ids, function(id) {
      sz <- SET_SIZES$size_in_background[SET_SIZES$collection == cl &
                                           SET_SIZES$term == id]
      g <- intersect(T2G[[cl]]$gene[T2G[[cl]]$term == id], arm_genes)
      e <- ENRICHED[ENRICHED$arm == arm & ENRICHED$ID == id, ]
      status <- if (!length(sz)) "absent from collection"
      else if (sz < MIN_GS) paste0("below min_gs_size (", sz, ")")
      else if (sz > MAX_GS) paste0("above max_gs_size (", sz, ")")
      else if (!nrow(e)) "tested, not enriched" else "tested, enriched"
      data.frame(arm = arm, collection = cl, term = id, status = status,
                 size_in_background = if (length(sz)) as.integer(sz) else NA_integer_,
                 n_arm_hits = length(g), hits = collapse_genes(g),
                 p.adjust = if (nrow(e)) e$p.adjust[1] else NA_real_,
                 p_matched = if (nrow(e)) e$p_matched[1] else NA_real_)
    }))
  }))
}))

HYP_OVERLAP <- dplyr::bind_rows(lapply(ARM_ORDER, function(arm) {
  h <- HYP[HYP$arm == arm & !is.na(HYP$term) & HYP$n_arm_hits > 0, ]
  if (nrow(h) < 2) return(NULL)
  cmb <- utils::combn(seq_len(nrow(h)), 2)
  dplyr::bind_rows(lapply(seq_len(ncol(cmb)), function(k) {
    i <- cmb[1, k]; j <- cmb[2, k]
    a <- strsplit(h$hits[i], "/")[[1]]; b <- strsplit(h$hits[j], "/")[[1]]
    data.frame(arm = arm, term_a = h$term[i], term_b = h$term[j],
               n_a = length(a), n_b = length(b), n_shared = length(intersect(a, b)),
               jaccard = jaccard(a, b))
  }))
}))

## ---------------------------------------------------------------------------
## 9. Concordance with 15_go_decomposition's GO BP
##
## Not the same test: that stage propagates org.Hs.eg.db with an IEA switch, this one uses
## MSigDB C5 GO:BP. Reported descriptively so a reader can see how much the mapping
## choice moves the answer.
## ---------------------------------------------------------------------------

GO15 <- file.path(stage_dir("15_go_decomposition", "tables"), "_overview",
                  "go_enriched_terms.csv")
CONCORD <- if (!file.exists(GO15)) {
  message("[9] 15_go_decomposition table absent; concordance skipped")
  data.frame()
} else {
  g15 <- utils::read.csv(GO15)
  dplyr::bind_rows(lapply(ARM_ORDER, function(arm) {
    a <- g15[g15$arm == arm & g15$ontology == "BP" & g15$iea_variant == "with_iea", ]
    b <- ENRICHED[ENRICHED$arm == arm & ENRICHED$collection == "GO_BP", ]
    # names differ by construction (GO id vs MSigDB stem), so compare on the term label
    lab15 <- toupper(gsub("[^A-Z0-9]+", "_", toupper(a$Description)))
    lab26 <- toupper(sub("^GOBP_", "", b$ID))
    data.frame(arm = arm,
               n_terms_15_go_decomposition = nrow(a), n_terms_here = nrow(b),
               n_label_matched = length(intersect(lab15, lab26)),
               jaccard_by_label = jaccard(lab15, lab26),
               n_genes_covered_15 = length(unique(unlist(strsplit(a$geneID, "/")))),
               n_genes_covered_here = length(unique(unlist(strsplit(b$hit_genes, "/")))))
  }))
}

## ---------------------------------------------------------------------------
## 10. Write
## ---------------------------------------------------------------------------

NULL_SUMMARY <- dplyr::bind_rows(lapply(names(NULLS), function(nm) {
  ne <- NULLS[[nm]]; if (is.null(ne)) return(NULL)
  parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
  s <- ne$summary
  s$arm <- parts[1]; s$collection <- parts[2]
  s$n_query_annotated <- ne$n_query_annotated
  s$n_terms_testable <- ne$n_terms_testable
  s$not_tested_reason <- NA_character_
  s
}))
# Untested cells ride in the same table so a reader counting rows finds every arm x
# collection pair, with the reason in place of the statistics.
if (length(NOT_TESTED))
  NULL_SUMMARY <- dplyr::bind_rows(NULL_SUMMARY, dplyr::bind_rows(NOT_TESTED))

# The permutation floor is a resolution limit and belongs on the record: p_matched cannot
# go below 1/(n_null+1), and where a selected set sits on that floor the ordering that
# follows it is what actually decided the assignment.
PM_FLOOR <- 1 / (N_NULL + 1)
SATURATION <- dplyr::bind_rows(lapply(names(SELECTED), function(key) {
  k <- SELECTED[[key]]
  if (!nrow(k)) return(NULL)
  parts <- strsplit(key, "|", fixed = TRUE)[[1]]
  data.frame(arm = parts[1], variant = parts[2], n_selected = nrow(k),
             p_matched_floor = PM_FLOOR,
             n_at_floor = sum(!is.na(k$p_matched) & k$p_matched <= PM_FLOOR + 1e-12),
             frac_at_floor = mean(!is.na(k$p_matched) & k$p_matched <= PM_FLOOR + 1e-12))
}))

w <- function(df, stem) {
  p <- file.path(TBL_OVW, paste0(stem, ".csv"))
  utils::write.csv(df, p, row.names = FALSE)
  message("    -> ", p, "  (", nrow(df), " rows)")
}

PROV <- data.frame(key = c(
  "script", "stage", "run_date", "R_version", "clusterProfiler_version",
  "genesets_object", "genesets_md5", "collections", "collection_sources",
  "background", "background_n", "p_cutoff", "q_cutoff", "p_adjust_method",
  "min_gs_size", "max_gs_size", "p_matched_cutoff", "max_null_recurrence",
  "prune_hit_jaccard", "top_n_per_collection", "accountings", "n_null",
  "depth_band_quantiles", "seed", "residual_label",
  "p_matched_floor", "winner_take_all_tiebreak",
  "seam_max_abs_dk", "seam_max_abs_dq"),
  value = c(
    SCRIPT, STAGE, as.character(Sys.Date()), R.version.string,
    as.character(utils::packageVersion("clusterProfiler")),
    AC$genesets_object, GS_MD5, paste(COLLECTIONS, collapse = ","),
    paste(sprintf("%s=%s", COLLECTIONS, SET_SOURCE), collapse = " | "),
    GO$background, length(BACKGROUND), P_CUT, Q_CUT, PADJ_METHOD,
    MIN_GS, MAX_GS, PM_CUT, MAX_RECUR, PRUNE_J, TOP_N,
    paste(ACCOUNTINGS, collapse = ","), N_NULL,
    paste(BAND_Q, collapse = ","), SEED, RESIDUAL,
    1 / (N_NULL + 1),
    "enriched, then p_matched, then raw hypergeometric pvalue, then fold_enrichment",
    max(SEAM$max_abs_dk, na.rm = TRUE), max(SEAM$max_abs_dq, na.rm = TRUE)))

# Every frame carrying a `variant` column must use the same vocabulary, or the directory
# disagrees with itself and a filter on the published label silently returns nothing.
for (.nm in c("ENRICHED", "ASSIGNMENT", "SHARES", "REMAINDER", "SATURATION")) {
  .f <- get(.nm)
  if (!is.data.frame(.f) || !"variant" %in% names(.f)) next
  .bad <- setdiff(unique(.f$variant), VARIANT_LEVELS)
  if (length(.bad))
    stop(sprintf("frame %s carries variant level(s) outside the vocabulary: %s",
                 .nm, paste(.bad, collapse = ", ")))
}
.aud <- dplyr::bind_rows(AUDIT)
if ("variant" %in% names(.aud)) {
  .bad <- setdiff(unique(.aud$variant), VARIANT_LEVELS)
  if (length(.bad))
    stop("the audit frame carries variant level(s) outside the vocabulary: ",
         paste(.bad, collapse = ", "))
}
message("[10] variant vocabulary consistent across frames: ",
        paste(VARIANT_LEVELS, collapse = ", "))

message("[10] writing")
w(ENRICHED, "composition_ora_terms")
w(dplyr::bind_rows(AUDIT), "composition_selected")
w(ASSIGNMENT, "composition_gene_assignment")
w(SHARES, "composition_shares")
w(REMAINDER, "composition_remainder")
w(HYP, "composition_hypoxia_sources")
w(HYP_OVERLAP, "composition_hypoxia_overlap")
w(NULL_SUMMARY, "composition_null_summary")
w(SATURATION, "composition_permutation_floor")
w(CONCORD, "composition_go_concordance")
w(PROV, "composition_provenance")

saveRDS(list(enriched = ENRICHED, audit = dplyr::bind_rows(AUDIT),
             assignment = ASSIGNMENT, shares = SHARES, remainder = REMAINDER,
             hypoxia = HYP, hypoxia_overlap = HYP_OVERLAP,
             null_summary = NULL_SUMMARY, saturation = SATURATION,
             concordance = CONCORD,
             provenance = PROV, arms = ARMS, residual = RESIDUAL),
        file.path(DIR_OBJECTS, "26_arm_composition.rds"))

message("26_arm_composition: COMPUTE DONE.")
