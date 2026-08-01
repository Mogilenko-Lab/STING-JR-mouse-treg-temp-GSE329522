#!/usr/bin/env Rscript
# 24_go_decomposition.R -- COMPUTE
# =============================================================================
# Ontology-wide over-representation composition of the projected up arms
# (stage 15_go_decomposition).
#
# Why this stage exists
#   Stage 12 and the human compartment's arm-decomposition both partition an up arm
#   against a hand-picked panel of curated programs. A panel that already contains
#   hypoxia can answer "is hypoxia in this arm?" and nothing else, and for WT_heat_up
#   two thirds of the arm (132 of 199) falls in no panel member at all. This stage asks
#   the ontology instead: over the WHOLE GO Biological Process space, and separately the
#   whole Molecular Function space, which terms does each arm's gene content actually
#   organise into? The candidate space is not curated by whoever is asking.
#
#   Over-representation is new machinery in this repository. Every other enrichment here
#   is pre-ranked GSEA; there is no enrichGO/enricher/groupGO call anywhere else in
#   either compartment. That is a reason for the seam checks below, not a reason to
#   trust the output more.
#
# What the background is, and why
#   03_results/human_projection/signatures/WT_heat/WT_heat_ranked.rnk column 1, 12,986
#   human symbols. This is the space the arms were SELECTED from, carried across the same
#   frozen ortholog map. A generic all-genes universe would credit the arm with
#   enrichment that is really the coverage of the mouse->human map.
#
# The four things this stage emits that a plain enrichGO call does not
#   1. Both IEA variants. enrichGO includes electronically inferred annotation;
#      13_semantic_decomp deliberately drops it. Both are run, both are reported, and
#      the primary is named in config rather than inherited from a default.
#   2. An annotation-depth-matched permutation null, on BOTH headline quantities. The gate
#      that produced these arms thresholds effect size, which preferentially admits
#      well-studied genes; the hypergeometric conditions on set size and universe only, so
#      it does not see this. The null counts how many terms a matched random draw takes to
#      significance AND how many of the drawn genes land in one, because the coverage
#      number is the one a reader takes away and it needs a comparison of its own.
#   3. A GO-term-versus-curated-lens gene-block comparison, so "GO says hypoxia" and
#      "Hallmark says hypoxia" can be read as the two different gene sets they are.
#   4. A searched-and-absent proteostasis probe, so the project's central negative
#      result is a row with a number instead of a silence in the reader's head.
#
# Inputs (frozen; read only)
#   03_results/human_projection/signatures/{WT_heat,KO_heat,Interaction}/*_up.txt
#   03_results/human_projection/signatures/WT_heat/WT_heat_ranked.rnk   (the background)
#   03_results/objects/17_signature_sets.rds + 03_results/objects/gene_universe.txt
#   ../human_treg_arthritis/00_data/references/{msigdb_hallmark,temp_hsr_lens}/*.txt
#   ../sting_positive_control/03_results/06_reference_axis/signatures/*_up.txt
#   org.Hs.eg.db / org.Mm.eg.db / GO.db  (versions recorded in tables/go_provenance.csv)
#
# Outputs -> 03_results/15_go_decomposition/tables/_overview/
#   go_enriched_terms.csv     go_term_summary.csv      go_iea_comparison.csv
#   go_depth_null.csv         go_null_design.csv       go_null_draws.csv
#   go_term_clusters.csv      go_vs_curated_lenses.csv go_gene_coverage.csv
#   go_coverage_summary.csv   go_proteostasis_probe.csv
#   go_mouse_replicate.csv    go_provenance.csv
#   plus 03_results/objects/24_go_decomposition.rds for the viz script.
#   Figures are 25_go_decomposition_viz.R's; this script never plots.
#
# Which GOSemSim build this runs on
#   simplify() and pairwise_termsim() both call the Wang measure. GOSemSim through
#   2.36.0 matches its edge weights against GO.db relationship strings without
#   normalising them: GO.db spells them "isa"/"part of", getSV keys on "is_a"/"part_of",
#   so every BP edge falls through to the "other" weight and the measure runs on a
#   uniform 0.7 -- which erases the one distinction the Wang measure is made of. 2.39.2
#   normalises first. This stage opts in to the project-local library and STOPS rather
#   than fall back to the system build; the version it used is in go_provenance.csv.
#
# Run from project root (~10 min):
#   Rscript 02_analysis/scripts/24_go_decomposition.R
# =============================================================================

source("02_analysis/config/config.R")
use_project_rlib(require_pkgs = "GOSemSim")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(GOSemSim)
  library(enrichplot)
  library(GO.db)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(Matrix)
  library(dplyr)
})
options(stringsAsFactors = FALSE)

STAGE   <- "15_go_decomposition"
SCRIPT  <- "02_analysis/scripts/24_go_decomposition.R"
TBL_DIR <- stage_dir(STAGE, "tables")
TBL_OVW <- file.path(TBL_DIR, "_overview")
dir.create(TBL_OVW, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_SCRATCH, recursive = TRUE, showWarnings = FALSE)

GO <- YAML_CONFIG$go_decomposition
stopifnot(!is.null(GO))

ONTS        <- as.character(unlist(GO$ontologies         %||% c("BP", "MF")))
ONT_PRIMARY <- as.character(GO$primary_ontology          %||% "BP")
IEA_VARS    <- as.character(unlist(GO$iea_variants       %||% c("with_iea", "no_iea")))
IEA_PRIMARY <- as.character(GO$primary_iea               %||% "with_iea")
DROP_EVID   <- as.character(unlist(GO$drop_evidence      %||% "IEA"))
P_CUT       <- as.numeric(GO$p_cutoff                    %||% 0.05)
Q_CUT       <- as.numeric(GO$q_cutoff                    %||% 0.2)
PADJ_METHOD <- as.character(GO$p_adjust_method           %||% "BH")
MIN_GS      <- as.integer(GO$min_gs_size                 %||% 10L)
MAX_GS      <- as.integer(GO$max_gs_size                 %||% 500L)
SIM_MEASURE <- as.character(GO$sim_measure               %||% "Wang")
SIMP_CUT    <- as.numeric(GO$simplify_cutoff             %||% 0.7)
SIMP_BY     <- as.character(GO$simplify_by               %||% "p.adjust")
CLUST_LINK  <- as.character(GO$cluster_linkage           %||% "average")
CLUST_H     <- as.numeric(GO$cluster_height              %||% 0.7)
SEED        <- as.integer(GO$seed                        %||% 20260730L)
N_NULL      <- as.integer(GO$n_null                      %||% 2000L)
NULL_ONTS   <- as.character(unlist(GO$null_ontologies    %||% "BP"))
BAND_Q      <- as.numeric(unlist(GO$depth_band_quantiles %||% c(.2, .4, .55, .7, .8, .875, .94, .98)))
PROBE_IDS   <- as.character(unlist(GO$proteostasis_probe))
PROBE_PAT   <- as.character(GO$probe_name_pattern %||% "chaperon|unfolded|heat")

GSS_VER <- as.character(utils::packageVersion("GOSemSim"))
message("=================================================================")
message("24_go_decomposition: whole-ontology over-representation (", STAGE, ")")
message("  ontologies ", paste(ONTS, collapse = "/"), " (primary ", ONT_PRIMARY, ")",
        "   IEA variants ", paste(IEA_VARS, collapse = "/"), " (primary ", IEA_PRIMARY, ")")
message("  GOSemSim ", GSS_VER, " from ", find.package("GOSemSim"),
        "   nulls ", N_NULL, "   seed ", SEED)
message("=================================================================")

## ---------------------------------------------------------------------------
## 0. Small helpers
## ---------------------------------------------------------------------------

read_symbols <- function(path) unique(trimws(readLines(path, warn = FALSE)))

zero_na <- function(x) { x[is.na(x)] <- 0L; as.integer(x) }

# p.adjust()'s default n counts NAs, which would silently deflate every adjusted
# p-value in a null replicate (most terms have zero hits and are not in the family).
# Adjust the non-NA entries only -- exactly the family enricher() adjusts over.
bh_col <- function(p) {
  ok <- !is.na(p)
  out <- rep(NA_real_, length(p))
  out[ok] <- stats::p.adjust(p[ok], method = PADJ_METHOD)
  out
}

go_attr <- function(ids, what = c("Term", "Ontology")) {
  what <- match.arg(what)
  f <- if (what == "Term") AnnotationDbi::Term else AnnotationDbi::Ontology
  vapply(ids, function(i) {
    t <- tryCatch(GO.db::GOTERM[[i]], error = function(e) NULL)
    if (is.null(t)) NA_character_ else f(t)
  }, character(1), USE.NAMES = FALSE)
}

collapse_genes <- function(g) if (!length(g)) "" else paste(sort(g), collapse = "/")

# One vocabulary for "why is there no test result in this row", shared by the curated-lens
# table and the proteostasis probe. A term that never entered the hypergeometric and a term
# that entered and came back null are different statements, and a bare FALSE says the
# second when it means the first. The branches are in the order the gates apply.
test_status <- function(ontology, term_size, k_arm) {
  if (is.na(ontology))                      return("absent_from_GO.db")
  if (!ontology %in% ONTS)                  return("ontology_not_run")
  if (is.na(term_size) || term_size == 0L)  return("term_absent_from_universe")
  if (term_size < MIN_GS)                   return("below_min_gs_size")
  if (term_size > MAX_GS)                   return("above_max_gs_size")
  if (is.na(k_arm) || k_arm == 0L)          return("no_arm_gene_in_term")
  "tested"
}

## ---------------------------------------------------------------------------
## 1. Inputs: arms, background, curated lenses
## ---------------------------------------------------------------------------

BG_PATH <- file.path(PROJECT_ROOT, GO$background)
stopifnot(file.exists(BG_PATH))
bg_tbl <- utils::read.delim(BG_PATH, header = FALSE, col.names = c("gene", "score"))
BACKGROUND <- unique(bg_tbl$gene)
message("[1] background: ", length(BACKGROUND), " human symbols from ", GO$background)

ARMS <- lapply(GO$arms, function(a) {
  p <- file.path(PROJECT_ROOT, a$file)
  stopifnot(file.exists(p))
  g <- read_symbols(p)
  list(name = as.character(a$name), role = as.character(a$role %||% "comparator"),
       file = as.character(a$file), genes = g,
       n_nominal = length(g), n_in_background = length(intersect(g, BACKGROUND)))
})
names(ARMS) <- vapply(ARMS, function(a) a$name, character(1))
for (a in ARMS)
  message("[1] arm ", a$name, ": ", a$n_nominal, " genes, ", a$n_in_background,
          " in the background")
ARM_ORDER <- names(ARMS)

LENSES <- lapply(GO$curated_lenses, function(l) {
  p <- file.path(PROJECT_ROOT, l$file)
  list(name = as.character(l$name), go_id = as.character(l$go_id),
       file = as.character(l$file), present = file.exists(p),
       genes = if (file.exists(p)) read_symbols(p) else character(0))
})
names(LENSES) <- vapply(LENSES, function(l) l$name, character(1))
for (l in LENSES)
  message("[1] lens ", l$name, ": ", length(l$genes), " genes",
          if (!l$present) "  [FILE MISSING]" else "")
stopifnot(all(vapply(LENSES, function(l) l$present, logical(1))))

## ---------------------------------------------------------------------------
## 2. The annotation map -- ONE propagated gene->term table, both IEA variants
##
## AnnotationDbi's GOALL/ONTOLOGYALL/EVIDENCEALL columns give the ancestor-propagated
## map with the evidence code attached, which is exactly the map clusterProfiler's
## get_GO_data() builds for enrichGO (it just never filters on evidence). Building it
## here once buys three things: the two IEA variants differ ONLY in an evidence filter
## and are otherwise the same code path; the permutation null uses the identical map the
## test used, so a null hit count and an observed hit count are the same quantity; and
## the map is cached, because the select() call is the slow step.
## ---------------------------------------------------------------------------

build_go_map <- function(orgdb_name, symbols, tag) {
  orgdb_v <- as.character(utils::packageVersion(orgdb_name))
  cache <- file.path(DIR_SCRATCH, sprintf("24_goall__%s__%s.rds", tag, orgdb_v))
  key <- list(symbols = sort(symbols), orgdb = orgdb_name, orgdb_v = orgdb_v,
              godb_v = as.character(utils::packageVersion("GO.db")))
  if (file.exists(cache)) {
    cached <- readRDS(cache)
    if (identical(cached$key, key)) {
      message("[2] ", tag, ": reusing cached GOALL map (", nrow(cached$map), " rows)")
      return(cached$map)
    }
  }
  message("[2] ", tag, ": building the propagated GOALL map from ", orgdb_name, " ...")
  db <- get(orgdb_name, envir = asNamespace(orgdb_name))
  m <- suppressMessages(AnnotationDbi::select(
    db, keys = sort(symbols), keytype = "SYMBOL",
    columns = c("GOALL", "ONTOLOGYALL", "EVIDENCEALL")))
  m <- m[!is.na(m$GOALL), c("SYMBOL", "GOALL", "ONTOLOGYALL", "EVIDENCEALL")]
  names(m) <- c("gene", "term", "ontology", "evidence")
  saveRDS(list(key = key, map = m), cache)
  message("[2] ", tag, ": ", nrow(m), " gene-term-evidence rows -> ", cache)
  m
}

GOMAP_HS <- build_go_map("org.Hs.eg.db", BACKGROUND, "human")

#' Term-to-gene table for one ontology and one IEA variant, restricted to a universe.
term2gene <- function(map, ont, variant, universe) {
  m <- map[map$ontology == ont & map$gene %in% universe, ]
  if (identical(variant, "no_iea")) m <- m[!m$evidence %in% DROP_EVID, ]
  unique(data.frame(term = m$term, gene = m$gene))
}

## ---------------------------------------------------------------------------
## 3. Over-representation: every arm x ontology x IEA variant
##
## enricher() rather than enrichGO() so both IEA variants run through one code path.
## The equivalence is not assumed: 3b re-runs the primary arm through enrichGO() and
## records the agreement in go_provenance.csv.
## ---------------------------------------------------------------------------

run_ora <- function(genes, t2g, universe) {
  ids <- unique(t2g$term)
  t2n <- data.frame(term = ids, name = go_attr(ids, "Term"))
  t2n$name[is.na(t2n$name)] <- t2n$term[is.na(t2n$name)]
  suppressMessages(clusterProfiler::enricher(
    gene = genes, TERM2GENE = t2g, TERM2NAME = t2n, universe = universe,
    pAdjustMethod = PADJ_METHOD, pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT,
    minGSSize = MIN_GS, maxGSSize = MAX_GS))
}

T2G <- list(); ORA <- list(); ORA_OBJ <- list(); MAP_STATS <- list()
for (ont in ONTS) for (v in IEA_VARS) {
  key <- paste(ont, v, sep = "|")
  T2G[[key]] <- term2gene(GOMAP_HS, ont, v, BACKGROUND)
  MAP_STATS[[key]] <- data.frame(
    ontology = ont, iea_variant = v,
    n_background = length(BACKGROUND),
    n_background_annotated = length(unique(T2G[[key]]$gene)),
    n_terms_in_map = length(unique(T2G[[key]]$term)),
    n_gene_term_pairs = nrow(T2G[[key]]))
  message("[3] map ", key, ": ", MAP_STATS[[key]]$n_background_annotated,
          " annotated background genes, ", MAP_STATS[[key]]$n_terms_in_map, " terms")
  for (arm in ARM_ORDER) {
    e <- run_ora(ARMS[[arm]]$genes, T2G[[key]], BACKGROUND)
    ORA_OBJ[[paste(arm, key, sep = "|")]] <- e
    df <- if (is.null(e)) data.frame() else as.data.frame(e)
    if (nrow(df)) {
      df$arm <- arm; df$ontology <- ont; df$iea_variant <- v
      df$n_query_annotated    <- as.integer(sub(".*/", "", df$GeneRatio))
      df$n_universe_annotated <- as.integer(sub(".*/", "", df$BgRatio))
      df$term_size            <- as.integer(sub("/.*", "", df$BgRatio))
      df$fold_enrichment <- (df$Count / df$n_query_annotated) /
        (df$term_size / df$n_universe_annotated)
    }
    ORA[[paste(arm, key, sep = "|")]] <- df
    message("      arm ", arm, " (", ont, ", ", v, "): ", nrow(df), " enriched terms")
  }
}

## 3b. Seam check -- does the hand-built map reproduce the canonical enrichGO call?
primary_arm <- ARM_ORDER[vapply(ARMS, function(a) identical(a$role, "primary"),
                                logical(1))][1]
eg <- suppressMessages(clusterProfiler::enrichGO(
  gene = ARMS[[primary_arm]]$genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
  ont = ONT_PRIMARY, universe = BACKGROUND, pAdjustMethod = PADJ_METHOD,
  pvalueCutoff = P_CUT, qvalueCutoff = Q_CUT, minGSSize = MIN_GS, maxGSSize = MAX_GS))
eg_df  <- as.data.frame(eg)
own_df <- ORA[[paste(primary_arm, ONT_PRIMARY, IEA_PRIMARY, sep = "|")]]
seam_same_ids <- setequal(eg_df$ID, own_df$ID)
seam_max_dp <- if (seam_same_ids && nrow(eg_df))
  max(abs(eg_df[order(eg_df$ID), "pvalue"] - own_df[order(own_df$ID), "pvalue"])) else NA_real_
message("[3b] enricher-vs-enrichGO seam on ", primary_arm, ": same term ids ",
        seam_same_ids, ", max |dp| ", format(seam_max_dp))

## ---------------------------------------------------------------------------
## 4. Redundancy collapse, two levels, both under the project GOSemSim
##      level 1  simplify()          -- drop a term whose Wang similarity to a MORE
##                                      significant term exceeds the cutoff
##      level 2  pairwise_termsim()  -- + average-linkage cut, giving every enriched
##                                      term a block id and a block representative
## ---------------------------------------------------------------------------

SEMDATA <- list()
for (ont in ONTS)
  SEMDATA[[ont]] <- GOSemSim::godata(annoDb = "org.Hs.eg.db", keytype = "SYMBOL",
                                     ont = ont, computeIC = FALSE)

simplify_keep <- list(); CLUSTERS <- list()
for (nm in names(ORA_OBJ)) {
  parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
  arm <- parts[1]; ont <- parts[2]; v <- parts[3]
  e  <- ORA_OBJ[[nm]]
  df <- ORA[[nm]]
  if (is.null(e) || nrow(df) < 2) {
    simplify_keep[[nm]] <- if (nrow(df)) df$ID else character(0)
    CLUSTERS[[nm]] <- data.frame()
    next
  }
  e@ontology <- ont          # enricher() leaves this UNKNOWN; simplify() dispatches on it
  s <- tryCatch(clusterProfiler::simplify(e, cutoff = SIMP_CUT, by = SIMP_BY,
                                          select_fun = min, measure = SIM_MEASURE,
                                          semData = SEMDATA[[ont]]),
                error = function(err) {
                  message("      simplify failed for ", nm, ": ", conditionMessage(err))
                  NULL })
  simplify_keep[[nm]] <- if (is.null(s)) character(0) else as.data.frame(s)$ID
  pt <- tryCatch(enrichplot::pairwise_termsim(e, method = SIM_MEASURE,
                                              semData = SEMDATA[[ont]],
                                              showCategory = nrow(df)),
                 error = function(err) NULL)
  if (!is.null(pt) && !is.null(pt@termsim) && nrow(pt@termsim) > 1) {
    S <- pt@termsim
    S[is.na(S)] <- 0
    S <- pmax(S, t(S)); diag(S) <- 1
    ## pairwise_termsim() labels its matrix by Description, since that is what the
    ## enrichplot figures print. The clustering is unaffected, but every join
    ## downstream keys on ID, so translate the labels here. Left as descriptions the
    ## join below matches nothing and NAs the entire block structure without raising,
    ## which is the same silent-seam failure the gene map is guarded against.
    if (!all(rownames(S) %in% df$ID)) {
      desc2id <- stats::setNames(df$ID, df$Description)
      relabelled <- unname(desc2id[rownames(S)])
      if (anyNA(relabelled))
        stop("termsim labels for ", nm, " are neither GO ids nor known Descriptions (",
             sum(is.na(relabelled)), " of ", length(relabelled), " unmatched)")
      dimnames(S) <- list(relabelled, relabelled)
    }
    hc <- stats::hclust(stats::as.dist(1 - S), method = CLUST_LINK)
    cl <- stats::cutree(hc, h = CLUST_H)
    cdf <- data.frame(arm = arm, ontology = ont, iea_variant = v,
                      ID = names(cl), cluster_id = as.integer(cl))
    cdf <- dplyr::left_join(cdf, df[, c("ID", "Description", "p.adjust", "Count")],
                            by = "ID")
    if (anyNA(cdf$Description))
      stop("the cluster join dropped ", sum(is.na(cdf$Description)), " of ", nrow(cdf),
           " terms for ", nm, "; termsim labels do not match the enriched ids")
    cdf <- cdf %>%
      dplyr::group_by(.data$cluster_id) %>%
      dplyr::mutate(cluster_size = dplyr::n(),
                    cluster_representative = .data$Description[which.min(.data$p.adjust)][1],
                    cluster_min_padj = min(.data$p.adjust)) %>%
      dplyr::ungroup() %>% as.data.frame()
    CLUSTERS[[nm]] <- cdf
  } else {
    CLUSTERS[[nm]] <- data.frame()
  }
  message("[4] ", nm, ": ", nrow(df), " enriched -> ", length(simplify_keep[[nm]]),
          " after simplify(", SIMP_CUT, "), ",
          length(unique(CLUSTERS[[nm]]$cluster_id)), " blocks at h=", CLUST_H)
}

for (nm in names(ORA)) {
  if (!nrow(ORA[[nm]])) next
  ORA[[nm]]$simplify_kept <- ORA[[nm]]$ID %in% simplify_keep[[nm]]
  cdf <- CLUSTERS[[nm]]
  if (nrow(cdf)) {
    ORA[[nm]] <- dplyr::left_join(
      ORA[[nm]],
      cdf[, c("ID", "cluster_id", "cluster_size", "cluster_representative")], by = "ID")
    if (anyNA(ORA[[nm]]$cluster_id))
      stop("cluster_id is NA for ", sum(is.na(ORA[[nm]]$cluster_id)), " of ",
           nrow(ORA[[nm]]), " enriched terms in ", nm)
  } else {
    ORA[[nm]]$cluster_id <- NA_integer_
    ORA[[nm]]$cluster_size <- NA_integer_
    ORA[[nm]]$cluster_representative <- NA_character_
  }
}

## ---------------------------------------------------------------------------
## 5. Annotation-depth-matched permutation null
##
## Depth is the number of terms a gene carries in the SAME propagated map the test uses,
## i.e. the gene's actual number of chances to be a hit. Bands are background quantiles;
## every replicate reproduces the query's per-band histogram exactly. Where a band pool
## is thinner than the query needs the shortfall is borrowed from the nearest bands and
## COUNTED, so every replicate is exactly the query's size -- a hit count and a share
## both scale with n, so a short replicate is not harmless the way a mean-of-pairwise
## statistic would be.
##
## enricher() is never called in the loop; the whole null is one sparse matrix product.
## ---------------------------------------------------------------------------

null_engine <- function(t2g, query, n_null, seed) {
  genes <- sort(unique(t2g$gene))
  terms <- sort(unique(t2g$term))
  M0 <- Matrix::sparseMatrix(i = match(t2g$term, terms), j = match(t2g$gene, genes),
                             x = 1L, dims = c(length(terms), length(genes)))
  depth <- stats::setNames(Matrix::colSums(M0), genes)
  sz0  <- Matrix::rowSums(M0)
  keep <- sz0 >= MIN_GS & sz0 <= MAX_GS
  M <- M0[keep, , drop = FALSE]; sz <- sz0[keep]; terms <- terms[keep]
  NG <- length(genes)
  qa <- intersect(query, genes)
  n  <- length(qa)
  if (n < 2L || !length(terms)) return(NULL)

  brk <- unique(c(-1, stats::quantile(depth, BAND_Q, names = FALSE), Inf))
  qband <- cut(depth[qa], brk)
  pool  <- setdiff(genes, query)
  pband <- cut(depth[pool], brk)
  need  <- table(qband)
  lev   <- levels(qband)
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
  uniform <- lapply(seq_len(n_null), function(i) list(genes = sample(pool, n),
                                                      borrowed = 0L))

  score <- function(sets) {
    idx <- lapply(sets, function(s) match(s$genes, genes))
    B <- Matrix::sparseMatrix(i = unlist(idx), j = rep(seq_along(idx), lengths(idx)),
                              x = 1L, dims = c(NG, length(idx)))
    K <- as.matrix(M %*% B)
    P <- stats::phyper(K - 1L, sz, NG - sz, n, lower.tail = FALSE)
    P[K == 0L] <- NA_real_        # enricher's family is the terms with >=1 query hit
    Q <- apply(P, 2, bh_col)
    ## THE COVERAGE NULL: per replicate, how many of the DRAWN genes fall in at least one
    ## term reaching significance. Without it the stage's headline coverage number has
    ## nothing to be compared against and cannot be read as evidence. Same statistic as
    ## 26_arm_composition.R's; the two engines must change together.
    covered <- vapply(seq_along(idx), function(j) {
      sig <- which(!is.na(Q[, j]) & Q[, j] < P_CUT)
      if (!length(sig)) return(0L)
      as.integer(sum(Matrix::colSums(M[sig, idx[[j]], drop = FALSE]) > 0))
    }, integer(1))
    list(K = K, Q = Q,
         n_signif = colSums(Q < P_CUT, na.rm = TRUE),
         family   = colSums(K > 0L),
         covered  = covered,
         borrowed = vapply(sets, function(s) s$borrowed, integer(1)))
  }
  rm_ <- score(matched); ru <- score(uniform)

  # Observed, recomputed on the same matrix so a null count and an observed count are
  # the same quantity. Cross-checked against enricher() by the caller.
  kobs <- as.integer(M %*% Matrix::sparseMatrix(i = match(qa, genes), j = rep(1L, n),
                                                x = 1L, dims = c(NG, 1L)))
  pobs <- stats::phyper(kobs - 1L, sz, NG - sz, n, lower.tail = FALSE)
  pobs[kobs == 0L] <- NA_real_
  qobs <- bh_col(pobs)
  sig_obs <- which(!is.na(qobs) & qobs < P_CUT)
  covered_obs <- if (!length(sig_obs)) 0L else as.integer(
    sum(Matrix::colSums(M[sig_obs, match(qa, genes), drop = FALSE]) > 0))
  n_signif_obs <- length(sig_obs)

  per_term <- data.frame(
    ID = terms, term_size = as.integer(sz), k_obs = kobs,
    p_obs_recomputed = pobs, q_obs_recomputed = qobs,
    k_null_matched_mean = rowMeans(rm_$K),
    k_null_matched_q95  = apply(rm_$K, 1, stats::quantile, probs = 0.95, names = FALSE),
    p_matched = (1 + rowSums(rm_$K >= kobs)) / (n_null + 1),
    frac_matched_reaching_q = rowMeans(rm_$Q < P_CUT, na.rm = TRUE),
    k_null_uniform_mean = rowMeans(ru$K),
    p_uniform = (1 + rowSums(ru$K >= kobs)) / (n_null + 1),
    frac_uniform_reaching_q = rowMeans(ru$Q < P_CUT, na.rm = TRUE))
  per_term$frac_matched_reaching_q[is.nan(per_term$frac_matched_reaching_q)] <- 0
  per_term$frac_uniform_reaching_q[is.nan(per_term$frac_uniform_reaching_q)] <- 0

  bands <- data.frame(
    band = lev, band_lower = utils::head(brk, -1), band_upper = utils::tail(brk, -1),
    n_query = as.integer(need),
    n_pool = as.integer(vapply(lev, function(b)
      length(pool_by_band[[b]] %||% character(0)), integer(1))))
  bands$pool_short_of_query <- bands$n_pool < bands$n_query

  ## Every replicate's two statistics, kept whole. Storing only mean/median/q95/max left
  ## the observed term count sitting between q95 and max with no p that could be quoted,
  ## and left the coverage figure with no substrate to draw the null as a distribution.
  draws <- rbind(
    data.frame(null = "depth_matched", replicate = seq_len(n_null),
               n_signif = as.integer(rm_$n_signif), covered = as.integer(rm_$covered),
               genes_borrowed = as.integer(rm_$borrowed)),
    data.frame(null = "uniform", replicate = seq_len(n_null),
               n_signif = as.integer(ru$n_signif), covered = as.integer(ru$covered),
               genes_borrowed = 0L))

  perm_p <- function(v, obs) (1 + sum(v >= obs)) / (n_null + 1)

  list(per_term = per_term, bands = bands, draws = draws,
       n_query_annotated = n, n_universe = NG, n_terms_testable = length(terms),
       n_observed_family = sum(kobs > 0L),
       n_observed_signif = n_signif_obs,
       n_covered_obs = covered_obs,
       summary = data.frame(
         null = c("depth_matched", "uniform"),
         n_replicates = n_null,
         mean_family_size = c(mean(rm_$family), mean(ru$family)),
         mean_n_signif    = c(mean(rm_$n_signif), mean(ru$n_signif)),
         median_n_signif  = c(stats::median(rm_$n_signif), stats::median(ru$n_signif)),
         q95_n_signif     = c(stats::quantile(rm_$n_signif, .95, names = FALSE),
                              stats::quantile(ru$n_signif, .95, names = FALSE)),
         max_n_signif     = c(max(rm_$n_signif), max(ru$n_signif)),
         n_observed_signif = n_signif_obs,
         p_n_signif       = c(perm_p(rm_$n_signif, n_signif_obs),
                              perm_p(ru$n_signif, n_signif_obs)),
         mean_covered     = c(mean(rm_$covered), mean(ru$covered)),
         median_covered   = c(stats::median(rm_$covered), stats::median(ru$covered)),
         q95_covered      = c(stats::quantile(rm_$covered, .95, names = FALSE),
                              stats::quantile(ru$covered, .95, names = FALSE)),
         max_covered      = c(max(rm_$covered), max(ru$covered)),
         n_covered_obs    = covered_obs,
         p_covered        = c(perm_p(rm_$covered, covered_obs),
                              perm_p(ru$covered, covered_obs)),
         mean_genes_borrowed = c(mean(rm_$borrowed), 0),
         max_genes_borrowed  = c(max(rm_$borrowed), 0)))
}

NULLS <- list()
for (ont in NULL_ONTS) {
  key <- paste(ont, IEA_PRIMARY, sep = "|")
  for (arm in ARM_ORDER) {
    message("[5] depth-matched null: ", arm, " (", ont, ", ", IEA_PRIMARY, ") x ", N_NULL)
    NULLS[[paste(arm, key, sep = "|")]] <- null_engine(T2G[[key]], ARMS[[arm]]$genes,
                                                       N_NULL, SEED)
  }
}

# Observed-count seam: the matrix recomputation must agree with enricher() term by term.
NULL_SEAM <- dplyr::bind_rows(lapply(names(NULLS), function(nm) {
  ne <- NULLS[[nm]]; if (is.null(ne)) return(NULL)
  obs <- ORA[[nm]]
  if (!nrow(obs))
    return(data.frame(key = nm, n_compared = 0L, max_abs_dk = NA_real_,
                      max_abs_dq = NA_real_))
  j <- dplyr::inner_join(obs[, c("ID", "Count", "p.adjust")], ne$per_term, by = "ID")
  data.frame(key = nm, n_compared = nrow(j),
             max_abs_dk = max(abs(j$Count - j$k_obs)),
             max_abs_dq = max(abs(j$p.adjust - j$q_obs_recomputed)))
}))
print(NULL_SEAM)

# Attach the null columns to the enriched-term table (primary ontology / IEA variant).
for (nm in names(NULLS)) {
  ne <- NULLS[[nm]]; if (is.null(ne) || !nrow(ORA[[nm]])) next
  ORA[[nm]] <- dplyr::left_join(
    ORA[[nm]],
    ne$per_term[, c("ID", "k_null_matched_mean", "k_null_matched_q95", "p_matched",
                    "frac_matched_reaching_q", "k_null_uniform_mean", "p_uniform",
                    "frac_uniform_reaching_q")], by = "ID")
}

ENRICHED <- dplyr::bind_rows(ORA)
if (nrow(ENRICHED)) {
  ENRICHED <- ENRICHED[order(ENRICHED$arm, ENRICHED$ontology, ENRICHED$iea_variant,
                             ENRICHED$p.adjust), ]
  ENRICHED$arm <- factor(ENRICHED$arm, levels = ARM_ORDER)
  ENRICHED <- ENRICHED[order(ENRICHED$arm), ]
  ENRICHED$arm <- as.character(ENRICHED$arm)
}

## ---------------------------------------------------------------------------
## 6. Per-gene coverage -- the three-way split
##      (a) the ontology cannot see the gene at all
##      (b) it can, but no enriched term contains it
##      (c) at least one enriched term contains it
## ---------------------------------------------------------------------------

GENE_COVERAGE <- dplyr::bind_rows(lapply(ARM_ORDER, function(arm) {
  g <- ARMS[[arm]]$genes
  row <- data.frame(arm = arm, gene = g, in_background = g %in% BACKGROUND)
  for (ont in ONTS) {
    t2g <- T2G[[paste(ont, IEA_PRIMARY, sep = "|")]]
    nt  <- table(t2g$gene)
    row[[paste0("n_terms_", tolower(ont))]] <- zero_na(nt[match(g, names(nt))])
    enr <- ORA[[paste(arm, ont, IEA_PRIMARY, sep = "|")]]
    hits <- if (nrow(enr)) t2g[t2g$term %in% enr$ID, ] else t2g[0, , drop = FALSE]
    ne <- table(hits$gene)
    row[[paste0("n_enriched_terms_", tolower(ont))]] <- zero_na(ne[match(g, names(ne))])
  }
  po <- tolower(ONT_PRIMARY)
  row$coverage_class <- ifelse(
    row[[paste0("n_terms_", po)]] == 0L, "no_annotation",
    ifelse(row[[paste0("n_enriched_terms_", po)]] == 0L, "annotated_no_enriched_term",
           "in_enriched_term"))
  # A gene the primary ontology cannot see may still be visible in another ontology.
  other <- setdiff(tolower(ONTS), po)
  if (length(other))
    row$visible_in_other_ontology <- rowSums(
      as.data.frame(row[, paste0("n_terms_", other), drop = FALSE])) > 0
  else row$visible_in_other_ontology <- NA
  row
}))

COVERAGE_SUMMARY <- GENE_COVERAGE %>%
  dplyr::group_by(.data$arm) %>%
  dplyr::summarise(
    n_genes = dplyr::n(),
    n_in_background = sum(.data$in_background),
    n_no_annotation = sum(.data$coverage_class == "no_annotation"),
    n_annotated_no_enriched_term = sum(.data$coverage_class == "annotated_no_enriched_term"),
    n_in_enriched_term = sum(.data$coverage_class == "in_enriched_term"),
    n_no_bp_but_other_ontology = sum(.data$coverage_class == "no_annotation" &
                                       .data$visible_in_other_ontology %in% TRUE),
    .groups = "drop") %>%
  dplyr::mutate(
    primary_ontology = ONT_PRIMARY, iea_variant = IEA_PRIMARY,
    frac_no_annotation = .data$n_no_annotation / .data$n_genes,
    frac_annotated_no_enriched_term = .data$n_annotated_no_enriched_term / .data$n_genes,
    frac_in_enriched_term = .data$n_in_enriched_term / .data$n_genes) %>%
  as.data.frame()
COVERAGE_SUMMARY <- COVERAGE_SUMMARY[match(ARM_ORDER, COVERAGE_SUMMARY$arm), ]

## The coverage null, carried beside the coverage count it has to be read against.
## n_in_enriched_term is out of the arm's nominal size; the null draws gene sets of the
## arm's ANNOTATED size, so covered_obs is out of that instead. Both are reported, and
## the p is the one-sided permutation p over the depth-matched draws.
cov_null <- dplyr::bind_rows(lapply(ARM_ORDER, function(arm) {
  ne <- NULLS[[paste(arm, ONT_PRIMARY, IEA_PRIMARY, sep = "|")]]
  if (is.null(ne)) return(data.frame(arm = arm))
  s <- ne$summary[ne$summary$null == "depth_matched", ]
  data.frame(arm = arm,
             n_arm_annotated = ne$n_query_annotated,
             n_covered_obs = ne$n_covered_obs,
             null_median_covered = s$median_covered,
             null_q95_covered = s$q95_covered,
             null_max_covered = s$max_covered,
             p_covered_matched = s$p_covered,
             uniform_median_covered =
               ne$summary$median_covered[ne$summary$null == "uniform"])
}))
COVERAGE_SUMMARY <- dplyr::left_join(COVERAGE_SUMMARY, cov_null, by = "arm")

## ---------------------------------------------------------------------------
## 7. GO term versus curated lens, gene block against gene block
##
## The headline: GO "response to hypoxia" and HALLMARK_HYPOXIA are two different gene
## sets, and how much of an arm each one claims is a property of the set, not of the
## word. Reported as three counts and three gene lists for every arm x lens pair.
##
## Every row also carries WHY it has the test result it has. A GO term can miss the
## hypergeometric for four reasons that have nothing to do with the arm -- it is outside
## the size window, it holds no arm gene, it is in an ontology this run does not cover, or
## the installed GO.db has never heard of it -- and all four used to print as
## go_term_enriched FALSE. GO:0006954 inflammatory response is the case that matters: 688
## genes in the with_iea universe puts it above max_gs_size, so it was never tested, while
## claiming 43 of WT_heat_up's 199 genes, the largest claim of any term in the panel. Under
## no_iea it falls under the cap, enters the test, and enriches, so the secondary variant
## is carried on the same row.
## ---------------------------------------------------------------------------

t2g_primary <- T2G[[paste(ONT_PRIMARY, IEA_PRIMARY, sep = "|")]]
t2g_by_ont <- stats::setNames(
  lapply(ONTS, function(o) T2G[[paste(o, IEA_PRIMARY, sep = "|")]]), ONTS)

IEA_SECONDARY <- setdiff(IEA_VARS, IEA_PRIMARY)[1]
t2g_by_ont_alt <- if (is.na(IEA_SECONDARY)) NULL else stats::setNames(
  lapply(ONTS, function(o) T2G[[paste(o, IEA_SECONDARY, sep = "|")]]), ONTS)

#' Everything one GO id has to say about one arm under one IEA variant: the term's size in
#' that universe, how many of the arm's genes it claims, why it did or did not enter the
#' hypergeometric, and -- only when it did -- the result. GO:0006954 inflammatory response
#' is the reason this exists: it holds 688 genes in the with_iea universe, above
#' max_gs_size, so it was never tested, and the earlier `!is.na(p.adjust)` encoding
#' reported it as a plain FALSE beside the largest curated claim in the arm.
go_term_claim <- function(arm, go_id, variant) {
  o <- go_attr(go_id, "Ontology")
  maps <- if (identical(variant, IEA_PRIMARY)) t2g_by_ont else t2g_by_ont_alt
  present <- !is.na(o) && !is.null(maps) && o %in% names(maps)
  gg <- if (!present) character(0) else {
    t <- maps[[o]]; t$gene[t$term == go_id]
  }
  arm_hits <- intersect(ARMS[[arm]]$genes, gg)
  K <- if (present) length(gg) else NA_integer_
  status <- test_status(o, K, length(arm_hits))
  row <- if (present) ORA[[paste(arm, o, variant, sep = "|")]] else data.frame()
  hit <- nrow(row) > 0L && go_id %in% row$ID
  list(ontology = o, in_go_db = !is.na(o), genes = arm_hits,
       n_term_genes = K, n_claimed = length(arm_hits),
       status = status, tested = identical(status, "tested"),
       enriched = if (identical(status, "tested")) hit else NA,
       p_adjust = if (hit) row$p.adjust[match(go_id, row$ID)] else NA_real_,
       qvalue   = if (hit) row$qvalue[match(go_id, row$ID)] else NA_real_,
       rank_in_arm = if (hit) match(go_id, row$ID[order(row$p.adjust)]) else NA_integer_,
       n_enriched_in_arm = nrow(row))
}

GO_VS_CURATED <- dplyr::bind_rows(lapply(ARM_ORDER, function(arm) {
  a <- ARMS[[arm]]$genes
  dplyr::bind_rows(lapply(names(LENSES), function(ln) {
    l <- LENSES[[ln]]
    pr  <- go_term_claim(arm, l$go_id, IEA_PRIMARY)
    go_arm  <- pr$genes
    len_arm <- intersect(a, l$genes)
    shared  <- intersect(go_arm, len_arm)
    uni     <- union(go_arm, len_arm)
    out <- data.frame(
      arm = arm, n_arm_genes = length(a),
      curated_lens = ln, n_lens_genes = length(l$genes),
      go_id = l$go_id, go_term = go_attr(l$go_id, "Term"),
      go_ontology = pr$ontology,
      go_id_in_go_db = pr$in_go_db,
      iea_variant = IEA_PRIMARY,
      n_go_term_genes_in_universe = pr$n_term_genes,
      n_claimed_by_go = pr$n_claimed,
      n_claimed_by_lens = length(len_arm),
      n_shared = length(shared),
      jaccard = if (length(uni)) length(shared) / length(uni) else NA_real_,
      genes_go_only = collapse_genes(setdiff(go_arm, len_arm)),
      genes_shared = collapse_genes(shared),
      genes_lens_only = collapse_genes(setdiff(len_arm, go_arm)),
      go_term_status = pr$status,
      go_term_tested = pr$tested,
      go_term_enriched = pr$enriched,
      go_term_p_adjust = pr$p_adjust,
      go_term_qvalue = pr$qvalue,
      go_term_rank_in_arm = pr$rank_in_arm,
      n_enriched_terms_in_arm = pr$n_enriched_in_arm)
    ## The other IEA variant carried alongside, because the evidence filter decides
    ## whether some of these terms are testable at all: under no_iea the inflammatory
    ## response term falls under max_gs_size, enters the test, and enriches.
    if (!is.na(IEA_SECONDARY)) {
      se <- go_term_claim(arm, l$go_id, IEA_SECONDARY)
      alt <- data.frame(
        se$n_term_genes, se$n_claimed, se$status, se$tested, se$enriched,
        se$p_adjust, se$rank_in_arm, se$n_enriched_in_arm)
      names(alt) <- paste0(IEA_SECONDARY, "_",
                           c("n_go_term_genes_in_universe", "n_claimed_by_go", "status",
                             "tested", "enriched", "p_adjust", "rank_in_arm",
                             "n_enriched_terms_in_arm"))
      out <- cbind(out, alt)
    }
    out
  }))
}))

## The depth-matched null p for the terms that were tested, so the claim size, the
## hypergeometric result and the null verdict all sit on one row.
if (!"p_matched" %in% names(ENRICHED)) ENRICHED$p_matched <- NA_real_
enr_primary <- ENRICHED[ENRICHED$ontology == ONT_PRIMARY &
                          ENRICHED$iea_variant == IEA_PRIMARY,
                        c("arm", "ID", "p_matched", "frac_matched_reaching_q")]
names(enr_primary) <- c("arm", "go_id", "go_term_p_matched",
                        "go_term_frac_matched_reaching_q")
GO_VS_CURATED <- dplyr::left_join(GO_VS_CURATED, enr_primary, by = c("arm", "go_id"))

## ---------------------------------------------------------------------------
## 8. Searched-and-absent: the proteostasis terms a heat-shock reading would predict
##
## Every probe id is looked up and reported with its actual gene count and q whether or
## not it passes anything, so the project's central negative result is a row with a
## number. The name sweep additionally catches a proteostasis-adjacent term nobody
## listed -- which is how "heat generation" surfaces.
## ---------------------------------------------------------------------------

probe_one <- function(arm, go_id) {
  o <- go_attr(go_id, "Ontology")
  base <- data.frame(arm = arm, go_id = go_id, go_term = go_attr(go_id, "Term"),
                     ontology = o, probe_source = "configured_id")
  if (is.na(o) || !o %in% ONTS) {
    return(cbind(base, data.frame(
      status = test_status(o, NA_integer_, NA_integer_),
      term_size_in_universe = NA_integer_, k_arm = NA_integer_, expected_k = NA_real_,
      pvalue = NA_real_, p_adjust = NA_real_, tested = FALSE, enriched = NA,
      genes_hit = "")))
  }
  ne <- NULLS[[paste(arm, o, IEA_PRIMARY, sep = "|")]]
  t2g <- t2g_by_ont[[o]]
  tg  <- t2g$gene[t2g$term == go_id]
  qa  <- intersect(ARMS[[arm]]$genes, unique(t2g$gene))
  n   <- length(qa)
  NGa <- length(unique(t2g$gene))
  K   <- length(tg)
  k   <- length(intersect(qa, tg))
  pv  <- if (K == 0L || n == 0L) NA_real_
         else stats::phyper(k - 1L, K, NGa - K, n, lower.tail = FALSE)
  row <- ORA[[paste(arm, o, IEA_PRIMARY, sep = "|")]]
  padj <- if (nrow(row) && go_id %in% row$ID) row$p.adjust[match(go_id, row$ID)] else NA_real_
  if (is.na(padj) && !is.null(ne)) {
    pt <- ne$per_term
    if (go_id %in% pt$ID) padj <- pt$q_obs_recomputed[match(go_id, pt$ID)]
  }
  # The depth-matched verdict carried on the probe row too: two of these terms clear the
  # hypergeometric, and how often a random depth-matched draw clears it is the number that
  # says whether that means anything.
  pm  <- NA_real_; frq <- NA_real_
  if (!is.null(ne) && go_id %in% ne$per_term$ID) {
    i <- match(go_id, ne$per_term$ID)
    pm  <- ne$per_term$p_matched[i]
    frq <- ne$per_term$frac_matched_reaching_q[i]
  }
  status <- test_status(o, K, k)
  tested <- identical(status, "tested")
  cbind(base, data.frame(
    status = status,
    term_size_in_universe = K, k_arm = k,
    expected_k = if (NGa > 0L) n * K / NGa else NA_real_,
    pvalue = pv, p_adjust = padj, tested = tested,
    # A term the hypergeometric never saw carries NA here, never FALSE.
    enriched = if (tested) nrow(row) > 0L && go_id %in% row$ID else NA,
    p_matched = pm, frac_matched_reaching_q = frq,
    genes_hit = collapse_genes(intersect(qa, tg))))
}

PROTEOSTASIS <- dplyr::bind_rows(lapply(ARM_ORDER, function(arm)
  dplyr::bind_rows(lapply(PROBE_IDS, function(id) probe_one(arm, id)))))

## The name sweep over the ENRICHED terms, so an absence cannot be an oversight.
name_sweep <- if (nrow(ENRICHED)) {
  hit <- grepl(PROBE_PAT, ENRICHED$Description, ignore.case = TRUE) &
    ENRICHED$iea_variant == IEA_PRIMARY
  if (any(hit)) {
    s <- ENRICHED[hit, ]
    data.frame(arm = s$arm, go_id = s$ID, go_term = s$Description,
               ontology = s$ontology, probe_source = "name_sweep",
               status = "tested", term_size_in_universe = s$term_size,
               k_arm = s$Count, expected_k = s$n_query_annotated * s$term_size /
                 s$n_universe_annotated,
               pvalue = s$pvalue, p_adjust = s$p.adjust, tested = TRUE, enriched = TRUE,
               p_matched = if ("p_matched" %in% names(s)) s$p_matched else NA_real_,
               frac_matched_reaching_q = if ("frac_matched_reaching_q" %in% names(s))
                 s$frac_matched_reaching_q else NA_real_,
               genes_hit = s$geneID)
  } else data.frame()
} else data.frame()
PROTEOSTASIS <- dplyr::bind_rows(PROTEOSTASIS, name_sweep)
PROTEOSTASIS$probe_pattern <- PROBE_PAT

## ---------------------------------------------------------------------------
## 9. IEA both ways -- what the evidence filter changes
## ---------------------------------------------------------------------------

IEA_COMPARISON <- dplyr::bind_rows(lapply(ARM_ORDER, function(arm)
  dplyr::bind_rows(lapply(ONTS, function(ont) {
    a <- ORA[[paste(arm, ont, "with_iea", sep = "|")]]
    b <- ORA[[paste(arm, ont, "no_iea", sep = "|")]]
    ia <- if (nrow(a)) a$ID else character(0)
    ib <- if (nrow(b)) b$ID else character(0)
    top_n <- min(10L, length(ia), length(ib))
    top_a <- if (nrow(a)) a$Description[order(a$p.adjust)][seq_len(min(10L, nrow(a)))] else character(0)
    top_b <- if (nrow(b)) b$Description[order(b$p.adjust)][seq_len(min(10L, nrow(b)))] else character(0)
    data.frame(
      arm = arm, ontology = ont,
      n_annotated_query_with_iea = if (nrow(a)) a$n_query_annotated[1] else NA_integer_,
      n_annotated_query_no_iea   = if (nrow(b)) b$n_query_annotated[1] else NA_integer_,
      n_annotated_bg_with_iea = MAP_STATS[[paste(ont, "with_iea", sep = "|")]]$n_background_annotated,
      n_annotated_bg_no_iea   = MAP_STATS[[paste(ont, "no_iea",   sep = "|")]]$n_background_annotated,
      n_terms_with_iea = length(ia), n_terms_no_iea = length(ib),
      n_shared = length(intersect(ia, ib)),
      jaccard = if (length(union(ia, ib))) length(intersect(ia, ib)) / length(union(ia, ib)) else NA_real_,
      n_top10_shared = length(intersect(top_a, top_b)),
      top10_with_iea = paste(top_a, collapse = " | "),
      top10_no_iea   = paste(top_b, collapse = " | "),
      top10_identical_order = identical(top_a[seq_len(top_n)], top_b[seq_len(top_n)]))
  }))))

## ---------------------------------------------------------------------------
## 10. Mouse-space replicate (secondary)
##
## The same question asked BEFORE the ortholog map, so a term that appears only in human
## space can be attributed to the projection rather than to the biology. Secondary by
## construction: the arms are defined in human projection space.
## ---------------------------------------------------------------------------

MR <- GO$mouse_replicate
MOUSE_REPLICATE <- data.frame(); MOUSE_TERMS <- data.frame()
mouse_sig_path <- file.path(PROJECT_ROOT, MR$signature_object)
mouse_uni_path <- file.path(PROJECT_ROOT, MR$universe)
if (file.exists(mouse_sig_path) && file.exists(mouse_uni_path)) {
  sig <- readRDS(mouse_sig_path)
  mq <- sig$sets[[MR$contrast]][[MR$direction]][[MR$gate]]
  mu <- read_symbols(mouse_uni_path)
  message("[10] mouse-space replicate: ", length(mq), " query, ", length(mu), " universe")
  GOMAP_MM <- build_go_map(as.character(MR$orgdb), mu, "mouse")
  for (ont in ONTS) {
    t2g_m <- term2gene(GOMAP_MM, ont, IEA_PRIMARY, mu)
    e <- run_ora(mq, t2g_m, mu)
    d <- if (is.null(e)) data.frame() else as.data.frame(e)
    if (nrow(d)) {
      d$arm <- paste0(MR$contrast, "_", MR$direction, "__mouse_space")
      d$ontology <- ont; d$iea_variant <- IEA_PRIMARY
      d$n_query_annotated <- as.integer(sub(".*/", "", d$GeneRatio))
      MOUSE_TERMS <- dplyr::bind_rows(MOUSE_TERMS, d)
    }
    hs <- ORA[[paste(primary_arm, ont, IEA_PRIMARY, sep = "|")]]
    ih <- if (nrow(hs)) hs$ID else character(0)
    im <- if (nrow(d)) d$ID else character(0)
    MOUSE_REPLICATE <- dplyr::bind_rows(MOUSE_REPLICATE, data.frame(
      ontology = ont,
      mouse_query_nominal = length(mq),
      mouse_query_annotated = length(intersect(mq, unique(t2g_m$gene))),
      mouse_universe = length(mu),
      mouse_universe_annotated = length(unique(t2g_m$gene)),
      n_terms_mouse_space = length(im),
      human_arm = primary_arm, n_terms_human_space = length(ih),
      n_shared = length(intersect(ih, im)),
      jaccard = if (length(union(ih, im))) length(intersect(ih, im)) / length(union(ih, im)) else NA_real_,
      top10_mouse = paste(if (nrow(d)) d$Description[order(d$p.adjust)][seq_len(min(10L, nrow(d)))] else character(0),
                          collapse = " | "),
      top10_human = paste(if (nrow(hs)) hs$Description[order(hs$p.adjust)][seq_len(min(10L, nrow(hs)))] else character(0),
                          collapse = " | ")))
  }
} else {
  message("[10] mouse-space replicate skipped: input missing")
}

## ---------------------------------------------------------------------------
## 11. Summaries, provenance, and the writes
## ---------------------------------------------------------------------------

TERM_SUMMARY <- dplyr::bind_rows(lapply(names(ORA), function(nm) {
  parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
  df <- ORA[[nm]]
  ne <- NULLS[[nm]]
  ms <- MAP_STATS[[paste(parts[2], parts[3], sep = "|")]]
  data.frame(
    arm = parts[1], ontology = parts[2], iea_variant = parts[3],
    is_primary = parts[2] == ONT_PRIMARY && parts[3] == IEA_PRIMARY,
    n_arm_genes = ARMS[[parts[1]]]$n_nominal,
    n_arm_annotated = if (nrow(df)) df$n_query_annotated[1]
                      else length(intersect(ARMS[[parts[1]]]$genes,
                                            unique(T2G[[paste(parts[2], parts[3], sep = "|")]]$gene))),
    n_background = ms$n_background, n_background_annotated = ms$n_background_annotated,
    n_terms_testable = if (!is.null(ne)) ne$n_terms_testable else NA_integer_,
    n_terms_in_family = if (!is.null(ne)) ne$n_observed_family else NA_integer_,
    n_enriched = nrow(df),
    n_after_simplify = length(simplify_keep[[nm]]),
    n_blocks = length(unique(CLUSTERS[[nm]]$cluster_id)),
    null_median_n_signif = if (!is.null(ne)) ne$summary$median_n_signif[ne$summary$null == "depth_matched"] else NA_real_,
    null_q95_n_signif = if (!is.null(ne)) ne$summary$q95_n_signif[ne$summary$null == "depth_matched"] else NA_real_,
    null_max_n_signif = if (!is.null(ne)) ne$summary$max_n_signif[ne$summary$null == "depth_matched"] else NA_real_,
    p_n_signif_matched = if (!is.null(ne)) ne$summary$p_n_signif[ne$summary$null == "depth_matched"] else NA_real_,
    uniform_median_n_signif = if (!is.null(ne)) ne$summary$median_n_signif[ne$summary$null == "uniform"] else NA_real_,
    n_genes_covered = if (!is.null(ne)) ne$n_covered_obs else NA_integer_,
    null_median_covered = if (!is.null(ne)) ne$summary$median_covered[ne$summary$null == "depth_matched"] else NA_real_,
    null_q95_covered = if (!is.null(ne)) ne$summary$q95_covered[ne$summary$null == "depth_matched"] else NA_real_,
    p_covered_matched = if (!is.null(ne)) ne$summary$p_covered[ne$summary$null == "depth_matched"] else NA_real_,
    uniform_median_covered = if (!is.null(ne)) ne$summary$median_covered[ne$summary$null == "uniform"] else NA_real_)
}))
TERM_SUMMARY$arm <- factor(TERM_SUMMARY$arm, levels = ARM_ORDER)
TERM_SUMMARY <- TERM_SUMMARY[order(TERM_SUMMARY$arm, TERM_SUMMARY$ontology,
                                   TERM_SUMMARY$iea_variant), ]
TERM_SUMMARY$arm <- as.character(TERM_SUMMARY$arm)

NULL_DESIGN <- dplyr::bind_rows(lapply(names(NULLS), function(nm) {
  ne <- NULLS[[nm]]; if (is.null(ne)) return(NULL)
  parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
  b <- ne$bands
  b$arm <- parts[1]; b$ontology <- parts[2]; b$iea_variant <- parts[3]
  b$n_replicates <- N_NULL; b$seed <- SEED
  b
}))
NULL_SUMMARY <- dplyr::bind_rows(lapply(names(NULLS), function(nm) {
  ne <- NULLS[[nm]]; if (is.null(ne)) return(NULL)
  parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
  s <- ne$summary
  s$arm <- parts[1]; s$ontology <- parts[2]; s$iea_variant <- parts[3]
  s$n_observed_enriched <- ne$n_observed_signif
  s$n_query_annotated <- ne$n_query_annotated
  s
}))
NULL_DESIGN <- dplyr::bind_rows(
  cbind(NULL_DESIGN, row_kind = "band"),
  cbind(NULL_SUMMARY, row_kind = "null_summary"))

## Every replicate of both nulls, so the two summary statistics on the figure can be
## checked against the distribution they came from rather than taken on trust.
NULL_DRAWS <- dplyr::bind_rows(lapply(names(NULLS), function(nm) {
  ne <- NULLS[[nm]]; if (is.null(ne)) return(NULL)
  parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
  d <- ne$draws
  d$arm <- parts[1]; d$ontology <- parts[2]; d$iea_variant <- parts[3]
  d$n_signif_obs <- ne$n_observed_signif
  d$n_covered_obs <- ne$n_covered_obs
  d$n_query_annotated <- ne$n_query_annotated
  d$seed <- SEED
  d
}))

DEPTH_NULL <- dplyr::bind_rows(lapply(names(NULLS), function(nm) {
  ne <- NULLS[[nm]]; if (is.null(ne)) return(NULL)
  parts <- strsplit(nm, "|", fixed = TRUE)[[1]]
  df <- ORA[[nm]]
  if (!nrow(df)) return(NULL)
  out <- ne$per_term[ne$per_term$ID %in% df$ID, ]
  out$arm <- parts[1]; out$ontology <- parts[2]; out$iea_variant <- parts[3]
  out$Description <- df$Description[match(out$ID, df$ID)]
  out$p_adjust_hypergeometric <- df$p.adjust[match(out$ID, df$ID)]
  out$simplify_kept <- df$simplify_kept[match(out$ID, df$ID)]
  out[order(out$p_adjust_hypergeometric), ]
}))

TERM_CLUSTERS <- dplyr::bind_rows(CLUSTERS)

PROVENANCE <- data.frame(
  key = c("script", "stage", "run_date", "R_version",
          "GOSemSim_version", "GOSemSim_path", "project_rlib",
          "GO.db_version", "org.Hs.eg.db_version", "org.Mm.eg.db_version",
          "clusterProfiler_version", "enrichplot_version", "DOSE_version",
          "background_file", "background_n",
          "ontologies", "primary_ontology", "iea_variants", "primary_iea",
          "drop_evidence", "p_cutoff", "q_cutoff", "p_adjust_method",
          "min_gs_size", "max_gs_size",
          "sim_measure", "simplify_cutoff", "cluster_linkage", "cluster_height",
          "n_null", "null_ontologies", "depth_band_quantiles", "seed",
          "seam_enricher_vs_enrichGO_same_ids", "seam_enricher_vs_enrichGO_max_abs_dp",
          "seam_null_matrix_max_abs_dk", "seam_null_matrix_max_abs_dq"),
  value = c(SCRIPT, STAGE, format(Sys.Date()), R.version.string,
            GSS_VER, find.package("GOSemSim"),
            normalizePath(DIR_R_LIBRARY, mustWork = FALSE),
            as.character(utils::packageVersion("GO.db")),
            as.character(utils::packageVersion("org.Hs.eg.db")),
            as.character(utils::packageVersion("org.Mm.eg.db")),
            as.character(utils::packageVersion("clusterProfiler")),
            as.character(utils::packageVersion("enrichplot")),
            as.character(utils::packageVersion("DOSE")),
            GO$background, length(BACKGROUND),
            paste(ONTS, collapse = ","), ONT_PRIMARY,
            paste(IEA_VARS, collapse = ","), IEA_PRIMARY,
            paste(DROP_EVID, collapse = ","), P_CUT, Q_CUT, PADJ_METHOD,
            MIN_GS, MAX_GS,
            SIM_MEASURE, SIMP_CUT, CLUST_LINK, CLUST_H,
            N_NULL, paste(NULL_ONTS, collapse = ","),
            paste(BAND_Q, collapse = ","), SEED,
            as.character(seam_same_ids), format(seam_max_dp),
            format(max(NULL_SEAM$max_abs_dk, na.rm = TRUE)),
            format(max(NULL_SEAM$max_abs_dq, na.rm = TRUE))))

w <- function(df, name) {
  p <- file.path(TBL_OVW, paste0(name, ".csv"))
  if (is.null(df) || !nrow(df))
    df <- data.frame(note = sprintf("%s produced no rows in this run", name))
  utils::write.csv(df, p, row.names = FALSE)
  message("  [write] ", p, "  (", nrow(df), " rows)")
}
w(ENRICHED,         "go_enriched_terms")
w(TERM_SUMMARY,     "go_term_summary")
w(IEA_COMPARISON,   "go_iea_comparison")
w(DEPTH_NULL,       "go_depth_null")
w(NULL_DESIGN,      "go_null_design")
w(NULL_DRAWS,       "go_null_draws")
w(TERM_CLUSTERS,    "go_term_clusters")
w(GO_VS_CURATED,    "go_vs_curated_lenses")
w(GENE_COVERAGE,    "go_gene_coverage")
w(COVERAGE_SUMMARY, "go_coverage_summary")
w(PROTEOSTASIS,     "go_proteostasis_probe")
w(MOUSE_REPLICATE,  "go_mouse_replicate")
w(PROVENANCE,       "go_provenance")

OBJ <- list(
  stage = STAGE, script = SCRIPT,
  arms = lapply(ARMS, function(a) a[c("name", "role", "file", "n_nominal",
                                      "n_in_background")]),
  arm_order = ARM_ORDER, primary_arm = primary_arm,
  background_n = length(BACKGROUND),
  ontology_primary = ONT_PRIMARY, iea_primary = IEA_PRIMARY,
  enriched = ENRICHED, term_summary = TERM_SUMMARY, iea_comparison = IEA_COMPARISON,
  depth_null = DEPTH_NULL, null_design = NULL_DESIGN, null_draws = NULL_DRAWS,
  term_clusters = TERM_CLUSTERS,
  go_vs_curated = GO_VS_CURATED, gene_coverage = GENE_COVERAGE,
  coverage_summary = COVERAGE_SUMMARY, proteostasis = PROTEOSTASIS,
  mouse_replicate = MOUSE_REPLICATE, mouse_terms = MOUSE_TERMS,
  provenance = PROVENANCE, gosemsim_version = GSS_VER,
  cutoffs = list(p = P_CUT, q = Q_CUT, padj = PADJ_METHOD,
                 min_gs = MIN_GS, max_gs = MAX_GS),
  top_n_terms = as.integer(GO$top_n_terms %||% 18L),
  top_n_blocks = as.integer(GO$top_n_blocks %||% 8L))
saveRDS(OBJ, file.path(DIR_OBJECTS, "24_go_decomposition.rds"))
message("[SAVE] ", file.path(DIR_OBJECTS, "24_go_decomposition.rds"))

message("=================================================================")
message("24_go_decomposition: done. ", nrow(ENRICHED), " enriched-term rows across ",
        length(ARM_ORDER), " arms x ", length(ONTS), " ontologies x ",
        length(IEA_VARS), " IEA variants.")
message("=================================================================")
