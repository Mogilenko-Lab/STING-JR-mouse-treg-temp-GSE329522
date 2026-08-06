#!/usr/bin/env Rscript
# 20_semantic_decomposition.R -- COMPUTE
# =============================================================================
# Ontology-side composition of the mouse 39 °C-derived up arm (stage 13_semantic_decomp).
#
# Purpose
#   Stage 12 asked how much of WT_heat_up falls inside curated lenses, a binary
#   membership question answering 22 of 213 in, 191 out. This stage asks three further
#   questions of the same 213 genes using GO BP semantic similarity (GOSemSim, Wang graph
#   measure, BMA gene-level combiner):
#
#     1. Can the ontology see the set at all? An annotation coverage ladder that names
#        the genes it cannot see.
#     2. Is the set one coherent process, or many? Mean pairwise semantic similarity
#        against annotation-depth-matched nulls and single-GO-term reference sets, all
#        size-matched.
#     3. Can it be split into modules? A silhouette sweep over linkage and k, run on every
#        set including reference sets coherent by construction, so a failure to partition
#        can be attributed to the method or to the set.
#
#   Then, for the genes no lens contains: how close each one's annotation sits to each
#   lens, graded, against a depth-stratified background.
#
#   Nothing here reads expression. The similarity structure is a property of the gene set
#   and the ontology, which keeps every partition computed here independent of enrichment.
#
# The confound this stage is built around
#   The export gate thresholds effect size, which selects well-studied genes. WT_heat_up
#   members carry a median 11 GO BP terms against the expressed background's 5 (Wilcoxon
#   p 5e-22). A better-annotated gene has more term pairs and so more chances to sit near
#   any reference set, and a uniform-random null would credit the set with coherence it
#   gets from annotation depth alone. Every null here matches the query's terms-per-gene
#   distribution band by band, and every percentile is taken within a gene's own band. The
#   uniform null is computed as well and reported beside the matched one to size the
#   confound.
#
# Inputs (frozen; read only)
#   03_results/objects/17_signature_sets.rds
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_mouse_lens.rds
#   00_data/references/gene_sets/tcr_activation_lens/tcr_activation_mouse.rds
#   00_data/references/gene_sets/lombardi2022_hif_consensus_mouse.rds
#   org.Mm.eg.db + GO.db (versions recorded in tables/semantic_provenance.csv)
#
# Outputs
#   03_results/13_semantic_decomp/tables/semantic_coverage.csv
#   03_results/13_semantic_decomp/tables/semantic_coherence_sets.csv
#   03_results/13_semantic_decomp/tables/lens_proximity_per_gene.csv
#   03_results/13_semantic_decomp/tables/semantic_provenance.csv
#   03_results/13_semantic_decomp/tables/_overview/wtheatup_lens_proximity.csv
#   03_results/objects/20_semantic_decomp.rds
#   03_results/objects/20_semantic_decomp__GOSemSim-<version>.rds   (engine-stamped copy)
#
# Honest ceiling
#   Semantic proximity to a lens says a gene's GO annotation resembles that of the lens
#   members. Participation in the named process and co-regulation are further claims that
#   need other evidence. GO BP coverage is also uneven in a way tied to the biology here:
#   P4ha2, one of the seven Lombardi2022_HIF genes in the set, carries CC and MF
#   annotation and no BP annotation in either species, so this stage leaves it unplaced.
#
# Which GOSemSim build this runs on, and why it is pinned
#   Wang similarity weights a parent edge by its relationship: is_a 0.8, part_of 0.6,
#   everything else 0.7. GOSemSim through 2.36.0 matches those weight names against the
#   relationship strings in its own packaged gotbl without normalising them, and gotbl
#   spells the same relationships "isa" and "part of". Every edge falls through to
#   "other", so on BP all 66,947 edges run at a uniform 0.7 — the distinction the Wang
#   measure is made of is absent. 2.39.2 (Bioconductor devel, upstream issue #51)
#   normalises the spellings first. This stage runs on 2.39.2 from the project-local R
#   library. Install it with:
#     R CMD INSTALL -l 01_modules/Rlib <GOSemSim source carrying the fix>
#
#   Two guards, since a version string is the weaker claim. use_project_rlib() stops if
#   GOSemSim resolves outside the project library, and section 0 recovers the effective
#   is_a weight from the measure and stops unless it is 0.8. The recovered value goes into
#   semantic_provenance.csv and into the term-matrix cache key, so the published table
#   shows the weighting that produced it. 20b_semantic_engine_validation.R measures how
#   much the difference is worth.
#
# Run from project root (expensive; one shared term matrix, then indexing):
#   Rscript 02_analysis/scripts/20_semantic_decomposition.R
#   PROJECT_RLIB=none Rscript 02_analysis/scripts/20_semantic_decomposition.R   # system build
# =============================================================================

source("02_analysis/config/config.R")
use_project_rlib(require_pkgs = "GOSemSim")

suppressPackageStartupMessages({
  library(GOSemSim)
  library(org.Mm.eg.db)
  library(GO.db)
  library(AnnotationDbi)
  library(parallel)
  library(cluster)
  library(dplyr)
})
options(stringsAsFactors = FALSE)

STAGE  <- "13_semantic_decomp"
SCRIPT <- "02_analysis/scripts/20_semantic_decomposition.R"
TBL_DIR <- stage_dir(STAGE, "tables")
TBL_OVW <- file.path(TBL_DIR, "_overview")
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_OVW, recursive = TRUE, showWarnings = FALSE)

SEM <- YAML_CONFIG$semantic
ONT       <- SEM$ontology          %||% "BP"
MEASURE   <- SEM$measure           %||% "Wang"
COMBINE   <- SEM$combine           %||% "BMA"
SEED      <- SEM$seed              %||% 20260728L
BANDS     <- c(as.numeric(unlist(SEM$depth_bands %||% c(0, 3, 6, 10, 15, 22, 32, 50))), Inf)
N_UNIFORM <- as.integer(SEM$n_null_uniform     %||% 20L)
N_MATCHED <- as.integer(SEM$n_null_matched     %||% 20L)
N_THEMESUB <- as.integer(SEM$n_theme_subsamples %||% 5L)
BG_N      <- as.integer(SEM$bg_sample_size     %||% 2000L)
MAX_K     <- as.integer(SEM$max_k_partition    %||% 15L)
THEMES    <- as.character(unlist(SEM$theme_references))
DROP_CODES <- as.character(unlist(SEM$drop_evidence %||% "IEA"))
NCORE     <- max(1L, min(60L, parallel::detectCores() - 4L))
GSS_VER   <- as.character(utils::packageVersion("GOSemSim"))
EXPECT_ISA <- as.numeric(SEM$expect_isa_weight %||% 0.8)
KNOWN_FIXED <- as.character(unlist(SEM$known_fixed_gosemsim %||% character(0)))

## ---------------------------------------------------------------------------
## 0. Which weighting is actually in force
## ---------------------------------------------------------------------------
# A version string records which build was loaded, and that is a weaker claim than the
# one this stage needs. 2.36.0 keys the Wang weight table on "is_a"/"part_of" and matches
# it against the packaged gotbl, which spells the same relationships "isa"/"part of", so
# every BP edge silently takes the 0.7 "other" weight. The stamped version would still
# read 2.36.0 and look self-consistent, and a future build could regress the same way
# under a version this script has never seen. So recover the weight from the arithmetic.
#
# For a BP term whose only gotbl parent is an is_a edge to the BP root, the two DAGs are
# {child: 1, root: w, all: 0} and {root: 1, all: 0}, the root and "all" are the shared
# terms, and Wang similarity reduces to (1 + w) / (2 + w). Inverting, w = (2s - 1)/(1 - s).
# Terms carrying extra parents only dilute w downward, so the maximum over several
# candidates recovers the true weight and needs no hardcoded GO id.
measured_isa_weight <- function(semdata) {
  e <- new.env()
  utils::data("gotbl", package = "GOSemSim", envir = e)
  g <- get("gotbl", envir = e)
  g <- unique(g[g$Ontology == ONT, c("go_id", "relationship", "parent")])
  root <- switch(ONT, BP = "GO:0008150", MF = "GO:0003674", CC = "GO:0005575")
  npar <- table(g$go_id)
  cand <- unique(g$go_id[g$parent == root & g$relationship == "isa" &
                         npar[g$go_id] == 1L])
  if (length(cand) < 1L) return(NA_real_)
  s <- vapply(utils::head(sort(cand), 8L), function(x)
    suppressWarnings(GOSemSim::goSim(x, root, semData = semdata, measure = MEASURE)),
    numeric(1))
  w <- (2 * s - 1) / (1 - s)
  w <- w[is.finite(w) & w > 0.4 & w < 1.05]
  if (!length(w)) return(NA_real_)
  max(w)
}

message("=================================================================")
message("20_semantic_decomposition: GO ", ONT, " composition (", STAGE, ")")
message("  measure: ", MEASURE, " / ", COMBINE, "   cores: ", NCORE, "   seed: ", SEED)
message("  GOSemSim ", GSS_VER, " from ", find.package("GOSemSim"))
message("=================================================================")

## ---------------------------------------------------------------------------
## 1. Annotation substrate
## ---------------------------------------------------------------------------
semdata <- godata(annoDb = "org.Mm.eg.db", keytype = "SYMBOL", ont = ONT,
                  computeIC = TRUE)

ISA_W <- if (identical(MEASURE, "Wang")) measured_isa_weight(semdata) else NA_real_
if (identical(MEASURE, "Wang")) {
  message("[0] effective Wang is_a edge weight, recovered from the measure: ",
          sprintf("%.3f", ISA_W), " (expected ", EXPECT_ISA, ")")
  if (!is.finite(ISA_W))
    stop("[20] could not recover the Wang is_a weight; the guard cannot be skipped ",
         "silently. Check GOSemSim's gotbl and the ", ONT, " root id.")
  isa_ok <- abs(ISA_W - EXPECT_ISA) <= 0.01
  if (!isa_ok) {
    msg <- sprintf(paste0("effective Wang is_a weight is %.3f, not the expected %.2f. ",
                          "GOSemSim %s is collapsing is_a/part_of onto the 0.7 'other' ",
                          "weight, so the measure has lost the distinction it is made of. ",
                          "Install a fixed build (known good: %s) into %s."),
                   ISA_W, EXPECT_ISA, GSS_VER,
                   paste(KNOWN_FIXED, collapse = ", "), DIR_R_LIBRARY)
    # PROJECT_RLIB is the deliberate escape hatch used to regenerate the pre-fix numbers
    # for 20b's comparison, so there it warns; on a default run it stops.
    if (nzchar(Sys.getenv("PROJECT_RLIB")))
      warning("[20] ", msg, " Continuing because PROJECT_RLIB is set.")
    else
      stop("[20] ", msg)
  }
  # Only meaningful once the weight check has passed: a build outside the advisory list
  # that satisfies the behavioural guard is fine, and the list is what should be updated.
  if (isa_ok && length(KNOWN_FIXED) && !GSS_VER %in% KNOWN_FIXED)
    message("[0] GOSemSim ", GSS_VER, " passed the behavioural guard and is not yet in ",
            "semantic.known_fixed_gosemsim (", paste(KNOWN_FIXED, collapse = ", "),
            "); add it there.")
}

IC <- semdata@IC
raw_anno <- semdata@geneAnno[!is.na(semdata@geneAnno$GO) &
                             semdata@geneAnno$ONTOLOGY == ONT, ]
# GOSemSim's own gene-level entry points drop electronically inferred annotations by
# default, and this stage follows that: an IEA term has never been seen by a curator.
# It is also the conservative choice for the confound this stage is built around,
# since poorly studied genes lean on IEA more heavily than well studied ones. The
# coverage table reports what the drop costs.
anno <- unique(raw_anno[!raw_anno$EVIDENCE %in% DROP_CODES, c("SYMBOL", "GO")])
anno_with_iea <- unique(raw_anno[, c("SYMBOL", "GO")])
GT <- split(anno$GO, anno$SYMBOL)
ANNOT <- names(GT)
NT <- vapply(GT, length, integer(1))
message("[1] genes with curated ", ONT, " annotation: ", length(ANNOT),
        " (", length(unique(anno_with_iea$SYMBOL)), " if ",
        paste(DROP_CODES, collapse = "/"), " kept)   distinct terms: ",
        length(unique(anno$GO)))

sig <- readRDS(file.path(DIR_OBJECTS, "17_signature_sets.rds"))
stopifnot(!is.null(sig$sets$WT_heat$up$fdr_logfc))
GATE_SETS <- list(
  WT_heat_up   = sig$sets$WT_heat$up$fdr_logfc,
  WT_heat_down = sig$sets$WT_heat$down$fdr_logfc,
  KO_heat_up   = sig$sets$KO_heat$up$fdr_logfc
)
query_nominal <- GATE_SETS$WT_heat_up
query <- sort(intersect(query_nominal, ANNOT))
background <- setdiff(intersect(sig$sets$WT_heat$ranked$gene_symbol, ANNOT),
                      query_nominal)
message("[1] query ", length(query_nominal), " nominal -> ", length(query),
        " annotated; background pool ", length(background))

hsr_path <- file.path("00_data/references/gene_sets/temp_hsr_lens",
                      "temp_hsr_mouse_lens.rds")
tcr_path <- file.path("00_data/references/gene_sets/tcr_activation_lens",
                      "tcr_activation_mouse.rds")
hif_path <- file.path("00_data/references/gene_sets",
                      "lombardi2022_hif_consensus_mouse.rds")
for (p in c(hsr_path, tcr_path, hif_path)) stopifnot(file.exists(p))
.flat <- function(x) unique(as.character(if (is.list(x)) unlist(x, use.names = FALSE) else x))
hsr <- readRDS(hsr_path); tcr <- readRDS(tcr_path); hif <- readRDS(hif_path)
LENSES <- list(HSR_core         = .flat(hsr$HSR_core),
               TCR_activation   = .flat(tcr$TCR_activation),
               Lombardi2022_HIF = .flat(hif))
for (n in names(LENSES))
  message("[1] lens ", n, ": ", length(LENSES[[n]]), " genes, ",
          length(intersect(LENSES[[n]], ANNOT)), " annotated")

## ---------------------------------------------------------------------------
## 2. Coverage ladder — name what the ontology cannot see
## ---------------------------------------------------------------------------
missing <- setdiff(query_nominal, ANNOT)
all_go <- suppressMessages(AnnotationDbi::select(
  org.Mm.eg.db, keys = missing, keytype = "SYMBOL",
  columns = c("ONTOLOGY", "GENENAME")))
onts_by_gene <- tapply(all_go$ONTOLOGY, all_go$SYMBOL,
                       function(x) paste(sort(unique(na.omit(x))), collapse = "/"))
gname <- tapply(all_go$GENENAME, all_go$SYMBOL, function(x) na.omit(x)[1])
class_missing <- function(g) {
  o <- onts_by_gene[[g]] %||% ""
  predicted <- grepl("^Gm[0-9]+$|Rik$|-ps[0-9]*$|pseudogene", g) ||
    grepl("pseudogene|predicted gene", gname[[g]] %||% "", ignore.case = TRUE)
  if (nzchar(o) && o != "") return("other_ontology_only")
  if (predicted) return("predicted_or_pseudogene")
  "no_GO_annotation"
}
miss_class <- vapply(missing, class_missing, character(1))
# A predicted gene that nonetheless carries CC/MF is reported by its annotation, not
# its name, so the classes stay mutually exclusive and sum to the shortfall.
coverage <- dplyr::bind_rows(
  data.frame(step = 1L, label = "WT_heat_up, nominal",
             n_genes = length(query_nominal), denominator = length(query_nominal),
             genes = NA_character_),
  data.frame(step = 2L, label = sprintf("with GO %s annotation", ONT),
             n_genes = length(query), denominator = length(query_nominal),
             genes = NA_character_),
  data.frame(step = 3L, label = "no BP: carries only CC or MF terms",
             n_genes = sum(miss_class == "other_ontology_only"),
             denominator = length(query_nominal),
             genes = paste(sort(names(miss_class)[miss_class == "other_ontology_only"]),
                           collapse = "; ")),
  data.frame(step = 4L, label = "no BP: predicted gene or pseudogene",
             n_genes = sum(miss_class == "predicted_or_pseudogene"),
             denominator = length(query_nominal),
             genes = paste(sort(names(miss_class)[miss_class == "predicted_or_pseudogene"]),
                           collapse = "; ")),
  data.frame(step = 5L, label = "no BP: no GO annotation at all",
             n_genes = sum(miss_class == "no_GO_annotation"),
             denominator = length(query_nominal),
             genes = paste(sort(names(miss_class)[miss_class == "no_GO_annotation"]),
                           collapse = "; ")),
  # What the evidence-code filter costs, so the choice is visible and reversible.
  data.frame(step = 6L,
             label = sprintf("recovered if %s annotation were kept",
                             paste(DROP_CODES, collapse = "/")),
             n_genes = length(intersect(setdiff(query_nominal, ANNOT),
                                        unique(anno_with_iea$SYMBOL))),
             denominator = length(query_nominal),
             genes = paste(sort(intersect(setdiff(query_nominal, ANNOT),
                                          unique(anno_with_iea$SYMBOL))),
                           collapse = "; ")))
coverage$fraction <- coverage$n_genes / coverage$denominator
# Which curated-lens members the ontology cannot place. This is the honest ceiling made
# concrete: a lens overlap the euler draws can contain genes this stage must drop.
lens_lost <- lapply(names(LENSES), function(ln) {
  ov <- intersect(query_nominal, LENSES[[ln]])
  data.frame(lens = ln, n_overlap_nominal = length(ov),
             n_overlap_annotated = length(intersect(ov, ANNOT)),
             lost_genes = paste(sort(setdiff(ov, ANNOT)), collapse = "; "))
})
lens_lost <- dplyr::bind_rows(lens_lost)
message("[2] coverage: ", length(query), " of ", length(query_nominal),
        " annotated; lens members lost: ",
        paste(sprintf("%s %d", lens_lost$lens,
                      lens_lost$n_overlap_nominal - lens_lost$n_overlap_annotated),
              collapse = ", "))

## ---------------------------------------------------------------------------
## 3. Reference sets: single-GO-term themes, and the two nulls
## ---------------------------------------------------------------------------
go_offspring_genes <- function(goid) {
  kids <- tryCatch(unlist(AnnotationDbi::mget(goid, GO.db::GOBPOFFSPRING,
                                              ifnotfound = NA)),
                   error = function(e) NA)
  ids <- unique(na.omit(c(goid, kids)))
  ANNOT[vapply(GT, function(v) any(v %in% ids), logical(1))]
}
go_term_name <- function(goid)
  tryCatch(AnnotationDbi::Term(GO.db::GOTERM[[goid]]), error = function(e) goid)

theme_full <- lapply(THEMES, go_offspring_genes)
names(theme_full) <- THEMES
for (g in THEMES)
  message("[3] theme ", g, " ", go_term_name(g), ": ", length(theme_full[[g]]), " genes")

N_MATCH_SIZE <- length(query)
q_band <- cut(NT[query], BANDS)
b_band <- cut(NT[background], BANDS)
band_need <- table(q_band)

set.seed(SEED)
uniform_sets <- lapply(seq_len(N_UNIFORM), function(i)
  sort(sample(background, N_MATCH_SIZE)))
matched_sets <- lapply(seq_len(N_MATCHED), function(i) {
  sort(unique(unlist(lapply(levels(q_band), function(bd) {
    k <- band_need[[bd]]; pool <- background[b_band == bd]
    if (k == 0L || !length(pool)) return(character(0))
    sample(pool, min(k, length(pool)))
  }))))
})
theme_sets <- list()
for (g in THEMES) {
  full <- theme_full[[g]]
  nm <- sprintf("%s %s", g, go_term_name(g))
  if (length(full) <= N_MATCH_SIZE) {
    theme_sets[[sprintf("%s [all %d]", nm, length(full))]] <- sort(full)
  } else {
    for (j in seq_len(N_THEMESUB))
      theme_sets[[sprintf("%s [%d of %d, draw %d]", nm, N_MATCH_SIZE,
                          length(full), j)]] <- sort(sample(full, N_MATCH_SIZE))
  }
}
bg_sample <- sort(sample(background, min(BG_N, length(background))))

## ---------------------------------------------------------------------------
## 4. One shared term-similarity matrix over the whole analysis universe
## ---------------------------------------------------------------------------
universe <- sort(unique(c(query, unlist(GATE_SETS), unlist(LENSES),
                          unlist(theme_sets), unlist(uniform_sets),
                          unlist(matched_sets), bg_sample)))
universe <- intersect(universe, ANNOT)
terms <- sort(unique(unlist(GT[universe])))
gb <- length(terms)^2 * 8 / 1e9
message("[4] universe ", length(universe), " genes -> ", length(terms),
        " terms (", sprintf("%.2f", gb), " GB matrix)")
if (gb > 24) stop("[20] term matrix would exceed 24 GB; lower bg_sample_size or n_null_*.")

# The term matrix is the one expensive step and depends on (terms, measure, ontology,
# GO.db version, GOSemSim version). Cache it in the sanctioned throwaway zone so a rerun
# after a downstream change costs seconds. Any mismatch recomputes.
# The GOSemSim version belongs in both the key and the filename: the 2.36.0 and 2.39.2
# builds return different Wang values for the same term pair, so a key that omitted the
# engine would hand a rerun the previous build's matrix without saying so, and stamping
# the filename lets the two builds be compared without each recomputing the other's.
CACHE <- file.path(DIR_SCRATCH, sprintf("20_termsim_cache__GOSemSim-%s.rds", GSS_VER))
dir.create(DIR_SCRATCH, recursive = TRUE, showWarnings = FALSE)
cache_key <- list(terms = terms, measure = MEASURE, ont = ONT,
                  godb = as.character(utils::packageVersion("GO.db")),
                  gosemsim = GSS_VER,
                  isa_weight = if (is.finite(ISA_W)) round(ISA_W, 3) else NA_real_)
TS <- NULL
if (file.exists(CACHE)) {
  cached <- readRDS(CACHE)
  if (identical(cached$key, cache_key)) {
    TS <- cached$TS
    message("[4] reused cached term matrix (", basename(CACHE), ")")
  } else {
    message("[4] cache present but the key differs; recomputing")
  }
}
if (is.null(TS)) {
  t0 <- Sys.time()
  chunks <- split(seq_along(terms), cut(seq_along(terms), NCORE, labels = FALSE))
  blocks <- mclapply(chunks, function(ix)
    GOSemSim::termSim(terms[ix], terms, semData = semdata, method = MEASURE),
    mc.cores = NCORE)
  failed <- vapply(blocks, function(b) !is.matrix(b) || !nrow(b), logical(1))
  if (any(failed)) stop("[20] termSim failed in ", sum(failed), " of ", length(blocks),
                        " chunks: ", paste(utils::head(unlist(blocks[failed]), 2),
                                           collapse = " | "))
  # NA is preserved as NA: pair_bma drops all-missing rows and columns
  # the way GOSemSim does, and a zero would instead be treated as a real low score.
  TS <- matrix(NA_real_, length(terms), length(terms), dimnames = list(terms, terms))
  for (b in blocks) TS[rownames(b), colnames(b)] <- b
  diag(TS) <- 1
  rm(blocks); invisible(gc())
  message("[4] termSim computed in ", sprintf("%.1f",
          as.numeric(difftime(Sys.time(), t0, units = "mins"))), " min")
  saveRDS(list(key = cache_key, TS = TS), CACHE)
}
n_dead <- sum(rowSums(!is.na(TS)) <= 1L)
message("[4] terms comparable to no other term: ", n_dead, " of ", length(terms))

TIDX <- setNames(seq_along(terms), terms)
gidx <- lapply(GT[universe], function(v) unname(TIDX[v[v %in% terms]]))
names(gidx) <- universe
gidx <- gidx[vapply(gidx, length, integer(1)) > 0L]

# Reproduces GOSemSim:::combineScores(..., "BMA") exactly, verified against
# GOSemSim::mgeneSim to within its 3-digit rounding (max abs difference 4.9e-4 over a
# 12-gene test). Two details matter and neither is optional. BMA is the term-count
# WEIGHTED mean of the best matches, sum(rowMax, colMax) / sum(dim), which differs from the average
# of the two directional means; the two coincide only when both genes carry the same
# number of terms, and the whole point of this stage is that they do not. Rows and
# columns that are entirely missing are dropped before combining, over being treated
# as zero, or a term the ontology cannot compare would drag the score down.
pair_bma <- function(ia, ib) {
  M <- TS[ia, ib, drop = FALSE]
  rk <- !apply(M, 1, function(r) all(is.na(r)))
  ck <- !apply(M, 2, function(c) all(is.na(c)))
  M <- M[rk, ck, drop = FALSE]
  if (!length(M)) return(NA_real_)
  if (nrow(M) == 1L || ncol(M) == 1L) return(max(M, na.rm = TRUE))
  (sum(apply(M, 1, max, na.rm = TRUE)) +
   sum(apply(M, 2, max, na.rm = TRUE))) / sum(dim(M))
}

## ---------------------------------------------------------------------------
## 5. Set coherence + partition diagnostics
## ---------------------------------------------------------------------------
sim_matrix <- function(genes) {
  gg <- intersect(genes, names(gidx))
  if (length(gg) < 2L) return(matrix(numeric(0), 0, 0))
  ix <- gidx[gg]
  rows <- mclapply(seq_along(gg), function(i) vapply(seq_along(gg), function(j)
    if (j < i) 0 else pair_bma(ix[[i]], ix[[j]]), numeric(1)), mc.cores = NCORE)
  bad <- vapply(rows, function(r) !is.numeric(r) || length(r) != length(gg), logical(1))
  if (any(bad)) stop("[20] sim_matrix: ", sum(bad), " of ", length(gg),
                     " rows came back malformed (a forked worker likely died)")
  S <- do.call(rbind, rows)
  S[lower.tri(S)] <- t(S)[lower.tri(S)]
  dimnames(S) <- list(gg, gg)
  S
}

set_stats <- function(genes, label, family) {
  S <- sim_matrix(genes)
  gg <- rownames(S)
  if (length(gg) < 10L) {
    message("[5] skipped ", label, ": only ", length(gg), " annotated genes")
    return(NULL)
  }
  od <- S[upper.tri(S)]
  D <- as.dist(1 - S / max(S))
  best <- list(sil = -Inf, cfg = NA_character_, k = NA_integer_, big = NA_real_)
  kmax <- min(MAX_K, length(gg) - 1L)
  for (lk in c("average", "ward.D2", "complete")) {
    hc <- hclust(D, method = lk)
    for (k in 2:kmax) {
      cl <- cutree(hc, k = k)
      sw <- mean(cluster::silhouette(cl, D)[, "sil_width"])
      if (sw > best$sil) best <- list(sil = sw, cfg = lk, k = k,
                                      big = max(table(cl)) / length(gg))
    }
  }
  for (k in 2:kmax) {
    pm <- cluster::pam(D, k = k)
    if (pm$silinfo$avg.width > best$sil)
      best <- list(sil = pm$silinfo$avg.width, cfg = "pam", k = k,
                   big = max(table(pm$clustering)) / length(gg))
  }
  data.frame(label = label, family = family, n_used = length(gg),
             med_terms_per_gene = as.numeric(median(NT[gg])),
             mean_sim = mean(od), sd_sim = sd(od),
             q10_sim = as.numeric(quantile(od, 0.10)),
             q90_sim = as.numeric(quantile(od, 0.90)),
             best_sil = best$sil, best_cfg = best$cfg, best_k = best$k,
             largest_cluster_frac = best$big)
}

message("[5] set coherence: query, gate siblings, curated lenses, ",
        length(theme_sets), " theme draws, ", N_UNIFORM, " uniform + ",
        N_MATCHED, " depth-matched nulls")
rows <- list()
rows[[length(rows) + 1L]] <- set_stats(query, "WT_heat_up", "query")
for (nm in setdiff(names(GATE_SETS), "WT_heat_up"))
  rows[[length(rows) + 1L]] <- set_stats(GATE_SETS[[nm]], nm, "gate_sibling")
for (nm in names(LENSES))
  rows[[length(rows) + 1L]] <- set_stats(LENSES[[nm]], nm, "curated_lens")
for (nm in names(theme_sets))
  rows[[length(rows) + 1L]] <- set_stats(theme_sets[[nm]], nm, "single_go_term")
for (i in seq_along(uniform_sets))
  rows[[length(rows) + 1L]] <- set_stats(uniform_sets[[i]],
                                         sprintf("uniform_null_%02d", i), "null_uniform")
for (i in seq_along(matched_sets))
  rows[[length(rows) + 1L]] <- set_stats(matched_sets[[i]],
                                         sprintf("matched_null_%02d", i), "null_depth_matched")
coh <- dplyr::bind_rows(Filter(Negate(is.null), rows))

q_mean <- coh$mean_sim[coh$family == "query"][1]
nu <- coh$mean_sim[coh$family == "null_uniform"]
nm_ <- coh$mean_sim[coh$family == "null_depth_matched"]
scale_summary <- data.frame(
  quantity = c("query_mean_sim",
               "uniform_null_mean", "uniform_null_sd", "z_vs_uniform", "p_vs_uniform",
               "matched_null_mean", "matched_null_sd", "z_vs_matched", "p_vs_matched",
               "query_med_terms_per_gene", "matched_null_med_terms_per_gene",
               "uniform_null_med_terms_per_gene"),
  value = c(q_mean,
            mean(nu), sd(nu), (q_mean - mean(nu)) / sd(nu),
            (1 + sum(nu >= q_mean)) / (1 + length(nu)),
            mean(nm_), sd(nm_), (q_mean - mean(nm_)) / sd(nm_),
            (1 + sum(nm_ >= q_mean)) / (1 + length(nm_)),
            median(NT[query]),
            mean(coh$med_terms_per_gene[coh$family == "null_depth_matched"]),
            mean(coh$med_terms_per_gene[coh$family == "null_uniform"])))
message("[5] query mean_sim ", sprintf("%.4f", q_mean),
        " | uniform null ", sprintf("%.4f", mean(nu)),
        " | depth-matched null ", sprintf("%.4f", mean(nm_)))
message("[5] best silhouette across all ", nrow(coh), " sets: median ",
        sprintf("%.3f", median(coh$best_sil)), ", max ",
        sprintf("%.3f", max(coh$best_sil)))

## ---------------------------------------------------------------------------
## 6. Graded lens proximity, percentiles taken within a gene's own depth band
## ---------------------------------------------------------------------------
# lens_ix is a NAMED list, so vapply over it yields a named vector and which.max()
# carries the winning gene's name through — which would rename the third column to
# which_max.<gene> and break every downstream lookup. unname() before assembling.
# A gene is never scored against itself. Without this a lens member's nearest member
# is itself at similarity 1, every member pins to the ceiling by construction, and the
# members then count towards the very excess the panel is testing for. which_max is
# returned as an index into the FULL lens so the caller can name the nearest member.
to_lens <- function(g, lens_ix) {
  ig <- gidx[[g]]
  na3 <- c(max = NA_real_, mean = NA_real_, which_max = NA_real_)
  if (is.null(ig)) return(na3)
  keep <- names(lens_ix) != g
  if (!any(keep)) return(na3)
  v <- rep(NA_real_, length(lens_ix))
  v[keep] <- vapply(lens_ix[keep], function(li) pair_bma(ig, li), numeric(1))
  ok <- is.finite(v)
  if (!any(ok)) return(na3)
  c(max = unname(max(v[ok])), mean = unname(mean(v[ok])),
    which_max = unname(which.max(replace(v, !ok, -Inf))))
}

prox_rows <- list()
lens_summary <- list()
for (ln in names(LENSES)) {
  lg <- intersect(LENSES[[ln]], names(gidx))
  lens_ix <- gidx[lg]
  qv <- do.call(rbind, mclapply(query, to_lens, lens_ix = lens_ix, mc.cores = NCORE))
  # The background stands for "a gene this lens does not contain", so lens members are
  # taken out of it. Leaving them in would raise the bar using the very genes the
  # query's members are being compared against.
  bgg <- setdiff(intersect(bg_sample, names(gidx)), lg)
  bv <- do.call(rbind, mclapply(bgg, to_lens, lens_ix = lens_ix, mc.cores = NCORE))
  rownames(qv) <- query; rownames(bv) <- bgg
  qb <- cut(NT[query], BANDS); bb <- cut(NT[bgg], BANDS)
  # Mid-rank percentile, over ecdf. A gene carrying one or two curated terms either
  # shares a term with the lens, scoring exactly 1, or does not, so the top of the
  # sparse bands is a pile of ties. An ecdf hands every tied gene the full 1.0 and an
  # exceedance test hands it nothing; the mid-rank splits the tied mass and is the only
  # one of the three that does not depend on which side of the tie you stand.
  midrank <- function(ref, v) {
    ref <- ref[is.finite(ref)]
    if (!length(ref)) return(rep(NA_real_, length(v)))
    vapply(v, function(x) if (!is.finite(x)) NA_real_ else
      (sum(ref < x) + 0.5 * sum(ref == x)) / length(ref), numeric(1))
  }
  pct_pool  <- midrank(bv[, "max"], qv[, "max"])
  pct_band  <- rep(NA_real_, length(query))
  n_band_bg <- rep(NA_integer_, length(query))
  # A band is saturated when its own background already reaches the ceiling at the
  # threshold, so no gene in it can be distinguished from the top of its band. Those
  # genes are still drawn, and are excluded from the counted excess rather than
  # silently scored.
  band_sat <- rep(NA, length(query))
  for (bd in levels(qb)) {
    qi <- which(qb == bd); bi <- which(bb == bd)
    if (!length(qi)) next
    n_band_bg[qi] <- length(bi)
    if (length(bi) < 20L) next          # too thin to quantile honestly; left NA
    pct_band[qi] <- midrank(bv[bi, "max"], qv[qi, "max"])
    band_sat[qi] <- as.numeric(quantile(bv[bi, "max"], 0.95, na.rm = TRUE)) >= 1 - 1e-9
  }
  scored <- !is.na(pct_band) & !is.na(band_sat) & !band_sat
  nearest <- lg[qv[, "which_max"]]
  is_mem <- query %in% lg
  prox_rows[[ln]] <- data.frame(
    gene = query, lens = ln,
    max_to_lens = qv[, "max"], mean_to_lens = qv[, "mean"],
    nearest_lens_member = nearest,
    pct_uniform_background = pct_pool, pct_depth_matched_background = pct_band,
    depth_band = as.character(qb), n_background_in_band = n_band_bg,
    band_saturated = band_sat, counted = scored,
    n_terms = NT[query], is_lens_member = is_mem)
  lens_summary[[ln]] <- data.frame(
    lens = ln, n_lens_annotated = length(lg),
    n_members_in_query = sum(is_mem),
    n_scored = sum(scored),
    n_excluded_saturated_band = sum(!is.na(band_sat) & band_sat),
    n_above_p95 = sum(scored & pct_band >= 0.95),
    n_expected_by_chance = 0.05 * sum(scored),
    n_member_above_p95 = sum(scored & pct_band >= 0.95 & is_mem),
    n_nonmember_above_p95 = sum(scored & pct_band >= 0.95 & !is_mem),
    n_members_scored = sum(scored & is_mem),
    mean_pct_uniform = mean(pct_pool[scored], na.rm = TRUE),
    mean_pct_matched = mean(pct_band[scored], na.rm = TRUE),
    rho_max_vs_n_terms = suppressWarnings(cor(qv[, "max"], NT[query],
                                              method = "spearman",
                                              use = "complete.obs")))
  ls_ <- lens_summary[[ln]]
  message("[6] ", ln, ": ", ls_$n_members_in_query, " members; ", ls_$n_above_p95,
          " of ", ls_$n_scored, " scored genes above the depth-matched 95th percentile",
          " (chance ", sprintf("%.1f", ls_$n_expected_by_chance), "); ",
          ls_$n_excluded_saturated_band, " left uncounted in saturated bands")
}
prox <- dplyr::bind_rows(prox_rows)
lens_sum <- dplyr::bind_rows(lens_summary)

## ---------------------------------------------------------------------------
## 7. Provenance and writes
## ---------------------------------------------------------------------------
prov <- data.frame(
  key = c("R", "GOSemSim", "wang_isa_weight_measured",
          "org.Mm.eg.db", "GO.db", "ontology", "measure", "combine",
          "drop_evidence", "seed", "n_universe_genes", "n_terms", "n_terms_isolated",
          "query_nominal", "query_annotated", "background_pool", "bg_sample",
          "n_null_uniform", "n_null_matched", "n_theme_draws", "depth_bands"),
  value = c(paste(R.version$major, R.version$minor, sep = "."),
            as.character(utils::packageVersion("GOSemSim")),
            if (is.finite(ISA_W)) sprintf("%.3f", ISA_W) else NA_character_,
            as.character(utils::packageVersion("org.Mm.eg.db")),
            as.character(utils::packageVersion("GO.db")),
            ONT, MEASURE, COMBINE, paste(DROP_CODES, collapse = "/"),
            as.character(SEED),
            length(universe), length(terms), n_dead,
            length(query_nominal), length(query), length(background),
            length(bg_sample), N_UNIFORM, N_MATCHED, length(theme_sets),
            paste(BANDS, collapse = ",")))

write.csv(coverage, file.path(TBL_DIR, "semantic_coverage.csv"), row.names = FALSE)
write.csv(lens_lost, file.path(TBL_DIR, "semantic_lens_coverage_loss.csv"),
          row.names = FALSE)
write.csv(coh, file.path(TBL_DIR, "semantic_coherence_sets.csv"), row.names = FALSE)
write.csv(prox, file.path(TBL_DIR, "lens_proximity_per_gene.csv"), row.names = FALSE)
write.csv(prov, file.path(TBL_DIR, "semantic_provenance.csv"), row.names = FALSE)
write.csv(scale_summary, file.path(TBL_DIR, "semantic_scale_summary.csv"),
          row.names = FALSE)

Sq <- sim_matrix(query)
per_gene <- data.frame(
  gene = rownames(Sq),
  mean_sim_to_set = (rowSums(Sq) - diag(Sq)) / (nrow(Sq) - 1L),
  n_terms = NT[rownames(Sq)])
per_gene <- per_gene[order(-per_gene$mean_sim_to_set), ]
write.csv(per_gene, file.path(TBL_DIR, "query_per_gene_coherence.csv"), row.names = FALSE)

state <- list(query_sim = Sq, coverage = coverage, coherence = coh,
              proximity = prox, lens_summary = lens_sum, scale = scale_summary,
              per_gene = per_gene, lens_lost = lens_lost, provenance = prov,
              theme_names = setNames(vapply(THEMES, go_term_name, character(1)), THEMES),
              built_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"), built_by = SCRIPT)
saveRDS(state, file.path(DIR_OBJECTS, "20_semantic_decomp.rds"))
# A second copy stamped with the engine, so a run on one GOSemSim build does not
# overwrite the record of the other and 20b can put the two side by side.
stamped <- sprintf("20_semantic_decomp__GOSemSim-%s.rds", GSS_VER)
saveRDS(state, file.path(DIR_OBJECTS, stamped))

message("[7] wrote ", length(list.files(TBL_DIR, pattern = "[.]csv$")),
        " tables + objects/20_semantic_decomp.rds + objects/", stamped)
message("20_semantic_decomposition: COMPUTE DONE.")
