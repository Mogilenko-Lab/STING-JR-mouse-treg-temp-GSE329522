#!/usr/bin/env Rscript
# 20b_semantic_engine_validation.R -- COMPUTE
# =============================================================================
# Puts the two GOSemSim builds side by side and measures what changed (stage
# 13_semantic_decomp).
#
# Why the stage was re-run at all
#   Wang similarity weights each parent edge of the GO DAG by its relationship type: is_a
#   0.8, part_of 0.6, anything else 0.7. GOSemSim through 2.36.0 keys that weight table on
#   the strings "is_a" and "part_of" and matches it against GO.db's own relationship column
#   without normalising the spelling. GO.db writes those same relationships as "isa" and
#   "part of". Every edge therefore fails both matches and is remapped to "other", so the
#   whole DAG is traversed at a uniform 0.7 weight and the is_a/part_of distinction the Wang
#   measure is built on is absent from the result. In the GO BP parent table that is 44,495
#   isa and 4,772 part-of edges collapsed onto one weight. GOSemSim 2.39.2 normalises the
#   spellings first.
#
# What this script measures
#   1. Term pairs. Wang similarity for a fixed sample of GO BP term pairs under each build:
#      correlation, mean shift, how many values move at all.
#   2. Gene pairs. The same for gene-level BMA over a fixed sample of WT_heat_up members,
#      which is the combiner stage 20 uses.
#   3. The stage's own pairwise matrix. Correlation between the two builds' 196 x 196
#      WT_heat_up similarity matrices.
#   4. Every headline number the stage publishes, old beside new.
#
# Inputs (both must exist; the driver below says how to make either one)
#   03_results/objects/20_semantic_decomp__GOSemSim-<old>.rds
#   03_results/objects/20_semantic_decomp__GOSemSim-<new>.rds
#
# Outputs
#   03_results/13_semantic_decomp/tables/semantic_engine_validation.csv
#   03_results/_scratch/20b_probe__GOSemSim-<version>.rds   (per-build probe, reused)
#
# Run from project root, after stage 20 has been run on both builds:
#   PROJECT_RLIB=none Rscript 02_analysis/scripts/20_semantic_decomposition.R
#   Rscript 02_analysis/scripts/20_semantic_decomposition.R
#   Rscript 02_analysis/scripts/20b_semantic_engine_validation.R
#
# The probe in sections 1-2 needs both builds loaded, which exceeds what one R session can
# hold. The script therefore re-invokes itself once per build as a worker (the
# SEMANTIC_ENGINE_WORKER environment variable selects that mode) and reads the two probe
# files back.
# =============================================================================

source("02_analysis/config/config.R")

STAGE  <- "13_semantic_decomp"
SCRIPT <- "02_analysis/scripts/20b_semantic_engine_validation.R"

SEM      <- YAML_CONFIG$semantic
ONT      <- SEM$ontology %||% "BP"
MEASURE  <- SEM$measure  %||% "Wang"
COMBINE  <- SEM$combine  %||% "BMA"
SEED     <- SEM$seed     %||% 20260728L
N_TERM_PAIRS <- 4000L    # term pairs probed per build
N_GENE_PROBE <- 40L      # WT_heat_up members probed per build -> 780 gene pairs

probe_path <- function(v)
  file.path(DIR_SCRATCH, sprintf("20b_probe__GOSemSim-%s.rds", v))
stage_path <- function(v)
  file.path(DIR_OBJECTS, sprintf("20_semantic_decomp__GOSemSim-%s.rds", v))

## ---------------------------------------------------------------------------
## WORKER MODE -- one build, one probe file
## ---------------------------------------------------------------------------
if (nzchar(Sys.getenv("SEMANTIC_ENGINE_WORKER"))) {
  use_project_rlib()
  suppressPackageStartupMessages({
    library(GOSemSim); library(org.Mm.eg.db); library(GO.db)
  })
  v <- as.character(utils::packageVersion("GOSemSim"))
  message("[worker] GOSemSim ", v, " from ", find.package("GOSemSim"))

  semdata <- godata(annoDb = "org.Mm.eg.db", keytype = "SYMBOL", ont = ONT,
                    computeIC = FALSE)
  # Both workers draw from the same GO.db and the same seed, so the two probes
  # compare the same pairs term for term and gene for gene.
  set.seed(SEED)
  pool <- sort(unique(semdata@geneAnno$GO[semdata@geneAnno$ONTOLOGY == ONT]))
  t1 <- sample(pool, N_TERM_PAIRS, replace = TRUE)
  t2 <- sample(pool, N_TERM_PAIRS, replace = TRUE)
  ts <- vapply(seq_len(N_TERM_PAIRS), function(i)
    GOSemSim::goSim(t1[i], t2[i], semData = semdata, measure = MEASURE), numeric(1))

  sig   <- readRDS(file.path(DIR_OBJECTS, "17_signature_sets.rds"))
  qgene <- sort(intersect(sig$sets$WT_heat$up$fdr_logfc,
                          unique(semdata@geneAnno$SYMBOL[semdata@geneAnno$ONTOLOGY == ONT])))
  set.seed(SEED + 1L)
  gsub_ <- sort(sample(qgene, min(N_GENE_PROBE, length(qgene))))
  gm <- GOSemSim::mgeneSim(gsub_, semData = semdata, measure = MEASURE,
                           combine = COMBINE, verbose = FALSE)

  saveRDS(list(version = v, t1 = t1, t2 = t2, term_sim = ts,
               genes = gsub_, gene_sim = gm),
          probe_path(v))
  message("[worker] wrote ", probe_path(v))
  quit(save = "no", status = 0)
}

## ---------------------------------------------------------------------------
## DRIVER -- resolve the two builds
## ---------------------------------------------------------------------------
ver_of_lib <- function(env) {
  out <- system2("Rscript", c("-e", shQuote(paste(
    'lp <- Sys.getenv("PROJECT_RLIB");',
    'if (nzchar(lp) && !identical(tolower(lp),"none")) .libPaths(c(lp, .libPaths()));',
    'cat(as.character(packageVersion("GOSemSim")))'))),
    env = env, stdout = TRUE, stderr = FALSE)
  trimws(paste(out, collapse = ""))
}
ENV_SYS <- "PROJECT_RLIB=none"
ENV_PRJ <- paste0("PROJECT_RLIB=", DIR_R_LIBRARY)
V_OLD <- ver_of_lib(ENV_SYS)
V_NEW <- ver_of_lib(ENV_PRJ)
message("=================================================================")
message("20b_semantic_engine_validation: GOSemSim ", V_OLD, " vs ", V_NEW)
message("=================================================================")
if (identical(V_OLD, V_NEW))
  stop("Both library paths resolve to GOSemSim ", V_OLD,
       "; install the second build into ", DIR_R_LIBRARY, " before running this.")

for (p in c(stage_path(V_OLD), stage_path(V_NEW))) {
  if (!file.exists(p))
    stop("Missing ", p, ".\n  Produce it by running stage 20 on that build:\n",
         "    PROJECT_RLIB=none Rscript 02_analysis/scripts/20_semantic_decomposition.R\n",
         "    Rscript 02_analysis/scripts/20_semantic_decomposition.R")
}

run_worker <- function(env, v) {
  if (file.exists(probe_path(v))) {
    message("[probe] reusing ", basename(probe_path(v)))
    return(invisible(NULL))
  }
  message("[probe] running the ", v, " worker")
  st <- system2("Rscript", c(SCRIPT), env = c(env, "SEMANTIC_ENGINE_WORKER=1"))
  if (st != 0L) stop("the ", v, " probe worker exited with status ", st)
}
run_worker(ENV_SYS, V_OLD)
run_worker(ENV_PRJ, V_NEW)
P_OLD <- readRDS(probe_path(V_OLD))
P_NEW <- readRDS(probe_path(V_NEW))
stopifnot(identical(P_OLD$t1, P_NEW$t1), identical(P_OLD$t2, P_NEW$t2),
          identical(P_OLD$genes, P_NEW$genes))

## ---------------------------------------------------------------------------
## Row builder
## ---------------------------------------------------------------------------
rows <- list()
add <- function(section, quantity, old, new, unit = "", reading = "") {
  old <- as.numeric(old); new <- as.numeric(new)
  d <- new - old
  rows[[length(rows) + 1L]] <<- data.frame(
    section = section, quantity = quantity,
    engine_old = V_OLD, engine_new = V_NEW,
    old_value = old, new_value = new, delta = d,
    rel_delta = if (is.finite(old) && abs(old) > 1e-12) d / abs(old) else NA_real_,
    unit = unit,
    changed = if (!is.finite(d)) NA else abs(d) > 1e-9,
    reading = reading)
  invisible(NULL)
}

## ---------------------------------------------------------------------------
## 1-2. Probe comparisons
## ---------------------------------------------------------------------------
cmp_vec <- function(o, n, section, what) {
  k <- is.finite(o) & is.finite(n)
  add(section, sprintf("%s compared", what), sum(k), sum(k), "count")
  add(section, sprintf("%s with an identical value", what),
      sum(o[k] == n[k]), sum(o[k] == n[k]), "count",
      "one count over both builds; the same value in both columns")
  add(section, sprintf("mean %s similarity", what), mean(o[k]), mean(n[k]), MEASURE)
  add(section, sprintf("median %s similarity", what), median(o[k]), median(n[k]), MEASURE)
  add(section, sprintf("Pearson r, %s, old vs new", what),
      cor(o[k], n[k]), cor(o[k], n[k]), "r",
      "one correlation between the builds; the same number in both columns")
  add(section, sprintf("Spearman rho, %s, old vs new", what),
      cor(o[k], n[k], method = "spearman"), cor(o[k], n[k], method = "spearman"), "rho",
      "one correlation between the builds; the same number in both columns")
  d <- n[k] - o[k]
  add(section, sprintf("largest downward shift, %s", what), min(d), min(d), MEASURE, "new minus old")
  add(section, sprintf("largest upward shift, %s", what), max(d), max(d), MEASURE, "new minus old")
  invisible(NULL)
}
cmp_vec(P_OLD$term_sim, P_NEW$term_sim, "term_pairs",
        sprintf("GO %s term pairs", ONT))
u <- upper.tri(P_OLD$gene_sim)
cmp_vec(P_OLD$gene_sim[u], P_NEW$gene_sim[u], "gene_pairs",
        sprintf("WT_heat_up gene pairs (%s)", COMBINE))

## ---------------------------------------------------------------------------
## 3-4. The stage itself
## ---------------------------------------------------------------------------
O <- readRDS(stage_path(V_OLD)); N <- readRDS(stage_path(V_NEW))
pv <- function(x, k) setNames(x$provenance$value, x$provenance$key)[[k]]
sv <- function(x, k) setNames(x$scale$value, x$scale$quantity)[[k]]
for (k in c("ontology", "measure", "combine", "GO.db", "org.Mm.eg.db", "seed",
            "n_terms", "query_annotated"))
  if (!identical(pv(O, k), pv(N, k)))
    warning("provenance key '", k, "' differs between the two runs: ",
            pv(O, k), " vs ", pv(N, k), " -- the comparison is not engine-only")

# The full pairwise matrix. Both runs place the same genes in the same order, so the
# upper triangles line up element for element.
stopifnot(identical(rownames(O$query_sim), rownames(N$query_sim)))
uq <- upper.tri(O$query_sim)
qo <- O$query_sim[uq]; qn <- N$query_sim[uq]
add("query_matrix", "WT_heat_up gene pairs in the matrix", length(qo), length(qn), "count")
add("query_matrix", "Pearson r over the full pairwise matrix",
    cor(qo, qn), cor(qo, qn), "r", "one correlation between the builds")
add("query_matrix", "Spearman rho over the full pairwise matrix",
    cor(qo, qn, method = "spearman"), cor(qo, qn, method = "spearman"), "rho",
    "one correlation between the builds")
add("query_matrix", "pairs whose value moved", sum(qo != qn), sum(qo != qn), "count")
add("query_matrix", "mean of the full pairwise matrix", mean(qo), mean(qn), MEASURE)

# Coverage is a property of the annotation, so it should not move with the engine.
add("coverage", "WT_heat_up genes, nominal",
    O$coverage$n_genes[O$coverage$step == 1L], N$coverage$n_genes[N$coverage$step == 1L], "genes")
add("coverage", sprintf("with GO %s annotation on curated evidence", ONT),
    O$coverage$n_genes[O$coverage$step == 2L], N$coverage$n_genes[N$coverage$step == 2L], "genes")
p4 <- function(x) {
  r <- x$coverage[x$coverage$step == 3L, ]
  as.integer(grepl("(^|; )P4ha2(;|$)", r$genes[1]))
}
add("coverage", sprintf("P4ha2 carries no GO %s term (1 = confirmed)", ONT),
    p4(O), p4(N), "flag", "engine-independent; read from the coverage ladder")

fam <- function(x, f, col) x$coherence[[col]][x$coherence$family == f]
lab <- function(x, l, col) x$coherence[[col]][x$coherence$label == l][1]
add("coherence", "WT_heat_up mean pairwise similarity",
    sv(O, "query_mean_sim"), sv(N, "query_mean_sim"), MEASURE)
add("coherence", "depth-matched null, mean",
    sv(O, "matched_null_mean"), sv(N, "matched_null_mean"), MEASURE)
add("coherence", "depth-matched null, SD",
    sv(O, "matched_null_sd"), sv(N, "matched_null_sd"), MEASURE)
add("coherence", "WT_heat_up vs depth-matched null, z",
    sv(O, "z_vs_matched"), sv(N, "z_vs_matched"), "SD")
add("coherence", "WT_heat_up vs depth-matched null, p",
    sv(O, "p_vs_matched"), sv(N, "p_vs_matched"), "p",
    sprintf("permutation p over %s draws; floor 1/(n+1)", pv(N, "n_null_matched")))
add("coherence", "uniform null, mean",
    sv(O, "uniform_null_mean"), sv(N, "uniform_null_mean"), MEASURE)
add("coherence", "WT_heat_up vs uniform null, z",
    sv(O, "z_vs_uniform"), sv(N, "z_vs_uniform"), "SD")
add("coherence", "WT_heat_up vs uniform null, p",
    sv(O, "p_vs_uniform"), sv(N, "p_vs_uniform"), "p")
add("coherence", "single-GO-term references, family mean",
    mean(fam(O, "single_go_term", "mean_sim")), mean(fam(N, "single_go_term", "mean_sim")), MEASURE)
# The published range 0.295-0.358 is over per-theme means, one per GO term, so the
# five draws of an oversized theme are averaged before the range is taken.
theme_means <- function(x) {
  d <- x$coherence[x$coherence$family == "single_go_term", ]
  g <- trimws(sub("^(GO:[0-9]+ [^\\[]+)\\[.*$", "\\1", d$label))
  tapply(d$mean_sim, g, mean)
}
add("coherence", "single-GO-term references, lowest per-theme mean",
    min(theme_means(O)), min(theme_means(N)), MEASURE)
add("coherence", "single-GO-term references, highest per-theme mean",
    max(theme_means(O)), max(theme_means(N)), MEASURE)
for (l in c("TCR_activation", "HSR_core", "Lombardi2022_HIF"))
  add("coherence", sprintf("%s lens, mean pairwise similarity", l),
      lab(O, l, "mean_sim"), lab(N, l, "mean_sim"), MEASURE)
for (l in c("WT_heat_down", "KO_heat_up"))
  add("coherence", sprintf("%s, mean pairwise similarity", l),
      lab(O, l, "mean_sim"), lab(N, l, "mean_sim"), MEASURE)

add("splittability", "WT_heat_up, genes in the largest cluster of the best split",
    lab(O, "WT_heat_up", "largest_cluster_frac"), lab(N, "WT_heat_up", "largest_cluster_frac"),
    "fraction")
add("splittability", "WT_heat_up, best average silhouette",
    lab(O, "WT_heat_up", "best_sil"), lab(N, "WT_heat_up", "best_sil"), "silhouette")
for (f in c("null_depth_matched", "single_go_term", "null_uniform"))
  add("splittability", sprintf("%s, mean genes in the largest cluster", f),
      mean(fam(O, f, "largest_cluster_frac")), mean(fam(N, f, "largest_cluster_frac")),
      "fraction")

ls_ <- function(x, l, col) x$lens_summary[[col]][x$lens_summary$lens == l]
for (l in c("TCR_activation", "HSR_core", "Lombardi2022_HIF")) {
  add("lens_proximity", sprintf("%s, genes scorable against a depth-matched background", l),
      ls_(O, l, "n_scored"), ls_(N, l, "n_scored"), "genes")
  add("lens_proximity", sprintf("%s, genes above the depth-matched 95th percentile", l),
      ls_(O, l, "n_above_p95"), ls_(N, l, "n_above_p95"), "genes")
  add("lens_proximity", sprintf("%s, count expected by chance", l),
      ls_(O, l, "n_expected_by_chance"), ls_(N, l, "n_expected_by_chance"), "genes")
  add("lens_proximity", sprintf("%s, genes left uncounted in a saturated band", l),
      ls_(O, l, "n_excluded_saturated_band"), ls_(N, l, "n_excluded_saturated_band"), "genes")
}

add("annotation_depth", "WT_heat_up, median GO BP terms per gene",
    sv(O, "query_med_terms_per_gene"), sv(N, "query_med_terms_per_gene"), "terms",
    "engine-independent; a property of the annotation")
add("annotation_depth", "uniform null, median GO BP terms per gene",
    sv(O, "uniform_null_med_terms_per_gene"), sv(N, "uniform_null_med_terms_per_gene"),
    "terms", "engine-independent; a property of the annotation")

val <- dplyr::bind_rows(rows)

## ---------------------------------------------------------------------------
## Write
## ---------------------------------------------------------------------------
TBL_DIR <- stage_dir(STAGE, "tables")
out <- file.path(TBL_DIR, "semantic_engine_validation.csv")
write.csv(val, out, row.names = FALSE)
message("[out] ", nrow(val), " rows -> ", out)
message("[out] rows that moved: ", sum(val$changed, na.rm = TRUE),
        " of ", sum(!is.na(val$changed)))
message("20b_semantic_engine_validation: COMPUTE DONE.")
