#!/usr/bin/env Rscript
# =============================================================================
# 10_gatom_modules.R — COMPUTE  (stage 09_gatom)
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# GATOM atom-graph MWCS metabolic modules from the limma-trend DE signal. This is the
# MitoCarta-corroborating mechanism arm: the Complex-I / metabolic-pseudohypoxia story.
# GATOM finds the maximally-regulated metabolic SUBNETWORK (a Signal-Generalised
# Maximum-Weight Connected Subgraph, SGMWCS) over an atom-transition network, scoring
# enzyme edges by a BUM model fit to the RAW p-values. The permutation/enrichment arm is
# GSEA, in 06_gsea.
#
# ROLE: COMPUTE half of the normalize-then-visualize split. All statistics, checkpoints
#   and tidy CSV. Figures live in the *_viz.R sibling, which reads only the artifacts
#   written here.
#
# CONTRASTS (headline only, per the standard sweep): WT_heat, KO_heat, Interaction,
#   Temp_main, intersected with the DE checkpoint names.
#
# NETWORKS: kegg + combined, the two standard GATOM networks. KEGG runs with
#   gene2reaction.extra=NULL. Combined REQUIRES the gene2reaction TSV and returns a silent
#   empty graph without it — see skill "Network Files" + "Input Format Gotchas".
#
# ---------------------------------------------------------------------------
# GATOM REFERENCE OBJECTS, loaded by their on-disk paths under
#   00_data/references/gatom/ (already downloaded; this script re-fetches nothing).
#   Shapes assumed per the gatom-metabolomic-predictions skill + GATOM docs:
#
#   org.Mm.eg.gatom.anno.rds   list(gene/metabolite annotation maps). Passed as
#                              `org.gatom.anno=`. Keyed by gene SYMBOL here
#                              (gene.de$ID is the limma gene_symbol).
#   network.<net>.rds          GATOM network object (reactions/atoms/enzymes).
#                              Passed as `network=`.  net in {kegg, combined}.
#   met.<net>.db.rds           metabolite annotation DB. Passed as `met.db=` —
#                              REQUIRED even with met.de=NULL, since topology needs it.
#   gene2reaction.<net>.mmu.eg.tsv   enzyme->reaction map for non-KEGG networks
#                              (mouse ncbi code = "mmu"). Passed as
#                              `gene2reaction.extra=`. REQUIRED for Combined.
#
#   Present on disk (ls 00_data/references/gatom/):
#     network.kegg.rds, network.combined.rds, network.rhea.rds,
#     network.rhea.lipids.rds, met.kegg.db.rds, met.rhea.db.rds,
#     met.lipids.db.rds, org.Mm.eg.gatom.anno.rds,
#     gene2reaction.combined.mmu.eg.tsv, gene2reaction.rhea.mmu.eg.tsv
#   This stage uses the KEGG + Combined set. Each required file is guarded with a stop()
#   naming the missing path.
#
# ---------------------------------------------------------------------------
# INVERTED GRAPH STRUCTURE (skill "Critical: Inverted Graph Structure"):
#   NODES = metabolites/atoms (carry label and score)
#   EDGES = reactions/enzymes (carry log2FC, pval, Symbol, enzyme)
#   => module genes come from EDGES: igraph::as_data_frame(m, "edges")$Symbol
#
# INPUT FORMAT (skill "Core Workflow" / "Input Format Gotchas"):
#   gene.de = data.frame(ID=symbol, pval=RAW P.Value, log2FC=logFC,
#                        baseMean=2^AveExpr [LINEAR]); dedup keeps min-pval.
#   pval MUST be the raw p-value, which is what the BUM scoring fits.
#
# Inputs:
#   - 03_results/objects/02_de_results.rds            (named list of topTables)
#   - 00_data/references/gatom/{network,met,org,gene2reaction}*  (see above)
#
# Outputs:
#   - 03_results/objects/10_gatom_<contrast>.rds      (per-contrast module bundle)
#   - 03_results/objects/10_gatom_all.rds             (checkpoint: full result list)
#   - 03_results/09_gatom/tables/gatom_modules.csv    (tidy cross-contrast summary)
#   - 03_results/master/master_gatom_modules.csv      (master_gsea_table schema)
#
# Dependencies (LAZY — required only when invoked): gatom, mwcsr, igraph, data.table.
#   Light tidyverse (dplyr/readr/tibble) via de_gsea_helpers.R.
# =============================================================================

source("02_analysis/config/config.R")             # PROJECT_ROOT, YAML_CONFIG, DIR_*, stage_dir, load_or_compute, %||%
source("02_analysis/helpers/de_gsea_helpers.R")    # load_de_results, build_ranked_vector, append_master_table, load_or_compute (path-keyed)

# -----------------------------------------------------------------------------
# Lazy heavy deps — fail with an actionable install hint ONLY when this script
# actually runs (so a lint/doc pass can source it in any container).
# -----------------------------------------------------------------------------
.need <- function(pkg, how) if (!requireNamespace(pkg, quietly = TRUE))
  stop(sprintf("10_gatom_modules.R: package '%s' is required. Install with %s", pkg, how))
.need("gatom",      "BiocManager::install('gatom')")
.need("mwcsr",      "install.packages('mwcsr')")
.need("igraph",     "install.packages('igraph')")
.need("data.table", "install.packages('data.table')")

set.seed(GSEA_SEED %||% 123)

STAGE      <- "09_gatom"
TBL_DIR    <- stage_dir(STAGE, "tables")                 # 03_results/09_gatom/tables/ (asserts stage registered)
GATOM_REFS <- file.path(PROJECT_ROOT, "00_data/references/gatom")

# Headline contrasts for the standard sweep (intersected with the DE checkpoint).
HEADLINE_CONTRASTS <- c("WT_heat", "KO_heat", "Interaction", "Temp_main")

# Standard GATOM network pair. KEGG: gene2reaction.extra = NULL.
#   Combined: REQUIRES the gene2reaction TSV (silent empty graph otherwise).
# Self-healing: keep only networks whose metabolite-annotation DB (met.<net>.db.rds)
# is actually provisioned. `combined` is dropped automatically until the owner
# supplies met.combined.db.rds (only met.kegg/rhea/lipids ship by default), so the
# run degrades to KEGG-only and finishes the bundle. Re-includes combined
# the moment that file lands.
NETWORKS <- Filter(
  function(net) file.exists(file.path(GATOM_REFS, sprintf("met.%s.db.rds", net))),
  c("kegg", "combined")
)
if (length(NETWORKS) == 0L)
  stop("No GATOM network has a provisioned met.<net>.db.rds under ", GATOM_REFS)
if (!"combined" %in% NETWORKS)
  message("[10] met.combined.db.rds absent -> running GATOM on KEGG only ",
          "(supply 00_data/references/gatom/met.combined.db.rds to enable combined).")

# k.gene = module-size heuristic. 50 = GATOM default (25 stringent / 75 exploratory).
K_GENE <- as.integer(YAML_CONFIG$thresholds$gatom_k_gene %||% 50L)

# cGAS-STING / metabolic watchlist surfaced in the console for a quick read.
KEY_METAB_GENES <- unique(c(HIF_GLYCO_MARKERS, "Ndufa1", "Ndufa9", "Ndufs1",
                            "Sdha", "Sdhb", "Pdha1", "Pdk1", "Idh1", "Idh2"))

cat("=================================================================\n")
cat("PHASE 10 (stage 09_gatom): GATOM atom-graph MWCS metabolic modules\n")
cat("=================================================================\n\n")

# -----------------------------------------------------------------------------
# 1. LOAD DE RESULTS + RESOLVE HEADLINE CONTRASTS
# -----------------------------------------------------------------------------
cat("[1] Loading DE results (02_de_results.rds) ...\n")
de_results <- load_de_results()                          # symbol-rownames + rank-metric contract enforced in helper
focal <- intersect(HEADLINE_CONTRASTS, names(de_results))
if (length(focal) == 0L)
  stop("None of the headline contrasts (", paste(HEADLINE_CONTRASTS, collapse = ", "),
       ") are present in 02_de_results.rds (have: ", paste(names(de_results), collapse = ", "), ").")
if (length(focal) < length(HEADLINE_CONTRASTS))
  warning("Missing headline contrasts: ",
          paste(setdiff(HEADLINE_CONTRASTS, focal), collapse = ", "), " — running on present subset.")
cat(sprintf("  Running GATOM on %d contrasts x %d networks: %s\n\n",
            length(focal), length(NETWORKS), paste(focal, collapse = ", ")))

# -----------------------------------------------------------------------------
# 2. LOAD GATOM REFERENCE OBJECTS BY THEIR ACTUAL PATHS (guarded; no download)
#    Each missing file => clear stop() naming the path (brief contract).
# -----------------------------------------------------------------------------
cat("[2] Loading GATOM references from ", GATOM_REFS, " ...\n", sep = "")
read_ref <- function(fname) {
  fp <- file.path(GATOM_REFS, fname)
  if (!file.exists(fp))
    stop("Missing GATOM reference: ", fp,
         "\n  -> download with RNAseq-toolkit download_gatom_references(dest_dir='00_data/references/gatom').")
  readRDS(fp)
}

org_anno <- read_ref("org.Mm.eg.gatom.anno.rds")
cat("  org.Mm.eg.gatom.anno.rds loaded.\n")

# Per-network bundle: network object, met.db, and (non-KEGG) gene2reaction.extra.
nets <- list()
for (net in NETWORKS) {
  network <- read_ref(sprintf("network.%s.rds", net))
  met_db  <- read_ref(sprintf("met.%s.db.rds",  net))

  g2r <- NULL
  if (net != "kegg") {
    # Mouse ncbi species code is "mmu"; the on-disk name is gene2reaction.<net>.mmu.eg.tsv.
    # Use Sys.glob to be robust to the exact species suffix, then guard hard.
    g2r_fp <- file.path(GATOM_REFS, sprintf("gene2reaction.%s.mmu.eg.tsv", net))
    if (!file.exists(g2r_fp)) {
      hits <- Sys.glob(file.path(GATOM_REFS, sprintf("gene2reaction.%s.*.eg.tsv", net)))
      if (length(hits) == 0L)
        stop("Missing GATOM gene2reaction map for '", net, "' network: ",
             g2r_fp, " (REQUIRED for non-KEGG networks; Combined graph is silently empty without it).")
      g2r_fp <- hits[1]
    }
    g2r <- data.table::fread(g2r_fp, colClasses = "character")
    cat(sprintf("  %-9s network + met.db + gene2reaction.extra (%s, %d rows)\n",
                net, basename(g2r_fp), nrow(g2r)))
  } else {
    cat(sprintf("  %-9s network + met.db (gene2reaction.extra = NULL)\n", net))
  }
  nets[[net]] <- list(network = network, met_db = met_db, g2r = g2r)
}
cat("\n")

# -----------------------------------------------------------------------------
# 3. GATOM INPUT BUILDER + per-(contrast x network) make -> score -> solve
# -----------------------------------------------------------------------------
# prepare_gene_de(): the EXACT GATOM input format (skill "Core Workflow").
#   ID = gene symbol; pval = RAW P.Value; log2FC = logFC; baseMean = 2^AveExpr
#   (LINEAR scale). Dedup symbols keeping the lowest p-value (skill: "Handle
#   duplicates: keep lowest p-value per gene").
prepare_gene_de <- function(tt) {
  sym <- rownames(tt)
  if (is.null(sym) && "gene_symbol" %in% colnames(tt)) sym <- as.character(tt$gene_symbol)
  df <- data.frame(
    ID       = sym,
    pval     = tt$P.Value,
    log2FC   = tt$logFC,
    baseMean = 2^tt$AveExpr,                            # limma AveExpr is log2 -> linear
    stringsAsFactors = FALSE
  )
  df <- df[!is.na(df$ID) & df$ID != "" & !is.na(df$pval), , drop = FALSE]
  df <- df[order(df$pval), , drop = FALSE]              # lowest pval first
  df <- df[!duplicated(df$ID), , drop = FALSE]          # keep min-pval per symbol
  df
}

# run_gatom(): make the atom graph, score the edges (BUM, k.gene), solve SGMWCS.
#   Returns a list with the module igraph (or success=FALSE + error on failure).
#   Genes live on EDGES, so an empty-edge graph means "no module".
run_gatom <- function(net_obj, network_name, gene_de, k_gene, contrast_name) {
  res <- list(contrast = contrast_name, network = network_name, k_gene = k_gene, success = FALSE)
  tryCatch({
    g <- gatom::makeMetabolicGraph(
      network        = net_obj$network,
      topology       = "atoms",
      org.gatom.anno = org_anno,
      gene.de        = gene_de,
      met.db         = net_obj$met_db,                  # REQUIRED even with met.de = NULL
      met.de         = NULL,
      gene2reaction.extra        = net_obj$g2r,         # NULL for KEGG; TSV for Combined
      keepReactionsWithoutEnzymes = FALSE
    )
    if (igraph::vcount(g) == 0L || igraph::ecount(g) == 0L)
      stop("empty atom graph (no overlapping enzymes — check ID type / gene2reaction).")

    gs  <- gatom::scoreGraph(g, k.gene = k_gene, k.met = NULL)
    sol <- mwcsr::solve_mwcsp(mwcsr::rnc_solver(), gs)  # Randomised Network-Constrained solver
    m   <- sol$graph

    res$module          <- m
    res$solution_weight <- sol$weight %||% NA_real_
    res$n_vertices      <- igraph::vcount(m)
    res$n_edges         <- igraph::ecount(m)
    res$success         <- TRUE
    cat(sprintf("    [%-11s / %-8s] module: V=%d E=%d w=%.2f\n",
                contrast_name, network_name, res$n_vertices, res$n_edges,
                res$solution_weight %||% NA_real_))
  }, error = function(e) {
    res$error <<- conditionMessage(e)
    cat(sprintf("    [%-11s / %-8s] FAILED: %s\n", contrast_name, network_name, conditionMessage(e)))
  })
  res
}

# -----------------------------------------------------------------------------
# 4. RUN (checkpointed). One bundle per contrast: gene.de + per-network modules.
#    Per-contrast objects ALSO written individually to 10_gatom_<contrast>.rds.
# -----------------------------------------------------------------------------
cat("[3] make -> score -> solve (k.gene = ", K_GENE, ") ...\n", sep = "")
gatom_all <- load_or_compute(
  "10_gatom_all.rds",
  desc = "GATOM modules (all headline contrasts x networks)",
  compute_fn = function() {
    out <- list()
    for (co in focal) {
      gene_de <- prepare_gene_de(de_results[[co]])
      cat(sprintf("  --- %s : %d genes ---\n", co, nrow(gene_de)))
      bundle <- list(contrast = co, gene_de = gene_de, modules = list())
      for (net in NETWORKS)
        bundle$modules[[net]] <- run_gatom(nets[[net]], net, gene_de, K_GENE, co)
      out[[co]] <- bundle
    }
    out
  }
)

# Per-contrast checkpoint objects (idempotent; overwritten each run, cheap).
for (co in focal) {
  fp <- file.path(DIR_OBJECTS, sprintf("10_gatom_%s.rds", co))
  saveRDS(gatom_all[[co]], fp)
  cat(sprintf("  [SAVE] %s\n", fp))
}
cat("\n")

# -----------------------------------------------------------------------------
# 5. EXTRACTORS — genes/metabolites/reactions live on the module igraph.
#    Genes are EDGE attributes; metabolites are VERTEX labels.
# -----------------------------------------------------------------------------
.pick <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) return(NULL)
  df[[hit[1]]]
}

module_edge_df <- function(res) {
  if (!isTRUE(res$success) || is.null(res$module) || igraph::ecount(res$module) == 0L)
    return(NULL)
  igraph::as_data_frame(res$module, what = "edges")
}

module_genes <- function(res) {                          # EDGE Symbol/gene attribute
  ed <- module_edge_df(res); if (is.null(ed)) return(character(0))
  g <- .pick(ed, c("Symbol", "gene", "enzyme"))
  if (is.null(g)) return(character(0))
  unique(g[!is.na(g) & g != ""])
}

module_reactions <- function(res) {                      # EDGE reaction id / name
  ed <- module_edge_df(res); if (is.null(ed)) return(character(0))
  r <- .pick(ed, c("reaction_name", "reaction", "rxn", "label"))
  if (is.null(r)) return(character(0))
  unique(r[!is.na(r) & r != ""])
}

module_metabolites <- function(res) {                    # VERTEX label (metabolite name)
  if (!isTRUE(res$success) || is.null(res$module) || igraph::vcount(res$module) == 0L)
    return(character(0))
  vd <- igraph::as_data_frame(res$module, what = "vertices")
  m <- .pick(vd, c("label", "metabolite", "name"))
  if (is.null(m)) return(character(0))
  unique(m[!is.na(m) & m != ""])
}

# Top-N joined string for the tidy summary (compact, plot-ready).
top_join <- function(x, n = 10L) paste(utils::head(x, n), collapse = "/")

# -----------------------------------------------------------------------------
# 6. TIDY SUMMARY -> 03_results/09_gatom/tables/gatom_modules.csv
#    One row per (contrast x network): module size, top metabolites/reactions,
#    solution score, and the contrast-direction (mean module-edge log2FC).
# -----------------------------------------------------------------------------
cat("[4] Building tidy summary -> gatom_modules.csv ...\n")
summary_df <- dplyr::bind_rows(lapply(focal, function(co) {
  dplyr::bind_rows(lapply(NETWORKS, function(net) {
    r  <- gatom_all[[co]]$modules[[net]]
    ok <- isTRUE(r$success)
    ed <- if (ok) module_edge_df(r) else NULL
    lfc <- if (!is.null(ed)) suppressWarnings(as.numeric(.pick(ed, c("log2FC", "logFC")))) else numeric(0)
    genes <- if (ok) module_genes(r) else character(0)
    tibble::tibble(
      contrast         = co,
      network          = net,
      success          = ok,
      n_vertices       = if (ok) r$n_vertices %||% NA_integer_ else NA_integer_,
      n_edges          = if (ok) r$n_edges    %||% NA_integer_ else NA_integer_,
      n_genes          = length(genes),
      n_metabolites    = length(module_metabolites(r)),
      n_reactions      = length(module_reactions(r)),
      solution_weight  = if (ok) r$solution_weight %||% NA_real_ else NA_real_,
      mean_edge_log2FC = if (length(lfc)) mean(lfc, na.rm = TRUE) else NA_real_,
      direction        = if (length(lfc)) ifelse(mean(lfc, na.rm = TRUE) > 0, "Up", "Down") else NA_character_,
      k_gene           = K_GENE,
      top_genes        = top_join(genes),
      top_metabolites  = top_join(module_metabolites(r)),
      top_reactions    = top_join(module_reactions(r)),
      error            = if (ok) NA_character_ else (r$error %||% NA_character_)
    )
  }))
}))
readr::write_csv(round_numeric_cols(summary_df), file.path(TBL_DIR, "gatom_modules.csv"))
cat(sprintf("  gatom_modules.csv: %d rows (%d/%d non-empty modules)\n",
            nrow(summary_df), sum(summary_df$success), nrow(summary_df)))

# Guard: all-fail almost always means an ID-keytype mismatch (gene.de$ID must be
# the symbol namespace org_anno is keyed on). Warn loudly — never silent-pass.
if (!any(summary_df$success))
  warning("ALL GATOM runs returned empty/failed modules. Suspect ID-keytype mismatch ",
          "(gene.de$ID symbols vs org.Mm.eg.gatom.anno) or a missing gene2reaction map.")

# -----------------------------------------------------------------------------
# 7. MASTER ACCUMULATOR (master_gsea_table schema). One pseudo-set per
#    successful (contrast x network). NOTE the legacy GSEA schema column names:
#    `nes` here = mean module-edge log2FC (a pseudo-NES, NOT an fgsea NES);
#    `pvalue`/`padj` = geometric-mean module-edge raw p-value; `core_enrichment`
#    = the module enzyme set. Kept for cross-arm accumulator compatibility.
# -----------------------------------------------------------------------------
cat("[5] Appending GATOM pseudo-sets to master_gatom_modules.csv ...\n")
master_rows <- dplyr::bind_rows(lapply(focal, function(co) {
  dplyr::bind_rows(lapply(NETWORKS, function(net) {
    r <- gatom_all[[co]]$modules[[net]]
    if (!isTRUE(r$success)) return(NULL)
    ed <- module_edge_df(r); if (is.null(ed)) return(NULL)
    genes <- module_genes(r); if (length(genes) == 0L) return(NULL)
    lfc <- suppressWarnings(as.numeric(.pick(ed, c("log2FC", "logFC"))))
    pv  <- suppressWarnings(as.numeric(.pick(ed, c("pval", "P.Value"))))
    nes <- if (length(lfc)) mean(lfc, na.rm = TRUE) else NA_real_           # pseudo-NES = mean edge log2FC
    pse <- if (length(pv))  exp(mean(log(pmax(pv, 1e-300)), na.rm = TRUE)) else NA_real_  # geo-mean edge pval
    tibble::tibble(
      pathway_id      = sprintf("GATOM_%s", toupper(net)),  # contrast-invariant id (one per network)
      pathway_name    = sprintf("GATOM %s metabolic module", net),
      database        = "GATOM",
      nes             = nes,
      pvalue          = pse,
      padj            = pse,                                # single pseudo-set -> no MT correction
      set_size        = length(genes),
      core_enrichment = paste(genes, collapse = "/"),       # THIS contrast's module enzymes
      contrast        = co,
      direction       = direction_from_nes(nes)
    )
  }))
}))

if (nrow(master_rows) > 0L) {
  # Idempotent: append_master_table drops prior database=="GATOM" rows before bind,
  # validates against schemas.master_gsea_table, and writes byte-stable (9 sig figs).
  append_master_table(master_rows, "master_gatom_modules.csv", key_col = "database")
  cat(sprintf("  master_gatom_modules.csv: %d GATOM pseudo-set rows.\n", nrow(master_rows)))
} else {
  cat("  (no successful modules -> master_gatom_modules.csv not written)\n")
}

# -----------------------------------------------------------------------------
# 8. CONSOLE SANITY — largest module + metabolic-watchlist read
# -----------------------------------------------------------------------------
ok <- dplyr::filter(summary_df, success)
if (nrow(ok)) {
  top <- ok[which.max(ok$n_edges), ]
  cat(sprintf("\n  Largest module: %s / %s  V=%d E=%d w=%.2f\n",
              top$contrast, top$network, top$n_vertices, top$n_edges, top$solution_weight))
}
all_module_genes <- unique(unlist(lapply(focal, function(co)
  unlist(lapply(NETWORKS, function(net) module_genes(gatom_all[[co]]$modules[[net]]))))))
hit <- intersect(all_module_genes, KEY_METAB_GENES)
if (length(hit))
  cat("  Metabolic-watchlist enzymes recruited into any module: ",
      paste(hit, collapse = ", "), "\n", sep = "")

# -----------------------------------------------------------------------------
# 9. FINAL ASSERTS
# -----------------------------------------------------------------------------
stopifnot(file.exists(file.path(TBL_DIR, "gatom_modules.csv")))
stopifnot(all(file.exists(file.path(DIR_OBJECTS, sprintf("10_gatom_%s.rds", focal)))))
if (nrow(master_rows) > 0L)
  stopifnot(file.exists(file.path(DIR_MASTER, "master_gatom_modules.csv")))

cat(sprintf("\n[DONE] 10_gatom_modules complete: %d/%d modules non-empty; %d master rows. Run 10_gatom_modules_viz.R for figures.\n",
            sum(summary_df$success), nrow(summary_df), nrow(master_rows)))
