#!/usr/bin/env Rscript
# 11_emit_universe.R — COMPUTE
## Emit the GSEA background gene universe and candidate-set frame for the STING standard
## sweep. Stage: 06_gsea. Run from project root AFTER 02_de_limma_trend.R and BEFORE the
## GSEA producers (04_gsea_set_prep.R, 05_gsea_msigdb_run.R, etc.).
##
## PURPOSE
##   Writes the two outputs every downstream GSEA arm consumes:
##     (1) 03_results/objects/gene_universe.txt  — one MGI symbol per line, sorted.
##         The contrast-invariant enrichment background: the union of all modelled symbols
##         from 02_de_results.rds (post dup-collapse, ~19,679 symbols). fgsea /
##         clusterProfiler callers pass this as their `universe` argument.
##     (2) 03_results/master/universe_frame.csv  — gene-level candidate-set frame.
##         One row per background gene, mirroring 14839's per-gene attribute semantics:
##         detected (appears in the modelled DE gene list, always TRUE here by
##         construction), expressed (AveExpr > log2(1 CPM) across >=1 contrast), plus
##         per-contrast significance and direction flags that downstream arms use to subset
##         the background or annotate enrichment results.
##
## HOW THIS DIFFERS FROM 14839's universe_frame.csv (a set-level frame)
##   14839/11_emit_universe.R writes a SET-level frame (pathway_id / entity_type /
##   genes_full_set) that feeds a set-level embedding. This compartment needs the
##   GENE-level background annotation frame, which is the natural "candidate-set frame"
##   for GSEA background annotation, one row per symbol.
##   Column-for-column mapping vs 14839:
##     gene_symbol   — present (the primary key)
##     detected      — present (always TRUE; synonym of being in the modelled universe)
##     expressed     — present (AveExpr > 0 in >=1 contrast, as a weak filter proxy)
##     in_universe   — present (always TRUE for every row in this file; explicit flag)
##     [per-contrast] <contrast>_logFC, <contrast>_t, <contrast>_padj, <contrast>_sig,
##                    <contrast>_direction — derived from 02_de_results.rds topTables.
##   The columns 14839 records from atlas/entity/set membership (entity_type,
##   genes_full_set) are set-level, and this frame is gene-level by design.
##
## ASSUMED SHAPE OF 02_de_results.rds (verify in-container):
##   A named list, names == 7 YAML contrast names
##   (WT_heat, KO_heat, Interaction, Geno_at_39, Geno_at_37, Temp_main, Geno_main).
##   Each element is a data.frame (limma topTable) with:
##     - Gene SYMBOL rownames (and a `gene_symbol` column).
##     - Columns: gene_symbol, ensembl, logFC, AveExpr, t, P.Value, adj.P.Val, B, contrast.
##   (Consistent with 02_de_limma_trend.R lines 118-127 and de_gsea_helpers.R load_de_results().)
##
## Idempotent: overwrites outputs on every run. No plotting. Byte-stable.

source("02_analysis/config/config.R")        # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS, DIR_MASTER, etc.
source("02_analysis/helpers/de_gsea_helpers.R")  # load_de_results(), gene_universe(), round_numeric_cols()

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})
options(stringsAsFactors = FALSE)

## ---------------------------------------------------------------------------
## Paths
## ---------------------------------------------------------------------------
UNIVERSE_TXT <- file.path(DIR_OBJECTS, "gene_universe.txt")
UNIVERSE_CSV <- file.path(DIR_MASTER, "universe_frame.csv")

dir.create(DIR_OBJECTS, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_MASTER,  recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------------------
## (1) Load DE results (the GSEA hub)
##     load_de_results() validates Symbol rownames and presence of the `t` column.
## ---------------------------------------------------------------------------
de <- load_de_results()   # named list of topTables; default path = DIR_OBJECTS/02_de_results.rds
contrasts <- names(de)
message(sprintf("DE results loaded: %d contrasts (%s).", length(contrasts), paste(contrasts, collapse = ", ")))

## ---------------------------------------------------------------------------
## (2) Compute the background gene universe — union of all contrast Symbol rownames.
##     gene_universe() de-duplicates and sorts. Output: character vector of MGI symbols.
## ---------------------------------------------------------------------------
universe_genes <- load_or_compute(
  "gene_universe_vec.rds",
  function() gene_universe(de),
  force = FALSE
)
message(sprintf("Background gene universe: %d symbols.", length(universe_genes)))

## ---------------------------------------------------------------------------
## (3) Write gene_universe.txt — one symbol per line (sorted, no header).
##     GSEA producers pass this to fgsea / clusterProfiler as their `universe`.
## ---------------------------------------------------------------------------
writeLines(universe_genes, UNIVERSE_TXT)
message(sprintf("Written: %s  (%d symbols, one per line)", UNIVERSE_TXT, length(universe_genes)))

## ---------------------------------------------------------------------------
## (4) Build the gene-level candidate-set frame (universe_frame.csv).
##
##     For each background gene, record:
##       gene_symbol  — the MGI symbol (primary key)
##       detected     — TRUE for every row (being in the modelled universe is the definition)
##       in_universe  — TRUE for every row (redundant but explicit; for downstream joins)
##       expressed    — AveExpr > 0 in >=1 contrast (limma AveExpr = mean log2 CPM;
##                      > 0 means mean CPM > 1, a weak but unambiguous "expressed" proxy)
##       For each contrast c:
##         <c>_logFC     — log2 fold-change
##         <c>_t         — limma-trend t-statistic (the GSEA rank metric)
##         <c>_padj      — BH-adjusted p-value (adj.P.Val)
##         <c>_sig       — logical: padj < DE_FDR & |logFC| >= DE_LOGFC
##         <c>_direction — "Up" / "Down" / "NS"
##
##     All numeric columns rounded to 9 significant figures (byte-stable; round_numeric_cols()).
##     Genes absent from a contrast's topTable (should not occur by design but handled
##     defensively with NA) do NOT break the join — they get NA in that contrast's columns.
## ---------------------------------------------------------------------------
build_universe_frame <- function(de_list, universe_syms, de_fdr, de_logfc) {
  ## Seed from the universe symbol vector
  frame <- tibble::tibble(
    gene_symbol = universe_syms,
    detected    = TRUE,
    in_universe = TRUE
  )

  ## Per-contrast columns via left-join on gene_symbol
  for (cn in names(de_list)) {
    tt <- de_list[[cn]]
    ## Ensure gene_symbol column exists (rownames are the canonical source)
    if (!"gene_symbol" %in% colnames(tt)) tt$gene_symbol <- rownames(tt)

    ## Select the columns we need (any_of() silently drops absent cols — defensive for
    ## non-standard topTable shapes). Then normalise adj.P.Val -> padj before prefixing.
    want <- intersect(c("gene_symbol", "logFC", "t", "adj.P.Val", "AveExpr"), colnames(tt))
    sub  <- tt[, want, drop = FALSE]
    ## Normalise column name: adj.P.Val -> padj (friendlier for downstream joins/prints)
    if ("adj.P.Val" %in% colnames(sub)) colnames(sub)[colnames(sub) == "adj.P.Val"] <- "padj"
    ## Prefix all non-key columns with the contrast name
    data_cols <- setdiff(colnames(sub), "gene_symbol")
    colnames(sub)[colnames(sub) %in% data_cols] <- paste0(cn, "_", data_cols)

    frame <- dplyr::left_join(frame, sub, by = "gene_symbol")
  }

  ## expressed: AveExpr > 0 in >=1 contrast (log2 CPM > 0 means mean CPM > 1).
  ## AveExpr columns are named <contrast>_AveExpr.
  aveexpr_cols <- grep("_AveExpr$", colnames(frame), value = TRUE)
  if (length(aveexpr_cols) > 0) {
    frame$expressed <- apply(
      frame[, aveexpr_cols, drop = FALSE], 1,
      function(row) any(!is.na(row) & row > 0)
    )
  } else {
    frame$expressed <- NA  # defensive: if AveExpr absent, mark NA (owner verifies in-container)
  }

  ## sig + direction columns (per-contrast)
  for (cn in names(de_list)) {
    padj_col  <- paste0(cn, "_padj")
    logfc_col <- paste0(cn, "_logFC")
    sig_col   <- paste0(cn, "_sig")
    dir_col   <- paste0(cn, "_direction")

    if (padj_col %in% colnames(frame) && logfc_col %in% colnames(frame)) {
      frame[[sig_col]] <- !is.na(frame[[padj_col]]) &
        !is.na(frame[[logfc_col]]) &
        frame[[padj_col]] < de_fdr &
        abs(frame[[logfc_col]]) >= de_logfc

      frame[[dir_col]] <- dplyr::case_when(
        is.na(frame[[logfc_col]])         ~ NA_character_,
        !frame[[sig_col]]                  ~ "NS",
        frame[[logfc_col]] > 0             ~ "Up",
        TRUE                               ~ "Down"
      )
    }
  }

  ## Drop the raw AveExpr columns (expressed is the derived summary; keep the per-contrast
  ## logFC/t/padj/sig/direction; dropping AveExpr keeps the frame compact).
  frame <- frame %>% dplyr::select(-dplyr::ends_with("_AveExpr"))

  ## Column order: gene_symbol, detected, in_universe, expressed, then contrast blocks
  meta_cols    <- c("gene_symbol", "detected", "in_universe", "expressed")
  contrast_cols <- setdiff(colnames(frame), meta_cols)
  frame <- frame %>% dplyr::select(dplyr::all_of(meta_cols), dplyr::all_of(contrast_cols))

  frame
}

universe_frame <- load_or_compute(
  "universe_frame_obj.rds",
  function() build_universe_frame(de, universe_genes, DE_FDR, DE_LOGFC),
  force = FALSE
)

## Validate
stopifnot("universe_frame must have a gene_symbol column" = "gene_symbol" %in% colnames(universe_frame))
stopifnot("universe_frame must not be empty"              = nrow(universe_frame) > 0)
stopifnot("universe_frame gene_symbol must be unique"     = !anyDuplicated(universe_frame$gene_symbol))
stopifnot("universe_frame nrow must match universe"       = nrow(universe_frame) == length(universe_genes))

## Byte-stable write (round all double columns to 9 sig figs)
readr::write_csv(round_numeric_cols(universe_frame, sig = 9L), UNIVERSE_CSV)
message(sprintf("Written: %s  (%d rows × %d cols)", UNIVERSE_CSV, nrow(universe_frame), ncol(universe_frame)))

## ---------------------------------------------------------------------------
## Summary report
## ---------------------------------------------------------------------------
sig_counts <- vapply(contrasts, function(cn) {
  col <- paste0(cn, "_sig")
  if (col %in% colnames(universe_frame)) sum(universe_frame[[col]], na.rm = TRUE) else NA_integer_
}, integer(1))

message("\n--- universe_frame.csv column summary ---")
message(sprintf("  Total background genes : %d", nrow(universe_frame)))
if ("expressed" %in% colnames(universe_frame) && !all(is.na(universe_frame$expressed))) {
  message(sprintf("  Expressed (AveExpr>0)  : %d", sum(universe_frame$expressed, na.rm = TRUE)))
}
message("  Significant DE per contrast (padj<", DE_FDR, ", |logFC|>=", DE_LOGFC, "):")
for (cn in contrasts) {
  n <- sig_counts[cn]
  message(sprintf("    %-15s : %s", cn, if (!is.na(n)) as.character(n) else "n/a"))
}

message("\n11_emit_universe complete.")
message("  gene_universe.txt  -> ", UNIVERSE_TXT)
message("  universe_frame.csv -> ", UNIVERSE_CSV)
