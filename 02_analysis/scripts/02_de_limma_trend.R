#!/usr/bin/env Rscript
# =============================================================================
# 02_de_limma_trend.R  --  Phase 2 COMPUTE: limma-trend 2x2 DE (statistics only)
# =============================================================================
# Phase:    2 (stage 03_de)
# Role:     COMPUTE half of the "normalize-then-visualize" split. Does ALL
#           statistics and emits checkpoints + plot-ready tidy tables. Contains
#           NO ggplot()/ggsave(); figures live in 02_de_limma_trend_viz.R.
#
# Inputs:   03_results/objects/01_eda.rds  (cpm_mat, logcpm_mat = log2(CPM+0.5),
#                                           metadata aligned 021..040, collapse_record)
#           02_analysis/config/analysis_config.yaml  (design.contrasts, 7 contrasts)
#           00_data/processed/GSE329522_normalized_counts_CPM_iTreg.csv (ensembl map)
# Contrasts (7, built generically from config:design.contrasts via makeContrasts):
#   WT_heat    = WT_39 - WT_37          heat response in the cGAS-competent genotype
#   KO_heat    = cGASKO_39 - cGASKO_37 heat response with cGAS removed
#   Interaction= (WT_39-cGASKO_39) - (WT_37-cGASKO_37)  cGAS-dependence test (1 df)
#   Geno_at_39 = WT_39 - cGASKO_39     genotype effect under heat
#   Geno_at_37 = WT_37 - cGASKO_37     genotype effect at baseline
#   Temp_main  = ½(WT_39+cGASKO_39) - ½(WT_37+cGASKO_37)  pooled average 39-vs-37 °C response
#   Geno_main  = ½(WT_37+WT_39) - ½(cGASKO_37+cGASKO_39)  pooled average genotype effect
#
# Outputs:  03_results/objects/02_de_results.rds
#             named list of 7 topTables; gene-symbol rownames;
#             cols: gene_symbol / ensembl / logFC / AveExpr / t /
#                   P.Value / adj.P.Val / B / contrast
#           03_results/master/master_de_genes.csv
#             7 contrasts stacked;
#             cols: gene_symbol / ensembl / logFC / t / P.Value /
#                   adj.P.Val / contrast / significant / direction
#           03_results/03_de/tables/marker_cgas_dependence.csv  (marker readout)
#           03_results/03_de/tables/fig2_marker_means.csv       (plot-ready, tidy)
#
# METHOD NOTE (CRITICAL):
#   The deposit is CPM (pre-normalized), NOT raw counts -> voom CANNOT be run
#   faithfully (no mean-variance trend from library sizes). We therefore use
#   limma-trend on log2(CPM+0.5):
#       design <- model.matrix(~0 + group)
#       fit    <- lmFit(logcpm_mat, design)
#       fit2   <- eBayes(contrasts.fit(fit, cmat), trend = TRUE, robust = TRUE)
#   trend=TRUE  -> intensity-dependent prior variance (replaces voom weights)
#   robust=TRUE -> robust empirical-Bayes hyperparameter estimation.
#   This is the recommended method for pre-normalized log-expression matrices
#   (Law et al. 2014; limma User's Guide sec. on "limma-trend").
#
# Dependencies: limma, dplyr, tidyr, yaml
# =============================================================================

source("02_analysis/config/config.R")
load_packages(c("tidyr"))
suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(tidyr)
})

set.seed(123)

# -----------------------------------------------------------------------------
# 1. Load Phase 1 checkpoint
# -----------------------------------------------------------------------------
eda <- readRDS(file.path(DIR_OBJECTS, "01_eda.rds"))
logcpm_mat <- eda$logcpm_mat
metadata   <- eda$metadata

stopifnot(identical(colnames(logcpm_mat), as.character(metadata$library_id)))

# Factor with reference levels WT, 37 already encoded in group; rebuild group
# factor with the YAML-declared level order so design columns are predictable.
group_levels <- unlist(YAML_CONFIG$design$groups)            # WT_37, cGASKO_37, WT_39, cGASKO_39
metadata$group <- factor(as.character(metadata$group), levels = group_levels)
stopifnot(!anyNA(metadata$group))
cat("Group sizes:\n"); print(table(metadata$group))

# -----------------------------------------------------------------------------
# 2. Ensembl-id map (gene_symbol -> winner Ensembl id after Phase-1 collapse)
# -----------------------------------------------------------------------------
# Phase 1 collapsed duplicate gene_name to the max-mean-CPM row. The CSV stores
# every (gene_id, gene_name); for duplicates we must keep the winner id recorded
# in collapse_record, otherwise the first occurrence.
.raw <- read.csv(FILE_CPM, check.names = FALSE, stringsAsFactors = FALSE)
ens_map <- setNames(.raw$gene_id, .raw$gene_name)            # first occurrence wins
# Override duplicates with the collapse winner id (max mean CPM).
cr <- eda$collapse_record
if (!is.null(cr) && nrow(cr) > 0) {
  ens_map[cr$gene_name] <- cr$winner_id
}
ensembl_for <- function(symbols) unname(ens_map[symbols])    # NA if absent

# -----------------------------------------------------------------------------
# 3. limma-trend fit
# -----------------------------------------------------------------------------
design <- model.matrix(~0 + group, data = metadata)
colnames(design) <- levels(metadata$group)                   # name columns by group level
cat("\nDesign columns:", paste(colnames(design), collapse = ", "), "\n")

fit <- lmFit(logcpm_mat, design)

# -----------------------------------------------------------------------------
# 4. Build contrast matrix from YAML (simple = num-denom; formula = verbatim)
# -----------------------------------------------------------------------------
contrast_defs <- YAML_CONFIG$design$contrasts
contrast_exprs <- vapply(contrast_defs, function(cd) {
  if (!is.null(cd$formula)) {
    cd$formula
  } else {
    sprintf("%s - %s", cd$numerator, cd$denominator)
  }
}, character(1))
names(contrast_exprs) <- vapply(contrast_defs, function(cd) cd$name, character(1))
cat("\nContrasts (", length(contrast_exprs), "):\n", sep = "")
for (nm in names(contrast_exprs)) cat(sprintf("  %-12s = %s\n", nm, contrast_exprs[[nm]]))

cmat <- makeContrasts(contrasts = unname(contrast_exprs), levels = design)
colnames(cmat) <- names(contrast_exprs)                      # name by YAML contrast name

fit2 <- contrasts.fit(fit, cmat)
fit2 <- eBayes(fit2, trend = TRUE, robust = TRUE)            # <-- trend+robust, NOT voom

# robust=TRUE => df.prior is a per-gene vector (down-weighted for outlier genes).
cat("\n[eBayes] df.prior: median =", round(median(fit2$df.prior), 3),
    "| range =", paste(round(range(fit2$df.prior), 3), collapse = " - "),
    "(per-gene; robust=TRUE)\n")
cat("[eBayes] s2.prior: median =", signif(median(fit2$s2.prior), 4), "\n")
cat("[eBayes] df.residual (first 6 genes):",
    paste(head(fit2$df.residual, 6), collapse = ", "),
    "(= n-#groups = 20-4 = 16 baseline)\n")
cat("[eBayes] df.total range:", paste(round(range(fit2$df.total), 2), collapse = " - "), "\n")

# -----------------------------------------------------------------------------
# 5. Per-contrast topTables -> named list keyed exactly by YAML contrast names
# -----------------------------------------------------------------------------
de_results <- lapply(colnames(cmat), function(cn) {
  tt <- topTable(fit2, coef = cn, number = Inf, sort.by = "none")
  tt$gene_symbol <- rownames(tt)
  tt$ensembl     <- ensembl_for(rownames(tt))
  tt$contrast    <- cn
  # keep t intact (decoupleR needs it) -- topTable already carries `t`
  tt[, c("gene_symbol", "ensembl", "logFC", "AveExpr", "t",
         "P.Value", "adj.P.Val", "B", "contrast")]
})
names(de_results) <- colnames(cmat)

saveRDS(de_results, file.path(DIR_OBJECTS, "02_de_results.rds"))
cat("\n[SAVE]", file.path(DIR_OBJECTS, "02_de_results.rds"),
    "-- ", length(de_results), "contrasts\n")

# -----------------------------------------------------------------------------
# 6. master_de_genes.csv  (YAML schema + ensembl + t)
# -----------------------------------------------------------------------------
master_de <- bind_rows(lapply(de_results, function(tt) {
  data.frame(
    gene_symbol = tt$gene_symbol,
    ensembl     = tt$ensembl,
    logFC       = tt$logFC,
    t           = tt$t,
    P.Value     = tt$P.Value,
    adj.P.Val   = tt$adj.P.Val,
    contrast    = tt$contrast,
    significant = tt$adj.P.Val < DE_FDR,
    direction   = ifelse(tt$logFC >= 0, "Up", "Down"),
    stringsAsFactors = FALSE
  )
}))
req_cols <- YAML_CONFIG$schemas$master_de_table$required_columns
stopifnot(all(req_cols %in% colnames(master_de)))
stopifnot(all(c("ensembl", "t") %in% colnames(master_de)))
dir.create(DIR_MASTER, recursive = TRUE, showWarnings = FALSE)
master_path <- file.path(DIR_MASTER, "master_de_genes.csv")
write.csv(master_de, master_path, row.names = FALSE)
cat("[SAVE]", master_path, "--", nrow(master_de), "rows,",
    length(unique(master_de$contrast)), "contrasts\n")

# -----------------------------------------------------------------------------
# 7. Sanity checks on contrast signs (independent group-mean recomputation)
# -----------------------------------------------------------------------------
grp_mean <- function(gene) {
  v <- logcpm_mat[gene, ]
  tapply(v, metadata$group, mean)
}
sanity <- function(gene, label) {
  gm <- grp_mean(gene)
  wt_heat_manual <- gm["WT_39"]   - gm["WT_37"]
  ko_heat_manual <- gm["cGASKO_39"] - gm["cGASKO_37"]
  inter_manual   <- wt_heat_manual - ko_heat_manual
  cat(sprintf("  %-8s (%s): meanWT37=%.2f WT39=%.2f KO37=%.2f KO39=%.2f | WT_heat=%.2f KO_heat=%.2f Inter=%.2f\n",
              gene, label, gm["WT_37"], gm["WT_39"], gm["cGASKO_37"], gm["cGASKO_39"],
              wt_heat_manual, ko_heat_manual, inter_manual))
  invisible(c(WT_heat = wt_heat_manual, KO_heat = ko_heat_manual, Interaction = inter_manual))
}
cat("\n--- Contrast-sign sanity (group means recomputed from logcpm_mat) ---\n")
if ("Hspa1b" %in% rownames(logcpm_mat)) sanity("Hspa1b", "heat-shock; WT_heat should be +")
for (g in intersect(c("Irf7", "Isg15", "Ifit1"), rownames(logcpm_mat))) sanity(g, "ISG")

# Cross-check manual vs limma for WT_heat on Hspa1b
if ("Hspa1b" %in% rownames(logcpm_mat)) {
  limma_wt <- de_results$WT_heat$logFC[de_results$WT_heat$gene_symbol == "Hspa1b"]
  man_wt   <- grp_mean("Hspa1b")["WT_39"] - grp_mean("Hspa1b")["WT_37"]
  cat(sprintf("  CHECK Hspa1b WT_heat: limma=%.3f manual=%.3f (match=%s)\n",
              limma_wt, man_wt, isTRUE(all.equal(unname(limma_wt), unname(man_wt), tolerance = 1e-6))))
}

# -----------------------------------------------------------------------------
# 8. KEY RESULT table: marker genes x {WT_heat, KO_heat, Interaction logFC + adjP}
# -----------------------------------------------------------------------------
get_stat <- function(contrast, gene, col) {
  tt <- de_results[[contrast]]
  v <- tt[[col]][tt$gene_symbol == gene]
  if (length(v) == 0) NA_real_ else v[1]
}
marker_genes <- c("Ifit1", "Isg15", "Irf7",                       # ISG arm
                  "Slc2a1", "Vegfa", "Egln3", "Bnip3",            # HIF-specific
                  "Pgk1", "Ldha")                                 # shared-glycolytic
marker_arm <- c(rep("ISG", 3), rep("HIF-specific", 4), rep("shared-glyco", 2))
key_tbl <- data.frame(
  gene        = marker_genes,
  arm         = marker_arm,
  WT_heat_lfc = sapply(marker_genes, get_stat, contrast = "WT_heat",     col = "logFC"),
  KO_heat_lfc = sapply(marker_genes, get_stat, contrast = "KO_heat",     col = "logFC"),
  Inter_lfc   = sapply(marker_genes, get_stat, contrast = "Interaction", col = "logFC"),
  Inter_adjP  = sapply(marker_genes, get_stat, contrast = "Interaction", col = "adj.P.Val"),
  row.names = NULL, stringsAsFactors = FALSE
)
key_tbl$Inter_sig <- key_tbl$Inter_adjP < DE_FDR
cat("\n=== KEY RESULT: marker-gene cGAS-dependence (Interaction = cGAS-dependence test) ===\n")
print(format(key_tbl, digits = 3, nsmall = 3))

# also save it as a stage table
key_path <- file.path(stage_dir("03_de", "tables"), "marker_cgas_dependence.csv")
write.csv(key_tbl, key_path, row.names = FALSE)
cat("[SAVE]", key_path, "\n")

# -----------------------------------------------------------------------------
# 9. PLOT-READY TIDY TABLE for Fig 2 (so the viz does ZERO statistics)
# -----------------------------------------------------------------------------
# Per-GROUP mean log2(CPM+0.5) for the Fig 2 marker set (ISG arm + HIF/glyco arm),
# in tidy long form, JOINED with each gene's Interaction logFC + adj.P.Val. The
# viz reads this verbatim -- no group-mean recompute, no topTable lookup.
isg <- intersect(ISG_MARKERS, rownames(logcpm_mat))
hif <- intersect(HIF_GLYCO_MARKERS, rownames(logcpm_mat))
hif_specific <- c("Slc2a1", "Vegfa", "Egln3", "Bnip3")
fig_genes <- c(isg, hif)

fig2_marker_means <- do.call(rbind, lapply(fig_genes, function(g) {
  gm <- grp_mean(g)
  data.frame(
    gene         = g,
    group        = names(gm),
    mean_log2cpm = as.numeric(gm),
    stringsAsFactors = FALSE
  )
})) %>%
  mutate(
    arm       = ifelse(gene %in% isg, "IFN / ISG arm", "HIF / glycolysis arm"),
    hif_class = case_when(
      gene %in% isg          ~ NA_character_,
      gene %in% hif_specific ~ "HIF-specific",
      TRUE                   ~ "shared-glycolytic"
    ),
    genotype  = ifelse(grepl("^WT", group), "WT", "cGASKO"),
    temp      = ifelse(grepl("_39$", group), "39C", "37C"),
    # join each gene's Interaction logFC + adj.P.Val so the figure annotation
    # (int.adjP=..) needs no recompute downstream.
    inter_logFC   = sapply(gene, get_stat, contrast = "Interaction", col = "logFC"),
    inter_adjP    = sapply(gene, get_stat, contrast = "Interaction", col = "adj.P.Val")
  ) %>%
  # column order per the documented schema
  select(gene, arm, hif_class, group, genotype, temp, mean_log2cpm,
         inter_logFC, inter_adjP)

fig2_means_path <- file.path(stage_dir("03_de", "tables"), "fig2_marker_means.csv")
write.csv(fig2_marker_means, fig2_means_path, row.names = FALSE)
cat("[SAVE]", fig2_means_path, "--", nrow(fig2_marker_means), "rows,",
    length(unique(fig2_marker_means$gene)), "genes\n")

cat("\n=== Phase 2 COMPUTE complete:", length(de_results),
    "contrasts | tidy tables emitted for viz ===\n")
