#!/usr/bin/env Rscript
# =============================================================================
# 01_mapping_qc.R  --  PHASE 1: Sample-mapping QC + scramble exhibit (COMPUTE)
# =============================================================================
# Phase:        1 (SCIENCE GATE -- go/no-go on the inferred sample mapping)
# Inputs:       00_data/processed/GSE329522_normalized_counts_CPM_iTreg.csv
#               03_results/objects/sample_metadata.rds  (Phase 0)
# Outputs:      03_results/objects/01_eda.rds  (cpm_mat, logcpm_mat, collapse
#                 record, aligned metadata, cgas_symbol, + derived plot objects)
#               03_results/01_qc/tables/mapping_verdict.csv
#               03_results/01_qc/README.md  (verdict paragraph)
#               -- PLOT-READY TIDY TABLES (consumed by 01_mapping_qc_viz.R) --
#               03_results/01_qc/tables/fig1a_thermometer_data.csv
#               03_results/01_qc/tables/fig1b_cgas_data.csv
#               03_results/02_eda/tables/fig1c_pca_data.csv
#               03_results/02_eda/tables/fig1c_pca_varexp.csv
#               03_results/01_qc/tables/fig1d_scramble_data.csv
# Dependencies: config.R; dplyr, tidyr (via load_packages)
#
# NORMALIZE-THEN-VISUALIZE: this script does ALL computation and emits tidy,
# plot-ready tables. It contains NO ggplot/ggsave. The companion viz script
# 01_mapping_qc_viz.R reads these tables and renders the 4 figures, doing zero
# statistics.
#
# Decisions made here once for the whole project:
#   * NA CPM cells are treated as 0 for matrix building (a missing CPM in a
#     normalized-counts deposit means "not detected"; 0 is the correct fill for
#     log2(CPM+0.5) and for variance/PCA). Reported to console.
#   * Duplicate gene_name -> collapse to the SINGLE row with the MAX MEAN-CPM
#     across the 20 samples (standard "most-expressed isoform wins"). The winning
#     Ensembl id per collapsed symbol is recorded.
# =============================================================================

source("02_analysis/config/config.R")
load_packages(extra = c("tidyr"))

set.seed(GSEA_SEED)

# -----------------------------------------------------------------------------
# A) LOAD + BUILD MATRICES (+ duplicate-symbol collapse decision)
# -----------------------------------------------------------------------------
message("\n==================  A) LOAD + BUILD MATRICES  ==================")

meta <- readRDS(file.path(DIR_OBJECTS, "sample_metadata.rds"))
lib_order <- meta$library_id                          # 021..040, temperature-major

raw <- read.csv(FILE_CPM, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot("gene_id" %in% names(raw), "gene_name" %in% names(raw))

sample_cols <- setdiff(names(raw), c("gene_id", "gene_name"))
stopifnot(setequal(sample_cols, lib_order))           # exactly the 20 libraries
# Reorder sample columns to metadata library_id order.
sample_cols <- lib_order

# Coerce the 20 sample columns to numeric.
for (cc in sample_cols) raw[[cc]] <- suppressWarnings(as.numeric(raw[[cc]]))

# --- NA report ---
n_na_cells <- sum(is.na(as.matrix(raw[, sample_cols])))
n_rows_with_na <- sum(apply(raw[, sample_cols], 1, function(x) any(is.na(x))))
message(sprintf("[NA] %d NA cells across %d gene rows (of %d). Policy: NA -> 0 for matrix building.",
                n_na_cells, n_rows_with_na, nrow(raw)))
raw[sample_cols][is.na(raw[sample_cols])] <- 0        # NA -> 0

# --- duplicate gene_name collapse: max mean-CPM wins ---
raw$.mean_cpm <- rowMeans(as.matrix(raw[, sample_cols]))
dup_symbols <- unique(raw$gene_name[duplicated(raw$gene_name)])
message(sprintf("[COLLAPSE] %d gene rows; %d gene_name symbols are duplicated (will collapse to max-mean-CPM row).",
                nrow(raw), length(dup_symbols)))

# Record, per duplicated symbol, all Ensembl ids and which won.
collapse_record <- raw %>%
  dplyr::filter(gene_name %in% dup_symbols) %>%
  dplyr::group_by(gene_name) %>%
  dplyr::arrange(dplyr::desc(.mean_cpm), .by_group = TRUE) %>%
  dplyr::summarise(
    n_ids        = dplyr::n(),
    winner_id    = gene_id[1],
    winner_meancpm = .mean_cpm[1],
    loser_ids    = paste(gene_id[-1], collapse = ";"),
    .groups = "drop"
  )

# Keep, for every symbol, the single max-mean-CPM row.
collapsed <- raw %>%
  dplyr::group_by(gene_name) %>%
  dplyr::slice_max(order_by = .mean_cpm, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

stopifnot(!any(duplicated(collapsed$gene_name)))

# --- marker-survival check ---
all_markers <- unique(c(THERMO_MARKERS, CGAS_GENE, ISG_MARKERS, HIF_GLYCO_MARKERS))
markers_present_pre  <- intersect(all_markers, raw$gene_name)
markers_present_post <- intersect(all_markers, collapsed$gene_name)
lost_markers <- setdiff(markers_present_pre, markers_present_post)
message(sprintf("[MARKER CHECK] %d/%d markers present pre-collapse, %d post-collapse.",
                length(markers_present_pre), length(all_markers), length(markers_present_post)))
if (length(lost_markers) > 0) {
  stop(sprintf("FATAL: marker(s) LOST in collapse: %s", paste(lost_markers, collapse = ", ")))
} else {
  message("[MARKER CHECK] PASS -- no marker gene lost in the collapse.")
}
# Explicit Aldoa note (a duplicated glycolysis symbol central to forensics).
if ("Aldoa" %in% dup_symbols) {
  aw <- collapse_record$winner_id[collapse_record$gene_name == "Aldoa"]
  message(sprintf("[MARKER CHECK] Aldoa is duplicated; kept Ensembl id %s (max mean-CPM).", aw))
}
cgas_symbol <- intersect(CGAS_GENE, collapsed$gene_name)[1]
message(sprintf("[MARKER CHECK] Cgas gene present as: %s", cgas_symbol))

# --- build matrices (genes x 20, rownames = gene_name) ---
cpm_mat <- as.matrix(collapsed[, sample_cols])
rownames(cpm_mat) <- collapsed$gene_name
colnames(cpm_mat) <- sample_cols                      # already lib_order
logcpm_mat <- log2(cpm_mat + 0.5)

# Verify column order matches metadata; align metadata rows to columns.
stopifnot(identical(colnames(cpm_mat), meta$library_id))
rownames(meta) <- meta$library_id
meta <- meta[colnames(cpm_mat), , drop = FALSE]
message(sprintf("[MATRIX] cpm_mat %d x %d ; columns aligned to metadata library_id order.",
                nrow(cpm_mat), ncol(cpm_mat)))

# Shared label scaffolding (used by derived tables; cosmetic only).
lib_short <- sub("^12630-RS-", "", meta$library_id)    # 021..040 short labels

# =============================================================================
# B) DERIVED: FIG 1a THERMOMETER DATA (label-blind temperature validation)
# =============================================================================
message("\n==================  B) FIG 1a THERMOMETER (compute)  ==================")

thermo_present <- intersect(THERMO_MARKERS, rownames(logcpm_mat))
thermo_long <- as.data.frame(t(logcpm_mat[thermo_present, , drop = FALSE]))
thermo_long$library_id <- rownames(thermo_long)
thermo_long <- thermo_long %>%
  tidyr::pivot_longer(dplyr::all_of(thermo_present),
                      names_to = "gene", values_to = "log2cpm") %>%
  dplyr::left_join(meta[, c("library_id", "temp", "genotype", "group")],
                   by = "library_id")
thermo_long$lib <- sub("^12630-RS-", "", thermo_long$library_id)

# Plot-ready tidy table: (library_id, marker, value, inferred_temp, inferred_group)
fig1a_data <- data.frame(
  library_id    = thermo_long$library_id,
  lib           = thermo_long$lib,
  marker        = as.character(thermo_long$gene),
  value         = thermo_long$log2cpm,
  inferred_temp = as.character(thermo_long$temp),
  inferred_group = as.character(thermo_long$group),
  stringsAsFactors = FALSE
)

# Quantify thermometer (acceptance #1): hot (031-040) vs cool (021-030) means.
hot_libs  <- meta$library_id[meta$temp == "39"]
cool_libs <- meta$library_id[meta$temp == "37"]
thermo_quant <- data.frame(
  gene = thermo_present,
  mean_cool_log2 = rowMeans(logcpm_mat[thermo_present, cool_libs, drop = FALSE]),
  mean_hot_log2  = rowMeans(logcpm_mat[thermo_present, hot_libs,  drop = FALSE])
)
thermo_quant$delta_hot_minus_cool <- thermo_quant$mean_hot_log2 - thermo_quant$mean_cool_log2
message("[ACCEPT #1] Thermometer hot vs cool (log2CPM means):")
print(thermo_quant, row.names = FALSE)

# =============================================================================
# C) DERIVED: FIG 1b CGAS GENOTYPE-CHECK DATA
# =============================================================================
message("\n==================  C) FIG 1b CGAS (compute)  ==================")

# Plot-ready tidy table: (library_id, cgas_value, inferred_genotype, inferred_temp)
fig1b_data <- data.frame(
  library_id       = meta$library_id,
  lib              = lib_short,
  cgas_value       = logcpm_mat[cgas_symbol, meta$library_id],
  inferred_genotype = as.character(meta$genotype),
  inferred_temp    = as.character(meta$temp),
  inferred_group   = as.character(meta$group),
  stringsAsFactors = FALSE
)

# Quantify Cgas WT vs KO within each temp half (acceptance #2).
cgas_quant <- fig1b_data %>%
  dplyr::rename(genotype = inferred_genotype, temp = inferred_temp,
                log2cpm = cgas_value) %>%
  dplyr::group_by(temp, genotype) %>%
  dplyr::summarise(mean_log2cpm = mean(log2cpm), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = genotype, values_from = mean_log2cpm) %>%
  dplyr::mutate(WT_minus_KO = WT - cGASKO)
message("[ACCEPT #2] Cgas WT vs cGASKO within each temperature half (log2CPM means):")
print(as.data.frame(cgas_quant), row.names = FALSE)

# =============================================================================
# D) DERIVED: FIG 1c PCA DATA (label-blind, overlay inferred labels)
# =============================================================================
message("\n==================  D) FIG 1c PCA 2x2 (compute)  ==================")

# The PCA runs on the SAME gene universe the DE stage models: every delivered
# symbol, unfiltered, minus the handful that are constant across all 20 libraries
# and so carry no variance to decompose. Selecting the most variable genes first
# would concentrate variance into PC1 by construction, which inflates the % var
# the figure reports without changing the axis it finds; the correspondence check
# below quantifies exactly that.
gene_var  <- apply(logcpm_mat, 1, var)
pca_genes <- names(gene_var)[gene_var > 0]
n_genes   <- length(pca_genes)
n_zerovar <- sum(gene_var == 0)
pca <- prcomp(t(logcpm_mat[pca_genes, ]), center = TRUE, scale. = FALSE)
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
message(sprintf("[PCA] %d genes entered (%d dropped for zero variance across all libraries).",
                n_genes, n_zerovar))

# Correspondence check: a top-2000-variable-gene PCA finds the same leading axis,
# at a much larger apparent % var. Carried as figure provenance so the caption can
# state the number instead of asserting the equivalence.
top2000    <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(2000, n_genes))]
pca_top    <- prcomp(t(logcpm_mat[top2000, ]), center = TRUE, scale. = FALSE)
pv_top     <- pca_top$sdev^2 / sum(pca_top$sdev^2) * 100
pc1_cor_top2000 <- abs(cor(pca$x[, 1], pca_top$x[, 1]))
pc2_cor_top2000 <- abs(cor(pca$x[, 2], pca_top$x[, 2]))
message(sprintf(
  "[ACCEPT #3b] top-2000 subset PCA: PC1 %.1f%% / PC2 %.1f%% vs %.1f%% / %.1f%% on all %d genes; |r| PC1 %.5f, PC2 %.3f.",
  pv_top[1], pv_top[2], percentVar[1], percentVar[2], n_genes,
  pc1_cor_top2000, pc2_cor_top2000))

# Plot-ready tidy table: (library_id, PC1, PC2, genotype, temp)
fig1c_data <- data.frame(
  library_id = rownames(pca$x),
  PC1 = pca$x[, 1], PC2 = pca$x[, 2],
  lib = lib_short,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(meta[, c("library_id", "temp", "genotype", "group")], by = "library_id")
fig1c_data$temp     <- as.character(fig1c_data$temp)
fig1c_data$genotype <- as.character(fig1c_data$genotype)
fig1c_data$group    <- as.character(fig1c_data$group)

# Small variance-explained table: (PC, pct_var). The gene count that entered the
# PCA and the top-2000 correspondence are carried as label columns on every row so
# the viz title and caption read them instead of hardcoding them.
fig1c_varexp <- data.frame(
  PC = paste0("PC", seq_along(percentVar)),
  pct_var = percentVar,
  n_genes = n_genes,
  n_zerovar_dropped = n_zerovar,
  pc1_pct_var_top2000 = pv_top[1],
  pc1_cor_vs_top2000 = pc1_cor_top2000,
  stringsAsFactors = FALSE
)

# Which axis tracks temperature vs genotype? (acceptance #3)
pc1_temp  <- summary(lm(PC1 ~ temp,     data = fig1c_data))$r.squared
pc1_geno  <- summary(lm(PC1 ~ genotype, data = fig1c_data))$r.squared
pc2_temp  <- summary(lm(PC2 ~ temp,     data = fig1c_data))$r.squared
pc2_geno  <- summary(lm(PC2 ~ genotype, data = fig1c_data))$r.squared
message(sprintf("[ACCEPT #3] PC1 = %.1f%% var; PC2 = %.1f%% var.", percentVar[1], percentVar[2]))
message(sprintf("[ACCEPT #3] R^2 of PC1 ~ temp = %.2f, PC1 ~ genotype = %.2f", pc1_temp, pc1_geno))
message(sprintf("[ACCEPT #3] R^2 of PC2 ~ temp = %.2f, PC2 ~ genotype = %.2f", pc2_temp, pc2_geno))

# =============================================================================
# E) DERIVED: FIG 1d SCRAMBLE EXHIBIT DATA (competence exhibit)
# =============================================================================
message("\n==================  E) FIG 1d SCRAMBLE (compute)  ==================")

# Map each GSM (true, marker-derived) to the group it actually belongs to, and
# each GSM (naive positional join) to the group a column-order GSM join assigns.
# We display, per library/column, the TRUE inferred condition vs the condition a
# naive accession-order join would assign; discordant libraries are highlighted.
gsm_to_truegroup <- setNames(as.character(meta$group), meta$gsm_id)

scram <- data.frame(
  library_id          = meta$library_id,
  lib                 = lib_short,
  inferred_group      = as.character(meta$group),               # marker-derived (TRUE column condition)
  gsm_positional      = meta$gsm_id_positional,                 # naive accession-order GSM
  thermo_score        = colMeans(logcpm_mat[thermo_present, meta$library_id, drop = FALSE]),
  stringsAsFactors    = FALSE
)
# Condition a naive positional GSM->column join would STAMP onto this column:
scram$naive_join_group <- gsm_to_truegroup[scram$gsm_positional]
scram$discordant <- scram$inferred_group != scram$naive_join_group

n_disc <- sum(scram$discordant)
disc_libs <- sub("^12630-RS-", "", scram$library_id[scram$discordant])
message(sprintf("[SCRAMBLE] %d/20 libraries are MISLABELED by a naive positional GSM->column join: %s",
                n_disc, paste(disc_libs, collapse = ", ")))

# Plot-ready tidy table: (library_id, inferred_condition, naive_positional_condition,
# discordant, thermo_score).
fig1d_data <- data.frame(
  library_id                 = scram$library_id,
  lib                        = scram$lib,
  inferred_condition         = scram$inferred_group,
  naive_positional_condition = scram$naive_join_group,
  discordant                 = scram$discordant,
  thermo_score               = scram$thermo_score,
  stringsAsFactors           = FALSE
)

# =============================================================================
# F) CHECKPOINT (cpm/logcpm + collapse record + metadata + derived plot objects)
# =============================================================================
message("\n==================  F) CHECKPOINT  ==================")

eda <- load_or_compute("01_eda.rds", function() {
  list(
    cpm_mat         = cpm_mat,
    logcpm_mat      = logcpm_mat,
    metadata        = meta,
    collapse_record = collapse_record,
    na_policy       = "NA->0",
    n_na_cells      = n_na_cells,
    cgas_symbol     = cgas_symbol
  )
}, force = TRUE, desc = "01_eda (cpm/logcpm matrices + collapse record + metadata)")

# Convenience handles (use the checkpoint so re-runs are consistent).
cpm_mat    <- eda$cpm_mat
logcpm_mat <- eda$logcpm_mat
meta       <- eda$metadata
cgas_symbol <- eda$cgas_symbol

# --- emit the plot-ready tidy tables ---
qc_tables  <- stage_dir("01_qc", "tables")
eda_tables <- stage_dir("02_eda", "tables")

write.csv(fig1a_data,   file.path(qc_tables,  "fig1a_thermometer_data.csv"), row.names = FALSE)
write.csv(fig1b_data,   file.path(qc_tables,  "fig1b_cgas_data.csv"),        row.names = FALSE)
write.csv(fig1c_data,   file.path(eda_tables, "fig1c_pca_data.csv"),         row.names = FALSE)
write.csv(fig1c_varexp, file.path(eda_tables, "fig1c_pca_varexp.csv"),       row.names = FALSE)
write.csv(fig1d_data,   file.path(qc_tables,  "fig1d_scramble_data.csv"),    row.names = FALSE)
message("[TABLES] wrote fig1a/b/d -> 01_qc/tables, fig1c (+varexp) -> 02_eda/tables")

# =============================================================================
# G) MACHINE-CHECKABLE MAPPING VERDICT
# =============================================================================
message("\n==================  G) MAPPING VERDICT  ==================")

# Data-driven predicted temp: thermometer score (mean log2CPM of thermo markers).
# Threshold = midpoint between the two thermometer clusters (k-means, k=2) so it
# is fully data-driven, not the inferred label.
thermo_score <- colMeans(logcpm_mat[thermo_present, meta$library_id, drop = FALSE])
km_t <- kmeans(matrix(thermo_score, ncol = 1), centers = 2, nstart = 25)
hot_cluster <- which.max(km_t$centers[, 1])
pred_temp <- ifelse(km_t$cluster == hot_cluster, "39", "37")
names(pred_temp) <- meta$library_id

# Data-driven predicted genotype: Cgas value, split WITHIN each predicted-temp
# half (KO has lower Cgas). Use median split within each half (n=5 vs 5).
cgas_val <- logcpm_mat[cgas_symbol, meta$library_id]
names(cgas_val) <- meta$library_id
pred_geno <- rep(NA_character_, length(cgas_val)); names(pred_geno) <- meta$library_id
for (tt in c("37", "39")) {
  in_half <- names(pred_temp)[pred_temp == tt]
  thr <- median(cgas_val[in_half])
  pred_geno[in_half] <- ifelse(cgas_val[in_half] >= thr, "WT", "cGASKO")
}

verdict <- data.frame(
  library_id         = meta$library_id,
  inferred_group     = as.character(meta$group),
  thermo_score       = round(thermo_score, 3),
  cgas_value         = round(cgas_val, 3),
  predicted_temp     = pred_temp[meta$library_id],
  predicted_genotype = pred_geno[meta$library_id],
  stringsAsFactors   = FALSE
)
verdict$inferred_temp     <- as.character(meta$temp)
verdict$inferred_genotype <- as.character(meta$genotype)
verdict$concordant <- with(verdict,
  predicted_temp == inferred_temp & predicted_genotype == inferred_genotype)
verdict$notes <- with(verdict, ifelse(concordant, "",
  paste0(
    ifelse(predicted_temp     != inferred_temp,     "temp_mismatch ", ""),
    ifelse(predicted_genotype != inferred_genotype, "genotype_mismatch", "")
  )))
verdict$notes <- trimws(verdict$notes)

verdict_out <- verdict[, c("library_id", "inferred_group", "thermo_score",
                           "cgas_value", "predicted_temp", "predicted_genotype",
                           "concordant", "notes")]
verdict_path <- file.path(stage_dir("01_qc", "tables"), "mapping_verdict.csv")
write.csv(verdict_out, verdict_path, row.names = FALSE)
message(sprintf("[TABLE] saved mapping_verdict.csv (%s)", verdict_path))

n_concordant <- sum(verdict$concordant)
disc_rows <- verdict_out$library_id[!verdict$concordant]
message(sprintf("[ACCEPT #4] %d/20 libraries concordant (data-derived label == inferred label).",
                n_concordant))
if (length(disc_rows) > 0) {
  message(sprintf("[ACCEPT #4] DISCORDANT libraries (FLAGGED, not smoothed): %s",
                  paste(disc_rows, collapse = ", ")))
}
message("\n--- mapping_verdict.csv ---")
print(verdict_out, row.names = FALSE)

# --- verdict paragraph to stage README ---
gate <- if (n_concordant == 20) "SUPPORTED (gate PASS)" else "CONCERN -- discordant libraries present"
thermo_mean_hot  <- mean(thermo_quant$mean_hot_log2)
thermo_mean_cool <- mean(thermo_quant$mean_cool_log2)

readme_path <- file.path(DIR_RESULTS, "01_qc", "README.md")
dir.create(dirname(readme_path), recursive = TRUE, showWarnings = FALSE)
readme_txt <- sprintf(
"# Sample-mapping QC of the iTreg libraries

**Generated:** %s by `02_analysis/scripts/01_mapping_qc.R`

## Verdict (mapping %s)

The temperature-major mapping under test (021-025 WT_37, 026-030 cGASKO_37, 031-035 WT_39, 036-040 cGASKO_39) is **%s**: **%d/20** libraries have a data-derived label (thermometer-predicted temperature + Cgas-predicted genotype) that matches it.%s

**Evidence.** The heat-shock thermometer (%s) is monotone with the assigned temperature: mean log2CPM = %.2f in the hot half (031-040) vs %.2f in the cool half (021-030), a +%.2f log2 shift. %s (Cgas) is higher in WT than cGAS-KO within both temperature halves (37 °C: WT-KO = %.2f; 39 °C: WT-KO = %.2f log2CPM). PCA on all %d genes carrying variance places %.1f%%/%.1f%% of variance on PC1/PC2. PC1 is temperature almost exactly (R2 %.2f against temperature, %.2f against genotype); genotype sits weakly on PC2 (R2 %.2f), so the 2x2 is recoverable from expression but its two factors are recoverable to very different degrees.

**Scramble caveat.** The deposited CPM column order is temperature-major; the GEO GSM accessions are genotype-major. A naive positional GSM->column join therefore mislabels **%d** libraries (%s) -- see `figures/fig1d_scramble.png`. This is why the mapping must be marker-derived, not accession-positional.

**Bottom line:** the label-blind mapping is **%s**. %s

See `tables/mapping_verdict.csv` for the per-library machine-checkable verdict.
",
  format(Sys.time(), "%Y-%m-%d %H:%M"),
  tolower(sample_mapping_status()),
  gate, n_concordant,
  if (length(disc_rows) > 0) sprintf(" Discordant (flagged): %s.", paste(disc_rows, collapse = ", ")) else " No discordant libraries.",
  paste(thermo_present, collapse = "/"),
  thermo_mean_hot, thermo_mean_cool, thermo_mean_hot - thermo_mean_cool,
  cgas_symbol,
  cgas_quant$WT_minus_KO[cgas_quant$temp == "37"],
  cgas_quant$WT_minus_KO[cgas_quant$temp == "39"],
  n_genes, percentVar[1], percentVar[2], pc1_temp, pc1_geno, pc2_geno,
  n_disc, paste(disc_libs, collapse = ", "),
  if (n_concordant == 20) "SUPPORTED" else "NOT fully supported",
  sample_mapping_caption()
)
writeLines(readme_txt, readme_path)
message(sprintf("[README] wrote verdict to %s", readme_path))

# =============================================================================
# SUMMARY
# =============================================================================
message("\n==================  ACCEPTANCE SUMMARY  ==================")
message(sprintf("1. Thermometer hot>cool: +%.2f log2 (hot %.2f vs cool %.2f) -- %s",
                thermo_mean_hot - thermo_mean_cool, thermo_mean_hot, thermo_mean_cool,
                if (all(thermo_quant$delta_hot_minus_cool > 0)) "MONOTONE PASS" else "CHECK"))
message(sprintf("2. Cgas WT>KO: 37C +%.2f, 39C +%.2f -- %s",
                cgas_quant$WT_minus_KO[cgas_quant$temp == "37"],
                cgas_quant$WT_minus_KO[cgas_quant$temp == "39"],
                if (all(cgas_quant$WT_minus_KO > 0)) "PASS both halves" else "CHECK"))
message(sprintf("3. PCA PC1=%.1f%% PC2=%.1f%%; PC1~temp R2=%.2f PC2~geno R2=%.2f",
                percentVar[1], percentVar[2], pc1_temp, pc2_geno))
message(sprintf("4. Concordance: %d/20 -- %s", n_concordant,
                if (n_concordant == 20) "PASS" else "CONCERN"))
message(sprintf("5. Markers lost in collapse: %d -- %s",
                length(lost_markers), if (length(lost_markers) == 0) "PASS" else "FAIL"))
message(sprintf("6. Plot-ready tidy tables emitted; run 01_mapping_qc_viz.R to render fig1a-d (sample mapping: %s)",
                sample_mapping_status()))
message(sprintf("\nGATE RECOMMENDATION: %s", gate))
message("\n[DONE] 01_mapping_qc.R (compute) complete.")
