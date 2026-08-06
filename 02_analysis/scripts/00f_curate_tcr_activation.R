#!/usr/bin/env Rscript
# 00f_curate_tcr_activation.R - Freeze the curated TCR/IEG T-cell activation lens
# Project: STING-cGAS-GSE329522
# Phase: 0 (infrastructure / curated reference asset)
# Inputs:
#   00_data/references/gene_sets/tcr_activation_lens/tcr_activation_panel.csv  (frozen HUMAN panel)
#   msigdbr offline package data, species="Mus musculus"  (native human->mouse conversion)
#   00_data/processed/GSE329522_normalized_counts_CPM_iTreg.csv  (this dataset's gene_name column)
# Outputs:
#   00_data/references/gene_sets/tcr_activation_lens/tcr_activation_human.rds
#   00_data/references/gene_sets/tcr_activation_lens/tcr_activation_mouse.rds
#   00_data/references/gene_sets/tcr_activation_lens/tcr_activation_ortholog_map.csv
#   00_data/references/gene_sets/tcr_activation_lens/tables/tcr_activation_setsizes.csv
# Dependencies: 02_analysis/config/config.R; msigdbr (offline MSigDB)
#
# ===========================================================================================
# METHOD & RATIONALE
# ===========================================================================================
#   (A) ACTIVATION POLE DEFINITION
#       The source panel is a tight, literature-grounded HUMAN T-cell activation panel
#       spanning TCR-proximal signaling, early costimulation, immediate-early transcription
#       factors and activation effector genes. The frozen CSV is curated input and is read
#       verbatim.
#
#   (B) MOUSE BASIS
#       Mouse symbols come from msigdbr(species="Mus musculus") on the same native
#       db_gene_symbol -> gene_symbol convention as the HSR lens. For this panel the
#       validated map is 1:1: 66 HUMAN symbols -> 66 unique mouse symbols, all present in
#       GSE329522.
#
#   (C) DISJOINT FROM HSR BY CONSTRUCTION
#       This activation pole is built separately from the frozen HSR lens, so the
#       three-lens decomposition can tell the two memberships apart.
#
#   (D) HONEST FRAMING
#       This is the ACTIVATION pole. Overlap of the empirical WT_heat_up with it indicates
#       activation; overlap with the HSR core indicates HSR-core membership. Fever causality
#       for WT_heat_up is a further claim in either case.
#
#   DETERMINISM: no RNG is used; all symbols and rows are sorted before writing.
# ===========================================================================================

source("02_analysis/config/config.R")
suppressPackageStartupMessages({ library(msigdbr) })

OUT_DIR <- file.path(PROJECT_ROOT, "00_data/references/gene_sets/tcr_activation_lens")
TBL_DIR <- file.path(OUT_DIR, "tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)

PANEL_FILE <- file.path(OUT_DIR, "tcr_activation_panel.csv")
EXPECTED_N <- 66L
EXPECTED_CATEGORIES <- c(
  "tcr_proximal",
  "costim_early",
  "immediate_early_tf",
  "activation_effector"
)
HONEST_NOTE <- paste(
  "This is the ACTIVATION pole; overlap of the empirical WT_heat_up with it indicates",
  "activation, overlap with the HSR core indicates HSR-core membership - neither makes",
  "WT_heat_up causal for fever."
)

stop_mismatch <- function(label, observed, expected) {
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(sprintf("%s mismatch: observed=%d expected=%d", label, observed, expected))
  }
}

stop_with_values <- function(prefix, values) {
  stop(sprintf("%s: %s", prefix, paste(sort(unique(values)), collapse = ", ")))
}

required_panel_cols <- c("gene", "category", "rationale", "key_citation", "curator_note", "source")
panel <- read.csv(PANEL_FILE, stringsAsFactors = FALSE, check.names = FALSE)
missing_panel_cols <- setdiff(required_panel_cols, names(panel))
if (length(missing_panel_cols)) {
  stop_with_values("Panel is missing required column(s)", missing_panel_cols)
}

stop_mismatch("Panel row count", nrow(panel), EXPECTED_N)
if (anyDuplicated(panel$gene)) {
  stop_with_values("Panel contains duplicated HUMAN gene symbol(s)", panel$gene[duplicated(panel$gene)])
}

observed_categories <- sort(unique(panel$category))
if (!identical(observed_categories, sort(EXPECTED_CATEGORIES))) {
  stop(sprintf("Panel categories mismatch: observed={%s}; expected={%s}",
               paste(observed_categories, collapse = ", "),
               paste(sort(EXPECTED_CATEGORIES), collapse = ", ")))
}

panel_genes <- sort(unique(panel$gene))
stop_mismatch("Unique HUMAN panel gene count", length(panel_genes), EXPECTED_N)

# ---------------------------------------------------------------------------
# 1) Build a GENERAL human->mouse map from msigdbr's mouse-native conversion,
#    then restrict to this frozen HUMAN panel.
# ---------------------------------------------------------------------------
mm <- msigdbr::msigdbr(species = "Mus musculus")
if (!all(c("db_gene_symbol", "gene_symbol") %in% names(mm))) {
  stop("Mouse msigdbr output must contain columns 'db_gene_symbol' and 'gene_symbol'.")
}

ortholog_map <- unique(as.data.frame(mm[, c("db_gene_symbol", "gene_symbol")], stringsAsFactors = FALSE))
names(ortholog_map) <- c("human_symbol", "mouse_symbol")
ortholog_map <- ortholog_map[
  ortholog_map$human_symbol %in% panel_genes &
    !is.na(ortholog_map$human_symbol) & ortholog_map$human_symbol != "" &
    !is.na(ortholog_map$mouse_symbol) & ortholog_map$mouse_symbol != "",
  ,
  drop = FALSE
]
ortholog_map <- ortholog_map[order(ortholog_map$human_symbol, ortholog_map$mouse_symbol), , drop = FALSE]
rownames(ortholog_map) <- NULL

missing_human_map <- setdiff(panel_genes, unique(ortholog_map$human_symbol))
if (length(missing_human_map)) {
  stop_with_values("Panel HUMAN gene(s) missing from msigdbr mouse-native map", missing_human_map)
}
stop_mismatch("HUMAN genes covered by msigdbr map", length(unique(ortholog_map$human_symbol)), EXPECTED_N)

human_map_counts <- table(ortholog_map$human_symbol)
multi_human <- names(human_map_counts)[human_map_counts > 1L]
if (length(multi_human)) {
  detail <- vapply(multi_human, function(h) {
    sprintf("%s -> %s", h, paste(ortholog_map$mouse_symbol[ortholog_map$human_symbol == h], collapse = ";"))
  }, character(1))
  stop_with_values("Panel has HUMAN symbols mapping to multiple mouse symbols", detail)
}

mouse_map_counts <- table(ortholog_map$mouse_symbol)
multi_mouse <- names(mouse_map_counts)[mouse_map_counts > 1L]
if (length(multi_mouse)) {
  detail <- vapply(multi_mouse, function(m) {
    sprintf("%s <- %s", m, paste(ortholog_map$human_symbol[ortholog_map$mouse_symbol == m], collapse = ";"))
  }, character(1))
  stop_with_values("Panel has mouse symbols mapped from multiple HUMAN symbols", detail)
}

mouse_syms <- sort(unique(ortholog_map$mouse_symbol))
stop_mismatch("Unique mouse ortholog count", length(mouse_syms), EXPECTED_N)

# ---------------------------------------------------------------------------
# 2) Validate mapped mouse symbols against THIS dataset's gene_name column.
# ---------------------------------------------------------------------------
cpm <- read.csv(FILE_CPM, check.names = FALSE)
if (!"gene_name" %in% names(cpm)) {
  stop(sprintf("FILE_CPM must contain a 'gene_name' column: %s", FILE_CPM))
}
dataset_symbols <- unique(cpm$gene_name)
ortholog_map$in_GSE329522 <- ortholog_map$mouse_symbol %in% dataset_symbols
missing_dataset <- ortholog_map$mouse_symbol[!ortholog_map$in_GSE329522]
if (length(missing_dataset)) {
  stop_with_values("Mapped mouse ortholog(s) absent from GSE329522 CPM gene_name", missing_dataset)
}
stop_mismatch("Mouse orthologs present in GSE329522", sum(ortholog_map$in_GSE329522), EXPECTED_N)

# ---------------------------------------------------------------------------
# 3) Save named-list assets with provenance attributes.
# ---------------------------------------------------------------------------
category_counts <- as.integer(table(factor(panel$category, levels = EXPECTED_CATEGORIES)))
names(category_counts) <- EXPECTED_CATEGORIES

provenance <- list(
  signature = "curated TCR/IEG T-cell activation lens",
  panel_source = "tcr_activation_panel.csv (agy_gemini_3.1_pro_web_2026-07-15 + curator: FOXP3 dropped)",
  ortholog_source = "msigdbr species=Mus musculus (native conversion; 1:1 for this panel)",
  msigdbr_version = as.character(utils::packageVersion("msigdbr")),
  n_human = EXPECTED_N,
  n_mouse = EXPECTED_N,
  n_mouse_in_dataset = EXPECTED_N,
  categories = category_counts,
  built_by = "02_analysis/scripts/00f_curate_tcr_activation.R",
  note = HONEST_NOTE
)

human_set <- structure(
  list(TCR_activation = structure(sort(unique(panel$gene)), provenance = provenance)),
  provenance = provenance
)
mouse_set <- structure(
  list(TCR_activation = structure(mouse_syms, provenance = provenance)),
  provenance = provenance
)

saveRDS(human_set, file.path(OUT_DIR, "tcr_activation_human.rds"))
saveRDS(mouse_set, file.path(OUT_DIR, "tcr_activation_mouse.rds"))

# ---------------------------------------------------------------------------
# 4) Write per-pair map and per-category accounting table.
# ---------------------------------------------------------------------------
write.csv(ortholog_map, file.path(OUT_DIR, "tcr_activation_ortholog_map.csv"), row.names = FALSE)

panel_map <- merge(
  panel[, c("gene", "category"), drop = FALSE],
  ortholog_map,
  by.x = "gene",
  by.y = "human_symbol",
  all.x = TRUE,
  sort = FALSE
)

category_rows <- do.call(rbind, lapply(EXPECTED_CATEGORIES, function(category_name) {
  rows <- panel_map[panel_map$category == category_name, , drop = FALSE]
  data.frame(
    category = category_name,
    n_human = length(unique(rows$gene)),
    n_mouse = length(unique(rows$mouse_symbol)),
    n_mouse_in_GSE329522 = length(unique(rows$mouse_symbol[rows$in_GSE329522])),
    stringsAsFactors = FALSE
  )
}))
total_row <- data.frame(
  category = "TCR_activation",
  n_human = length(unique(panel_map$gene)),
  n_mouse = length(unique(panel_map$mouse_symbol)),
  n_mouse_in_GSE329522 = length(unique(panel_map$mouse_symbol[panel_map$in_GSE329522])),
  stringsAsFactors = FALSE
)
setsizes <- rbind(category_rows, total_row)
write.csv(setsizes, file.path(TBL_DIR, "tcr_activation_setsizes.csv"), row.names = FALSE)

cat("\n================ TCR/IEG activation lens ================\n")
cat(sprintf("Panel HUMAN genes:           %d\n", length(panel_genes)))
cat(sprintf("Mouse orthologs:             %d\n", length(mouse_syms)))
cat(sprintf("Mouse orthologs in GSE329522: %d\n", sum(ortholog_map$in_GSE329522)))
cat("\nPer-category counts:\n")
for (i in seq_len(nrow(category_rows))) {
  cat(sprintf("  %-20s human=%2d mouse=%2d mouse_in_GSE329522=%2d\n",
              category_rows$category[[i]],
              category_rows$n_human[[i]],
              category_rows$n_mouse[[i]],
              category_rows$n_mouse_in_GSE329522[[i]]))
}
cat(sprintf("\nHonest note: %s\n", HONEST_NOTE))
cat("\nWrote:\n")
cat("  ", file.path(OUT_DIR, "tcr_activation_human.rds"), "\n")
cat("  ", file.path(OUT_DIR, "tcr_activation_mouse.rds"), "\n")
cat("  ", file.path(OUT_DIR, "tcr_activation_ortholog_map.csv"), "\n")
cat("  ", file.path(TBL_DIR, "tcr_activation_setsizes.csv"), "\n")
cat("==========================================================\n")
