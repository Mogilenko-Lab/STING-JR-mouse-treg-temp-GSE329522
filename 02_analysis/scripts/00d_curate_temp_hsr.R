#!/usr/bin/env Rscript
# 00d_curate_temp_hsr.R - Freeze the curated heat-shock-response (HSR) lens
# Project: STING-cGAS-GSE329522
# Phase: 0 (infrastructure / curated reference asset)
# Inputs:
#   msigdbr offline package data, species="Homo sapiens", MSigDB v2026.1.Hs
#   00_data/processed/GSE329522_normalized_counts_CPM_iTreg.csv  (this dataset's gene_name column)
# Outputs:
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_human.rds
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_mouse.rds
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_mouse_sensitivity.rds
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_curation_log.csv
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_ortholog_map.csv
#   00_data/references/gene_sets/temp_hsr_lens/tables/temp_hsr_setsizes.csv
# Dependencies: 02_analysis/config/config.R; msigdbr (offline MSigDB)
#
# ===========================================================================================
# METHOD & RATIONALE
# ===========================================================================================
#   (A) HSR CORE DEFINITION
#       The curated HSR core is the UNION of three public MSigDB v2026.1.Hs sets:
#       REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS,
#       REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE, and GOBP_RESPONSE_TO_HEAT.
#       The shared HSF1 -> HSPA* -> DNAJ* -> BAG* -> HSPB* proteostasis core is the
#       temperature-specific signal this lens measures, so it is kept in full.
#
#   (B) HONEST CEILING
#       This clean core is proteotoxic-stress-general: HSF1 also fires on oxidative,
#       proteasome and metal stress, so fever specificity is a further claim. The mouse 37/39
#       contrast is what measures thermal-ness here.
#
#   (C) TWO FROZEN TIERS
#       sensitivity = the full mapped mouse union, before GSE329522 feature filtering.
#       core        = the anchor's expressed-filtered mouse set, the union intersected with
#                     genes present in the GSE329522 CPM gene_name column.
#       The human asset is the portable sensitivity union in HUMAN symbols, and each human
#       compartment derives its own expressed core from that union locally.
#
#   DETERMINISM: no RNG is used; all symbols are sorted before writing.
# ===========================================================================================

source("02_analysis/config/config.R")
suppressPackageStartupMessages({ library(msigdbr) })

OUT_DIR <- file.path(PROJECT_ROOT, "00_data/references/gene_sets/temp_hsr_lens")
TBL_DIR <- file.path(OUT_DIR, "tables")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)

TARGET_SETS <- c(
  "REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS",
  "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE",
  "GOBP_RESPONSE_TO_HEAT"
)
EXPECTED_N <- c(
  "REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS" = 101L,
  "REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE" = 82L,
  "GOBP_RESPONSE_TO_HEAT" = 104L
)
EXCLUDED_SETS <- c(
  "GOBP_DETECTION_OF_TEMPERATURE_STIMULUS",
  "GOBP_DETECTION_OF_TEMPERATURE_STIMULUS_INVOLVED_IN_THERMOCEPTION",
  "HP_FEVER"
)

HONEST_CEILING <- paste(
  "Even this clean core is proteotoxic-stress-general: HSF1 also fires on oxidative,",
  "proteasome, and metal stress, so it is not fever-specific; the mouse 37/39 contrast",
  "is the only thing that can measure thermal-ness."
)

# ---------------------------------------------------------------------------
# 1) Pull and validate the three HUMAN source sets from offline msigdbr data.
# ---------------------------------------------------------------------------
msig <- msigdbr::msigdbr(species = "Homo sapiens")
if (!all(c("gs_name", "gene_symbol") %in% names(msig))) {
  stop("msigdbr output must contain columns 'gs_name' and 'gene_symbol'.")
}

hsr_df <- msig[msig$gs_name %in% TARGET_SETS, c("gs_name", "gene_symbol"), drop = FALSE]
found_sets <- sort(unique(hsr_df$gs_name))
missing_sets <- setdiff(TARGET_SETS, found_sets)
if (length(missing_sets)) {
  stop(sprintf("Missing required MSigDB set(s): %s", paste(missing_sets, collapse = ", ")))
}

source_genes <- lapply(TARGET_SETS, function(set_name) {
  sort(unique(hsr_df$gene_symbol[hsr_df$gs_name == set_name]))
})
names(source_genes) <- TARGET_SETS

observed_n <- vapply(source_genes, length, integer(1))
bad_n <- names(observed_n)[observed_n != EXPECTED_N[names(observed_n)]]
if (length(bad_n)) {
  detail <- sprintf("%s observed=%d expected=%d",
                    bad_n, observed_n[bad_n], EXPECTED_N[bad_n])
  stop(sprintf(
    "MSigDB HSR source-set size drift detected; freeze no longer matches validated msigdbr release: %s",
    paste(detail, collapse = "; ")
  ))
}

human_union <- sort(unique(unlist(source_genes, use.names = FALSE)))
if (length(human_union) != 176L) {
  stop(sprintf("Human HSR union size drift detected: observed=%d expected=176", length(human_union)))
}

# ---------------------------------------------------------------------------
# 2) Per-gene membership flags across the three source sets.
# ---------------------------------------------------------------------------
membership <- data.frame(
  human_symbol = human_union,
  in_reactome_crhs = human_union %in% source_genes[["REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS"]],
  in_reactome_hsf1 = human_union %in% source_genes[["REACTOME_REGULATION_OF_HSF1_MEDIATED_HEAT_SHOCK_RESPONSE"]],
  in_gobp_response_to_heat = human_union %in% source_genes[["GOBP_RESPONSE_TO_HEAT"]],
  stringsAsFactors = FALSE
)
membership$n_sets <- rowSums(membership[, c(
  "in_reactome_crhs", "in_reactome_hsf1", "in_gobp_response_to_heat"
)])

# ---------------------------------------------------------------------------
# 3) Mouse-native msigdbr orthology, preserving 1:many paralog mappings.
# ---------------------------------------------------------------------------
mm <- msigdbr::msigdbr(species = "Mus musculus")
if (!all(c("gs_name", "db_gene_symbol", "gene_symbol") %in% names(mm))) {
  stop("Mouse msigdbr output must contain columns 'gs_name', 'db_gene_symbol', and 'gene_symbol'.")
}
mm_hsr <- mm[mm$gs_name %in% TARGET_SETS, c("db_gene_symbol", "gene_symbol"), drop = FALSE]
pair <- unique(mm_hsr)
names(pair) <- c("human_symbol", "mouse_symbol")
pair <- pair[pair$human_symbol %in% human_union &
               !is.na(pair$human_symbol) & pair$human_symbol != "" &
               !is.na(pair$mouse_symbol) & pair$mouse_symbol != "", , drop = FALSE]
pair <- pair[order(pair$human_symbol, pair$mouse_symbol), , drop = FALSE]
rownames(pair) <- NULL

mouse_union <- sort(unique(pair$mouse_symbol))
if (length(mouse_union) != 170L) {
  detail <- sprintf(
    "observed=%d expected=170; missing sanity genes: %s",
    length(mouse_union),
    paste(setdiff(c("Hspa1a", "Hspa1b", "Gml", "Gml2"), mouse_union), collapse = ", ")
  )
  stop(sprintf("Mouse msigdbr-native HSR union size drift detected: %s", detail))
}
required_mouse <- c("Hspa1a", "Hspa1b", "Gml", "Gml2")
missing_required_mouse <- setdiff(required_mouse, mouse_union)
if (length(missing_required_mouse)) {
  stop(sprintf(
    "Mouse msigdbr-native HSR union is missing required paralog sanity gene(s): %s",
    paste(missing_required_mouse, collapse = ", ")
  ))
}

# ---------------------------------------------------------------------------
# 4) Validate mapped mouse symbols against THIS dataset's gene_name column.
# ---------------------------------------------------------------------------
cpm <- read.csv(FILE_CPM, check.names = FALSE)
if (!"gene_name" %in% names(cpm)) {
  stop(sprintf("FILE_CPM must contain a 'gene_name' column: %s", FILE_CPM))
}
dataset_symbols <- unique(cpm$gene_name)
pair$in_GSE329522 <- pair$mouse_symbol %in% dataset_symbols
mouse_core <- sort(unique(pair$mouse_symbol[pair$in_GSE329522]))

# ---------------------------------------------------------------------------
# 5) Full per-gene audit log, one row per HUMAN union gene.
# ---------------------------------------------------------------------------
pair_split <- split(pair, pair$human_symbol)
orth_summary <- data.frame(
  human_symbol = human_union,
  n_mouse_orthologs = vapply(human_union, function(x) {
    if (!x %in% names(pair_split)) {
      return(0L)
    }
    length(unique(pair_split[[x]]$mouse_symbol))
  }, integer(1)),
  mouse_symbols = vapply(human_union, function(x) {
    if (!x %in% names(pair_split)) {
      return("")
    }
    paste(sort(unique(pair_split[[x]]$mouse_symbol)), collapse = ";")
  }, character(1)),
  any_mouse_in_GSE329522 = vapply(human_union, function(x) {
    if (!x %in% names(pair_split)) {
      return(FALSE)
    }
    any(pair_split[[x]]$in_GSE329522)
  }, logical(1)),
  stringsAsFactors = FALSE
)
log_df <- merge(membership, orth_summary, by = "human_symbol", all.x = TRUE, sort = FALSE)
log_df <- log_df[order(-log_df$n_sets, log_df$human_symbol), c(
  "human_symbol",
  "in_reactome_crhs",
  "in_reactome_hsf1",
  "in_gobp_response_to_heat",
  "n_sets",
  "n_mouse_orthologs",
  "mouse_symbols",
  "any_mouse_in_GSE329522"
)]

# ---------------------------------------------------------------------------
# 6) Tidy set-size table, including staged candidates deliberately excluded.
# ---------------------------------------------------------------------------
setsizes <- data.frame(
  set = c(TARGET_SETS, "HSR_thermal_core_union", EXCLUDED_SETS),
  role = c(rep("core_component", length(TARGET_SETS)), "union", rep("excluded", length(EXCLUDED_SETS))),
  n_human = c(unname(observed_n[TARGET_SETS]), length(human_union), rep(NA_integer_, length(EXCLUDED_SETS))),
  note = c(
    "Reactome cellular heat-stress response; included as HSR-core component; mouse basis is msigdbr-native.",
    "Reactome HSF1-mediated heat-shock-response regulation; included as HSR-core component; mouse basis is msigdbr-native.",
    "GO biological process response to heat; included as HSR-core component; mouse basis is msigdbr-native.",
    "Union of the three included HSR-core source sets; mouse basis is msigdbr-native, paralog-complete.",
    "Thermosensory-neuron detection program; excluded as irrelevant to T-cell transcriptional heat response.",
    "Thermoception-specific thermosensory-neuron detection program; excluded as irrelevant to T cells.",
    "Human Phenotype Ontology febrile-syndrome etiology panel; excluded because mutation etiology is not a transcriptional heat-response signature."
  ),
  stringsAsFactors = FALSE
)

# ---------------------------------------------------------------------------
# 7) Save assets with provenance attributes.
# ---------------------------------------------------------------------------
component_sets <- data.frame(
  gs_name = TARGET_SETS,
  n_human = unname(observed_n[TARGET_SETS]),
  stringsAsFactors = FALSE
)
provenance <- list(
  signature = "Curated HSR core (HSF1 + co-chaperone categories)",
  citation = "MSigDB v2026.1.Hs (Reactome, GO:BP)",
  source = "msigdbr",
  msigdbr_version = as.character(utils::packageVersion("msigdbr")),
  component_sets = component_sets,
  union_rule = "union of the three HSR-core sets",
  n_human_union = length(human_union),
  mouse_basis = "msigdbr species=Mus musculus (source-native ortholog conversion; paralog-complete: HSPA1A/HSPA1B->Hspa1a/Hspa1b, GML->Gml/Gml2)",
  n_mouse_union = length(mouse_union),
  n_mouse_in_dataset = length(mouse_core),
  built_by = "02_analysis/scripts/00d_curate_temp_hsr.R",
  note = HONEST_CEILING
)

human_set <- structure(list(HSR_thermal_core = human_union), provenance = provenance)
mouse_set <- structure(list(HSR_thermal_core = mouse_core), provenance = provenance)
mouse_sensitivity_set <- structure(list(HSR_thermal_core = mouse_union), provenance = provenance)

saveRDS(human_set, file.path(OUT_DIR, "temp_hsr_human.rds"))
saveRDS(mouse_set, file.path(OUT_DIR, "temp_hsr_mouse.rds"))
saveRDS(mouse_sensitivity_set, file.path(OUT_DIR, "temp_hsr_mouse_sensitivity.rds"))
write.csv(log_df, file.path(OUT_DIR, "temp_hsr_curation_log.csv"), row.names = FALSE)
write.csv(pair, file.path(OUT_DIR, "temp_hsr_ortholog_map.csv"), row.names = FALSE)
write.csv(setsizes, file.path(TBL_DIR, "temp_hsr_setsizes.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 8) Report.
# ---------------------------------------------------------------------------
cat("\n================ Curated HSR thermal proteostasis lens ================\n")
cat(sprintf("msigdbr version:           %s\n", as.character(utils::packageVersion("msigdbr"))))
cat("Source sets:\n")
for (set_name in TARGET_SETS) {
  cat(sprintf("  %-64s %4d genes\n", set_name, observed_n[[set_name]]))
}
cat(sprintf("Human union:               %d genes\n", length(human_union)))
cat(sprintf("Mouse union:               %d genes\n", length(mouse_union)))
cat(sprintf("In GSE329522 (mouse core): %d genes\n", length(mouse_core)))
cat(sprintf("Hspa1b present:            %s\n", "Hspa1b" %in% mouse_union))
cat(sprintf("Gml2 present:              %s\n", "Gml2" %in% mouse_union))
cat(sprintf("Honest ceiling:            %s\n", HONEST_CEILING))
cat("\nWrote:\n")
cat("  ", file.path(OUT_DIR, "temp_hsr_human.rds"), "\n")
cat("  ", file.path(OUT_DIR, "temp_hsr_mouse.rds"), "\n")
cat("  ", file.path(OUT_DIR, "temp_hsr_mouse_sensitivity.rds"), "\n")
cat("  ", file.path(OUT_DIR, "temp_hsr_curation_log.csv"), "\n")
cat("  ", file.path(OUT_DIR, "temp_hsr_ortholog_map.csv"), "\n")
cat("  ", file.path(TBL_DIR, "temp_hsr_setsizes.csv"), "\n")
cat("======================================================================\n")
