#!/usr/bin/env Rscript
# 00e_curate_temp_hsr_lens.R - Freeze the taxonomy-refined HSR lens
# Project: STING-cGAS-GSE329522
# Phase: 0 (infrastructure / curated reference asset)
# Inputs:
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_gene_taxonomy.csv   (frozen; agy + curator)
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_ortholog_map.csv    (authoritative 1:many human<->mouse)
# Outputs:
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_human_lens.rds
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_mouse_lens.rds
#   00_data/references/gene_sets/temp_hsr_lens/tables/temp_hsr_taxonomy_summary.csv
# Dependencies: 02_analysis/config/config.R
#
# ===========================================================================================
# METHOD & RATIONALE (decisions are auditable, not arbitrary)
# ===========================================================================================
#   (A) CLEANED CYTOSOLIC CORE
#       The scored downstream lens is the frozen taxonomy rows where cytosolic_core == TRUE
#       (categories {hsf1_core_hsr, co_chaperone}). npc_transport / upr_er / thermosensory /
#       generic_stress are bleed or separate-axis and are excluded from the scored core.
#
#   (B) SENSITIVITY TIER
#       The sensitivity tier is the full 176-gene union, carried for robustness only.
#
#   (C) FROZEN INPUTS ONLY; MOUSE VIA THE AUTHORITATIVE ORTHOLOG MAP
#       No MSigDB re-pull is performed. Mouse symbols are joined from
#       temp_hsr_ortholog_map.csv (msigdbr species="Mus musculus" native conversion, 1:many,
#       paralog-complete: HSPA1A/HSPA1B -> Hspa1a/Hspa1b, GML -> Gml/Gml2). A human gene may
#       map to >1 mouse symbol, so mouse set sizes are counts of UNIQUE mouse symbols, not of
#       human rows.
#
#   (D) HONEST CEILING
#       Even the cleaned cytosolic core is proteotoxic-stress-general, not fever-specific;
#       only the mouse 37/39 contrast can measure thermal-ness.
#
#   DETERMINISM: no RNG is used; all symbols are sorted before writing.
# ===========================================================================================

source("02_analysis/config/config.R")

OUT_DIR <- file.path(PROJECT_ROOT, "00_data/references/gene_sets/temp_hsr_lens")
TBL_DIR <- file.path(OUT_DIR, "tables")
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)

TAXONOMY_FILE     <- file.path(OUT_DIR, "temp_hsr_gene_taxonomy.csv")
ORTHOLOG_MAP_FILE <- file.path(OUT_DIR, "temp_hsr_ortholog_map.csv")

EXPECTED_CATEGORIES <- c(
  "co_chaperone",
  "generic_stress",
  "hsf1_core_hsr",
  "npc_transport",
  "thermosensory",
  "upr_er"
)

HONEST_CEILING <- paste(
  "Even the cleaned cytosolic core is proteotoxic-stress-general, not fever-specific;",
  "HSF1 also fires on oxidative, proteasome, and heavy-metal stress, so only the",
  "mouse 37/39 contrast can measure thermal-ness."
)

stop_mismatch <- function(label, observed, expected) {
  if (!identical(as.integer(observed), as.integer(expected))) {
    stop(sprintf("%s size mismatch: observed=%d expected=%d", label, observed, expected))
  }
}

required_taxonomy_cols <- c("human_symbol", "category", "cytosolic_core", "n_sets",
                            "rationale", "source")
required_map_cols      <- c("human_symbol", "mouse_symbol", "in_GSE329522")

taxonomy     <- read.csv(TAXONOMY_FILE, stringsAsFactors = FALSE, check.names = FALSE)
ortholog_map <- read.csv(ORTHOLOG_MAP_FILE, stringsAsFactors = FALSE, check.names = FALSE)

for (chk in list(list("Taxonomy", required_taxonomy_cols, names(taxonomy)),
                 list("Ortholog map", required_map_cols, names(ortholog_map)))) {
  miss <- setdiff(chk[[2]], chk[[3]])
  if (length(miss)) stop(sprintf("%s is missing required column(s): %s",
                                 chk[[1]], paste(miss, collapse = ", ")))
}

stop_mismatch("Taxonomy row count", nrow(taxonomy), 176L)

observed_categories <- sort(unique(taxonomy$category))
if (!identical(observed_categories, EXPECTED_CATEGORIES)) {
  stop(sprintf("Taxonomy categories mismatch: observed={%s}; expected={%s}",
               paste(observed_categories, collapse = ", "),
               paste(EXPECTED_CATEGORIES, collapse = ", ")))
}

if (anyDuplicated(taxonomy$human_symbol)) stop("Taxonomy contains duplicated human_symbol values.")

taxonomy$cytosolic_core   <- as.logical(taxonomy$cytosolic_core)
ortholog_map$in_GSE329522 <- as.logical(ortholog_map$in_GSE329522)

# Every ortholog-map human gene must be in the union; every mouse-mapped taxonomy gene must
# appear in the map. (Some human genes have no mouse ortholog and simply never appear in the map.)
if (!all(ortholog_map$human_symbol %in% taxonomy$human_symbol)) {
  stop("Ortholog map contains human symbols absent from the taxonomy union.")
}

stop_mismatch("Taxonomy cytosolic_core TRUE count", sum(taxonomy$cytosolic_core, na.rm = TRUE), 56L)

# ---------------------------------------------------------------------------
# Derive the four gene vectors. Human from the taxonomy; mouse from the join to the
# authoritative 1:many ortholog map (unique mouse symbols).
# ---------------------------------------------------------------------------
core_human_symbols <- taxonomy$human_symbol[taxonomy$cytosolic_core]

human_core <- sort(unique(core_human_symbols))
human_sens <- sort(unique(taxonomy$human_symbol))

map_core <- ortholog_map[ortholog_map$human_symbol %in% core_human_symbols, , drop = FALSE]
mouse_core <- sort(unique(map_core$mouse_symbol[map_core$in_GSE329522]))   # expressed cytosolic core
mouse_sens <- sort(unique(ortholog_map$mouse_symbol))                      # full mapped union

stop_mismatch("human_core", length(human_core), 56L)
stop_mismatch("human_sens", length(human_sens), 176L)
# 56 cytosolic-core human genes -> 47 unique expressed mouse symbols. The msigdbr-native map
# keeps both HSPA1A->Hspa1a and HSPA1B->Hspa1a/Hspa1b, so the inducible paralog Hspa1b is
# retained in the mouse core (the babelgene 1:1 collapse had dropped it).
stop_mismatch("mouse_core", length(mouse_core), 47L)
stop_mismatch("mouse_sens", length(mouse_sens), 170L)

provenance <- list(
  signature        = "curated HSR lens (taxonomy-refined)",
  taxonomy_source  = "temp_hsr_gene_taxonomy.csv (agy_gemini_3.1_pro_2026-07-14 + curator npc_transport split)",
  ortholog_source  = "temp_hsr_ortholog_map.csv (msigdbr species=Mus musculus; 1:many, paralog-complete)",
  core_rule        = "categories {hsf1_core_hsr, co_chaperone}; cytosolic_core==TRUE",
  sensitivity_rule = "full 176-gene union",
  n_core_human     = 56L,
  n_sens_human     = 176L,
  n_core_mouse     = 47L,
  n_sens_mouse     = 170L,
  built_by         = "02_analysis/scripts/00e_curate_temp_hsr_lens.R",
  note             = HONEST_CEILING
)

human_lens <- structure(
  list(HSR_core        = structure(human_core, provenance = provenance),
       HSR_sensitivity = structure(human_sens, provenance = provenance)),
  provenance = provenance
)
mouse_lens <- structure(
  list(HSR_core        = structure(mouse_core, provenance = provenance),
       HSR_sensitivity = structure(mouse_sens, provenance = provenance)),
  provenance = provenance
)

saveRDS(human_lens, file.path(OUT_DIR, "temp_hsr_human_lens.rds"))
saveRDS(mouse_lens, file.path(OUT_DIR, "temp_hsr_mouse_lens.rds"))

# ---------------------------------------------------------------------------
# Per-category summary + the two derived-term rows. Mouse counts are unique mouse symbols
# from the ortholog-map join, so paralog expansion is reflected honestly.
# ---------------------------------------------------------------------------
mouse_for <- function(human_syms, expressed_only = FALSE) {
  m <- ortholog_map[ortholog_map$human_symbol %in% human_syms, , drop = FALSE]
  if (expressed_only) m <- m[m$in_GSE329522, , drop = FALSE]
  length(unique(m$mouse_symbol))
}

category_rows <- do.call(rbind, lapply(EXPECTED_CATEGORIES, function(category_name) {
  hs <- taxonomy$human_symbol[taxonomy$category == category_name]
  data.frame(
    category             = category_name,
    n_human              = length(unique(hs)),
    n_mouse_mapped       = mouse_for(hs, expressed_only = FALSE),
    n_mouse_in_GSE329522 = mouse_for(hs, expressed_only = TRUE),
    in_cytosolic_core    = all(taxonomy$cytosolic_core[taxonomy$category == category_name]),
    stringsAsFactors     = FALSE
  )
}))

summary_rows <- data.frame(
  category             = c("HSR_core", "HSR_sensitivity"),
  n_human              = c(length(human_core), length(human_sens)),
  n_mouse_mapped       = c(mouse_for(core_human_symbols, FALSE), length(mouse_sens)),
  n_mouse_in_GSE329522 = c(length(mouse_core), mouse_for(taxonomy$human_symbol, TRUE)),
  in_cytosolic_core    = c(TRUE, FALSE),
  stringsAsFactors     = FALSE
)

taxonomy_summary <- rbind(category_rows, summary_rows)
write.csv(taxonomy_summary, file.path(TBL_DIR, "temp_hsr_taxonomy_summary.csv"), row.names = FALSE)

cat("\n================ Taxonomy-refined HSR lens ================\n")
cat(sprintf("Human HSR_core:        %d genes\n", length(human_core)))
cat(sprintf("Human HSR_sensitivity: %d genes\n", length(human_sens)))
cat(sprintf("Mouse HSR_core:        %d genes (unique mouse symbols; incl. Hspa1b)\n", length(mouse_core)))
cat(sprintf("Mouse HSR_sensitivity: %d genes\n", length(mouse_sens)))
cat("\nPer-category counts:\n")
for (i in seq_len(nrow(category_rows))) {
  cat(sprintf("  %-16s human=%3d mouse_mapped=%3d mouse_in_GSE329522=%3d cytosolic_core=%s\n",
              category_rows$category[[i]], category_rows$n_human[[i]],
              category_rows$n_mouse_mapped[[i]], category_rows$n_mouse_in_GSE329522[[i]],
              category_rows$in_cytosolic_core[[i]]))
}
cat(sprintf("\nHonest ceiling: %s\n", HONEST_CEILING))
cat("\nWrote:\n")
cat("  ", file.path(OUT_DIR, "temp_hsr_human_lens.rds"), "\n")
cat("  ", file.path(OUT_DIR, "temp_hsr_mouse_lens.rds"), "\n")
cat("  ", file.path(TBL_DIR, "temp_hsr_taxonomy_summary.csv"), "\n")
cat("===========================================================\n")
