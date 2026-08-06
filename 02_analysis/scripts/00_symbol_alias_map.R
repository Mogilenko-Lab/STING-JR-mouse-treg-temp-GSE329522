#!/usr/bin/env Rscript
# 00_symbol_alias_map.R — ASSET BUILDER. The committed reference-to-matrix symbol map.
# =============================================================================
# The collaborators' delivered CPM table was quantified against GRCm38 + GENCODE vM25, so
# this project's vocabulary is frozen to that build's MGI vintage: 2,341 of the 19,679
# modelled symbols have since moved off the official MGI list. Reference gene sets ship
# current symbols and every match here is an exact string match, so a renamed gene leaves a
# set silently and the loss reads as biological absence. The worst case is stark: the matrix
# carries the whole legacy ATP-synthase block Atp5a1 Atp5b Atp5c1 Atp5d Atp5e Atp5g1..g3
# Atp5h Atp5j Atp5j2 Atp5k Atp5l Atp5o Atpif1 while MitoCarta 3.0 ships Atp5f1a..Atp5po, so
# MITOPATHWAYS_OXPHOS.Complex_V matched 7 of its 22 genes, fell under gsea_min_size = 15,
# and had no row in master_gsea_table.csv at all. This script resolves that once and writes
# the answer down.
#
# WHY A COMMITTED CSV. The map is a property of (matrix vocabulary x org.Mm.eg.db release)
# and of nothing else, so it is computed once here. Hardcoding the pairs in config would go
# stale without a diff, and resolving live in every consumer would load the annotation
# database repeatedly and leave nothing to inspect. A committed CSV is auditable and moves
# when someone regenerates it on purpose.
#
# WHAT THIS SCRIPT COVERS. The reference side alone: reference symbol DOWN into the matrix
# vocabulary. The ortholog query side runs the other way — matrix symbol UP to current, which
# is what babelgene keys on — and lives in helpers/ortholog_utils.R::normalise_mouse_query(),
# applied inside build_ortholog_map() at 18_projection_export.R and published there as its
# own ledger. The one-to-one safety condition means something different in each direction, so
# the two stay separate.
#
# Reads, read-only:
#   02_analysis/config/analysis_config.yaml            symbol_alias + databases
#   03_results/objects/gene_universe.txt               the vocabulary resolved INTO
#   the 8 MSigDB collections and the 5 custom databases named under `databases:`
#
# Writes (committed assets, NOT stage results):
#   00_data/references/symbol_alias/symbol_alias_map.csv
#   00_data/references/symbol_alias/symbol_alias_provenance.csv
#
# Run from the project root:
#   Rscript 02_analysis/scripts/00_symbol_alias_map.R

source("02_analysis/config/config.R")             # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS
source("02_analysis/helpers/de_gsea_helpers.R")   # load_msigdb_collection, load_custom_geneset
source("02_analysis/helpers/symbol_alias.R")      # build_alias_map, ALIAS_RESOLUTIONS
source("02_analysis/helpers/source_hash_manifest.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})
options(stringsAsFactors = FALSE)

SA <- YAML_CONFIG$symbol_alias
if (is.null(SA))
  stop("[00_alias] analysis_config.yaml has no `symbol_alias:` block — add it before running.")
if (is.null(YAML_CONFIG$databases))
  stop("[00_alias] analysis_config.yaml has no `databases:` block.")

MAP_OUT   <- file.path(PROJECT_ROOT, SA$map_path)
PROV_OUT  <- file.path(PROJECT_ROOT, SA$provenance_path)
VOCAB_TXT <- file.path(PROJECT_ROOT, SA$matrix_vocabulary)
FLAGGED   <- unlist(SA$flagged_for_review) %||% character()
dir.create(dirname(MAP_OUT), recursive = TRUE, showWarnings = FALSE)

message("=================================================================")
message("00_symbol_alias_map — reference MGI symbols -> this matrix's vocabulary")
message("=================================================================")

# ============================================================================
# 1. The vocabulary the map resolves INTO
# ============================================================================
# gene_universe.txt is the modelled DE background, and in this project it is also exactly
# the set of symbols the delivered CPM table carries: the duplicate-symbol collapse at
# 01_mapping_qc.R drops Ensembl ids and keeps every symbol, so 19,679 unique gene_name values
# survive as 19,679 modelled symbols. That single layer is worth stating, because it means
# an unmatched reference gene here is either a vocabulary result or genuinely not in the
# delivered quantification — there is no expression filter in between to blame.

if (!file.exists(VOCAB_TXT))
  stop("[00_alias] matrix vocabulary not found at ", VOCAB_TXT,
       " — run 11_emit_universe.R first.")
MATRIX_SYMBOLS <- unique(trimws(readLines(VOCAB_TXT, warn = FALSE)))
MATRIX_SYMBOLS <- MATRIX_SYMBOLS[nzchar(MATRIX_SYMBOLS)]
message(sprintf("[1] matrix vocabulary: %d symbols from %s",
                length(MATRIX_SYMBOLS), SA$matrix_vocabulary))
for (p in c("Atp5a1", "Atp5f1a", "Cgas", "Mb21d1", "Sting1", "Ddx58", "Rigi"))
  message(sprintf("    %-8s in matrix vocabulary: %s", p, p %in% MATRIX_SYMBOLS))

# ============================================================================
# 2. Every reference symbol this project consumes on the gene-set side
# ============================================================================
# The union is deliberately wider than any single collection needs, because the map has to
# answer for every consumer without being rebuilt per collection.

sources <- list()
note <- function(label, syms) {
  syms <- unique(syms[!is.na(syms) & nzchar(syms)])
  sources[[label]] <<- syms
  message(sprintf("    %-28s %6d unique symbols", label, length(syms)))
  invisible(syms)
}

message("[2] collecting reference symbols ...")
for (m in YAML_CONFIG$databases$msigdb)
  note(paste0("msigdb:", m$name),
       unlist(load_msigdb_collection(category = m$category,
                                     subcategory = m$subcategory %||% "",
                                     species = SPECIES), use.names = FALSE))
for (db in YAML_CONFIG$databases$custom) {
  p <- file.path(PROJECT_ROOT, db$path)
  if (!file.exists(p)) { message("    missing custom DB: ", db$path); next }
  note(paste0("custom:", db$name), unlist(load_custom_geneset(p), use.names = FALSE))
}

REFERENCE_SYMBOLS <- sort(unique(unlist(sources, use.names = FALSE)))
message(sprintf("[2] reference universe: %d unique symbols over %d sources",
                length(REFERENCE_SYMBOLS), length(sources)))

# ============================================================================
# 3. Build the map
# ============================================================================

message("[3] resolving against org.Mm.eg.db ",
        as.character(packageVersion("org.Mm.eg.db")), " ...")
MAP <- build_alias_map(REFERENCE_SYMBOLS, MATRIX_SYMBOLS,
                       db = org.Mm.eg.db::org.Mm.eg.db, flagged_pairs = FLAGGED)
SUMM <- attr(MAP, "summary")

# Every flagged pair must be a pair the guard would otherwise have accepted, or the
# exclusion list is stale and reads as a decision that is no longer doing anything.
missed_flags <- setdiff(FLAGGED, paste0(MAP$reference_symbol, "->", MAP$matrix_symbol))
if (length(missed_flags))
  stop(sprintf(paste0("[00_alias] flagged_for_review names %d pair(s) this build does not ",
                      "produce as a candidate at all (%s). Either the pair is stale or the ",
                      "reference universe no longer contains it — decide, do not leave it ",
                      "silently inert."),
               length(missed_flags), paste(missed_flags, collapse = ", ")))

MAP <- MAP %>% arrange(match(.data$resolution, ALIAS_RESOLUTIONS), .data$reference_symbol)
readr::write_csv(MAP, MAP_OUT)
message(sprintf("  [SAVE] %s  %d candidate pairs", SA$map_path, nrow(MAP)))
print(as.data.frame(MAP %>% count(resolution, name = "n_pairs")), row.names = FALSE)

message("  the ATP-synthase block, the reason this asset exists:")
print(as.data.frame(MAP %>% filter(grepl("^Atp5", .data$reference_symbol)) %>%
                      transmute(pair = paste0(.data$reference_symbol, "->",
                                              .data$matrix_symbol), .data$resolution)),
      row.names = FALSE)
message("  every pair the guards withheld:")
print(as.data.frame(MAP %>% filter(.data$resolution != "accepted") %>%
                      transmute(pair = paste0(.data$reference_symbol, "->",
                                              .data$matrix_symbol), .data$resolution)),
      row.names = FALSE)

# ============================================================================
# 4. Provenance — what the map was built from, so a stale map is visible
# ============================================================================
# The no-candidate counts live here, outside the map's rows. A reference symbol with
# no alias in this vocabulary is not a decision about a pair, and carrying thousands of such
# rows would defeat the one review step this asset most needs: reading the map cold and
# recognising every accepted pair as a nomenclature update.

PROV <- tibble(
  built_by = "02_analysis/scripts/00_symbol_alias_map.R",
  org_mm_eg_db_version = as.character(packageVersion("org.Mm.eg.db")),
  annotationdbi_version = as.character(packageVersion("AnnotationDbi")),
  msigdbr_version = as.character(packageVersion("msigdbr")),
  species = SPECIES,
  matrix_vocabulary_path = SA$matrix_vocabulary,
  matrix_vocabulary_sha256 = source_sha256(VOCAB_TXT),
  n_matrix_symbols = length(MATRIX_SYMBOLS),
  n_reference_symbols = SUMM$n_reference_symbols,
  n_reference_sources = length(sources),
  reference_sources = paste(names(sources), collapse = "; "),
  n_reference_absent_from_vocabulary = SUMM$n_absent_from_vocabulary,
  n_reference_not_in_org_db = SUMM$n_not_in_org_db,
  n_reference_no_alias_in_vocabulary = SUMM$n_no_alias_in_vocabulary,
  n_candidate_pairs = SUMM$n_candidates,
  n_accepted = SUMM$n_accepted,
  n_flagged_for_review = SUMM$n_flagged,
  n_rejected = SUMM$n_rejected,
  flagged_for_review = paste(FLAGGED, collapse = "; "))
readr::write_csv(PROV, PROV_OUT)
message(sprintf("  [SAVE] %s", SA$provenance_path))
print(as.data.frame(t(PROV)))

# The map may only ever ADD a symbol to a set, so nothing it accepts may already be a
# reference symbol that matched, and every accepted target must be in the vocabulary.
stopifnot(
  "an accepted target must be in the matrix vocabulary" =
    all(MAP$matrix_symbol[MAP$resolution == "accepted"] %in% MATRIX_SYMBOLS),
  "an accepted reference symbol must be ABSENT from the matrix vocabulary" =
    !any(MAP$reference_symbol[MAP$resolution == "accepted"] %in% MATRIX_SYMBOLS),
  "one reference symbol resolves to at most one matrix symbol" =
    !any(duplicated(MAP$reference_symbol[MAP$resolution == "accepted"])),
  "the ATP-synthase block must resolve" =
    all(c("Atp5f1a", "Atp5f1b", "Atp5po", "Atp5if1") %in%
          MAP$reference_symbol[MAP$resolution == "accepted"]))

message("\n[DONE] symbol alias map built. 04_gsea_set_prep.R consumes it.")
