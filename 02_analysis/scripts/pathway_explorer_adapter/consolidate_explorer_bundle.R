#!/usr/bin/env Rscript
# consolidate_explorer_bundle.R — validate and build the pathway-explorer input bundle
# Run from project root.

source("02_analysis/config/config.R")          # CONFIG, proj_path, cfg_path
if (file.exists("02_analysis/helpers/de_gsea_helpers.R")) {
  source("02_analysis/helpers/de_gsea_helpers.R")
}
`%||%` <- if (exists("%||%")) `%||%` else function(a, b) if (is.null(a) || length(a) == 0) b else a
if (!exists("round_numeric_cols")) round_numeric_cols <- function(df, sig = 9L) {
  for (nm in names(df)) { col <- df[[nm]]; if (is.double(col)) df[[nm]] <- signif(col, sig) }; df
}
suppressPackageStartupMessages({ library(readr); library(dplyr) })
options(stringsAsFactors = FALSE)

PID      <- YAML_CONFIG$project$id
MST_DIR  <- Sys.getenv("EXPLORER_MASTER_DIR", unset = "")
if (!nzchar(MST_DIR)) MST_DIR <- DIR_MASTER
MST_DIR  <- sub("/+$", "", MST_DIR)

UNIFIED  <- file.path(MST_DIR, "master_unified.csv")
DE_MST   <- file.path(MST_DIR, "master_de_table.csv")
ATLAS_FP <- file.path(MST_DIR, "atlas_gene_universe.txt")
MANIFEST <- file.path(MST_DIR, "explorer_manifest.json")

# ===========================================================================
# 0. Copy and Reconstruct Files (STING-specific step)
# ===========================================================================
message("[adapter] Consolidating inputs in master directory ...")

# 0a. Copy gene universe list to atlas_gene_universe.txt
univ_obj <- file.path(DIR_OBJECTS, "gene_universe.txt")
if (file.exists(univ_obj)) {
  file.copy(univ_obj, ATLAS_FP, overwrite = TRUE)
  message("[adapter] Copied gene_universe.txt to ", ATLAS_FP)
} else {
  stop("[adapter] gene_universe.txt missing from objects dir — cannot build atlas universe sidecar.")
}

# 0b. Copy master_de_genes.csv to master_de_table.csv
de_genes_fp <- file.path(MST_DIR, "master_de_genes.csv")
if (file.exists(de_genes_fp)) {
  file.copy(de_genes_fp, DE_MST, overwrite = TRUE)
  message("[adapter] Copied master_de_genes.csv to ", DE_MST)
} else {
  warning("[adapter] master_de_genes.csv missing from master dir — running-sum panel will be empty.")
}

# 0c. Build master_unified.csv
gsea_mst_fp <- file.path(MST_DIR, "master_gsea_table.csv")
if (!file.exists(gsea_mst_fp)) {
  stop("[adapter] master_gsea_table.csv is missing at ", gsea_mst_fp)
}
gsea_mst <- readr::read_csv(gsea_mst_fp, show_col_types = FALSE)

# Add entity_type
gsea_mst$entity_type <- "Pathway"

# Intersect with atlas to populate genes_full_set
atlas_universe <- readLines(ATLAS_FP)
atlas_universe <- unique(atlas_universe[!is.na(atlas_universe) & nzchar(atlas_universe)])

obj_dir <- DIR_OBJECTS
geneset_files <- list.files(obj_dir, pattern = "^(geneset_.*|coresh_derived_sets)\\.rds$", full.names = TRUE)
all_gene_sets <- list()
for (f in geneset_files) {
  obj <- readRDS(f)
  if (is.list(obj)) {
    if (all(c("T2G") %in% names(obj)) && is.data.frame(obj$T2G)) {
      sets <- split(as.character(obj$T2G$gene_symbol), as.character(obj$T2G$gs_name))
      all_gene_sets <- c(all_gene_sets, sets)
    } else {
      names_obj <- names(obj)
      if (!is.null(names_obj)) {
        sets <- lapply(obj, function(g) unique(as.character(g)))
        all_gene_sets <- c(all_gene_sets, sets)
      }
    }
  }
}

gsea_mst$genes_full_set <- vapply(
  seq_len(nrow(gsea_mst)),
  function(i) {
    pid <- gsea_mst$pathway_id[i]
    members <- all_gene_sets[[pid]]
    if (is.null(members)) {
      # Fallback to core_enrichment if not found in gene sets
      return(gsea_mst$core_enrichment[i])
    }
    paste(intersect(members, atlas_universe), collapse = "/")
  },
  character(1)
)

# Load TF activities
tf_mst_fp <- file.path(MST_DIR, "master_tf_activities.csv")
tf_mst <- if (file.exists(tf_mst_fp)) {
  readr::read_csv(tf_mst_fp, show_col_types = FALSE) %>%
    dplyr::mutate(entity_type = "TF", genes_full_set = core_enrichment)
} else {
  NULL
}

# Load PROGENy activities
progeny_mst_fp <- file.path(MST_DIR, "master_progeny_activities.csv")
progeny_mst <- if (file.exists(progeny_mst_fp)) {
  readr::read_csv(progeny_mst_fp, show_col_types = FALSE) %>%
    dplyr::mutate(entity_type = "PROGENy", genes_full_set = core_enrichment)
} else {
  NULL
}

# Load GATOM modules
gatom_mst_fp <- file.path(MST_DIR, "master_gatom_modules.csv")
gatom_mst <- if (file.exists(gatom_mst_fp) && file.size(gatom_mst_fp) > 10L) {
  gatom_df <- readr::read_csv(gatom_mst_fp, show_col_types = FALSE)
  if (nrow(gatom_df) > 0) {
    gatom_gfs <- gatom_df %>%
      dplyr::group_by(pathway_id) %>%
      dplyr::summarise(
        genes_full_set = paste(unique(unlist(strsplit(core_enrichment, "/", fixed = TRUE))), collapse = "/"),
        .groups = "drop"
      )
    gatom_df %>%
      dplyr::left_join(gatom_gfs, by = "pathway_id") %>%
      dplyr::mutate(entity_type = "GATOM")
  } else {
    NULL
  }
} else {
  NULL
}

# Stack all together
unified_list <- list(gsea_mst)
if (!is.null(tf_mst)) unified_list[[length(unified_list) + 1]] <- tf_mst
if (!is.null(progeny_mst)) unified_list[[length(unified_list) + 1]] <- progeny_mst
if (!is.null(gatom_mst)) unified_list[[length(unified_list) + 1]] <- gatom_mst

unified_df <- dplyr::bind_rows(unified_list)
readr::write_csv(unified_df, UNIFIED)
message("[adapter] Built master_unified.csv containing GSEA, TF, PROGENy, and GATOM modules.")

# ===========================================================================
# 1. Atlas sidecar — HARD REQUIREMENT
# ===========================================================================
if (!file.exists(ATLAS_FP))
  stop("[adapter] REQUIRED: atlas_gene_universe.txt absent at ", ATLAS_FP)

atlas_universe <- readLines(ATLAS_FP)
atlas_n        <- length(atlas_universe)
message(sprintf("[adapter] Atlas loaded: %d symbols from %s.", atlas_n, ATLAS_FP))

# ===========================================================================
# 2. master_unified.csv — validate and apply idempotent schema repair only
# ===========================================================================
REQ_UNI <- c("entity_type","pathway_id","pathway_name","database","nes","padj",
             "set_size","core_enrichment","genes_full_set","contrast")

if (!file.exists(UNIFIED)) {
  stop("[adapter] master_unified.csv absent (", UNIFIED, ")")
}

u <- readr::read_csv(UNIFIED, show_col_types = FALSE)
changed <- FALSE

# Idempotent schema repair: NES -> nes
if ("NES" %in% names(u) && !"nes" %in% names(u)) { u <- dplyr::rename(u, nes = NES); changed <- TRUE }
if ("NES" %in% names(u) &&  "nes" %in% names(u)) { u <- dplyr::select(u, -NES);       changed <- TRUE }

# Idempotent schema repair: adj.P.Val -> padj
if (!"padj" %in% names(u) && "adj.P.Val" %in% names(u)) { u <- dplyr::rename(u, padj = adj.P.Val); changed <- TRUE }

if (!"entity_type" %in% names(u))
  stop("[adapter] master_unified.csv lacks 'entity_type'")

if (!"genes_full_set" %in% names(u))
  stop("[adapter] genes_full_set column absent from master_unified.csv")

empty_gfs <- is.na(u$genes_full_set) | u$genes_full_set == ""
if (any(empty_gfs)) {
  # Clean up rows with empty genes_full_set or report error
  u <- u[!empty_gfs, ]
  changed <- TRUE
  message(sprintf("[adapter] Cleaned up %d rows with empty genes_full_set", sum(empty_gfs)))
}

miss <- setdiff(REQ_UNI, names(u))
if (length(miss)) stop("[adapter] master_unified.csv missing required columns: ", paste(miss, collapse = ", "))

# ===========================================================================
# 3. Contrast-invariance assertion
# ===========================================================================
u$gfs_card <- vapply(u$genes_full_set, function(x) {
  if (is.na(x) || x == "") NA_integer_ else length(strsplit(x, "/", fixed = TRUE)[[1]])
}, integer(1))

drift_check <- u %>%
  dplyr::group_by(pathway_id) %>%
  dplyr::summarize(n_contrasts  = dplyr::n(),
                   n_card       = dplyr::n_distinct(gfs_card, na.rm = TRUE),
                   .groups      = "drop") %>%
  dplyr::filter(n_contrasts >= 2L, n_card > 1L)

if (nrow(drift_check) > 0) {
  bad <- paste(head(drift_check$pathway_id, 5), collapse = ", ")
  stop("[adapter] genes_full_set cardinality DRIFTS across contrasts for ", nrow(drift_check), " pathway_id(s): ", bad)
}
message("[adapter] Contrast-invariance OK: ", length(unique(u$pathway_id)), " unique pathways.")

# ===========================================================================
# 4. Per-type membership assertion
# ===========================================================================
.parse_genes <- function(x) {
  if (is.na(x) || x == "") character(0) else strsplit(x, "/", fixed = TRUE)[[1]]
}

subset_rows <- u[u$entity_type %in% c("Pathway", "GATOM"), ]
if (nrow(subset_rows) > 0) {
  subset_violations <- vapply(seq_len(nrow(subset_rows)), function(i) {
    ce  <- .parse_genes(subset_rows$core_enrichment[i])
    gfs <- .parse_genes(subset_rows$genes_full_set[i])
    length(ce) > 0 && !all(ce %in% gfs)
  }, logical(1))
  if (any(subset_violations)) {
    n_bad <- sum(subset_violations)
    ex <- subset_rows$pathway_id[which(subset_violations)[1]]
    stop("[adapter] ", n_bad, " Pathway/GATOM row(s) violate core_enrichment ⊆ genes_full_set (ex: ", ex, ")")
  }
  message("[adapter] Pathway/GATOM membership OK: core_enrichment ⊆ genes_full_set for all rows.")
}

nonpw_rows <- u[u$entity_type %in% c("TF", "PROGENy"), ]
if (nrow(nonpw_rows) > 0) {
  nonpw_violations <- vapply(seq_len(nrow(nonpw_rows)), function(i) {
    ce  <- sort(.parse_genes(nonpw_rows$core_enrichment[i]))
    gfs <- sort(.parse_genes(nonpw_rows$genes_full_set[i]))
    !identical(ce, gfs)
  }, logical(1))
  if (any(nonpw_violations)) {
    n_bad <- sum(nonpw_violations)
    ex <- nonpw_rows$pathway_id[which(nonpw_violations)[1]]
    stop("[adapter] ", n_bad, " TF/PROGENy row(s) violate genes_full_set == core_enrichment (ex: ", ex, ")")
  }
  message("[adapter] TF/PROGENy membership OK: genes_full_set == core_enrichment for all rows.")
}

# Write repaired/reconstructed CSV
readr::write_csv(round_numeric_cols(u %>% dplyr::select(-gfs_card)), UNIFIED)
message("[adapter] master_unified.csv written (conformant).")

# ===========================================================================
# 5. master_de_table.csv — validate and repair
# ===========================================================================
REQ_DE <- c("gene_symbol","t","logFC","adj.P.Val","contrast")
if (!file.exists(DE_MST)) {
  message("[adapter] master_de_table.csv absent — running-sum panel skipped.")
} else {
  de <- readr::read_csv(DE_MST, show_col_types = FALSE)
  # Ensure 't' is present
  if (!"t" %in% names(de)) {
    pcol <- intersect(c("P.Value","PValue","pvalue","adj.P.Val","padj"), names(de))
    if (length(pcol) && "logFC" %in% names(de)) {
      de$t <- sign(de$logFC) * qnorm(1 - pmin(pmax(de[[pcol[1]]], 1e-300), 1) / 2)
      message("[adapter] master_de_table.csv had no 't' — derived t-statistic proxy.")
    } else {
      stop("[adapter] master_de_table.csv has neither 't' nor logFC+p-value.")
    }
  }
  # Ensure adj.P.Val is present
  if (!"adj.P.Val" %in% names(de) && "padj" %in% names(de)) {
    de$adj.P.Val <- de$padj
  }
  miss <- setdiff(REQ_DE, names(de))
  if (length(miss)) stop("[adapter] master_de_table.csv missing columns: ", paste(miss, collapse = ", "))
  if (any(is.na(de$gene_symbol))) de <- de[!is.na(de$gene_symbol), ]
  readr::write_csv(round_numeric_cols(de), DE_MST)
  message("[adapter] master_de_table.csv conformant.")
}

# ===========================================================================
# 6. Emit explorer_manifest.json
# ===========================================================================
atlas_sha256 <- tryCatch({
  digest::digest(paste(atlas_universe, collapse = "\n"), algo = "sha256", serialize = FALSE)
}, error = function(e) {
  tryCatch(tools::md5sum(ATLAS_FP), error = function(e2) "unavailable")
})

git_commit <- tryCatch(
  trimws(system("git rev-parse --short HEAD 2>/dev/null", intern = TRUE, ignore.stderr = TRUE)),
  error = function(e) "unknown"
)
if (length(git_commit) == 0) git_commit <- "unknown"

manifest <- list(
  produced_by          = "consolidate_explorer_bundle.R",
  produced_at          = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  git_commit           = git_commit,
  atlas_file           = ATLAS_FP,
  atlas_sha256         = atlas_sha256,
  atlas_size           = atlas_n,
  namespace            = "mouse_symbol_GRCm39",
  genes_full_set_basis = "atlas_intersected",
  n_rows               = nrow(u),
  contrasts            = sort(unique(u$contrast)),
  entity_types         = sort(unique(u$entity_type))
)

if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(manifest, MANIFEST, pretty = TRUE, auto_unbox = TRUE)
} else {
  to_json_str <- function(x) {
    if (is.null(x))       return("null")
    if (is.logical(x))    return(tolower(as.character(x)))
    if (is.numeric(x))    return(as.character(x))
    if (is.character(x) && length(x) == 1) return(sprintf('"%s"', gsub('"', '\\"', x)))
    if (is.character(x))  return(sprintf('[%s]', paste(sprintf('"%s"', gsub('"', '\\"', x)), collapse = ",")))
    return(sprintf('"%s"', as.character(x)[1]))
  }
  lines <- c("{")
  nms <- names(manifest)
  for (i in seq_along(nms)) {
    sep <- if (i < length(nms)) "," else ""
    lines <- c(lines, sprintf('  "%s": %s%s', nms[i], to_json_str(manifest[[nms[i]]]), sep))
  }
  lines <- c(lines, "}")
  writeLines(lines, MANIFEST)
}
message("[adapter] explorer_manifest.json written.")

# ===========================================================================
# 7. Document the bundle in interactive/README.md
# ===========================================================================
INT_DIR <- Sys.getenv("EXPLORER_INTERACTIVE_DIR", unset = "")
if (!nzchar(INT_DIR)) INT_DIR <- DIR_INTERACTIVE
INT_DIR <- sub("/+$", "", INT_DIR)
dir.create(INT_DIR, recursive = TRUE, showWarnings = FALSE)
readme_path <- file.path(INT_DIR, "README.md")

caption_block <- paste0(
"## master_unified.csv (explorer bundle view)\n",
"**The consolidated pathway-explorer input accumulator (entity_type-tagged Pathway/TF/PROGENy/GATOM rows, lowercase nes, genes_full_set = membership ∩ atlas, contrast-invariant): the cross-stratum enrichment/activity universe the dashboards lay out.**\n",
"| | |\n",
"|---|---|\n",
"| Script   | `02_analysis/scripts/pathway_explorer_adapter/consolidate_explorer_bundle.R` |\n",
"| Function | (validation + manifest) |\n",
"| Config   | `paths.master = 03_results/master/` |\n",
"| Input    | `03_results/master/master_unified.csv; 03_results/master/atlas_gene_universe.txt` |\n\n",
"## master_de_table.csv (explorer schema view)\n",
"**The explorer-schema DE view (gene_symbol, t, logFC, adj.P.Val, contrast): the t-ranked per-contrast gene list that drives the dashboards' running-sum panel.**\n",
"| | |\n",
"|---|---|\n",
"| Script   | `02_analysis/scripts/pathway_explorer_adapter/consolidate_explorer_bundle.R` |\n",
"| Function | (validation) |\n",
"| Config   | `paths.master = 03_results/master/` |\n",
"| Input    | `03_results/master/master_de_table.csv` |\n"
)

if (!file.exists(readme_path)) {
  writeLines("# 03_results/interactive — Pathway Explorer Dashboards\n", readme_path)
}
existing <- readLines(readme_path)
if (!any(grepl("## master_unified.csv \\(explorer bundle view\\)", existing))) {
  cat(caption_block, file = readme_path, append = TRUE)
  message("[adapter] README.md bundle captions written to ", readme_path)
}

message("[adapter] consolidate_explorer_bundle.R complete — all validations passed.")
