#!/usr/bin/env Rscript
# emit_explorer_universe.R — emit the set-level candidate-pathway universe for pathway-explorer
# Adapt to STING-cGAS-GSE329522 RDS checkpoints and folders.
# Run from project root.

source("02_analysis/config/config.R")           # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS, DIR_MASTER, SPECIES, %||%
if (file.exists("02_analysis/helpers/de_gsea_helpers.R")) {
  source("02_analysis/helpers/de_gsea_helpers.R")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(purrr)
})
options(stringsAsFactors = FALSE)

obj_dir <- DIR_OBJECTS
mst_dir <- DIR_MASTER
dir.create(mst_dir, recursive = TRUE, showWarnings = FALSE)

MINSZ   <- as.integer(YAML_CONFIG$thresholds$gsea_min_size %||% 15L)
MAXSZ   <- as.integer(YAML_CONFIG$thresholds$gsea_max_size %||% 500L)
SPECIES <- YAML_CONFIG$project$species %||% "Mus musculus"
OUT_FP  <- file.path(mst_dir, "explorer_universe.csv")

# ---------------------------------------------------------------------------
# Load modeled gene universe (atlas_gene_universe equivalent)
# ---------------------------------------------------------------------------
atlas_fp <- file.path(obj_dir, "gene_universe.txt")
if (!file.exists(atlas_fp)) {
  stop("gene_universe.txt absent — run 11_emit_universe.R first to produce the modeled gene list.")
}
atlas_universe <- readLines(atlas_fp)
atlas_universe <- unique(atlas_universe[!is.na(atlas_universe) & nzchar(atlas_universe)])
message(sprintf("[universe] Loaded %d symbols from gene_universe.txt", length(atlas_universe)))

# Helper builder: size-filter on original size, intersect with atlas, format output
.build_rows <- function(members_by_set, entity_type, size_filter = TRUE) {
  ids <- names(members_by_set)
  out <- purrr::map_dfr(ids, function(id) {
    mem  <- unique(as.character(members_by_set[[id]]))
    mem  <- mem[!is.na(mem) & nzchar(mem)]
    if (size_filter && (length(mem) < MINSZ || length(mem) > MAXSZ)) return(NULL)
    keep <- intersect(mem, atlas_universe)
    if (length(keep) == 0L) return(NULL)
    tibble::tibble(pathway_id = id, entity_type = entity_type,
                   genes_full_set = paste(keep, collapse = "/"))
  })
  out
}

# ===========================================================================
# (1) Pathway — MSigDB (8 collections) + custom (4) + CoReSh GMT
# ===========================================================================
suppressPackageStartupMessages({ library(msigdbr) })

.msigdbr_sets <- function(collection, subcollection) {
  call_msig <- function(coll, sub) {
    params      <- names(formals(msigdbr::msigdbr))
    use_new_api <- "collection" %in% params
    if (use_new_api) {
      if (nzchar(sub)) msigdbr(species = SPECIES, collection = coll, subcollection = sub)
      else             msigdbr(species = SPECIES, collection = coll)
    } else {
      if (nzchar(sub)) msigdbr(species = SPECIES, category = coll, subcategory = sub)
      else             msigdbr(species = SPECIES, category = coll)
    }
  }
  df <- tryCatch(call_msig(collection, subcollection), error = function(e) {
    # CP:KEGG -> CP:KEGG_LEGACY fallback for modern msigdbr versions
    if (nzchar(subcollection) && grepl("Unknown|unknown", conditionMessage(e))) {
      alt <- sub("CP:KEGG$", "CP:KEGG_LEGACY", subcollection)
      if (alt != subcollection) {
        message("  [msigdbr] retry with subcollection='", alt, "' (v26 rename)")
        return(call_msig(collection, alt))
      }
    }
    stop(e)
  })
  if (is.null(df) || nrow(df) == 0) return(list())
  split(as.character(df$gene_symbol), as.character(df$gs_name))
}

pathway_rows <- list()
for (db in YAML_CONFIG$databases$msigdb) {
  sub_ <- db$subcategory %||% ""
  sets <- .msigdbr_sets(db$category, sub_)
  rows <- .build_rows(sets, "Pathway", size_filter = TRUE)
  message(sprintf("[universe] MSigDB %-15s: %d sets -> %d kept.",
                  db$name, length(sets), nrow(rows)))
  pathway_rows[[db$name]] <- rows
}

# Custom DB sets
custom_outfile <- c(
  MitoPathways      = "geneset_custom_MitoPathways.rds",
  MitoXplorer       = "geneset_custom_MitoXplorer.rds",
  TransportDB       = "geneset_custom_TransportDB.rds",
  Lombardi2022_HIF  = "geneset_custom_Lombardi2022_HIF.rds"
)
for (dbn in names(custom_outfile)) {
  fp <- file.path(obj_dir, custom_outfile[[dbn]])
  if (!file.exists(fp)) {
    warning("Custom collection RDS missing: ", fp)
    next
  }
  db <- readRDS(fp)
  # Check if in standard list or list(T2G) shape
  if (is.list(db) && all(c("T2G") %in% names(db)) && is.data.frame(db$T2G)) {
    sets <- split(as.character(db$T2G$gene_symbol), as.character(db$T2G$gs_name))
  } else if (is.list(db) && !is.null(names(db))) {
    sets <- lapply(db, function(g) unique(as.character(g)))
  } else {
    warning("Unexpected custom RDS format for ", dbn)
    next
  }
  rows <- .build_rows(sets, "Pathway", size_filter = TRUE)
  message(sprintf("[universe] Custom %-15s: %d sets -> %d kept.", dbn, length(sets), nrow(rows)))
  pathway_rows[[dbn]] <- rows
}

# CoReSh derived sets
coresh_fp <- file.path(obj_dir, "coresh_derived_sets.rds")
if (file.exists(coresh_fp)) {
  coresh_sets <- readRDS(coresh_fp)
  rows <- .build_rows(coresh_sets, "Pathway", size_filter = TRUE)
  message(sprintf("[universe] CoReSh %-15s: %d sets -> %d kept.", "derived", length(coresh_sets), nrow(rows)))
  pathway_rows[["CoReSh"]] <- rows
}

pathway_df <- dplyr::bind_rows(pathway_rows)

# ===========================================================================
# (2) TF — full CollecTRI network
# ===========================================================================
tf_fp <- file.path(obj_dir, "net_collectri_mouse.rds")
if (!file.exists(tf_fp)) {
  stop("net_collectri_mouse.rds absent — run TF decoupling first.")
}
collectri_net <- readRDS(tf_fp)
tf_sets <- split(as.character(collectri_net$target), as.character(collectri_net$source))
names(tf_sets) <- paste0("TF_", names(tf_sets))
tf_df <- .build_rows(tf_sets, "TF", size_filter = FALSE)
message(sprintf("[universe] TF (CollecTRI): %d regulons -> %d kept.", length(tf_sets), nrow(tf_df)))

# ===========================================================================
# (3) PROGENy — footprints
# ===========================================================================
pg_fp <- file.path(obj_dir, "net_progeny_mouse.rds")
if (!file.exists(pg_fp)) {
  stop("net_progeny_mouse.rds absent — run PROGENy first.")
}
progeny_net <- readRDS(pg_fp)
pg_sets <- split(as.character(progeny_net$target), as.character(progeny_net$source))
names(pg_sets) <- paste0("PROGENY_", names(pg_sets))
pg_df <- .build_rows(pg_sets, "PROGENy", size_filter = FALSE)
message(sprintf("[universe] PROGENy: %d footprints -> %d kept.", length(pg_sets), nrow(pg_df)))

# ===========================================================================
# (4) GATOM modules
# ===========================================================================
gatom_df <- tibble::tibble(pathway_id = character(), entity_type = character(), genes_full_set = character())
gatom_fp <- file.path(obj_dir, "10_gatom_all.rds")
if (file.exists(gatom_fp)) {
  suppressPackageStartupMessages({ library(igraph) })
  gatom_results <- readRDS(gatom_fp)
  edge_genes <- function(res) {
    if (!isTRUE(res$success) || is.null(res$module) || igraph::ecount(res$module) == 0) return(character(0))
    ed <- igraph::as_data_frame(res$module, what = "edges")
    g  <- if ("Symbol" %in% names(ed)) ed$Symbol else ed$gene
    unique(g[!is.na(g) & g != ""])
  }
  
  # Determine focal contrasts
  focal <- if (!is.null(YAML_CONFIG$design$contrasts)) {
    intersect(vapply(YAML_CONFIG$design$contrasts, function(x) x$name, ""), names(gatom_results))
  } else {
    names(gatom_results)
  }
  
  # Get networks present (e.g. kegg, combined)
  networks <- names(gatom_results[[1]]$modules)
  gatom_sets <- list()
  for (net in networks) {
    g_all <- unique(unlist(lapply(focal, function(co) {
      r <- gatom_results[[co]]$modules[[net]]
      if (is.null(r)) character(0) else edge_genes(r)
    })))
    gatom_sets[[sprintf("GATOM_%s", toupper(net))]] <- g_all
  }
  gatom_df <- .build_rows(gatom_sets, "GATOM", size_filter = FALSE)
  message(sprintf("[universe] GATOM: %d net-unions -> %d kept.", length(gatom_sets), nrow(gatom_df)))
} else {
  message("[universe] 10_gatom_all.rds absent — skipping GATOM universe.")
}

# ===========================================================================
# Assemble + write
# ===========================================================================
universe <- dplyr::bind_rows(pathway_df, tf_df, pg_df, gatom_df) %>%
  dplyr::distinct(pathway_id, entity_type, .keep_all = TRUE) %>%
  dplyr::arrange(entity_type, pathway_id) %>%
  dplyr::select(pathway_id, entity_type, genes_full_set)

stopifnot(nrow(universe) > 0)
stopifnot(!any(is.na(universe$genes_full_set) | universe$genes_full_set == ""))
stopifnot(!any(duplicated(paste(universe$pathway_id, universe$entity_type, sep = "\r"))))

readr::write_csv(universe, OUT_FP)
message(sprintf("[universe] Successfully wrote set-level explorer_universe.csv (%d rows) to %s",
                nrow(universe), OUT_FP))
by_type <- universe %>% dplyr::count(entity_type, name = "n")
for (i in seq_len(nrow(by_type))) {
  message(sprintf("  %-8s: %d rows", by_type$entity_type[i], by_type$n[i]))
}
