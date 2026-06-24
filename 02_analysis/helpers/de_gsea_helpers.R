## de_gsea_helpers.R — STING-cGAS-GSE329522 DE/GSEA COMPUTE toolbox (NO plotting)
## ===========================================================================
## The keystone compute helper for the standard bulk-RNA-seq sweep. Every new
## arm (GSEA / CoReSh / PROGENy / GATOM) sources this AFTER the project config:
##
##   source("02_analysis/config/config.R")        # PROJECT_ROOT, YAML_CONFIG, stage_dir, DIR_*, etc.
##   source("02_analysis/helpers/de_gsea_helpers.R")
##
## Ported from 14839-DM-cGAS/02_analysis/helpers/de_gsea_helpers.R and adapted to
## STING's config API + existing compute idiom (03_decoupler_tf.R: load_or_compute +
## schema-validated master table + tidy CSVs; rank metric = the limma-trend t-statistic).
##
## SCOPE — COMPUTE ONLY. No ggplot/ggsave, no DE re-fitting, no TMM/voom. The DE is
## already done (02_de_limma_trend.R → 03_results/objects/02_de_results.rds); we only
## consume its topTables and rank on the existing `t`. Captions/figures live in the
## figure-style viz shim (02_analysis/helpers/figure-style/figure_helpers.R), NOT here.
##
## CONFIG KEYS READ (all via the project config.R accessors / YAML_CONFIG — never hardcoded):
##   thresholds.gsea_min_size  -> GSEA_MIN_SIZE   (15)        min gene-set size
##   thresholds.gsea_max_size  -> GSEA_MAX_SIZE   (500)       max gene-set size
##   thresholds.gsea_seed      -> GSEA_SEED       (123)       fgsea RNG seed (reproducibility)
##   thresholds.gsea_nperm     -> GSEA_NPERM      (100000)    permutations (fgsea simple path only)
##   thresholds.gsea_fdr       -> GSEA_FDR_CUTOFF (0.05)      display FDR (callers, not run)
##   RANK_METRIC               -> "t"             (config.R constant; NEVER logFC)
##   design.contrasts          -> YAML_CONFIG$design$contrasts (per-contrast `name`)
##   databases.msigdb          -> 8 collections (category + subcategory + name)
##   databases.custom          -> custom DBs (name + path RDS, already mouse symbols)
##   schemas.master_gsea_table -> the master GSEA schema (required_columns)
##   paths.objects / paths.master / project.species  (via DIR_OBJECTS / DIR_MASTER / SPECIES)
##
## COMPOSES WITH config.R — this file requires the project config to have been sourced
## already. It uses config.R's accessors/constants (PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS,
## DIR_MASTER, SPECIES, GSEA_MIN_SIZE/MAX_SIZE/SEED/NPERM, RANK_METRIC, the `%||%` op).
## It deliberately re-defines `load_or_compute` with a path-keyed signature (see note at
## that function) so the new sweep arms get an RDS cache keyed on an explicit cache_path,
## matching the 14839 idiom and the per-arm checkpoint pattern; the config.R variant
## (basename-under-DIR_OBJECTS) keeps working for the legacy 00–03 scripts that call it
## before this file is sourced.
##
## LAZY-DEPS DESIGN — heavy enrichment packages (fgsea, msigdbr, decoupleR) are NOT
## attached at source() time. Every function that needs one calls requireNamespace()
## and stops with an actionable install hint only when actually invoked. This lets the
## file source() cleanly in any container (e.g. a doc/lint pass) without those installed.
## Light tidyverse deps (dplyr/tibble/readr) ARE required up front — they are part of
## the project's base stack (config.R::load_packages) and every compute script needs them.
## ===========================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(readr)
})

## Guard: this helper must be sourced AFTER the project config.R (which sets PROJECT_ROOT,
## YAML_CONFIG, DIR_OBJECTS/DIR_MASTER, SPECIES, the GSEA_* constants and `%||%`).
if (!exists("YAML_CONFIG"))
  stop("de_gsea_helpers.R: source 02_analysis/config/config.R first (need YAML_CONFIG, DIR_OBJECTS, ...).")

## Null/empty-coalesce — config.R already defines `%||%`, but re-define defensively
## (idempotent) so this file is robust if sourced standalone after a partial config.
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
}

## ---------------------------------------------------------------------------
## (0) GSEA parameter accessors — read straight from config.R constants when present,
##     else fall back to YAML_CONFIG$thresholds with the documented defaults. Centralised
##     so callers never re-read thresholds (AGENTS.md rule 1: config, not hardcoding).
## ---------------------------------------------------------------------------
.gsea_min_size <- function() if (exists("GSEA_MIN_SIZE")) GSEA_MIN_SIZE else (YAML_CONFIG$thresholds$gsea_min_size %||% 15L)
.gsea_max_size <- function() if (exists("GSEA_MAX_SIZE")) GSEA_MAX_SIZE else (YAML_CONFIG$thresholds$gsea_max_size %||% 500L)
.gsea_seed     <- function() if (exists("GSEA_SEED"))     GSEA_SEED     else (YAML_CONFIG$thresholds$gsea_seed     %||% 123L)
.gsea_nperm    <- function() if (exists("GSEA_NPERM"))    GSEA_NPERM    else (YAML_CONFIG$thresholds$gsea_nperm    %||% 100000L)
.rank_metric   <- function() if (exists("RANK_METRIC"))   RANK_METRIC   else (YAML_CONFIG$gsea$rank_metric        %||% "t")

## ---------------------------------------------------------------------------
## (1) load_or_compute(cache_path, compute_fn, force = FALSE) — RDS cache keyed on PATH.
##     Matches STING's existing load_or_compute SEMANTICS (cache-or-compute-and-save into
##     03_results/objects/) but keyed on an explicit cache_path rather than a bare
##     basename, mirroring the 14839 per-arm checkpoint idiom. A bare filename (no dir)
##     is resolved under DIR_OBJECTS so callers can pass either "gsea_msigdb_WT_heat.rds"
##     or an absolute/relative path. Idempotent + provenance-friendly: a cache hit logs and
##     returns; a miss computes, saves, returns. Re-defines the config.R name on purpose
##     (see header) — the path-keyed form is what the new sweep arms use.
## ---------------------------------------------------------------------------
load_or_compute <- function(cache_path, compute_fn, force = FALSE) {
  ## bare filename -> place under the project objects dir (03_results/objects/)
  if (identical(basename(cache_path), cache_path)) {
    odir <- if (exists("DIR_OBJECTS")) DIR_OBJECTS
            else file.path(PROJECT_ROOT, YAML_CONFIG$paths$objects %||% "03_results/objects")
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)
    cache_path <- file.path(odir, cache_path)
  } else {
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
  }
  if (file.exists(cache_path) && !force) {
    message("[load_or_compute] cache hit: ", cache_path)
    return(readRDS(cache_path))
  }
  message("[load_or_compute] computing: ", cache_path)
  obj <- compute_fn()
  saveRDS(obj, cache_path)
  obj
}

## ---------------------------------------------------------------------------
## (2) load_de_results(path = ...) — load the limma-trend topTables (the GSEA hub).
##     Default path = 03_results/objects/02_de_results.rds (produced by 02_de_limma_trend.R;
##     NEVER re-fit DE here). Returns the per-contrast NAMED LIST of topTables.
##
##     ASSUMED SHAPE of 02_de_results.rds (verified against 02_de_limma_trend.R lines
##     118-127 — owner should re-confirm in-container):
##       * a named list, names == the 7 YAML contrast names
##         (WT_heat, KO_heat, Interaction, Geno_at_39, Geno_at_37, Temp_main, Geno_main);
##       * each element is a data.frame (topTable) with gene-SYMBOL ROWNAMES and columns:
##           gene_symbol (== rownames), ensembl, logFC, AveExpr, t, P.Value, adj.P.Val, B, contrast.
##     The rank metric `t`, plus logFC / P.Value / AveExpr, are all present per element.
##     We assert the Symbol-rownames contract (the load-bearing invariant every consumer
##     relies on) and the presence of the rank metric column.
## ---------------------------------------------------------------------------
load_de_results <- function(path = NULL) {
  if (is.null(path)) {
    odir <- if (exists("DIR_OBJECTS")) DIR_OBJECTS
            else file.path(PROJECT_ROOT, YAML_CONFIG$paths$objects %||% "03_results/objects")
    path <- file.path(odir, "02_de_results.rds")
  }
  if (!file.exists(path)) stop("load_de_results: missing DE checkpoint: ", path,
                               " — run 02_de_limma_trend.R first (do NOT re-fit DE here).")
  de <- readRDS(path)
  if (!is.list(de) || is.null(names(de)) || any(names(de) == ""))
    stop("load_de_results: expected a NAMED list of per-contrast topTables; got ", class(de)[1])

  rm <- .rank_metric()
  for (cn in names(de)) {
    tt <- de[[cn]]
    if (!is.data.frame(tt)) stop("load_de_results: contrast '", cn, "' is not a data.frame.")
    ## Symbol-rownames contract: if rownames are missing but a gene_symbol column exists,
    ## restore them (02_de_limma_trend.R keeps both, but be defensive for re-keyed objects).
    if (is.null(rownames(tt)) && "gene_symbol" %in% colnames(tt)) {
      rownames(tt) <- make.unique(as.character(tt$gene_symbol)); de[[cn]] <- tt
    }
    rn <- rownames(tt)
    if (mean(grepl("^ENSMUSG", rn)) > 0.5)
      stop("load_de_results: contrast '", cn, "' has Ensembl-like rownames (", rn[1],
           ") — enrichment needs gene Symbols. Re-key the topTable to Symbol rownames.")
    if (!rm %in% colnames(tt))
      stop("load_de_results: contrast '", cn, "' lacks the rank metric column '", rm, "'.")
  }
  de
}

## ---------------------------------------------------------------------------
## (3) build_ranked_vector(topTable, metric = "t") — the GSEA ranking input.
##     Returns a named numeric vector (gene SYMBOL -> metric), de-duplicated and sorted
##     DECREASING. Names come from rownames (the Symbol contract); ties on duplicate symbols
##     keep the larger-|metric| representative (the more extreme rank, never silently
##     averaged). NA metrics, empty/NA names dropped. Infinities dropped (fgsea rejects them).
##     Fallback ladder (if `metric` absent): t -> stat -> sign(logFC)*-log10(P.Value).
## ---------------------------------------------------------------------------
build_ranked_vector <- function(topTable, metric = NULL) {
  metric <- metric %||% .rank_metric()
  cols <- colnames(topTable)
  if (metric %in% cols) {
    v <- topTable[[metric]]
  } else if ("t" %in% cols) {
    v <- topTable[["t"]]
  } else if ("stat" %in% cols) {
    v <- topTable[["stat"]]
  } else if (all(c("logFC", "P.Value") %in% cols)) {
    v <- sign(topTable[["logFC"]]) * -log10(pmax(topTable[["P.Value"]], 1e-300))
  } else {
    stop("build_ranked_vector: no usable rank column (", metric, "/t/stat/logFC+P.Value) in topTable.")
  }
  nm <- rownames(topTable)
  if (is.null(nm) && "gene_symbol" %in% cols) nm <- as.character(topTable[["gene_symbol"]])
  names(v) <- nm
  keep <- !is.na(v) & is.finite(v) & !is.na(names(v)) & names(v) != ""
  v <- v[keep]
  ## de-duplicate symbols: keep the most extreme (largest |metric|) per symbol
  if (anyDuplicated(names(v))) {
    ord <- order(abs(v), decreasing = TRUE)
    v   <- v[ord]
    v   <- v[!duplicated(names(v))]
  }
  sort(v, decreasing = TRUE)
}

## ---------------------------------------------------------------------------
## (4) gene_universe(de_results) — the background symbol universe for GSEA/CoReSh.
##     Union of every contrast's topTable symbols (rownames), de-duplicated. This is the
##     Symbol-deduped DE namespace (the modeled gene list, ~19,679 symbols) used as the
##     contrast-invariant background. Accepts either the named list from load_de_results()
##     or a single topTable data.frame.
## ---------------------------------------------------------------------------
gene_universe <- function(de_results) {
  syms_of <- function(tt) {
    rn <- rownames(tt)
    if (is.null(rn) && "gene_symbol" %in% colnames(tt)) rn <- as.character(tt$gene_symbol)
    rn
  }
  syms <- if (is.data.frame(de_results)) syms_of(de_results)
          else unlist(lapply(de_results, syms_of), use.names = FALSE)
  syms <- syms[!is.na(syms) & syms != ""]
  sort(unique(syms))
}

## ---------------------------------------------------------------------------
## (5) load_msigdb_collection(category, subcategory, species = "Mus musculus")
##     -> a NAMED LIST of gene-symbol sets (gs_name -> character vector), the fgsea
##     `pathways` format. Used for the 8 config MSigDB collections (databases.msigdb).
##
##     msigdbr API note (handled here): v8/v26 renamed the args category->collection and
##     subcategory->subcollection and rejects the legacy db_species="MM"; we detect the
##     installed API via formals() and pass species=<common name> ("Mus musculus") only.
##     v26 also renamed CP:KEGG -> CP:KEGG_LEGACY; we retry with that suffix on an
##     "Unknown subcollection" error. Returns list() (empty) if the collection yields no sets.
## ---------------------------------------------------------------------------
load_msigdb_collection <- function(category, subcategory = "", species = NULL) {
  if (!requireNamespace("msigdbr", quietly = TRUE))
    stop("load_msigdb_collection: package 'msigdbr' is required. ",
         "Install with BiocManager::install('msigdbr').")
  species <- species %||% (if (exists("SPECIES")) SPECIES else (YAML_CONFIG$project$species %||% "Mus musculus"))
  sub_ <- subcategory %||% ""
  has_new <- "collection" %in% names(formals(msigdbr::msigdbr))   # v8+/v26 vs legacy v7.5

  fetch <- function(sub_use) {
    if (has_new) {
      msigdbr::msigdbr(species = species, collection = category,
                       subcollection = if (nzchar(sub_use)) sub_use else NULL)
    } else {
      msigdbr::msigdbr(species = species, category = category,
                       subcategory = if (nzchar(sub_use)) sub_use else NULL)
    }
  }
  df <- tryCatch(fetch(sub_), error = function(e) {
    ## v26 CP:KEGG -> CP:KEGG_LEGACY rename retry
    if (nzchar(sub_) && grepl("[Uu]nknown subcollection", conditionMessage(e))) {
      alt <- sub("CP:KEGG$", "CP:KEGG_LEGACY", sub_)
      if (alt != sub_) {
        message("  [msigdb ", category, "/", sub_, "] retrying subcollection='", alt, "' (msigdbr v26 rename)")
        return(tryCatch(fetch(alt), error = function(e2) {
          warning("load_msigdb_collection: ", category, "/", sub_, " failed: ", conditionMessage(e2)); NULL
        }))
      }
    }
    warning("load_msigdb_collection: ", category, "/", sub_, " failed: ", conditionMessage(e)); NULL
  })
  if (is.null(df) || nrow(df) == 0) return(list())
  ## msigdbr column for the symbol differs slightly across versions: prefer gene_symbol.
  sym_col <- if ("gene_symbol" %in% colnames(df)) "gene_symbol" else "human_gene_symbol"
  split(as.character(df[[sym_col]]), as.character(df$gs_name))
}

## ---------------------------------------------------------------------------
## (6) load_custom_geneset(path) — read a custom-DB .rds into a NAMED LIST of mouse-symbol
##     sets (the fgsea `pathways` format). The config databases.custom RDS come in TWO shapes;
##     both are normalised here:
##       (a) toolkit T2G/T2N list  -> list(T2G = data.frame(gs_name, gene_symbol), T2N, ...)
##           (e.g. transportdb_genesets.rds, mito_*.rds) -> split T2G$gene_symbol by gs_name.
##       (b) named list of symbol vectors -> list(<SetName> = c("Gene1","Gene2",...), ...)
##           (e.g. lombardi2022_hif_consensus_mouse.rds = list(Lombardi2022_HIF = <genes>),
##            with a `provenance` attribute) -> returned as-is (coerced to character sets).
##     Already mouse symbols per config — NO species conversion here.
## ---------------------------------------------------------------------------
load_custom_geneset <- function(path) {
  if (!file.exists(path)) stop("load_custom_geneset: file not found: ", path)
  obj <- readRDS(path)
  ## (a) toolkit T2G/T2N shape
  if (is.list(obj) && all(c("T2G") %in% names(obj)) &&
      is.data.frame(obj$T2G) && all(c("gs_name", "gene_symbol") %in% colnames(obj$T2G))) {
    return(split(as.character(obj$T2G$gene_symbol), as.character(obj$T2G$gs_name)))
  }
  ## (b) named list of gene-symbol vectors. lapply() keeps names but drops list-level
  ##     attrs (e.g. `provenance`), giving a clean fgsea pathways list.
  if (is.list(obj) && !is.null(names(obj)) && all(nzchar(names(obj)))) {
    return(lapply(obj, function(g) unique(as.character(g))))
  }
  stop("load_custom_geneset: unrecognised RDS shape at ", path,
       " (expected toolkit {T2G,T2N} list or a named list of gene-symbol vectors).")
}

## ---------------------------------------------------------------------------
## (7a) direction_from_nes(nes) — small sign helper: NES > 0 -> "Up", else "Down".
##      NA NES -> NA_character_. The single source of truth for the master `direction`
##      column across every enrichment arm (keeps the schema value consistent).
## ---------------------------------------------------------------------------
direction_from_nes <- function(nes) {
  ifelse(is.na(nes), NA_character_, ifelse(nes > 0, "Up", "Down"))
}

## ---------------------------------------------------------------------------
## (7b) run_fgsea(ranked, gene_sets, database, contrast, ...) — fgsea wrapper returning a
##      TIDY data.frame in the master_gsea_table schema:
##        pathway_id, pathway_name, database, nes, pvalue, padj, set_size,
##        core_enrichment, contrast, direction
##      (lowercase `nes`, per the master/explorer contract). `database` and `contrast` are
##      filled from the args; `direction` from sign(nes) via direction_from_nes().
##
##      Engine: fgsea::fgseaMultilevel (the modern, accurate-tail default) when available,
##      else fgsea::fgsea(nperm=) (older API). Size bounds from config (gsea_min_size /
##      gsea_max_size); seeded from gsea_seed for reproducibility. `eps` exposes the
##      multilevel p-value floor (0 = estimate to machine precision; the run-time default).
##      `core_enrichment` = the leadingEdge genes, slash-joined (matches the toolkit /
##      master convention so downstream viz reads one format).
##
##      pathway_name defaults to a cleaned pathway_id (prefix stripped, "_"->" "); callers
##      with a T2N description table can override post-hoc. Returns a 0-row schema-shaped
##      data.frame if fgsea yields nothing (so binds/append stay well-typed).
## ---------------------------------------------------------------------------
run_fgsea <- function(ranked, gene_sets, database, contrast,
                      minSize = NULL, maxSize = NULL,
                      eps = 0, nperm = NULL, seed = NULL,
                      pathway_names = NULL) {
  if (!requireNamespace("fgsea", quietly = TRUE))
    stop("run_fgsea: package 'fgsea' is required. Install with BiocManager::install('fgsea').")

  minSize <- minSize %||% .gsea_min_size()
  maxSize <- maxSize %||% .gsea_max_size()
  nperm   <- nperm   %||% .gsea_nperm()
  seed    <- seed    %||% .gsea_seed()

  ## drop empty gene sets up front (fgsea errors on zero-length pathways)
  gene_sets <- gene_sets[vapply(gene_sets, function(g) length(g) > 0L, logical(1))]
  if (length(gene_sets) == 0L || length(ranked) == 0L) return(.empty_gsea_df())

  set.seed(seed)
  res <- tryCatch({
    if ("fgseaMultilevel" %in% getNamespaceExports("fgsea")) {
      fgsea::fgseaMultilevel(pathways = gene_sets, stats = ranked,
                             minSize = minSize, maxSize = maxSize, eps = eps)
    } else {
      ## older fgsea API: permutation-based call (nperm) — eps not supported there
      fgsea::fgsea(pathways = gene_sets, stats = ranked,
                   minSize = minSize, maxSize = maxSize, nperm = nperm)
    }
  }, error = function(e) {
    warning("run_fgsea: fgsea failed for database='", database, "', contrast='", contrast,
            "': ", conditionMessage(e)); NULL
  })
  if (is.null(res) || nrow(res) == 0L) return(.empty_gsea_df())

  res <- as.data.frame(res, stringsAsFactors = FALSE)
  ## leadingEdge is a list-column -> slash-joined string (the master core_enrichment format)
  core <- vapply(res$leadingEdge, function(le) paste(le, collapse = "/"), character(1))

  pid <- as.character(res$pathway)
  pname <- if (!is.null(pathway_names) && all(pid %in% names(pathway_names))) {
    unname(pathway_names[pid])
  } else {
    .clean_pathway_name(pid)
  }

  out <- data.frame(
    pathway_id      = pid,
    pathway_name    = pname,
    database        = database,
    nes             = as.numeric(res$NES),
    pvalue          = as.numeric(res$pval),
    padj            = as.numeric(res$padj),
    set_size        = as.integer(res$size),
    core_enrichment = core,
    contrast        = contrast,
    direction       = direction_from_nes(as.numeric(res$NES)),
    stringsAsFactors = FALSE
  )
  ## stable order: most significant first (ties by |NES|)
  out[order(out$padj, -abs(out$nes)), , drop = FALSE]
}

## Empty schema-shaped GSEA data.frame (well-typed; used on no-result so binds stay clean).
.empty_gsea_df <- function() {
  data.frame(
    pathway_id = character(), pathway_name = character(), database = character(),
    nes = numeric(), pvalue = numeric(), padj = numeric(), set_size = integer(),
    core_enrichment = character(), contrast = character(), direction = character(),
    stringsAsFactors = FALSE
  )
}

## Basic pathway-name cleaner (prefix strip + underscores->spaces + title case).
## Internal — viz owns pretty formatting; this just gives a readable default `pathway_name`.
.clean_pathway_name <- function(ids) {
  x <- as.character(ids)
  prefixes <- c("HALLMARK_", "KEGG_", "KEGG_MEDICUS_", "REACTOME_", "WP_",
                "GOBP_", "GOCC_", "GOMF_", "GO_", "MITOPATHWAYS_", "MITOXPLORER_",
                "MITOCARTA_", "TRANSPORTDB_", "GTRD_", "TFT_")
  for (p in prefixes) x <- sub(paste0("^", p), "", x)
  x <- gsub("_", " ", x)
  vapply(x, function(s) {
    w <- strsplit(tolower(s), " ", fixed = TRUE)[[1]]
    w <- ifelse(nzchar(w), paste0(toupper(substr(w, 1, 1)), substr(w, 2, nchar(w))), w)
    paste(w, collapse = " ")
  }, character(1), USE.NAMES = FALSE)
}

## ---------------------------------------------------------------------------
## (8a) round_numeric_cols(df, sig = 9L) — byte-stability helper for master CSV writes.
##      signif()-rounds every is.double() column to `sig` significant figures so that
##      re-running fgsea/decoupleR (which can perturb doubles at the ULP level) produces NO
##      meaningless diff in the master CSVs. 9 sig figs preserves all biological content
##      (p/NES/logFC/t carry <=4 meaningful digits) while discarding IEEE-754 noise, and is
##      idempotent. Integer-valued doubles (set_size) are <=9 digits so signif() keeps them
##      exact. Character/logical cols (core_enrichment, pathway names, ...) are NEVER touched.
## ---------------------------------------------------------------------------
round_numeric_cols <- function(df, sig = 9L) {
  for (nm in names(df)) {
    col <- df[[nm]]
    if (is.double(col)) df[[nm]] <- signif(col, sig)
  }
  df
}

## ---------------------------------------------------------------------------
## (8b) append_master_table(df, master_path, key_col = "database") — the canonical master
##      writer for compute scripts. IDEMPOTENT append keyed on key_col (default "database"):
##        1. validate df against the master_gsea_table schema (required_columns);
##        2. enforce lowercase `nes` (rename NES->nes if needed) — the explorer contract;
##        3. read any existing CSV, DROP every row whose key_col value is in df[[key_col]],
##           bind_rows the new df, round to 9 sig figs (byte-stable), write.
##      Re-running an arm therefore replaces exactly its own rows and produces no byte-diff.
##      master_path may be a bare filename (resolved under 03_results/master/) or a full path.
##      The schema is read from YAML_CONFIG$schemas$master_gsea_table (config, not hardcoded).
## ---------------------------------------------------------------------------
append_master_table <- function(df, master_path, key_col = "database") {
  ## enforce lowercase nes BEFORE schema validation (explorer/master contract)
  if ("NES" %in% colnames(df) && !"nes" %in% colnames(df))
    df <- dplyr::rename(df, nes = NES)

  req <- YAML_CONFIG$schemas$master_gsea_table$required_columns
  if (is.null(req))
    stop("append_master_table: schemas.master_gsea_table.required_columns missing from config.")
  missing <- setdiff(req, colnames(df))
  if (length(missing))
    stop("append_master_table: df missing required schema cols: ", paste(missing, collapse = ", "))
  if (!key_col %in% colnames(df))
    stop("append_master_table: key_col '", key_col, "' not present in df.")

  ## resolve master path (bare filename -> 03_results/master/)
  if (identical(basename(master_path), master_path)) {
    mdir <- if (exists("DIR_MASTER")) DIR_MASTER
            else file.path(PROJECT_ROOT, YAML_CONFIG$paths$master %||% "03_results/master")
    dir.create(mdir, recursive = TRUE, showWarnings = FALSE)
    master_path <- file.path(mdir, master_path)
  } else {
    dir.create(dirname(master_path), recursive = TRUE, showWarnings = FALSE)
  }

  new_keys <- unique(df[[key_col]])
  if (file.exists(master_path)) {
    old <- readr::read_csv(master_path, show_col_types = FALSE, progress = FALSE)
    if (key_col %in% colnames(old)) old <- old[!(old[[key_col]] %in% new_keys), , drop = FALSE]
    combined <- dplyr::bind_rows(old, df)             # bind_rows unions cols, NA-filling
  } else {
    combined <- df
  }
  readr::write_csv(round_numeric_cols(combined), master_path)
  invisible(master_path)
}

## ---------------------------------------------------------------------------
## End of de_gsea_helpers.R — STING-cGAS-GSE329522 COMPUTE toolbox.
## Provides: load_or_compute, load_de_results, build_ranked_vector, gene_universe,
##   load_msigdb_collection, load_custom_geneset, run_fgsea, direction_from_nes,
##   round_numeric_cols, append_master_table  (+ internal .gsea_*/.clean_pathway_name/.empty_gsea_df).
## ---------------------------------------------------------------------------
message("[de_gsea_helpers] loaded (STING compute toolbox; lazy deps: fgsea/msigdbr/decoupleR).")
