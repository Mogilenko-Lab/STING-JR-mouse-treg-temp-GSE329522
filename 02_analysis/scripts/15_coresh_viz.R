#!/usr/bin/env Rscript
# 15_coresh_viz.R — VIZ
## CoReSh compendium-ranking + derived-GSEA visualisation (stage 08_coresh).
## Run from project root: Rscript 02_analysis/scripts/15_coresh_viz.R
##
## VIZ-ONLY. Reads objects/coresh_ranked.rds, objects/coresh_derived_sets.rds,
## objects/02_de_results.rds (ranked vectors), master/master_gsea_table.csv rows with
## database="CoReSh_derived", and 08_coresh/tables/_overview/coresh_dataset_annotation.csv
## (dataset identity from 08b). Never re-runs coresh_batch/build_coresh_gmt/fgsea.
##
## GATED: the CoReSh arm requires the ~20 GB mmu Synapse compendium (syn66227307), consumed
## read-only from the shared reference cache, before 07/08 can be run. This script GUARDS on
## the first compute output and stops cleanly (status 0) if the arm has not been run.
##
## Figures produced
## ----------------
## _overview/ (CoReSh-specific — no analogue in the standard GSEA sweep):
##   coresh_pctvar_overview — faceted bars of top-ranked compendium datasets per query.
##   coresh_nes_dotplot     — CoReSh-derived set NES x contrast dotplot (query-origin strip).
##
## by_contrast/<c>/CoReSh_derived/ (SAME four toolkit panels every standard database gets in
## the 06_gsea sweep, so CoReSh reads at parity with GO/Hallmark/KEGG/... — via a gseaResult
## reconstructed viz-side and the RNAseq-toolkit plotters):
##   dotplot / facet / barplot / running_sum
##
## Derived-set labels carry what each recovered public GEO dataset actually is (from 08b's
## annotation), so a reader sees "GSE89069 · viral infection · embryonic brain", not a raw id.
##
## Sample mapping is owner-confirmed (GSE329522 2x2 genotype x temperature); captions carry
## that stamp via sample_mapping_stamp().

# =============================================================================
# 0. Style contract + toolkit GSEA plotters (MANDATORY FIRST — no inline theme/hex)
# =============================================================================

source("02_analysis/helpers/figure_style.R")    # FIG_CFG, project_theme, save_figure,
                                                 # save_overview, style_series, contrast_path,
                                                 # overview_path, contrast_label, direction_cue, %||%
source("02_analysis/config/config.R")            # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS,
                                                 # DIR_MASTER, DIR_RESULTS, SPECIES, stage_dir,
                                                 # sample_mapping_stamp, RANK_METRIC,
                                                 # GSEA_FDR_CUTOFF, GSEA_MIN_SIZE/MAX_SIZE, %||%
source("02_analysis/helpers/de_gsea_helpers.R")  # load_de_results(), build_ranked_vector()

# RNAseq-toolkit GSEA plotters (format_pathway_names BEFORE the plotters that call it).
.RTK <- "01_modules/RNAseq-toolkit"
.GP  <- file.path(.RTK, "scripts", "GSEA", "GSEA_plotting")
source(file.path(.RTK, "scripts", "custom_minimal_theme.R"))  # custom_minimal_theme_with_grid()
source(file.path(.GP, "format_pathway_names.R"))              # format_pathway_name()
source(file.path(.GP, "gsea_plotting_utils.R"))               # smart_wrap(), get_db_plot_params()
source(file.path(.GP, "gsea_dotplot.R"))
source(file.path(.GP, "gsea_dotplot_facet.R"))
source(file.path(.GP, "gsea_barplot.R"))
source(file.path(.GP, "gsea_running_sum_plot.R"))

options(bitmapType = "cairo")     # headless PNG
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(tibble); library(readr)
  library(scales); library(patchwork); library(stringr)
  library(methods); library(enrichplot); library(DOSE)   # gseaResult S4 + deterministic ES curve
})
options(stringsAsFactors = FALSE)

## The toolkit plotters hardcode format_pathway_name(use_formatting=TRUE) on @result$Description,
## which is built for MSigDB-style IDs and title-cases accession tokens (GSE89069 -> "Gse89069").
## This stage only ever plots CoReSh sets whose Description is ALREADY a final annotated label,
## so shadow format_pathway_name with an identity pass in this script's scope: labels render
## verbatim, smart_wrap() still wraps long ones inside the plotters.
format_pathway_name <- function(text, use_formatting = TRUE, strip_prefix = TRUE) as.character(text)

# =============================================================================
# 1. Parameters (config, not hardcoded)
# =============================================================================

STAGE       <- "08_coresh"
SCRIPT_PATH <- "02_analysis/scripts/15_coresh_viz.R"

tbl_dir <- stage_dir(STAGE, "tables")
fig_dir <- stage_dir(STAGE, "figures")
BY_CONTRAST_DIR <- YAML_CONFIG$figures$by_contrast_dir %||% "by_contrast"
OVERVIEW_DIR    <- YAML_CONFIG$figures$overview_dir    %||% "_overview"

coresh_cfg <- YAML_CONFIG$coresh %||% list()
TOP_N   <- as.integer(coresh_cfg$top_n_hits        %||% 5L)   # top datasets / query (pctVar bars)
SHOWCAT <- as.integer(YAML_CONFIG$figures$top_pathways %||% YAML_CONFIG$figures$top_n %||% 20L)
FDR     <- as.numeric(GSEA_FDR_CUTOFF %||% 0.05)
NES_CAP <- as.numeric(YAML_CONFIG$figures$nes_cap  %||% 3.5)
RSTOP   <- as.integer(FIG_CFG$figures$running_sum_top %||% 5L)
RSYLIM  <- as.numeric(unlist(FIG_CFG$figures$running_sum_ylim %||% c(-1, 1)))
MINGS   <- as.integer(GSEA_MIN_SIZE %||% 15L)
MAXGS   <- as.integer(GSEA_MAX_SIZE %||% 500L)
CAP     <- sample_mapping_stamp()     # owner-confirmed sample-mapping stamp
DB_NAME <- "CoReSh_derived"          # the single "database" this stage contributes

qsigs_cfg <- coresh_cfg$query_signatures %||% list()

## Colours from FIG_CFG (never inline hex)
NEG <- FIG_CFG$colors$diverging$down    %||% "steelblue4"
MID <- FIG_CFG$colors$diverging$neutral %||% "grey97"
POS <- FIG_CFG$colors$diverging$up      %||% "sienna"
OI  <- unlist(FIG_CFG$colors$okabe_ito  %||% list(), use.names = FALSE)
if (length(OI) == 0) OI <- c("darkorange", "deepskyblue3", "mediumseagreen", "gold",
                             "steelblue", "firebrick2", "orchid3", "black")

CONTRASTS <- vapply(YAML_CONFIG$design$contrasts, function(x) x$name %||% "", character(1))
CONTRASTS <- CONTRASTS[nzchar(CONTRASTS)]

# =============================================================================
# 2. GUARD — compute outputs must exist (NEVER fabricate)
# =============================================================================

rk_fp <- file.path(DIR_OBJECTS, "coresh_ranked.rds")
if (!file.exists(rk_fp)) {
  message("CoReSh outputs not found — run 07/08 after provisioning the mmu compendium; ",
          "exiting cleanly with no figures produced.")
  quit(save = "no", status = 0)
}

# =============================================================================
# 3. Query specs + dataset annotation -> per-set display labels & origins
# =============================================================================

.query_specs <- function(cfg) {
  specs <- cfg$sets %||% list()
  if (is.null(cfg$object) || !nzchar(as.character(cfg$object)))
    stop("15_coresh_viz: coresh.query_signatures.object is required.")
  if (!is.list(specs) || length(specs) == 0L)
    stop("15_coresh_viz: coresh.query_signatures.sets must be a non-empty list.")
  out <- lapply(seq_along(specs), function(i) {
    spec <- specs[[i]]
    required <- c("contrast", "direction", "gate")
    missing <- required[!nzchar(vapply(required, function(k) spec[[k]] %||% "", character(1)))]
    if (length(missing) > 0L)
      stop("15_coresh_viz: coresh.query_signatures.sets[[", i, "]] missing: ",
           paste(missing, collapse = ", "))
    data.frame(contrast = as.character(spec$contrast),
               direction = as.character(spec$direction),
               gate = as.character(spec$gate),
               origin = as.character(spec$origin %||% paste0(spec$contrast, "_", spec$direction)),
               stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(out) %>%
    mutate(query_name = paste0("Q_sig_", contrast, "_", direction, "_", gate))
}
QUERY_SPECS <- .query_specs(qsigs_cfg)
if (anyDuplicated(QUERY_SPECS$query_name))
  stop("15_coresh_viz: duplicate query name(s) in coresh.query_signatures.")
EXPECTED_QUERIES <- QUERY_SPECS$query_name

## Origin label per distinct (contrast, direction) — data-driven, not row-order.
ORIGIN_SPECS <- QUERY_SPECS %>%
  distinct(contrast, direction, .keep_all = TRUE) %>%
  transmute(query_prefix = paste0("Q_sig_", contrast, "_", direction, "_"),
            origin_label = origin)

## Dataset annotation (from 08b). Absent -> labels fall back to bare GSE (no hard fail).
ann_fp <- file.path(tbl_dir, OVERVIEW_DIR, "coresh_dataset_annotation.csv")
ANN <- if (file.exists(ann_fp)) readr::read_csv(ann_fp, show_col_types = FALSE) else NULL
if (is.null(ANN)) message("[15] NOTE: dataset annotation absent (run 08b) — using bare GSE labels.")

## Parse the trailing GSE from a derived-set id "CORESH_<query>_<GSE>".
.gse_of <- function(sid) sub(".*_(GSE[0-9]+)$", "\\1", sid)

## Build a readable, format-safe label for a derived-set id (FULL text, never truncated):
##   "GSE89069 · viral infection · embryonic brain (E17.5)"  (context + tissue from 08b)
## Long labels are shown in full by wrapping to multiple lines (.wrap / plotter wrap_width),
## not by cutting them.
.set_label <- function(sid) {
  gse <- .gse_of(sid)
  if (!is.null(ANN)) {
    row <- ANN[ANN$gse == gse, , drop = FALSE]
    if (nrow(row) == 1L) {
      ctx <- gsub("_", " ", row$context_class[1])
      tis <- row$tissue_cell[1]
      tis <- if (!is.na(tis) && nzchar(tis)) tis else NULL
      return(paste(c(gse, ctx, tis), collapse = " · "))
    }
  }
  gse
}

## Wrap a label to multiple lines so long names display in full (used for axis/legend text
## that the toolkit plotters do not wrap themselves).
WRAP_W <- 40L
.wrap <- function(x) vapply(x, function(s) paste(strwrap(s, width = WRAP_W), collapse = "\n"),
                            character(1), USE.NAMES = FALSE)

## Query-origin family for a derived-set id (for the nes_dotplot colour strip).
.query_origin <- function(sid) {
  for (i in seq_len(nrow(ORIGIN_SPECS))) {
    if (grepl(ORIGIN_SPECS$query_prefix[i], sid, fixed = TRUE))
      return(ORIGIN_SPECS$origin_label[i])
  }
  "Other CoReSh query"
}

# =============================================================================
# 4. Load the compendium ranking (pctVar sweep)
# =============================================================================

ranked <- as.data.frame(readRDS(rk_fp), stringsAsFactors = FALSE)
stopifnot(all(c("query_name", "gse", "pctVar", "size", "rank") %in% colnames(ranked)))
if (any(EXPECTED_QUERIES %in% ranked$query_name))
  ranked <- ranked[ranked$query_name %in% EXPECTED_QUERIES, , drop = FALSE]
message(sprintf("[15] ranked compendium: %d rows, %d queries, %d unique GSEs",
                nrow(ranked), length(unique(ranked$query_name)), length(unique(ranked$gse))))

.query_label <- function(q) gsub("_", " ", sub("^Q_sig_", "Signature: ", q))

# =============================================================================
# 5. FIGURE — pctVar overview (CoReSh-specific: top compendium datasets per query)
# =============================================================================

create_pctvar_overview <- function(ranked, top_n = TOP_N) {
  if (nrow(ranked) == 0L) return(invisible(NULL))
  ## Compact facet label, e.g. "Interaction ↑ (fdr_logfc)" — the long "Signature:" form
  ## clips inside the free-scale panels.
  .pv_facet <- function(q) {
    gate <- sub("^Q_sig_.*_(fdr_logfc|fdr_only)$", "\\1", q)
    base <- sub("^Q_sig_(.*)_up_(fdr_logfc|fdr_only)$", "\\1", q)
    sprintf("%s ↑ (%s)", base, gate)
  }
  df <- ranked %>%
    group_by(query_name) %>% slice_max(pctVar, n = top_n, with_ties = FALSE) %>% ungroup() %>%
    mutate(query_label = .pv_facet(query_name),
           gse_lab = sprintf("%s (k=%d)", gse, size))
  if (nrow(df) == 0L) return(invisible(NULL))
  ## The GSE accession IS the y-axis tick (no floating labels to clip). gse_lab carries k
  ## and is unique per facet (k differs across gates), so a single factor ordered by
  ## (query, pctVar) lands each facet's bars in rank order — highest pctVar at top.
  df <- df %>% arrange(query_name, pctVar) %>%
    mutate(gse_lab = factor(gse_lab, levels = unique(gse_lab)))

  ggplot(df, aes(x = pctVar, y = gse_lab)) +
    geom_col(width = 0.7, alpha = 0.88) +
    facet_wrap(~ query_label, scales = "free") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.06)),
                       labels = function(x) paste0(round(x, 0), "%")) +
    labs(title = "CoReSh compendium ranking — project signature queries",
         subtitle = sprintf("top %d public mouse GEO datasets per query by pctVar (rank 1 = top of each panel)", top_n),
         x = "Variance explained by query signature in public GEO dataset (pctVar, %)", y = NULL,
         caption = "Each bar = one GEO dataset · y = GSE accession (k = mapped query genes) · length = pctVar\npctVar = PCA-inspired co-regulation score; higher = query set co-moves more across that dataset") +
    project_theme(config = FIG_CFG)
}

p_pctvar <- tryCatch(create_pctvar_overview(ranked),
                     error = function(e) { message("  pctvar_overview: ", e$message); NULL })
top_ds_ov <- if (nrow(ranked) > 0L) {
  ranked %>% group_by(query_name) %>% slice_max(pctVar, n = TOP_N, with_ties = FALSE) %>% ungroup() %>%
    arrange(query_name, rank) %>% mutate(query_label = .query_label(query_name)) %>%
    select(query_name, query_label, gse, gpl, pctVar, pval, size, rank)
} else NULL

if (!is.null(p_pctvar) && !is.null(top_ds_ov) && nrow(top_ds_ov) > 0L) {
  tryCatch(save_overview(
    plot = p_pctvar, stage = STAGE, name = "coresh_pctvar_overview", table = top_ds_ov,
    finding = "Top public mouse GEO datasets co-regulating each project signature query, ranked by CoReSh pctVar.",
    script = SCRIPT_PATH, fn = "create_pctvar_overview",
    config_kv = sprintf("coresh.top_n_hits=%d; figures.top_pathways=%d; thresholds.gsea_fdr=%.2f", TOP_N, SHOWCAT, FDR),
    input = "03_results/objects/coresh_ranked.rds",
    how_to_read = paste0(
      "Each bar = one public GEO dataset. Bar length = pctVar (% of that dataset's variance explained by the query), ",
      "a PCA-inspired co-regulation score (higher = stronger co-regulation); label = GSE, mapped query size k, rank. ",
      "Facets separate the project signature queries. Claim tier: L3-DE (compendium co-regulation score). ",
      "pctVar >= 0 (variance fraction, unsigned)."),
    config = FIG_CFG),
    error = function(e) message("  coresh_pctvar_overview: ", e$message))
}

# =============================================================================
# 6. Reconstruct a gseaResult for CoReSh_derived (the toolkit-plotter input)
# =============================================================================
## Mirrors the standard sweep's viz-side reconstruction: @result taken VERBATIM from
## master_gsea_table.csv (exact NES/padj — zero engine drift), @geneList the ranked DE
## vector, @geneSets the CoReSh-derived sets restricted to the ranked universe, and
## @params$exponent = 1 so the running-ES is deterministic (no permutation re-run).

sets_fp <- file.path(DIR_OBJECTS, "coresh_derived_sets.rds")
CORESH_SETS <- if (file.exists(sets_fp)) readRDS(sets_fp) else NULL
DE_RESULTS  <- tryCatch(load_de_results(), error = function(e) NULL)

mg_fp <- file.path(DIR_MASTER, "master_gsea_table.csv")
MG <- if (file.exists(mg_fp)) readr::read_csv(mg_fp, show_col_types = FALSE) else NULL
MG_CORESH <- if (!is.null(MG)) MG[!is.na(MG$database) & MG$database == DB_NAME, , drop = FALSE] else NULL

.RANKED_CACHE <- new.env(parent = emptyenv())
ranked_for <- function(co) {
  if (is.null(.RANKED_CACHE[[co]])) {
    tt <- DE_RESULTS[[co]]
    if (is.null(tt)) return(NULL)
    .RANKED_CACHE[[co]] <- build_ranked_vector(tt, RANK_METRIC)
  }
  .RANKED_CACHE[[co]]
}

as_coresh_gsearesult <- function(co) {
  if (is.null(MG_CORESH) || is.null(CORESH_SETS)) return(NULL)
  rows <- MG_CORESH %>% dplyr::filter(contrast == co, !is.na(nes))
  if (nrow(rows) == 0L) return(NULL)
  rk <- ranked_for(co); if (is.null(rk)) return(NULL)

  res <- data.frame(
    ID = rows$pathway_id,
    Description = vapply(rows$pathway_id, .set_label, character(1)),  # final annotated label
    setSize = as.integer(rows$set_size),
    enrichmentScore = NA_real_,
    NES = as.numeric(rows$nes),
    pvalue = as.numeric(rows$pvalue),
    p.adjust = as.numeric(rows$padj),
    qvalue = as.numeric(rows$padj),
    rank = NA_integer_, leading_edge = NA_character_,
    core_enrichment = as.character(rows$core_enrichment),
    stringsAsFactors = FALSE)
  rownames(res) <- res$ID

  sets <- lapply(CORESH_SETS, function(g) intersect(as.character(g), names(rk)))
  common <- intersect(res$ID, names(sets)); sets <- sets[common]

  methods::new("gseaResult",
    result = res, organism = (if (exists("SPECIES")) SPECIES else "Mus musculus"),
    setType = "CoReSh_derived", geneSets = sets, geneList = rk, keytype = "SYMBOL",
    permScores = matrix(),
    params = list(pvalueCutoff = 1, eps = 0, pAdjustMethod = "fdr",
                  exponent = 1, minGSSize = MINGS, maxGSSize = MAXGS),
    readable = FALSE)
}

# =============================================================================
# 7. Per-contrast toolkit battery: dotplot / facet / barplot / running_sum
#    -> by_contrast/<c>/CoReSh_derived/  (parity with every standard database)
# =============================================================================

emit_coresh_cell <- function(co) {
  g <- tryCatch(as_coresh_gsearesult(co),
                error = function(e) { message(sprintf("  [15] %s build error: %s", co, conditionMessage(e))); NULL })
  if (is.null(g) || nrow(g@result) == 0L) {
    message(sprintf("  [15] SKIP (%s): no gseaResult (0 CoReSh_derived rows or no gene sets)", co)); return(invisible(NULL))
  }
  n_sig <- sum(g@result$p.adjust < FDR, na.rm = TRUE)
  message(sprintf("  [15] %s / %s: %d sets, %d sig (padj < %.2g)", co, DB_NAME, nrow(g@result), n_sig, FDR))

  ttl <- sprintf("%s — %s", DB_NAME, contrast_label(co, short = TRUE))
  sub <- "NES > 0 = up in numerator, NES < 0 = down · sample mapping owner-confirmed"

  fig_db_dir <- file.path(contrast_path(STAGE, co, "figures", config = FIG_CFG), DB_NAME)
  dir.create(fig_db_dir, recursive = TRUE, showWarnings = FALSE)
  stale <- list.files(fig_db_dir, pattern = "\\.(png|pdf)$", full.names = TRUE)
  if (length(stale)) file.remove(stale)

  # 1. dotplot
  p_dot <- gsea_dotplot(g, filterBy = "p.adjust", showCategory = SHOWCAT, padj_cutoff = FDR,
                        title = ttl, neg_color = NEG, mid_color = MID, pos_color = POS,
                        nes_limits = c(-NES_CAP, NES_CAP), strip_prefix = FALSE) +
    labs(subtitle = sub) + project_theme(config = FIG_CFG)
  if (n_sig == 0)
    p_dot <- p_dot + annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
             label = sprintf("No significant sets (FDR < %.2g)", FDR),
             size = (FIG_CFG$figures$label_size %||% 5) + 1, fontface = "bold",
             colour = "white", fill = NEG, label.size = 0)
  save_figure(p_dot, STAGE, file.path(DB_NAME, "dotplot"), contrast = co, config = FIG_CFG, height = 9)

  # 2. up/down faceted dotplot
  p_fac <- gsea_dotplot_facet(g, showCategory = SHOWCAT, padj_cutoff = FDR, title = ttl,
                              neg_color = NEG, mid_color = MID, pos_color = POS,
                              nes_limits = c(-NES_CAP, NES_CAP), strip_prefix = FALSE) +
    labs(subtitle = sub) + project_theme(config = FIG_CFG)
  save_figure(p_fac, STAGE, file.path(DB_NAME, "facet"), contrast = co, config = FIG_CFG, height = 9, wide = TRUE)

  # 3. NES barplot (significant sets only, top-N by |NES|)
  p_bar <- gsea_barplot(g, padj_cutoff = FDR, top_n = SHOWCAT, title = ttl,
                        neg_color = NEG, mid_color = MID, pos_color = POS,
                        nes_limits = c(-NES_CAP, NES_CAP), strip_prefix = FALSE) +
    labs(subtitle = sub) + project_theme(config = FIG_CFG)
  save_figure(p_bar, STAGE, file.path(DB_NAME, "barplot"), contrast = co, config = FIG_CFG, height = 9)

  # 4. running-sum — REAL three-panel ES curves (Description already = annotated label)
  have_set <- g@result$ID %in% names(g@geneSets)
  cand <- g@result[have_set, , drop = FALSE]
  if (nrow(cand) > 0L) {
    top_ids <- cand$ID[order(abs(cand$NES), decreasing = TRUE)][seq_len(min(RSTOP, nrow(cand)))]
    g_rs <- g                                   # wrap legend labels so long names show in full
    g_rs@result$Description <- .wrap(g_rs@result$Description)
    p_rs <- tryCatch(gsea_running_sum_plot(g_rs, gene_set_ids = top_ids, title = ttl, max_name_length = 200),
                     error = function(e) { message(sprintf("  [15] running_sum skipped (%s): %s", co, conditionMessage(e))); NULL })
    if (!is.null(p_rs)) {
      # x becomes rank/n_ranked, taken from this contrast's own ranked vector.
      #
      # Each curve also gets its own named tick row, and the plotter stacks the rows from the
      # bottom up in the order the ids arrive, so the labels go in that order. The row key is
      # the source accession, which is the one token that separates these sets from each other
      # and is the tail of the legend entry beside it; the rest of the derived-set name is the
      # same words on every curve and would cost the panel width it does not repay.
      rs_labs <- sub("^.*_(GSE[0-9]+)$", "\\1", top_ids)
      p_rs <- style_series(p_rs, ylim = RSYLIM, n_ranked = length(g_rs@geneList),
                           trace_labels = rs_labs, config = FIG_CFG)
      save_figure(p_rs, STAGE, file.path(DB_NAME, "running_sum"), contrast = co, config = FIG_CFG, wide = TRUE)
    }
  } else {
    message(sprintf("  [15] running_sum skipped (%s): no set has a gene-set list", co))
  }
  invisible(TRUE)
}

message(sprintf("[15] Generating CoReSh_derived battery: %d contrasts x {dotplot,facet,barplot,running_sum}",
                length(CONTRASTS)))
for (co in CONTRASTS) emit_coresh_cell(co)

# =============================================================================
# 8. FIGURE — cross-contrast NES dotplot (CoReSh-derived sets x contrasts, _overview/)
# =============================================================================
## Pooled dotplot: y = derived set (annotated), x = contrast, fill = NES,
## size = -log10(padj), outline = FDR sig; left strip = query origin.

pool <- if (!is.null(MG_CORESH) && nrow(MG_CORESH) > 0L) {
  tibble::tibble(set_id = MG_CORESH$pathway_id,
                 nes = MG_CORESH$nes, padj = MG_CORESH$padj,
                 set_size = MG_CORESH$set_size %||% NA_integer_, contrast = MG_CORESH$contrast)
} else NULL

create_nes_dotplot <- function(pool, top_n = SHOWCAT) {
  if (is.null(pool) || nrow(pool) == 0L) return(invisible(NULL))
  rank_sets <- pool %>% group_by(set_id) %>%
    summarise(m = max(abs(nes), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(m)) %>% slice_head(n = top_n) %>% pull(set_id)
  df <- pool %>% filter(set_id %in% rank_sets) %>%
    mutate(set_lab = .wrap(vapply(set_id, .set_label, character(1))),
           origin  = vapply(set_id, .query_origin, character(1)),
           sig = !is.na(padj) & padj < FDR,
           neglogp = pmin(-log10(pmax(padj, 1e-50)), 16),
           nes_capped = pmax(pmin(nes, NES_CAP), -NES_CAP),
           contrast = factor(contrast, levels = intersect(CONTRASTS, unique(contrast))))
  ## Set order: highest median NES at top (within each origin band via space="free_y").
  set_order <- df %>% group_by(set_lab) %>%
    summarise(med = stats::median(nes, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med)) %>% pull(set_lab)
  df$set_lab <- factor(df$set_lab, levels = rev(set_order))
  df$origin  <- factor(df$origin, levels = c(ORIGIN_SPECS$origin_label, "Other CoReSh query"))
  ## x-axis: short contrast labels (flattened to one line) in biological (config) order.
  flat <- function(x) gsub("\\s*\n\\s*", " ", contrast_label(x, short = TRUE))
  xlev <- flat(levels(df$contrast))
  df$contrast_lab <- factor(flat(as.character(df$contrast)), levels = xlev)

  ## Query origin is a left-side facet band (single panel; handles the long dataset
  ## labels the patchwork strip could not without pushing the matrix off-page).
  ggplot(df, aes(x = contrast_lab, y = set_lab)) +
    geom_point(aes(size = neglogp, fill = nes_capped, color = sig), shape = 21, stroke = 0.8) +
    facet_grid(rows = vars(origin), scales = "free_y", space = "free_y", switch = "y",
               labeller = labeller(origin = label_wrap_gen(16))) +
    scale_fill_gradient2(low = NEG, mid = MID, high = POS, midpoint = 0, name = "NES",
                         limits = c(-NES_CAP, NES_CAP), oob = scales::squish) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "transparent"), name = sprintf("FDR < %.2f", FDR)) +
    scale_size_continuous(range = c(2.5, 9), name = "−log₁₀(padj)") +
    labs(title = "CoReSh-derived set enrichment across contrasts",
         subtitle = sprintf("fill = NES (orange up / blue down) · size = −log₁₀(padj) · outline = FDR<%.2f · rows = query origin", FDR),
         x = NULL, y = NULL) +
    project_theme(config = FIG_CFG) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.placement = "outside")
}

p_nes <- tryCatch(create_nes_dotplot(pool), error = function(e) { message("  nes_dotplot: ", e$message); NULL })
nes_tbl <- if (!is.null(pool) && nrow(pool) > 0L) {
  pool %>% mutate(set_label = vapply(set_id, .set_label, character(1)),
                  origin = vapply(set_id, .query_origin, character(1)),
                  sig = !is.na(padj) & padj < FDR) %>%
    arrange(contrast, padj) %>% select(set_id, set_label, origin, contrast, nes, padj, set_size, sig)
} else NULL

if (!is.null(p_nes)) {
  tryCatch(save_overview(
    plot = p_nes, stage = STAGE, name = "coresh_nes_dotplot", table = nes_tbl,
    finding = "Cross-contrast enrichment of CoReSh-derived co-regulation sets (NES fill, -log10(padj) size, FDR outline); rows grouped by seeding query origin; labels carry each recovered dataset's identity.",
    script = SCRIPT_PATH, fn = "create_nes_dotplot",
    config_kv = sprintf("figures.top_pathways=%d; figures.nes_cap=%.1f; thresholds.gsea_fdr=%.2f", SHOWCAT, NES_CAP, FDR),
    input = "03_results/master/master_gsea_table.csv rows database=CoReSh_derived",
    how_to_read = paste0(
      "Each circle = one CoReSh-derived gene set (row) in one contrast (column). Row label = the recovered public ",
      "GEO dataset (GSE · context · tissue, from the dataset annotation). Fill = NES: orange positive (up in numerator), ",
      "blue negative; clamped at ±", NES_CAP, ". Size = -log10(padj). Black outline = FDR < ", FDR, ". ",
      "Rows are grouped into left-side facet bands by the project signature query that seeded each set; ",
      "within a band, sets are ordered by median NES (highest at top). Claim tier: L3-DE (fgsea, BH-FDR)."),
    config = FIG_CFG, width = 10, height = 11),
    error = function(e) message("  coresh_nes_dotplot: ", e$message))
}

# =============================================================================
# 9. Stage README header (write-if-absent; full captions live per-figure)
# =============================================================================

readme_fp <- file.path(DIR_RESULTS, STAGE, "README.md")
if (!file.exists(readme_fp)) {
  dir.create(dirname(readme_fp), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    sprintf("# %s — CoReSh compendium ranking + derived-set GSEA", STAGE), "",
    sprintf("Stage: `%s` | Scripts: `07_coresh_search.R` (search) · `08_coresh_derived_gsea.R` (GSEA) · `08b_coresh_annotate.R` (dataset annotation) · `%s` (viz)", STAGE, SCRIPT_PATH), "",
    "**What CoReSh does here.** It tests where the project's mouse-native signatures co-regulate across the public mouse mmu GEO compendium, turns the top public datasets into derived gene sets, and runs those sets through the same GSEA as every other database.", "",
    "**Query sets used.** Four project signatures seed the search: `WT_heat_up` (generic thermal, fdr_logfc + fdr_only gates) and `Interaction_up` (cGAS-dependent, fdr_logfc + fdr_only gates), exported from the signature-derivation stage. The resulting `CoReSh_derived` sets are data-driven co-regulation modules that do not exist in any curated collection (Hallmark/GO/KEGG/Reactome/...), so they are unique to this stage and populate `08_coresh` rather than the standard GSEA sweep.", "",
    "**Dataset annotation.** Each recovered dataset is labelled with what it actually is (title/organism/tissue/perturbation + context flags) in `tables/_overview/coresh_dataset_annotation.csv`; figure labels carry that identity.", "",
    "**Sample mapping is owner-confirmed** (GSE329522 iTreg 2×2 genotype × temperature).", "",
    "**GATED:** requires the ~20 GB mmu Synapse compendium (syn66227307), mounted read-only, before `07` can run.", "", "---", ""
  ), readme_fp)
  message(sprintf("[15] wrote stage README: %s", readme_fp))
} else {
  message(sprintf("[15] stage README exists: %s", readme_fp))
}

# =============================================================================
# 10. Uniqueness check — CoReSh_derived sets are disjoint from the GSEA collections
# =============================================================================
## The CoReSh-derived sets are data-driven co-regulation modules; by construction they
## carry the "CORESH_" id prefix and database="CoReSh_derived". Assert their pathway_ids
## share nothing with any curated collection in the master table, so it is explicit that
## this stage contributes UNIQUE sets (they populate 08_coresh, not the 06_gsea sweep).

if (!is.null(MG) && !is.null(MG_CORESH) && nrow(MG_CORESH) > 0L) {
  coresh_ids <- unique(MG_CORESH$pathway_id)
  other_ids  <- unique(MG$pathway_id[MG$database != DB_NAME])
  overlap <- intersect(coresh_ids, other_ids)
  if (length(overlap) > 0L)
    warning(sprintf("[15] uniqueness: %d CoReSh_derived id(s) also appear under another database: %s",
                    length(overlap), paste(utils::head(overlap, 5), collapse = ", ")))
  else
    message(sprintf("[15] uniqueness OK: all %d CoReSh_derived set ids are disjoint from the %d GSEA-collection ids.",
                    length(coresh_ids), length(other_ids)))
}

# =============================================================================
# 11. Structural asserts (fail loudly; overviews need a same-stem sidecar table)
# =============================================================================

ov_fig_dir <- overview_path(STAGE, "figures", config = FIG_CFG)
ov_tbl_dir <- overview_path(STAGE, "tables",  config = FIG_CFG)
ov_stems <- unique(sub("\\.(pdf|png)$", "", list.files(ov_fig_dir, pattern = "\\.(pdf|png)$")))
ov_orphans <- ov_stems[!file.exists(file.path(ov_tbl_dir, paste0(ov_stems, ".csv")))]
if (length(ov_orphans) > 0L)
  warning("Overview figure(s) without a same-stem sidecar CSV: ", paste(ov_orphans, collapse = ", "))

n_battery <- sum(vapply(CONTRASTS, function(co)
  file.exists(file.path(contrast_path(STAGE, co, "figures", config = FIG_CFG), DB_NAME, "dotplot.png")),
  logical(1)))
message(sprintf("[15] COMPLETE — %d _overview stem(s), %d/%d contrasts with a CoReSh_derived battery. | %s",
                length(ov_stems), n_battery, length(CONTRASTS), CAP))
