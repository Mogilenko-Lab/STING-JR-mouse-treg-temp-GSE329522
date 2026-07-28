# 12_gsea_viz.R — VIZ
## GSEA visualization for the STING-cGAS standard sweep (stage 06_gsea).
##
## Per (contrast x database) cell the four canonical toolkit panels:
##   running_sum.{print.pdf,screen.png}  — REAL three-panel GSEA enrichment curve
##                                         (running ES / gene-hit ticks / ranked metric)
##                                         via enrichplot::gseaplot2() through the toolkit
##                                         gsea_running_sum_plot() (NOT a lollipop proxy).
##   dotplot.{print.pdf,screen.png}      — gsea_dotplot()        (x = GeneRatio, fill = NES)
##   facet.{print.pdf,screen.png}        — gsea_dotplot_facet()  (up/down split)
##   barplot.{print.pdf,screen.png}      — gsea_barplot()        (NES bars, significant only)
##
## Cross-contrast overview (03_results/06_gsea/figures/_overview/):
##   gsea_asymmetry_panel           — cGAS-dependence asymmetry (IFN/ISG vs HIF/glycolysis)
##   gsea_lombardi_vs_bespoke_hif   — Lombardi-2022 48-gene vs Hallmark Hypoxia HIF check
##
## Run from project root AFTER 05_gsea_msigdb_run.R and 06_gsea_custom_run.R:
##   Rscript 02_analysis/scripts/12_gsea_viz.R
##
## VIZ-ONLY. Reads master_gsea_table.csv + cached DE topTables (02_de_results.rds) +
## cached gene-set lists (geneset_*_<DB>.rds). NEVER re-runs GSEA, never re-fits DE,
## never writes 03_results/master/.
##
## How the toolkit plotters are fed (the headline fix)
## ---------------------------------------------------
## The toolkit GSEA plotters require a clusterProfiler `gseaResult` S4 object whose
## @result carries ID/Description/NES/p.adjust/setSize/core_enrichment AND whose
## @geneList (the ranked stat vector) + @geneSets (the gene-set membership lists)
## drive the running-sum ES curve. The compute master CSV is correct but tidy; it
## discards the S4 + the ranked vector. We REBUILD the S4 viz-side from inputs that
## already exist on disk (as_gsearesult()):
##   * @result   <- the master_gsea_table.csv rows for this (contrast x DB)  -> the
##                  *quotable* NES/padj are reproduced EXACTLY (zero engine drift).
##   * @geneList <- build_ranked_vector(02_de_results.rds[[contrast]], "t") -> ranked
##                  decreasing, gene-symbol named (the exact compute ranking input).
##   * @geneSets <- the cached geneset_<src>_<DB>.rds list (gs_name -> symbols),
##                  restricted to the ranked universe.
##   * @params$exponent = 1 -> the running ES is computed DETERMINISTICALLY by
##                  enrichplot::gseaScores(geneList, geneSet, exponent); NO permutation
##                  re-run, so the curve never disagrees with the master NES.
## The fgsea(compute) vs clusterProfiler(plot) NES engines were checked to agree on the
## one cached gseaResult (gsea_msigdb_WT_heat.rds): max |dNES| = 0.04, r = 0.99999,
## 100%% sign agreement (logged below). Because @result is taken verbatim from the
## master, the figures carry the master's own statistics regardless.
##
## House constraints (non-negotiable in all captions):
##   * Never "cGAS-independent" -> say "no detectable cGAS-dependence at n=5"
##   * Never crown HIF1alpha/HIF2alpha as the driver
##   * Lombardi-2022 (48-gene) is the orthogonal HIF check; contrasted with hand-made 16-gene
##   * Claims floor at L3 (DE/enrichment stats); stamp PROVISIONAL sample labels where relevant
##   * NES sign: positive = enriched in NUMERATOR (hot/WT side); negative = denominator

# ============================================================================
# 0. Setup + style contract + RNAseq-toolkit GSEA plotters
# ============================================================================

source("02_analysis/helpers/figure_style.R")    # -> FIG_CFG, project_theme(), save_figure(),
                                                 #    save_overview(), style_series(), contrast_path(),
                                                 #    overview_path(), purge_figures(), write_caption()
source("02_analysis/config/config.R")            # -> YAML_CONFIG, DIR_OBJECTS, DIR_MASTER,
                                                 #    GSEA_FDR_CUTOFF, GSEA_MIN_SIZE/MAX_SIZE, RANK_METRIC
source("02_analysis/helpers/de_gsea_helpers.R")  # -> load_de_results(), build_ranked_vector()

# RNAseq-toolkit GSEA plotters (source format_pathway_names BEFORE the plotters that call it).
.RTK <- "01_modules/RNAseq-toolkit"
.GP  <- file.path(.RTK, "scripts", "GSEA", "GSEA_plotting")
source(file.path(.RTK, "scripts", "custom_minimal_theme.R"))  # custom_minimal_theme_with_grid()
source(file.path(.GP, "format_pathway_names.R"))              # format_pathway_name()  <- BEFORE plotters
source(file.path(.GP, "gsea_plotting_utils.R"))               # smart_wrap(), get_db_plot_params()
source(file.path(.GP, "gsea_dotplot.R"))                      # gsea_dotplot()
source(file.path(.GP, "gsea_dotplot_facet.R"))                # gsea_dotplot_facet()
source(file.path(.GP, "gsea_barplot.R"))                      # gsea_barplot()
source(file.path(.GP, "gsea_running_sum_plot.R"))             # gsea_running_sum_plot()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(methods)
  library(enrichplot)   # gseaplot2() + gseaScores() (running-ES, deterministic)
  library(DOSE)         # gseaResult S4 class definition
})
options(stringsAsFactors = FALSE)

# ============================================================================
# 1. Constants from config (NEVER hardcoded)
# ============================================================================

STAGE    <- "06_gsea"
SCRIPT   <- "02_analysis/scripts/12_gsea_viz.R"
PID      <- YAML_CONFIG$project$id

FDR      <- GSEA_FDR_CUTOFF                                    # 0.05
TOPN     <- FIG_CFG$figures$top_pathways   %||% 20L            # toolkit showCategory / top_n
RSTOP    <- FIG_CFG$figures$running_sum_top %||% 5L            # # curves per running-sum panel
RSYLIM   <- as.numeric(unlist(FIG_CFG$figures$running_sum_ylim %||% c(-1, 1)))  # shared ES y-range
NESCAP   <- FIG_CFG$figures$nes_cap        %||% 3.5
MINGS    <- GSEA_MIN_SIZE                                      # 15
MAXGS    <- GSEA_MAX_SIZE                                      # 500

NEG      <- FIG_CFG$colors$diverging$down       %||% "steelblue4"    # blue  (NES < 0)
MID      <- FIG_CFG$colors$diverging$neutral    %||% "grey97"        # white (NES ~ 0)
POS      <- FIG_CFG$colors$diverging$up         %||% "sienna"        # orange (NES > 0)
COL_PURPLE <- FIG_CFG$colors$okabe_ito$reddish_purple %||% "orchid3" # Hallmark Hypoxia trace

CONTRASTS     <- vapply(YAML_CONFIG$design$contrasts, function(x) x$name, character(1))

MSIGDB_NAMES <- vapply(YAML_CONFIG$databases$msigdb, function(d) d$name, character(1))
CUSTOM_NAMES <- vapply(YAML_CONFIG$databases$custom,  function(d) d$name, character(1))
ALL_DBS      <- c(MSIGDB_NAMES, CUSTOM_NAMES)

# ============================================================================
# 2. Guard: master table + DE checkpoint must exist (compute scripts ran first)
# ============================================================================

mg_fp <- file.path(DIR_MASTER, "master_gsea_table.csv")
if (!file.exists(mg_fp)) {
  stop("12_gsea_viz: master_gsea_table.csv not found at ", mg_fp,
       "\n  Run 05_gsea_msigdb_run.R and 06_gsea_custom_run.R first.")
}

mg <- readr::read_csv(mg_fp, show_col_types = FALSE)
req_cols <- c("pathway_id", "pathway_name", "database", "nes", "pvalue",
              "padj", "set_size", "core_enrichment", "contrast", "direction")
missing_cols <- setdiff(req_cols, colnames(mg))
if (length(missing_cols) > 0)
  stop("12_gsea_viz: master_gsea_table.csv is missing required columns: ",
       paste(missing_cols, collapse = ", "))

message(sprintf("[12] Master GSEA table loaded: %d rows, %d databases, %d contrasts",
                nrow(mg), dplyr::n_distinct(mg$database), dplyr::n_distinct(mg$contrast)))

CONTRASTS_AVAIL <- intersect(CONTRASTS, unique(mg$contrast))
if (length(CONTRASTS_AVAIL) == 0)
  stop("12_gsea_viz: no config contrasts found in master_gsea_table.csv. ",
       "Config: ", paste(CONTRASTS, collapse = ","),
       " Master: ", paste(unique(mg$contrast), collapse = ","))
if (!identical(sort(CONTRASTS), sort(CONTRASTS_AVAIL)))
  message(sprintf("[12] NOTE: %d/%d config contrasts in master (missing: %s)",
                  length(CONTRASTS_AVAIL), length(CONTRASTS),
                  paste(setdiff(CONTRASTS, CONTRASTS_AVAIL), collapse = ", ")))
CONTRASTS <- CONTRASTS_AVAIL

# DE topTables (ranked-vector source) + gene-set manifest (gene-set lists).
DE_RESULTS <- load_de_results()                                       # named list of 7 topTables
GS_MANIFEST <- {
  mp <- file.path(DIR_OBJECTS, "geneset_manifest.csv")
  if (!file.exists(mp)) stop("12_gsea_viz: geneset_manifest.csv not found at ", mp,
                             " — run 05/06_gsea_*_run.R first.")
  readr::read_csv(mp, show_col_types = FALSE)
}

# ============================================================================
# 3. Per-(contrast x DB) gseaResult reconstruction (the toolkit-plotter input)
# ============================================================================

# Cache the ranked stat vector per contrast (build_ranked_vector is cheap but
# called once per DB; memoise across the inner DB loop).
.RANKED_CACHE <- new.env(parent = emptyenv())
ranked_for <- function(co) {
  if (is.null(.RANKED_CACHE[[co]])) {
    tt <- DE_RESULTS[[co]]
    if (is.null(tt)) stop("ranked_for: contrast '", co, "' absent from 02_de_results.rds")
    .RANKED_CACHE[[co]] <- build_ranked_vector(tt, RANK_METRIC)       # named, sorted decreasing
  }
  .RANKED_CACHE[[co]]
}

# Cache the gene-set list per DB (a named list gs_name -> symbol character vector).
.GENESET_CACHE <- new.env(parent = emptyenv())
geneset_for <- function(db) {
  if (is.null(.GENESET_CACHE[[db]])) {
    fp <- GS_MANIFEST$cache_path[GS_MANIFEST$database == db]
    if (length(fp) == 0 || is.na(fp[1]) || !file.exists(fp[1])) {
      .GENESET_CACHE[[db]] <- list(.absent = TRUE)                    # sentinel -> NULL on read
    } else {
      sets <- readRDS(fp[1])
      if (!is.list(sets) || is.null(names(sets))) sets <- list(.absent = TRUE)
      .GENESET_CACHE[[db]] <- sets
    }
  }
  out <- .GENESET_CACHE[[db]]
  if (length(out) == 1L && isTRUE(out$.absent)) return(NULL)
  out
}

#' Build a clusterProfiler `gseaResult` for one (contrast x DB) from cached inputs.
#'
#' @result is filled VERBATIM from master_gsea_table.csv (exact NES/padj/setSize/
#'   core_enrichment — the quotable statistics, no re-run). @geneList is the ranked
#'   vector, @geneSets the gene-set lists restricted to the ranked universe, and
#'   @params$exponent = 1 so the running-ES is computed deterministically by gseaScores.
#' Returns NULL when the cell has no master rows or no gene-set list on disk.
as_gsearesult <- function(co, db) {
  rows <- mg %>% dplyr::filter(contrast == co, database == db, !is.na(nes))
  if (nrow(rows) == 0) return(NULL)
  sets <- geneset_for(db)
  if (is.null(sets)) return(NULL)
  ranked <- ranked_for(co)

  # @result in the clusterProfiler gseaResult schema. Description = raw pathway_id so the
  # toolkit format_pathway_name()/strip_prefix renders the SAME display name everywhere.
  res <- data.frame(
    ID              = rows$pathway_id,
    Description     = rows$pathway_id,
    setSize         = as.integer(rows$set_size),
    enrichmentScore = NA_real_,
    NES             = as.numeric(rows$nes),
    pvalue          = as.numeric(rows$pvalue),
    p.adjust        = as.numeric(rows$padj),
    qvalue          = as.numeric(rows$padj),
    rank            = NA_integer_,
    leading_edge    = NA_character_,
    core_enrichment = as.character(rows$core_enrichment),
    stringsAsFactors = FALSE
  )
  rownames(res) <- res$ID

  # Restrict each gene set to the ranked universe; key by the @result IDs (running-sum
  # needs @geneSets[[ID]] for every plotted pathway). Sets absent from the cached list
  # (e.g. size-filtered out) simply won't be available for the running-sum, which is fine.
  sets <- lapply(sets, function(g) intersect(as.character(g), names(ranked)))
  common <- intersect(res$ID, names(sets))
  sets <- sets[common]

  methods::new("gseaResult",
    result      = res,
    organism    = (if (exists("SPECIES")) SPECIES else "Mus musculus"),
    setType     = "UNKNOWN",
    geneSets    = sets,
    geneList    = ranked,
    keytype     = "UNKNOWN",
    permScores  = matrix(),
    params      = list(pvalueCutoff = 1, eps = 0, pAdjustMethod = "fdr",
                       exponent = 1, minGSSize = MINGS, maxGSSize = MAXGS),
    gene2Symbol = character(),
    readable    = FALSE)
}

# ---- NES-agreement sanity gate (logged once) -------------------------------
# If a cached clusterProfiler gseaResult is on disk (gsea_msigdb_WT_heat.rds), confirm
# the plot engine (clusterProfiler) agrees with the master (fgsea) before trusting the
# rebuilt curves. We never alter the master; this only flags a material engine drift.
.nes_agreement_gate <- function() {
  cand <- file.path(DIR_OBJECTS, "gsea_msigdb_WT_heat.rds")
  if (!file.exists(cand)) {
    message("[12] NES-agreement gate: no cached gseaResult to cross-check (skipping; ",
            "@result is taken verbatim from the master either way).")
    return(invisible(NULL))
  }
  cached <- tryCatch(readRDS(cand), error = function(e) NULL)
  if (is.null(cached) || is.null(cached[["Hallmark"]])) return(invisible(NULL))
  h <- cached[["Hallmark"]]
  cdf <- data.frame(pathway_id = h@result$ID, nes_engine = h@result$NES,
                    stringsAsFactors = FALSE)
  m <- mg %>% dplyr::filter(contrast == "WT_heat", database == "Hallmark") %>%
    dplyr::transmute(pathway_id, nes_master = nes)
  j <- dplyr::inner_join(cdf, m, by = "pathway_id")
  if (nrow(j) == 0) return(invisible(NULL))
  dmax <- max(abs(j$nes_engine - j$nes_master), na.rm = TRUE)
  sgn  <- mean(sign(j$nes_engine) == sign(j$nes_master), na.rm = TRUE)
  rho  <- suppressWarnings(stats::cor(j$nes_engine, j$nes_master))
  message(sprintf("[12] NES-agreement gate (WT_heat/Hallmark, n=%d): max|dNES|=%.3f, r=%.5f, sign=%.2f",
                  nrow(j), dmax, rho, sgn))
  if (dmax > 0.5 || sgn < 0.98)
    warning("[12] NES-agreement gate FLAGGED material engine drift (max|dNES|=", round(dmax, 3),
            ", sign=", round(sgn, 3), "). Curves use the rebuilt geneList but @result NES is ",
            "the master's own — inspect before quoting running-sum leading edges.")
  invisible(j)
}
.nes_agreement_gate()

# ============================================================================
# 4. Per-contrast subtitle + caption metadata
# ============================================================================

#' Standard per-contrast subtitle line for every panel.
contrast_subtitle <- function(co) {
  sprintf("%s — NES>0 enriched in numerator/hot/WT side, NES<0 in denominator — PROVISIONAL sample labels",
          contrast_label(co))
}

CFG_KV <- sprintf(
  "thresholds.gsea_fdr=%.2g; figures.top_pathways=%d; figures.running_sum_top=%d; figures.running_sum_ylim=[%.1f,%.1f]; figures.nes_cap=%.1f; colors.diverging",
  FDR, TOPN, RSTOP, RSYLIM[1], RSYLIM[2], NESCAP)

# ============================================================================
# 5. Per-cell emitter — the four toolkit panels for one (contrast x DB)
# ============================================================================

emit_cell <- function(co, db_name, src) {
  g <- tryCatch(as_gsearesult(co, db_name),
                error = function(e) { message(sprintf("  [12] %s/%s build error: %s",
                                                      co, db_name, conditionMessage(e))); NULL })
  if (is.null(g) || nrow(g@result) == 0) {
    message(sprintf("  [12] SKIP (%s / %s): no gseaResult (0 master rows or no gene-set list)",
                    co, db_name))
    return(invisible(NULL))
  }

  n_sig <- sum(g@result$p.adjust < FDR, na.rm = TRUE)
  message(sprintf("  [12] %s / %s: %d pathways, %d sig (padj < %.2g)",
                  co, db_name, nrow(g@result), n_sig, FDR))

  ttl <- sprintf("%s — %s", db_name, contrast_label(co))
  sub <- contrast_subtitle(co)

  # Pre-create + purge the per-(contrast x DB) figure subdir so the run owns its namespace
  # (save_figure's own prefix-purge does not reach into the <DB>/ subdir).
  fig_db_dir <- file.path(contrast_path(STAGE, co, "figures", config = FIG_CFG), db_name)
  dir.create(fig_db_dir, recursive = TRUE, showWarnings = FALSE)
  stale <- list.files(fig_db_dir, pattern = "\\.(png|pdf)$", full.names = TRUE)
  if (length(stale)) file.remove(stale)

  # ---- 1. dotplot (show top-N by p.adjust; black outline = FDR-significant) ----
  p_dot <- gsea_dotplot(g, filterBy = "p.adjust", showCategory = TOPN, padj_cutoff = FDR,
                        title = ttl, neg_color = NEG, mid_color = MID, pos_color = POS,
                        nes_limits = c(-NESCAP, NESCAP)) +
    labs(subtitle = sub) +
    project_theme(config = FIG_CFG)
  if (n_sig == 0) {
    p_dot <- p_dot +
      annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4,
               label = sprintf("No significant pathways (FDR < %.2g)", FDR),
               size = (FIG_CFG$figures$label_size %||% 5) + 1, fontface = "bold",
               colour = "white", fill = NEG, label.size = 0, label.r = unit(2, "pt"),
               label.padding = unit(0.4, "lines"))
  }
  save_figure(p_dot, STAGE, file.path(db_name, "dotplot"), contrast = co, config = FIG_CFG)

  # ---- 2. up/down faceted dotplot ----
  p_fac <- gsea_dotplot_facet(g, showCategory = TOPN, padj_cutoff = FDR, title = ttl,
                              neg_color = NEG, mid_color = MID, pos_color = POS,
                              nes_limits = c(-NESCAP, NESCAP)) +
    labs(subtitle = sub) +
    project_theme(config = FIG_CFG)
  save_figure(p_fac, STAGE, file.path(db_name, "facet"), contrast = co, config = FIG_CFG)

  # ---- 3. NES barplot (significant sets only, top-N by |NES|) ----
  p_bar <- gsea_barplot(g, padj_cutoff = FDR, top_n = TOPN, title = ttl,
                        neg_color = NEG, mid_color = MID, pos_color = POS,
                        nes_limits = c(-NESCAP, NESCAP)) +
    labs(subtitle = sub) +
    project_theme(config = FIG_CFG)
  save_figure(p_bar, STAGE, file.path(db_name, "barplot"), contrast = co, config = FIG_CFG)

  # ---- 4. running-sum — REAL three-panel ES curves (the headline fix) ----
  # Top-RSTOP pathways by |NES| that have a gene-set in @geneSets (running-sum needs membership).
  # gseaplot2() draws its legend from @result$Description, so pre-format Description on a COPY of
  # the object (the dotplot/barplot/facet already ran above on the raw-id Description and
  # re-format it themselves — format_pathway_name() is not idempotent, so we must not mutate the
  # shared object before those calls). This yields a clean legend ("TNF-alpha Signaling via
  # NF-kappaB" not "HALLMARK_TNFA_...").
  have_set <- g@result$ID %in% names(g@geneSets)
  cand <- g@result[have_set, , drop = FALSE]
  if (nrow(cand) > 0) {
    top_ids <- cand$ID[order(abs(cand$NES), decreasing = TRUE)][seq_len(min(RSTOP, nrow(cand)))]
    g_rs <- g
    g_rs@result$Description <- format_pathway_name(g_rs@result$ID,
                                                   use_formatting = TRUE, strip_prefix = TRUE)
    p_rs <- tryCatch(
      gsea_running_sum_plot(g_rs, gene_set_ids = top_ids, title = ttl, max_name_length = 40),
      error = function(e) { message(sprintf("  [12] running_sum skipped (%s/%s): %s",
                                            co, db_name, conditionMessage(e))); NULL })
    if (!is.null(p_rs)) {
      # Pin the running-ES y-range to [-1,1] (config running_sum_ylim) + inside legend so
      # the per-DB curves stay directly comparable (style_series = the ported style_running_sum).
      p_rs <- style_series(p_rs, ylim = RSYLIM, config = FIG_CFG)
      save_figure(p_rs, STAGE, file.path(db_name, "running_sum"), contrast = co, config = FIG_CFG)
    }
  } else {
    message(sprintf("  [12] running_sum skipped (%s/%s): no pathway has a gene-set list", co, db_name))
  }

  invisible(TRUE)
}

# ============================================================================
# 6. Main per-contrast x per-database loop
# ============================================================================

message(sprintf("[12] Generating per-contrast panels: %d contrasts x %d databases",
                length(CONTRASTS), length(ALL_DBS)))

for (co in CONTRASTS) {
  message(sprintf("[12] == contrast: %s ==", co))
  for (db_name in MSIGDB_NAMES) emit_cell(co, db_name, "msigdb")
  for (db_name in CUSTOM_NAMES)  emit_cell(co, db_name, "custom")
}

# ============================================================================
# 7. Cross-contrast coverage metadata (db_in_master, for sections 9 & 11)
# ============================================================================

db_in_master <- intersect(ALL_DBS, unique(mg$database))

# ============================================================================
# 9. cGAS-dependence ASYMMETRY panel (_overview/)
# ============================================================================
ASYMM_IFN_PATTERN  <- "INTERFERON|INNATE_IMMUNE|INFLAMMATORY|TNF|IL6|IL2|ALLOGRAFT"
ASYMM_HIF_PATTERN  <- "HYPOXIA|GLYCOLYSIS|MTORC1|OXIDATIVE_PHOSPHORYLATION"
FOCAL_CONTRASTS    <- intersect(c("WT_heat", "KO_heat", "Interaction", "Temp_main"), CONTRASTS)

if ("Hallmark" %in% db_in_master && length(FOCAL_CONTRASTS) >= 2) {
  hall_mg <- mg %>% dplyr::filter(database == "Hallmark", contrast %in% FOCAL_CONTRASTS)

  asym_df <- hall_mg %>%
    dplyr::mutate(
      arm = dplyr::case_when(
        grepl(ASYMM_IFN_PATTERN, pathway_id, ignore.case = TRUE) ~ "IFN / Immune",
        grepl(ASYMM_HIF_PATTERN, pathway_id, ignore.case = TRUE) ~ "HIF / Glycolysis",
        TRUE ~ NA_character_)) %>%
    dplyr::filter(!is.na(arm)) %>%
    dplyr::mutate(
      label_y     = format_pathway_name(pathway_name, use_formatting = TRUE, strip_prefix = TRUE),
      significant = padj < FDR,
      contrast_lab = factor(contrast_label(contrast, short = TRUE), levels = contrast_label(FOCAL_CONTRASTS, short = TRUE)),
      arm         = factor(arm, levels = c("IFN / Immune", "HIF / Glycolysis")))

  if (nrow(asym_df) > 0) {
    top_ids_asym <- asym_df %>%
      dplyr::group_by(arm, pathway_id) %>%
      dplyr::summarise(max_abs_nes = max(abs(nes), na.rm = TRUE), .groups = "drop") %>%
      dplyr::group_by(arm) %>%
      dplyr::slice_max(order_by = max_abs_nes, n = min(TOPN, 12L), with_ties = FALSE) %>%
      dplyr::ungroup()
    asym_df <- asym_df %>% dplyr::semi_join(top_ids_asym, by = c("arm", "pathway_id"))

    p_asym <- ggplot(asym_df, aes(x = contrast_lab, y = reorder(label_y, nes), fill = nes)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(data = dplyr::filter(asym_df, significant),
                aes(label = "*"), colour = "grey15", size = 5, vjust = 0.8,
                show.legend = FALSE) +
      scale_fill_gradient2(low = NEG, mid = MID, high = POS, midpoint = 0,
                           limits = c(-NESCAP, NESCAP), oob = scales::squish, name = "NES") +
      facet_grid(arm ~ ., scales = "free_y", space = "free_y",
                 labeller = ggplot2::labeller(arm = ggplot2::label_wrap_gen(22))) +
      labs(
        title    = "cGAS-dependence asymmetry: IFN vs HIF arms (Hallmark)",
        subtitle = sprintf("* = padj < %.2g", FDR),
        x = NULL, y = NULL, caption = NULL) +
      project_theme(config = FIG_CFG) +
      ggplot2::theme(axis.text.x  = ggplot2::element_text(angle = 35, hjust = 1),
                     strip.text.y = ggplot2::element_text(angle = 0, face = "bold"),
                     panel.spacing = ggplot2::unit(0.5, "lines"))

    save_overview(
      plot      = p_asym,
      stage     = STAGE,
      name      = "gsea_asymmetry_panel",
      table     = asym_df %>%
                    dplyr::select(arm, pathway_id, pathway_name, contrast, nes, padj, set_size) %>%
                    dplyr::mutate(across(where(is.factor), as.character)),
      finding   = sprintf(
        "cGAS-dependence asymmetry in Hallmark GSEA (focal contrasts: %s). IFN/immune sets (top block) are expected to show positive NES in WT_heat and positive Interaction NES (cGAS-dependent heat-induced IFN arm). HIF/glycolysis sets (bottom block) are expected positive in WT_heat/KO_heat/Temp_main but NOT significantly enriched in Interaction — interpreted as no detectable cGAS-dependence at n=5, NOT as proven cGAS-independence (the 1-df Interaction term is underpowered).",
        paste(FOCAL_CONTRASTS, collapse = ", ")),
      script    = SCRIPT,
      fn        = "geom_tile",
      config_kv = CFG_KV,
      input     = "03_results/master/master_gsea_table.csv",
      how_to_read = paste0(
        "Two-row heatmap: top block = IFN/immune Hallmark sets; bottom block = HIF/glycolysis sets. ",
        "Color = NES (orange = NES>0 / enriched in hot/WT numerator; blue = NES<0 / enriched in cold/KO denominator). ",
        "* = padj < ", FDR, ". The Interaction contrast is the cGAS-dependence read-out: ",
        "positive/significant Interaction NES for IFN sets = cGAS-dependent; non-significant Interaction NES ",
        "for HIF sets = no detectable cGAS-dependence at n=5 (NEVER 'cGAS-independent'). Claim tier: L3. PROVISIONAL sample labels."),
      config    = FIG_CFG, wide = TRUE
    )
    message("[12] Asymmetry panel emitted: gsea_asymmetry_panel")
  } else {
    message("[12] SKIP asymmetry panel: no Hallmark IFN/HIF rows matched the patterns.")
  }
} else {
  message("[12] SKIP asymmetry panel: Hallmark not in master or <2 focal contrasts available.")
}

# ============================================================================
# 10. Lombardi-2022 vs Hallmark-Hypoxia HIF comparison panel (_overview/)
# ============================================================================
LOMBARDI_DB <- "Lombardi2022_HIF"

if (LOMBARDI_DB %in% db_in_master) {
  lomb_df <- mg %>%
    dplyr::filter(database == LOMBARDI_DB, contrast %in% CONTRASTS) %>%
    dplyr::mutate(source_label = "Lombardi-2022\n(48-gene published HIF consensus)",
                  contrast_lab = factor(contrast_label(contrast, short = TRUE), levels = contrast_label(CONTRASTS, short = TRUE)),
                  significant  = padj < FDR)

  hall_hypoxia <- mg %>%
    dplyr::filter(database == "Hallmark",
                  grepl("HYPOXIA", pathway_id, ignore.case = TRUE),
                  contrast %in% CONTRASTS) %>%
    dplyr::mutate(source_label = "Hallmark HYPOXIA\n(MSigDB canonical)",
                  contrast_lab = factor(contrast_label(contrast, short = TRUE), levels = contrast_label(CONTRASTS, short = TRUE)),
                  significant  = padj < FDR)

  combo_df <- dplyr::bind_rows(lomb_df, hall_hypoxia)

  if (nrow(combo_df) > 0) {
    p_lomb <- ggplot(combo_df,
                     aes(x = contrast_lab, y = nes, colour = source_label, shape = significant,
                         group = interaction(pathway_name, source_label))) +
      geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
      geom_line(aes(group = interaction(database, pathway_id)), linewidth = 0.7, alpha = 0.6) +
      geom_point(size = 4, stroke = 1.2) +
      scale_colour_manual(
        values = c("Lombardi-2022\n(48-gene published HIF consensus)" = POS,
                   "Hallmark HYPOXIA\n(MSigDB canonical)"             = COL_PURPLE),
        name = "HIF signature") +
      scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 1),
                         name   = sprintf("padj < %.2g", FDR),
                         labels = c("TRUE" = "Significant", "FALSE" = "Not significant")) +
      scale_y_continuous(limits = c(-NESCAP, NESCAP), oob = scales::squish) +
      labs(
        title    = "HIF orthogonal check: Lombardi-2022 vs Hallmark Hypoxia",
        subtitle = NULL,
        x = NULL,
        y = sprintf("NES (capped at ±%.1f)", NESCAP),
        caption = NULL) +
      project_theme(config = FIG_CFG) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
                     legend.position = "right")

    save_overview(
      plot      = p_lomb,
      stage     = STAGE,
      name      = "gsea_lombardi_vs_bespoke_hif",
      table     = combo_df %>%
                    dplyr::select(source_label, pathway_id, pathway_name, contrast, nes, padj, set_size) %>%
                    dplyr::mutate(across(where(is.factor), as.character)),
      finding   = sprintf(
        "Comparison of HIF-signature GSEA NES across contrasts: published Lombardi-2022 48-gene consensus (orange) vs Hallmark HYPOXIA (purple). Both are expected to track the heat response in WT_heat/KO_heat/Temp_main but NOT show significant Interaction NES — consistent with no detectable cGAS-dependence of the HIF/glycolysis arm at n=5. The Lombardi-2022 set provides an orthogonal published benchmark absent from the original hand-curated 16-gene list (~92%% heat-shock/glycolytic; the hypoxia-diagnostic core repressed, not induced)."),
      script    = SCRIPT,
      fn        = "geom_point + geom_line",
      config_kv = CFG_KV,
      input     = "03_results/master/master_gsea_table.csv",
      how_to_read = paste0(
        "Line-point traces show NES across contrasts for each HIF signature. ",
        "Filled points = padj < ", FDR, "; open points = not significant. ",
        "Positive NES = enriched in numerator (hot/WT); negative = denominator. ",
        "The Interaction contrast (if present) is the cGAS-dependence read-out: ",
        "non-significant NES for the HIF arm = no detectable cGAS-dependence at n=5 ",
        "(NEVER 'cGAS-independent'). House rule: never crown HIF1α/HIF2α as driver; ",
        "NES is an enrichment statistic (L3 tier). PROVISIONAL sample labels."),
      config    = FIG_CFG, wide = TRUE
    )
    message("[12] Lombardi-vs-bespoke panel emitted: gsea_lombardi_vs_bespoke_hif")
  } else {
    message("[12] SKIP Lombardi panel: no rows for Lombardi2022_HIF or Hallmark Hypoxia.")
  }
} else {
  message(sprintf("[12] SKIP Lombardi panel: '%s' not present in master_gsea_table.csv", LOMBARDI_DB))
}

# ============================================================================
# 11. Per-DB README captions for the per-cell toolkit panels
# ============================================================================
# One caption block per database family covering the four per-cell stems. The
# running-sum now describes a REAL three-panel enrichment curve, not a lollipop.
for (db_name in ALL_DBS) {
  if (!(db_name %in% db_in_master)) next
  src <- if (db_name %in% MSIGDB_NAMES) "msigdb" else "custom"
  db_body <- if (src == "msigdb" || db_name %in% c("Lombardi2022_HIF"))
    "IFN/ISG and immune sets enrich positively in WT-time contrasts and reverse sign in the STING/IFNAR1 interaction contrasts"
  else
    "per-contrast set enrichment for this metabolic/transport database — read each set's direction on its own, NOT as an IFN-axis signature (these DBs carry no IFN/ISG sets)"
  write_caption(
    stage    = STAGE,
    filename = sprintf("figures/by_contrast/<contrast>/%s/*.png", db_name),
    finding  = sprintf(
      "%s GSEA per contrast (dotplot/facet/barplot/running_sum) via the RNAseq-toolkit plotters on a viz-side-reconstructed gseaResult (NES/padj taken verbatim from master_gsea_table.csv; ranked vector from 02_de_results.rds; gene sets from geneset_%s_%s.rds). %s.",
      db_name, src, db_name, db_body),
    script   = SCRIPT,
    fn       = "gsea_dotplot / gsea_dotplot_facet / gsea_barplot / gsea_running_sum_plot",
    config_kv = CFG_KV,
    input    = sprintf("03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_%s_%s.rds}", src, db_name),
    how_to_read = sprintf(
      "dotplot: x = GeneRatio (leading-edge / set size), point size = -log10(FDR), fill = NES (orange %s up / blue %s down), black outline = padj < %.2g. facet: same dotplot split into Up (NES>0) vs Down (NES<0). barplot: NES bars for FDR-significant sets only, top %d by |NES|. running_sum: a REAL three-panel GSEA enrichment curve (top = running enrichment score with the leading-edge peak; middle = gene-hit ticks at each set member's rank; bottom = the ranked t-statistic) for the top %d sets by |NES|; ES y-range pinned to [%.1f,%.1f] so curves are comparable across databases. NES>0 = enriched in numerator (hot/WT side). Claim tier: L3. PROVISIONAL sample labels.",
      POS, NEG, FDR, TOPN, RSTOP, RSYLIM[1], RSYLIM[2])
  )
}

# ============================================================================
# 12. Stage README overview
# ============================================================================
write_caption(
  stage      = STAGE,
  filename   = "README overview",
  finding    = sprintf(
    "06_gsea GSEA stage: MSigDB 8 collections + %d custom databases (TransportDB, MitoPathways, MitoXplorer, Lombardi2022_HIF) across %d contrasts. Per-contrast figures in by_contrast/<c>/<DB>/ are built with the RNAseq-toolkit plotters (gsea_dotplot/facet/barplot/running_sum) on a gseaResult reconstructed viz-side from the master table + cached DE ranks + gene-set lists; the running_sum is a real enrichplot::gseaplot2 three-panel ES curve. Cross-contrast overviews in _overview/. Produced by 12_gsea_viz.R (VIZ-ONLY; GSEA computed by 05/06_gsea_*_run.R). Key biology: IFN/ISG arm is expected cGAS-dependent (positive WT_heat + Interaction NES); HIF/glycolysis arm shows no detectable cGAS-dependence at n=5 (non-significant Interaction NES). NEVER write 'cGAS-independent' — always 'no detectable cGAS-dependence at n=5'. Lombardi-2022 48-gene consensus (doi:10.1016/j.celrep.2022.111652) is the published HIF orthogonal check. Claim tier: L3. Sample labels: PROVISIONAL.",
    length(CUSTOM_NAMES), length(CONTRASTS)),
  script     = SCRIPT,
  fn         = "save_overview / save_figure",
  config_kv  = CFG_KV,
  input      = "03_results/master/master_gsea_table.csv",
  how_to_read = paste0(
    "NES > 0 = enriched in numerator side of the contrast (39 degrees C or WT genotype). ",
    "NES < 0 = enriched in denominator side (cold=37 degrees C, cGASKO). ",
    "Glyph convention: orange = up/positive, blue = down/negative (consistent with colors.diverging in config). ",
    "* or filled point = padj < ", FDR, ". ",
    "Interaction contrast: (WT_heat) - (KO_heat); positive Interaction NES = cGAS-dependent induction. ",
    "PROVISIONAL sample labels: group identity inferred from marker genes (Hspa1b thermometer + Cgas), pending collaborator sample sheet. ",
    "Claim tier: L3 (DE/enrichment statistics). Never L7 (mechanism) in a figure title."),
  config     = FIG_CFG
)

# ============================================================================
# 13. Final structural asserts
# ============================================================================
ov_fig  <- overview_path(STAGE, "figures", config = FIG_CFG)
ov_tbl  <- overview_path(STAGE, "tables",  config = FIG_CFG)

# At least one per-cell running_sum must be a real multi-panel patchwork. We re-derive
# one representative cell and assert the rebuilt object had a non-empty @geneList +
# at least one pathway with a gene-set list (so the lollipop regression cannot return).
.assert_rs <- function() {
  for (co in CONTRASTS) for (db in ALL_DBS) {
    g <- tryCatch(as_gsearesult(co, db), error = function(e) NULL)
    if (!is.null(g) && length(g@geneList) > 0 && any(g@result$ID %in% names(g@geneSets)))
      return(invisible(TRUE))
  }
  stop("12_gsea_viz: no (contrast x DB) cell yielded a gseaResult with a ranked vector AND a ",
       "gene-set list — running-sum curves cannot be real. Check 02_de_results.rds + geneset_*.rds.")
}
.assert_rs()

n_pdf <- length(list.files(file.path("03_results", STAGE, "figures"),
                           pattern = "\\.(pdf|png)$", recursive = TRUE))
n_csv <- length(list.files(file.path("03_results", STAGE, "tables"),
                           pattern = "\\.csv$", recursive = TRUE))

message(sprintf(
  "\n12_gsea_viz complete.\n  Figures: %d PDF/PNG files under 03_results/%s/figures/\n  Tables:  %d CSV files under 03_results/%s/tables/\n  Per-cell panels: dotplot/facet/barplot/running_sum (toolkit) per (contrast x DB)\n  Overview panels: gsea_asymmetry_panel, gsea_lombardi_vs_bespoke_hif",
  n_pdf, STAGE, n_csv, STAGE))
stopifnot("No output figures produced" = n_pdf >= 1)
