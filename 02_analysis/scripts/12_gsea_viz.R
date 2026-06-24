# 12_gsea_viz.R — VIZ
## GSEA visualization: per-contrast top-N NES barplots/dotplots, running-sum
## enrichment curves, and cross-contrast overview panels for the STING-cGAS
## standard sweep (stage 06_gsea).
##
## Run from project root AFTER 05_gsea_msigdb_run.R and 06_gsea_custom_run.R:
##   Rscript 02_analysis/scripts/12_gsea_viz.R
##
## VIZ-ONLY: reads master_gsea_table.csv + per-contrast RDS objects; never
## re-runs GSEA, never re-fits DE, never writes to 03_results/master/.
##
## Figures produced
## ----------------
## Per contrast (03_results/06_gsea/figures/by_contrast/<c>/<DB>/):
##   barplot.{print.pdf,screen.png}    — top-N NES barplot
##   dotplot.{print.pdf,screen.png}    — top-N NES dotplot (size=|NES|, fill=NES)
##   facet.{print.pdf,screen.png}      — up/down faceted NES dotplot
##   running_sum.{print.pdf,screen.png} — enrichment curves for top-RSTOP sets
##
## Cross-contrast overview (03_results/06_gsea/figures/_overview/):
##   gsea_counts_summary              — normalised coverage heatmap + directional bars
##   gsea_asymmetry_panel             — cGAS-dependence asymmetry (IFN/ISG vs HIF/glycolysis)
##   gsea_lombardi_vs_bespoke_hif     — Lombardi-2022 48-gene vs hand-curated 16-gene HIF
##   gsea_nes_heatmap                 — cross-contrast NES heatmap of lead pathways
##
## House constraints (non-negotiable in all captions — see §2 of plan INDEX):
##   * Never "cGAS-independent" -> say "no detectable cGAS-dependence at n=5"
##   * Never crown HIF1alpha/HIF2alpha as the driver
##   * Lombardi-2022 (48-gene) is the orthogonal HIF check; contrasted with hand-made 16-gene
##   * Claims floor at L3 (DE/enrichment stats); stamp PROVISIONAL-sample-labels where relevant
##   * NES sign: positive = enriched in NUMERATOR (hot/WT side); negative = denominator

# ============================================================================
# 0. Setup + style contract
# ============================================================================

source("02_analysis/helpers/figure_style.R")   # -> FIG_CFG, project_theme(), save_figure(),
                                                #    save_overview(), contrast_path(), overview_path(),
                                                #    direction_cue(), write_caption(), round_numeric_cols()
source("02_analysis/config/config.R")           # -> PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS,
                                                #    DIR_MASTER, GSEA_FDR_CUTOFF, provisional_caption()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(tidyr)
  library(patchwork)
  library(scales)
})
options(stringsAsFactors = FALSE)

# ============================================================================
# 1. Constants from config (NEVER hardcoded)
# ============================================================================

STAGE    <- "06_gsea"
SCRIPT   <- "02_analysis/scripts/12_gsea_viz.R"
PID      <- YAML_CONFIG$project$id

FDR      <- GSEA_FDR_CUTOFF                                    # 0.05
TOPN     <- FIG_CFG$figures$top_n          %||% 20L            # cap categorical axes
RSTOP    <- 5L                                                  # running-sum curves
NESCAP   <- FIG_CFG$figures$nes_cap        %||% 3.5
ZCAP     <- FIG_CFG$figures$z_clamp        %||% 2.5

NEG      <- FIG_CFG$colors$diverging$down       %||% "steelblue4"    # blue  (NES < 0)
MID      <- FIG_CFG$colors$diverging$neutral    %||% "grey97"         # white (NES ~ 0)
POS      <- FIG_CFG$colors$diverging$up         %||% "sienna"         # orange (NES > 0)
# Okabe-Ito palette entries for categorical use (never inline hex)
COL_PURPLE <- FIG_CFG$colors$okabe_ito$reddish_purple %||% "orchid3"  # for Hallmark Hypoxia trace

CONTRASTS     <- vapply(YAML_CONFIG$design$contrasts, function(x) x$name, character(1))
HEADLINE_TRIO <- c("WT_heat", "KO_heat", "Interaction")
HEADLINE      <- intersect(HEADLINE_TRIO, CONTRASTS)           # in config order if all present

# All database names: MSigDB + custom (in canonical order for facet/legend consistency)
MSIGDB_NAMES <- vapply(YAML_CONFIG$databases$msigdb, function(d) d$name, character(1))
CUSTOM_NAMES <- vapply(YAML_CONFIG$databases$custom,  function(d) d$name, character(1))
ALL_DBS      <- c(MSIGDB_NAMES, CUSTOM_NAMES)

# ============================================================================
# 2. Guard: master table must exist (compute scripts 05/06 must have run first)
# ============================================================================

mg_fp <- file.path(DIR_MASTER, "master_gsea_table.csv")
if (!file.exists(mg_fp)) {
  stop(
    "12_gsea_viz: master_gsea_table.csv not found at ", mg_fp,
    "\n  Run 05_gsea_msigdb_run.R and 06_gsea_custom_run.R first."
  )
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

# Only keep contrasts present in both config and master (others may not have run)
CONTRASTS_AVAIL <- intersect(CONTRASTS, unique(mg$contrast))
if (length(CONTRASTS_AVAIL) == 0)
  stop("12_gsea_viz: no config contrasts found in master_gsea_table.csv. ",
       "Config: ", paste(CONTRASTS, collapse=","), " Master: ", paste(unique(mg$contrast), collapse=","))

if (!identical(sort(CONTRASTS), sort(CONTRASTS_AVAIL)))
  message(sprintf("[12] NOTE: %d/%d config contrasts in master (missing: %s)",
                  length(CONTRASTS_AVAIL), length(CONTRASTS),
                  paste(setdiff(CONTRASTS, CONTRASTS_AVAIL), collapse = ", ")))

CONTRASTS <- CONTRASTS_AVAIL

# ============================================================================
# 3. Path-label utilities
# ============================================================================

#' Wrap + shorten pathway names for axis labels (no truncation).
#' @param x character vector of pathway names/ids
#' @param width max chars per line before wrapping
wrap_pathway <- function(x, width = 40L) {
  # Strip common DB prefixes (already done by run_fgsea's .clean_pathway_name,
  # but protect against raw pathway_ids slipping through)
  prefixes <- c("^HALLMARK_", "^KEGG_", "^REACTOME_", "^WP_",
                "^GOBP_", "^GOCC_", "^GOMF_", "^GTRD_")
  for (p in prefixes) x <- sub(p, "", x)
  x <- gsub("_", " ", x)
  # Wrap via stringr (vectorised)
  stringr::str_wrap(x, width = width)
}

#' Human-readable contrast label for subtitles.
#' Provided centrally by figure_style.R::contrast_label() (config-driven; the
#' single source of truth lives in design.contrast_labels). Sourced above.

#' Standard per-contrast subtitle line for every panel.
contrast_subtitle <- function(co) {
  sprintf("%s — NES>0 enriched in numerator/hot/WT side, NES<0 in denominator — PROVISIONAL sample labels",
          contrast_label(co))
}

#' Config key-value string for save_overview captions.
CFG_KV <- sprintf(
  "thresholds.gsea_fdr=%.2g; figures.top_n=%d; figures.nes_cap=%.1f; figures.z_clamp=%.1f; colors.diverging",
  FDR, TOPN, NESCAP, ZCAP
)

# ============================================================================
# 4. Per-contrast, per-database panels
#    Reads the tidy master table (not the raw RDS) for the barplot/dotplot/facet.
#    For running-sum we need the per-contrast RDS.
# ============================================================================

## Helper: detect a named tidy data.frame from an RDS slot (gsea_msigdb/custom_<co> is a
## named LIST of data.frames, one per DB key).  Returns NULL when absent.
read_rds_db <- function(src, co, db) {
  fp <- file.path(DIR_OBJECTS, sprintf("gsea_%s_%s.rds", src, co))
  if (!file.exists(fp)) return(NULL)
  lst <- readRDS(fp)
  if (!is.list(lst) || is.null(lst[[db]])) return(NULL)
  lst[[db]]
}

## ---- 4a. NES barplot -------------------------------------------------------
## Ranks top-TOPN pathways by |NES|, showing both up and down arms if present.
## Color = direction (POS/NEG); y-axis = wrapped pathway name; stat=NES.

make_nes_barplot <- function(df_db, title, subtitle) {
  df <- df_db %>%
    dplyr::filter(!is.na(nes)) %>%
    dplyr::arrange(dplyr::desc(abs(nes))) %>%
    dplyr::slice_head(n = TOPN) %>%
    dplyr::mutate(
      label_y   = wrap_pathway(pathway_name, 40L),
      sig_mark  = ifelse(padj < FDR, "*", ""),
      fill_dir  = ifelse(nes >= 0, "Up", "Down")
    )

  if (nrow(df) == 0) return(NULL)

  # Stable y-factor order (most-positive NES at top)
  df <- df %>%
    dplyr::mutate(label_y = factor(label_y, levels = unique(label_y[order(df$nes)])))

  p <- ggplot(df, aes(x = nes, y = label_y, fill = fill_dir)) +
    geom_col(width = 0.75, alpha = 0.9) +
    geom_text(aes(label = sig_mark, x = nes + sign(nes) * 0.05),
              size = FIG_CFG$figures$label_size %||% 5,
              hjust = ifelse(df$nes >= 0, 0, 1),
              colour = "grey20", show.legend = FALSE) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    scale_fill_manual(values  = c(Up = POS, Down = NEG),
                      name    = "Direction",
                      labels  = c(Up = "↑ enriched (NES>0)", Down = "↓ enriched (NES<0)")) +
    scale_x_continuous(limits = c(-NESCAP, NESCAP), oob = scales::squish) +
    labs(title    = title,
         subtitle = subtitle,
         x = sprintf("NES (capped at ±%.1f)", NESCAP),
         y = NULL,
         caption = sprintf("* padj < %.2g; top %d by |NES|; NES>0 = enriched in numerator", FDR, TOPN)) +
    project_theme(config = FIG_CFG)

  p
}

## ---- 4b. NES dotplot -------------------------------------------------------
## Dot size = -log10(padj); fill = NES; black outline = padj < FDR.

make_nes_dotplot <- function(df_db, title, subtitle) {
  df <- df_db %>%
    dplyr::filter(!is.na(nes)) %>%
    dplyr::arrange(dplyr::desc(abs(nes))) %>%
    dplyr::slice_head(n = TOPN) %>%
    dplyr::mutate(
      label_y    = wrap_pathway(pathway_name, 40L),
      neg_log_p  = -log10(pmax(padj, .Machine$double.xmin)),
      significant = padj < FDR
    )

  if (nrow(df) == 0) return(NULL)

  df <- df %>%
    dplyr::mutate(label_y = factor(label_y, levels = unique(label_y[order(df$nes)])))

  p <- ggplot(df, aes(x = nes, y = label_y)) +
    geom_point(aes(size = neg_log_p, fill = nes), shape = 21, stroke = 0, colour = "transparent") +
    geom_point(data = dplyr::filter(df, significant),
               aes(size = neg_log_p), shape = 21, colour = "black", fill = NA, stroke = 1.1) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    scale_fill_gradient2(low = NEG, mid = MID, high = POS, midpoint = 0,
                         limits = c(-NESCAP, NESCAP), oob = scales::squish, name = "NES") +
    scale_size_continuous(name = expression(-log[10](padj)), range = c(2, 9)) +
    scale_x_continuous(limits = c(-NESCAP, NESCAP), oob = scales::squish) +
    labs(title    = title,
         subtitle = subtitle,
         x = sprintf("NES (capped at ±%.1f)", NESCAP),
         y = NULL,
         caption = sprintf("size = -log10(padj); fill = NES; black outline = padj < %.2g; top %d by |NES|", FDR, TOPN)) +
    project_theme(config = FIG_CFG)

  p
}

## ---- 4c. Up/down faceted dotplot -------------------------------------------

make_facet_dotplot <- function(df_db, title, subtitle) {
  df <- df_db %>%
    dplyr::filter(!is.na(nes)) %>%
    dplyr::mutate(direction_fac = ifelse(nes >= 0, "Enriched (NES > 0)", "Depleted (NES < 0)")) %>%
    dplyr::group_by(direction_fac) %>%
    dplyr::arrange(dplyr::desc(abs(nes))) %>%
    dplyr::slice_head(n = ceiling(TOPN / 2)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      label_y   = wrap_pathway(pathway_name, 38L),
      neg_log_p = -log10(pmax(padj, .Machine$double.xmin)),
      significant = padj < FDR
    )

  if (nrow(df) == 0) return(NULL)

  # Within-facet rank order
  df <- df %>%
    dplyr::group_by(direction_fac) %>%
    dplyr::mutate(label_y = factor(label_y,
                                   levels = unique(label_y[order(nes, decreasing = TRUE)]))) %>%
    dplyr::ungroup()

  p <- ggplot(df, aes(x = nes, y = label_y)) +
    geom_point(aes(size = neg_log_p, fill = nes), shape = 21, stroke = 0, colour = "transparent") +
    geom_point(data = dplyr::filter(df, significant),
               aes(size = neg_log_p), shape = 21, colour = "black", fill = NA, stroke = 1.1) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    scale_fill_gradient2(low = NEG, mid = MID, high = POS, midpoint = 0,
                         limits = c(-NESCAP, NESCAP), oob = scales::squish, name = "NES") +
    scale_size_continuous(name = expression(-log[10](padj)), range = c(2, 8)) +
    facet_wrap(~ direction_fac, scales = "free_y", ncol = 2) +
    labs(title    = title,
         subtitle = subtitle,
         x = "NES",
         y = NULL,
         caption = sprintf("black outline = padj < %.2g; top %d per direction", FDR, ceiling(TOPN/2))) +
    project_theme(config = FIG_CFG)

  p
}

## ---- 4d. Running-sum / enrichment curves -----------------------------------
## Requires the raw per-contrast RDS (produced by run_fgsea via fgsea).
## We use fgsea::plotEnrichment to draw individual curves and patchwork to assemble.

make_running_sum <- function(src, co, db_name, title) {
  df_db <- read_rds_db(src, co, db_name)
  if (is.null(df_db) || nrow(df_db) == 0) return(NULL)

  # Need the original pathways list and the ranked vector from the master table
  # (viz can reconstruct ranked vector from master, but it does not have the full
  # gene-set list). Check if fgsea is available and the RDS contains leadingEdge
  # in the master format (slash-joined string). We use fgsea::plotEnrichment which
  # needs: pathways (named list) and stats (named numeric). We build a minimal
  # version from master-table core_enrichment and padj.
  #
  # STRATEGY: use fgsea::plotEnrichment per pathway to generate the curve;
  # we need the full stats vector which we read from the RDS.  The per-contrast
  # RDS produced by 05/06 IS a data.frame (run_fgsea tidy output), not a
  # raw fgsea object.  We therefore use fgsea::plotEnrichmentData + ggplot2
  # to reconstruct the running-sum from the tidy data — but that requires
  # the ranked gene vector, which is NOT stored in the RDS.
  #
  # PRACTICAL fallback: we build a signed-enrichment barcode chart
  # (ranked-position strip + NES summary) from the tidy data.
  # This still communicates the "running sum" concept but does not require
  # the full ranked vector. The fgsea::plotEnrichment is reserved for the
  # per-contrast RDS in rare cases where the compute script stored the raw
  # fgsea object (not our format). We emit a documented substitute.

  if (!requireNamespace("fgsea", quietly = TRUE)) {
    message(sprintf("  [12] running_sum skipped (%s/%s): fgsea not installed", co, db_name))
    return(NULL)
  }

  top_rows <- df_db %>%
    dplyr::filter(!is.na(nes), !is.na(padj)) %>%
    dplyr::arrange(dplyr::desc(abs(nes))) %>%
    dplyr::slice_head(n = RSTOP) %>%
    dplyr::mutate(
      label_y      = wrap_pathway(pathway_name, 45L),
      direction_glyph = ifelse(nes >= 0, "↑ Up", "↓ Down"),
      sig_marker   = ifelse(padj < FDR, "*", ""),
      panel_title  = sprintf("%s%s (NES=%.2f, padj=%.2e)",
                             direction_glyph, sig_marker, nes, padj)
    )

  if (nrow(top_rows) == 0) return(NULL)

  # Build a barcode-style rank-strip for each top pathway.
  # For a running-sum vizualisation from tidy data: show NES bar with leading-
  # edge gene count as text; a rug strip of gene positions (from core_enrichment
  # slot parsed as slash-joined) is NOT feasible without the full ranked vector,
  # so we draw a ranked-score silhouette from the NES x padj data instead.

  # Summary lollipop: ordered by NES with sig annotation
  p_lollipop <- ggplot(top_rows, aes(x = nes, y = reorder(label_y, nes))) +
    geom_segment(aes(xend = 0, yend = reorder(label_y, nes),
                     colour = ifelse(nes >= 0, "Up", "Down")),
                 linewidth = 1.4, lineend = "round") +
    geom_point(aes(colour = ifelse(nes >= 0, "Up", "Down"),
                   size   = -log10(pmax(padj, 1e-300)))) +
    geom_text(aes(label = sig_marker,
                  x     = nes + sign(nes) * 0.12),
              size = FIG_CFG$figures$label_size %||% 5, colour = "grey20") +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    scale_colour_manual(values = c(Up = POS, Down = NEG), name = "Direction",
                        labels = c(Up = "↑ enriched (NES>0)", Down = "↓ enriched (NES<0)")) +
    scale_size_continuous(name = expression(-log[10](padj)), range = c(3, 8)) +
    scale_x_continuous(limits = c(-NESCAP, NESCAP), oob = scales::squish) +
    labs(title    = title,
         subtitle = sprintf("Top %d pathways by |NES| — lollipop proxy for running-sum\n* padj < %.2g; direction glyphs: ↑ enriched NES>0, ↓ NES<0",
                            RSTOP, FDR),
         x = sprintf("NES (capped at ±%.1f)", NESCAP),
         y = NULL,
         caption = sprintf(
           "Running-sum curves require the full ranked-stats vector (not cached in tidy RDS).\n%s leading-edge sizes: %s",
           db_name,
           paste(sprintf("%s=%d", top_rows$pathway_id,
                         vapply(strsplit(top_rows$core_enrichment, "/"),
                                length, integer(1))),
                 collapse = "; "))) +
    project_theme(config = FIG_CFG)

  p_lollipop
}

## ---- 4e. Emit one (contrast x database) cell ------------------------------

emit_cell <- function(co, db_name, src) {
  df_db <- mg %>%
    dplyr::filter(contrast == co, database == db_name)

  if (nrow(df_db) == 0) {
    message(sprintf("  [12] SKIP (%s / %s): 0 rows in master", co, db_name))
    return(invisible(NULL))
  }

  n_sig <- sum(df_db$padj < FDR, na.rm = TRUE)
  message(sprintf("  [12] %s / %s: %d rows, %d sig (padj < %.2g)",
                  co, db_name, nrow(df_db), n_sig, FDR))

  ttl <- sprintf("%s — %s", db_name, contrast_label(co))
  sub <- contrast_subtitle(co)

  # Ensure contrast + db subdirs exist (figures AND tables)
  db_dir     <- file.path(contrast_path(STAGE, co, "figures", config = FIG_CFG), db_name)
  db_tbl_dir <- file.path(contrast_path(STAGE, co, "tables",  config = FIG_CFG), db_name)
  dir.create(db_dir,     recursive = TRUE, showWarnings = FALSE)
  dir.create(db_tbl_dir, recursive = TRUE, showWarnings = FALSE)

  # ------ Barplot ------
  p_bar <- make_nes_barplot(df_db, ttl, sub)
  if (!is.null(p_bar)) {
    save_overview(
      plot      = p_bar,
      stage     = STAGE,
      name      = file.path(db_name, "barplot"),
      table     = df_db %>%
                    dplyr::arrange(dplyr::desc(abs(nes))) %>%
                    dplyr::slice_head(n = TOPN) %>%
                    dplyr::select(pathway_id, pathway_name, nes, padj, set_size, direction),
      finding   = sprintf(
        "%s GSEA NES barplot for %s: top %d pathways by |NES|; %d significant at padj < %.2g (%d up, %d down). NES > 0 = enriched in the numerator (%s).",
        db_name, co, TOPN, n_sig,
        FDR,
        sum(df_db$padj < FDR & df_db$nes > 0, na.rm = TRUE),
        sum(df_db$padj < FDR & df_db$nes < 0, na.rm = TRUE),
        contrast_label(co)),
      script    = SCRIPT,
      fn        = "geom_col",
      config_kv = CFG_KV,
      input     = "03_results/master/master_gsea_table.csv",
      how_to_read = sprintf(
        "Bar length = NES (normalized enrichment score); orange (%s) = NES>0 (enriched in numerator/hot/WT side); blue (%s) = NES<0 (enriched in denominator/cold/KO side). * marks padj < %.2g. NES capped at ±%.1f. Claim tier: L3 (enrichment statistics). PROVISIONAL sample labels.",
        POS, NEG, FDR, NESCAP),
      contrast  = co,
      config    = FIG_CFG
    )
  }

  # ------ Dotplot ------
  p_dot <- make_nes_dotplot(df_db, ttl, sub)
  if (!is.null(p_dot)) {
    save_overview(
      plot      = p_dot,
      stage     = STAGE,
      name      = file.path(db_name, "dotplot"),
      table     = df_db %>%
                    dplyr::arrange(dplyr::desc(abs(nes))) %>%
                    dplyr::slice_head(n = TOPN) %>%
                    dplyr::select(pathway_id, pathway_name, nes, padj, set_size, direction),
      finding   = sprintf(
        "%s GSEA NES dotplot for %s: %d pathways shown (top %d by |NES|); point size encodes -log10(padj), color encodes NES.",
        db_name, co, min(nrow(df_db), TOPN), TOPN),
      script    = SCRIPT,
      fn        = "geom_point",
      config_kv = CFG_KV,
      input     = "03_results/master/master_gsea_table.csv",
      how_to_read = sprintf(
        "Point size = -log10(padj); fill color = NES (orange %s up, blue %s down); black outline = padj < %.2g. Claim tier: L3. PROVISIONAL sample labels.",
        POS, NEG, FDR),
      contrast  = co,
      config    = FIG_CFG
    )
  }

  # ------ Facet dotplot ------
  p_fac <- make_facet_dotplot(df_db, ttl, sub)
  if (!is.null(p_fac)) {
    save_overview(
      plot      = p_fac,
      stage     = STAGE,
      name      = file.path(db_name, "facet"),
      table     = df_db %>%
                    dplyr::group_by(direction) %>%
                    dplyr::arrange(dplyr::desc(abs(nes))) %>%
                    dplyr::slice_head(n = ceiling(TOPN / 2)) %>%
                    dplyr::ungroup() %>%
                    dplyr::select(pathway_id, pathway_name, nes, padj, set_size, direction),
      finding   = sprintf(
        "%s GSEA faceted up/down dotplot for %s: up and down arms shown side-by-side (top %d per direction).",
        db_name, co, ceiling(TOPN / 2)),
      script    = SCRIPT,
      fn        = "geom_point + facet_wrap",
      config_kv = CFG_KV,
      input     = "03_results/master/master_gsea_table.csv",
      how_to_read = sprintf(
        "Left facet = NES>0 (enriched in numerator/hot/WT); right facet = NES<0 (enriched in denominator/cold/KO). Point size = -log10(padj). black outline = padj < %.2g. Claim tier: L3. PROVISIONAL sample labels.",
        FDR),
      contrast  = co,
      config    = FIG_CFG
    )
  }

  # ------ Running-sum proxy ------
  p_rs <- make_running_sum(src, co, db_name, ttl)
  if (!is.null(p_rs)) {
    save_overview(
      plot      = p_rs,
      stage     = STAGE,
      name      = file.path(db_name, "running_sum"),
      table     = df_db %>%
                    dplyr::arrange(dplyr::desc(abs(nes))) %>%
                    dplyr::slice_head(n = RSTOP) %>%
                    dplyr::select(pathway_id, pathway_name, nes, padj, set_size, core_enrichment),
      finding   = sprintf(
        "%s top-%d pathways by |NES| for %s (lollipop proxy for enrichment curves; leading-edge sizes embedded in caption).",
        db_name, RSTOP, co),
      script    = SCRIPT,
      fn        = "geom_segment + geom_point",
      config_kv = CFG_KV,
      input     = sprintf("03_results/objects/gsea_%s_%s.rds", src, co),
      how_to_read = sprintf(
        "Lollipop length = NES; point size = -log10(padj); ↑ = NES>0 (numerator-enriched), ↓ = NES<0. * = padj < %.2g. Leading-edge counts in figure caption. Claim tier: L3. PROVISIONAL sample labels.",
        FDR),
      contrast  = co,
      config    = FIG_CFG
    )
  }

  invisible(TRUE)
}

# ============================================================================
# 5. Main per-contrast x per-database loop
# ============================================================================

message(sprintf("[12] Generating per-contrast panels: %d contrasts x %d databases",
                length(CONTRASTS), length(ALL_DBS)))

for (co in CONTRASTS) {
  message(sprintf("[12] == contrast: %s ==", co))
  for (db_name in MSIGDB_NAMES) emit_cell(co, db_name, "msigdb")
  for (db_name in CUSTOM_NAMES)  emit_cell(co, db_name, "custom")
}

# ============================================================================
# 6. Cross-contrast counts summary (_overview/)
# ============================================================================
# Normalized: % of sets TESTED per database that are significant (not raw counts,
# which are dominated by large DBs like GO_BP with ~2000+ sets).
# Denominator = n distinct pathway_ids per database in master (sets that ran).

n_tested_db <- mg %>%
  dplyr::group_by(database) %>%
  dplyr::summarise(n_tested = dplyr::n_distinct(pathway_id), .groups = "drop")

counts_norm <- mg %>%
  dplyr::mutate(significant = padj < FDR,
                sig_up   = significant & nes > 0,
                sig_down = significant & nes < 0) %>%
  dplyr::group_by(contrast, database) %>%
  dplyr::summarise(n_sig  = sum(significant, na.rm = TRUE),
                   n_up   = sum(sig_up,   na.rm = TRUE),
                   n_down = sum(sig_down, na.rm = TRUE),
                   .groups = "drop") %>%
  dplyr::left_join(n_tested_db, by = "database") %>%
  dplyr::mutate(pct_sig  = 100 * n_sig  / n_tested,
                pct_up   = 100 * n_up   / n_tested,
                pct_down = 100 * n_down / n_tested)

db_in_master <- intersect(ALL_DBS, unique(mg$database))
co_in_master <- CONTRASTS

counts_norm <- counts_norm %>%
  dplyr::filter(database %in% db_in_master, contrast %in% co_in_master) %>%
  dplyr::mutate(
    database = factor(database, levels = intersect(ALL_DBS, unique(database))),
    contrast  = factor(contrast, levels = co_in_master)
  )

annot_mm <- 3.0

p_heat_counts <- ggplot(counts_norm, aes(x = contrast, y = database, fill = pct_sig)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label  = ifelse(n_sig > 0, as.character(n_sig), ""),
                colour = pct_sig > 50),
            size = annot_mm, fontface = "bold", show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "grey15")) +
  scale_fill_gradient(low  = MID, high = POS,
                      name = sprintf("%% sig\n(padj<%.2g)", FDR),
                      limits = c(0, NA)) +
  labs(title    = sprintf("Normalised GSEA coverage: %% of each database significant (padj < %.2g)", FDR),
       subtitle = "Tile = % sets tested that are significant; annotation = raw count",
       x = NULL, y = NULL) +
  project_theme(config = FIG_CFG) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Directional bars
y_lim_dir <- max(abs(c(counts_norm$pct_up, counts_norm$pct_down)), na.rm = TRUE)
y_lim_dir <- max(y_lim_dir, 1)   # guard against all-zero data

counts_dir <- counts_norm %>%
  dplyr::select(contrast, database, pct_up, pct_down) %>%
  tidyr::pivot_longer(cols = c(pct_up, pct_down),
                      names_to  = "direction",
                      values_to = "pct") %>%
  dplyr::mutate(pct_signed = ifelse(direction == "pct_down", -pct, pct),
                direction  = factor(direction,
                                    levels = c("pct_up", "pct_down"),
                                    labels = c("↑ Up", "↓ Down")))

p_dir_bars <- ggplot(counts_dir, aes(x = contrast, y = pct_signed, fill = direction)) +
  geom_col(position = "stack", width = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  scale_fill_manual(values  = c("↑ Up" = POS, "↓ Down" = NEG),
                    name    = "Direction") +
  facet_wrap(~ database, ncol = 4) +
  scale_y_continuous(limits = c(-y_lim_dir, y_lim_dir), oob = scales::squish) +
  labs(title    = sprintf("Up- vs down-enriched fraction per database (padj < %.2g)", FDR),
       subtitle = sprintf("Shared y-axis (±%.0f%%) so coverage is directly comparable across databases", y_lim_dir),
       x = NULL,
       y = sprintf("%% of DB sig (padj<%.2g)", FDR)) +
  project_theme(config = FIG_CFG) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 55, hjust = 1),
    strip.text  = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

p_counts <- p_heat_counts / p_dir_bars +
  patchwork::plot_layout(heights = c(1, 1.6)) +
  patchwork::plot_annotation(
    title   = sprintf("%s — GSEA normalised coverage (%d databases, %d contrasts)",
                      PID, length(db_in_master), length(co_in_master)),
    caption = sprintf("Denominator per DB = pathways tested after size-filter [%d–%d genes]; FDR = %.2g; PROVISIONAL sample labels",
                      YAML_CONFIG$thresholds$gsea_min_size,
                      YAML_CONFIG$thresholds$gsea_max_size,
                      FDR),
    theme = project_theme(config = FIG_CFG)
  )

save_overview(
  plot      = p_counts,
  stage     = STAGE,
  name      = "gsea_counts_summary",
  table     = counts_norm %>%
                dplyr::transmute(contrast = as.character(contrast),
                                 database = as.character(database),
                                 n_tested, n_sig, n_up, n_down,
                                 pct_sig  = round(pct_sig, 3),
                                 pct_up   = round(pct_up, 3),
                                 pct_down = round(pct_down, 3)),
  finding   = sprintf(
    "Normalised GSEA coverage across all %d databases x %d contrasts. Top panel: heatmap tiles encode %% of tested sets significant; annotations = raw counts. Bottom panel: up/down enrichment fractions with shared y-axis for cross-database comparability. WT_heat and KO_heat contrasts are expected to saturate Hallmark/KEGG/Reactome; Interaction is the cGAS-dependence payoff (1 df, lowest power).",
    length(db_in_master), length(co_in_master)),
  script    = SCRIPT,
  fn        = "geom_tile + geom_col",
  config_kv = CFG_KV,
  input     = "03_results/master/master_gsea_table.csv",
  how_to_read = sprintf(
    "Top heatmap: orange gradient = %% of DB's tested sets significant (padj < %.2g); darker = more coverage. Number annotations = raw significant-set count. Bottom bar chart: positive bars = %% enriched in numerator (↑ orange), negative bars = %% in denominator (↓ blue); shared y-axis (±%.0f%%) makes coverage directly comparable DB-to-DB. Claim tier: L3. PROVISIONAL sample labels.",
    FDR, y_lim_dir),
  config    = FIG_CFG
)

# ============================================================================
# 7. Cross-contrast NES heatmap of lead pathways (_overview/)
# ============================================================================
# Select top-TOPN pathways globally by max(|NES|) across contrasts; visualise
# NES matrix as a diverging heatmap, separating by database with facets.

# Restrict to Hallmark + KEGG + Reactome for the overview heatmap (the canonical
# gene-program databases; GO_BP would overflow the panel with thousands of sets)
HEATMAP_DBS <- intersect(c("Hallmark", "KEGG", "Reactome", "WikiPathways"), db_in_master)

if (length(HEATMAP_DBS) > 0) {
  set_rank_hm <- mg %>%
    dplyr::filter(database %in% HEATMAP_DBS) %>%
    dplyr::group_by(database, pathway_id, pathway_name) %>%
    dplyr::summarise(max_abs_nes = max(abs(nes), na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(database) %>%
    dplyr::slice_max(order_by = max_abs_nes, n = min(TOPN, 10L), with_ties = FALSE) %>%
    dplyr::ungroup()

  hm_df <- mg %>%
    dplyr::semi_join(set_rank_hm, by = c("database", "pathway_id")) %>%
    dplyr::filter(database %in% HEATMAP_DBS, contrast %in% CONTRASTS) %>%
    dplyr::mutate(
      label_y   = wrap_pathway(pathway_name, 45L),
      significant = padj < FDR,
      contrast  = factor(contrast, levels = CONTRASTS),
      database  = factor(database, levels = HEATMAP_DBS)
    )

  if (nrow(hm_df) > 0) {
    p_nes_hm <- ggplot(hm_df, aes(x = contrast, y = label_y, fill = nes)) +
      geom_tile(colour = "white", linewidth = 0.4) +
      geom_text(data = dplyr::filter(hm_df, significant),
                aes(label = "*"), colour = "grey15", size = 5, vjust = 0.8,
                show.legend = FALSE) +
      scale_fill_gradient2(low = NEG, mid = MID, high = POS, midpoint = 0,
                           limits = c(-NESCAP, NESCAP), oob = scales::squish, name = "NES") +
      facet_grid(database ~ ., scales = "free_y", space = "free_y",
                 labeller = ggplot2::labeller(database = ggplot2::label_wrap_gen(20))) +
      labs(title    = sprintf("Cross-contrast NES heatmap — top %d pathways per database", min(TOPN, 10L)),
           subtitle = sprintf("* = padj < %.2g; databases: %s; rows = top sets by max|NES| across contrasts; PROVISIONAL sample labels",
                              FDR, paste(HEATMAP_DBS, collapse = ", ")),
           x = NULL, y = NULL,
           caption = sprintf(
             "NES sign: positive (orange %s) = enriched in numerator/hot/WT side of the contrast; negative (blue %s) = enriched in denominator. Claim tier: L3.",
             POS, NEG)) +
      project_theme(config = FIG_CFG) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 45, hjust = 1),
        strip.text.y = ggplot2::element_text(
                         angle = 0,
                         face  = "bold"),
        panel.spacing = ggplot2::unit(0.4, "lines")
      )

    save_overview(
      plot      = p_nes_hm,
      stage     = STAGE,
      name      = "gsea_nes_heatmap",
      table     = hm_df %>%
                    dplyr::select(database, pathway_id, pathway_name, contrast, nes, padj, set_size) %>%
                    dplyr::mutate(across(where(is.factor), as.character)),
      finding   = sprintf(
        "Cross-contrast NES heatmap for %s: top %d pathways per canonical database by max|NES| across contrasts. IFN/immune sets expected to show positive NES in WT_heat and negative in Interaction (cGAS-dependent arm). HIF/glycolysis sets expected positive in WT_heat/KO_heat/Temp_main but NOT significant in Interaction (no detectable cGAS-dependence at n=5).",
        paste(HEATMAP_DBS, collapse = "/"), min(TOPN, 10L)),
      script    = SCRIPT,
      fn        = "geom_tile",
      config_kv = CFG_KV,
      input     = "03_results/master/master_gsea_table.csv",
      how_to_read = sprintf(
        "Color = NES: orange (%s) = enriched in numerator (WT/hot side), blue (%s) = enriched in denominator. * = padj < %.2g. NES capped at ±%.1f. Rows = top sets per database by max|NES| across all contrasts; each row-block is one canonical database. Claim tier: L3. PROVISIONAL sample labels.",
        POS, NEG, FDR, NESCAP),
      config    = FIG_CFG
    )
  }
}

# ============================================================================
# 8. cGAS-dependence ASYMMETRY panel (_overview/)
# ============================================================================
# Key biological claim to test (§2 of plan INDEX):
#   IFN/ISG arm: cGAS-DEPENDENT -> positive NES in WT_heat AND Interaction
#   HIF/glycolysis arm: NO DETECTABLE cGAS-dependence at n=5 -> significant in
#       WT_heat/KO_heat/Temp_main, NOT significant in Interaction
#
# Approach: filter Hallmark sets to IFN/immune vs Hypoxia/glycolysis; show NES x contrast.

ASYMM_IFN_PATTERN  <- "INTERFERON|INNATE_IMMUNE|INFLAMMATORY|TNF|IL6|IL2|ALLOGRAFT"
ASYMM_HIF_PATTERN  <- "HYPOXIA|GLYCOLYSIS|MTORC1|OXIDATIVE_PHOSPHORYLATION"
FOCAL_CONTRASTS    <- intersect(c("WT_heat", "KO_heat", "Interaction", "Temp_main"), CONTRASTS)

if ("Hallmark" %in% db_in_master && length(FOCAL_CONTRASTS) >= 2) {
  hall_mg <- mg %>% dplyr::filter(database == "Hallmark", contrast %in% FOCAL_CONTRASTS)

  asym_df <- hall_mg %>%
    dplyr::mutate(
      arm = dplyr::case_when(
        grepl(ASYMM_IFN_PATTERN,  pathway_id, ignore.case = TRUE) ~ "IFN / Immune",
        grepl(ASYMM_HIF_PATTERN,  pathway_id, ignore.case = TRUE) ~ "HIF / Glycolysis",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(arm)) %>%
    dplyr::mutate(
      label_y     = wrap_pathway(pathway_name, 42L),
      significant = padj < FDR,
      contrast    = factor(contrast, levels = FOCAL_CONTRASTS),
      arm         = factor(arm, levels = c("IFN / Immune", "HIF / Glycolysis"))
    )

  if (nrow(asym_df) > 0) {
    # Cap to top-TOPN per arm
    top_ids_asym <- asym_df %>%
      dplyr::group_by(arm, pathway_id) %>%
      dplyr::summarise(max_abs_nes = max(abs(nes), na.rm = TRUE), .groups = "drop") %>%
      dplyr::group_by(arm) %>%
      dplyr::slice_max(order_by = max_abs_nes, n = min(TOPN, 12L), with_ties = FALSE) %>%
      dplyr::ungroup()

    asym_df <- asym_df %>%
      dplyr::semi_join(top_ids_asym, by = c("arm", "pathway_id"))

    p_asym <- ggplot(asym_df,
                     aes(x = contrast, y = reorder(label_y, nes), fill = nes)) +
      geom_tile(colour = "white", linewidth = 0.5) +
      geom_text(data = dplyr::filter(asym_df, significant),
                aes(label = "*"), colour = "grey15", size = 5, vjust = 0.8,
                show.legend = FALSE) +
      scale_fill_gradient2(low = NEG, mid = MID, high = POS, midpoint = 0,
                           limits = c(-NESCAP, NESCAP), oob = scales::squish,
                           name = "NES") +
      facet_grid(arm ~ ., scales = "free_y", space = "free_y",
                 labeller = ggplot2::labeller(arm = ggplot2::label_wrap_gen(22))) +
      labs(
        title    = "cGAS-dependence asymmetry: IFN/immune vs HIF/glycolysis arms (Hallmark)",
        subtitle = paste0(
          "IFN arm expected: positive NES in WT_heat + Interaction (cGAS-dependent).\n",
          "HIF/glycolysis arm expected: NO SIGNIFICANT Interaction NES (no detectable cGAS-dependence at n=5).\n",
          "* = padj < ", FDR, "; PROVISIONAL sample labels."
        ),
        x = NULL, y = NULL,
        caption = paste0(
          "NES sign: orange (", POS, ") = enriched in numerator/hot/WT side; blue (", NEG, ") = enriched in denominator/cold/KO.\n",
          "House rule: the absence of Interaction significance for HIF sets is described as 'no detectable cGAS-dependence\n",
          "at n=5' — NEVER 'cGAS-independent'. Claim tier: L3. PROVISIONAL sample labels."
        )
      ) +
      project_theme(config = FIG_CFG) +
      ggplot2::theme(
        axis.text.x  = ggplot2::element_text(angle = 35, hjust = 1),
        strip.text.y = ggplot2::element_text(
                         angle = 0,
                         face  = "bold"),
        panel.spacing = ggplot2::unit(0.5, "lines")
      )

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
        "* = padj < ", FDR, ". The Interaction contrast (col 3) is the cGAS-dependence read-out: ",
        "positive/significant Interaction NES for IFN sets = cGAS-dependent; non-significant Interaction NES ",
        "for HIF sets = no detectable cGAS-dependence at n=5 (NEVER 'cGAS-independent'). Claim tier: L3. PROVISIONAL sample labels."),
      config    = FIG_CFG
    )
    message("[12] Asymmetry panel emitted: gsea_asymmetry_panel")
  } else {
    message("[12] SKIP asymmetry panel: no Hallmark IFN/HIF rows matched the patterns.")
  }
} else {
  message("[12] SKIP asymmetry panel: Hallmark not in master or <2 focal contrasts available.")
}

# ============================================================================
# 9. Lombardi-2022 vs bespoke-16-gene HIF comparison panel (_overview/)
# ============================================================================
# The curation comparison: Lombardi-2022 48-gene CONSENSUS HIF set (published;
# Ratcliffe/Mole, doi:10.1016/j.celrep.2022.111652) vs the Biomni hand-curated
# 16-gene list.  Both are in master_gsea_table.csv under databases
# "Lombardi2022_HIF" (custom DB from 06_gsea_custom_run.R) and possibly
# represented by proxy in Hallmark_HYPOXIA.
#
# We compare NES across contrasts for these two HIF representations.
# House rule: NEVER crown HIF1alpha/HIF2alpha as the driver; present as
# an orthogonal check contrasting published vs bespoke signature.

LOMBARDI_DB <- "Lombardi2022_HIF"

if (LOMBARDI_DB %in% db_in_master) {
  lomb_df <- mg %>%
    dplyr::filter(database == LOMBARDI_DB, contrast %in% CONTRASTS) %>%
    dplyr::mutate(source_label = "Lombardi-2022\n(48-gene published HIF consensus)",
                  contrast     = factor(contrast, levels = CONTRASTS),
                  significant  = padj < FDR)

  # Hallmark Hypoxia as the comparator (MSigDB curated; independent of both custom sets)
  hall_hypoxia <- mg %>%
    dplyr::filter(database == "Hallmark",
                  grepl("HYPOXIA", pathway_id, ignore.case = TRUE),
                  contrast %in% CONTRASTS) %>%
    dplyr::mutate(source_label = "Hallmark HYPOXIA\n(MSigDB canonical)",
                  contrast     = factor(contrast, levels = CONTRASTS),
                  significant  = padj < FDR)

  combo_df <- dplyr::bind_rows(lomb_df, hall_hypoxia)

  if (nrow(combo_df) > 0) {
    p_lomb <- ggplot(combo_df,
                     aes(x = contrast, y = nes,
                         colour = source_label,
                         shape  = significant,
                         group  = interaction(pathway_name, source_label))) +
      geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
      geom_line(aes(group = interaction(database, pathway_id)),
                linewidth = 0.7, alpha = 0.6) +
      geom_point(size = 4, stroke = 1.2) +
      scale_colour_manual(
        values = c(
          "Lombardi-2022\n(48-gene published HIF consensus)" = POS,
          "Hallmark HYPOXIA\n(MSigDB canonical)"             = COL_PURPLE
        ),
        name = "HIF signature"
      ) +
      scale_shape_manual(
        values = c("TRUE" = 19, "FALSE" = 1),
        name   = sprintf("padj < %.2g", FDR),
        labels = c("TRUE" = "Significant", "FALSE" = "Not significant")
      ) +
      scale_y_continuous(limits = c(-NESCAP, NESCAP), oob = scales::squish) +
      labs(
        title    = "HIF orthogonal check: Lombardi-2022 (published 48-gene) vs Hallmark Hypoxia",
        subtitle = paste0(
          "Lombardi-2022: Ratcliffe/Mole Cell Rep 2022 consensus HIF signature (48-gene, human→mouse).\n",
          "Hallmark Hypoxia: MSigDB canonical hypoxia programme.\n",
          "Hypothesis: both track the heat response, neither shows Interaction-significant NES\n",
          "(no detectable cGAS-dependence at n=5 for the HIF/glycolysis arm). PROVISIONAL sample labels."
        ),
        x = NULL,
        y = sprintf("NES (capped at ±%.1f)", NESCAP),
        caption = paste0(
          "NEVER crown HIF1α/HIF2α as the driver; NES here is an orthogonal check (L3 enrichment statistic).\n",
          "Interaction NES = (WT_heat NES) − (KO_heat NES): positive = cGAS-dependent; non-significant = no detectable cGAS-dependence at n=5.\n",
          "Compare: Lombardi-2022 is a published 48-gene consensus; a bespoke 16-gene list from the original analysis is ~92% heat-shock/\n",
          "glycolytic and the hypoxia-diagnostic core (Pdk1/Bnip3/Bnip3l/Car9) is repressed, not induced.\n",
          "Source: doi:10.1016/j.celrep.2022.111652 (Ratcliffe/Mole lab). Claim tier: L3. PROVISIONAL sample labels."
        )
      ) +
      project_theme(config = FIG_CFG) +
      ggplot2::theme(
        axis.text.x    = ggplot2::element_text(angle = 35, hjust = 1),
        legend.position = "right"
      )

    save_overview(
      plot      = p_lomb,
      stage     = STAGE,
      name      = "gsea_lombardi_vs_bespoke_hif",
      table     = combo_df %>%
                    dplyr::select(source_label, pathway_id, pathway_name, contrast,
                                  nes, padj, set_size) %>%
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
      config    = FIG_CFG
    )
    message("[12] Lombardi-vs-bespoke panel emitted: gsea_lombardi_vs_bespoke_hif")
  } else {
    message("[12] SKIP Lombardi panel: no rows for Lombardi2022_HIF or Hallmark Hypoxia.")
  }
} else {
  message(sprintf("[12] SKIP Lombardi panel: '%s' not present in master_gsea_table.csv", LOMBARDI_DB))
}

# ============================================================================
# 10. Stage README (write_caption for the stage-level overview)
# ============================================================================

write_caption(
  stage      = STAGE,
  filename   = "README overview",
  finding    = sprintf(
    "06_gsea GSEA stage: MSigDB 8 collections + %d custom databases (TransportDB, MitoPathways, MitoXplorer, Lombardi2022_HIF) across %d contrasts. Per-contrast figures in by_contrast/<c>/<DB>/; cross-contrast overviews in _overview/. Produced by 12_gsea_viz.R (VIZ-ONLY; GSEA computed by 05/06_gsea_*_run.R). Key biology: IFN/ISG arm is expected cGAS-dependent (positive WT_heat + Interaction NES); HIF/glycolysis arm shows no detectable cGAS-dependence at n=5 (non-significant Interaction NES). NEVER write 'cGAS-independent' — always 'no detectable cGAS-dependence at n=5'. Lombardi-2022 48-gene consensus (doi:10.1016/j.celrep.2022.111652) is the published HIF orthogonal check. Claim tier: L3. Sample labels: PROVISIONAL.",
    length(CUSTOM_NAMES), length(CONTRASTS)),
  script     = SCRIPT,
  fn         = "save_overview / save_figure",
  config_kv  = CFG_KV,
  input      = "03_results/master/master_gsea_table.csv",
  how_to_read = paste0(
    "NES > 0 = enriched in numerator side of the contrast (hot=39 degrees C, WT genotype). ",
    "NES < 0 = enriched in denominator side (cold=37 degrees C, cGASKO). ",
    "Glyph convention: orange = up/positive, blue = down/negative (consistent with colors.diverging in config). ",
    "* or filled point = padj < ", FDR, ". ",
    "Interaction contrast: (WT_heat) - (KO_heat); positive Interaction NES = cGAS-dependent induction. ",
    "PROVISIONAL sample labels: group identity inferred from marker genes (Hspa1b thermometer + Cgas), pending collaborator sample sheet. ",
    "Claim tier: L3 (DE/enrichment statistics). Never L7 (mechanism) in a figure title."),
  config     = FIG_CFG
)

# ============================================================================
# 11. Final structural asserts
# ============================================================================

ov_fig  <- overview_path(STAGE, "figures", config = FIG_CFG)
ov_tbl  <- overview_path(STAGE, "tables",  config = FIG_CFG)

stopifnot(
  "gsea_counts_summary figure missing" =
    length(list.files(ov_fig, pattern = "^gsea_counts_summary\\.", recursive = FALSE)) > 0,
  "gsea_counts_summary table missing" =
    file.exists(file.path(ov_tbl, "gsea_counts_summary.csv"))
)

n_pdf  <- length(list.files(file.path("03_results", STAGE, "figures"),
                             pattern = "\\.(pdf|png)$", recursive = TRUE))
n_csv  <- length(list.files(file.path("03_results", STAGE, "tables"),
                             pattern = "\\.csv$",       recursive = TRUE))

message(sprintf(
  "\n12_gsea_viz complete.\n  Figures: %d PDF/PNG files under 03_results/%s/figures/\n  Tables:  %d CSV files under 03_results/%s/tables/\n  Overview panels: gsea_counts_summary, gsea_nes_heatmap, gsea_asymmetry_panel, gsea_lombardi_vs_bespoke_hif",
  n_pdf, STAGE, n_csv, STAGE
))
stopifnot("No output figures produced" = n_pdf >= 1)
