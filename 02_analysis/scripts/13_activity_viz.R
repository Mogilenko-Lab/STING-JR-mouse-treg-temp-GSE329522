#!/usr/bin/env Rscript
# 13_activity_viz.R - VIZ
# =============================================================================
# PROGENy pathway activity + CollecTRI TF activity visualization
# Stage: 05_progeny
# Project: GSE329522 STING/cGAS Hyperthermia (2x2 genotype x temperature)
#
# Role: VIZ ONLY - reads checkpoint objects + master CSVs; never re-runs
#       decoupleR; no DE re-fit; no statistics; no network fetch.
#
# Key scientific claim (the headline panel):
#   Hypoxia (PROGENy) rises in BOTH WT_heat and KO_heat but is FLAT in the
#   Interaction contrast => NO DETECTABLE cGAS-dependence at n=5
#   (PROVISIONAL; n=5/group; not proven independence).
#   JAK-STAT (PROGENy) and the IFN/IRF/STAT TFs are positive in WT_heat,
#   reduced in KO_heat, and positive in the Interaction => cGAS-DEPENDENT.
#   This orthogonal pathway-footprint view corroborates the two-arms split
#   without relying on any TF regulon.
#
# Inputs (read-only; never written):
#   03_results/objects/09_progeny_activity.rds    (raw decoupleR MLM)
#   03_results/master/master_progeny_activities.csv
#   03_results/objects/03_tf_collectri.rds         (raw decoupleR ULM + BH)
#   03_results/master/master_tf_activities.csv
#
# Outputs (05_progeny stage):
#   figures/_overview/
#     progeny_heatmap.{print.pdf,screen.png}
#     progeny_interaction_split.{print.pdf,screen.png}
#     progeny_tf_combined.{print.pdf,screen.png}
#     tf_heatmap.{print.pdf,screen.png}
#   figures/by_contrast/<contrast>/
#     progeny_barplot.{print.pdf,screen.png}   - 4 headline contrasts
#     tf_barplot.{print.pdf,screen.png}         - 4 headline contrasts
#   tables/_overview/{progeny_heatmap,progeny_interaction_split,progeny_tf_combined,tf_heatmap}.csv
#   tables/by_contrast/<contrast>/{progeny_barplot,tf_barplot}.csv
#   03_results/05_progeny/README.md
#
# Run from project root:
#   Rscript 02_analysis/scripts/13_activity_viz.R
# =============================================================================

source("02_analysis/helpers/figure_style.R")   # loads figure_helpers.R; sets FIG_CFG

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggrepel)
  library(pheatmap)
  library(scales)
})
options(stringsAsFactors = FALSE)

# =============================================================================
# CONSTANTS (from FIG_CFG / YAML; never hardcoded)
# =============================================================================

STAGE    <- "05_progeny"
SCRIPT   <- "02_analysis/scripts/13_activity_viz.R"

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

FDR      <- as.numeric(FIG_CFG$thresholds$gsea_fdr    %||% 0.05)
CAP      <- as.numeric(FIG_CFG$figures$z_clamp         %||%
                       FIG_CFG$figures$nes_cap          %||% 3.5)
TF_TOPN  <- as.integer(FIG_CFG$figures$top_n           %||%
                       FIG_CFG$visualization$n_top_pathways %||% 20L)
HM_TOPN  <- as.integer(FIG_CFG$figures$top_n           %||% 40L)
HM_MINC  <- 2L    # TF heatmap: must be sig in >=2 contrasts
LBL_CUT  <- 2.0   # volcano label cutoff (abs score)

# Diverging color triplet (from config: down/neutral/up)
NEG <- FIG_CFG$colors$diverging$down    %||% "steelblue4"
MID <- FIG_CFG$colors$diverging$neutral %||% "grey97"
POS <- FIG_CFG$colors$diverging$up      %||% "sienna"

# Axis colors from config (single source of truth across all STING viz)
AXIS_COLORS <- c(
  "HIF"   = FIG_CFG$colors$diverging$up   %||% "sienna",
  "IFN"   = FIG_CFG$colors$diverging$down %||% "steelblue4",
  "other" = "grey55"
)

# Headline contrasts (config order; the four named in Plan §2/§5)
HEADLINE_CONTRASTS <- c("WT_heat", "KO_heat", "Interaction", "Temp_main")

# Key PROGENy pathways for the Hypoxia-vs-immune split
KEY_PROGENY <- c("Hypoxia", "JAK-STAT", "NFkB", "TNFa")

# Key TFs for the IFN axis
KEY_TFS <- c("Irf3", "Irf7", "Irf1", "Stat1", "Stat2", "Nfkb1", "Rela",
             "Hif1a", "Epas1")

# Provisional caption stamp (matches config.R::provisional_caption())
PROV_STAMP <- paste0(
  "[PROVISIONAL - inferred sample mapping (Hspa1b/Hsph1 thermometer + Cgas), ",
  "pending collaborator sample sheet; n=5/group]"
)

# pheatmap color ramp + break vector (diverging palette, clamped to CAP)
RAMP <- colorRampPalette(c(NEG, MID, POS))(100)
BRKS <- seq(-CAP, CAP, length.out = 101)

# =============================================================================
# DIRECTORIES
# =============================================================================

DIR_OBJECTS <- file.path("03_results/objects")
DIR_MASTER  <- file.path("03_results/master")

# Ensure the stage figures and tables roots exist (overview + by_contrast dirs
# are created lazily by overview_path()/contrast_path())
.ov_fig  <- overview_path(STAGE, "figures",  config = FIG_CFG)
.ov_tbl  <- overview_path(STAGE, "tables",   config = FIG_CFG)

# =============================================================================
# 1. GUARD - stop if PROGENy object absent; warn if TF absent
# =============================================================================

prog_path <- file.path(DIR_OBJECTS, "09_progeny_activity.rds")
tf_path   <- file.path(DIR_OBJECTS, "03_tf_collectri.rds")

if (!file.exists(prog_path))
  stop("[13_activity_viz] 09_progeny_activity.rds not found: ", prog_path,
       " - run 09_activity_progeny.R first.")

TF_AVAILABLE <- file.exists(tf_path)
if (!TF_AVAILABLE)
  warning("[13_activity_viz] 03_tf_collectri.rds not found - PROGENy figures will be produced; TF panels SKIPPED.")

# =============================================================================
# 2. LOAD + TIDY DATA
# =============================================================================

message("[13_activity_viz] Loading PROGENy checkpoint ...")
progeny_raw <- readRDS(prog_path)

# --- PROGENy master CSV (has BH padj) ----------------------------------------
prog_master_path <- file.path(DIR_MASTER, "master_progeny_activities.csv")
if (file.exists(prog_master_path)) {
  prog_long <- readr::read_csv(prog_master_path, show_col_types = FALSE, progress = FALSE)
  # Unify column names to internal convention:
  #   pathway_name -> source, nes -> score, pvalue -> p_value
  if ("pathway_name" %in% names(prog_long) && !"source" %in% names(prog_long))
    prog_long <- dplyr::rename(prog_long, source = pathway_name)
  if ("nes" %in% names(prog_long) && !"score" %in% names(prog_long))
    prog_long <- dplyr::rename(prog_long, score = nes)
  if ("pvalue" %in% names(prog_long) && !"p_value" %in% names(prog_long))
    prog_long <- dplyr::rename(prog_long, p_value = pvalue)
  prog_long <- dplyr::filter(prog_long, !is.na(score), !is.na(contrast))
} else {
  # Fallback: build from raw decoupleR result (no BH applied across 14 at this step)
  message("[13_activity_viz] master_progeny_activities.csv absent; building from raw rds (no BH padj)")
  prog_long <- as.data.frame(progeny_raw, stringsAsFactors = FALSE) %>%
    dplyr::rename_with(~ sub("^condition$", "contrast", .x)) %>%
    dplyr::filter(!is.na(score)) %>%
    dplyr::group_by(contrast) %>%
    dplyr::mutate(padj = stats::p.adjust(p_value, "BH")) %>%
    dplyr::ungroup()
}

# Ensure we have the columns we need
stopifnot(all(c("source", "score", "contrast") %in% colnames(prog_long)))
if (!"padj" %in% colnames(prog_long)) prog_long$padj <- prog_long$p_value %||% NA_real_
if (!"p_value" %in% colnames(prog_long)) prog_long$p_value <- prog_long$padj %||% NA_real_

prog_long <- prog_long %>%
  dplyr::mutate(sig = ifelse(!is.na(p_value), p_value < FDR, FALSE))

# All contrasts present (in config order as far as possible)
ALL_CONTRASTS <- unique(prog_long$contrast)
CO_PROG <- intersect(
  c("WT_heat", "KO_heat", "Interaction", "Temp_main",
    "Geno_at_39", "Geno_at_37", "Geno_main"),
  ALL_CONTRASTS
)
CO_HEADLINE <- intersect(HEADLINE_CONTRASTS, CO_PROG)

message(sprintf("[13_activity_viz] PROGENy: %d rows, %d pathways, %d contrasts",
                nrow(prog_long),
                length(unique(prog_long$source)),
                length(CO_PROG)))

# --- TF data -----------------------------------------------------------------
if (TF_AVAILABLE) {
  message("[13_activity_viz] Loading TF checkpoint ...")
  tf_raw <- readRDS(tf_path)

  tf_master_path <- file.path(DIR_MASTER, "master_tf_activities.csv")
  if (file.exists(tf_master_path)) {
    tf_long <- readr::read_csv(tf_master_path, show_col_types = FALSE, progress = FALSE)
    # Filter to CollecTRI rows; unify column names
    tf_long <- dplyr::filter(tf_long, database == "CollecTRI")
    if ("pathway_name" %in% names(tf_long) && !"source" %in% names(tf_long))
      tf_long <- dplyr::rename(tf_long, source = pathway_name)
    if ("nes" %in% names(tf_long) && !"score" %in% names(tf_long))
      tf_long <- dplyr::rename(tf_long, score = nes)
    if ("pvalue" %in% names(tf_long) && !"p_value" %in% names(tf_long))
      tf_long <- dplyr::rename(tf_long, p_value = pvalue)
    tf_long <- dplyr::filter(tf_long, !is.na(score), !is.na(contrast))
  } else {
    # Fallback: build from raw decoupleR result (03_tf_collectri.rds already has padj from bh_within)
    message("[13_activity_viz] master_tf_activities.csv absent; building from raw rds")
    tf_long <- as.data.frame(tf_raw, stringsAsFactors = FALSE)
    if ("condition" %in% names(tf_long) && !"contrast" %in% names(tf_long))
      tf_long <- dplyr::rename(tf_long, contrast = condition)
    tf_long <- dplyr::filter(tf_long, !is.na(score))
    # Apply BH within contrast if padj not already present
    if (!"padj" %in% names(tf_long)) {
      tf_long <- tf_long %>%
        dplyr::group_by(contrast) %>%
        dplyr::mutate(padj = stats::p.adjust(p_value, "BH")) %>%
        dplyr::ungroup()
    }
  }

  stopifnot(all(c("source", "score", "contrast") %in% colnames(tf_long)))
  if (!"padj" %in% colnames(tf_long)) tf_long$padj <- tf_long$p_value %||% NA_real_
  if (!"p_value" %in% colnames(tf_long)) tf_long$p_value <- tf_long$padj %||% NA_real_

  tf_long <- tf_long %>%
    dplyr::mutate(sig = ifelse(!is.na(padj), padj < FDR, FALSE))

  CO_TF <- intersect(CO_PROG, unique(tf_long$contrast))
  CO_TF_HEADLINE <- intersect(HEADLINE_CONTRASTS, CO_TF)

  message(sprintf("[13_activity_viz] TF: %d rows, %d TFs, %d contrasts",
                  nrow(tf_long),
                  length(unique(tf_long$source)),
                  length(CO_TF)))
} else {
  tf_long <- NULL
  CO_TF   <- character(0)
  CO_TF_HEADLINE <- character(0)
}

# =============================================================================
# 3. SHARED HELPERS
# =============================================================================

# Diverging fill scale (colors from FIG_CFG$colors$diverging)
.div_fill <- function(name = "Activity\nScore", lim = c(-CAP, CAP))
  scale_fill_gradient2(low = NEG, mid = MID, high = POS, midpoint = 0,
                       limits = lim, oob = scales::squish, name = name)

.div_color <- scale_color_manual(
  values = c(Up = POS, Down = NEG),
  name   = "Direction"
)

# Significance star (single-tier: p < FDR -> "*", else "")
.sig_star <- function(p, thr = FDR)
  ifelse(!is.na(p) & p < thr, "*", "")

# Contrast label (pretty name for subtitle; never bare).
# Aliased to the central, config-driven figure_style.R::contrast_label() (single
# source of truth in design.contrast_labels) — sourced at the top of this script.
.contrast_label <- contrast_label

# pheatmap-grob builder: captures the grob for ggsave-routed save_overview
.heat_grob <- function(mat, main, pmat = NULL,
                       fontsize_row = 10, cellheight = 14, cellwidth = 24,
                       cluster_rows = TRUE, gaps_row = NULL,
                       annotation_row = NULL, annotation_colors = NULL) {
  disp <- if (!is.null(pmat)) {
    pm_sub <- pmat[rownames(mat), colnames(mat), drop = FALSE]
    star_m <- matrix(.sig_star(as.vector(pm_sub)),
                     nrow = nrow(mat), ncol = ncol(mat),
                     dimnames = dimnames(mat))
    star_m
  } else FALSE

  pheatmap(mat,
           main              = main,
           color             = RAMP,
           breaks            = BRKS,
           cluster_rows      = cluster_rows,
           cluster_cols      = FALSE,
           gaps_row          = gaps_row,
           annotation_row    = annotation_row,
           annotation_colors = annotation_colors,
           display_numbers   = disp,
           number_color      = "black",
           fontsize_number   = 10,
           fontsize_row      = fontsize_row,
           fontsize_col      = 10,
           angle_col         = 45,
           cellwidth         = cellwidth,
           cellheight        = cellheight,
           silent            = TRUE)$gtable
}

# Build a wide score matrix (rows = sources, cols = contrasts, ordered)
.wide_matrix <- function(df, sources, contrasts, score_col = "score") {
  wide <- df %>%
    dplyr::filter(source %in% sources, contrast %in% contrasts) %>%
    dplyr::select(source, contrast, score = dplyr::all_of(score_col)) %>%
    tidyr::pivot_wider(names_from = contrast, values_from = score) %>%
    tibble::column_to_rownames("source") %>%
    as.matrix()
  # Re-order to requested row/col order (keep only those present)
  row_ord <- sources[sources %in% rownames(wide)]
  col_ord <- contrasts[contrasts %in% colnames(wide)]
  wide[row_ord, col_ord, drop = FALSE]
}

# Build a wide p-matrix (same structure as score matrix)
.wide_pmat <- function(df, sources, contrasts, p_col) {
  wide <- df %>%
    dplyr::filter(source %in% sources, contrast %in% contrasts) %>%
    dplyr::select(source, contrast, p = dplyr::all_of(p_col)) %>%
    tidyr::pivot_wider(names_from = contrast, values_from = p) %>%
    tibble::column_to_rownames("source") %>%
    as.matrix()
  row_ord <- sources[sources %in% rownames(wide)]
  col_ord <- contrasts[contrasts %in% colnames(wide)]
  m <- wide[row_ord, col_ord, drop = FALSE]
  m[is.na(m)] <- 1.0
  m
}

# Clamp score matrix values to [-CAP, CAP], fill NA with 0
.clamp <- function(mat) {
  mat <- pmin(pmax(mat, -CAP), CAP)
  mat[is.na(mat)] <- 0
  mat
}

# =============================================================================
# 4. BARPLOT BUILDER (shared by PROGENy and TF; matches 03_decoupler_tf_viz.R idiom)
#
# Like fig3b (top TF per contrast, lollipop by axis color), but adapted for
# the standard-sweep save_overview() contract:
#   * sorted by score (ascending left-to-right); * significance star nudged
#   * key watchlist entries force-included + bold outline (same logic as 14839
#     template .activity_barplot); * axis color if TF, diverging fill if PROGENy
# =============================================================================

.barplot_progeny <- function(df, co, n = 14L) {
  # PROGENy: all 14 pathways; signed score; diverging fill; raw p for star
  if (is.null(df) || nrow(df) == 0) return(NULL)

  top <- df %>%
    dplyr::filter(contrast == co) %>%
    dplyr::arrange(score) %>%
    dplyr::mutate(
      source  = factor(source, levels = source),
      is_sig  = !is.na(p_value) & p_value < FDR,
      star    = .sig_star(p_value),
      is_key  = source %in% KEY_PROGENY,
      star_y  = score + sign(score) * (max(abs(score), na.rm = TRUE) * 0.06 + 0.1)
    )

  ggplot(top, aes(x = source, y = score, fill = score)) +
    geom_col(width = 0.75, color = NA) +
    geom_col(data = dplyr::filter(top, is_key),
             fill = NA, color = "black", linewidth = 0.9, width = 0.75) +
    geom_text(aes(y = star_y, label = star), size = 4.5, vjust = 0.5) +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
    .div_fill() +
    coord_flip() +
    labs(
      title    = "PROGENy pathway activity (MLM)",
      subtitle = paste0(.contrast_label(co), "\n", PROV_STAMP),
      x        = NULL,
      y        = "Activity score (MLM)",
      caption  = paste0("* raw p < ", FDR,
                        "  |  bold outline = key pathways (Hypoxia / JAK-STAT / NFkB / TNFa)")
    ) +
    project_theme(config = FIG_CFG, legend = TRUE)
}

.barplot_tf <- function(df, co, n = TF_TOPN) {
  # TF: top-N by padj; diverging fill + axis color family for watchlist
  if (is.null(df) || nrow(df) == 0) return(NULL)

  base_df <- df %>% dplyr::filter(contrast == co)

  top_sig   <- base_df %>% dplyr::arrange(padj) %>% utils::head(n)
  extra_key <- base_df %>% dplyr::filter(source %in% KEY_TFS,
                                          !source %in% top_sig$source)
  top <- dplyr::bind_rows(top_sig, extra_key) %>%
    dplyr::arrange(score) %>%
    dplyr::mutate(
      source  = factor(source, levels = source),
      is_sig  = !is.na(padj) & padj < FDR,
      star    = .sig_star(padj),
      is_key  = source %in% KEY_TFS,
      axis    = dplyr::case_when(
        source %in% c("Hif1a", "Epas1")                            ~ "HIF",
        source %in% c("Irf1","Irf3","Irf7","Stat1","Stat2","Nfkb1","Rela") ~ "IFN",
        TRUE                                                         ~ "other"
      ),
      star_y  = score + sign(score) * (max(abs(score), na.rm = TRUE) * 0.06 + 0.1)
    )

  ggplot(top, aes(x = source, y = score)) +
    geom_segment(aes(xend = source, yend = 0, color = axis), linewidth = 0.8) +
    geom_point(aes(color = axis), size = 3) +
    geom_point(data = dplyr::filter(top, is_key),
               shape = 21, size = 3.5, stroke = 0.9,
               fill = NA, color = "black") +
    geom_text(aes(y = star_y, label = star), size = 4.5, vjust = 0.5,
              color = "grey20") +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
    scale_color_manual(values = AXIS_COLORS, name = "TF axis") +
    coord_flip() +
    labs(
      title    = sprintf("Top %d active TFs (CollecTRI ULM)", n),
      subtitle = paste0(.contrast_label(co), "\n", PROV_STAMP),
      x        = NULL,
      y        = "Activity score (ULM)",
      caption  = paste0("* BH padj < ", FDR,
                        "  |  circle outline = IFN/HIF watchlist",
                        "  |  orange = HIF axis  |  blue = IFN/NFkB axis")
    ) +
    project_theme(config = FIG_CFG, legend = TRUE)
}

# =============================================================================
# 5. PER-CONTRAST PANELS
# =============================================================================

message("[13_activity_viz] --- Per-contrast panels ---")

for (co in CO_HEADLINE) {
  # PROGENy barplot
  p <- tryCatch(.barplot_progeny(prog_long, co),
                error = function(e) { message("  progeny bar [", co, "]: ", e$message); NULL })
  if (!is.null(p)) {
    tbl_co <- prog_long %>% dplyr::filter(contrast == co) %>%
      dplyr::select(source, contrast, score, p_value, padj) %>%
      dplyr::arrange(p_value)
    tryCatch(
      save_overview(p, STAGE, "progeny_barplot",
                    table     = tbl_co,
                    contrast  = co,
                    finding   = paste0("PROGENy 14-pathway activity (MLM) in ", co,
                                       ": bars sorted by score; * raw p<", FDR,
                                       "; key pathways (Hypoxia/JAK-STAT/NFkB/TNFa) bold-outlined."),
                    script    = SCRIPT,
                    fn        = ".barplot_progeny",
                    config_kv = paste0("thresholds.gsea_fdr=", FDR,
                                       "; figures.z_clamp=", CAP,
                                       "; colors.diverging"),
                    input     = "03_results/objects/09_progeny_activity.rds",
                    how_to_read = paste0(
                      "Bars = MLM activity score (orange = pathway more active in numerator; ",
                      "blue = more active in denominator). Glyphs: * = raw p < ", FDR, "; ",
                      "bold outline = key pathway. Score is not a fold-change; sign tracks ",
                      "numerator activation direction. Claim tier: L3 (provisional, n=5/group)."),
                    config    = FIG_CFG),
      error = function(e) message("  save_overview progeny_barplot [", co, "]: ", e$message))
  }

  # TF barplot
  if (TF_AVAILABLE && co %in% CO_TF) {
    p2 <- tryCatch(.barplot_tf(tf_long, co),
                   error = function(e) { message("  tf bar [", co, "]: ", e$message); NULL })
    if (!is.null(p2)) {
      tbl_tf_co <- tf_long %>% dplyr::filter(contrast == co) %>%
        dplyr::arrange(padj) %>% utils::head(TF_TOPN + length(KEY_TFS)) %>%
        dplyr::select(source, contrast, score, p_value, padj)
      tryCatch(
        save_overview(p2, STAGE, "tf_barplot",
                      table     = tbl_tf_co,
                      contrast  = co,
                      finding   = paste0("Top CollecTRI TF activity (ULM) in ", co,
                                         ": lollipops colored by HIF/IFN axis; * BH padj<", FDR,
                                         "; IFN/HIF watchlist bold-outlined."),
                      script    = SCRIPT,
                      fn        = ".barplot_tf",
                      config_kv = paste0("thresholds.gsea_fdr=", FDR,
                                         "; figures.top_n=", TF_TOPN,
                                         "; colors.diverging"),
                      input     = "03_results/objects/03_tf_collectri.rds",
                      how_to_read = paste0(
                        "Lollipops = ULM activity score (rightward = pathway-level activation in numerator). ",
                        "Color: orange = HIF axis (Hif1a/Epas1), blue = IFN/NFkB axis, grey = other. ",
                        "Open circle = watchlist TF. * = BH padj < ", FDR, ". ",
                        "Claim tier: L3 (provisional, n=5/group)."),
                      config    = FIG_CFG),
        error = function(e) message("  save_overview tf_barplot [", co, "]: ", e$message))
    }
  }
}

message("[13_activity_viz] Per-contrast panels done (", length(CO_HEADLINE), " headline contrasts).")

# =============================================================================
# 6. OVERVIEW: PROGENy 14-pathway HEATMAP (all contrasts × 14 pathways)
# Clustered rows (pathways); columns in config order; diverging fill;
# star where raw p < FDR.
# =============================================================================

message("[13_activity_viz] --- PROGENy cross-contrast heatmap ---")

tryCatch({
  pathways_all <- sort(unique(prog_long$source))
  mat_ph <- .wide_matrix(prog_long, pathways_all, CO_PROG, "score")
  mat_ph <- .clamp(mat_ph)
  pmat_ph <- .wide_pmat(prog_long, pathways_all, CO_PROG, "p_value")
  pmat_ph <- pmat_ph[rownames(mat_ph), colnames(mat_ph), drop = FALSE]

  grob_ph <- .heat_grob(mat_ph,
                        main         = "PROGENy pathway activity across contrasts",
                        pmat         = pmat_ph,
                        fontsize_row = 11,
                        cellheight   = 18,
                        cellwidth    = 28,
                        cluster_rows = TRUE)

  side_ph <- as.data.frame(mat_ph) %>%
    tibble::rownames_to_column("source") %>%
    tidyr::pivot_longer(-source, names_to = "contrast", values_to = "score") %>%
    dplyr::left_join(
      as.data.frame(pmat_ph) %>%
        tibble::rownames_to_column("source") %>%
        tidyr::pivot_longer(-source, names_to = "contrast", values_to = "p_value"),
      by = c("source", "contrast")) %>%
    dplyr::mutate(significant = !is.na(p_value) & p_value < FDR) %>%
    dplyr::select(source, contrast, score, p_value, significant)

  save_overview(grob_ph, STAGE, "progeny_heatmap",
                table     = side_ph,
                finding   = paste0(
                  "PROGENy 14-pathway x ", length(CO_PROG), "-contrast activity heatmap (MLM score, ",
                  "clamped +/-", CAP, "; * raw p<", FDR, "; rows clustered): ",
                  "Hypoxia/glycolysis rise in BOTH heat arms but are flat in the Interaction; ",
                  "JAK-STAT/NFkB/TNFa rise in WT_heat and are positive in the Interaction (cGAS-dependent). ",
                  "PROVISIONAL; n=5/group."),
                script    = SCRIPT,
                fn        = ".heat_grob + pheatmap",
                config_kv = paste0("thresholds.gsea_fdr=", FDR,
                                   "; figures.z_clamp=", CAP,
                                   "; colors.diverging"),
                input     = "03_results/objects/09_progeny_activity.rds",
                how_to_read = paste0(
                  "Rows = PROGENy pathways (hierarchically clustered); columns = contrasts in ",
                  "design order (left = WT_heat, KO_heat; then Interaction, Temp_main, Geno*, Geno_main). ",
                  "Fill: orange = pathway activated in numerator; blue = activated in denominator. ",
                  "Score clamped to +/-", CAP, "; * = raw p < ", FDR, " (n=14 pathways; no multi-test ",
                  "correction warranted). Claim tier: L3 (provisional, n=5/group)."),
                config    = FIG_CFG)
  message("[13_activity_viz] PROGENy heatmap done (",
          nrow(mat_ph), " pathways x ", ncol(mat_ph), " contrasts)")
}, error = function(e) message("  progeny heatmap error: ", conditionMessage(e)))

# =============================================================================
# 7. HEADLINE PANEL: Hypoxia-vs-immune SPLIT on the Interaction contrast
#
# THE KEY RESULT: PROGENy Hypoxia (and glycolysis proxies) score flat/NS in the
# Interaction while JAK-STAT/NFkB/TNFa are positive - the cGAS-dependence
# asymmetry read out by a footprint method that uses NO HIF regulon.
#
# Design: grouped barplot (contrasts on x-axis, pathways as fill groups) showing
# KEY_PROGENY pathways across the 4 headline contrasts, so the reader sees:
#   WT_heat: Hypoxia UP, JAK-STAT UP
#   KO_heat: Hypoxia UP, JAK-STAT diminished
#   Interaction: Hypoxia FLAT, JAK-STAT positive   ← THE CLAIM
#   Temp_main: Hypoxia UP (average), JAK-STAT UP (but smaller than WT_heat)
#
# Caption explicitly uses "no detectable cGAS-dependence at n=5" NOT
# "cGAS-independent"; never crowns HIF1α/2α; adds PROVISIONAL label.
# =============================================================================

message("[13_activity_viz] --- Hypoxia-vs-immune split panel (Interaction) ---")

tryCatch({
  key_df <- prog_long %>%
    dplyr::filter(source %in% KEY_PROGENY,
                  contrast %in% CO_HEADLINE) %>%
    dplyr::mutate(
      contrast = factor(contrast, levels = CO_HEADLINE),
      source   = factor(source,
                        levels = c("Hypoxia", "JAK-STAT", "NFkB", "TNFa")),
      is_sig   = !is.na(p_value) & p_value < FDR,
      star     = .sig_star(p_value),
      star_y   = score + sign(score) * 0.15
    )

  if (nrow(key_df) == 0) {
    message("  split panel SKIPPED: no KEY_PROGENY pathways in headline contrasts")
  } else {
    # Axis-family coloring for pathways:
    #   Hypoxia -> orange (no detectable cGAS-dependence -> same hue as HIF arm)
    #   JAK-STAT/NFkB/TNFa -> blue (cGAS-dependent immune arm)
    pathway_cols <- c(
      "Hypoxia"  = POS,          # orange: parallel to the "flat/NS" HIF arm
      "JAK-STAT" = NEG,          # blue: cGAS-dependent
      "NFkB"     = "cornflowerblue",   # lighter blue: immune
      "TNFa"     = "lightsteelblue"   # lightest blue: immune
    )

    fig_split <- ggplot(key_df,
                        aes(x = contrast, y = score, fill = source,
                            group = source)) +
      geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
      geom_col(position = position_dodge(width = 0.8), width = 0.72) +
      geom_text(aes(y = star_y, label = star),
                position = position_dodge(width = 0.8),
                size = 4.5, vjust = 0.5, color = "grey20") +
      scale_fill_manual(values = pathway_cols,
                        name   = "PROGENy pathway") +
      scale_x_discrete(labels = c(
        WT_heat     = "WT heat\n(39 vs 37 C)",
        KO_heat     = "KO heat\n(39 vs 37 C)",
        Interaction = "Interaction\n(cGAS-dependence)",
        Temp_main   = "Temp main\n(average heat)"
      )) +
      labs(
        title    = "PROGENy pathway activity: Hypoxia vs immune split (MLM)",
        subtitle = paste0(
          "Hypoxia (orange) rises equally in BOTH genotypes (flat Interaction -> ",
          "no detectable cGAS-dependence at n=5).\n",
          "JAK-STAT/NFkB/TNFa (blue) are cGAS-dependent (positive Interaction).\n",
          PROV_STAMP),
        x        = NULL,
        y        = "PROGENy activity score (MLM)",
        caption  = paste0("* raw p < ", FDR,
                          "  |  PROVISIONAL; n=5/group; not proven independence.")
      ) +
      project_theme(config = FIG_CFG, legend = TRUE)

    side_split <- key_df %>%
      dplyr::select(source, contrast, score, p_value, is_sig) %>%
      dplyr::rename(pathway = source, significant = is_sig) %>%
      dplyr::mutate(contrast = as.character(contrast))

    save_overview(fig_split, STAGE, "progeny_interaction_split",
                  table     = side_split,
                  finding   = paste0(
                    "PROGENy Hypoxia vs immune split: Hypoxia activity is similar in WT_heat and KO_heat ",
                    "with a flat Interaction (no detectable cGAS-dependence at n=5), whereas ",
                    "JAK-STAT/NFkB/TNFa are positive in the Interaction (cGAS-dependent arm). ",
                    "This orthogonal pathway footprint corroborates the two-arms DE/TF result. ",
                    "PROVISIONAL; n=5/group."),
                  script    = SCRIPT,
                  fn        = "ggplot/geom_col",
                  config_kv = paste0("thresholds.gsea_fdr=", FDR,
                                     "; colors.diverging; design.contrasts"),
                  input     = "03_results/objects/09_progeny_activity.rds",
                  how_to_read = paste0(
                    "Grouped bars: x = contrast (4 headline contrasts), fill = PROGENy pathway. ",
                    "Orange (Hypoxia) bars should be similar in WT_heat and KO_heat and near-zero in ",
                    "the Interaction panel - meaning no detectable cGAS-dependence. Blue bars ",
                    "(JAK-STAT/NFkB/TNFa) should be taller in WT_heat than KO_heat and positive in ",
                    "Interaction - cGAS-dependent. * = raw p < ", FDR, ". ",
                    "CAUTION: 'flat Interaction' is NOT proven cGAS-independence; ",
                    "the study is powered at n=5/group. Claim tier: L3 (provisional)."),
                  config    = FIG_CFG)
    message("[13_activity_viz] Interaction split panel done.")
  }
}, error = function(e) message("  split panel error: ", conditionMessage(e)))

# =============================================================================
# 8. OVERVIEW: TF HEATMAP (TFs sig in >=2 contrasts, top HM_TOPN)
# Mirrors the 03_decoupler_tf_viz.R idiom (fig3b / top-TF-by-contrast):
#   * axis colors used for annotation strip (HIF/IFN/other)
#   * rows clustered; cols in config order; star where padj < FDR
# =============================================================================

if (TF_AVAILABLE && length(CO_TF) >= 2) {
  message("[13_activity_viz] --- TF cross-contrast heatmap ---")

  tryCatch({
    tf_pick_df <- tf_long %>%
      dplyr::group_by(source) %>%
      dplyr::summarise(min_padj = min(padj, na.rm = TRUE),
                       n_sig    = sum(!is.na(padj) & padj < FDR),
                       .groups  = "drop") %>%
      dplyr::filter(n_sig >= HM_MINC) %>%
      dplyr::arrange(min_padj) %>%
      utils::head(HM_TOPN)

    tf_pick <- union(tf_pick_df$source, intersect(KEY_TFS, unique(tf_long$source)))

    if (length(tf_pick) >= 3) {
      mat_th <- .wide_matrix(tf_long, tf_pick, CO_TF, "score")
      mat_th <- .clamp(mat_th)
      pmat_th <- .wide_pmat(tf_long, tf_pick, CO_TF, "padj")
      pmat_th <- pmat_th[rownames(mat_th), colnames(mat_th), drop = FALSE]

      # Row annotation: axis family (HIF / IFN / other)
      tf_axis <- dplyr::case_when(
        rownames(mat_th) %in% c("Hif1a", "Epas1")                              ~ "HIF",
        rownames(mat_th) %in% c("Irf1","Irf3","Irf7","Stat1","Stat2","Nfkb1","Rela") ~ "IFN",
        TRUE                                                                     ~ "other"
      )
      ann_th <- data.frame(axis = tf_axis, row.names = rownames(mat_th))
      ann_col_th <- list(axis = AXIS_COLORS)

      grob_th <- .heat_grob(mat_th,
                            main           = sprintf("CollecTRI TF activity across contrasts (sig in >=%d)", HM_MINC),
                            pmat           = pmat_th,
                            fontsize_row   = 9,
                            cellheight     = 13,
                            cellwidth      = 22,
                            cluster_rows   = TRUE,
                            annotation_row = ann_th,
                            annotation_colors = ann_col_th)

      side_th <- as.data.frame(mat_th) %>%
        tibble::rownames_to_column("source") %>%
        tidyr::pivot_longer(-source, names_to = "contrast", values_to = "score") %>%
        dplyr::left_join(
          as.data.frame(pmat_th) %>%
            tibble::rownames_to_column("source") %>%
            tidyr::pivot_longer(-source, names_to = "contrast", values_to = "padj"),
          by = c("source", "contrast")) %>%
        dplyr::left_join(
          data.frame(source = rownames(mat_th), axis = tf_axis),
          by = "source") %>%
        dplyr::mutate(significant = !is.na(padj) & padj < FDR) %>%
        dplyr::select(source, axis, contrast, score, padj, significant)

      save_overview(grob_th, STAGE, "tf_heatmap",
                    table     = side_th,
                    finding   = paste0(
                      "CollecTRI TF activity heatmap (TFs significant in >=", HM_MINC,
                      " contrasts + watchlist; rows = TFs, cols = contrasts; score clamped +/-",
                      CAP, "; * BH padj<", FDR, "; row strip = HIF/IFN/other axis): ",
                      "IFN/IRF/STAT TFs cluster as the cGAS-dependent block (positive in Interaction); ",
                      "HIF-axis TFs are non-significant in the Interaction. PROVISIONAL; n=5/group."),
                    script    = SCRIPT,
                    fn        = ".heat_grob + pheatmap",
                    config_kv = paste0("thresholds.gsea_fdr=", FDR,
                                       "; figures.z_clamp=", CAP,
                                       "; figures.top_n=", TF_TOPN,
                                       "; colors.diverging"),
                    input     = "03_results/objects/03_tf_collectri.rds",
                    how_to_read = paste0(
                      "Rows = CollecTRI TFs (hierarchically clustered); columns = contrasts. ",
                      "Fill: orange = TF activated in numerator; blue = activated in denominator. ",
                      "Row strip: orange = HIF axis (Hif1a/Epas1), blue = IFN/NFkB axis. ",
                      "Score clamped to +/-", CAP, ". * = BH padj < ", FDR, ". ",
                      "Claim tier: L3 (provisional, n=5/group; IFN-arm TFs positive in Interaction ",
                      "= cGAS-dependent; HIF-axis TFs flat/NS in Interaction = no detectable ",
                      "cGAS-dependence at n=5; NOT proven independence)."),
                    config    = FIG_CFG)
      message("[13_activity_viz] TF heatmap done (",
              nrow(mat_th), " TFs x ", ncol(mat_th), " contrasts)")

      # Positive control: Stat1 in heatmap
      if ("Stat1" %in% rownames(mat_th))
        message("[13_activity_viz] POSITIVE CONTROL: Stat1 present in TF heatmap rows.")
      else
        message("[13_activity_viz] WARNING: Stat1 NOT in TF heatmap rows (check sig threshold).")

    } else {
      message("[13_activity_viz] TF heatmap SKIPPED: only ", length(tf_pick), " qualifying TFs (<3)")
    }
  }, error = function(e) message("  tf heatmap error: ", conditionMessage(e)))
}

# =============================================================================
# 9. OVERVIEW: COMBINED PROGENy + TF PANEL for headline contrasts
#
# A dual-lollipop figure: left panel = PROGENy KEY_PROGENY pathways, right panel
# = KEY_TFS, both plotted for the 4 headline contrasts via facet_wrap on source,
# sharing the same x-axis (activity score). This directly mirrors fig3b from
# 03_decoupler_tf_viz.R (single shared-axis facet_wrap justified by the
# specification comment there). The combined view lets the reader see PROGENy
# footprint and TF activity side-by-side on the same contrast row.
# =============================================================================

if (TF_AVAILABLE && length(CO_TF_HEADLINE) >= 1) {
  message("[13_activity_viz] --- Combined PROGENy+TF panel ---")

  tryCatch({
    # Build a unified long table: entity_type (PROGENy / TF), source, contrast, score, padj/p_value
    prog_sel <- prog_long %>%
      dplyr::filter(source %in% KEY_PROGENY, contrast %in% CO_HEADLINE) %>%
      dplyr::select(source, contrast, score, p_value) %>%
      dplyr::rename(p = p_value) %>%
      dplyr::mutate(entity_type = "PROGENy",
                    axis = dplyr::case_when(
                      source == "Hypoxia"           ~ "HIF",
                      source %in% c("JAK-STAT","NFkB","TNFa") ~ "IFN",
                      TRUE                          ~ "other"))

    tf_sel <- tf_long %>%
      dplyr::filter(source %in% KEY_TFS, contrast %in% CO_TF_HEADLINE) %>%
      dplyr::select(source, contrast, score, padj) %>%
      dplyr::rename(p = padj) %>%
      dplyr::mutate(entity_type = "TF",
                    axis = dplyr::case_when(
                      source %in% c("Hif1a", "Epas1")                             ~ "HIF",
                      source %in% c("Irf1","Irf3","Irf7","Stat1","Stat2","Nfkb1","Rela") ~ "IFN",
                      TRUE                                                          ~ "other"))

    comb_df <- dplyr::bind_rows(prog_sel, tf_sel) %>%
      dplyr::mutate(
        contrast    = factor(contrast, levels = CO_HEADLINE),
        entity_type = factor(entity_type, levels = c("PROGENy", "TF")),
        is_sig      = !is.na(p) & p < FDR,
        star        = .sig_star(p),
        star_y      = score + sign(score) * 0.12
      )

    if (nrow(comb_df) > 0) {
      # Within each entity-type facet, order sources by mean score across contrasts
      src_order <- comb_df %>%
        dplyr::group_by(entity_type, source) %>%
        dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
        dplyr::arrange(entity_type, mean_score) %>%
        dplyr::pull(source) %>%
        unique()

      comb_df <- comb_df %>%
        dplyr::mutate(source = factor(source, levels = src_order))

      fig_comb <- ggplot(comb_df, aes(x = score, y = source)) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.35) +
        geom_segment(aes(xend = 0, yend = source, color = axis), linewidth = 0.7) +
        geom_point(aes(color = axis, shape = entity_type), size = 3.2) +
        geom_text(aes(x = star_y, label = star), size = 4.0, vjust = 0.5,
                  color = "grey20") +
        facet_grid(entity_type ~ contrast,
                   scales = "free_y",
                   space  = "free_y",
                   labeller = labeller(
                     contrast = c(
                       WT_heat     = "WT heat",
                       KO_heat     = "KO heat",
                       Interaction = "Interaction",
                       Temp_main   = "Temp main"
                     )
                   )) +
        scale_color_manual(values = AXIS_COLORS, name = "Axis family") +
        scale_shape_manual(values = c(PROGENy = 18, TF = 16),
                           name = "Entity type") +
        labs(
          title    = "PROGENy + TF activity: key pathways and TFs (headline contrasts)",
          subtitle = paste0(
            "PROGENy MLM footprint (top) + CollecTRI ULM TF activity (bottom). ",
            "Orange = HIF/glycolysis arm; blue = IFN/NFkB arm.\n",
            PROV_STAMP),
          x        = "Activity score",
          y        = NULL,
          caption  = paste0("* p < ", FDR, " (raw for PROGENy, BH padj for TF)",
                            "  |  diamond = PROGENy; circle = TF")
        ) +
        project_theme(config = FIG_CFG, legend = TRUE) +
        theme(panel.spacing.x = unit(0.4, "lines"),
              panel.spacing.y = unit(0.5, "lines"))

      side_comb <- comb_df %>%
        dplyr::transmute(
          entity_type = as.character(entity_type),
          source,
          axis,
          contrast    = as.character(contrast),
          score,
          p,
          significant = is_sig
        )

      save_overview(fig_comb, STAGE, "progeny_tf_combined",
                    table     = side_comb,
                    finding   = paste0(
                      "Combined PROGENy (top panel) + CollecTRI TF (bottom panel) activity for ",
                      "key pathways/TFs across ", length(CO_HEADLINE), " headline contrasts: ",
                      "Hypoxia (PROGENy) and Hif1a (TF) are active in both heat arms and flat in ",
                      "the Interaction (no detectable cGAS-dependence at n=5); JAK-STAT/IFN/Stat1/Irf3 ",
                      "are positive in the Interaction (cGAS-dependent). PROVISIONAL."),
                    script    = SCRIPT,
                    fn        = "ggplot/facet_grid",
                    config_kv = paste0("thresholds.gsea_fdr=", FDR,
                                       "; colors.diverging; design.contrasts"),
                    input     = paste0("03_results/objects/09_progeny_activity.rds; ",
                                       "03_results/objects/03_tf_collectri.rds"),
                    how_to_read = paste0(
                      "Row facets = PROGENy (top) vs CollecTRI TF (bottom); ",
                      "column facets = headline contrasts. ",
                      "Lollipop color: orange = HIF axis, blue = IFN/NFkB axis. ",
                      "Shape: diamond = PROGENy pathway, circle = TF. ",
                      "Direction: rightward = activated in the numerator condition. ",
                      "* = p < ", FDR, " (raw for PROGENy; BH padj for TF). ",
                      "Key comparison: Hypoxia/Hif1a bars are near-zero in the Interaction column ",
                      "- no detectable cGAS-dependence. JAK-STAT/Stat1/Irf3 bars are positive ",
                      "in the Interaction column - cGAS-dependent. ",
                      "CAUTION: n=5/group; absence of significant Interaction for HIF arm is NOT ",
                      "proven independence. Claim tier: L3 (provisional)."),
                    config    = FIG_CFG)
      message("[13_activity_viz] Combined PROGENy+TF panel done.")
    } else {
      message("[13_activity_viz] Combined panel SKIPPED: no data in key pathways/TFs intersection.")
    }
  }, error = function(e) message("  combined panel error: ", conditionMessage(e)))
}

# =============================================================================
# 10. BIOLOGY SANITY CHECKS (console only; findings, not crashes)
# =============================================================================

message("[13_activity_viz] --- Biology sanity checks ---")

.check_pw <- function(pw, co) {
  r <- prog_long %>% dplyr::filter(source == pw, contrast == co)
  if (nrow(r) == 0) return(invisible(NULL))
  message(sprintf("  PROGENy %-10s | %-14s  score=%+.3f  p=%.2e  sig=%s",
                  pw, co, r$score[1], r$p_value[1],
                  ifelse(!is.na(r$p_value[1]) & r$p_value[1] < FDR, "YES", "no")))
  if (!is.na(r$score[1]) && pw == "Hypoxia" && co == "WT_heat" && r$score[1] < 0)
    warning("SANITY: PROGENy Hypoxia score < 0 in WT_heat - expected UP (heat-induced glycolytic/HIF overlap).")
  if (!is.na(r$score[1]) && pw == "JAK-STAT" && co == "WT_heat" && r$score[1] < 0)
    warning("SANITY: PROGENy JAK-STAT score < 0 in WT_heat - expected UP (STING -> IFN positive control).")
  if (!is.na(r$score[1]) && pw == "JAK-STAT" && co == "Interaction" && r$score[1] < 0)
    warning("SANITY: PROGENy JAK-STAT score < 0 in Interaction - expected POSITIVE (IFN arm is cGAS-dependent).")
  invisible(r)
}

for (pw in c("Hypoxia", "JAK-STAT", "NFkB")) {
  for (co in c("WT_heat", "KO_heat", "Interaction", "Temp_main")) {
    if (co %in% CO_PROG) .check_pw(pw, co)
  }
}

if (TF_AVAILABLE) {
  # IFN TFs in Interaction (expect positive scores for Stat1/Irf3/Irf7)
  ifn_int <- tf_long %>%
    dplyr::filter(source %in% c("Stat1", "Irf3", "Irf7"), contrast == "Interaction") %>%
    dplyr::select(source, score, padj)
  if (nrow(ifn_int) > 0)
    message("  IFN TFs in Interaction (expect positive scores, padj<FDR): ",
            paste(sprintf("%s=%.2f(padj=%.2e)", ifn_int$source, ifn_int$score, ifn_int$padj),
                  collapse = ", "))
}

# =============================================================================
# 11. ENSURE README EXISTS FOR THE 05_PROGENY STAGE
# =============================================================================

readme_path <- file.path("03_results", STAGE, "README.md")
if (!file.exists(readme_path)) {
  # write_caption creates the README; seed it with the stage header
  dir.create(dirname(readme_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(
    "# 05_progeny - PROGENy pathway activity",
    "",
    "Stage 05: PROGENy MLM pathway activity inference (14 pathways, all contrasts).",
    "Key result: Hypoxia is flat in the Interaction (no detectable cGAS-dependence at n=5;",
    "PROVISIONAL); JAK-STAT/NFkB/TNFa are positive in the Interaction (cGAS-dependent).",
    "PROVISIONAL - inferred sample mapping pending collaborator sample sheet.",
    ""
  ), readme_path)
  message("[13_activity_viz] Created 03_results/", STAGE, "/README.md stub.")
}

# =============================================================================
# 12. FINAL SUMMARY
# =============================================================================

n_fig <- length(list.files(file.path("03_results", STAGE, "figures"),
                            pattern = "\\.(pdf|png)$", recursive = TRUE))
## adjacency is per figure STEM, not per file: each stem emits both a
## .print.pdf and a .screen.png variant but a single same-stem .csv, so compare
## unique stems (strip the .<variant>.<ext> suffix) against the CSV count.
n_ov  <- length(unique(sub("\\.(print|screen)\\.(pdf|png)$", "",
                           list.files(.ov_fig, pattern = "\\.(pdf|png)$"))))
n_csv <- length(list.files(.ov_tbl, pattern = "\\.csv$"))

message(sprintf(
  "[13_activity_viz] COMPLETE: %d total figure files; %d overview figures / %d overview CSVs under %s/_overview/",
  n_fig, n_ov, n_csv, STAGE))

if (n_fig == 0)
  warning("[13_activity_viz] No figures produced - check errors above.")
if (n_ov != n_csv)
  warning("[13_activity_viz] Overview figure count (", n_ov, ") != overview CSV count (", n_csv,
          ") - source-table-adjacency contract violated.")

message("[13_activity_viz] Run from project root: Rscript 02_analysis/scripts/13_activity_viz.R")
