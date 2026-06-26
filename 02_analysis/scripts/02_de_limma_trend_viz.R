#!/usr/bin/env Rscript
# =============================================================================
# 02_de_limma_trend_viz.R  --  Phase 2 VIZ: cGAS-dependence figure + volcanoes
# =============================================================================
# Phase:    2 (stage 03_de)
# Role:     VISUALIZE half of the "normalize-then-visualize" split. Reads the
#           plot-ready tidy table + the DE checkpoint emitted by
#           02_de_limma_trend.R and renders figures. Performs NO statistics
#           (no lmFit/eBayes/topTable/makeContrasts/p.adjust); it only plots
#           already-computed columns. Runs STANDALONE after the compute script.
#
#           Plot objects come from the RNAseq-toolkit DE plotters
#           (create_standard_volcano, create_MD_plot); every save routes through
#           the figure-style contract (project_theme / save_figure / save_overview).
#           Contrast display labels come from the contract's VECTORIZED
#           contrast_label() (figure_style.R) -- never hardcoded here.
#
# Inputs:   03_results/03_de/tables/fig2_marker_means.csv   (per-group means +
#                                           Interaction logFC/adjP, tidy; also the
#                                           HIGHLIGHT watchlist source for volcano/MD)
#           03_results/objects/02_de_results.rds             (7 topTables, for
#                                           volcano + MD -- existing columns only)
#           03_results/master/master_de_genes.csv            (cross-contrast
#                                           DE accumulator for overview panel)
# Outputs:  03_results/03_de/figures/fig2_cgas_dependence_markers.{print.pdf,screen.png}
#           03_results/03_de/figures/by_contrast/<c>/volcano.{print.pdf,screen.png}  (x7)
#           03_results/03_de/figures/by_contrast/<c>/md.{print.pdf,screen.png}       (x7)
#           03_results/03_de/figures/_overview/de_counts_summary.{print.pdf,screen.png}
#           03_results/03_de/tables/by_contrast/<c>/{volcano,md}.csv
#           03_results/03_de/tables/_overview/de_counts_summary.csv
#           03_results/03_de/README.md                       (captions, idempotent)
#
# NOTE on PCA: the reference (11_de_viz.R §5.6) emits an _overview/pca panel, but
#           that needs a filtered DGEList / logCPM matrix. Only the 7 topTables
#           (02_de_results.rds) exist in 03_results/objects/ -- NO DGEList and NO
#           MArrayLM fit. PCA is therefore SKIPPED (it would force a compute
#           re-run of 02_de_limma_trend.R, which this VIZ-ONLY script must not do).
#
# Dependencies: ggplot2, dplyr, ggrepel, readr, scales, tidyr
# =============================================================================

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme(), save_figure(),
                                                # save_overview(), contrast_path(),
                                                # overview_path(), write_caption(),
                                                # purge_figures(), round_numeric_cols(),
                                                # direction_cue(), VECTORIZED contrast_label()
source("02_analysis/config/config.R")          # YAML_CONFIG, DIR_OBJECTS, DIR_MASTER,
                                                # DE_FDR, DE_LOGFC, provisional_caption(),
                                                # stage_dir()

## RNAseq-toolkit DE plotters (source-only library; theme NOT required here --
## project_theme() is appended after each call as the single house-style entry point).
source("01_modules/RNAseq-toolkit/scripts/DE/plot_standard_volcano.R")  # create_standard_volcano()
source("01_modules/RNAseq-toolkit/scripts/DE/create_MD_plot.R")         # create_MD_plot()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(readr)
  library(scales)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# 1. Style-contract constants (NEVER hardcode; always from FIG_CFG / config.R)
# -----------------------------------------------------------------------------
STAGE   <- "03_de"
SCRIPT  <- "02_analysis/scripts/02_de_limma_trend_viz.R"

FDR     <- as.numeric(FIG_CFG$thresholds$de_fdr   %||% DE_FDR)     # 0.05
LFC     <- as.numeric(FIG_CFG$thresholds$de_logfc %||% DE_LOGFC)   # 1.0
LBL_TOP <- as.integer(FIG_CFG$figures$volcano_label_top %||% 10L)  # 10
CUE     <- as.numeric(FIG_CFG$figures$cue_size %||% 4)             # in-panel cue text size

NEG <- FIG_CFG$colors$diverging$down    %||% "#2166AC"   # blue  (down in numerator)
MID <- FIG_CFG$colors$diverging$neutral %||% "#F7F7F7"   # white
POS <- FIG_CFG$colors$diverging$up      %||% "#B35806"   # orange (up in numerator)

# All 7 contrasts in config order — NEVER hardcode names
CONTRASTS <- vapply(YAML_CONFIG$design$contrasts, function(x) x$name, character(1))

CFG_KV <- sprintf(
  "thresholds.de_fdr=%.2g; thresholds.de_logfc=%.1f; figures.volcano_label_top=%d; colors.diverging",
  FDR, LFC, LBL_TOP
)

# -----------------------------------------------------------------------------
# 2. LOCAL helpers (interpretive caption notes + dynamic x-tick step)
#    NOTE: contrast DISPLAY labels come from the contract's VECTORIZED
#    contrast_label() (figure_style.R). We do NOT redefine it here. The scalar
#    shadow that used to live in this block (with a hardcoded inline LUT) was a
#    direct contract violation and has been DELETED.
# -----------------------------------------------------------------------------

#' Contrast-specific interpretive note for caption `finding` strings ONLY
#' (never a plot title). Preserves the two-arms framing constraint:
#' a non-significant Interaction gene = "no detectable cGAS-dependence at n=5",
#' NEVER "cGAS-independent".
interp_note <- function(co) {
  switch(co,
    Interaction = paste0(
      " Interaction is THE cGAS-dependence test (1 df, n=5, lowest-powered term); ",
      "a non-significant gene = no detectable cGAS-dependence at n=5, NOT independence."),
    Temp_main = paste0(
      " Pooled/marginal heat effect (collapses genotype); ",
      "NOT the cGAS-dependence read-out (that is the Interaction)."),
    Geno_main = paste0(
      " Pooled/marginal genotype effect (collapses temperature); ",
      "NOT the cGAS-dependence read-out (that is the Interaction)."),
    ""
  )
}

#' Per-panel subtitle: contract label + direction key (+ the Interaction caveat).
#' Direction text is generic (no n-badge) -- the sig counts are annotated inside
#' the panel via direction_cue glyphs (mirrors reference 11_de_viz.R §5.4).
panel_subtitle <- function(co) {
  paste0(contrast_label(co),
         "  |  logFC>0 = up in numerator; logFC<0 = up in denominator",
         if (co == "Interaction")
           "  |  non-sig = no detectable cGAS-dependence at n=5"
         else "")
}

#' Pick a "nice" x-break step (1-2-5 ladder) targeting ~10 ticks over the full
#' symmetric ±range so volcano x-axis labels never crowd. LOCAL (ref §5.3a).
.volcano_x_step <- function(x_range) {
  n_target_half <- 5L
  raw_step   <- x_range / n_target_half
  candidates <- c(0.5, 1, 2, 5, 10, 20, 50)
  step <- candidates[candidates >= raw_step][1]
  if (is.na(step)) step <- ceiling(raw_step)
  step
}

# -----------------------------------------------------------------------------
# 3. Read plot-ready inputs (NO recomputation)
# -----------------------------------------------------------------------------
fig2_means_path <- file.path(stage_dir("03_de", "tables"), "fig2_marker_means.csv")
if (!file.exists(fig2_means_path)) {
  stop(sprintf("Missing %s -- run 02_de_limma_trend.R (compute) first.", fig2_means_path))
}
gm_long <- read.csv(fig2_means_path, stringsAsFactors = FALSE, check.names = FALSE)

de_path <- file.path(DIR_OBJECTS, "02_de_results.rds")
if (!file.exists(de_path)) {
  stop(sprintf("Missing %s -- run 02_de_limma_trend.R (compute) first.", de_path))
}
de_results <- readRDS(de_path)

# HIGHLIGHT watchlist: config has no key_genes block, so derive the cGAS-STING /
# two-arms biology watchlist from the fig2 marker genes (7 ISG-arm + 8 HIF/glyco-arm).
# These are exactly the genes the reference would always label on every volcano/MD.
HIGHLIGHT <- unique(gm_long$gene)
stopifnot(length(HIGHLIGHT) >= 1)
message(sprintf("[02_viz] HIGHLIGHT watchlist (%d two-arms marker genes): %s",
                length(HIGHLIGHT), paste(HIGHLIGHT, collapse = ", ")))

# -----------------------------------------------------------------------------
# 5. Read-only source-table captions (fig2_marker_means.csv + marker_cgas_dependence.csv)
# -----------------------------------------------------------------------------
write_caption(
  stage      = STAGE,
  filename   = "tables/fig2_marker_means.csv",
  finding    = paste0(
    "Per-group mean log2(CPM+0.5) tidy table for 15 cGAS-dependence marker genes ",
    "(7 ISG-arm: Ifit1, Isg15, Irf7, Oasl2, Mx1, Stat1, Cxcl10; ",
    "8 HIF/glycolysis-arm: Slc2a1, Vegfa, Egln3, Bnip3, Pgk1, Ldha, Aldoa, Hk2) ",
    "across 4 groups (WT_37, WT_39, cGASKO_37, cGASKO_39). ",
    "Also carries Interaction logFC and adj.P.Val from limma-trend (no recomputation in viz). ",
    "READ-ONLY in this script; schema-frozen (also consumed by 03_decoupler_tf.R:687). ",
    "Also the HIGHLIGHT watchlist source for the volcano + MD sweeps. ",
    "Claim tier: L3. PROVISIONAL sample labels."),
  script     = SCRIPT,
  fn         = "read.csv (read-only; produced by 02_de_limma_trend.R)",
  config_kv  = CFG_KV,
  input      = "03_results/objects/02_de_results.rds + fig2_marker_means.csv",
  how_to_read = paste0(
    "Columns: gene, arm (IFN/ISG arm | HIF/glycolysis arm), hif_class (HIF-specific | shared-glycolytic | NA for ISG), ",
    "genotype (WT | cGASKO), temp (37C | 39C), mean_log2cpm, inter_logFC, inter_adjP. ",
    "One row per gene x group (4 rows per gene). inter_adjP is constant within gene. ",
    "Schema frozen: do NOT rename/relocate; consumed by downstream scripts. Claim tier: L3. PROVISIONAL."),
  config     = FIG_CFG
)

write_caption(
  stage      = STAGE,
  filename   = "tables/marker_cgas_dependence.csv",
  finding    = paste0(
    "Per-marker Interaction contrast statistics table (subset of 02_de_results.rds, Interaction slot) ",
    "for the cGAS-dependence marker genes. Lists logFC, P.Value, adj.P.Val, AveExpr for each marker. ",
    "Used as supplemental evidence for the two-arms framing. Claim tier: L3. PROVISIONAL."),
  script     = SCRIPT,
  fn         = "produced by 02_de_limma_trend.R (read-only in viz)",
  config_kv  = CFG_KV,
  input      = "03_results/objects/02_de_results.rds",
  how_to_read = paste0(
    "Rows = marker genes; columns = gene_symbol, ensembl, logFC, AveExpr, t, P.Value, adj.P.Val, contrast. ",
    "All values from the Interaction contrast (1 df; logFC > 0 = up in WT relative to cGASKO after removing shared heat effect). ",
    "Non-significant adj.P = no detectable cGAS-dependence at n=5 (NEVER 'cGAS-independent'). ",
    "Claim tier: L3. PROVISIONAL sample labels."),
  config     = FIG_CFG
)

# -----------------------------------------------------------------------------
# 6. Purge stale flat figures (the redundant flat volcano trio + any old fig2)
#    The flat volcano_{WT_heat,KO_heat,Interaction}.pdf are now SUPERSEDED by the
#    full by_contrast/<c>/volcano sweep (§7) for all 7 contrasts. Remove them so
#    the run owns its figure namespace (avoids stale off-contract leftovers).
# -----------------------------------------------------------------------------
purge_figures(STAGE, "volcano_", config = FIG_CFG)

# -----------------------------------------------------------------------------
# 7. Per-contrast volcano sweep over ALL 7 contrasts (toolkit + figure-style)
# -----------------------------------------------------------------------------
# Plot object built by the RNAseq-toolkit create_standard_volcano() (decision-by-
# FDR; raw-p on y, FDR for the color/label decision), then styled via project_theme
# and saved through save_overview (figure + same-stem table + README caption).
# The toolkit helper does its own thresholding for COLORING/LABELING only -- no DE
# statistics are (re)computed here. Symbol rownames are required (labels die on
# Ensembl ids); we set them via make.unique(gene_symbol).

message(sprintf("[02_viz] Per-contrast volcano sweep: %d contrasts", length(CONTRASTS)))

for (co in CONTRASTS) {
  tt <- de_results[[co]]
  if (is.null(tt) || nrow(tt) == 0) {
    message(sprintf("[02_viz] SKIP %s (NULL or 0-row topTable)", co))
    next
  }

  # Symbol rownames so create_standard_volcano() labels render (not Ensembl ids).
  rownames(tt) <- make.unique(as.character(tt$gene_symbol))

  # Shared significance rule (matches reference sig_of): FDR AND |logFC|.
  sig    <- (tt$adj.P.Val < FDR) & (abs(tt$logFC) >= LFC)
  n_up   <- sum(sig & tt$logFC >  0, na.rm = TRUE)
  n_down <- sum(sig & tt$logFC <  0, na.rm = TRUE)
  n_sig  <- n_up + n_down

  message(sprintf("[02_viz] %s: %d sig (%d up-in-numerator, %d up-in-denominator)",
                  co, n_sig, n_up, n_down))

  # Dynamic x-tick step: ~10 labelled ticks across the symmetric axis range.
  x_step <- .volcano_x_step(max(abs(tt$logFC), na.rm = TRUE))

  p_volc <- create_standard_volcano(
    de_results    = tt[, c("logFC", "P.Value", "adj.P.Val")],
    decision_by   = "fdr",
    p_cutoff      = FDR,
    fc_cutoff     = LFC,
    top_n         = LBL_TOP,
    highlight_gene = HIGHLIGHT,         # two-arms ISG/HIF watchlist always labelled
    title         = contrast_label(co),
    subtitle      = panel_subtitle(co),
    x_breaks      = x_step,
    max.overlaps  = 15,
    annotate_counts = TRUE             # toolkit appends sig ↑/↓ counts under the
                                       # highest-priority populated legend line
  ) +
    labs(caption = NULL) +
    project_theme(config = FIG_CFG) +
    # Give the appended count line ("↑ n  ↓ m") a little air below its sig line —
    # the legend label's second line was sitting too tight against the first.
    # Volcano-only (merges onto project_theme's legend.text); leaves single-line
    # and running-sum wrapped legends unaffected.
    theme(legend.text = element_text(lineheight = 1.5))

  # Source table: key columns, ordered by significance then effect size.
  tbl_volc <- tt[order(tt$adj.P.Val, -abs(tt$logFC)),
                 intersect(c("gene_symbol", "ensembl", "logFC", "AveExpr",
                             "t", "P.Value", "adj.P.Val", "contrast"),
                           names(tt)),
                 drop = FALSE]

  save_overview(
    plot      = p_volc,
    stage     = STAGE,
    name      = "volcano",
    table     = tbl_volc,
    contrast  = co,
    finding   = sprintf(
      "limma-trend volcano for %s (%s): %d genes pass adj.P < %.2g & |log2FC| >= %.1f (%d up in numerator, %d up in denominator). cGAS-STING/ISG + HIF watchlist genes always labelled; sig up/down counts appended (toolkit) as a second line under the highest-priority populated significance line in the legend; x-tick step chosen dynamically to avoid label crowding.%s",
      co, contrast_label(co), n_sig, FDR, LFC, n_up, n_down, interp_note(co)
    ),
    script    = SCRIPT,
    fn        = "create_standard_volcano (toolkit)",
    config_kv = CFG_KV,
    input     = "03_results/objects/02_de_results.rds",
    how_to_read = sprintf(
      paste0(
        "x = log2FC (>0 = up in numerator; <0 = up in denominator); ",
        "y = -log10(P.Value) (raw p on y; FDR for the color decision). ",
        "Dashed lines = adj.P < %.2g boundary (horizontal) & |log2FC| >= %.1f (vertical). ",
        "Colored points = significant by the toolkit decision-by-FDR rule. ",
        "Labelled = top %d genes by significance per side, PLUS the always-on two-arms ",
        "watchlist (ISG: Ifit1/Isg15/Irf7/...; HIF/glyco: Slc2a1/Vegfa/Egln3/...). ",
        "Data: CPM / limma-trend on log2(CPM+0.5) (no voom). Claim tier: L3. ",
        "PROVISIONAL sample labels; n=5/group."
      ),
      FDR, LFC, LBL_TOP
    ),
    config    = FIG_CFG,
    # Volcano carries ~30 repel labels (top-N per side + 15-gene watchlist); the
    # default 3.5x3in print column is too small for ggrepel (zero-dimension viewport).
    # Give both variants a roomier square canvas so labels fit and stay legible.
    width = 9, height = 8
  )
}

# -----------------------------------------------------------------------------
# 8. Per-contrast MD (mean-difference) sweep over ALL 7 contrasts (NEW)
# -----------------------------------------------------------------------------
# Mirrors reference 11_de_viz.R §5.5 (create_MD_plot). The reference passes the
# MArrayLM `fit`, but NO fit object exists in 03_results/objects/ (VIZ-ONLY; we
# must not re-run the compute script to emit one). The topTables already carry
# AveExpr + logFC + adj.P.Val -- exactly what create_MD_plot needs -- so we drive
# it from the topTable directly with fit = NULL (the fit is only dereferenced for
# a numeric coef or the AveExpr fallback, neither of which we hit: coef is a name
# and AveExpr is present). Same toolkit look (NS cloud + colored sig dots + LOESS
# trend + median guide + quadrant counts), styled via project_theme, saved via
# save_overview.

message(sprintf("[02_viz] Per-contrast MD sweep: %d contrasts", length(CONTRASTS)))

for (co in CONTRASTS) {
  tt <- de_results[[co]]
  if (is.null(tt) || nrow(tt) == 0) {
    message(sprintf("[02_viz] SKIP MD %s (NULL or 0-row topTable)", co))
    next
  }
  rownames(tt) <- make.unique(as.character(tt$gene_symbol))

  p_md <- create_MD_plot(
    fit            = NULL,              # no MArrayLM on disk; topTable carries AveExpr
    coef           = co,               # character -> fit never dereferenced
    de_results     = tt[, c("logFC", "AveExpr", "P.Value", "adj.P.Val")],
    fc_cutoff      = LFC,
    fdr_cutoff     = FDR,
    top_n          = 5,
    highlight_gene = HIGHLIGHT,
    label_method   = "top",
    show_quadrant_counts = TRUE,
    title          = contrast_label(co)
  ) +
    labs(subtitle = panel_subtitle(co), caption = NULL) +
    project_theme(config = FIG_CFG)

  # Source table: same key columns as the volcano (significance-ordered).
  tbl_md <- tt[order(tt$adj.P.Val, -abs(tt$logFC)),
               intersect(c("gene_symbol", "ensembl", "logFC", "AveExpr",
                           "t", "P.Value", "adj.P.Val", "contrast"),
                         names(tt)),
               drop = FALSE]

  save_overview(
    plot      = p_md,
    stage     = STAGE,
    name      = "md",
    table     = tbl_md,
    contrast  = co,
    finding   = sprintf(
      "limma-trend mean-difference (MD) plot for %s (%s): logFC vs average expression. LOESS trend (red dashed) + median-expression guide expose any expression-range normalization bias; quadrant counts quantify directionality. cGAS-STING/ISG + HIF watchlist genes always labelled.%s",
      co, contrast_label(co), interp_note(co)
    ),
    script    = SCRIPT,
    fn        = "create_MD_plot (toolkit; driven from topTable AveExpr/logFC)",
    config_kv = CFG_KV,
    input     = "03_results/objects/02_de_results.rds",
    how_to_read = sprintf(
      paste0(
        "x = average expression (AveExpr, log2CPM scale); y = log2FC. ",
        "Grey cloud = NS genes; orange = up (FDR < %.2g); blue = down (FDR < %.2g). ",
        "Red dashed = LOESS trend (should hug logFC~0; a tilt flags a normalization confound); ",
        "dotted vertical = median expression; dashed horizontals = |log2FC| >= %.1f. ",
        "Quadrant numbers = significant-gene counts per quadrant. ",
        "Labelled = top 5 sig genes per side, PLUS the always-on two-arms watchlist. ",
        "Driven from 02_de_results.rds topTable columns (no MArrayLM fit on disk). ",
        "Claim tier: L3. PROVISIONAL sample labels; n=5/group."
      ),
      FDR, FDR, LFC
    ),
    config    = FIG_CFG,
    # MD also repel-labels the watchlist + top genes; widen the print canvas.
    width = 9, height = 8
  )
}

# -----------------------------------------------------------------------------
# 9. Cross-contrast DE counts summary (_overview/)
# -----------------------------------------------------------------------------
# Reads master_de_genes.csv (significant flag = adj.P < 0.05 ONLY, no logFC gate
# -- different from the volcano combined threshold; noted explicitly in caption).
# Signed bar chart: n_up positive/orange, -n_down negative/blue.

master_de_path <- file.path(DIR_MASTER, "master_de_genes.csv")

if (!file.exists(master_de_path)) {
  message(sprintf(
    "[02_viz] SKIP de_counts_summary: master_de_genes.csv not found at %s\n  Run 02_de_limma_trend.R first.",
    master_de_path
  ))
} else {
  mde <- readr::read_csv(master_de_path, show_col_types = FALSE)

  req_cols_mde <- c("significant", "direction", "contrast")
  missing_mde  <- setdiff(req_cols_mde, colnames(mde))
  if (length(missing_mde) > 0) {
    message(sprintf(
      "[02_viz] SKIP de_counts_summary: master_de_genes.csv missing columns: %s",
      paste(missing_mde, collapse = ", ")
    ))
  } else {
    # Restrict to config contrasts present in master (others may not have run)
    contrasts_avail <- intersect(CONTRASTS, unique(mde$contrast))
    if (length(contrasts_avail) == 0) {
      message("[02_viz] SKIP de_counts_summary: no config contrasts in master_de_genes.csv")
    } else {
      # NOTE: master `significant` = adj.P < 0.05 ONLY (no logFC gate).
      # This is DELIBERATELY different from the volcano combined threshold
      # (adj.P < FDR AND |logFC| >= LFC); the discrepancy is explained in captions.
      de_counts <- mde %>%
        dplyr::filter(contrast %in% contrasts_avail) %>%
        dplyr::group_by(contrast) %>%
        dplyr::summarise(
          n_tested = dplyr::n(),
          n_sig    = sum(significant, na.rm = TRUE),
          n_up     = sum(significant & direction == "Up",   na.rm = TRUE),
          n_down   = sum(significant & direction == "Down", na.rm = TRUE),
          .groups  = "drop"
        ) %>%
        # Factor in config order for stable facet/axis ordering
        dplyr::mutate(contrast = factor(contrast, levels = contrasts_avail))

      # Signed bar chart: up positive (orange), -down negative (blue)
      de_counts_long <- de_counts %>%
        dplyr::select(contrast, n_up, n_down) %>%
        tidyr::pivot_longer(
          cols      = c(n_up, n_down),
          names_to  = "direction",
          values_to = "n"
        ) %>%
        dplyr::mutate(
          n_signed  = ifelse(direction == "n_down", -n, n),
          direction = factor(direction,
                             levels = c("n_up", "n_down"),
                             labels = c("Up in numerator", "Down in numerator"))
        )

      y_abs <- max(abs(de_counts_long$n_signed), na.rm = TRUE)
      y_abs <- max(y_abs, 1L)   # guard against all-zero data

      p_de_counts <- ggplot(de_counts_long,
                            aes(x = contrast, y = n_signed, fill = direction)) +
        geom_col(position = "stack", width = 0.8) +
        geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
        scale_x_discrete(labels = function(x) contrast_label(x, short = TRUE)) +
        scale_fill_manual(
          values = c("Up in numerator"   = POS,
                     "Down in numerator" = NEG),
          name   = "Direction"
        ) +
        scale_y_continuous(
          limits = c(-y_abs, y_abs),
          oob    = scales::squish,
          labels = function(x) as.character(abs(x))  # display absolute values on axis
        ) +
        labs(
          title    = sprintf("DE gene counts per contrast (adj.P < %.2g; no |logFC| gate)", FDR),
          subtitle = NULL,
          x        = NULL,
          y        = sprintf("Gene count (adj.P < %.2g; |bars| = absolute count)", FDR),
          caption  = NULL
        ) +
        project_theme(config = FIG_CFG) +
        ggplot2::theme(
          # Short labels are already two-line wrapped (contrast_labels_short carries
          # an embedded newline); keep them HORIZONTAL so the 7 wrapped ticks do not
          # collide (a 45deg rotation of two-line labels overlaps badly).
          axis.text.x     = ggplot2::element_text(angle = 0, hjust = 0.5, size = 10),
          legend.position = "bottom"
        )

      save_overview(
        plot      = p_de_counts,
        stage     = STAGE,
        name      = "de_counts_summary",
        table     = de_counts,
        # no contrast= argument -> routes to _overview/
        finding   = sprintf(
          paste0(
            "DE gene counts per contrast (adj.P < %.2g; NO logFC gate) across %d contrasts. ",
            "Signed bar chart: positive orange bars = genes up in numerator; negative blue bars = genes up in denominator. ",
            "IMPORTANT: the master significant flag = adj.P < %.2g ONLY (no |logFC| gate), unlike the ",
            "per-contrast volcano panels which combine adj.P < %.2g AND |logFC| >= %.1f; ",
            "volcano counts will therefore be lower than the summary bars here. ",
            "Interaction contrast (1 df, n=5) is expected to have the fewest significant genes ",
            "(lowest power): non-significant gene = no detectable cGAS-dependence at n=5, NOT independence."
          ),
          FDR, length(contrasts_avail), FDR, FDR, LFC
        ),
        script    = SCRIPT,
        fn        = "geom_col (signed bar chart)",
        config_kv = CFG_KV,
        input     = "03_results/master/master_de_genes.csv",
        how_to_read = sprintf(
          paste0(
            "x = contrast (config order; short display labels); y = gene count (absolute value on axis). ",
            "Orange bars above zero = up in numerator side; blue bars below zero = up in denominator side. ",
            "Threshold for this panel: adj.P < %.2g (no logFC filter). ",
            "Threshold for by_contrast/*/volcano panels: adj.P < %.2g AND |logFC| >= %.1f -- ",
            "the panels are complementary, NOT contradictory. ",
            "Interaction contrast is the cGAS-dependence payoff: 1 df, n=5, lowest power. ",
            "Claim tier: L3. PROVISIONAL sample labels."
          ),
          FDR, FDR, LFC
        ),
        config    = FIG_CFG,
        # Wide canvas so the 7 two-line short contrast labels sit side-by-side
        # without colliding (print gets a roomier square too).
        wide = TRUE, width = 12, height = 8
      )

      message("[02_viz] de_counts_summary overview emitted")
    }
  }
}

# -----------------------------------------------------------------------------
# 10. Final structural asserts (volcano + MD per contrast; LOUD, no tryCatch)
# -----------------------------------------------------------------------------

bc_dir <- FIG_CFG$figures$by_contrast_dir %||% "by_contrast"

expected_volc_pngs <- vapply(CONTRASTS, function(co) {
  file.path("03_results", STAGE, "figures", bc_dir, co, "volcano.png")
}, character(1))
expected_md_pngs <- vapply(CONTRASTS, function(co) {
  file.path("03_results", STAGE, "figures", bc_dir, co, "md.png")
}, character(1))

missing_volc <- expected_volc_pngs[!file.exists(expected_volc_pngs)]
missing_md   <- expected_md_pngs[!file.exists(expected_md_pngs)]

if (length(missing_volc) > 0) {
  warning(sprintf("[02_viz] %d expected by_contrast volcano PNG(s) missing:\n  %s",
                  length(missing_volc), paste(missing_volc, collapse = "\n  ")))
} else {
  message(sprintf("[02_viz] All %d by_contrast volcano.png present.", length(CONTRASTS)))
}
if (length(missing_md) > 0) {
  warning(sprintf("[02_viz] %d expected by_contrast md PNG(s) missing:\n  %s",
                  length(missing_md), paste(missing_md, collapse = "\n  ")))
} else {
  message(sprintf("[02_viz] All %d by_contrast md.png present.", length(CONTRASTS)))
}

n_pdf <- length(list.files(file.path("03_results", STAGE, "figures"),
                            pattern = "\\.(pdf|png)$", recursive = TRUE))
n_csv <- length(list.files(file.path("03_results", STAGE, "tables"),
                            pattern = "\\.csv$",       recursive = TRUE))

message(sprintf(
  "\n02_de_limma_trend_viz.R complete.\n  Figures: %d PDF/PNG files under 03_results/%s/figures/\n  Tables:  %d CSV files under 03_results/%s/tables/\n  Sweep:   by_contrast/<c>/{volcano,md}.{pdf,png} for %d contrasts\n  Overview: _overview/de_counts_summary.{pdf,png}\n  PCA: SKIPPED (no DGEList/MArrayLM in 03_results/objects/; viz-only, no compute re-run).",
  n_pdf, STAGE, n_csv, STAGE, length(CONTRASTS)
))

stopifnot("No output figures produced" = n_pdf >= 1)
