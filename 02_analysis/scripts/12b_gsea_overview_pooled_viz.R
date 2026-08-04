# 12b_gsea_overview_pooled_viz.R — VIZ
## Faceted lollipop overview of the GSEA sweep under the pooled Benjamini-Hochberg
## FDR family computed by 06b_gsea_pooled_padj.R. One facet per displayed database,
## the top-N sets inside each facet, NES on x, pooled-FDR status in the glyph shape.
##
## Nothing here recomputes a statistic: the pooled family, its padj_pooled column and
## the per-database significance counts all arrive from 06b. This script selects what
## to draw, orders it, and renders it.
##
## The facet list, the per-facet cap and the canvas geometry are read from
## analysis_config.yaml (`gsea_pooled_overview:`) — see that block for why
## exclude_databases drops a database from the panel while keeping it in the pooled family.
##
## Reads:
##   03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv
##
## Emits via save_overview():
##   03_results/06_gsea/figures/_overview/gsea_pooled_overview_WT_heat.{pdf,png}
##   03_results/06_gsea/figures/_overview/gsea_pooled_overview_Interaction.{pdf,png}
##   03_results/06_gsea/figures/_overview/gsea_pooled_overview.{pdf,png}
##   and sibling CSV tables in 03_results/06_gsea/tables/_overview/

source("02_analysis/helpers/figure_style.R")
source("02_analysis/config/config.R")

.RTK <- "01_modules/RNAseq-toolkit"
.GP  <- file.path(.RTK, "scripts", "GSEA", "GSEA_plotting")
source(file.path(.GP, "format_pathway_names.R"))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(grid)
})

STAGE  <- "06_gsea"
SCRIPT <- "02_analysis/scripts/12b_gsea_overview_pooled_viz.R"
FDR    <- GSEA_FDR_CUTOFF %||% 0.05

# --- Panel contract (analysis_config.yaml :: gsea_pooled_overview) ------------
OVR        <- FIG_CFG$gsea_pooled_overview %||% list()
TOP_N      <- as.integer(OVR$top_n_per_db      %||% 5)
DB_EXCLUDE <- as.character(unlist(OVR$exclude_databases %||% character(0)))
LABEL_WRAP <- as.integer(OVR$label_wrap        %||% 46)
FACET_NCOL <- as.integer(OVR$facet_ncol        %||% 2)
ROW_H      <- as.numeric(OVR$facet_row_height  %||% 0.30)
STRIP_H    <- as.numeric(OVR$facet_strip_height %||% 0.34)
CHROME_H   <- as.numeric(OVR$chrome_height     %||% 3.4)

CANVAS_W <- as.numeric(FIG_CFG$figures$width_wide %||% 13)   # save_overview(wide = TRUE)

NEG <- FIG_CFG$colors$diverging$down    %||% "steelblue4"
MID <- FIG_CFG$colors$diverging$neutral %||% "grey97"
POS <- FIG_CFG$colors$diverging$up      %||% "sienna"

# Database display order = the order the databases are declared in the config
# `databases:` block, minus the panel exclusions. No facet list is hardcoded here.
.db_names <- function(x) vapply(x, function(d) as.character(d$name), character(1))
DB_ALL    <- c(.db_names(FIG_CFG$databases$msigdb), .db_names(FIG_CFG$databases$custom))
DB_SHOWN  <- setdiff(DB_ALL, DB_EXCLUDE)

# Pooling-FDR status dictionary. All four levels always appear in the legend (a
# stable glyph dictionary across contrasts), which is why the plot carries a
# transparent one-row-per-level key layer: since ggplot2 3.5 a legend key whose
# level has no rows in any layer is drawn as an empty slot with no glyph.
SIG_LEVELS   <- c("Significant (both)", "Lost on pooling", "Gained on pooling", "Not significant")
SIG_SHAPES   <- stats::setNames(c(21, 23, 24, 1), SIG_LEVELS)      # circle / diamond / triangle / open circle
SIG_KEY_FILL <- c("grey35", "grey35", "grey35", NA)                # neutral stand-in for the NES fill

# The ONLY prose kept on the canvas: it decodes two glyph conventions the legend
# cannot carry by itself. Everything else that used to sit in a footnote block now
# lives in the stage README caption.
FACE_NOTE <- paste0(
  "Facet header (n/N): n of the N sets in that database pass the pooled FDR, counted over the whole database. ",
  "Panel fill is the NES colour; the legend keys use a neutral fill."
)
KEY_SEP      <- "\u001F"   # unit separator: makes the y key unique per facet

# Load master pooled GSEA table from compute step (06b)
fp_pooled <- "03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv"
if (!file.exists(fp_pooled)) {
  stop("Master pooled table missing at ", fp_pooled, ". Run 06b_gsea_pooled_padj.R first.")
}

df_all <- readr::read_csv(fp_pooled, show_col_types = FALSE)

missing_db <- setdiff(DB_SHOWN, unique(df_all$database))
if (length(missing_db)) {
  warning("Configured databases absent from the pooled table: ", paste(missing_db, collapse = ", "))
  DB_SHOWN <- intersect(DB_SHOWN, unique(df_all$database))
}

#' Measure the rendered title against the canvas and warn if it will run off.
#' Measured on the built gtable, not assumed from the character count.
.title_fits <- function(p, canvas_w) {
  gt  <- ggplot2::ggplotGrob(p)
  idx <- which(gt$layout$name == "title")
  if (!length(idx)) return(invisible(NA_real_))
  w <- grid::convertWidth(grid::grobWidth(gt$grobs[[idx[1]]]), "in", valueOnly = TRUE)
  if (w > canvas_w - 0.4) {
    warning(sprintf("Title is %.2f in wide on a %.2f in canvas and will be clipped.", w, canvas_w))
  } else {
    message(sprintf("  [12b] title width %.2f in of %.2f in canvas", w, canvas_w))
  }
  invisible(w)
}

#' Panel block height (inches) for one faceted contrast panel.
.panel_block_h <- function(n_db) ceiling(n_db / FACET_NCOL) * (STRIP_H + TOP_N * ROW_H)

# ---------------------------------------------------------------------------
# Per-contrast panel
# ---------------------------------------------------------------------------
make_contrast_overview_plot <- function(df_sub, contrast_name) {
  # Short, gloss-free contrast name from the config label map (newline -> space).
  co_short    <- gsub("\n", " ", contrast_label(contrast_name, short = TRUE), fixed = TRUE)
  n_db_pooled <- dplyr::n_distinct(df_sub$database)   # the BH family, all databases

  # Significance category:
  #   "Significant (both)" : padj < FDR & padj_pooled < FDR
  #   "Lost on pooling"    : padj < FDR & padj_pooled >= FDR
  #   "Gained on pooling"  : padj >= FDR & padj_pooled < FDR
  #   "Not significant"    : neither
  df_mod <- df_sub %>%
    dplyr::filter(database %in% DB_SHOWN) %>%
    dplyr::mutate(
      sig_status = dplyr::case_when(
        padj <  FDR & padj_pooled <  FDR ~ SIG_LEVELS[1],
        padj <  FDR & padj_pooled >= FDR ~ SIG_LEVELS[2],
        padj >= FDR & padj_pooled <  FDR ~ SIG_LEVELS[3],
        TRUE                             ~ SIG_LEVELS[4]
      ),
      sig_status     = factor(sig_status, levels = SIG_LEVELS),
      database       = factor(database, levels = DB_SHOWN),
      formatted_name = format_pathway_name(pathway_name, use_formatting = TRUE, strip_prefix = TRUE)
    )
  n_db_shown <- dplyr::n_distinct(df_mod$database)

  # SELECTION RULE, in full: the top TOP_N sets inside each database ranked by
  # padj_pooled, ties broken by the raw p-value and then by |NES| so the draw is
  # deterministic. Nothing else is added — no thematic pin, no manual rescue.
  plot_data <- df_mod %>%
    dplyr::group_by(database) %>%
    dplyr::arrange(padj_pooled, pvalue, dplyr::desc(abs(nes)), .by_group = TRUE) %>%
    dplyr::slice_head(n = TOP_N) %>%
    dplyr::ungroup()

  # Facet header carries each database's pooled-significant count over the WHOLE
  # database, so a facet whose drawn sets are all non-significant says so in its header
  # instead of costing a facet to say nothing.
  db_meta <- df_mod %>%
    dplyr::distinct(database, n_tests_in_db, sig_after_in_db) %>%
    dplyr::mutate(db_label = sprintf("%s (%d/%d)",
                                     as.character(database),
                                     as.integer(sig_after_in_db),
                                     as.integer(n_tests_in_db))) %>%
    dplyr::arrange(database)
  db_levels <- db_meta$db_label

  # y key must be unique ACROSS facets: two databases can format to the same set
  # name (e.g. "Inflammatory Response" in Hallmark and GO_BP), and a shared factor
  # level would draw one row in both facets.
  plot_data <- plot_data %>%
    dplyr::left_join(db_meta %>% dplyr::select(database, db_label), by = "database") %>%
    dplyr::mutate(db_label = factor(db_label, levels = db_levels)) %>%
    dplyr::arrange(database, nes) %>%
    dplyr::mutate(row_key = paste(as.integer(database), formatted_name, sep = KEY_SEP),
                  row_key = factor(row_key, levels = unique(row_key)))

  # One transparent point per status level so every legend key has a glyph.
  key_stub <- data.frame(
    sig_status = factor(SIG_LEVELS, levels = SIG_LEVELS),
    db_label   = factor(db_levels[1], levels = db_levels),
    row_key    = factor(levels(plot_data$row_key)[1], levels = levels(plot_data$row_key)),
    nes        = 0
  )

  n_total  <- df_sub$n_tests_pooled[1]
  n_before <- sum(df_sub$padj < FDR, na.rm = TRUE)
  n_after  <- sum(df_sub$padj_pooled < FDR, na.rm = TRUE)
  n_lost   <- sum(df_sub$padj < FDR & df_sub$padj_pooled >= FDR, na.rm = TRUE)
  n_gained <- sum(df_sub$padj >= FDR & df_sub$padj_pooled < FDR, na.rm = TRUE)

  sub_text <- sprintf(
    paste0("Top %d sets per database, ranked by pooled FDR (ties broken by raw p, then |NES|); ",
           "%d of the %d pooled databases shown.\n",
           "Pooled family: %s tests. Significant sets: %d at per-database FDR < %s, %d at pooled FDR (%d lost, %d gained)."),
    TOP_N, n_db_shown, n_db_pooled,
    format(n_total, big.mark = ","), n_before, format(FDR), n_after, n_lost, n_gained
  )

  p <- ggplot(plot_data, aes(x = nes, y = row_key)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    geom_segment(aes(x = 0, xend = nes, y = row_key, yend = row_key),
                 colour = "grey70", linewidth = 0.6) +
    geom_point(aes(fill = nes, shape = sig_status), size = 3.5, stroke = 0.9, colour = "grey20") +
    geom_point(data = key_stub, aes(x = nes, y = row_key, shape = sig_status),
               inherit.aes = FALSE, alpha = 0, size = 3.5, stroke = 0.9,
               colour = "grey20", show.legend = TRUE) +
    scale_fill_gradient2(
      low = NEG, mid = MID, high = POS, midpoint = 0,
      limits = c(-1, 1) * as.numeric(FIG_CFG$figures$nes_cap %||% 3.5),
      oob = scales::squish, name = "NES"
    ) +
    scale_shape_manual(values = SIG_SHAPES,
                       name = sprintf("Set status when the per-database FDR is replaced by one pooled FDR (both at %s)",
                                      format(FDR)),
                       drop = FALSE) +
    scale_y_discrete(labels = function(x)
      stringr::str_wrap(sub(paste0("^[0-9]+", KEY_SEP), "", x), width = LABEL_WRAP)) +
    guides(
      # Titles ride ABOVE their keys: with both titles inline the guide row was
      # wider than the canvas and the trailing guide was clipped by the device.
      shape = guide_legend(
        order = 1, nrow = 1, title.position = "top",
        override.aes = list(alpha = 1, size = 4, stroke = 0.9,
                            colour = "grey20", fill = SIG_KEY_FILL)),
      fill = guide_colourbar(
        order = 2, title.position = "top",
        theme = theme(legend.key.width  = unit(3, "cm"),
                      legend.key.height = unit(0.5, "cm")))
    ) +
    facet_wrap(~ db_label, scales = "free_y", ncol = FACET_NCOL) +
    labs(
      title    = sprintf("GSEA overview across %d databases: %s", n_db_shown, co_short),
      subtitle = sub_text,
      x        = "Normalized Enrichment Score (NES)",
      y        = NULL,
      caption  = FACE_NOTE
    ) +
    project_theme(config = FIG_CFG, base_size = FIG_CFG$figures$base_size %||% 14) +
    theme(
      axis.text.y                 = element_text(lineheight = 0.85),
      # Anchor the caption and the guide box to the CANVAS, not to the panel
      # region: panel-anchored prose starts after the y-axis labels and runs
      # off the right edge on a two-column facet grid.
      plot.caption.position       = "plot",
      legend.position             = "bottom",
      legend.box                  = "vertical",
      legend.box.just             = "left",
      legend.justification.bottom = "left",
      panel.spacing.x             = unit(1.2, "lines")
    )

  attr(p, "n_db_shown")  <- n_db_shown
  attr(p, "n_db_pooled") <- n_db_pooled
  p
}

# ---------------------------------------------------------------------------
# Caption text, computed from the same numbers the panel draws
# ---------------------------------------------------------------------------
contrast_finding <- function(df_sub, contrast_name, n_db_shown, n_db_pooled) {
  co_long  <- contrast_label(contrast_name)
  n_total  <- df_sub$n_tests_pooled[1]
  n_before <- sum(df_sub$padj < FDR, na.rm = TRUE)
  n_after  <- sum(df_sub$padj_pooled < FDR, na.rm = TRUE)
  n_lost   <- sum(df_sub$padj < FDR & df_sub$padj_pooled >= FDR, na.rm = TRUE)
  n_gained <- sum(df_sub$padj >= FDR & df_sub$padj_pooled < FDR, na.rm = TRUE)

  db_null <- df_sub %>%
    dplyr::filter(database %in% DB_SHOWN) %>%
    dplyr::distinct(database, sig_after_in_db) %>%
    dplyr::filter(sig_after_in_db == 0) %>%
    dplyr::pull(database)
  null_txt <- if (length(db_null) == 1L) {
    sprintf("%s carries no pooled-significant set at all and keeps its facet so the null stays visible.",
            as.character(db_null))
  } else if (length(db_null)) {
    sprintf("%d of them carry no pooled-significant set at all (%s) and keep their facets so the nulls stay visible.",
            length(db_null), paste(sort(as.character(db_null)), collapse = ", "))
  } else {
    "Every displayed database carries at least one pooled-significant set."
  }

  sprintf(paste0(
    "Cross-database GSEA overview for the %s contrast, labelled %s, under the pooled Benjamini-Hochberg family ",
    "of %s tests across %d databases. Per-database FDR < %s gives %d significant sets; the pooled FDR gives ",
    "%d (%d lost on pooling, %d gained). The panel draws the top %d sets per database for %d of the %d pooled ",
    "databases. %s"),
    contrast_name, co_long, format(n_total, big.mark = ","), n_db_pooled, format(FDR),
    n_before, n_after, n_lost, n_gained, TOP_N, n_db_shown, n_db_pooled, null_txt)
}

HOW_TO_READ_COMMON <- paste0(
  "Lollipop panel, one facet per database, x = Normalized Enrichment Score (NES). Orange fill = positive NES ",
  "(enriched in the contrast numerator), blue = negative, squished at the configured NES cap. ",
  "GLYPH = what becomes of a set when the per-database Benjamini-Hochberg correction is replaced by one pooled ",
  "correction over every database: filled circle = significant under both, diamond = per-database only (lost on ",
  "pooling), triangle = pooled only (gained on pooling), open circle = neither. The legend keys carry a neutral ",
  "grey fill because the panel spends fill on NES. ",
  "SELECTION RULE: top N per database by padj_pooled, ties broken by the raw p-value then by |NES|; N is ",
  "gsea_pooled_overview.top_n_per_db. Nothing is pinned on top of that rank. An earlier revision force-included ",
  "the Hallmark HYPOXIA and INTERFERON sets and the three project-curated lenses whatever their rank, which let ",
  "the panel argue for a conclusion; that pin is gone. ",
  "FACET HEADER (n/N): n of the N sets in that database pass the pooled FDR, counted over the whole database. ",
  "A header reading (0/23) marks a real null, and those facets stay on the canvas, since dropping a database ",
  "because it came out empty would leave a panel selected for positives. ",
  "EXCLUDED FROM THE PANEL: whatever gsea_pooled_overview.exclude_databases names, which is ",
  "currently empty, so the facet count and the pooled-family database count agree. An excluded ",
  "database stays inside the pooled correction and in every 06_gsea table. ",
  "CAVEAT: padj_pooled is a comparability device across databases rather than a calibrated error rate, since GO ",
  "terms and pathway sets share genes. ",
  "Claim tier: L3 (enrichment statistics). Sample mapping owner-confirmed."
)

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------
CONTRASTS <- c("WT_heat", "Interaction")
panels <- list()

for (co in CONTRASTS) {
  df_co <- df_all %>% dplyr::filter(contrast == co)
  p_co  <- make_contrast_overview_plot(df_co, co)
  panels[[co]] <- p_co

  n_shown  <- attr(p_co, "n_db_shown")
  n_pooled <- attr(p_co, "n_db_pooled")
  .title_fits(p_co, CANVAS_W)

  save_overview(
    plot        = p_co,
    stage       = STAGE,
    name        = sprintf("gsea_pooled_overview_%s", co),
    table       = df_co,
    finding     = contrast_finding(df_co, co, n_shown, n_pooled),
    script      = SCRIPT,
    fn          = "ggplot2 / geom_point + facet_wrap",
    config_kv   = sprintf("thresholds.gsea_fdr=%s; gsea_pooled_overview.top_n_per_db=%d; gsea_pooled_overview.exclude_databases=[%s]; gsea_pooled_overview.facet_ncol=%d; figures.nes_cap=%s",
                          format(FDR), TOP_N, paste(DB_EXCLUDE, collapse = ","), FACET_NCOL,
                          format(FIG_CFG$figures$nes_cap %||% 3.5)),
    input       = "03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv",
    how_to_read = HOW_TO_READ_COMMON,
    config      = FIG_CFG,
    wide        = TRUE,
    height      = CHROME_H + .panel_block_h(n_shown)
  )
}

# Combined two-panel figure: the same two panels, same selection rule, SIDE BY SIDE.
# Side by side rather than stacked because each panel is already ~15 in tall at N sets
# per database across a 7-row facet grid; stacking them produced a canvas over 30 in
# tall that no reader can take in, while placing them abreast keeps a landscape aspect
# and lines the database blocks up row for row across the seam.
n_shown_all <- attr(panels[["WT_heat"]], "n_db_shown")
# The per-panel legend and caption are identical, but patchwork cannot merge them
# (each panel's key layer carries its own data, so the guides do not hash equal).
# Drop both from the right panel and hoist the note to the figure level instead.
p_combined <- ((panels[["WT_heat"]] + labs(caption = NULL)) |
               (panels[["Interaction"]] + labs(caption = NULL) + theme(legend.position = "none"))) +
  plot_annotation(
    title   = sprintf("Pooled-FDR GSEA overview across %d databases: WT heat beside Interaction",
                      n_shown_all),
    caption = FACE_NOTE,
    theme   = project_theme(config = FIG_CFG, base_size = FIG_CFG$figures$base_size %||% 14)
  )

save_overview(
  plot        = p_combined,
  stage       = STAGE,
  name        = "gsea_pooled_overview",
  table       = df_all,
  finding     = sprintf(paste0(
    "The two per-contrast pooled-FDR overviews abreast on one canvas, WT heat on the left and Interaction on the ",
    "right, so the %d displayed databases line up block for block between the contrasts. Under the shared pooled ",
    "family of %s tests across %d databases, WT heat carries %d significant sets before pooling and %d after, ",
    "Interaction %d and %d. The per-contrast panels hold the full numbers."),
    n_shown_all, format(df_all$n_tests_pooled[1], big.mark = ","),
    attr(panels[["WT_heat"]], "n_db_pooled"),
    sum(df_all$contrast == "WT_heat"     & df_all$padj        < FDR, na.rm = TRUE),
    sum(df_all$contrast == "WT_heat"     & df_all$padj_pooled < FDR, na.rm = TRUE),
    sum(df_all$contrast == "Interaction" & df_all$padj        < FDR, na.rm = TRUE),
    sum(df_all$contrast == "Interaction" & df_all$padj_pooled < FDR, na.rm = TRUE)),
  script      = SCRIPT,
  fn          = "patchwork / p_wt beside p_int",
  config_kv   = sprintf("thresholds.gsea_fdr=%s; gsea_pooled_overview.top_n_per_db=%d; gsea_pooled_overview.exclude_databases=[%s]",
                        format(FDR), TOP_N, paste(DB_EXCLUDE, collapse = ",")),
  input       = "03_results/06_gsea/tables/_overview/gsea_pooled_overview.csv",
  how_to_read = paste0(
    "Two copies of the per-contrast panel abreast: WT heat on the left, Interaction on the right, sharing the ",
    "single legend at the foot. The facet order is identical on both sides, so a database sits at the same height ",
    "in each; each panel's own subtitle carries its pooled-family counts, which differ between contrasts. ",
    HOW_TO_READ_COMMON,
    " This canvas is double width by construction (two full facet grids at N sets per database); for reading one ",
    "contrast, prefer gsea_pooled_overview_WT_heat.png or gsea_pooled_overview_Interaction.png."
  ),
  config      = FIG_CFG,
  width       = 2 * CANVAS_W,
  height      = CHROME_H + 0.8 + .panel_block_h(n_shown_all)
)

cat("[12b] All GSEA pooled overview figures and tables generated successfully.\n")
