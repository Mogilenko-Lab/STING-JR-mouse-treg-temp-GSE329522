## 02_analysis/helpers/figure_style_unified.R — Wave-0 UNIFIED figure-style prototype.
## =====================================================================================
## Proves the new single-variant, dual-FORMAT contract IN-PROJECT before it is promoted
## up into SciAgent-toolkit. Sources the existing shim (toolkit path/caption/table helpers
## + FIG_CFG + contrast_label) and then SHADOWS the styling + export functions:
##
##   * project_theme()  — ONE legible theme (no print/screen variant), palette-aware.
##   * save_figure()    — ONE themed plot -> <name>.pdf AND <name>.png (same geometry,
##                        no .print/.screen suffix). cairo_pdf for Unicode glyphs.
##   * save_overview()  — same, plus sibling table + README caption (caption -> <name>.png).
##   * style_series()   — patchwork-CORRECT running-sum normalizer (uses `&`, not `+`),
##                        fixed ES y-range + inside legend; the running-sum headline fix.
##
## SOURCED BY figure_style.R (at its end), AFTER the toolkit lib + FIG_CFG + contrast_label
## are already in scope. Viz scripts source figure_style.R (the project entry point); they do
## NOT source this file directly. These definitions SHADOW the toolkit lib's project_theme /
## save_figure / save_overview / style_series with the unified single-variant, dual-format style.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## ---------------------------------------------------------------------------
## 1. THEME — single unified tier (compact-but-legible: ~14pt @ 8.5x6.5in).
##    No `variant`. Sizes from FIG_CFG$figures; legible shrunk to a column AND projected.
## ---------------------------------------------------------------------------
project_theme <- function(base_size = NULL, legend = TRUE, config = NULL, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("project_theme() needs ggplot2.")
  cfg <- config %||% FIG_CFG
  f   <- cfg$figures %||% list()
  bs  <- as.numeric(base_size %||% f$base_size %||% 14)
  ggplot2::theme_minimal(base_size = bs) +
    ggplot2::theme(
      text             = ggplot2::element_text(size = bs),
      plot.title       = ggplot2::element_text(size = f$title_size    %||% 16, face = "bold"),
      plot.subtitle    = ggplot2::element_text(size = f$subtitle_size %||% 11, colour = "grey25"),
      plot.caption     = ggplot2::element_text(size = f$caption_size   %||% 9,  colour = "grey45",
                                               hjust = 0, lineheight = 1.05),
      axis.title       = ggplot2::element_text(size = f$axis_title_size %||% 13),  # plain (not bold)
      axis.text        = ggplot2::element_text(size = f$axis_text_size  %||% 11),
      legend.text      = ggplot2::element_text(size = f$legend_text_size %||% 11),
      legend.title     = ggplot2::element_text(size = f$legend_text_size %||% 11, face = "bold"),
      strip.text       = ggplot2::element_text(size = f$strip_size %||% 12, face = "bold"),
      legend.key.spacing.y = ggplot2::unit(3, "pt"),     # a little air between legend rows
      legend.position  = if (isTRUE(legend)) "right" else "none",
      panel.grid.minor = ggplot2::element_blank(),
      axis.line        = ggplot2::element_line(linewidth = 0.4),
      panel.border     = ggplot2::element_blank(),
      plot.title.position = "plot",
      plot.margin      = ggplot2::margin(8, 12, 8, 8))
}
set_paper_style <- function(...) project_theme(...)

## ---------------------------------------------------------------------------
## 2. EXPORT — ONE plot -> <name>.pdf + <name>.png, ONE geometry, ONE theme.
##    `variant` accepted but IGNORED (drop-in compat with existing call sites).
##    `name` may carry a subdir (e.g. "Hallmark/dotplot"); the subdir is created.
## ---------------------------------------------------------------------------
.fig_geom <- function(config, width, height, wide) {
  f <- (config %||% FIG_CFG)$figures %||% list()
  w <- as.numeric(width  %||% (if (isTRUE(wide)) (f$width_wide %||% 13) else (f$width %||% 8.5)))
  h <- as.numeric(height %||% (f$height %||% 6.5))
  c(w, h)
}

.purge_stem <- function(out_dir, stem) {
  if (!dir.exists(out_dir)) return(invisible(0L))
  all <- list.files(out_dir, full.names = FALSE)
  hit <- startsWith(all, paste0(stem, ".")) &
         (endsWith(all, ".png") | endsWith(all, ".pdf"))
  if (any(hit)) file.remove(file.path(out_dir, all[hit]))
  invisible(sum(hit))
}

save_figure <- function(plot, stage, name, variant = NULL, contrast = NULL,
                        overview = FALSE, config = NULL,
                        width = NULL, height = NULL, wide = FALSE, void = FALSE) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("save_figure() needs ggplot2.")
  cfg <- config %||% FIG_CFG
  f   <- cfg$figures %||% list()
  geo <- .fig_geom(cfg, width, height, wide)

  base_dir <- if (!is.null(contrast)) contrast_path(stage, contrast, "figures", cfg)
              else if (isTRUE(overview)) overview_path(stage, "figures", cfg)
              else { d <- file.path(cfg$paths$results %||% "03_results/", stage, "figures")
                     dir.create(d, recursive = TRUE, showWarnings = FALSE); d }
  sub  <- dirname(name); stem <- basename(name)
  out_dir <- if (identical(sub, ".")) base_dir else file.path(base_dir, sub)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  .purge_stem(out_dir, stem)

  ## Theming is entirely the CALLER's job (every viz call site already adds
  ## project_theme()/style_series()). save_figure must NOT re-theme: project_theme is a
  ## COMPLETE theme (theme_minimal base), so a second application RESETS every element the
  ## caller set AFTER their own project_theme() — plot.tag.position, legend.position="bottom",
  ## axis.text.x angle, style_series's inside-legend pin — silently clobbering per-figure
  ## tweaks. This matches the reference contract (callers theme; the exporter only writes).
  styled <- plot
  if (isTRUE(void) && exists(".void_overlay")) styled <- styled + .void_overlay()

  written <- list()
  for (ext in c("pdf", "png")) {
    out <- file.path(out_dir, paste0(stem, ".", ext))
    dev <- if (identical(ext, "pdf") && isTRUE(capabilities()[["cairo"]]))
             grDevices::cairo_pdf else NULL
    ggplot2::ggsave(out, styled, width = geo[1], height = geo[2],
                    dpi = as.numeric(f$dpi %||% 300), device = dev)
    written[[ext]] <- out
  }
  message(sprintf("  [unified] save_figure: %s.{pdf,png} (%gx%gin)", basename(name), geo[1], geo[2]))
  invisible(written)
}

## ---------------------------------------------------------------------------
## 3. OVERVIEW — figure (both formats) + sibling table + README caption, atomic.
##    Caption fig path now points to <name>.png (no .screen suffix).
## ---------------------------------------------------------------------------
save_overview <- function(plot, stage, name, table, finding, script, fn, config_kv, input,
                          how_to_read, contrast = NULL, config = NULL,
                          width = NULL, height = NULL, wide = FALSE, void = FALSE) {
  cfg <- config %||% FIG_CFG
  overview <- is.null(contrast)
  figs <- save_figure(plot, stage, name, contrast = contrast, overview = overview,
                      config = cfg, width = width, height = height, wide = wide, void = void)
  bcd <- cfg$figures$by_contrast_dir %||% "by_contrast"
  ovd <- cfg$figures$overview_dir    %||% "_overview"
  if (!is.null(contrast)) {
    tdir    <- contrast_path(stage, contrast, "tables", cfg)
    fig_rel <- file.path("figures", bcd, contrast, sprintf("%s.png", name))
  } else {
    tdir    <- overview_path(stage, "tables", cfg)
    fig_rel <- file.path("figures", ovd, sprintf("%s.png", name))
  }
  table_path <- file.path(tdir, sprintf("%s.csv", basename(name)))
  if (!is.null(table)) utils::write.csv(round_numeric_cols(table), table_path, row.names = FALSE)
  readme <- write_caption(stage, fig_rel, finding = finding, script = script, fn = fn,
                          config_kv = config_kv, input = input, how_to_read = how_to_read,
                          config = cfg)
  invisible(list(figures = figs, table = table_path, readme = readme))
}

## ---------------------------------------------------------------------------
## 4. RUNNING-SUM / SERIES normalizer — now a THIN PASS-THROUGH to the toolkit.
##    The toolkit gsea_running_sum_plot() owns panel construction AND alignment;
##    it returns a patchwork carrying a `grs_restyle` closure (the clean extension
##    interface). We re-skin the running-sum via NAMED knobs — ES y clamped to
##    `ylim`, a SINGLE legend collected OUTSIDE on the right, x ticks ONLY on the
##    bottom panel, hidden rug y-index labels, project panel_heights, and our
##    project_theme as the base — WITHOUT ever indexing styled[[i]] (the old
##    desync-prone path). The composer applies project_theme first, then re-asserts
##    the per-panel chrome, so alignment + tick/label suppression hold by construction.
##    Per-panel aspect.ratio is never touched (it would break width alignment).
## ---------------------------------------------------------------------------
style_series <- function(plot, ylim = NULL, config = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("style_series() needs ggplot2.")
  cfg  <- config %||% FIG_CFG
  ylim <- as.numeric(unlist(ylim %||% (cfg$figures$running_sum_ylim %||% c(-1, 1))))
  stopifnot(length(ylim) == 2, all(is.finite(ylim)))
  ph <- as.numeric(unlist(cfg$figures$running_sum_heights %||% c(2.4, 0.7, 0.9)))

  ## Clean toolkit interface: re-skin via the attached composer closure (no indexing).
  restyle <- attr(plot, "grs_restyle")
  if (is.function(restyle)) {
    styled <- restyle(
      es_ylim         = ylim,                       # clamp ES y for comparability
      legend_position = "right",                    # ONE legend, outside-right
      xticks          = "bottom",                   # x ticks only on the bottom panel
      rug_ylabels     = FALSE,                       # hide rug y-index labels
      panel_heights   = ph,                          # ES : rug : metric proportions
      base_theme      = project_theme(config = cfg)) # project base; chrome re-asserted on top
    ## Top-align the collected outside-right legend so it sits at the level of the
    ## top (ES) panel rather than vertically centred across all three panels.
    return(styled & ggplot2::theme(
      legend.justification.right = "top",
      legend.justification       = "top"))
  }

  if (!inherits(plot, "patchwork")) {
    ## non-patchwork series figure: simple shared-y + inside legend
    styled <- tryCatch(plot + project_theme(config = cfg), error = function(e) plot)
    return(styled + ggplot2::coord_cartesian(ylim = ylim) +
             ggplot2::theme(legend.position = "inside",
                            legend.position.inside = c(0.98, 0.98),
                            legend.justification = c(1, 1),
                            legend.background = ggplot2::element_rect(fill = "white", colour = "grey90")))
  }

  ## Fallback: a patchwork WITHOUT the toolkit closure (e.g. a non-toolkit series
  ## figure). Theme + collect a single right legend at the figure level via `&`
  ## (no panel indexing). We deliberately do NOT clamp y here — a global `&` ylim
  ## would wrongly squash the rug/metric panels; ES clamping is the toolkit's job.
  styled <- tryCatch(plot & project_theme(config = cfg), error = function(e) plot)
  styled <- tryCatch(styled + patchwork::plot_layout(heights = ph, guides = "collect"),
                     error = function(e) styled)
  styled & ggplot2::theme(
    legend.position      = "right",
    legend.key.spacing.y = ggplot2::unit(3, "pt"),
    legend.background    = ggplot2::element_rect(fill = "white", colour = "grey90"),
    legend.key.size      = ggplot2::unit(0.8, "lines"))
}
## Canonical name for the same operation (so a viz script can call either).
style_running_sum <- function(p, ylim = NULL, config = NULL) style_series(p, ylim = ylim, config = config)

## ---------------------------------------------------------------------------
## 5. PALETTE — colorblind-safe scales sourced from FIG_CFG$colors (semantic categories).
## ---------------------------------------------------------------------------
.okabe <- function(config = NULL) {
  oi <- (config %||% FIG_CFG)$colors$okabe_ito %||% list()
  unname(unlist(oi))
}
scale_color_okabe <- function(..., config = NULL)
  ggplot2::scale_color_manual(values = .okabe(config), ...)
scale_fill_okabe  <- function(..., config = NULL)
  ggplot2::scale_fill_manual(values = .okabe(config), ...)

invisible(TRUE)
