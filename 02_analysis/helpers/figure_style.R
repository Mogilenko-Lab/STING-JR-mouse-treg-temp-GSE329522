## 02_analysis/helpers/figure_style.R — per-project figure-style shim.
## =====================================================================
## ONE import per viz script. Delegates to the SciAgent-toolkit contract
## lib (02_analysis/helpers/figure-style/figure_helpers.R, symlinked by
## `sciagent activate`) and loads the project analysis_config.yaml.
##
## Usage in any viz script (run from the project root):
##   source("02_analysis/helpers/figure_style.R")
##   p <- ggplot(...) + project_theme(config = FIG_CFG)
##   save_overview(p, "06_gsea", "name", table = df, ..., config = FIG_CFG)
##
## Provides (from the toolkit lib): load_figure_config, project_theme/
## set_paper_style, save_figure (print+screen variants), save_overview
## (atomic figure+table+caption), contrast_path/overview_path, style_series,
## purge_figures, write_caption, append_master_table, round_numeric_cols,
## direction_cue. FIG_CFG is the project-wide config handle.

# Portable null-coalesce, defined BEFORE first use.
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

## ---------------------------------------------------------------------------
## 1. Source the symlinked contract lib (graceful fallback if not linked).
##    Scripts are run from the project root, so the relative path is stable;
##    `sciagent activate` maintains the figure-style symlink.
## ---------------------------------------------------------------------------
.HELPERS_LIB <- "02_analysis/helpers/figure-style/figure_helpers.R"

.FIGURE_STYLE_LOADED <- FALSE
if (file.exists(.HELPERS_LIB)) {
    source(.HELPERS_LIB)
    .FIGURE_STYLE_LOADED <- TRUE
} else {
    warning(
        "[figure_style] toolkit lib not found at: ", .HELPERS_LIB,
        "\nRun `sciagent activate` in this repo to link lib/figure-style/. ",
        "Falling back to minimal stubs.", call. = FALSE
    )
    # --- MINIMAL FALLBACK stubs (keep scripts from hard-crashing) -----------
    load_figure_config <- function(path = "02_analysis/config/analysis_config.yaml") {
        if (requireNamespace("yaml", quietly = TRUE)) yaml::read_yaml(path)
        else { warning("[figure_style] yaml unavailable; using empty config."); list() }
    }
    project_theme <- function(base_size = NULL, legend = TRUE, variant = "screen",
                               config = NULL) {
        if (!requireNamespace("ggplot2", quietly = TRUE))
            stop("[figure_style] ggplot2 required for project_theme().")
        ggplot2::theme_minimal(base_size = base_size %||% 16)
    }
    set_paper_style <- function(base_size = NULL, legend = TRUE, variant = "screen",
                                 config = NULL) {
        project_theme(base_size = base_size, legend = legend, variant = variant,
                      config = config)
    }
}

## ---------------------------------------------------------------------------
## 2. Load the project config once (FIG_CFG is the stable project-wide handle).
## ---------------------------------------------------------------------------
FIG_CFG <- tryCatch(
    load_figure_config("02_analysis/config/analysis_config.yaml"),
    error = function(e) {
        warning("[figure_style] Could not load analysis_config.yaml: ", conditionMessage(e),
                call. = FALSE)
        list()
    }
)

## ---------------------------------------------------------------------------
## 3. Contrast DISPLAY labels — single source of truth for every viz script.
##    Read from design.contrast_labels{,_short} in the project config so all of
##    12/13/14/16 render the SAME human-readable contrast names (the analytic
##    keys WT_heat/Interaction/... stay the join keys; these are labels only).
##    contrast_label() is VECTORIZED (safe inside dplyr::mutate AND scalar in
##    per-contrast loops); falls back to the key when no label is configured.
## ---------------------------------------------------------------------------
.cl_to_vec <- function(x) {
    if (is.null(x) || length(x) == 0) return(character(0))
    v <- vapply(x, function(s) as.character(s)[1], character(1))
    stats::setNames(v, names(x))
}
CONTRAST_LABELS       <- .cl_to_vec(FIG_CFG$design$contrast_labels)
CONTRAST_LABELS_SHORT <- .cl_to_vec(FIG_CFG$design$contrast_labels_short)

contrast_label <- function(co, short = FALSE) {
    lut <- if (isTRUE(short)) CONTRAST_LABELS_SHORT else CONTRAST_LABELS
    lab <- unname(lut[as.character(co)])
    ifelse(is.na(lab), as.character(co), lab)
}

## ---------------------------------------------------------------------------
## 3b. SAMPLE-MAPPING PROVENANCE — single source of truth for every viz script.
##     Read from design.sample_mapping in the project config. No viz script may
##     hard-code the mapping status: the 2026-07-22 owner sample sheet flipped it
##     INFERRED -> CONFIRMED, and figures regenerated after that date kept printing
##     the old hedge because the string lived in nine scripts instead of one key.
##
##     config.R defines the same four accessors off the same config key, for the
##     compute scripts that source it instead of this shim. Definitions are guarded
##     by exists() so whichever loads first wins; both read design$sample_mapping,
##     so the two cannot disagree on a value. Use sample_mapping_stamp() for a
##     figure canvas (short) and sample_mapping_caption() for a README (names the
##     evidence). Full per-library account: 00_data/processed/PROVENANCE.md.
## ---------------------------------------------------------------------------
.SAMPLE_MAPPING <- FIG_CFG$design$sample_mapping %||% list()

if (!exists("sample_mapping_status", mode = "function")) {
    sample_mapping_status <- function() toupper(.SAMPLE_MAPPING$status %||% "INFERRED")
}
if (!exists("sample_mapping_confirmed", mode = "function")) {
    sample_mapping_confirmed <- function() identical(sample_mapping_status(), "CONFIRMED")
}
if (!exists("sample_mapping_stamp", mode = "function")) {
    sample_mapping_stamp <- function() {
        if (sample_mapping_confirmed()) {
            .SAMPLE_MAPPING$stamp %||% "Sample mapping owner-confirmed"
        } else {
            .SAMPLE_MAPPING$stamp_unconfirmed %||%
                "Sample mapping inferred, pending the owner sample sheet"
        }
    }
}
if (!exists("sample_mapping_caption", mode = "function")) {
    sample_mapping_caption <- function() {
        if (sample_mapping_confirmed()) {
            .SAMPLE_MAPPING$caption %||% sample_mapping_stamp()
        } else {
            .SAMPLE_MAPPING$caption_unconfirmed %||% sample_mapping_stamp()
        }
    }
}
## Deprecated alias; new code should call sample_mapping_stamp()/_caption().
if (!exists("provisional_caption", mode = "function")) {
    provisional_caption <- function() sample_mapping_stamp()
}

## ---------------------------------------------------------------------------
## 3c. CANVAS-FIT GUARD — make an over-wide title/subtitle line a hard error.
##     ggsave() draws a line that is wider than the canvas straight off the edge
##     and still exits 0, so a clipped title is invisible to every text-level
##     check and survives review. On 2026-08-05 six of the seventeen 04_tf
##     overview panels were shipping with their title or subtitle cut off at the
##     right edge for exactly this reason. Measure the drawn line and stop.
##
##     Deliberately OPT-IN — call it from a viz script next to its labs(), do not
##     wire it into save_figure(). A silent global guard would begin failing
##     renders across every stage at once, which is a migration, not a fix.
##
##     Usage:
##       fits_canvas(my_title,    FIG_CFG$figures$title_size,    "bold",  W, "title")
##       fits_canvas(my_subtitle, FIG_CFG$figures$subtitle_size, "plain", W, "subtitle")
##     `width_in` must be the SAME width later passed to save_figure/save_overview
##     (or FIG_CFG$figures$width when the call site does not override it).
## ---------------------------------------------------------------------------
##     MEASURE ON THE OUTPUT DEVICE, NOT pdf(NULL). The base pdf device carries
##     Helvetica metrics; save_figure() renders through cairo_pdf and ragg/cairo
##     png, whose default sans is ~16% wider. Measuring on pdf(NULL) passes
##     strings that then ship clipped -- that mistake was made and caught here on
##     2026-08-05, so the device below is chosen to match the real output.
fits_canvas <- function(text, fontsize, fontface, width_in, what = "text",
                        margin_in = 0.45) {
    if (is.null(text) || !nzchar(as.character(text))) return(invisible(TRUE))
    tmp <- tempfile(fileext = ".png")
    if (requireNamespace("ragg", quietly = TRUE)) {
        ragg::agg_png(tmp, width = width_in, height = 4, units = "in", res = 300)
    } else if (isTRUE(grDevices::capabilities()[["cairo"]])) {
        grDevices::png(tmp, width = width_in, height = 4, units = "in",
                       res = 300, type = "cairo")
    } else {
        grDevices::pdf(NULL, width = width_in, height = 4)
    }
    on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)
    avail <- width_in - margin_in
    for (ln in strsplit(as.character(text), "\n", fixed = TRUE)[[1]]) {
        if (!nzchar(ln)) next
        w <- grid::convertWidth(
            grid::grobWidth(grid::textGrob(
                ln, gp = grid::gpar(fontsize = fontsize, fontface = fontface))),
            "inches", valueOnly = TRUE)
        if (w > avail) {
            stop(sprintf(paste0("%s line is %.2fin wide on a %.2fin canvas ",
                                "(%.2fin usable) and would be CLIPPED: %s"),
                         what, w, width_in, avail, ln), call. = FALSE)
        }
    }
    invisible(TRUE)
}

## ---------------------------------------------------------------------------
## 4. UNIFIED STYLE OVERRIDES — single-variant theme + dual-format (pdf+png) export.
##    Sourced LAST so it SHADOWS the toolkit lib's project_theme / save_figure /
##    save_overview / style_series with the unified contract (no .print/.screen split;
##    compact-but-legible 14pt @ 8.5x6.5in; patchwork-correct running-sum normalizer).
##    Proven in-project (Wave-0, 2026-06-25); promoted to the toolkit in Wave 2.
## ---------------------------------------------------------------------------
.UNIFIED_OVERRIDES <- "02_analysis/helpers/figure_style_unified.R"
if (file.exists(.UNIFIED_OVERRIDES)) {
    source(.UNIFIED_OVERRIDES)
} else {
    warning("[figure_style] unified overrides not found at: ", .UNIFIED_OVERRIDES,
            "; falling back to the toolkit lib's dual-variant style.", call. = FALSE)
}
