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
