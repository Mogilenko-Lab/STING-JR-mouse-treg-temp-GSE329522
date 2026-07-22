# config.R - Configuration loader for STING-cGAS-GSE329522
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2: genotype x temperature)
# Species: Mus musculus
#
# Stage-based layout (NOT the flat .ref layout):
#   03_results/<stage_id>/{figures,tables}/   -- per-stage deliverables
#   03_results/objects/                       -- reusable state (.rds checkpoints)
#   03_results/master/                        -- cross-stage accumulator tables
#   03_results/interactive/                   -- standalone HTML dashboards
#   03_results/_scratch/                      -- throwaway artifacts
#
# Sourcing this file is intentionally side-effect-light: it reads the YAML and
# defines constants/helpers but does NOT library() heavy packages or create
# directories. Scripts call load_packages() and stage_dir() explicitly.

# ============================================================================
# NULL-COALESCE
# ============================================================================

if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
}

# ============================================================================
# PROJECT ROOT + YAML
# ============================================================================

PROJECT_ROOT <- getwd()

if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Package 'yaml' is required but not installed.")
}

.config_path <- file.path(PROJECT_ROOT, "02_analysis/config/analysis_config.yaml")
if (!file.exists(.config_path)) {
  # Fallback: relative to a sourced script location
  if (!is.null(sys.frame(1)$ofile)) {
    .config_path <- file.path(dirname(sys.frame(1)$ofile), "analysis_config.yaml")
  }
}
if (!file.exists(.config_path)) {
  stop(sprintf("analysis_config.yaml not found (looked at %s). Run from project root.", .config_path))
}

YAML_CONFIG <- yaml::read_yaml(.config_path)

# ============================================================================
# PROJECT CONSTANTS
# ============================================================================

PROJECT_ID <- YAML_CONFIG$project$id
SPECIES    <- YAML_CONFIG$project$species

# ============================================================================
# DIRECTORIES (stage-based layout; read roots from YAML$paths)
# ============================================================================

DIR_RESULTS     <- file.path(PROJECT_ROOT, sub("/+$", "", YAML_CONFIG$paths$results %||% "03_results"))
DIR_OBJECTS     <- file.path(PROJECT_ROOT, sub("/+$", "", YAML_CONFIG$paths$objects %||% "03_results/objects"))
DIR_MASTER      <- file.path(PROJECT_ROOT, sub("/+$", "", YAML_CONFIG$paths$master %||% "03_results/master"))
DIR_INTERACTIVE <- file.path(PROJECT_ROOT, sub("/+$", "", YAML_CONFIG$paths$interactive %||% "03_results/interactive"))
DIR_SCRATCH     <- file.path(PROJECT_ROOT, sub("/+$", "", YAML_CONFIG$paths$scratch %||% "03_results/_scratch"))

# Per-stage subdir names (figures/tables) from YAML, with sane defaults.
.STAGE_FIGURES_SUBDIR <- YAML_CONFIG$paths$stage_figures_subdir %||% "figures"
.STAGE_TABLES_SUBDIR  <- YAML_CONFIG$paths$stage_tables_subdir  %||% "tables"

#' Resolve (and create) a per-stage output directory.
#'
#' @param id    Stage id, e.g. "03_de". Must be registered in YAML_CONFIG$stages.
#' @param kind  "figures" (default) or "tables".
#' @return Absolute path to 03_results/<id>/<figures|tables>, created on call.
stage_dir <- function(id, kind = c("figures", "tables")) {
  kind <- match.arg(kind)
  registered <- vapply(YAML_CONFIG$stages, function(s) s$id, character(1))
  if (!id %in% registered) {
    warning(sprintf("stage_dir(): stage '%s' is NOT registered in analysis_config.yaml (stages: %s). Register it before writing artifacts.",
                    id, paste(registered, collapse = ", ")))
  }
  sub <- if (kind == "figures") .STAGE_FIGURES_SUBDIR else .STAGE_TABLES_SUBDIR
  path <- file.path(DIR_RESULTS, id, sub)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

# ============================================================================
# INPUT PATHS
# ============================================================================

FILE_CPM <- file.path(PROJECT_ROOT, "00_data/processed/GSE329522_normalized_counts_CPM_iTreg.csv")

# ============================================================================
# GSEA / RANKING PARAMETERS (from YAML thresholds; do NOT redefine in scripts)
# ============================================================================

GSEA_NPERM      <- YAML_CONFIG$thresholds$gsea_nperm    %||% 100000
GSEA_SEED       <- YAML_CONFIG$thresholds$gsea_seed     %||% 123
GSEA_MIN_SIZE   <- YAML_CONFIG$thresholds$gsea_min_size %||% 15
GSEA_MAX_SIZE   <- YAML_CONFIG$thresholds$gsea_max_size %||% 500
GSEA_FDR_CUTOFF <- YAML_CONFIG$thresholds$gsea_fdr      %||% 0.05
RANK_METRIC     <- "t"   # t-statistic from limma; NEVER logFC

# DE thresholds
DE_FDR   <- YAML_CONFIG$thresholds$de_fdr   %||% 0.05
DE_LOGFC <- YAML_CONFIG$thresholds$de_logfc %||% 1.0

# ============================================================================
# COLORS
# ============================================================================

# Diverging palette named for the reference 3.3_tf_viz.R expectations
# ($negative / $neutral / $positive). YAML stores them as down/neutral/up.
DIVERGING_COLORS <- list(
  negative = YAML_CONFIG$colors$diverging$down    %||% "#2166AC",  # blue (down)
  neutral  = YAML_CONFIG$colors$diverging$neutral %||% "#F7F7F7",  # white
  positive = YAML_CONFIG$colors$diverging$up      %||% "#B35806"   # orange (up)
)

# Group colors from YAML if present (named vector).
GROUP_COLORS <- if (!is.null(YAML_CONFIG$colors$groups)) {
  unlist(YAML_CONFIG$colors$groups)
} else {
  NULL
}

# TF axis/family color palette -- single source of truth for ALL lollipop figures.
#   HIF axis  (Hif1a, Epas1)                        -> orange (UP diverging color; "flat/NS" arm)
#   IFN/NFkB  (Irf1, Irf3, Irf7, Stat1, Stat2,
#               Nfkb1, Rela)                         -> blue   (interaction-positive arm)
#   other     (everything else)                      -> grey
# Used by: 03_decoupler_tf_viz.R (fig3a, fig3b) and
#          03b_decoupler_method_comparison_viz.R (fig3j).
AXIS_COLORS <- c(
  "HIF"   = "#B35806",   # orange -- matches DIVERGING_COLORS$positive
  "IFN"   = "#2166AC",   # blue   -- matches DIVERGING_COLORS$negative
  "other" = "grey55"
)

# Human-readable legend labels for each axis level (reused across figures).
AXIS_LABELS <- c(
  "HIF"   = "HIF axis (Hif1a/Epas1)",
  "IFN"   = "IFN/NFkB axis",
  "other" = "other"
)

# ----------------------------------------------------------------------------
# MODULE-BUCKET PALETTE -- the attribution arc (fig3g-expanded, fig3l, fig3m)
# Single source of truth for the biological re-bucketing of Hif1a's regulon.
# HUE-SEPARATED (audit item 5): the central diagnostic contrast is "repressed
# hypoxic core (DOWN, diagnostic)" vs "shared/glycolytic (UP)". The earlier
# three-shades-of-orange palette made exactly that contrast the LEAST separable,
# so the repressed core now gets a distinct cool hue (dark teal) that is clearly
# unlike the orange UP-set and the purple heat-shock contamination -- and is NOT
# the IFN navy (#2166AC), so it cannot be confused with the IFN axis. Heat-shock
# stays purple (out-of-family contamination).
MODULE_COLORS <- c(
  "heatshock_stress"     = "#762A83",  # purple -- the contamination (Hspa1a/Timp1/Sdc1/Cdkn1a/Serpine1/Eno2)
  "shared_angio_glucose" = "#E08214",  # orange -- Vegfa/Slc2a1 (shared HIF1/HIF2 angio-glucose; UP, least diagnostic)
  "autoreg_feedback"     = "#FDB863",  # pale orange -- Egln3/PHD3 (HIF-induced feedback brake; UP)
  "hif1a_hypoxic_core"   = "#01665E",  # dark teal -- Pdk1/Bnip3/Bnip3l/Car9 (diagnostic core; REPRESSED/DOWN)
  "other_unclassified"   = "grey80"    # faint context -- the ~333 uncurated "other" members
)

# fig3n heat-MAIN regulator axes: HIF/IFN reuse AXIS_COLORS; the heat-shock
# axis (Hsf1) reuses the module-bucket stress purple so HSF1 connects visually
# back to the regulon contaminant.
HEAT_AXIS_COLORS <- c(
  AXIS_COLORS,
  "heatshock" = "#762A83"
)

# ============================================================================
# ANALYSIS CONSTANTS (gene vectors consumed by downstream phases)
# ============================================================================

# Interferon-stimulated genes (cGAS-DEPENDENT arm).
ISG_MARKERS <- c("Ifit1", "Isg15", "Irf7", "Oasl2", "Mx1", "Stat1", "Cxcl10")

# HIF / glycolysis arm (shows NO DETECTABLE cGAS-dependence at n=5 -- an
# asymmetry vs the cGAS-dependent IFN arm, NOT proven independence).
#   HIF-specific:      Slc2a1, Vegfa, Egln3, Bnip3
#   shared-glycolytic: Pgk1,  Ldha,  Aldoa, Hk2   (carried by many programs, not HIF-exclusive)
HIF_GLYCO_MARKERS <- c(
  "Slc2a1", "Vegfa", "Egln3", "Bnip3",   # HIF-specific
  "Pgk1",   "Ldha",  "Aldoa", "Hk2"      # shared-glycolytic
)

# Heat-shock thermometer markers (used for label-blind temperature inference).
THERMO_MARKERS <- c("Hspa1b", "Hsph1", "Hspa1a", "Dnajb1")

# cGAS gene + its alias (Mb21d1 is the legacy symbol for Cgas).
CGAS_GENE <- c("Cgas", "Mb21d1")

# Key transcription factors for the TF-activity forensics.
KEY_TFS <- c("Irf3", "Irf7", "Irf1", "Stat1", "Stat2", "Hif1a", "Epas1", "Nfkb1", "Rela")

# Key PROGENy pathways for the orthogonal pathway split.
KEY_PROGENY <- c("Hypoxia", "JAK-STAT", "NFkB", "TNFa")

# ============================================================================
# SAMPLE-MAPPING CAPTION (single source of truth; every figure uses this)
# ============================================================================

#' Standard sample-provenance stamp for every figure caption/subtitle.
#' The 12630-RS-021..040 -> genotype x temperature mapping is owner-confirmed
#' (2026-07-22 sample sheet; see 00_data/processed/PROVENANCE.md). It matched the
#' Hspa1b/Hsph1 thermometer + Cgas inference exactly (20/20).
#' @return Single character string.
provisional_caption <- function() {
  "Sample mapping owner-confirmed (GSE329522 iTreg 2x2 genotype x temperature)"
}

# ============================================================================
# CHECKPOINT CACHING  (writes into DIR_OBJECTS, NOT a flat checkpoints dir)
# ============================================================================

#' Load a cached checkpoint or compute and save it.
#'
#' @param checkpoint_file  Filename (basename) under DIR_OBJECTS, e.g. "02_de_results.rds".
#' @param compute_fn       Zero-arg function producing the result.
#' @param force            If TRUE, recompute even if the cache exists.
#' @param desc             Human-readable description for log messages.
#' @return The cached or freshly computed result.
load_or_compute <- function(checkpoint_file, compute_fn, force = FALSE, desc = "Result") {
  dir.create(DIR_OBJECTS, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(DIR_OBJECTS, checkpoint_file)
  if (file.exists(path) && !force) {
    message(sprintf("[CACHE] Loading %s ...", desc))
    return(readRDS(path))
  }
  message(sprintf("[COMPUTE] Computing %s ...", desc))
  res <- compute_fn()
  saveRDS(res, path)
  message(sprintf("[SAVE] Saved %s", path))
  res
}

# ============================================================================
# PACKAGE LOADER (called explicitly by scripts; not at source time)
# ============================================================================

#' Attach the analysis package stack. Call from scripts, not config.R.
#' @param extra Optional character vector of additional packages to attach.
load_packages <- function(extra = character(0)) {
  pkgs <- unique(c("limma", "edgeR", "ggplot2", "dplyr", "yaml", extra))
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Required package '%s' is not installed.", p))
    }
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
  invisible(TRUE)
}
