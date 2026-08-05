# 28_hypoxia_focus_viz.R (VIZ)
## Two hypoxia-focused running-enrichment panels, written to
## 03_results/06_gsea/figures/_overview/:
##
##   hypoxia_running_sum_wt_heat   HALLMARK_HYPOXIA alone in the WT_heat contrast, drawn
##                                 as the canonical three-panel RNAseq-toolkit running-sum
##                                 figure (running ES / gene-hit rug / ranked metric) via
##                                 gsea_running_sum_plot() -> enrichplot::gseaplot2().
##   hypoxia_routes_by_contrast    Four hypoxia-named gene sets, one per database, drawn as
##                                 running-ES curves faceted by the three focal contrasts
##                                 WT_heat / KO_heat / Interaction.
##
## Viz only. Reads master_gsea_table.csv for the quotable NES and padj, the cached DE
## topTables (02_de_results.rds) for the ranked vector, and the cached gene-set lists
## (geneset_{msigdb,custom}_<DB>.rds). It re-runs no GSEA, re-fits no DE, and writes nothing
## under 03_results/master/.
##
## Why these panels exist alongside the general ones. The general running sums of this stage
## draw the top figures.running_sum_top sets per (contrast x database) cell by |NES|, and
## HALLMARK_HYPOXIA is 6th by |NES| in the WT_heat Hallmark cell (+1.91) while sitting 4th by
## adjusted p (4.2e-06), so it falls outside them by construction. The general panels pin no
## set by name, so a figure that does pin one owes that provenance line in every caption it
## writes, which is what BY_NAME_NOTE and BY_NAME_NOTE_ONE below carry.
##
## The fourth set earns its curve. REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA runs negative and
## non-significant in all three contrasts, and is drawn at the same weight with its own
## colour key, because whether two sets that share a word in their names behave alike is
## something the panel shows rather than assumes.
##
## Run from the compartment root:
##   Rscript 02_analysis/scripts/28_hypoxia_focus_viz.R

# ============================================================================
# 0. Setup + style contract + RNAseq-toolkit GSEA plotters
# ============================================================================

source("02_analysis/helpers/figure_style.R")     # -> FIG_CFG, project_theme(), save_figure(),
                                                  #    save_overview(), style_series(), overview_path(),
                                                  #    write_caption(), contrast_label()
source("02_analysis/config/config.R")             # -> YAML_CONFIG, DIR_OBJECTS, DIR_MASTER,
                                                  #    GSEA_FDR_CUTOFF, GSEA_MIN_SIZE/MAX_SIZE, RANK_METRIC
source("02_analysis/helpers/de_gsea_helpers.R")   # -> load_de_results(), build_ranked_vector()

# RNAseq-toolkit GSEA plotters, with format_pathway_names ahead of the plotters that call it.
.RTK <- "01_modules/RNAseq-toolkit"
.GP  <- file.path(.RTK, "scripts", "GSEA", "GSEA_plotting")
source(file.path(.RTK, "scripts", "custom_minimal_theme.R"))  # custom_minimal_theme_with_grid()
source(file.path(.GP, "format_pathway_names.R"))              # format_pathway_name()  <- BEFORE plotters
source(file.path(.GP, "gsea_plotting_utils.R"))               # smart_wrap(), get_db_plot_params()
source(file.path(.GP, "gsea_running_sum_plot.R"))             # gsea_running_sum_plot()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(tidyr)
  library(methods)
  library(patchwork)
  library(enrichplot)   # gseaplot2() (running-ES geometry, deterministic at exponent = 1)
  library(DOSE)         # gseaResult S4 class definition
})
options(stringsAsFactors = FALSE)

# ============================================================================
# 1. Constants, all read from config
# ============================================================================

STAGE  <- "06_gsea"
SCRIPT <- "02_analysis/scripts/28_hypoxia_focus_viz.R"

FDR    <- GSEA_FDR_CUTOFF                                                     # 0.05
RSTOP  <- FIG_CFG$figures$running_sum_top %||% 5L                             # general-panel |NES| cap
RSYLIM <- as.numeric(unlist(FIG_CFG$figures$running_sum_ylim %||% c(-1, 1)))  # shared ES y-range
TOPN   <- FIG_CFG$figures$top_pathways   %||% 20L
MINGS  <- GSEA_MIN_SIZE
MAXGS  <- GSEA_MAX_SIZE
LBLSZ  <- FIG_CFG$figures$label_size %||% 4
LINEW  <- FIG_CFG$figures$line_width %||% 1.0

OI     <- FIG_CFG$colors$okabe_ito %||% list()

# One colour per drawn set, from the config's colourblind-safe palette. Four separated hues,
# so the null-running Reactome curve reads as clearly as the rest.
SET_COLOR <- c(
  HALLMARK_HYPOXIA                                   = OI$vermillion    %||% "#D55E00",
  GOMF_2_OXOGLUTARATE_DEPENDENT_DIOXYGENASE_ACTIVITY = OI$blue          %||% "#0072B2",
  GOBP_RESPONSE_TO_OXYGEN_LEVELS                     = OI$bluish_green  %||% "#009E73",
  REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA              = OI$reddish_purple %||% "#CC79A7")

# The four sets, selected by name, one per database, in legend order. `tag` is the in-panel
# short key, and the database name alone identifies a set here, since each comes from its own.
PICKS <- tibble::tribble(
  ~pathway_id,                                          ~database,   ~tag,
  "HALLMARK_HYPOXIA",                                   "Hallmark",  "Hallmark",
  "GOMF_2_OXOGLUTARATE_DEPENDENT_DIOXYGENASE_ACTIVITY", "GO_MF",     "GO MF",
  "GOBP_RESPONSE_TO_OXYGEN_LEVELS",                     "GO_BP",     "GO BP",
  "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA",              "Reactome",  "Reactome")

FOCAL <- c("WT_heat", "KO_heat", "Interaction")   # the three contrasts in figure 2
FOCUS_SET      <- "HALLMARK_HYPOXIA"              # the single set in figure 1
FOCUS_DATABASE <- "Hallmark"

CFG_KV <- sprintf(
  "thresholds.gsea_fdr=%.2g; figures.running_sum_ylim=[%.1f,%.1f]; figures.running_sum_heights; running_sum_x=rank/n_ranked; figures.running_sum_top=%d; figures.top_pathways=%d; colors.okabe_ito",
  FDR, RSYLIM[1], RSYLIM[2], RSTOP, TOPN)

# ============================================================================
# 2. Guards: master table + DE checkpoint + gene-set manifest must exist
# ============================================================================

mg_fp <- file.path(DIR_MASTER, "master_gsea_table.csv")
if (!file.exists(mg_fp))
  stop("28_hypoxia_focus_viz: master_gsea_table.csv not found at ", mg_fp,
       "\n  Run 05_gsea_msigdb_run.R and 06_gsea_custom_run.R first.")
mg <- readr::read_csv(mg_fp, show_col_types = FALSE)

req_cols <- c("pathway_id", "pathway_name", "database", "nes", "pvalue",
              "padj", "set_size", "core_enrichment", "contrast")
missing_cols <- setdiff(req_cols, colnames(mg))
if (length(missing_cols) > 0)
  stop("28_hypoxia_focus_viz: master_gsea_table.csv is missing required columns: ",
       paste(missing_cols, collapse = ", "))

# Every (set x contrast) cell this script draws has to be in the master, and a missing one is
# a hard stop, since a silently dropped curve would read as a null.
want <- tidyr::crossing(PICKS, contrast = FOCAL)
have <- mg %>% dplyr::semi_join(want, by = c("pathway_id", "database", "contrast"))
if (nrow(have) != nrow(want)) {
  gap <- want %>% dplyr::anti_join(mg, by = c("pathway_id", "database", "contrast"))
  stop("28_hypoxia_focus_viz: ", nrow(gap), " of ", nrow(want),
       " (set x contrast) cells are absent from master_gsea_table.csv: ",
       paste(sprintf("%s/%s/%s", gap$database, gap$pathway_id, gap$contrast), collapse = "; "))
}
message(sprintf("[28] Master GSEA table loaded: %d rows; all %d hypoxia (set x contrast) cells present",
                nrow(mg), nrow(want)))

DE_RESULTS <- load_de_results()                     # named list of per-contrast topTables
absent_co  <- setdiff(FOCAL, names(DE_RESULTS))
if (length(absent_co) > 0)
  stop("28_hypoxia_focus_viz: contrast(s) absent from 02_de_results.rds: ",
       paste(absent_co, collapse = ", "))

GS_MANIFEST <- {
  mp <- file.path(DIR_OBJECTS, "geneset_manifest.csv")
  if (!file.exists(mp)) stop("28_hypoxia_focus_viz: geneset_manifest.csv not found at ", mp,
                             "; run 05/06_gsea_*_run.R first.")
  readr::read_csv(mp, show_col_types = FALSE)
}

# ============================================================================
# 3. gseaResult rehydration (the toolkit-plotter input)
# ============================================================================
# The contract of 12_gsea_viz.R's as_gsearesult(contrast, database), widened from one whole
# database to an explicit (pathway_id, database) pick list, since the four sets here live in
# four databases while one gseaResult carries a single @geneSets list. @result is filled
# verbatim from master_gsea_table.csv, so the figures carry the master's own NES and padj with
# no engine drift. @geneList is the ranked vector from the cached DE topTable, @geneSets the
# cached lists restricted to the ranked universe, and @params$exponent = 1, so
# enrichplot::gseaScores() computes the running ES deterministically with no permutation re-run.

.RANKED_CACHE <- new.env(parent = emptyenv())
ranked_for <- function(co) {
  if (is.null(.RANKED_CACHE[[co]]))
    .RANKED_CACHE[[co]] <- build_ranked_vector(DE_RESULTS[[co]], RANK_METRIC)
  .RANKED_CACHE[[co]]
}

.GENESET_CACHE <- new.env(parent = emptyenv())
geneset_for <- function(db) {
  if (is.null(.GENESET_CACHE[[db]])) {
    fp <- GS_MANIFEST$cache_path[GS_MANIFEST$database == db]
    if (length(fp) == 0 || is.na(fp[1]) || !file.exists(fp[1]))
      stop("geneset_for: no cached gene-set list on disk for database '", db, "'")
    sets <- readRDS(fp[1])
    if (!is.list(sets) || is.null(names(sets)))
      stop("geneset_for: cached gene-set list for '", db, "' is not a named list")
    .GENESET_CACHE[[db]] <- sets
  }
  .GENESET_CACHE[[db]]
}

#' Build a clusterProfiler `gseaResult` for one contrast over an explicit pick list.
#' @param co       contrast name (a key of DE_RESULTS and of mg$contrast)
#' @param picks    data frame with columns pathway_id + database, in draw order
as_gsearesult <- function(co, picks) {
  rows <- mg %>%
    dplyr::semi_join(picks, by = c("pathway_id", "database")) %>%
    dplyr::filter(contrast == co, !is.na(nes))
  rows <- rows[match(picks$pathway_id, rows$pathway_id), , drop = FALSE]
  if (anyNA(rows$pathway_id))
    stop("as_gsearesult: contrast '", co, "' is missing master rows for: ",
         paste(picks$pathway_id[is.na(rows$pathway_id)], collapse = ", "))

  ranked <- ranked_for(co)

  # Gene sets, pulled from each pick's own database cache and restricted to the ranked
  # universe. The running sum needs @geneSets[[ID]] for every plotted pathway, so a set
  # missing from its cache is a hard stop ahead of a quietly undrawn curve.
  sets <- list()
  for (i in seq_len(nrow(picks))) {
    id <- picks$pathway_id[i]
    g  <- geneset_for(picks$database[i])[[id]]
    if (is.null(g))
      stop("as_gsearesult: set '", id, "' absent from the cached '", picks$database[i],
           "' gene-set list")
    sets[[id]] <- intersect(as.character(g), names(ranked))
  }

  # @result in the clusterProfiler gseaResult schema. Description holds the raw pathway_id, so
  # the toolkit's format_pathway_name() and strip_prefix render one display name everywhere.
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
    stringsAsFactors = FALSE)
  rownames(res) <- res$ID

  methods::new("gseaResult",
    result      = res,
    organism    = (if (exists("SPECIES")) SPECIES else "Mus musculus"),
    setType     = "UNKNOWN",
    geneSets    = sets[res$ID],
    geneList    = ranked,
    keytype     = "UNKNOWN",
    permScores  = matrix(),
    params      = list(pvalueCutoff = 1, eps = 0, pAdjustMethod = "fdr",
                       exponent = 1, minGSSize = MINGS, maxGSSize = MAXGS),
    gene2Symbol = character(),
    readable    = FALSE)
}

# ============================================================================
# 4. Small display helpers
# ============================================================================

#' Adjusted p for the figure face: scientific below 1e-3, else three significant digits.
fmt_p <- function(p, digits = 2) {
  vapply(p, function(x) {
    if (is.na(x)) return("NA")
    if (x < 1e-3) sprintf(paste0("%.", digits, "e"), x) else sprintf("%.3g", x)
  }, character(1))
}

#' Number of leading-edge genes behind a master row's `core_enrichment` field.
n_leading <- function(core) {
  vapply(core, function(s) {
    if (is.na(s) || !nzchar(s)) return(0L)
    length(strsplit(s, "/", fixed = TRUE)[[1]])
  }, integer(1), USE.NAMES = FALSE)
}

wrap_text <- function(x, width) paste(strwrap(as.character(x), width = width), collapse = "\n")

#' Warn if a title or subtitle will run off the canvas. Three things this check needs, each
#' of which was a false pass in the ported form of it:
#'
#' 1. Font metrics from the output device. grid::convertWidth() answers in the units of
#'    whatever device is open, and the headless default understates cairo by about 14%: the
#'    subtitle that visibly overran here measured 12.30 in on the default device under a
#'    13 in canvas and 13.99 in on the cairo device save_figure() writes to. So the
#'    measurement happens inside a throwaway cairo_pdf of the saved figure's geometry.
#' 2. The subtitle measured alongside the title, since the subtitle is what overran.
#' 3. The strings passed in, because patchwork flattens its child plots into one top-level
#'    gtable whose "title" and "subtitle" cells are zeroGrobs, so a running-sum figure whose
#'    title lives on an inner panel measures 0.00 in and reads as comfortable.
#' Both arms are measured and the larger governs, so neither can hide a clipped string.
.text_fits <- function(p, canvas_w, title = NULL, subtitle = NULL,
                       canvas_h = NULL, margin_in = 0.4) {
  fg <- FIG_CFG$figures %||% list()
  canvas_h <- as.numeric(canvas_h %||% fg$height %||% 6.5)
  spec <- list(
    list(nm = "title",    txt = title,    pt = as.numeric(fg$title_size    %||% 16), face = "bold"),
    list(nm = "subtitle", txt = subtitle, pt = as.numeric(fg$subtitle_size %||% 11), face = "plain"))
  spec <- Filter(function(s) !is.null(s$txt) && nzchar(as.character(s$txt)[1]), spec)
  if (!length(spec)) return(invisible(numeric(0)))

  # The measurement device: same family and geometry as the saved artifact.
  tmp <- tempfile(fileext = ".pdf")
  dev_ok <- tryCatch({ grDevices::cairo_pdf(tmp, width = canvas_w, height = canvas_h); TRUE },
                     error = function(e) FALSE)
  if (!dev_ok) {
    warning("[28] text-fit check could not open a cairo measurement device; ",
            "widths would be measured on the wrong font metrics, so the check is SKIPPED.")
    return(invisible(numeric(0)))
  }
  on.exit({ grDevices::dev.off(); unlink(tmp) }, add = TRUE)

  str_w <- function(txt, pt, face) {
    lines <- strsplit(as.character(txt)[1], "\n", fixed = TRUE)[[1]]
    max(vapply(lines, function(s)
      grid::convertWidth(
        grid::grobWidth(grid::textGrob(s, gp = grid::gpar(fontsize = pt, fontface = face))),
        "in", valueOnly = TRUE), numeric(1)))
  }

  gt_w <- c(title = NA_real_, subtitle = NA_real_)
  gt <- tryCatch(ggplot2::ggplotGrob(p), error = function(e) NULL)
  if (!is.null(gt) && !is.null(gt$layout)) {
    for (nm in names(gt_w)) {
      i <- which(gt$layout$name == nm)
      if (length(i))
        gt_w[[nm]] <- tryCatch(
          grid::convertWidth(grid::grobWidth(gt$grobs[[i[1]]]), "in", valueOnly = TRUE),
          error = function(e) NA_real_)
    }
  }

  out <- numeric(0)
  for (s in spec) {
    w_str <- str_w(s$txt, s$pt, s$face)
    w_gt  <- gt_w[[s$nm]]
    w     <- max(c(w_str, if (is.finite(w_gt)) w_gt else numeric(0)))
    if (w > canvas_w - margin_in)
      warning(sprintf("[28] %s is %.2f in wide on a %.2f in canvas and will be clipped.",
                      s$nm, w, canvas_w))
    else
      message(sprintf("  [28] %s %.2f in of %.2f in canvas (string %.2f, gtable %s)",
                      s$nm, w, canvas_w, w_str,
                      if (is.finite(w_gt)) sprintf("%.2f", w_gt) else "absent"))
    out[[s$nm]] <- w
  }
  invisible(out)
}

CANVAS_W      <- as.numeric(FIG_CFG$figures$width      %||% 8.5)
CANVAS_W_WIDE <- as.numeric(FIG_CFG$figures$width_wide %||% 13)

#' Display name for a database key ("GO_MF" -> "GO MF").
db_display <- function(db) gsub("_", " ", db, fixed = TRUE)

#' The stage-standard display name for a set, plus its database, wrapped for a legend key.
set_legend_label <- function(id, db, width = 34L) {
  pretty <- format_pathway_name(id, use_formatting = TRUE, strip_prefix = TRUE)
  wrap_text(sprintf("%s: %s", db_display(db), pretty), width)
}

#' Thin a running-ES polyline to its exact vertices.
#' Between two consecutive set members the running score decays by a constant step, so the
#' curve is piecewise-linear with vertices at each member's rank and the rank just before it.
#' Keeping those ranks plus the two endpoints reproduces the curve exactly while cutting
#' ~19,700 vertices per curve to a few hundred, which keeps the PDF openable in a vector editor.
thin_running_sum <- function(d) {
  hits <- which(d$position == 1L)
  keep <- sort(unique(c(1L, nrow(d), hits, pmax(hits - 1L, 1L), pmin(hits + 1L, nrow(d)))))
  d[keep, , drop = FALSE]
}

#' Running-ES geometry for one gseaResult, as a tidy frame (pathway_id, x, running_score).
#' The geometry comes from the public plotter the stage uses elsewhere:
#' enrichplot::gseaplot2(subplots = 1) returns the running-ES panel and its `$data` is the
#' fortified running-sum frame, so these curves and figure 1's share one code path. At
#' @params$exponent = 1 the running score is deterministic, so it holds to the master NES
#' printed beside it.
es_curve <- function(g, ids, thin = TRUE) {
  p_es <- enrichplot::gseaplot2(g, geneSetID = ids, subplots = 1,
                                color = unname(SET_COLOR[ids]))   # palette UNNAMED, in `ids` order
  d <- p_es$data
  need <- c("x", "runningScore", "position", "Description")
  if (!all(need %in% colnames(d)))
    stop("es_curve: the running-sum frame from enrichplot::gseaplot2() lacks ",
         paste(setdiff(need, colnames(d)), collapse = ", "),
         ". The ES-curve extraction needs updating for this enrichplot version.")
  out <- if (isTRUE(thin)) {
    thinned <- d %>%
      dplyr::group_by(Description) %>%
      dplyr::group_modify(~ thin_running_sum(.x)) %>%
      dplyr::ungroup()
    # The thinning is lossless, so assert the per-set extrema survive it.
    chk <- dplyr::full_join(
      d       %>% dplyr::group_by(Description) %>%
                  dplyr::summarise(lo = min(runningScore), hi = max(runningScore), .groups = "drop"),
      thinned %>% dplyr::group_by(Description) %>%
                  dplyr::summarise(lo2 = min(runningScore), hi2 = max(runningScore), .groups = "drop"),
      by = "Description")
    if (any(abs(chk$lo - chk$lo2) > 1e-12) || any(abs(chk$hi - chk$hi2) > 1e-12))
      stop("es_curve: running-ES thinning changed a curve's range")
    thinned
  } else d
  ## `x` is the raw rank the plotter returned; `rank_fraction` is what gets drawn, so the
  ## curves of two rankings of different length share one axis.
  out %>% dplyr::transmute(pathway_id = as.character(Description),
                           x = x, running_score = runningScore,
                           n_ranked = nrow(d) / length(ids),
                           rank_fraction = x / n_ranked)
}

# The master rows behind BOTH figures, in draw order, with the display fields attached.
FOCUS_ROWS <- mg %>%
  dplyr::semi_join(PICKS, by = c("pathway_id", "database")) %>%
  dplyr::filter(contrast %in% FOCAL) %>%
  dplyr::left_join(PICKS, by = c("pathway_id", "database")) %>%
  dplyr::mutate(
    leading_edge_n = n_leading(core_enrichment),
    set_label      = vapply(seq_along(pathway_id),
                            function(i) set_legend_label(pathway_id[i], database[i]),
                            character(1)),
    contrast       = factor(contrast, levels = FOCAL),
    pathway_id     = factor(pathway_id, levels = PICKS$pathway_id)) %>%
  dplyr::arrange(contrast, pathway_id)

SET_LABEL <- FOCUS_ROWS %>%
  dplyr::distinct(pathway_id, set_label) %>%
  { stats::setNames(.$set_label, as.character(.$pathway_id)) }
SET_LABEL <- SET_LABEL[PICKS$pathway_id]     # legend key order = PICKS order

# How many sets the focus database offered the sweep, so the single-set caption can name the
# family its one set was drawn from.
N_FOCUS_DB_SETS <- dplyr::n_distinct(
  mg$pathway_id[mg$database == FOCUS_DATABASE & mg$contrast == "WT_heat"])

# The provenance sentence that travels with every artifact this script writes: these figures
# name their sets, where the general panels of this GSEA sweep pin nothing.
BY_NAME_NOTE <- paste0(
  "SELECTION RULE: these ", nrow(PICKS), " sets were chosen BY NAME, one per database, with no ",
  "input from their adjusted p or |NES|. The general per-database and pooled panels of this GSEA ",
  "sweep pin no set by name, and no pin from this figure is carried back into them. So this is a ",
  "named-set read-out and says nothing about how these four rank among the thousands of sets they ",
  "came from.")

BY_NAME_NOTE_ONE <- paste0(
  "SELECTION RULE: this set was chosen BY NAME, and its adjusted p and |NES| played no part ",
  "in that choice. The general per-database and pooled panels of this GSEA sweep pin no set ",
  "by name, and no pin from this figure is carried back into them. This panel draws one set ",
  "out of the ", N_FOCUS_DB_SETS, " that the ", db_display(FOCUS_DATABASE), " collection ",
  "contributed to this contrast, so read the two ranks quoted above as the whole of what it ",
  "says about where this set stands among them.")

COMPOSITION_NOTE <- paste0(
  "A set's enrichment is not evidence that the program its name invokes is present; ",
  "composition would have to be established separately. The four names overlap in wording ",
  "while their gene content differs, which is why each is drawn on its own curve.")

COMPOSITION_NOTE_ONE <- paste0(
  "A set's enrichment is not evidence that the program its name invokes is present; ",
  "composition would have to be established separately. What this curve locates is where ",
  "the set's member genes sit in one ranking.")

TIER_NOTE <- paste0("Claim tier: L3 (DE and enrichment statistics). ",
                    sample_mapping_caption())

SIGN_NOTE <- paste0(
  "NES > 0 = enriched in the contrast numerator (39 °C or WT). NES < 0 = enriched in the ",
  "denominator (37 °C or cGAS-KO).")

# ============================================================================
# 5. Figure 1: HALLMARK_HYPOXIA alone, WT_heat, three-panel running sum
# ============================================================================

f1_pick <- PICKS %>% dplyr::filter(pathway_id == FOCUS_SET, database == FOCUS_DATABASE)
g1      <- as_gsearesult("WT_heat", f1_pick)
r1      <- g1@result[FOCUS_SET, ]
r1_le   <- n_leading(r1$core_enrichment)

message(sprintf("[28] figure 1 (WT_heat / %s): NES %+.4f, p %.4g, padj %.4g, setSize %d, leading edge %d",
                FOCUS_SET, r1$NES, r1$pvalue, r1$p.adjust, r1$setSize, r1_le))

# Where the curve peaks, read off the drawn geometry so the caption cannot drift from it.
f1_curve <- es_curve(g1, FOCUS_SET, thin = FALSE)
f1_peak  <- f1_curve[which.max(abs(f1_curve$running_score)), ]
f1_n     <- f1_curve$n_ranked[1]

# The legend is this plotter's one channel for NES and adjusted p: gsea_running_sum_plot()
# calls enrichplot::gseaplot2() with pvalue_table = FALSE hardcoded, and the legend text comes
# from @result$Description, which `labels` overwrites. Two lines keep the legend column off the
# panel width, and max_name_length sits past the longest line so the toolkit leaves the wrap
# where this code put it.
f1_labs <- stats::setNames(
  sprintf("%s\n%d genes ranked\nNES %+.2f,  FDR %s",
          set_legend_label(FOCUS_SET, FOCUS_DATABASE, width = 60L),
          r1$setSize, r1$NES, fmt_p(r1$p.adjust)),
  FOCUS_SET)

# Held in a variable so .text_fits() measures the same string that gets drawn. The plotter
# offers no subtitle channel.
F1_TITLE <- sprintf("%s %s running enrichment, WT heat 39 vs 37 °C",
                    db_display(FOCUS_DATABASE),
                    format_pathway_name(FOCUS_SET, use_formatting = TRUE, strip_prefix = TRUE))

# The palette stays unnamed: the vendored gsea_running_sum_plot.R warns in its @note that a
# named palette breaks enrichplot::gseaplot2() for databases whose pathway ids differ from the
# Description field.
p1 <- gsea_running_sum_plot(
  g1,
  gene_set_ids    = FOCUS_SET,
  palette         = unname(SET_COLOR[[FOCUS_SET]]),
  labels          = f1_labs,
  max_name_length = max(nchar(f1_labs)) + 1L,
  title           = F1_TITLE)

# The required post-styling step: pins the ES y-range to figures.running_sum_ylim, puts rank
# fraction on x, collects one legend outside-right, keeps x ticks on the bottom panel, and
# applies project_theme.
p1 <- style_series(p1, ylim = RSYLIM, n_ranked = f1_n, config = FIG_CFG)
.text_fits(p1, CANVAS_W, title = F1_TITLE)

f1_tbl <- FOCUS_ROWS %>%
  dplyr::filter(pathway_id == FOCUS_SET, contrast == "WT_heat") %>%
  dplyr::transmute(pathway_id = as.character(pathway_id), pathway_name, database,
                   contrast = as.character(contrast), nes, pvalue, padj,
                   set_size, leading_edge_n, n_ranked = f1_n, panel_key = tag)

save_overview(
  plot      = p1,
  stage     = STAGE,
  name      = "hypoxia_running_sum_wt_heat",
  table     = f1_tbl,
  finding   = sprintf(
    "Running enrichment of the single Hallmark set HALLMARK_HYPOXIA against the WT 39-vs-37 °C ranked %s-statistic. NES %+.4f, p %.4g, adjusted p %.4g, %d of the set's genes present in the %s-gene ranked universe, %d of them in the leading edge. The curve peaks at an enrichment score of %+.3f at rank %s, which is the top %.0f%% of the ranking, so the set's members are concentrated toward the 39 °C end of the WT ordering. The statistics on the figure are taken verbatim from master_gsea_table.csv, and the curve geometry is recomputed deterministically from the same ranked vector the sweep used, at exponent 1, with no permutation re-run. This set sits 6th by |NES| and 4th by adjusted p in the WT_heat Hallmark cell, so it falls outside the general running-sum panels of this GSEA sweep, which draw the top %d sets per cell by |NES|. Its enrichment locates where these %d genes sit in one ranking, and the gene content of the set is a separate question from the name the set carries.",
    RANK_METRIC, r1$NES, r1$pvalue, r1$p.adjust, r1$setSize,
    format(f1_n, big.mark = ","), r1_le,
    f1_peak$running_score, format(f1_peak$x, big.mark = ","),
    100 * f1_peak$x / f1_n, RSTOP, r1$setSize),
  script    = SCRIPT,
  fn        = "gsea_running_sum_plot",
  config_kv = CFG_KV,
  input     = "03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_Hallmark.rds}",
  how_to_read = paste0(
    "Three stacked panels sharing one x axis: each gene's position in the WT_heat ranked list as a ",
    "FRACTION of its length, most 39 °C-shifted at 0, most 37 °C-shifted at 1, because ranked ",
    "universes differ in length between compartments; the axis title carries this one's size. ",
    "Top panel: the running enrichment score, which steps up at each gene belonging to the set ",
    "and decays between them. Its peak is the enrichment score, and the set members left of the ",
    "peak are the leading edge. ",
    "The y range is pinned to [", RSYLIM[1], ", ", RSYLIM[2], "] so the curve stays comparable to ",
    "every other running-sum figure in this GSEA sweep. Middle panel: one tick per set member at ",
    "that member's rank, over a band showing where the ranking crosses zero. Bottom panel: the ",
    "ranked ", RANK_METRIC, "-statistic, which shows how much signal each rank carries. The ",
    "legend carries the set name, its genes present in the ranked universe, its NES and its ",
    "adjusted p. ", SIGN_NOTE, " This set is 6th by |NES| and 4th by adjusted p in its cell, so ",
    "the general running-sum panels of this GSEA sweep, which draw the top ", RSTOP,
    " per cell by |NES|, leave it out by construction. ",
    BY_NAME_NOTE_ONE, " ", COMPOSITION_NOTE_ONE, " ", TIER_NOTE),
  config    = FIG_CFG)

# ============================================================================
# 6. Figure 2: four hypoxia-named sets, running ES faceted by contrast
# ============================================================================
# A gseaResult carries one ranked list, so three objects are built, one per contrast, and
# their running-ES geometry is pooled into a single frame that facets by contrast.

es_frames <- list()
for (co in FOCAL) {
  g   <- as_gsearesult(co, PICKS)
  ids <- as.character(PICKS$pathway_id)
  cur <- es_curve(g, ids)
  es_frames[[co]] <- cur %>% dplyr::mutate(contrast = co)
  message(sprintf("[28] figure 2 / %s: %d curve vertices retained from %d ranked positions x %d sets",
                  co, nrow(cur), cur$n_ranked[1], length(ids)))
}

es_df <- dplyr::bind_rows(es_frames) %>%
  dplyr::mutate(
    contrast_lab = factor(contrast_label(contrast, short = TRUE),
                          levels = contrast_label(FOCAL, short = TRUE)),
    set_label    = factor(unname(SET_LABEL[pathway_id]), levels = unname(SET_LABEL)))

# Every contrast's own ranked-universe size, so the axis states what the fraction hides. One
# value per distinct length, and the three contrasts here rank the same genes.
N_RANKED <- sort(unique(es_df$n_ranked))

# Per-facet statistics block. A shared legend carries one label per set, so the per-contrast
# NES and adjusted p are drawn inside each facet in the set's own colour, in the legend's
# top-to-bottom order, keyed by database name. They sit in the lower band of the panel, which
# the ES range of all twelve curves stays clear of.
stat_rows <- FOCUS_ROWS %>%
  dplyr::mutate(
    contrast_lab = factor(contrast_label(as.character(contrast), short = TRUE),
                          levels = contrast_label(FOCAL, short = TRUE)),
    set_label    = factor(unname(SET_LABEL[as.character(pathway_id)]),
                          levels = unname(SET_LABEL)),
    row          = as.integer(pathway_id),
    label        = sprintf("%s   NES %+.2f,  FDR %s", tag, nes, fmt_p(padj, digits = 1)),
    x            = 0.03,
    y            = RSYLIM[1] * (0.42 + 0.155 * (row - 1L)))

es_lo <- min(es_df$running_score)
if (min(stat_rows$y) < RSYLIM[1] || max(stat_rows$y) > es_lo)
  message(sprintf("[28] NOTE: statistics block spans y [%.2f, %.2f] against a curve minimum of %.2f",
                  min(stat_rows$y), max(stat_rows$y), es_lo))

# Held in variables so the strings that are drawn are the strings that get measured.
F2_TITLE <- "Running enrichment of four hypoxia-named sets by contrast"
F2_SUB   <- SIGN_NOTE   # the sign convention, and nothing else; the rest is in the caption

p2 <- ggplot(es_df, aes(x = rank_fraction, y = running_score, colour = set_label)) +
  geom_line(linewidth = LINEW) +
  geom_text(data = stat_rows,
            aes(x = x, y = y, label = label, colour = set_label),
            hjust = 0, vjust = 0.5, size = LBLSZ, fontface = "bold",
            show.legend = FALSE, inherit.aes = FALSE) +
  scale_colour_manual(values = unname(SET_COLOR[PICKS$pathway_id]), name = NULL,
                      drop = FALSE) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
  coord_cartesian(ylim = RSYLIM) +
  facet_wrap(~ contrast_lab, nrow = 1) +
  labs(
    title    = F2_TITLE,
    subtitle = F2_SUB,
    x        = rank_fraction_lab(
      N_RANKED, of = sprintf("each contrast's ordered %s-statistic", RANK_METRIC)),
    y        = "Running enrichment score") +
  project_theme(config = FIG_CFG) +
  ggplot2::theme(legend.position = "bottom",
                 legend.direction = "horizontal",
                 # Wide enough a gutter that one facet's 1.0 tick clears the next facet's 0.0.
                 panel.spacing = ggplot2::unit(1.8, "lines")) +
  ggplot2::guides(colour = ggplot2::guide_legend(ncol = 2, byrow = TRUE))

.text_fits(p2, CANVAS_W_WIDE, title = F2_TITLE, subtitle = F2_SUB)

# Each contrast's ranked-universe size, so the axis title's count is readable from the table.
f2_n <- es_df %>% dplyr::distinct(contrast, n_ranked)

f2_tbl <- FOCUS_ROWS %>%
  dplyr::transmute(pathway_id = as.character(pathway_id), pathway_name, database,
                   contrast = as.character(contrast), nes, pvalue, padj,
                   set_size, leading_edge_n,
                   n_ranked = f2_n$n_ranked[match(as.character(contrast), f2_n$contrast)],
                   panel_key = tag)

# The per-contrast sentence, built from the drawn numbers so the caption holds to the figure.
f2_sig <- FOCUS_ROWS %>% dplyr::filter(padj < FDR)
f2_line <- function(co) {
  r <- FOCUS_ROWS %>% dplyr::filter(contrast == co) %>% dplyr::arrange(pathway_id)
  sprintf("%s: %s", contrast_label(co),
          paste(sprintf("%s NES %+.2f (padj %s)", r$tag, r$nes, fmt_p(r$padj)),
                collapse = ", "))
}

save_overview(
  plot      = p2,
  stage     = STAGE,
  name      = "hypoxia_routes_by_contrast",
  table     = f2_tbl,
  finding   = sprintf(
    "Running-enrichment curves for four hypoxia-named gene sets, one per database, across the three focal contrasts. %s. %s. %s. Three of the four sets carry a positive NES and reach FDR < %.2g in both per-genotype heat contrasts, and the fourth, REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA, carries a negative NES and no significance in any of the three; %d of the 12 (set x contrast) cells reach FDR < %.2g. Every Interaction NES is non-significant, which reads as no detectable cGAS-dependence at n=5 and never as proven cGAS-independence: the 1-df interaction term is the least-powered contrast in this design. All four sets were chosen by name, so the panel reports these four curves and ranks nothing.",
    f2_line("WT_heat"), f2_line("KO_heat"), f2_line("Interaction"),
    FDR, nrow(f2_sig), FDR),
  script    = SCRIPT,
  fn        = "geom_line / facet_wrap",
  config_kv = CFG_KV,
  input     = "03_results/master/master_gsea_table.csv + 03_results/objects/{02_de_results.rds, geneset_msigdb_{Hallmark,GO_MF,GO_BP,Reactome}.rds}",
  how_to_read = paste0(
    "One facet per contrast, and within each facet four overlaid running-enrichment curves, one ",
    "per gene set, keyed by colour in the shared legend below the panels. The x axis is each gene's ",
    "position in THAT contrast's ranked list as a FRACTION of its length, most numerator-shifted at ",
    "0, most denominator-shifted at 1, so the facets share an axis and each carries its own ",
    "ordering; a fraction rather than a rank because ranked universes differ in length between ",
    "compartments, and the axis title carries these rankings' size. A curve steps up at each gene ",
    "belonging to its set and decays between them: an early high peak means the members are packed ",
    "at the numerator end, and a curve near zero ",
    "means they are spread through the list. The y range is pinned to [",
    RSYLIM[1], ", ", RSYLIM[2], "] so every curve stays comparable to every other running-sum ",
    "figure in this GSEA sweep. Inside each facet the NES and adjusted p for that contrast are ",
    "printed in each set's own colour, in the legend's own top-to-bottom order, keyed by database ",
    "name because each of the four sets comes from a different database. The member ticks and the ",
    "ranked-metric panel are left off to keep four curves per facet legible; the companion figure ",
    "hypoxia_running_sum_wt_heat.png draws all three panels for one set. ", SIGN_NOTE, " ",
    "The Interaction facet is the cGAS-dependence read-out: a positive significant Interaction NES ",
    "is consistent with cGAS-dependent induction, and a non-significant one means no detectable ",
    "cGAS-dependence at n=5. ", BY_NAME_NOTE, " ",
    "Read the fourth curve as carefully as the other three: three of these hypoxia-named sets ",
    "carry a positive NES in both heat contrasts and the Reactome set carries a negative one, so ",
    "shared wording in two set names leaves shared behaviour an open question. ",
    COMPOSITION_NOTE, " ", TIER_NOTE),
  config    = FIG_CFG, wide = TRUE)

# ============================================================================
# 7. Final asserts
# ============================================================================
ov_fig <- overview_path(STAGE, "figures", config = FIG_CFG)
ov_tbl <- overview_path(STAGE, "tables",  config = FIG_CFG)
expect <- c("hypoxia_running_sum_wt_heat", "hypoxia_routes_by_contrast")
for (nm in expect) {
  for (ext in c("png", "pdf")) {
    fp <- file.path(ov_fig, paste0(nm, ".", ext))
    if (!file.exists(fp)) stop("28_hypoxia_focus_viz: expected figure not written: ", fp)
  }
  fp <- file.path(ov_tbl, paste0(nm, ".csv"))
  if (!file.exists(fp)) stop("28_hypoxia_focus_viz: expected source table not written: ", fp)
}

message(sprintf(
  "\n28_hypoxia_focus_viz complete.\n  Figures: %s/{%s}.{pdf,png}\n  Tables:  %s/{%s}.csv\n  Caption: 03_results/%s/README.md",
  ov_fig, paste(expect, collapse = ","), ov_tbl, paste(expect, collapse = ","), STAGE))
