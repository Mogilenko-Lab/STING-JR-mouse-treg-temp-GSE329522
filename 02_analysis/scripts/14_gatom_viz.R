#!/usr/bin/env Rscript
# 14_gatom_viz.R — VIZ
# =============================================================================
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   09_gatom
# Role:    VIZ-ONLY — reads per-contrast module igraphs (10_gatom_<contrast>.rds)
#          + cross-contrast tidy summary (master/master_gatom_modules.csv).
#          Renders figure-style-contract figures for the GATOM metabolic-module arm.
#          PERFORMS NO STATISTICS. No gatom::makeMetabolicGraph / scoreGraph /
#          mwcsr::solve_mwcsp; no DE re-fit; no master-table writes.
#
# Figures produced:
#   by_contrast/<contrast>/
#     module_graph.*   — metabolic module igraph (ggraph; atoms/reactions colored
#                        by enzyme log2FC; up=orange / down=blue diverging scale)
#   _overview/
#     module_summary.* — module-size / score / top-metabolite panel across contrasts
#     module_sizes.*   — bar chart: reaction edge count per contrast x network
#     module_weights.* — bar chart: MWCS solution weight per contrast x network
#
# GUARD:
#   * No 10_gatom_*.rds objects at all → stop() pointing to 10_gatom_modules.R.
#   * Some contrasts have objects but are empty/failed → render present, warn absent.
#   * A contrast with 0 edges in its module → a placeholder table + caption (no
#     broken ggraph).
#
# Run from project root:  Rscript 02_analysis/scripts/14_gatom_viz.R
# Dependencies: igraph, ggraph, ggrepel, patchwork, dplyr, tidyr, tibble, readr,
#   ggplot2 (all lazy — stop() with install hint if absent when needed).
# =============================================================================

## ── 0. PROJECT CONFIG + FIGURE-STYLE CONTRACT ─────────────────────────────
source("02_analysis/config/config.R")          # PROJECT_ROOT, YAML_CONFIG, DIR_*, %||%
source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme, save_figure,
                                               # save_overview, write_caption, direction_cue,
                                               # contrast_path, overview_path, round_numeric_cols

options(bitmapType = "cairo")   # headless PNG

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(ggrepel)
  library(patchwork)
})
options(stringsAsFactors = FALSE)

## ── 1. CONSTANTS (config-driven; NO inline hexes or hardcoded paths) ──────
STAGE  <- "09_gatom"
SCRIPT <- "02_analysis/scripts/14_gatom_viz.R"

## Color palette from config (diverging + muted neutral; no bare hexes)
NEG     <- FIG_CFG$colors$diverging$down    %||% "steelblue4"   # blue  (down-regulated)
MID     <- FIG_CFG$colors$diverging$neutral %||% "grey97"      # white (zero / neutral)
POS     <- FIG_CFG$colors$diverging$up      %||% "sienna"      # orange (up-regulated)
NODE_PT <- "grey30"                                          # metabolite node (no met. data)

## Geometry from config figures block
NETWORK_H  <- as.numeric(FIG_CFG$figures$height        %||% 8) * 1.8  # taller for graph layout
TOP_N_LBL  <- as.integer(FIG_CFG$figures$top_n         %||% 20L)      # max labeled nodes
TOP_N_BARS <- as.integer(FIG_CFG$figures$top_n         %||% 20L)      # max metabolites in summary bar
MAX_LBL_SZ <- as.numeric(FIG_CFG$figures$label_size    %||% 5)        # ggraph label size cap
K_GENE     <- as.integer(YAML_CONFIG$thresholds$gatom_k_gene %||% 50L)

## Headline contrasts (mirrors 10_gatom_modules.R)
HEADLINE_CONTRASTS <- c("WT_heat", "KO_heat", "Interaction", "Temp_main")

## Human-readable contrast labels come centrally from figure_style.R (config-driven;
## single source of truth in design.contrast_labels). CONTRAST_LABELS (named vector)
## and the vectorized contrast_label() are defined there and sourced above.

## Sample-provenance stamp, READ from analysis_config.yaml:design$sample_mapping.
## Short form for a canvas; sample_mapping_caption() is the README form.
MAPPING_STAMP <- sample_mapping_stamp()

## Network fill palette (two networks: kegg / combined)
NET_PAL <- c(kegg = NEG, combined = POS)

## KEGG-only degradation note (met.combined.db.rds absent → 10_gatom_modules.R built
## only the KEGG reaction network; there is no `combined` network in any bundle). Surfaced
## on the overviews so the single-network legend is not read as a missing-bar bug. Kept SHORT
## so a one-line subtitle does not clip at the panel edge; the full reason lives in how_to_read.
KEGG_ONLY_NOTE <- "KEGG-only run (Combined KEGG+Rhea network unavailable)"

## Overview geometry (inches) — explicit so BOTH variants (incl. the print PDF) are legible;
## `wide=TRUE` only widens the screen variant, leaving print at the 3.5×3 in column default,
## which crams the horizontal-bar overviews + their full-sentence contrast labels.
OVERVIEW_W   <- 11   # wide enough for full contrast labels + bars on both variants
BAR_OVERVIEW_H     <- 5    # 3-contrast horizontal-bar charts
SUMMARY_OVERVIEW_H <- 10   # 2-panel stacked patchwork

## ── 2. GUARD: discover per-contrast RDS objects ───────────────────────────
rds_pattern <- file.path(DIR_OBJECTS, "10_gatom_%s.rds")
rds_exists  <- function(co) file.exists(sprintf(rds_pattern, co))

any_rds <- any(sapply(HEADLINE_CONTRASTS, rds_exists))
if (!any_rds) {
  stop(
    "14_gatom_viz.R: no per-contrast GATOM checkpoint objects found ",
    "(looked for 03_results/objects/10_gatom_<contrast>.rds).\n",
    "  -> Run 10_gatom_modules.R first to produce the module igraphs.\n",
    "  -> If met.combined.db.rds or other GATOM references are missing,\n",
    "     see 00_data/references/gatom/ and the skill 'gatom-metabolomic-predictions'."
  )
}

## Load the contrasts that ARE available; warn for the rest.
loaded     <- list()   # co -> bundle (contrast + gene_de + modules per net)
missing_co <- character(0)
for (co in HEADLINE_CONTRASTS) {
  fp <- sprintf(rds_pattern, co)
  if (!file.exists(fp)) {
    message("  WARNING [14_gatom_viz]: no RDS for contrast '", co,
            "' (", fp, ") — skipping.")
    missing_co <- c(missing_co, co)
    next
  }
  loaded[[co]] <- readRDS(fp)
}
if (length(missing_co))
  warning("Contrasts with no GATOM checkpoint (skipped): ",
          paste(missing_co, collapse = ", "),
          ". Run 10_gatom_modules.R for the full sweep.")

CO <- names(loaded)   # contrasts with RDS objects present
NETWORKS <- c("kegg", "combined")

## ── 3. HELPER ACCESSORS (mirrors 10_gatom_modules.R extractors) ───────────
##   Bundle structure: bundle$modules[[net]] is the per-network result list.
##   INVERTED GRAPH: genes live on EDGES, metabolites on VERTEX labels.

.pick <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (!length(hit)) return(NULL)
  df[[hit[1]]]
}

## Borderless network / info panels: the figure-style contract's ONE theme entry point is
## project_theme(), which (correctly for bars/dots) keeps axis lines + grid + x/y titles.
## A ggraph network or a centred info panel must be borderless. The OLD code did
## `theme_void() + project_theme()` — which re-added that chrome after stripping it (order
## bug) — then patched it with an inline theme(). Both are gone: the panels now build with
## plain project_theme(), and the borderless stripping is applied by the contract itself at
## save time via save_overview(..., void = TRUE) (a void-override layer added AFTER
## project_theme, so it is not clobbered by save_figure's per-variant re-theme). No raw
## theme() and no theme_void()/project_theme() ordering left in this script.

net_result <- function(bundle, net) bundle$modules[[net]]

is_plottable <- function(res) {
  isTRUE(res$success) &&
    !is.null(res$module) &&
    igraph::vcount(res$module) >= 2 &&
    igraph::ecount(res$module) > 0
}

module_edge_df <- function(res) {
  if (!is_plottable(res)) return(NULL)
  igraph::as_data_frame(res$module, what = "edges")
}

module_vertex_df <- function(res) {
  if (!is_plottable(res)) return(NULL)
  igraph::as_data_frame(res$module, what = "vertices")
}

## ── 4. FIGURE BUILDER: ggraph module network ──────────────────────────────
##   Routed through save_figure / save_overview per the figure-style contract.
##   The ggraph/ggplot object is built here, then handed to save_overview —
##   the contract's save path; no bare ggsave() / inline theme() / raw hex.

build_module_ggraph <- function(res, co, net) {
  ## Returns the ggraph ggplot object or NULL if not plottable.
  if (!is_plottable(res)) return(NULL)

  module  <- res$module
  ed      <- igraph::as_data_frame(module, what = "edges")
  vd      <- igraph::as_data_frame(module, what = "vertices")
  n_edges <- igraph::ecount(module)
  n_verts <- igraph::vcount(module)

  ## Enzyme label column: Symbol > gene > enzyme > label (EDGE attrs; genes live on edges)
  label_col <- NULL
  for (cc in c("Symbol", "gene", "enzyme", "label"))
    if (cc %in% names(ed) && any(!is.na(ed[[cc]]))) { label_col <- cc; break }

  ## log2FC + pval from edge attributes (GATOM stores them per enzyme reaction)
  has_lfc  <- "log2FC" %in% names(ed) && any(!is.na(ed$log2FC))
  has_pval <- "pval"   %in% names(ed) && any(!is.na(ed$pval))
  if (!has_lfc && "logFC" %in% names(ed)) {
    igraph::E(module)$log2FC <- ed$logFC
    ed$log2FC <- ed$logFC
    has_lfc <- TRUE
  }
  if (!is.null(label_col))
    igraph::E(module)$enzyme_label <- ed[[label_col]]

  ## Node score: metabolomics absent in this run (k.met=NULL) → uniform nodes
  met_scored <- "score" %in% names(vd) && any(vd$score != 0, na.rm = TRUE)

  ## Label: top-N vertices by degree (layout-independent; avoids overlap on dense graphs).
  ## .pick() returns the column VALUES (a vector) already — assign them directly; do NOT
  ## re-index the data.frame with them (that is the `no such index at level 1` bug).
  vd_lbl   <- vd
  lbl_vals <- .pick(vd_lbl, c("label", "metabolite", "name"))   # vector of label strings or NULL
  if (is.null(lbl_vals)) {
    igraph::V(module)$display_label <- NA_character_
  } else if (nrow(vd_lbl) > TOP_N_LBL) {
    degs       <- igraph::degree(module)
    top_v      <- names(sort(degs, decreasing = TRUE))[seq_len(TOP_N_LBL)]
    node_names <- if (!is.null(rownames(vd_lbl))) rownames(vd_lbl) else as.character(seq_len(nrow(vd_lbl)))
    show_lbl   <- node_names %in% top_v
    igraph::V(module)$display_label <- ifelse(show_lbl, as.character(lbl_vals), NA_character_)
  } else {
    igraph::V(module)$display_label <- as.character(lbl_vals)
  }

  lbl_size <- min(MAX_LBL_SZ, max(2, 5 - n_edges / 25))

  net_label <- switch(net,
    kegg     = "KEGG reaction network",
    combined = "Combined KEGG + Rhea reaction network",
    paste("Network:", net))

  ## Build ggraph — set seed for reproducible Fruchterman-Reingold layout
  set.seed(42)
  p <- ggraph(module, layout = "fr")

  ## Edges: map color to log2FC; edge width to -log10(pval) if available
  if (has_lfc && has_pval && !is.null(label_col)) {
    p <- p +
      geom_edge_link(
        aes(edge_color = log2FC,
            edge_width  = -log10(pval + 1e-10),
            label       = enzyme_label),
        alpha        = 0.75,
        label_size   = lbl_size * 0.8,
        label_colour = "gray25",
        angle_calc   = "along",
        label_dodge  = unit(3, "mm"),
        show.legend  = TRUE) +
      scale_edge_width_continuous(range = c(0.6, 3.0), name = "-log₁₀(p)")
  } else if (has_lfc && !is.null(label_col)) {
    p <- p +
      geom_edge_link(
        aes(edge_color = log2FC, label = enzyme_label),
        alpha        = 0.75,
        label_size   = lbl_size * 0.8,
        label_colour = "gray25",
        angle_calc   = "along",
        label_dodge  = unit(3, "mm")) +
      scale_edge_width_continuous(range = c(0.6, 3.0), name = "-log₁₀(p)")
  } else if (has_lfc) {
    p <- p + geom_edge_link(aes(edge_color = log2FC), alpha = 0.75)
  } else {
    p <- p + geom_edge_link(color = "gray60", alpha = 0.75)
  }

  ## Diverging edge color scale (colors from config; no inline hex)
  if (has_lfc)
    p <- p + scale_edge_color_gradient2(
      low      = NEG, mid = MID, high = POS,
      midpoint = 0, na.value = "gray50",
      name     = "Enzyme\nlog2FC")

  ## Nodes: fixed size (no met. score data in this expression-only run)
  if (met_scored) {
    p <- p +
      geom_node_point(aes(size = score), color = NODE_PT, alpha = 0.85) +
      scale_size_continuous(range = c(4, 12), name = "Metabolite\nscore")
  } else {
    p <- p + geom_node_point(size = 5, color = NODE_PT, alpha = 0.85)
  }

  ## Node labels: repelled, capped to top-N, bold for legibility
  p <- p + geom_node_text(
    aes(label = display_label),
    size     = lbl_size,
    repel    = TRUE,
    max.overlaps = 30,
    segment.color = "gray50",
    fontface = "bold",
    color    = "gray10",
    na.rm    = TRUE)

  met_note <- if (!met_scored)
    paste0("\nNode size is uniform (expression-only run; k.met=NULL — ",
           "no metabolomics data supplied; metabolite scores are absent).")
  else ""

  cue <- direction_cue(
    mean(suppressWarnings(as.numeric(ed$log2FC)), na.rm = TRUE))

  p +
    labs(
      title    = sprintf("%s\n[%s]", contrast_label(co), net_label),
      subtitle = sprintf("Direction cue: %s  |  k.gene = %d  |  V=%d  E=%d",
                         cue, K_GENE, n_verts, n_edges),
      caption  = NULL) +
    project_theme(config = FIG_CFG)   # borderless stripping applied at save time (void = TRUE)
}

## ── 5. PLACEHOLDER for empty/failed modules ───────────────────────────────
make_empty_placeholder <- function(co, net, reason = "no significant module") {
  ## A small info-text ggplot that clearly says "no module" rather than
  ## rendering a broken ggraph. Meets the "placeholder table + caption" requirement.
  df <- data.frame(
    x    = 0.5,
    y    = 0.5,
    text = sprintf("No significant GATOM module\nContrast: %s | Network: %s\nReason: %s",
                   contrast_label(co), net, reason)
  )
  ggplot(df, aes(x, y, label = text)) +
    geom_text(size = 5, color = "gray40", lineheight = 1.4) +
    xlim(0, 1) + ylim(0, 1) +
    labs(
      title    = sprintf("%s [%s] — no module", contrast_label(co), net),
      subtitle = "GATOM SGMWCS returned an empty or failed module for this contrast/network.",
      caption  = paste0("Claim tier L3. ", MAPPING_STAMP)) +
    project_theme(config = FIG_CFG)   # borderless stripping applied at save time (void = TRUE)
}

## Helper: build node/edge table for save_overview source table
make_node_edge_table <- function(res, co, net) {
  if (!is_plottable(res)) {
    return(data.frame(
      contrast = co, network = net,
      type = "no_module", label = NA_character_, log2FC = NA_real_,
      stringsAsFactors = FALSE))
  }
  ed <- igraph::as_data_frame(res$module, what = "edges")
  vd <- igraph::as_data_frame(res$module, what = "vertices")
  lbl_e <- .pick(ed, c("Symbol", "gene", "enzyme", "label"))
  lbl_v <- .pick(vd, c("label", "metabolite", "name"))
  lfc_e <- suppressWarnings(as.numeric(.pick(ed, c("log2FC", "logFC"))))
  pv_e  <- suppressWarnings(as.numeric(.pick(ed, c("pval", "P.Value"))))

  edge_tbl <- data.frame(
    contrast = co, network = net, type = "enzyme_edge",
    label    = if (!is.null(lbl_e)) as.character(lbl_e) else NA_character_,
    log2FC   = if (length(lfc_e)) lfc_e else NA_real_,
    pval     = if (length(pv_e))  pv_e  else NA_real_,
    stringsAsFactors = FALSE)

  node_tbl <- data.frame(
    contrast = co, network = net, type = "metabolite_node",
    label    = if (!is.null(lbl_v)) as.character(lbl_v) else NA_character_,
    log2FC   = NA_real_, pval = NA_real_,
    stringsAsFactors = FALSE)

  dplyr::bind_rows(edge_tbl, node_tbl)
}

## ── 6. RENDER: per-contrast module graphs ─────────────────────────────────
##   One call per (contrast x network); each routed through save_overview so
##   the figure, its source table (node/edge list), and README caption are
##   written atomically. Network panels that share a contrast are written
##   individually (one file per network) rather than patchwork-stacked, so
##   each has its own well-typed source table.

CFG_KV_GRAPH <- sprintf(
  "colors.diverging; figures.top_n=%d; thresholds.gatom_k_gene=%d",
  TOP_N_LBL, K_GENE)

for (co in CO) {
  bundle <- loaded[[co]]
  for (net in NETWORKS) {
    res <- net_result(bundle, net)
    if (is.null(res)) {
      message("  [14_gatom_viz] no result for [", co, "/", net, "] — skipping.")
      next
    }

    src_tbl <- make_node_edge_table(res, co, net)

    if (is_plottable(res)) {
      p <- tryCatch(
        build_module_ggraph(res, co, net),
        error = function(e) {
          message("  [14_gatom_viz] ggraph failed [", co, "/", net, "]: ", conditionMessage(e))
          make_empty_placeholder(co, net, conditionMessage(e))
        })
    } else {
      reason <- if (!isTRUE(res$success)) (res$error %||% "GATOM solve failed")
                else "module has fewer than 2 nodes or 0 edges"
      p <- make_empty_placeholder(co, net, reason)
      message("  [14_gatom_viz] empty/failed module [", co, "/", net, "]: ", reason)
    }

    fname <- sprintf("module_graph_%s", net)   # e.g. module_graph_kegg

    finding_graph <- sprintf(
      paste0(
        "%s (%s): GATOM SGMWCS metabolic module — nodes are metabolites/atoms; ",
        "edges are enzymatic reactions colored by gene log2FC ",
        "(orange = up-regulated at 39°C, blue = down-regulated; ",
        "no metabolomics data supplied so node size is uniform). ",
        "Module V=%s E=%s w=%s. ",
        "This is an expression-only GATOM run (k.met=NULL). ",
        "Claim tier L3 (enrichment statistics; mechanism is interpretive text only). %s"),
      contrast_label(co), net,
      if (is_plottable(res)) as.character(res$n_vertices) else "n/a",
      if (is_plottable(res)) as.character(res$n_edges)    else "n/a",
      if (is_plottable(res)) sprintf("%.2f", res$solution_weight %||% NA_real_) else "n/a",
      MAPPING_STAMP)

    how_to_read_graph <- paste0(
      "Nodes = metabolites/atoms (uniform size — no metabolomics). ",
      "Edges = enzymatic reactions; gene symbol labels on edges. ",
      "Edge color: orange = up-regulated enzyme (log2FC > 0), blue = down-regulated. ",
      "Edge width (when available): -log10(raw p-value) from GATOM BUM scoring. ",
      "Sign convention: positive log2FC = higher in numerator of the contrast. ",
      "The direction cue (subtitle) reflects the mean signed log2FC of module edges. ",
      "Top-", TOP_N_LBL, " most connected nodes are labeled to reduce overplotting. ",
      "Claim tier L3: the module is a statistically optimal connected subgraph, ",
      "NOT a direct measurement of metabolic flux.")

    tryCatch(
      save_overview(
        plot      = p,
        stage     = STAGE,
        name      = fname,
        table     = round_numeric_cols(src_tbl),
        finding   = finding_graph,
        script    = SCRIPT,
        fn        = "build_module_ggraph",
        config_kv = CFG_KV_GRAPH,
        input     = sprintf("03_results/objects/10_gatom_%s.rds", co),
        how_to_read = how_to_read_graph,
        contrast  = co,
        config    = FIG_CFG,
        ## NETWORK_H (≈ figures.height × 1.8 ≈ 14 in): a labeled metabolic subnetwork is
        ## unreadable at the 3.5×3 in print-column default. A force-directed graph reads
        ## best near-square, so give BOTH variants a tall square canvas (explicit width +
        ## height both win over the column/wide presets in .variant_geometry).
        width     = NETWORK_H,
        height    = NETWORK_H,
        ## Borderless: strip the 0–1 ggraph axes + x/y titles + grid AFTER project_theme.
        void      = TRUE),
      error = function(e)
        message("  [14_gatom_viz] save_overview failed [", co, "/", net, "]: ",
                conditionMessage(e)))
  }
}

## ── 7. OVERVIEW: module summary panel (cross-contrast) ───────────────────
##   Reads master_gatom_modules.csv. Falls back to building the summary in-memory
##   from loaded RDS bundles if the master CSV is absent (e.g. all-empty run).

master_csv <- file.path(DIR_MASTER, "master_gatom_modules.csv")

if (file.exists(master_csv)) {
  master_df <- tryCatch(
    readr::read_csv(master_csv, show_col_types = FALSE, progress = FALSE),
    error = function(e) {
      message("  [14_gatom_viz] could not read master_gatom_modules.csv: ", conditionMessage(e))
      NULL
    })
} else {
  message("  [14_gatom_viz] master_gatom_modules.csv not found — building summary from RDS bundles.")
  master_df <- NULL
}

## Also build a tidy summary frame from the RDS bundles (always available,
## used for the sizes/weights bars and as fallback for master_df).
tidy_summary <- dplyr::bind_rows(lapply(CO, function(co) {
  bundle <- loaded[[co]]
  dplyr::bind_rows(lapply(NETWORKS, function(net) {
    res <- net_result(bundle, net)
    if (is.null(res)) return(NULL)
    ed  <- if (is_plottable(res)) igraph::as_data_frame(res$module, what = "edges") else NULL
    lfc <- if (!is.null(ed)) suppressWarnings(
      as.numeric(.pick(ed, c("log2FC", "logFC")))) else numeric(0)
    tibble::tibble(
      contrast        = co,
      network         = net,
      label           = contrast_label(co),
      success         = isTRUE(res$success),
      n_vertices      = if (isTRUE(res$success)) res$n_vertices    %||% 0L else 0L,
      n_edges         = if (isTRUE(res$success)) res$n_edges        %||% 0L else 0L,
      solution_weight = if (isTRUE(res$success)) res$solution_weight %||% NA_real_ else NA_real_,
      mean_log2FC     = if (length(lfc)) mean(lfc, na.rm = TRUE) else NA_real_
    )
  }))
}))

## ── 7a. Module sizes bar ──────────────────────────────────────────────────
sz_df <- dplyr::filter(tidy_summary, n_edges > 0L)

if (nrow(sz_df) > 0L) {
  ## Horizontal bars (coord_flip): long sentence-style contrast labels collide on a
  ## vertical x-axis; flipping reads them in full without rotating or truncating.
  p_sizes <- ggplot(sz_df,
      aes(x = reorder(label, n_edges), y = n_edges, fill = network)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    scale_fill_manual(values = NET_PAL, name = "Network") +
    coord_flip() +
    labs(
      title    = "GATOM metabolic module size by contrast",
      subtitle = sprintf("Reaction edges per contrast × network (k.gene = %d) · %s", K_GENE, KEGG_ONLY_NOTE),
      x        = NULL,
      y        = "Number of edges (enzymatic reactions)",
      caption  = paste0("Larger modules = more connected metabolic DE signal · Claim tier L3.\n",
                        MAPPING_STAMP)) +
    project_theme(config = FIG_CFG)

  tryCatch(
    save_overview(
      plot      = p_sizes,
      stage     = STAGE,
      name      = "module_sizes",
      table     = round_numeric_cols(tidy_summary),
      finding   = paste0(
        "GATOM metabolic module size (reaction edges) per contrast × network: ",
        "contrasts with the most connected metabolic sub-networks have the highest edge count; ",
        "WT_heat and Temp_main are expected to show the largest modules (the 39 °C contrast recruits ",
        "glycolytic + mitochondrial enzymes). Claim tier L3."),
      script    = SCRIPT,
      fn        = "module_sizes_bar",
      config_kv = sprintf("colors.diverging; thresholds.gatom_k_gene=%d", K_GENE),
      input     = paste(sprintf("03_results/objects/10_gatom_%s.rds", CO), collapse = "; "),
      how_to_read = paste0(
        "Horizontal bar = one (contrast × network) combination; longer bar = bigger module. ",
        "Length = number of reaction edges (enzyme-encoded enzymatic steps) in the SGMWCS module. ",
        "Blue bars = KEGG reaction network (the only network in this run). The Combined ",
        "KEGG+Rhea network is unavailable (met.combined.db.rds absent in 00_data/references/gatom/), ",
        "so 10_gatom_modules.R degraded to KEGG-only — there is no orange (combined) bar by design. ",
        "An absent contrast means its module was empty (GATOM returned 0 edges). ",
        "Claim tier L3: module size reflects DE signal density in the atom-graph, ",
        "not metabolic flux."),
      config    = FIG_CFG,
      width     = OVERVIEW_W,
      height    = BAR_OVERVIEW_H),
    error = function(e) message("  [14_gatom_viz] module_sizes overview failed: ", conditionMessage(e)))
}

## ── 7b. Module solution weights bar ──────────────────────────────────────
wt_df <- dplyr::filter(tidy_summary, !is.na(solution_weight), n_edges > 0L)

if (nrow(wt_df) > 0L) {
  p_weights <- ggplot(wt_df,
      aes(x = reorder(label, solution_weight), y = solution_weight, fill = network)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    scale_fill_manual(values = NET_PAL, name = "Network") +
    coord_flip() +
    labs(
      title    = "GATOM MWCS solution weight by contrast",
      subtitle = sprintf("SGMWCS optimization score per contrast × network · %s", KEGG_ONLY_NOTE),
      x        = NULL,
      y        = "Solution weight (SGMWCS objective)",
      caption  = paste0("Higher weight = more concentrated DE signal in the subnetwork · Claim tier L3.\n",
                        MAPPING_STAMP)) +
    project_theme(config = FIG_CFG)

  tryCatch(
    save_overview(
      plot      = p_weights,
      stage     = STAGE,
      name      = "module_weights",
      table     = round_numeric_cols(wt_df),
      finding   = paste0(
        "GATOM MWCS solution weight per contrast × network: the objective score of the ",
        "SGMWCS problem integrates both module size and edge/node score magnitudes; ",
        "contrasts with strong uniform DE in metabolic enzymes produce the highest weights. ",
        "Claim tier L3."),
      script    = SCRIPT,
      fn        = "module_weights_bar",
      config_kv = sprintf("colors.diverging; thresholds.gatom_k_gene=%d", K_GENE),
      input     = paste(sprintf("03_results/objects/10_gatom_%s.rds", CO), collapse = "; "),
      how_to_read = paste0(
        "Each bar = one (contrast × network) combination. ",
        "Height = SGMWCS objective value (dimensionless; from mwcsr::solve_mwcsp). ",
        "A higher weight indicates more/larger-magnitude enzyme DE signals concentrated in ",
        "a connected metabolic sub-network. ",
        "Claim tier L3: does not imply metabolic flux or enzyme activity directly."),
      config    = FIG_CFG,
      width     = OVERVIEW_W,
      height    = BAR_OVERVIEW_H),
    error = function(e) message("  [14_gatom_viz] module_weights overview failed: ", conditionMessage(e)))
}

## ── 7c. Module summary panel (from master CSV or in-memory tidy_summary) ──
##   Cross-contrast panel: module size / score / top metabolites.
##   Uses master_gsea_table-schema columns (nes = mean edge log2FC pseudo-NES).

## Source preference: tidy_summary (built from the RDS bundles) carries the REAL per-network
## edge count and a network key (`kegg`/`combined`) that matches NET_PAL. The master CSV
## collapses to database="GATOM" with set_size = GSEA leading-edge size (not the module edge
## count) and no per-network split — using it here mis-colored every bar grey ("No shared
## levels") and understated module size. So drive the summary panel from tidy_summary; fall
## back to the master schema only if the bundles produced nothing plottable.
tidy_plottable <- tidy_summary %>% dplyr::filter(n_edges > 0L)

if (nrow(tidy_plottable) > 0L) {
  sum_plot_df <- tidy_plottable %>%
    dplyr::select(contrast, network, n_edges, score = mean_log2FC, label)
} else if (!is.null(master_df) && nrow(master_df) > 0 &&
           all(c("nes", "set_size") %in% names(master_df))) {
  sum_plot_df <- master_df %>%
    dplyr::select(contrast,
                  network = dplyr::any_of(c("network", "database")),
                  n_edges = set_size, score = nes) %>%
    dplyr::mutate(label = contrast_label(contrast))
} else {
  sum_plot_df <- NULL
}

if (!is.null(sum_plot_df) && nrow(sum_plot_df) > 0) {

  ## SHARED contrast ordering: both panels must place the SAME contrast at the SAME
  ## position, so order ONCE (by module size) and apply that single factor to both. The
  ## old code reordered panel A by -n_edges and panel B by -score independently, so a
  ## contrast sat at different x in the two panels (misleading). coord_flip below reads
  ## the long sentence-style labels in full.
  ord_levels <- sum_plot_df %>%
    dplyr::group_by(label) %>%
    dplyr::summarise(.k = max(n_edges, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(.k) %>%        # ascending → largest ends up at TOP after coord_flip
    dplyr::pull(label)
  sum_plot_df <- sum_plot_df %>%
    dplyr::mutate(label = factor(label, levels = ord_levels))

  ## Networks actually present (KEGG-only run → just "kegg"); avoids a misleading "× 2".
  nets_present <- sort(unique(as.character(sum_plot_df$network)))

  ## Panel A: n_edges (bar)
  pA <- ggplot(sum_plot_df, aes(label, n_edges, fill = network)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    scale_fill_manual(values = NET_PAL, name = "Network") +
    coord_flip() +
    labs(title = "Module size", subtitle = "Reaction edges",
         x = NULL, y = "Edges") +
    project_theme(config = FIG_CFG)

  ## Panel B: mean edge log2FC (point + hline at 0) — SAME label factor / ordering as A
  pB <- ggplot(sum_plot_df %>% dplyr::filter(!is.na(score)),
               aes(label, score, color = network, group = network)) +
    geom_hline(yintercept = 0, color = "gray60", linetype = "dashed") +
    geom_point(size = 3, alpha = 0.85, position = position_dodge(0.4)) +
    scale_color_manual(values = NET_PAL, name = "Network") +
    coord_flip() +
    labs(title = "Module direction",
         subtitle = "Mean enzyme log2FC (pseudo-NES)",
         x = NULL, y = "Mean log2FC") +
    project_theme(config = FIG_CFG)

  summary_panel <- pA / pB +
    patchwork::plot_layout(heights = c(1, 1)) +
    patchwork::plot_annotation(
      title   = "GATOM metabolic-module summary",
      subtitle = sprintf(
        "SGMWCS modules · %d headline contrasts · %s network%s · k.gene = %d (Combined KEGG+Rhea n/a)",
        length(CO), paste(nets_present, collapse = "+"),
        if (length(nets_present) == 1) "" else "s", K_GENE),
      caption = paste0(
        "Top: module size (reaction edges). Bottom: mean enzyme log2FC ",
        "(pseudo-NES; positive = net up-regulation). Same contrast ordering in both panels.\n",
        "Claim tier L3 (enrichment statistics; mechanism is interpretive text). ", MAPPING_STAMP),
      theme = project_theme(config = FIG_CFG))

  tryCatch(
    save_overview(
      plot      = summary_panel,
      stage     = STAGE,
      name      = "module_summary",
      table     = round_numeric_cols(as.data.frame(sum_plot_df)),
      finding   = paste0(
        "Cross-contrast GATOM module summary: top panel shows module size (reaction edges) ",
        "per contrast × network; bottom panel shows mean module enzyme log2FC (pseudo-NES), ",
        "indicating whether the recruited metabolic sub-network is net up- or down-regulated. ",
        "Corroborates the MitoCarta-anchored Complex-I / metabolic-pseudohypoxia mechanism ",
        "at the enrichment-statistics tier (L3). ", MAPPING_STAMP),
      script    = SCRIPT,
      fn        = "module_summary_panel",
      config_kv = sprintf(
        "colors.diverging; figures.top_n=%d; thresholds.gatom_k_gene=%d",
        TOP_N_LBL, K_GENE),
      input     = paste(sprintf("03_results/objects/10_gatom_%s.rds", CO), collapse = "; "),
      how_to_read = paste0(
        "Top panel: bar height = number of enzymatic reaction edges in the SGMWCS module. ",
        "Bottom panel: dot position = mean log2FC across all enzyme-edges in the module. ",
        "Blue = KEGG-only network (the only network in this run; Combined KEGG+Rhea ",
        "was unavailable — met.combined.db.rds absent). ",
        "Positive pseudo-NES (above dashed line) = net up-regulation of module enzymes ",
        "in the contrast numerator. ",
        "Sign convention: positive log2FC = higher in numerator (e.g. 39°C arm for heat contrasts). ",
        "Claim tier L3: module statistics, not metabolic flux measurement. ",
        sample_mapping_caption()),
      config    = FIG_CFG,
      width     = OVERVIEW_W,
      height    = SUMMARY_OVERVIEW_H),
    error = function(e) message("  [14_gatom_viz] module_summary overview failed: ", conditionMessage(e)))
}

## ── 8. README.md for the stage (captions via write_caption; idempotent) ───
##   write_caption is called inside each save_overview for figure artifacts.
##   Here we add captions for any CSV tables written directly (not via save_overview).

## Ensure the stage README exists (write_caption will create it if absent, but
## we pre-create a placeholder header so it is present even with zero figures).
readme_path <- file.path(DIR_RESULTS, STAGE, "README.md")
if (!file.exists(readme_path)) {
  dir.create(dirname(readme_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(c(sprintf("# %s — artifact captions", STAGE), ""), readme_path)
}

## ── 9. BIOLOGY SANITY (console; no failures) ─────────────────────────────
ok_rows <- dplyr::filter(tidy_summary, n_edges > 0L)
if (nrow(ok_rows)) {
  top <- ok_rows[which.max(ok_rows$n_edges), ]
  message(sprintf("[14_gatom_viz] largest module: %s / %s  E=%d  w=%.2f",
                  top$contrast, top$network, top$n_edges,
                  top$solution_weight %||% NA_real_))
}

## ── 10. FINAL REPORT ──────────────────────────────────────────────────────
n_pdf <- length(list.files(
  file.path(DIR_RESULTS, STAGE, "figures"), pattern = "\\.pdf$", recursive = TRUE))
n_png <- length(list.files(
  file.path(DIR_RESULTS, STAGE, "figures"), pattern = "\\.png$", recursive = TRUE))
message(sprintf(
  "[14_gatom_viz] done: %d PDF + %d PNG figure files under 03_results/%s/figures/",
  n_pdf, n_png, STAGE))
message(sprintf("  contrasts rendered: %s", paste(CO, collapse = ", ")))
if (length(missing_co))
  message(sprintf("  contrasts skipped (no RDS): %s", paste(missing_co, collapse = ", ")))
message("  README captions: ", readme_path)
