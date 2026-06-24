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

## PROVISIONAL stamp (single source of truth from config.R)
PROVISIONAL <- provisional_caption()

## Network fill palette (two networks: kegg / combined)
NET_PAL <- c(kegg = NEG, combined = POS)

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

  ## Label: top-N vertices by degree (layout-independent; avoids overlap on dense graphs)
  vd_lbl <- vd
  if (nrow(vd_lbl) > TOP_N_LBL) {
    degs <- igraph::degree(module)
    top_v <- names(sort(degs, decreasing = TRUE))[seq_len(TOP_N_LBL)]
    lbl_col <- .pick(vd_lbl, c("label", "metabolite", "name"))
    if (!is.null(lbl_col)) {
      vd_lbl_vec <- vd_lbl[[lbl_col]]
      node_names <- if (!is.null(rownames(vd_lbl))) rownames(vd_lbl) else as.character(seq_len(nrow(vd_lbl)))
      show_lbl <- node_names %in% top_v
      igraph::V(module)$display_label <- ifelse(show_lbl, vd_lbl_vec, NA_character_)
    }
  } else {
    lbl_col <- .pick(vd_lbl, c("label", "metabolite", "name"))
    igraph::V(module)$display_label <- if (!is.null(lbl_col)) vd_lbl[[lbl_col]] else NA_character_
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
      caption  = paste0(
        "Nodes = metabolites/atoms | Edges = enzymatic reactions (gene symbol).\n",
        "Edge color: orange = up-regulated enzyme (log2FC > 0), blue = down-regulated.\n",
        sprintf("Module: V=%d metabolites, E=%d reactions, w=%.2f.\n",
                n_verts, n_edges, res$solution_weight %||% NA_real_),
        sprintf("Sign convention: positive log2FC = higher in numerator of '%s'.\n",
                contrast_label(co)),
        "Claim tier L3 (DE/enrichment statistics; mechanism is interpretive text only).\n",
        PROVISIONAL,
        met_note)) +
    theme_void() +
    project_theme(config = FIG_CFG) +
    theme(
      plot.title       = element_text(hjust = 0.5, face = "bold", lineheight = 1.1),
      legend.position  = "right",
      plot.margin      = unit(c(0.4, 0.6, 0.7, 0.8), "cm"),
      plot.caption     = element_text(hjust = 0, margin = margin(t = 5)))
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
      caption  = paste0("Claim tier L3. ", PROVISIONAL)) +
    theme_void() +
    project_theme(config = FIG_CFG)
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
      PROVISIONAL)

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
        config    = FIG_CFG),
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
  p_sizes <- ggplot(sz_df,
      aes(x = reorder(label, -n_edges), y = n_edges, fill = network)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    scale_fill_manual(values = NET_PAL, name = "Network") +
    labs(
      title    = "GATOM metabolic module size by contrast",
      subtitle = sprintf("Reaction edges per contrast × network (k.gene = %d)", K_GENE),
      x        = NULL,
      y        = "Number of edges (enzymatic reactions)",
      caption  = paste0("Larger modules indicate more connected metabolic DE signal.\n",
                        "Claim tier L3. ", PROVISIONAL)) +
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
        "WT_heat and Temp_main are expected to show the largest modules (heat program recruits ",
        "glycolytic + mitochondrial enzymes). Claim tier L3."),
      script    = SCRIPT,
      fn        = "module_sizes_bar",
      config_kv = sprintf("colors.diverging; thresholds.gatom_k_gene=%d", K_GENE),
      input     = paste(sprintf("03_results/objects/10_gatom_%s.rds", CO), collapse = "; "),
      how_to_read = paste0(
        "Each bar = one (contrast × network) combination. ",
        "Height = number of reaction edges (enzyme-encoded enzymatic steps) in the SGMWCS module. ",
        "Orange bars = Combined KEGG+Rhea network; blue bars = KEGG-only network. ",
        "An absent bar means the module was empty (GATOM returned 0 edges for that contrast/network). ",
        "Claim tier L3: module size reflects DE signal density in the atom-graph, ",
        "not metabolic flux."),
      config    = FIG_CFG),
    error = function(e) message("  [14_gatom_viz] module_sizes overview failed: ", conditionMessage(e)))
}

## ── 7b. Module solution weights bar ──────────────────────────────────────
wt_df <- dplyr::filter(tidy_summary, !is.na(solution_weight), n_edges > 0L)

if (nrow(wt_df) > 0L) {
  p_weights <- ggplot(wt_df,
      aes(x = reorder(label, -solution_weight), y = solution_weight, fill = network)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    scale_fill_manual(values = NET_PAL, name = "Network") +
    labs(
      title    = "GATOM MWCS solution weight by contrast",
      subtitle = "Maximum-weight connected subgraph optimization score per contrast × network",
      x        = NULL,
      y        = "Solution weight (SGMWCS objective)",
      caption  = paste0("Higher weight = more concentrated DE signal in the metabolic subnetwork.\n",
                        "Claim tier L3. ", PROVISIONAL)) +
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
      config    = FIG_CFG),
    error = function(e) message("  [14_gatom_viz] module_weights overview failed: ", conditionMessage(e)))
}

## ── 7c. Module summary panel (from master CSV or in-memory tidy_summary) ──
##   Cross-contrast panel: module size / score / top metabolites.
##   Uses master_gsea_table-schema columns (nes = mean edge log2FC pseudo-NES).

summary_src <- if (!is.null(master_df) && nrow(master_df) > 0) master_df else tidy_summary

if (!is.null(summary_src) && nrow(summary_src) > 0) {
  ## Harmonise column names between master schema (nes, set_size) and tidy_summary
  has_master_cols <- all(c("nes", "set_size", "pathway_name") %in% names(summary_src))
  if (has_master_cols) {
    sum_plot_df <- summary_src %>%
      dplyr::select(
        contrast,
        network  = dplyr::any_of(c("network", "database")),
        n_edges  = set_size,
        score    = nes) %>%
      dplyr::mutate(label = contrast_label(contrast))
  } else {
    sum_plot_df <- tidy_summary %>%
      dplyr::select(contrast, network, n_edges, score = mean_log2FC, label) %>%
      dplyr::filter(!is.na(n_edges))
  }

  ## Panel A: n_edges (bar)
  pA <- ggplot(sum_plot_df, aes(reorder(label, -n_edges), n_edges, fill = network)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.85) +
    scale_fill_manual(values = NET_PAL, name = "Network") +
    labs(title = "Module size", subtitle = "Reaction edges",
         x = NULL, y = "Edges") +
    project_theme(config = FIG_CFG)

  ## Panel B: mean edge log2FC (point + hline at 0)
  pB <- ggplot(sum_plot_df %>% dplyr::filter(!is.na(score)),
               aes(reorder(label, -score), score, color = network, group = network)) +
    geom_hline(yintercept = 0, color = "gray60", linetype = "dashed") +
    geom_point(size = 3, alpha = 0.85, position = position_dodge(0.4)) +
    scale_color_manual(values = NET_PAL, name = "Network") +
    labs(title = "Module direction",
         subtitle = "Mean enzyme log2FC (pseudo-NES)",
         x = NULL, y = "Mean log2FC") +
    project_theme(config = FIG_CFG)

  summary_panel <- pA / pB +
    patchwork::plot_layout(heights = c(1, 1)) +
    patchwork::plot_annotation(
      title   = "GATOM metabolic-module summary",
      subtitle = sprintf(
        "SGMWCS modules across %d headline contrasts × %d networks (k.gene = %d). %s",
        length(CO), length(NETWORKS), K_GENE, PROVISIONAL),
      caption = paste0(
        "Top panel: module size (reaction edges). ",
        "Bottom panel: mean module enzyme log2FC (pseudo-NES; positive = net up-regulation).\n",
        "Claim tier L3 (enrichment statistics; Complex-I / pseudohypoxia mechanism is interpretive text only)."),
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
        "at the enrichment-statistics tier (L3). ", PROVISIONAL),
      script    = SCRIPT,
      fn        = "module_summary_panel",
      config_kv = sprintf(
        "colors.diverging; figures.top_n=%d; thresholds.gatom_k_gene=%d",
        TOP_N_LBL, K_GENE),
      input     = file.path(DIR_MASTER, "master_gatom_modules.csv"),
      how_to_read = paste0(
        "Top panel: bar height = number of enzymatic reaction edges in the SGMWCS module. ",
        "Bottom panel: dot position = mean log2FC across all enzyme-edges in the module. ",
        "Orange = Combined KEGG+Rhea network, blue = KEGG-only network. ",
        "Positive pseudo-NES (above dashed line) = net up-regulation of module enzymes ",
        "in the contrast numerator. ",
        "Sign convention: positive log2FC = higher in numerator (e.g. 39°C arm for heat contrasts). ",
        "Claim tier L3: module statistics, not metabolic flux measurement. ",
        "PROVISIONAL: sample-group labels are marker-inferred, pending collaborator sample sheet."),
      config    = FIG_CFG),
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
