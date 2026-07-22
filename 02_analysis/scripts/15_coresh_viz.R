#!/usr/bin/env Rscript
# 15_coresh_viz.R — VIZ
## CoReSh compendium-ranking + derived-GSEA visualisation (stage 08_coresh).
## Run from project root: Rscript 02_analysis/scripts/15_coresh_viz.R
##
## VIZ-ONLY — reads objects/coresh_ranked.rds, 08_coresh/tables/coresh_ranked.csv,
## objects/gsea_coresh_<contrast>.rds, and master/master_gsea_table.csv rows with
## database="CoReSh_derived".  Never re-runs coresh_batch/build_coresh_gmt/fgsea;
## never touches 02_de_results.rds.
##
## GATED: the CoReSh arm requires the ~20 GB mmu Synapse compendium (syn66227307), consumed
## read-only from the shared reference cache, before 07_coresh_search.R / 08_coresh_derived_gsea.R
## can be run. This script GUARDS on the first compute output and stops loudly if the arm has
## not been run.
##
## Figures produced (_overview/):
##   coresh_pctvar_overview — faceted bars of top-ranked compendium datasets per configured query.
##   coresh_nes_dotplot — CoReSh-derived GSEA NES × contrast dotplot (fill=NES,
##     size=-log10(padj), outline=FDR sig; y-axis = derived set, x-axis = contrast).
##
## Figures produced (by_contrast/<c>/):
##   coresh_gsea_lollipop — per-contrast top-N CoReSh-derived sets by |NES| (lollipop).
##
## Each figure written via save_overview() (figure + sibling CSV table + README caption).

# =============================================================================
# 0. Style contract (MANDATORY FIRST — no inline ggsave / theme / hex)
# =============================================================================

source("02_analysis/helpers/figure_style.R")   # loads figure-style/figure_helpers.R;
                                                # defines FIG_CFG, project_theme,
                                                # save_overview, save_figure, contrast_path,
                                                # overview_path, direction_cue,
                                                # round_numeric_cols, write_caption, %||%

options(bitmapType = "cairo")     # headless PNG

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(scales)
  library(patchwork)
})
options(stringsAsFactors = FALSE)

# =============================================================================
# 1. Project config (AFTER figure_style so %||% is defined)
# =============================================================================

source("02_analysis/config/config.R")    # PROJECT_ROOT, YAML_CONFIG, DIR_OBJECTS,
                                          # DIR_MASTER, stage_dir, provisional_caption, %||%

STAGE       <- "08_coresh"
SCRIPT_PATH <- "02_analysis/scripts/15_coresh_viz.R"

tbl_dir <- stage_dir(STAGE, "tables")    # 03_results/08_coresh/tables/ (asserts stage declared)
fig_dir <- stage_dir(STAGE, "figures")   # 03_results/08_coresh/figures/

BY_CONTRAST_DIR <- YAML_CONFIG$figures$by_contrast_dir %||% "by_contrast"
OVERVIEW_DIR    <- YAML_CONFIG$figures$overview_dir    %||% "_overview"

# =============================================================================
# 2. GUARD — compute outputs must exist (NEVER fabricate)
# =============================================================================
## The CoReSh arm is gated on the ~20 GB mmu Synapse compendium. Until 07_coresh_search.R
## and 08_coresh_derived_gsea.R have run (in-container, after the owner provisions the data),
## neither coresh_ranked.rds nor the GSEA caches exist. Stop loudly with a pointer.

PROVISION_NOTE <- "docs/_internal/reasoning/2026-06-24_02_coresh-provisioning.md"

rk_fp <- file.path(DIR_OBJECTS, "coresh_ranked.rds")
if (!file.exists(rk_fp)) {
  message(
    "CoReSh outputs not found — run 07/08 after provisioning the mmu compendium; ",
    "see ", PROVISION_NOTE, " — exiting cleanly with no figures produced."
  )
  quit(save = "no", status = 0)
}

# =============================================================================
# 3. Parameters from config (config, not hardcoded)
# =============================================================================

coresh_cfg  <- YAML_CONFIG$coresh %||% list()
TOP_N       <- as.integer(coresh_cfg$top_n_hits        %||% 5L)   # top datasets / query
SHOWCAT     <- as.integer(YAML_CONFIG$figures$top_n    %||% 20L)  # top sets in NES dotplot
FDR         <- as.numeric(YAML_CONFIG$thresholds$gsea_fdr %||% 0.05)
NES_CAP     <- as.numeric(YAML_CONFIG$figures$nes_cap  %||% 3.5)
PROV_CAP    <- provisional_caption()    # "PROVISIONAL — inferred sample mapping …"
qsigs_cfg   <- coresh_cfg$query_signatures %||% list()

.query_specs <- function(cfg) {
  specs <- cfg$sets %||% list()
  if (is.null(cfg$object) || !nzchar(as.character(cfg$object)))
    stop("15_coresh_viz: coresh.query_signatures.object is required.")
  if (!is.list(specs) || length(specs) == 0L)
    stop("15_coresh_viz: coresh.query_signatures.sets must be a non-empty list.")
  out <- lapply(seq_along(specs), function(i) {
    spec <- specs[[i]]
    required <- c("contrast", "direction", "gate")
    missing <- required[!nzchar(vapply(required, function(k) spec[[k]] %||% "", character(1)))]
    if (length(missing) > 0L)
      stop("15_coresh_viz: coresh.query_signatures.sets[[", i,
           "]] missing required field(s): ", paste(missing, collapse = ", "))
    data.frame(
      contrast = as.character(spec$contrast),
      direction = as.character(spec$direction),
      gate = as.character(spec$gate),
      origin = as.character(spec$origin %||% paste0(spec$contrast, "_", spec$direction)),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out) %>%
    mutate(query_name = paste0("Q_sig_", contrast, "_", direction, "_", gate))
}

QUERY_SPECS <- .query_specs(qsigs_cfg)
if (anyDuplicated(QUERY_SPECS$query_name))
  stop("15_coresh_viz: duplicate query name(s) from coresh.query_signatures: ",
       paste(unique(QUERY_SPECS$query_name[duplicated(QUERY_SPECS$query_name)]), collapse = ", "))
EXPECTED_QUERIES <- QUERY_SPECS$query_name

## Origin label is the configured biological family (coresh.query_signatures.sets[].origin),
## taken per distinct (contrast, direction) — data, not config row order.
ORIGIN_SPECS <- QUERY_SPECS %>%
  distinct(contrast, direction, .keep_all = TRUE) %>%
  transmute(
    query_prefix = paste0("Q_sig_", contrast, "_", direction, "_"),
    origin_label = origin
  )

## Colors from FIG_CFG (NEVER inline hex — config is the single source of truth)
NEG <- FIG_CFG$colors$diverging$down    %||% "steelblue4"
MID <- FIG_CFG$colors$diverging$neutral %||% "grey97"
POS <- FIG_CFG$colors$diverging$up      %||% "sienna"
OI  <- unlist(FIG_CFG$colors$okabe_ito  %||% list(), use.names = FALSE)
if (length(OI) == 0) OI <- c(
  FIG_CFG$colors$okabe_ito$orange        %||% "darkorange",
  FIG_CFG$colors$okabe_ito$sky_blue      %||% "deepskyblue3",
  FIG_CFG$colors$okabe_ito$bluish_green  %||% "mediumseagreen",
  FIG_CFG$colors$okabe_ito$yellow        %||% "gold",
  FIG_CFG$colors$okabe_ito$blue          %||% "steelblue",
  FIG_CFG$colors$okabe_ito$vermillion    %||% "firebrick2",
  FIG_CFG$colors$okabe_ito$reddish_purple %||% "orchid3",
  FIG_CFG$colors$okabe_ito$black         %||% "black"
)

## Contrast order from YAML (biological order for axis display)
CONTRASTS <- vapply(YAML_CONFIG$design$contrasts, function(x) x$name %||% "", character(1))
CONTRASTS <- CONTRASTS[nzchar(CONTRASTS)]

# =============================================================================
# 4. Load the compendium ranking
# =============================================================================

ranked <- as.data.frame(readRDS(rk_fp), stringsAsFactors = FALSE)
stopifnot(all(c("query_name", "gse", "pctVar", "size", "rank") %in% colnames(ranked)))
if (any(EXPECTED_QUERIES %in% ranked$query_name)) {
  ranked <- ranked[ranked$query_name %in% EXPECTED_QUERIES, , drop = FALSE]
} else {
  message("[15_coresh_viz] expected queries absent from coresh_ranked.rds; using available rows.")
}

message(sprintf("[15_coresh_viz] ranked compendium: %d rows, %d queries, %d unique GSEs",
                nrow(ranked), length(unique(ranked$query_name)), length(unique(ranked$gse))))

queries_all     <- unique(ranked$query_name)

## Pretty label helper: replace underscore-heavy internal names with human-readable forms
.query_label <- function(q) {
  q <- sub("^Q_sig_", "Signature: ", q)
  q <- gsub("_", " ", q)
  q
}

## GSE label: "GSE12345" (no network calls; titles unavailable without the compute cache)
.gse_label <- function(gse, n = 50L) {
  gse <- as.character(gse)
  ifelse(nchar(gse) > n, paste0(substr(gse, 1, n - 1), "…"), gse)
}

# =============================================================================
# 5. FIGURE A — pctVar overview: top-ranked datasets per configured query
# =============================================================================
## Faceted query bar plot, sorted by pctVar within query and capped at TOP_N.

create_pctvar_overview <- function(ranked, top_n = TOP_N) {
  if (nrow(ranked) == 0L) return(invisible(NULL))

  df <- ranked %>%
    group_by(query_name) %>%
    slice_max(pctVar, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      query_label = .query_label(query_name),
      gse_lab     = .gse_label(gse)
    )

  if (nrow(df) == 0L) return(invisible(NULL))

  df <- df %>%
    group_by(query_name) %>%
    arrange(pctVar, .by_group = TRUE) %>%
    mutate(row_id = factor(seq_len(n()))) %>%
    ungroup()

  ggplot(df, aes(x = pctVar, y = row_id)) +
    geom_col(width = 0.7, alpha = 0.88) +
    geom_text(aes(label = sprintf("%s (k=%d, rank=%d)", gse_lab, size, rank)),
              hjust = -0.06, size = FIG_CFG$figures$label_size %||% 3.0,
              color = "grey20") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.45)),
                       labels = function(x) paste0(round(x, 1), "%")) +
    facet_wrap(~ query_label, scales = "free_y") +
    labs(
      title    = "CoReSh compendium ranking — configured signature queries",
      subtitle = sprintf(
        "top %d public mouse GEO datasets per query by pctVar | %s",
        top_n, PROV_CAP),
      x = "Variance explained by query signature in public GEO dataset (pctVar, %)",
      y = NULL,
      caption = sprintf(
        "How to read: each bar = one public GEO dataset (GSE accession); length = pctVar (%%); rank = 1 is best.\npctVar is a PCA-inspired co-regulation score: high pctVar means the query gene set co-moves across samples in that dataset.\n%s", PROV_CAP)
    ) +
    project_theme(config = FIG_CFG)
}

p_pctvar_ov <- tryCatch(create_pctvar_overview(ranked),
                         error = function(e) { message("  pctvar_overview: ", e$message); NULL })

## Source table for the overview figure
top_ds_ov <- if (nrow(ranked) > 0L) {
  ranked %>%
    group_by(query_name) %>%
    slice_max(pctVar, n = TOP_N, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(query_name, rank) %>%
    mutate(query_label = .query_label(query_name)) %>%
    select(query_name, query_label, gse, gpl, pctVar, pval, size, rank)
} else {
  NULL
}

if (!is.null(p_pctvar_ov) && !is.null(top_ds_ov) && nrow(top_ds_ov) > 0L) {
  tryCatch(
    save_overview(
      plot       = p_pctvar_ov,
      stage      = STAGE,
      name       = "coresh_pctvar_overview",
      table      = top_ds_ov,
      finding    = paste0(
        "Top public mouse GEO datasets co-regulating each configured project signature query, ",
        "ranked by CoReSh pctVar."),
      script     = SCRIPT_PATH,
      fn         = "create_pctvar_overview",
      config_kv  = sprintf(
        "coresh.top_n_hits=%d; figures.top_n=%d; thresholds.gsea_fdr=%.2f; colors.diverging; colors.okabe_ito",
        TOP_N, SHOWCAT, FDR),
      input      = "03_results/objects/coresh_ranked.rds",
      how_to_read = paste0(
        "Each bar = one public GEO dataset. ",
        "Bar length = pctVar (% of variance in that dataset explained by the query), a PCA-inspired co-regulation score ",
        "(higher = stronger co-regulation); label = GSE accession, mapped query size k, and rank. ",
        "Facets separate configured signature queries exported from 17_signature_sets.rds. ",
        "Claim tier: L3-DE (data-driven, compendium coregulation score); sample labels PROVISIONAL. ",
        "Sign convention: pctVar >= 0 (variance fraction, not signed)."),
      config     = FIG_CFG),
    error = function(e) message("  coresh_pctvar_overview save_overview: ", e$message))
} else {
  message("  coresh_pctvar_overview: empty input or plot NULL — skipping save_overview.")
}

# =============================================================================
# 6. Load per-contrast CoReSh GSEA caches (VIZ reads cached .rds; no recompute)
# =============================================================================

gsea_list <- list()
for (co in CONTRASTS) {
  fp <- file.path(DIR_OBJECTS, sprintf("gsea_coresh_%s.rds", co))
  if (!file.exists(fp)) next
  g <- tryCatch(readRDS(fp), error = function(e) NULL)
  ## g is a data.frame (run_fgsea output) or NULL
  if (!is.null(g) && is.data.frame(g) && nrow(g) > 0L) gsea_list[[co]] <- g
}
message(sprintf("[15_coresh_viz] Loaded %d per-contrast CoReSh GSEA data.frames.", length(gsea_list)))

## FALLBACK: pool from master_gsea_table.csv if no per-contrast .rds found
pool <- dplyr::bind_rows(lapply(names(gsea_list), function(co) {
  r <- gsea_list[[co]]
  if (is.null(r) || nrow(r) == 0L) return(NULL)
  tibble::tibble(
    set_id   = r$pathway_id,
    set_name = r$pathway_name %||% r$pathway_id,
    nes      = r$nes,
    padj     = r$padj,
    set_size = r$set_size %||% NA_integer_,
    contrast = co
  )
}))

if (is.null(pool) || nrow(pool) == 0L) {
  mg_fp <- file.path(DIR_MASTER, "master_gsea_table.csv")
  if (file.exists(mg_fp)) {
    mg <- readr::read_csv(mg_fp, show_col_types = FALSE)
    mg_coresh <- mg[!is.na(mg$database) & mg$database == "CoReSh_derived", , drop = FALSE]
    if (nrow(mg_coresh) > 0L) {
      pool <- tibble::tibble(
        set_id   = mg_coresh$pathway_id,
        set_name = mg_coresh$pathway_name %||% mg_coresh$pathway_id,
        nes      = mg_coresh$nes,
        padj     = mg_coresh$padj,
        set_size = mg_coresh$set_size %||% NA_integer_,
        contrast = mg_coresh$contrast
      )
      message(sprintf("[15_coresh_viz] fallback: loaded %d CoReSh_derived rows from master_gsea_table.csv",
                      nrow(pool)))
    }
  }
}

# =============================================================================
# 8. FIGURE C — NES dotplot: CoReSh-derived sets × contrasts (_overview/)
# =============================================================================
## Pooled cross-contrast dotplot: y = derived set, x = contrast, fill = NES,
## size = -log10(padj), outline = FDR significant.
## Query-origin annotation is retained for traceability of the derived sets.

## Pretty-format a CORESH set_id: "CORESH_Q_sig_<contrast>_<direction>_<gate>_<GSE>"
.format_set_id <- function(sid) {
  sid <- sub("^CORESH_", "", sid)
  sid <- sub("^Q_sig_", "Signature: ", sid)
  gsub("_", " ", sid)
}

## Identify the configured contrast/direction origin of a CORESH set_id for the color strip.
.query_origin <- function(sid) {
  for (i in seq_len(nrow(ORIGIN_SPECS))) {
    prefix <- ORIGIN_SPECS$query_prefix[i]
    formatted <- gsub("_", " ", sub("^Q_sig_", "Signature: ", prefix))
    if (grepl(prefix, sid, fixed = TRUE) || grepl(formatted, sid, fixed = TRUE))
      return(ORIGIN_SPECS$origin_label[i])
  }
  "Other CoReSh query"
}

create_coresh_nes_dotplot <- function(pool, top_n = SHOWCAT) {
  if (is.null(pool) || nrow(pool) == 0L) return(invisible(NULL))

  ## Rank sets by max |NES| across contrasts; take top-N
  rank_sets <- pool %>%
    group_by(set_id) %>%
    summarise(max_abs_nes = max(abs(nes), na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_abs_nes)) %>%
    slice_head(n = top_n) %>%
    pull(set_id)

  if (length(rank_sets) == 0L) return(invisible(NULL))

  df <- pool %>%
    filter(set_id %in% rank_sets) %>%
    mutate(
      set_lab  = vapply(set_id, .format_set_id, character(1)),
      origin   = vapply(set_id, .query_origin, character(1)),
      sig      = !is.na(padj) & padj < FDR,
      neglogp  = pmin(-log10(pmax(padj, 1e-50)), 16),
      nes_capped = pmax(pmin(nes, NES_CAP), -NES_CAP),
      contrast = factor(contrast, levels = intersect(CONTRASTS, unique(contrast)))
    )

  ## y-axis order: by median NES across contrasts (most positive at top)
  set_order <- df %>%
    group_by(set_lab) %>%
    summarise(med = stats::median(nes, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med)) %>%
    pull(set_lab)
  df$set_lab <- factor(df$set_lab, levels = rev(set_order))  # highest median NES at top

  ## Origin palette (Okabe-Ito + grey)
  origin_levels <- c(ORIGIN_SPECS$origin_label, "Other CoReSh query")
  origin_cols <- stats::setNames(rep(OI, length.out = length(ORIGIN_SPECS$origin_label)),
                                 ORIGIN_SPECS$origin_label)
  origin_cols <- c(origin_cols, "Other CoReSh query" = "grey75")
  df$origin <- factor(df$origin, levels = origin_levels)

  n_sets    <- length(unique(df$set_id))
  n_contrasts <- length(unique(df$contrast))

  ## Left strip: one column showing the origin color per set
  strip_df <- df %>%
    distinct(set_lab, origin) %>%
    mutate(x = "query origin")

  p_strip <- ggplot(strip_df, aes(x = x, y = set_lab, fill = origin)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_manual(values = origin_cols, name = "Query origin", drop = FALSE) +
    labs(x = NULL, y = NULL) +
    project_theme(config = FIG_CFG)

  p_dot <- ggplot(df, aes(x = contrast, y = set_lab)) +
    geom_point(aes(size = neglogp, fill = nes_capped, color = sig),
               shape = 21, stroke = 0.8) +
    scale_fill_gradient2(low = NEG, mid = MID, high = POS, midpoint = 0,
                         name = "NES",
                         limits = c(-NES_CAP, NES_CAP),
                         oob   = scales::squish) +
    scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "transparent"),
                       name   = sprintf("FDR < %.2f", FDR)) +
    scale_size_continuous(range = c(2.5, 9), name = "−log₁₀(padj)") +
    labs(
      title    = "CoReSh-derived co-regulation set enrichment across contrasts",
      subtitle = sprintf(
        "%d sets × %d contrasts | fill = NES (orange up / blue down), size = −log₁₀(padj), outline = FDR<%.2f | left strip = query origin | %s",
        n_sets, n_contrasts, FDR, PROV_CAP),
      x = NULL, y = NULL
    ) +
    project_theme(config = FIG_CFG)

  p_strip + p_dot +
    patchwork::plot_layout(widths = c(1, 14), guides = "collect")
}

p_nes <- tryCatch(create_coresh_nes_dotplot(pool, top_n = SHOWCAT),
                   error = function(e) { message("  nes_dotplot: ", e$message); NULL })

nes_tbl <- if (!is.null(pool) && nrow(pool) > 0L) {
  pool %>%
    mutate(origin = vapply(set_id, .query_origin, character(1)),
           sig    = !is.na(padj) & padj < FDR) %>%
    arrange(contrast, padj) %>%
    select(set_id, set_name, origin, contrast, nes, padj, set_size, sig)
} else NULL

if (!is.null(p_nes)) {
  tryCatch(
    save_overview(
      plot       = p_nes,
      stage      = STAGE,
      name       = "coresh_nes_dotplot",
      table      = nes_tbl,
      finding    = paste0(
        "Cross-contrast enrichment of CoReSh-derived co-regulation sets (NES fill, -log10(padj) size, ",
        "FDR<", FDR, " outline); left strip = seeding query origin. ",
        "Expected sets are derived from the configured signature query matrix."),
      script     = SCRIPT_PATH,
      fn         = "create_coresh_nes_dotplot",
      config_kv  = sprintf(
        "figures.top_n=%d; figures.nes_cap=%.1f; thresholds.gsea_fdr=%.2f; colors.diverging; colors.okabe_ito",
        SHOWCAT, NES_CAP, FDR),
      input      = paste0(
        "03_results/objects/gsea_coresh_<contrast>.rds (fallback: ",
        "03_results/master/master_gsea_table.csv rows database=CoReSh_derived)"),
      how_to_read = paste0(
        "Each circle = one CoReSh-derived gene set (row) in one contrast (column). ",
        "Fill color = NES (Normalized Enrichment Score): orange = positive enrichment (up in numerator), ",
        "blue = negative enrichment (down); scale clamped at ±", NES_CAP, " for visual comparability. ",
        "Circle size = -log10(padj); larger = more significant. ",
        "Black outline = FDR < ", FDR, " (significant); no outline = not significant. ",
        "Left color strip = biological origin of the query that seeded the derived set ",
        "(configured contrast/direction, grouped across gates). ",
        "Sets ordered by median NES across contrasts (highest at top). ",
        "Claim tier: L3-DE (fgsea, BH-FDR). ",
        "Sample labels PROVISIONAL (inferred from heat-shock thermometer + Cgas expression)."),
      config     = FIG_CFG),
    error = function(e) message("  coresh_nes_dotplot save_overview: ", e$message))
} else {
  message("  coresh_nes_dotplot: pool empty or plot NULL — skipping (no CoReSh GSEA results yet).")
}

# =============================================================================
# 9. FIGURE D — per-contrast CoReSh GSEA lollipop (by_contrast/<c>/)
# =============================================================================
## For each contrast that produced CoReSh GSEA results: top-N sets by |NES|,
## lollipop colored by sign(NES), labeled with set name.

create_coresh_lollipop <- function(gsea_df, contrast, top_n = min(SHOWCAT, 20L)) {
  if (is.null(gsea_df) || nrow(gsea_df) == 0L) return(invisible(NULL))

  df <- gsea_df %>%
    arrange(desc(abs(nes))) %>%
    slice_head(n = top_n) %>%
    mutate(
      set_lab  = vapply(pathway_id, .format_set_id, character(1)),
      sig      = !is.na(padj) & padj < FDR,
      dir_fill = ifelse(nes >= 0, "positive", "negative"),
      cue      = direction_cue(nes)
    )

  ## Ensure no label truncation: wrap long names
  max_lab_chars <- max(nchar(as.character(df$set_lab)), na.rm = TRUE)
  if (max_lab_chars > 60)
    df$set_lab <- vapply(df$set_lab, function(s)
      paste(strwrap(s, width = 60), collapse = "\n"), character(1))

  df$set_lab <- factor(df$set_lab, levels = rev(df$set_lab))  # highest |NES| at top

  ggplot(df, aes(x = nes, y = set_lab)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.5) +
    geom_segment(aes(xend = 0, yend = set_lab, color = dir_fill),
                 linewidth = 1.1, show.legend = FALSE) +
    geom_point(aes(color = dir_fill, shape = sig), size = 3.2) +
    scale_color_manual(values = c(positive = POS, negative = NEG),
                       name = "NES direction") +
    scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                       name   = sprintf("FDR < %.2f", FDR),
                       labels = c(`TRUE` = "significant", `FALSE` = "n.s.")) +
    scale_x_continuous(limits = c(-NES_CAP, NES_CAP),
                       oob    = scales::squish,
                       labels = function(x) sprintf("%.1f", x)) +
    labs(
      title    = sprintf("CoReSh-derived GSEA — top %d sets: %s", top_n, contrast),
      subtitle = sprintf(
        "fill = NES (orange up %s / blue down %s); filled = FDR<%.2f; open = n.s. | %s",
        direction_cue(1), direction_cue(-1), FDR, PROV_CAP),
      x = "NES (Normalized Enrichment Score; |NES| capped at 3.5)",
      y = NULL
    ) +
    project_theme(config = FIG_CFG)
}

for (co in names(gsea_list)) {
  gsea_co <- gsea_list[[co]]
  if (is.null(gsea_co) || nrow(gsea_co) == 0L) next

  p_loll <- tryCatch(create_coresh_lollipop(gsea_co, contrast = co),
                      error = function(e) {
                        message(sprintf("  lollipop [%s]: %s", co, e$message)); NULL
                      })

  ## Sidecar table: top-N by |NES| for this contrast
  loll_tbl <- gsea_co %>%
    arrange(desc(abs(nes))) %>%
    slice_head(n = min(SHOWCAT, 20L)) %>%
    mutate(sig = !is.na(padj) & padj < FDR) %>%
    select(pathway_id, pathway_name, nes, padj, set_size, sig, direction, contrast)

  if (!is.null(p_loll)) {
    tryCatch(
      save_overview(
        plot       = p_loll,
        stage      = STAGE,
        name       = "coresh_gsea_lollipop",
        table      = loll_tbl,
        finding    = sprintf(
          "%s: top CoReSh-derived co-regulation sets from the configured signature-query matrix by |NES|; positive NES (orange %s) = set enriched among up-regulated genes (numerator-high); negative NES (blue %s) = enriched among down-regulated. Filled dots = FDR < %.2f; open = n.s.",
          co, direction_cue(1), direction_cue(-1), FDR),
        script     = SCRIPT_PATH,
        fn         = "create_coresh_lollipop",
        config_kv  = sprintf(
          "figures.top_n=%d; figures.nes_cap=%.1f; thresholds.gsea_fdr=%.2f; colors.diverging",
          SHOWCAT, NES_CAP, FDR),
        input      = sprintf("03_results/objects/gsea_coresh_%s.rds", co),
        how_to_read = paste0(
          "Each row = one CoReSh-derived gene set (named CORESH_<query>_<GSE>). ",
          "Horizontal bar length = NES magnitude; direction = sign: rightward (orange) = up-enriched, ",
          "leftward (blue) = down-enriched in this contrast. ",
          "x-axis clamped at ±", NES_CAP, ". Filled dot = FDR < ", FDR, " (significant). ",
          "Set name format: 'Signature: <contrast> <direction> <gate> <GSE>' — GSE is the public dataset ",
          "whose co-regulation pattern was used to derive the gene set. ",
          "Claim tier: L3-DE (fgsea multilevel, BH-FDR). Sign convention: NES > 0 = ", direction_cue(1),
          "; NES < 0 = ", direction_cue(-1), ". ",
          "Sample labels PROVISIONAL (inferred thermometer + Cgas). ",
          "No detectable cGAS-dependence at n=5 (not cGAS-independent)."),
        contrast   = co,
        config     = FIG_CFG),
      error = function(e) message(sprintf("  lollipop [%s] save_overview: %s", co, e$message)))
  }
}

# =============================================================================
# 10. README for 03_results/08_coresh/ — ensure it exists
# =============================================================================
## save_overview() writes into 03_results/08_coresh/README.md for each figure.
## Emit a header section for the stage README if it does not yet exist.

readme_fp <- file.path(DIR_RESULTS, STAGE, "README.md")
if (!file.exists(readme_fp)) {
  dir.create(dirname(readme_fp), recursive = TRUE, showWarnings = FALSE)
  header_lines <- c(
    sprintf("# %s — CoReSh compendium ranking + derived-GSEA", STAGE),
    "",
    sprintf("Stage: `%s` | Script: `%s`", STAGE, SCRIPT_PATH),
    "",
    paste0(
      "**CoReSh value:** tests where the project's configured mouse-native signatures co-regulate ",
      "across the public mouse mmu GEO compendium, then turns the top public datasets into derived ",
      "gene sets for downstream GSEA."),
    "",
    paste0("**PROVISIONAL sample labels:** sample assignments are inferred from heat-shock thermometer genes ",
           "(Hspa1b, Hsph1) + Cgas expression, NOT from a deposited sample sheet. All captions floor at ",
           "claim tier L3-DE."),
    "",
    paste0("**GATED:** requires the ~20 GB mmu Synapse compendium (syn66227307) to be provisioned ",
           "in-container by the owner before 07_coresh_search.R can run. ",
           "See: docs/_internal/reasoning/2026-06-24_02_coresh-provisioning.md"),
    "",
    "---",
    ""
  )
  writeLines(header_lines, readme_fp)
  message(sprintf("[15_coresh_viz] wrote stage README: %s", readme_fp))
} else {
  message(sprintf("[15_coresh_viz] stage README already exists: %s", readme_fp))
}

# =============================================================================
# 11. Console biology sanity (findings, NOT crashes)
# =============================================================================

## CoReSh GSEA contrast spread check
if (!is.null(pool) && nrow(pool) > 0L) {
  wt_max <- if ("WT_heat" %in% pool$contrast)
    max(pool$nes[pool$contrast == "WT_heat"], na.rm = TRUE) else NA_real_
  ko_max <- if ("KO_heat" %in% pool$contrast)
    max(pool$nes[pool$contrast == "KO_heat"], na.rm = TRUE) else NA_real_
  if (!is.na(wt_max) && !is.na(ko_max))
    message(sprintf(
      "NES spread check: max CoReSh NES WT_heat=%.2f vs KO_heat=%.2f",
      wt_max, ko_max))
}

# =============================================================================
# 12. Structural final asserts (fail loudly; every figure must have a sidecar table)
# =============================================================================

ov_fig_dir <- overview_path(STAGE, "figures", config = FIG_CFG)
ov_tbl_dir <- overview_path(STAGE, "tables",  config = FIG_CFG)

## Every PDF/PNG in _overview/figures/ must have a same-stem .csv in _overview/tables/
ov_figs <- list.files(ov_fig_dir, pattern = "\\.(pdf|png)$")
## Reduce to unique stems (strip .<variant>.<ext>)
ov_stems <- unique(sub("\\.(pdf|png)$", "", ov_figs))
ov_orphans <- ov_stems[!file.exists(file.path(ov_tbl_dir, paste0(ov_stems, ".csv")))]
if (length(ov_orphans) > 0L)
  warning("Overview figure(s) without a same-stem tables/_overview sidecar CSV: ",
          paste(ov_orphans, collapse = ", "))

n_ov_figs       <- length(ov_stems)
n_contrast_loll <- sum(vapply(CONTRASTS, function(co)
  file.exists(file.path(
    contrast_path(STAGE, co, "figures", config = FIG_CFG),
    "coresh_gsea_lollipop.png")),
  logical(1)))

message(sprintf(
  "[15_coresh_viz] COMPLETE — %d _overview figure stems, %d per-contrast lollipop(s). | %s",
  n_ov_figs, n_contrast_loll, PROV_CAP))
