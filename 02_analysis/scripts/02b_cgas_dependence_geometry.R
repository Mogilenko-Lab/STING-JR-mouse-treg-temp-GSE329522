#!/usr/bin/env Rscript
# 02b_cgas_dependence_geometry.R -- COMPUTE
# =============================================================================
# One wide per-gene table for reading the 39 C heat response of WT against the
# 39 C heat response of cGAS-KO iTregs, plus the handful of scalars the two
# panels print in-panel (stage 03_de).
#
# Why this table exists
#   The three heat contrasts are linearly dependent by construction:
#     WT_heat = KO_heat + Interaction, and Temp_main is their shared average.
#   A reader with no linear-model vocabulary can see that relation as GEOMETRY
#   if the contrasts sit side by side per gene: plot WT_heat against KO_heat and
#   the interaction becomes distance from the identity line. That needs the
#   contrasts joined on a stable per-gene key, which is what this script does.
#   No DE model is re-fit and no threshold is re-chosen.
#
# Inputs (frozen; read only)
#   03_results/objects/02_de_results.rds       (7 limma-trend topTables)
#     The per-contrast tables/by_contrast/<c>/md.csv files carry the SAME numbers
#     in CSV form; the checkpoint is read instead so there is one parse and one
#     rounding step rather than seven.
#
# Outputs
#   03_results/03_de/tables/_overview/cgas_dependence_wide.csv    (per gene)
#   03_results/03_de/tables/_overview/cgas_dependence_stats.csv   (scalars)
#
# Method
#   Join the four heat-relevant contrasts on `ensembl` (never on gene_symbol --
#   symbols are neither unique nor stable). Significance is the FDR-ONLY gate
#   (adj.P.Val < thresholds.de_fdr, no fold-change cut-off), matching
#   tables/_overview/de_counts_summary.csv; the frozen human-projection export
#   applies a stricter |logFC| gate and therefore reports far fewer genes.
#   Correlations, the WT-on-KO regression slope, the label ranking, and the
#   panel axis limits are ALL computed here so the viz sibling only draws.
#
# Honest ceiling
#   The interaction is a 1 df term at n=5/group -- the least-powered contrast in
#   the design. A gene that fails it has NO DETECTABLE cGAS-dependence at n=5,
#   which is not the same claim as cGAS-independence. Its 23 significant genes
#   are a floor on the cGAS-dependent arm, not a measurement of its full size.
#   Claim tier: L3 (association).
#
# Run from project root:
#   Rscript 02_analysis/scripts/02b_cgas_dependence_geometry.R
# =============================================================================

source("02_analysis/config/config.R")            # YAML_CONFIG, DE_FDR, stage_dir()
source("02_analysis/helpers/de_gsea_helpers.R")  # load_de_results(), round_numeric_cols()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})
options(stringsAsFactors = FALSE)

STAGE  <- "03_de"
SCRIPT <- "02_analysis/scripts/02b_cgas_dependence_geometry.R"
OVD    <- YAML_CONFIG$figures$overview_dir %||% "_overview"
OUT_DIR <- file.path(stage_dir(STAGE, "tables"), OVD)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

FDR <- DE_FDR   # FDR-only gate; the |logFC| cut-off is deliberately NOT applied here

# The four contrasts the geometry needs, as config KEYS (display labels are the
# viz sibling's job via contrast_label()).
CO_WT    <- "WT_heat"
CO_KO    <- "KO_heat"
CO_INT   <- "Interaction"
CO_SHARE <- "Temp_main"
CONTRASTS <- c(CO_WT, CO_KO, CO_INT, CO_SHARE)

message("=================================================================")
message("02b_cgas_dependence_geometry: wide per-gene heat-contrast table")
message("  gate: adj.P.Val < ", FDR, " (FDR only, no fold-change cut-off)")
message("  join key: ensembl")
message("=================================================================")

# -----------------------------------------------------------------------------
# 1. Read the frozen topTables and reshape each to a per-contrast slice
# -----------------------------------------------------------------------------
de <- load_de_results()
missing_co <- setdiff(CONTRASTS, names(de))
if (length(missing_co) > 0L)
  stop("[02b] Missing DE contrasts in the checkpoint: ", paste(missing_co, collapse = ", "))

slice_of <- function(co) {
  tt <- de[[co]]
  need <- c("gene_symbol", "ensembl", "logFC", "adj.P.Val")
  miss <- setdiff(need, colnames(tt))
  if (length(miss) > 0L)
    stop("[02b] Contrast '", co, "' lacks column(s): ", paste(miss, collapse = ", "))
  if (anyDuplicated(tt$ensembl) > 0L)
    stop("[02b] Contrast '", co, "' has duplicated ensembl ids -- the join key must be unique.")
  out <- tibble::tibble(
    ensembl     = as.character(tt$ensembl),
    gene_symbol = as.character(tt$gene_symbol),
    logFC       = as.numeric(tt$logFC),
    adjP        = as.numeric(tt$adj.P.Val)
  )
  out$sig <- !is.na(out$adjP) & out$adjP < FDR
  names(out)[names(out) == "logFC"] <- paste0("logFC_", co)
  names(out)[names(out) == "adjP"]  <- paste0("adjP_",  co)
  names(out)[names(out) == "sig"]   <- paste0("sig_",   co)
  out
}

slices <- lapply(CONTRASTS, slice_of)
names(slices) <- CONTRASTS

# -----------------------------------------------------------------------------
# 2. Join on ensembl (inner join; every contrast is tested on the same universe)
# -----------------------------------------------------------------------------
wide <- slices[[1]]
for (co in CONTRASTS[-1]) {
  nxt <- slices[[co]]
  # gene_symbol is carried once, from the first slice; assert the symbol map agrees
  # so a silently re-keyed checkpoint cannot mislabel a point.
  chk <- dplyr::inner_join(dplyr::select(wide, ensembl, gene_symbol),
                           dplyr::select(nxt, ensembl, gene_symbol),
                           by = "ensembl", suffix = c("_a", "_b"))
  if (!identical(chk$gene_symbol_a, chk$gene_symbol_b))
    stop("[02b] ensembl->symbol map disagrees between contrasts -- refusing to join.")
  wide <- dplyr::inner_join(wide, dplyr::select(nxt, -gene_symbol), by = "ensembl")
}

n_universe <- vapply(slices, nrow, integer(1))
if (length(unique(n_universe)) != 1L || nrow(wide) != n_universe[[1]])
  stop("[02b] Join lost genes: slices had ", paste(n_universe, collapse = "/"),
       " rows, joined table has ", nrow(wide), ".")

# -----------------------------------------------------------------------------
# 3. Derived per-gene flags (plot-ready; the viz sibling adds nothing)
# -----------------------------------------------------------------------------
wide <- wide %>%
  dplyr::mutate(
    # Responds to 39 C in at least one genotype -- the denominator for "how much of
    # the heat response carries a detectable cGAS-dependence".
    heat_responsive = .data[[paste0("sig_", CO_WT)]] | .data[[paste0("sig_", CO_KO)]],
    # The 1 df cGAS-dependence test. Named for what it TESTS, never for a result.
    cgas_dependent  = .data[[paste0("sig_", CO_INT)]],
    # Signed distance from the identity line in the WT-vs-KO panel IS the interaction
    # effect (WT_heat - KO_heat = Interaction), materialized so the panel's geometric
    # claim is checkable in the table.
    wt_minus_ko     = .data[[paste0("logFC_", CO_WT)]] - .data[[paste0("logFC_", CO_KO)]]
  )

# Label ranking: significant interaction genes ordered by evidence then effect size.
# The viz caps this with figures.volcano_label_top; ranking is a compute decision.
int_order <- order(wide[[paste0("adjP_", CO_INT)]],
                   -abs(wide[[paste0("logFC_", CO_INT)]]))
wide$interaction_rank <- NA_integer_
sig_rows <- int_order[wide$cgas_dependent[int_order]]
wide$interaction_rank[sig_rows] <- seq_along(sig_rows)

wide <- wide %>% dplyr::arrange(.data[[paste0("adjP_", CO_INT)]])

# -----------------------------------------------------------------------------
# 4. Scalars the panels print (correlation, slope, counts, axis limits)
# -----------------------------------------------------------------------------
x_wt <- wide[[paste0("logFC_", CO_WT)]]
y_ko <- wide[[paste0("logFC_", CO_KO)]]
x_sh <- wide[[paste0("logFC_", CO_SHARE)]]
y_in <- wide[[paste0("logFC_", CO_INT)]]

fit_ko_on_wt <- stats::lm(y_ko ~ x_wt)

# Axis limits: symmetric, rounded up to a whole log2 unit so ticks land on integers.
sym_lim <- function(...) ceiling(max(abs(c(...)), na.rm = TRUE))
lim_wt_ko <- sym_lim(x_wt, y_ko)
lim_share <- sym_lim(x_sh)
lim_inter <- sym_lim(y_in)

# The panel's geometric claim -- distance from the identity line IS the interaction
# effect -- is an arithmetic identity in a balanced 2x2. Assert it rather than assume it.
ident_resid <- max(abs(wide$wt_minus_ko - y_in), na.rm = TRUE)
if (!(ident_resid < 1e-6))
  stop("[02b] WT_heat - KO_heat != Interaction (max residual ", ident_resid,
       ") -- the contrast set is not the balanced 2x2 the panels assume.")

count_of <- function(co) {
  s <- wide[[paste0("sig_", co)]]
  l <- wide[[paste0("logFC_", co)]]
  c(n_sig = sum(s, na.rm = TRUE),
    n_up  = sum(s & l > 0, na.rm = TRUE),
    n_down = sum(s & l < 0, na.rm = TRUE))
}
cnt <- vapply(CONTRASTS, count_of, numeric(3))

stat_row <- function(metric, scope, value, note)
  tibble::tibble(metric = metric, scope = scope, value = value, note = note)

stats_tbl <- dplyr::bind_rows(
  stat_row("n_genes", "all", nrow(wide),
           "Genes tested in every contrast (shared universe after the ensembl join)."),
  stat_row("pearson_r", paste(CO_WT, CO_KO, sep = "_vs_"),
           stats::cor(x_wt, y_ko, method = "pearson", use = "complete.obs"),
           "Correlation of the two per-genotype heat responses, one point per gene."),
  stat_row("spearman_rho", paste(CO_WT, CO_KO, sep = "_vs_"),
           stats::cor(x_wt, y_ko, method = "spearman", use = "complete.obs"),
           "Rank correlation of the two per-genotype heat responses."),
  stat_row("ols_slope", paste(CO_KO, "on", CO_WT),
           unname(stats::coef(fit_ko_on_wt)[2]),
           "Slope of KO heat logFC regressed on WT heat logFC (1 = identical response size)."),
  stat_row("ols_intercept", paste(CO_KO, "on", CO_WT),
           unname(stats::coef(fit_ko_on_wt)[1]), "Intercept of the same regression."),
  stat_row("pearson_r", paste(CO_SHARE, CO_INT, sep = "_vs_"),
           stats::cor(x_sh, y_in, method = "pearson", use = "complete.obs"),
           "Correlation of the shared temperature axis with the cGAS-dependence axis."),
  stat_row("median_abs_logfc", CO_WT, stats::median(abs(x_wt), na.rm = TRUE),
           "Typical size of the heat response in WT."),
  stat_row("median_abs_logfc", CO_INT, stats::median(abs(y_in), na.rm = TRUE),
           "Typical size of the cGAS-dependent component."),
  stat_row("n_heat_responsive", paste(CO_WT, CO_KO, sep = "_or_"),
           sum(wide$heat_responsive, na.rm = TRUE),
           "Genes significant for heat in at least one genotype."),
  stat_row("n_heat_responsive_cgas_dependent", CO_INT,
           sum(wide$heat_responsive & wide$cgas_dependent, na.rm = TRUE),
           "Of those, the ones that also pass the cGAS-dependence test."),
  stat_row("axis_lim", paste(CO_WT, CO_KO, sep = "_vs_"), lim_wt_ko,
           "Symmetric axis half-range for the WT-vs-KO panel (both axes, equal scale)."),
  stat_row("axis_lim", CO_SHARE, lim_share, "Symmetric x half-range for the decomposition panel."),
  stat_row("axis_lim", CO_INT, lim_inter, "Symmetric y half-range for the decomposition panel."),
  stat_row("y_expansion", paste(CO_SHARE, CO_INT, sep = "_over_"),
           lim_share / lim_inter,
           "How much the cGAS-dependence axis is expanded relative to the shared temperature axis."),
  dplyr::bind_rows(lapply(CONTRASTS, function(co) dplyr::bind_rows(
    stat_row("n_sig",  co, cnt["n_sig",  co], "Genes at adj.P.Val < the configured FDR, no fold-change cut-off."),
    stat_row("n_up",   co, cnt["n_up",   co], "Of those, logFC > 0."),
    stat_row("n_down", co, cnt["n_down", co], "Of those, logFC < 0.")
  ))),
  stat_row("de_fdr", "gate", FDR, "The FDR-only gate used for every count in this table."),
  stat_row("max_abs_identity_residual", "WT_heat_minus_KO_heat_vs_Interaction", ident_resid,
           "Largest per-gene departure from WT_heat - KO_heat = Interaction (0 = exact).")
)

# -----------------------------------------------------------------------------
# 5. Write
# -----------------------------------------------------------------------------
wide_path  <- file.path(OUT_DIR, "cgas_dependence_wide.csv")
stats_path <- file.path(OUT_DIR, "cgas_dependence_stats.csv")

readr::write_csv(round_numeric_cols(as.data.frame(wide)), wide_path)
readr::write_csv(round_numeric_cols(as.data.frame(stats_tbl)), stats_path)

message(sprintf("[02b] %s  (%d genes x %d columns)", wide_path, nrow(wide), ncol(wide)))
message(sprintf("[02b] %s  (%d scalars)", stats_path, nrow(stats_tbl)))
message(sprintf("[02b] WT vs KO heat logFC: Pearson r = %.4f, OLS slope = %.4f",
                stats::cor(x_wt, y_ko), unname(stats::coef(fit_ko_on_wt)[2])))
message(sprintf("[02b] %s: %d sig (%d up, %d down) -- 1 df at n=5, a floor",
                CO_INT, cnt["n_sig", CO_INT], cnt["n_up", CO_INT], cnt["n_down", CO_INT]))

stopifnot("wide table is empty"  = nrow(wide) > 0L,
          "stats table is empty" = nrow(stats_tbl) > 0L)

message("02b_cgas_dependence_geometry complete.")
