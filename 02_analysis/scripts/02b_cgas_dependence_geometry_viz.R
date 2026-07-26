#!/usr/bin/env Rscript
# 02b_cgas_dependence_geometry_viz.R -- VIZ
# =============================================================================
# Two panels that carry the cGAS-dependence argument to a reader who has never
# seen a linear-model contrast (stage 03_de).
#
# The argument
#   Warming iTregs to 39 C changes thousands of genes, and it changes almost the
#   same genes by almost the same amount whether or not cGAS is present. A small,
#   one-sided arm behaves differently: 23 genes respond more strongly to heat in
#   WT than in cGAS-KO, none the other way round.
#
#   Panel A puts the WT heat response on x and the cGAS-KO heat response on y.
#   Both axes are the SAME familiar quantity measured in two genotypes, so no
#   modelling vocabulary is needed: a gene on the diagonal responded to heat the
#   same way in both, and distance from the diagonal is exactly how much removing
#   cGAS changed its response (asserted as an identity by the compute sibling).
#   Panel B splits the same data into the shared temperature axis (x) and the
#   cGAS-dependence axis (y), so the near-independence of the two threads is
#   visible as a horizontal band.
#
# Role: DRAWS ONLY. Every number on these panels -- correlations, regression
#   slope, counts, axis limits, which genes get labelled -- is read from the
#   tables emitted by 02b_cgas_dependence_geometry.R. No statistic, join, or
#   count is computed here.
#
# Inputs (read only)
#   03_results/03_de/tables/_overview/cgas_dependence_wide.csv
#   03_results/03_de/tables/_overview/cgas_dependence_stats.csv
#
# Outputs
#   03_results/03_de/figures/_overview/heat_response_wt_vs_ko.{pdf,png}
#   03_results/03_de/figures/_overview/heat_response_shared_vs_cgas_arm.{pdf,png}
#   03_results/03_de/tables/_overview/heat_response_wt_vs_ko.csv
#   03_results/03_de/tables/_overview/heat_response_shared_vs_cgas_arm.csv
#   03_results/03_de/README.md (captions, idempotent)
#
# Honest framing (both panels)
#   The interaction is LABELLED by what it tests -- cGAS-dependence of the heat
#   response -- never by a result, and its claim floors at L3. It is a 1 df term
#   at n=5/group, the least-powered contrast in the design, so its 23 genes are a
#   FLOOR on the cGAS-dependent arm. A gene on the diagonal has no detectable
#   cGAS-dependence at n=5, which is not a claim of cGAS-independence. The arm is
#   also one-sided (23 up, 0 down); nothing here may be drawn symmetrically.
#
# Dependencies: ggplot2, dplyr, readr, ggrepel
# Run from project root (after the compute sibling):
#   Rscript 02_analysis/scripts/02b_cgas_dependence_geometry_viz.R
# =============================================================================

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme(), save_overview(),
                                                # overview_path(), write_caption(),
                                                # rasterize_axes(), contrast_label()
source("02_analysis/config/config.R")          # YAML_CONFIG, GSEA_SEED, stage_dir()

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(ggrepel)
})

# -----------------------------------------------------------------------------
# 1. Style-contract constants (all from FIG_CFG / config.R; nothing hardcoded)
# -----------------------------------------------------------------------------
STAGE  <- "03_de"
SCRIPT <- "02_analysis/scripts/02b_cgas_dependence_geometry_viz.R"

FDR     <- as.numeric(FIG_CFG$thresholds$de_fdr %||% DE_FDR)
LBL_TOP <- as.integer(FIG_CFG$figures$volcano_label_top %||% 10L)
CUE     <- as.numeric(FIG_CFG$figures$cue_size   %||% 4)
LBL     <- as.numeric(FIG_CFG$figures$label_size %||% 4)
PT      <- as.numeric(FIG_CFG$figures$point_size %||% 2.4)
SEED    <- GSEA_SEED    # deterministic ggrepel layout, so re-runs are byte-stable

## DECLARED PALETTE — the only colour literals in this script.
HL <- FIG_CFG$colors$okabe_ito$vermillion %||% "#D55E00"  # the cGAS-dependent arm
BG <- "grey72"                                            # every other gene
AX <- "grey35"                                            # reference lines / rules
HALO <- "white"                                           # gene-label outline

CO_WT    <- "WT_heat"
CO_KO    <- "KO_heat"
CO_INT   <- "Interaction"
CO_SHARE <- "Temp_main"

CFG_KV <- sprintf(
  "thresholds.de_fdr=%.2g; figures.volcano_label_top=%d; figures.point_size=%.1f; colors.okabe_ito",
  FDR, LBL_TOP, PT
)

# -----------------------------------------------------------------------------
# 2. Read the compute sibling's tables (NO recomputation)
# -----------------------------------------------------------------------------
OVD      <- FIG_CFG$figures$overview_dir %||% "_overview"
TBL_OVW  <- overview_path(STAGE, "tables", FIG_CFG)
wide_path  <- file.path(TBL_OVW, "cgas_dependence_wide.csv")
stats_path <- file.path(TBL_OVW, "cgas_dependence_stats.csv")

for (p in c(wide_path, stats_path)) {
  if (!file.exists(p))
    stop(sprintf("Missing %s -- run 02_analysis/scripts/02b_cgas_dependence_geometry.R first.", p))
}

gw <- readr::read_csv(wide_path,  show_col_types = FALSE)
st <- readr::read_csv(stats_path, show_col_types = FALSE)

#' Look up one precomputed scalar. A row/column read, never a calculation. `subset`
#' names the gene set the statistic was measured on -- agreement between the two
#' genotypes is reported on heat-responsive genes, not on the whole universe.
S <- function(metric, scope, subset = "all_genes") {
  v <- st$value[st$metric == metric & st$scope == scope & st$subset == subset]
  if (length(v) != 1L)
    stop(sprintf("[02b_viz] expected exactly one '%s' / '%s' / '%s' row in the stats table, found %d.",
                 metric, scope, subset, length(v)))
  v
}
#' Same, formatted as a thousands-separated whole number for in-panel text.
N <- function(metric, scope, subset = "all_genes")
  format(as.integer(round(S(metric, scope, subset))), big.mark = ",", trim = TRUE)

HR   <- "heat_responsive"       # the scope every agreement statistic is reported on
ARM  <- "cgas_dependent_arm"    # the scope the arm anatomy is reported on
KO_ON_WT <- paste(CO_KO, "on", CO_WT)

WT_KO   <- paste(CO_WT, CO_KO, sep = "_vs_")
LIM_A   <- S("axis_lim", WT_KO)          # symmetric, equal on both axes (Panel A)
LIM_BX  <- S("axis_lim", CO_SHARE)       # Panel B x half-range
LIM_BY  <- S("axis_lim", CO_INT)         # Panel B y half-range

# Point cloud vs highlighted arm: one is 19,679 points, the other 23. Both sizes
# derive from figures.point_size so the pair scales with the contract.
PT_BG <- PT * 0.28
PT_HL <- PT * 0.95

# Legend keys are sentences, not contrast algebra. Three keys, because the arm splits
# into genes that keep the same direction and genes that REVERSE without cGAS -- the
# split is the plainest reading of the panel, so it gets its own glyph and is countable.
KEY_BG  <- "no detectable difference between genotypes at n=5"
KEY_HL  <- sprintf("heat response differs between genotypes (adj.P < %.2g)", FDR)
KEY_REV <- "reverses direction: up with heat in WT, down without cGAS"
KEYS    <- c(KEY_BG, KEY_HL, KEY_REV)
# arm_class comes from the compute sibling; the viz only attaches display labels.
CLASS_KEY <- c(no_detectable_difference = KEY_BG,
               differs_same_direction   = KEY_HL,
               differs_and_reverses     = KEY_REV)
SHAPES <- setNames(c(16, 16, 17), KEYS)   # filled circle / circle / triangle
FILLS  <- setNames(c(BG, HL, HL), KEYS)

gw <- gw %>%
  dplyr::mutate(arm = factor(unname(CLASS_KEY[arm_class]), levels = KEYS))

gw_bg   <- gw %>% dplyr::filter(arm_class == "no_detectable_difference")
gw_same <- gw %>% dplyr::filter(arm_class == "differs_same_direction")
gw_rev  <- gw %>% dplyr::filter(arm_class == "differs_and_reverses")
gw_lab  <- gw %>% dplyr::filter(!is.na(label_rank), label_rank <= LBL_TOP)

message(sprintf("[02b_viz] %s genes; %d in the arm (%d reverse); %d labelled (cap %d)",
                N("n_genes", "all"), nrow(gw_same) + nrow(gw_rev), nrow(gw_rev),
                nrow(gw_lab), LBL_TOP))

# The one-sidedness of the arm is a finding, not a rounding: state it wherever the
# arm is described, and never let a panel read as symmetric.
ARM_SENTENCE <- sprintf(
  "All %s genes in the arm respond to heat more strongly in WT than in cGAS-KO; %s do the reverse.",
  N("n_up", CO_INT), N("n_down", CO_INT))
REVERSE_SENTENCE <- sprintf(
  paste0("Every one of them falls with heat once cGAS is gone (%s of %s individually ",
         "significant), and %s rise with heat in WT instead -- an interferon response that heat ",
         "sustains only when cGAS is present."),
  N("n_down_with_heat_in_ko_significant", CO_INT, ARM),
  N("n_down_with_heat_in_ko", CO_INT, ARM),
  N("n_reverses_without_cgas", CO_INT, ARM))
FLOOR_SENTENCE <- sprintf(
  paste0("The genotype comparison of the heat response is a single-degree-of-freedom test at ",
         "n=5 per group, the least-powered contrast in this design, so %s genes is a floor on ",
         "the cGAS-dependent arm rather than its full size. Genes on the shared axis have no ",
         "detectable cGAS-dependence at n=5, which is a weaker statement than independence."),
  N("n_sig", CO_INT))

# Wrap reader-facing panel text at the contract's caption column so a sentence can
# never run off the right edge of the canvas (a truncated label is unreadable).
WRAP <- as.integer(FIG_CFG$figures$caption_wrap_column %||% 70)
wrap_at <- function(txt, width = WRAP) paste(strwrap(txt, width = width), collapse = "\n")

# -----------------------------------------------------------------------------
# 3. PANEL A — the heat response in WT against the heat response in cGAS-KO
# -----------------------------------------------------------------------------
# Agreement is quoted on the HEAT-RESPONSIVE genes, so the number cannot be carried by
# the mass of unchanged genes sitting at the origin. The regression slope leads: it is
# the direct "same response size" statistic and is not inflated by sample size.
stats_box_a <- paste(
  sprintf("%s genes, one dot each", N("n_genes", "all")),
  sprintf("39 °C changes %s genes in WT, %s in cGAS-KO",
          N("n_sig", CO_WT), N("n_sig", CO_KO)),
  sprintf("among heat-responsive genes the cGAS-KO response"),
  sprintf("is %.2f× the WT response (r = %.2f, n = %s)",
          S("ols_slope", KO_ON_WT, HR), S("pearson_r", WT_KO, HR),
          N("n_genes", "all", HR)),
  sprintf("%s differ by genotype; %s of them reverse (triangles)",
          N("n_sig", CO_INT), N("n_reverses_without_cgas", CO_INT, ARM)),
  sep = "\n")

p_a <- ggplot(mapping = aes(x = .data[[paste0("logFC_", CO_WT)]],
                            y = .data[[paste0("logFC_", CO_KO)]])) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = BG) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = BG) +
  geom_point(data = gw_bg, aes(colour = arm, shape = arm), size = PT_BG, alpha = 0.35) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, colour = AX) +
  geom_point(data = gw_same, aes(colour = arm, shape = arm), size = PT_HL) +
  geom_point(data = gw_rev,  aes(colour = arm, shape = arm), size = PT_HL * 1.25) +
  ggrepel::geom_text_repel(
    data = gw_lab, aes(label = gene_symbol), colour = HL, size = LBL,
    seed = SEED, max.overlaps = Inf, force = 8, force_pull = 0.4,
    box.padding = 0.9, point.padding = 0.6, min.segment.length = 0,
    segment.colour = HL, segment.linewidth = 0.4,
    # White halo: gene names stay readable where a leader line or a neighbouring
    # dot passes under them, without a filled box hiding data.
    bg.colour = HALO, bg.r = 0.14, show.legend = FALSE) +
  # Plain-language reading instructions, placed as fractions of the axis range.
  annotate("text", x = -0.97 * LIM_A, y = 0.97 * LIM_A, hjust = 0, vjust = 1,
           size = CUE, label = stats_box_a) +
  annotate("text", x = 0.40 * LIM_A, y = 0.40 * LIM_A, angle = 45, vjust = -0.9,
           size = CUE, colour = AX, label = "on the line: same response in both genotypes") +
  # The quadrant reading in plain language. Triangles sit right of 0 and below 0 by
  # construction, and the legend key names what that means, so this block only has to
  # point at the region -- the count lives in the stats box, never implied as all 23.
  annotate("text", x = 0.98 * LIM_A, y = -0.78 * LIM_A, hjust = 1, vjust = 1,
           size = CUE, colour = HL,
           label = paste("triangles, right of 0 and below 0:",
                         "heat raises these in WT, lowers them without cGAS",
                         sep = "\n")) +
  annotate("text", x = 0.98 * LIM_A, y = -0.97 * LIM_A, hjust = 1, vjust = 1,
           size = CUE, colour = AX,
           label = "distance from the line = WT minus cGAS-KO") +
  scale_colour_manual(values = FILLS,  breaks = KEYS, name = NULL, drop = FALSE) +
  scale_shape_manual(values  = SHAPES, breaks = KEYS, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = PT_HL * 2, alpha = 1)),
         shape  = guide_legend(ncol = 1, override.aes = list(size = PT_HL * 2, alpha = 1))) +
  coord_fixed(ratio = 1, xlim = c(-LIM_A, LIM_A), ylim = c(-LIM_A, LIM_A)) +
  labs(
    title    = paste("Warming to 39 °C changes the same genes",
                      "with and without cGAS — apart from one interferon arm", sep = "\n"),
    subtitle = wrap_at(paste0(
      "Each dot is one gene. Both axes are its log2 fold change at 39 versus 37 °C, ",
      "measured in WT (x) and in cGAS-KO (y). The dashed line marks an identical ",
      "response in the two genotypes.")),
    x = sprintf("%s — log2 fold change", contrast_label(CO_WT)),
    y = sprintf("%s — log2 fold change", contrast_label(CO_KO)),
    caption = wrap_at(paste(ARM_SENTENCE, REVERSE_SENTENCE, FLOOR_SENTENCE))
  ) +
  project_theme(config = FIG_CFG) +
  ggplot2::theme(legend.position = "bottom")

p_a <- rasterize_axes(p_a, config = FIG_CFG)   # dense dot layer -> raster inside the vector PDF

tbl_a <- gw %>%
  dplyr::select(gene_symbol, ensembl,
                dplyr::all_of(c(paste0("logFC_", CO_WT), paste0("adjP_", CO_WT),
                                paste0("logFC_", CO_KO), paste0("adjP_", CO_KO),
                                paste0("logFC_", CO_INT), paste0("adjP_", CO_INT))),
                cgas_dependent, interaction_rank)

save_overview(
  plot      = p_a,
  stage     = STAGE,
  name      = "heat_response_wt_vs_ko",
  table     = tbl_a,
  finding   = sprintf(
    paste0("Warming iTregs to 39 °C changes %s genes in WT and %s in cGAS-KO, and it changes ",
           "them in step: among the %s heat-responsive genes the cGAS-KO response is %.2f× the ",
           "WT response (r = %.2f, Spearman %.2f), so the dots pile onto the identity line. ",
           "%s %s Claim tier: L3."),
    N("n_sig", CO_WT), N("n_sig", CO_KO), N("n_genes", "all", HR),
    S("ols_slope", KO_ON_WT, HR), S("pearson_r", WT_KO, HR), S("spearman_rho", WT_KO, HR),
    ARM_SENTENCE, REVERSE_SENTENCE
  ),
  script    = SCRIPT,
  fn        = "geom_point + geom_abline (effect-versus-effect scatter, equal axes)",
  config_kv = CFG_KV,
  input     = "03_results/03_de/tables/_overview/cgas_dependence_wide.csv + cgas_dependence_stats.csv",
  how_to_read = sprintf(
    paste0("One dot per gene. x = log2 fold change at 39 vs 37 °C in WT, y = the same in ",
           "cGAS-KO, equal scales so the dashed identity line runs at 45°. Grey circles = no ",
           "detectable difference between genotypes at n=5. Vermillion circles = the heat ",
           "response differs (adj.P < %.2g); triangles = it also reverses sign, positive in WT ",
           "and negative in cGAS-KO, so they sit right of 0 and below 0. %s of those %s ",
           "reversals are individually significant in WT on their own. Labels put reversing ",
           "genes first, then evidence, capped at %d. Distance from the line = WT minus cGAS-KO. ",
           "Gate: FDR only. Claim tier: L3. PROVISIONAL sample labels; n=5/group."),
    FDR, N("n_reverses_wt_significant", CO_INT, ARM),
    N("n_reverses_without_cgas", CO_INT, ARM), LBL_TOP),
  config    = FIG_CFG,
  # Square canvas: coord_fixed keeps the identity line at 45°, which is the whole point.
  width = 9.5, height = 9.5
)

# -----------------------------------------------------------------------------
# 4. PANEL B — the shared temperature axis against the cGAS-dependence axis
# -----------------------------------------------------------------------------
stats_box_b <- paste(
  sprintf("%s genes, one dot each", N("n_genes", "all")),
  sprintf("the two axes are near-independent: r = %.2f",
          S("pearson_r", paste(CO_SHARE, CO_INT, sep = "_vs_"))),
  sprintf("%s of %s heat-responsive genes differ by genotype",
          N("n_heat_responsive_cgas_dependent", CO_INT),
          N("n_heat_responsive", paste(CO_WT, CO_KO, sep = "_or_"))),
  sprintf("%s of those %s reverse direction (triangles)",
          N("n_reverses_without_cgas", CO_INT, ARM), N("n_sig", CO_INT)),
  sep = "\n")

p_b <- ggplot(mapping = aes(x = .data[[paste0("logFC_", CO_SHARE)]],
                            y = .data[[paste0("logFC_", CO_INT)]])) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = BG) +
  geom_point(data = gw_bg, aes(colour = arm, shape = arm), size = PT_BG, alpha = 0.35) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, colour = AX) +
  geom_point(data = gw_same, aes(colour = arm, shape = arm), size = PT_HL) +
  geom_point(data = gw_rev,  aes(colour = arm, shape = arm), size = PT_HL * 1.25) +
  ggrepel::geom_text_repel(
    data = gw_lab, aes(label = gene_symbol), colour = HL, size = LBL,
    seed = SEED, max.overlaps = Inf, force = 8, force_pull = 0.4,
    box.padding = 0.9, point.padding = 0.6, min.segment.length = 0,
    segment.colour = HL, segment.linewidth = 0.4,
    # White halo: gene names stay readable where a leader line or a neighbouring
    # dot passes under them, without a filled box hiding data.
    bg.colour = HALO, bg.r = 0.14, show.legend = FALSE) +
  # Placed on separate rows of the empty lower half so no two blocks can collide.
  annotate("text", x = 0.98 * LIM_BX, y = -0.46 * LIM_BY, hjust = 1, vjust = 1,
           size = CUE, label = stats_box_b) +
  annotate("text", x = 0.97 * LIM_BX, y = 0.97 * LIM_BY, hjust = 1, vjust = 1,
           size = CUE, colour = HL,
           label = paste("above the line: heat response stronger in WT",
                         "— the cGAS-dependent arm, with no mirror below",
                         sep = "\n")) +
  annotate("text", x = -0.98 * LIM_BX, y = -0.30 * LIM_BY, hjust = 0, vjust = 1,
           size = CUE, colour = AX,
           label = paste("the dense band on the dashed line:",
                         "heat response shared by both genotypes",
                         sep = "\n")) +
  scale_colour_manual(values = FILLS,  breaks = KEYS, name = NULL, drop = FALSE) +
  scale_shape_manual(values  = SHAPES, breaks = KEYS, name = NULL, drop = FALSE) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = PT_HL * 2, alpha = 1)),
         shape  = guide_legend(ncol = 1, override.aes = list(size = PT_HL * 2, alpha = 1))) +
  coord_cartesian(xlim = c(-LIM_BX, LIM_BX), ylim = c(-LIM_BY, LIM_BY)) +
  labs(
    title    = paste("Two separate threads: a shared temperature response",
                      "and a small cGAS-dependent arm", sep = "\n"),
    subtitle = wrap_at(sprintf(
      paste0("Each dot is one gene, in log2 fold change. x = its heat response averaged ",
             "over both genotypes, y = how much stronger that response is in WT than in ",
             "cGAS-KO. Typical sizes: %.2f on x against %.2f on y, and the vertical axis is ",
             "drawn %.1f× the horizontal so the small arm is visible at all."),
      S("median_abs_logfc", CO_WT), S("median_abs_logfc", CO_INT),
      S("y_expansion", paste(CO_SHARE, CO_INT, sep = "_over_")))),
    x = sprintf("Shared heat response — %s", contrast_label(CO_SHARE)),
    y = "WT minus cGAS-KO heat response",
    caption = wrap_at(paste(ARM_SENTENCE, REVERSE_SENTENCE, FLOOR_SENTENCE))
  ) +
  project_theme(config = FIG_CFG) +
  ggplot2::theme(legend.position = "bottom")

p_b <- rasterize_axes(p_b, config = FIG_CFG)

tbl_b <- gw %>%
  dplyr::select(gene_symbol, ensembl,
                dplyr::all_of(c(paste0("logFC_", CO_SHARE), paste0("adjP_", CO_SHARE),
                                paste0("logFC_", CO_INT), paste0("adjP_", CO_INT))),
                heat_responsive, cgas_dependent, interaction_rank)

save_overview(
  plot      = p_b,
  stage     = STAGE,
  name      = "heat_response_shared_vs_cgas_arm",
  table     = tbl_b,
  finding   = sprintf(
    paste0("Splitting the same fold changes into a shared temperature axis and a ",
           "cGAS-dependence axis separates two threads that barely overlap: the two are ",
           "near-independent (r = %.2f), the shared response is the larger of the two ",
           "(median |log2FC| %.2f against %.2f), and %s of the %s heat-responsive genes also ",
           "differ by genotype. Every one of them sits above the line and none below. %s ",
           "Claim tier: L3."),  # %s = the reversal reading
    S("pearson_r", paste(CO_SHARE, CO_INT, sep = "_vs_")),
    S("median_abs_logfc", CO_WT), S("median_abs_logfc", CO_INT),
    N("n_heat_responsive_cgas_dependent", CO_INT),
    N("n_heat_responsive", paste(CO_WT, CO_KO, sep = "_or_")),
    REVERSE_SENTENCE
  ),
  script    = SCRIPT,
  fn        = "geom_point + geom_hline (shared-axis versus cGAS-dependence-axis scatter)",
  config_kv = CFG_KV,
  input     = "03_results/03_de/tables/_overview/cgas_dependence_wide.csv + cgas_dependence_stats.csv",
  how_to_read = sprintf(
    paste0("One dot per gene. x = heat response pooled over genotypes, y = the WT heat ",
           "response minus the cGAS-KO heat response, both log2. The dashed line is y = 0: a ",
           "gene on it responded to heat identically in the two genotypes. Grey circles = no ",
           "detectable difference at n=5. Vermillion circles = differs by genotype ",
           "(adj.P < %.2g); triangles = also reverses sign, up with heat in WT and down without ",
           "cGAS (%s of the %s, %s significant in WT on their own). Top %d labelled, reversing ",
           "genes first. The y axis is drawn %.1f× the x axis, so vertical spread is magnified ",
           "on purpose. Gate: FDR only. Claim tier: L3. PROVISIONAL sample labels."),
    FDR, N("n_reverses_without_cgas", CO_INT, ARM), N("n_sig", CO_INT),
    N("n_reverses_wt_significant", CO_INT, ARM), LBL_TOP,
    S("y_expansion", paste(CO_SHARE, CO_INT, sep = "_over_"))),
  config    = FIG_CFG,
  wide = TRUE, width = 11, height = 8.5
)

# -----------------------------------------------------------------------------
# 5. Captions for the sibling artifacts save_overview does not key on
#    (it captions the .png; the vector companion and the plotted source table
#    each need their own section so every file under 03_results/ is captioned).
# -----------------------------------------------------------------------------
panel_meta <- list(
  list(stem = "heat_response_wt_vs_ko",
       what = "the WT heat response plotted against the cGAS-KO heat response",
       read = paste0("Grey = no detectable difference between genotypes at n=5, vermillion = ",
                     "passes the genotype comparison of the heat response. On the dashed ",
                     "identity line = same response with and without cGAS; distance from it = ",
                     "the WT effect minus the cGAS-KO effect."),
       fn   = "geom_point + geom_abline (effect-versus-effect scatter, equal axes)"),
  list(stem = "heat_response_shared_vs_cgas_arm",
       what = "the shared temperature axis plotted against the cGAS-dependence axis",
       read = paste0("Grey = no detectable difference between genotypes at n=5, vermillion = ",
                     "passes the genotype comparison of the heat response. y = 0 marks an ",
                     "identical response in both genotypes, and every highlighted gene lies ",
                     "above it."),
       fn   = "geom_point + geom_hline (shared-axis versus cGAS-dependence-axis scatter)")
)

for (m in panel_meta) {
  write_caption(
    stage    = STAGE,
    filename = file.path("figures", OVD, paste0(m$stem, ".pdf")),
    finding  = sprintf(paste0("Vector companion to figures/%s/%s.png, %s. Same plot object, ",
                              "with the dense dot layer rasterised so the file stays small and ",
                              "the axes and text stay editable. Claim tier: L3."),
                       OVD, m$stem, m$what),
    script      = SCRIPT, fn = m$fn, config_kv = CFG_KV,
    input       = "03_results/03_de/tables/_overview/cgas_dependence_wide.csv + cgas_dependence_stats.csv",
    how_to_read = sprintf("%s Gate: adj.P < %.2g, no fold-change cut-off. %s Claim tier: L3.",
                          m$read, FDR, ARM_SENTENCE),
    config = FIG_CFG
  )
  write_caption(
    stage    = STAGE,
    filename = file.path("tables", OVD, paste0(m$stem, ".csv")),
    finding  = sprintf(paste0("Plotted values behind figures/%s/%s.png: one row per gene, the ",
                              "columns the panel draws, sliced from cgas_dependence_wide.csv. ",
                              "%s genes. Claim tier: L3."),
                       OVD, m$stem, N("n_genes", "all")),
    script      = SCRIPT, fn = "save_overview (writes the figure's same-stem source table)",
    config_kv   = CFG_KV,
    input       = "03_results/03_de/tables/_overview/cgas_dependence_wide.csv",
    how_to_read = sprintf(paste0("logFC_ / adjP_ columns are the log2 fold change and BH-adjusted ",
                                 "p-value per contrast. cgas_dependent = adj.P < %.2g on the genotype ",
                                 "comparison of the heat response (FDR only, no fold-change cut-off). ",
                                 "interaction_rank orders that arm by evidence and is NA elsewhere; ",
                                 "the panel labels its top %d. Claim tier: L3."), FDR, LBL_TOP),
    config = FIG_CFG
  )
}

# -----------------------------------------------------------------------------
# 6. Read-only captions for the compute sibling's two tables
# -----------------------------------------------------------------------------
write_caption(
  stage    = STAGE,
  filename = file.path("tables", OVD, "cgas_dependence_wide.csv"),
  finding  = sprintf(
    paste0("One row per gene with the four heat-relevant contrasts side by side, joined on ",
           "ensembl: the WT and cGAS-KO heat responses, their shared average, and the genotype ",
           "comparison of the heat response. %s genes, all tested in every contrast. This is the ",
           "single source for both scatter panels and the table in which the panels' geometric ",
           "claim is checkable: wt_minus_ko equals logFC_Interaction gene by gene. Claim tier: L3."),
    N("n_genes", "all")),
  script      = "02_analysis/scripts/02b_cgas_dependence_geometry.R",
  fn          = "inner_join on ensembl (read-only here; produced by the compute sibling)",
  config_kv   = CFG_KV,
  input       = "03_results/objects/02_de_results.rds",
  how_to_read = sprintf(
    paste0("Columns: ensembl, gene_symbol, then logFC_/adjP_/sig_ per contrast (%s, %s, %s, %s). ",
           "sig_ = adj.P.Val < %.2g, no fold-change cut-off. heat_responsive = significant for ",
           "heat in at least one genotype. cgas_dependent = passes the genotype comparison of the ",
           "heat response. up_with_heat_in_wt and down_with_heat_in_ko are logFC signs, ",
           "reverses_without_cgas is both at once, and arm_class folds these into the three ",
           "classes the panels draw. wt_minus_ko = the WT effect minus the cGAS-KO effect. ",
           "interaction_rank orders the arm by evidence and label_rank puts reversing genes ",
           "first; both are NA off the arm. Positive logFC = higher at 39 °C, except the ",
           "genotype comparison, where positive = a larger heat response in WT. Claim tier: L3."),
    CO_WT, CO_KO, CO_INT, CO_SHARE, FDR),
  config = FIG_CFG
)

write_caption(
  stage    = STAGE,
  filename = file.path("tables", OVD, "cgas_dependence_stats.csv"),
  finding  = paste0(
    "The scalars both scatter panels print, so no number on a figure is computed at draw ",
    "time. Correlations and the regression slope between the two per-genotype heat responses ",
    "are reported twice, over all genes and over the heat-responsive ones only, because ",
    "agreement measured on the whole universe could be carried by unchanged genes sitting at ",
    "the origin. The restricted scope is the one the panels quote, and it is the higher of the ",
    "two. Also carries per-contrast counts, the anatomy of the cGAS-dependent arm, typical ",
    "effect sizes, the panel axis ranges, and the identity residual. Claim tier: L3."),
  script      = "02_analysis/scripts/02b_cgas_dependence_geometry.R",
  fn          = "cor + lm + tallies (read-only here; produced by the compute sibling)",
  config_kv   = CFG_KV,
  input       = "03_results/objects/02_de_results.rds",
  how_to_read = paste0(
    "Long format: metric, scope, subset, value, note, so a lookup is (metric, scope, subset). ",
    "subset is all_genes, heat_responsive or cgas_dependent_arm and must always be read -- ",
    "pearson_r appears under two of them. The arm rows count how many of the significant genes ",
    "fall with heat in cGAS-KO, how many reverse sign, and how many of those hold up ",
    "individually in each genotype. axis_lim gives each panel's symmetric half-range and ",
    "y_expansion the ratio between them. max_abs_identity_residual is 0 to numerical precision, ",
    "which is what licenses reading distance from the identity line as the genotype difference. ",
    "Every count uses adj.P.Val < de_fdr with no fold-change cut-off. Claim tier: L3."),
  config = FIG_CFG
)

# -----------------------------------------------------------------------------
# 7. Structural asserts (LOUD)
# -----------------------------------------------------------------------------
fig_ovw <- overview_path(STAGE, "figures", FIG_CFG)
expected <- c(file.path(fig_ovw, paste0("heat_response_wt_vs_ko.", c("pdf", "png"))),
              file.path(fig_ovw, paste0("heat_response_shared_vs_cgas_arm.", c("pdf", "png"))),
              file.path(TBL_OVW, "heat_response_wt_vs_ko.csv"),
              file.path(TBL_OVW, "heat_response_shared_vs_cgas_arm.csv"))
missing <- expected[!file.exists(expected)]
if (length(missing) > 0L)
  stop("[02b_viz] expected artifact(s) missing:\n  ", paste(missing, collapse = "\n  "))

message("02b_cgas_dependence_geometry_viz.R complete.")
message("  figures/_overview/heat_response_wt_vs_ko.{pdf,png}")
message("  figures/_overview/heat_response_shared_vs_cgas_arm.{pdf,png}")
