#!/usr/bin/env Rscript
# 02b_cgas_dependence_geometry_viz.R -- VIZ
# =============================================================================
# One panel that carries the cGAS-dependence argument to a reader who has never
# seen a linear-model contrast (stage 03_de).
#
# The argument
#   Warming iTregs to 39 C changes thousands of genes, and it changes almost the
#   same genes by almost the same amount whether or not cGAS is present. A small,
#   one-sided arm behaves differently: 23 genes respond more strongly to heat in
#   WT than in cGAS-KO, none the other way round.
#
#   The panel puts the WT heat response on x and the cGAS-KO heat response on y.
#   Both axes are the SAME familiar quantity measured in two genotypes, so no
#   modelling vocabulary is needed: a gene on the diagonal responded to heat the
#   same way in both, and distance from the diagonal is exactly how much removing
#   cGAS changed its response (asserted as an identity by the compute sibling).
#
# Why there is no second, rotated panel
#   A companion panel used to plot Temp_main (x) against Interaction (y). Because
#   Temp_main = 1/2(WT_heat + KO_heat) and Interaction = WT_heat - KO_heat, that
#   view is a 45-degree ROTATION of this one, not a second measurement, and its
#   two headline numbers are artifacts of the rotation rather than results:
#   cov(x,y) = 1/2(var(WT_heat) - var(KO_heat)), so its r = -0.08 is fixed by the
#   two genotypes having near-equal response variance -- the same fact this panel
#   reports as slope 0.99 / r = 0.95 -- and "the shared axis carries the larger
#   effects" restates that a cloud hugging the identity line spreads further along
#   it than across it. Reporting either as independent threads over-read the
#   geometry, so the rotated view was dropped rather than re-captioned.
#
# Encoding (value first, hue second, shape third)
#   The three classes are ordered by LIGHTNESS so they survive a greyscale print
#   and deuteranopia/protanopia: pale grey cloud -> vermillion -> near-black, each
#   with its own glyph (small dot -> ringed circle -> ringed triangle). Hue only
#   reinforces a separation the value ladder already carries.
#
# The stringent gate is drawn as GEOMETRY, not as a fourth colour
#   Two subsets of the arm both number nine and are NOT the same genes: the
#   |logFC| >= de_logfc set, and the set that reverses sign without cGAS. They
#   share four members. Drawing the gate as dotted lines parallel to the identity
#   line makes membership exact and readable off the page: highlighted glyph
#   beyond the line = gate, triangle = reverses, so beyond+triangle = both.
#   The gate is NOT "the genes that travel to human". The frozen
#   03_results/human_projection/ contract exports this arm at BOTH gates -- the
#   nine at fdr_logfc and all 23 at fdr_only, the sensitivity read a 1 df term at
#   n=5 needs -- so the canvas states the geometry and leaves the export to the
#   README and to human_projection/manifest.csv.
#
# Role: DRAWS ONLY. Every number on these panels -- correlations, regression
#   slope, counts, class membership, the gate offset, axis limits, which genes get
#   labelled -- is read from the tables emitted by 02b_cgas_dependence_geometry.R.
#   No statistic, join, threshold or count is computed here.
#
# Inputs (read only)
#   03_results/03_de/tables/_overview/cgas_dependence_wide.csv
#   03_results/03_de/tables/_overview/cgas_dependence_stats.csv
#
# Outputs
#   03_results/03_de/figures/_overview/heat_response_wt_vs_ko.{pdf,png}
#   03_results/03_de/tables/_overview/heat_response_wt_vs_ko.csv
#   03_results/03_de/README.md (captions, idempotent)
#
# Honest framing
#   The interaction is LABELLED by what it tests -- cGAS-dependence of the heat
#   response -- never by a result, and its claim floors at L3. It is a 1 df term
#   at n=5/group, the least-powered contrast in the design, so its 23 genes are a
#   FLOOR on the cGAS-dependent arm. A gene on the diagonal has no detectable
#   cGAS-dependence at n=5, which is not a claim of cGAS-independence. The arm is
#   also one-sided (23 up, 0 down); nothing here may be drawn symmetrically.
#   One-sidedness is a property of the SIGNIFICANT arm ONLY. The pale cloud
#   straddles the diagonal in both directions -- 9,541 of the 19,679 genes have a
#   nominally weaker heat response in WT, and one of them (the pseudogene
#   Spcs2-ps, adj.P = 1) clears the gate offset on that side. So every "none the
#   other way" sentence, on the canvas or in a caption, must name the arm it is
#   scoped to; unscoped, it is false and the picture contradicts it.
#   The reader-facing prose lives in the README caption sections, not on the canvas.
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

## DECLARED PALETTE — the only colour constants in this script, all from FIG_CFG$colors.
## The three classes are ordered by VALUE, not hue, so the panel separates in pure
## greyscale as well as in colour. Approximate WCAG relative luminance in brackets:
##   VAL_LO  pale grey   [Y ~ 0.60]  the 19,656-gene cloud
##   VAL_MID vermillion  [Y ~ 0.22]  differs between genotypes, same direction
##   VAL_HI  black       [Y ~ 0.00]  differs AND reverses direction without cGAS
## Successive greyscale contrast ratios are 2.6 and 5.4, so the ladder is legible
## before hue is consulted; vermillion-versus-neutral is the deuteranopia-safe cue.
OI      <- FIG_CFG$colors$okabe_ito %||% list()
VAL_LO  <- "grey80"
VAL_MID <- OI$vermillion %||% "#D55E00"
VAL_HI  <- OI$black      %||% "#000000"
RING    <- "grey20"     # outline on the two highlighted glyphs (lifts them off the cloud)
AX      <- "grey35"     # reference lines: identity + gate
GRID    <- "grey72"     # zero rules
TXT     <- "grey15"     # in-panel reader text
LBL_INK <- "grey20"     # gene labels
HALO    <- "white"      # gene-label outline

CO_WT    <- "WT_heat"
CO_KO    <- "KO_heat"
CO_INT   <- "Interaction"
CO_SHARE <- "Temp_main"

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
LIM_A   <- S("axis_lim", WT_KO)          # symmetric, equal on both axes
GATE    <- S("de_logfc", "gate")         # the |logFC| cut-off drawn as dotted lines

CFG_KV <- sprintf(
  paste0("thresholds.de_fdr=%.2g; thresholds.de_logfc=%g; figures.volcano_label_top=%d; ",
         "figures.point_size=%.1f; figures.cue_size=%.1f; colors.okabe_ito"),
  FDR, GATE, LBL_TOP, PT, CUE)

# Point cloud vs highlighted arm: one is 19,679 points, the other 23. All sizes
# derive from figures.point_size so the family scales with the contract.
PT_BG <- PT * 0.28
PT_HL <- PT * 0.95

# -----------------------------------------------------------------------------
# 3. Class keying — legend keys name the GEOMETRY the reader can see
# -----------------------------------------------------------------------------
# The old keys ("heat response differs between genotypes", "reverses direction")
# described the test, not the picture: a reader saw two spatially separate clouds
# under one generic label. These keys say where each class sits and what that
# means biologically, in words a wet-lab reader already owns. arm_class comes from
# the compute sibling; the viz only attaches display labels.
KEY_BG   <- "no detectable difference between genotypes at n=5"
KEY_SAME <- "falls with heat in both genotypes, and further without cGAS"
KEY_REV  <- "rises with heat in WT, falls without cGAS"
KEYS     <- c(KEY_BG, KEY_SAME, KEY_REV)
CLASS_KEY <- c(no_detectable_difference = KEY_BG,
               differs_same_direction   = KEY_SAME,
               differs_and_reverses     = KEY_REV)
# Ordered light -> dark, one glyph each: small dot, ringed circle, ringed triangle.
SHAPES  <- setNames(c(16, 21, 24),                KEYS)
FILLS   <- setNames(c(VAL_LO, VAL_MID, VAL_HI),   KEYS)
COLOURS <- setNames(c(VAL_LO, RING, RING),        KEYS)

gw <- gw %>%
  dplyr::mutate(arm = factor(unname(CLASS_KEY[arm_class]), levels = KEYS))

gw_bg   <- gw %>% dplyr::filter(arm_class == "no_detectable_difference")
gw_same <- gw %>% dplyr::filter(arm_class == "differs_same_direction")
gw_rev  <- gw %>% dplyr::filter(arm_class == "differs_and_reverses")

# Labels: the compute sibling's evidence ranking capped at the contract's label
# budget, UNION the stringent-gate set. The gate set is the narrated one -- the genes
# the frozen mouse->human signature carries -- and it is itself smaller than the cap,
# so neither label source exceeds figures.volcano_label_top. Naming both subsets is
# the point: the reader must be able to check which of the two nines a gene is in.
# Membership is READ from in_stringent_gate / label_rank, never derived here.
# Gene names of gate members are set BOLD. Four of the nine sit within a hundredth of
# a log2 unit of the gate line -- a hard threshold cutting a continuum -- so the line
# alone cannot resolve them at glyph size. Bold is exact, costs no colour and no legend
# row, and survives greyscale. gate_face is a display attribute of a read column.
gw_lab <- gw %>%
  dplyr::filter((!is.na(label_rank) & label_rank <= LBL_TOP) | in_stringent_gate) %>%
  dplyr::mutate(gate_face = ifelse(in_stringent_gate, "bold", "plain"))

message(sprintf("[02b_viz] %s genes; arm %s (%s reverse, %s clear the |logFC| >= %g gate, %s in both)",
                N("n_genes", "all"), N("n_sig", CO_INT),
                N("n_reverses_without_cgas", CO_INT, ARM),
                N("n_stringent_gate", CO_INT, ARM), GATE,
                N("n_stringent_and_reverses", CO_INT, ARM)))

# -----------------------------------------------------------------------------
# 4. In-panel text — short blocks only; the paragraphs live in the README
# -----------------------------------------------------------------------------
# Wrap widths are fractions of the contract's caption column, chosen so each block
# fits the VERIFIED-EMPTY region it is placed in (see the placement comments below).
# A block that overruns its region lands on data, which is the failure mode here.
WRAP     <- as.integer(FIG_CFG$figures$caption_wrap_column %||% 70)
WRAP_BOX <- as.integer(round(WRAP * 0.66))   # the two corner blocks

#' Wrap one long string for a title/subtitle.
wrap_at <- function(txt, width = WRAP) paste(strwrap(txt, width = width), collapse = "\n")

# The arm's symbols read as interferon-stimulated, which is a look, not a result. State
# the CURATED composition instead: membership of MSigDB Hallmark interferon-alpha/gamma,
# with the unassigned remainder reported and the reason it is not a negative.
COMPOSITION_SENTENCE <- sprintf(
  paste0("Curated composition: %s of the %s are Hallmark interferon-alpha or -gamma members ",
         "and %s are unassigned, several of them mouse-specific paralogs the human-derived ",
         "sets omit."),
  N("n_hallmark_ifn", CO_INT, ARM), N("n_sig", CO_INT),
  N("n_hallmark_ifn_unassigned", CO_INT, ARM))
#' Assemble an in-panel block. Each element is wrapped SEPARATELY, so the line count
#' is predictable: short elements are one line each and long ones break cleanly.
box_text <- function(..., width = WRAP_BOX)
  paste(unlist(lapply(c(...), strwrap, width = width)), collapse = "\n")
#' Guard the block against the region it was placed in. A block that outgrows its
#' budget lands on data or runs off the canvas, which is this figure's failure mode,
#' so it fails the run loudly instead of shipping a collision.
fits <- function(txt, max_lines, max_chars, where) {
  ln <- strsplit(txt, "\n", fixed = TRUE)[[1]]
  if (length(ln) > max_lines || max(nchar(ln)) > max_chars)
    stop(sprintf("[02b_viz] the %s block is %d lines x %d chars; budget is %d x %d.",
                 where, length(ln), max(nchar(ln)), max_lines, max_chars))
  invisible(txt)
}

# -----------------------------------------------------------------------------
# 5. THE PANEL — the heat response in WT against the heat response in cGAS-KO
# -----------------------------------------------------------------------------
# The counts block is NUMBERS, not sentences: one row per quantity the geometry is built
# from, so a reader can check the picture against the arithmetic without parsing prose. It
# carries counts ONLY. Δ and the significance cut-off are defined in the subtitle, and the
# agreement statistics (r, slope) are left to the README caption — the dots piling onto the
# identity line already show the agreement, so printing r beside it spends canvas on a
# number the picture makes. The per-direction split uses bare arrow glyphs rather than
# direction_cue()'s "arrow + word" form because the two counts sit adjacent on one line and
# disambiguate each other (bare-glyph precedent: 15_coresh_viz.R facet labels).
box_a_stats <- fits(box_text(
  sprintf("%s genes", N("n_genes", "all")),
  sprintf("ΔWT: %s sig  ↑ %s  ↓ %s",
          N("n_sig", CO_WT), N("n_up", CO_WT), N("n_down", CO_WT)),
  sprintf("ΔcGAS-KO: %s sig  ↑ %s  ↓ %s",
          N("n_sig", CO_KO), N("n_up", CO_KO), N("n_down", CO_KO)),
  sprintf("ΔWT − ΔcGAS-KO: %s sig  ↑ %s  ↓ %s",
          N("n_sig", CO_INT), N("n_up", CO_INT), N("n_down", CO_INT))),
  max_lines = 5, max_chars = WRAP_BOX - 1, where = "upper-left counts")

# The key states the two reference lines as EQUATIONS on the same Δ vocabulary as the
# counts block, so the geometry is self-defining. Both claims are scoped to the
# HIGHLIGHTED arm: the pale cloud straddles the dashed line in both directions, so an
# unscoped "none below" would be contradicted by the picture it annotates. Where the gated
# genes go downstream is deliberately NOT here — the frozen export carries the arm at two
# gates, and one canvas line cannot say that without misleading (see the header note).
box_a_read <- fits(box_text(
  "dashed line: ΔWT = ΔcGAS-KO",
  sprintf("dotted lines: ΔWT − ΔcGAS-KO = ±%g", GATE),
  sprintf("all %s highlighted genes lie below the dashed line; the %s in bold also clear ±%g",
          N("n_sig", CO_INT), N("n_stringent_gate", CO_INT, ARM), GATE)),
  max_lines = 5, max_chars = WRAP_BOX - 1, where = "lower-right key")

p_a <- ggplot(mapping = aes(x = .data[[paste0("logFC_", CO_WT)]],
                            y = .data[[paste0("logFC_", CO_KO)]],
                            fill = arm, colour = arm, shape = arm)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = GRID) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = GRID) +
  geom_point(data = gw_bg, size = PT_BG, alpha = 0.35) +
  # Identity line, then the stringent gate as two lines PARALLEL to it: the gate is
  # |WT - cGAS-KO| >= de_logfc, so it is exactly an offset of the diagonal. Both
  # offsets are drawn; nothing highlighted lies beyond the upper one, which is the
  # one-sidedness of the arm shown rather than asserted.
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.6, colour = AX) +
  geom_abline(slope = 1, intercept = c(-GATE, GATE), linetype = "dotted",
              linewidth = 0.5, colour = AX) +
  geom_point(data = gw_same, size = PT_HL) +
  geom_point(data = gw_rev,  size = PT_HL * 1.25) +
  ggrepel::geom_text_repel(
    data = gw_lab, aes(label = gene_symbol, fontface = gate_face),
    colour = LBL_INK, size = LBL,
    seed = SEED, max.overlaps = Inf, force = 20, force_pull = 0.12,
    box.padding = 1.2, point.padding = 0.75, min.segment.length = 0,
    segment.colour = LBL_INK, segment.size = 0.4,
    # Confine labels to the region left free by the two annotation blocks (verified
    # empty of both blocks and of data at these bounds), so no label can collide.
    xlim = c(-0.97 * LIM_A, 0.95 * LIM_A), ylim = c(-0.52 * LIM_A, 0.35 * LIM_A),
    # White halo: gene names stay readable where a leader line or a neighbouring
    # dot passes under them, without a filled box hiding data.
    bg.colour = HALO, bg.r = 0.14, show.legend = FALSE) +
  # Both blocks sit in regions with zero data points at these bounds: the upper-left
  # triangle above the cloud, and the lower-right strip under it.
  annotate("text", x = -0.97 * LIM_A, y = 0.97 * LIM_A, hjust = 0, vjust = 1,
           size = CUE, colour = TXT, label = box_a_stats) +
  annotate("text", x = 0.97 * LIM_A, y = -0.97 * LIM_A, hjust = 1, vjust = 0,
           size = CUE, colour = TXT, label = box_a_read) +
  scale_fill_manual(values   = FILLS,   breaks = KEYS, name = NULL, drop = FALSE) +
  scale_colour_manual(values = COLOURS, breaks = KEYS, name = NULL, drop = FALSE) +
  scale_shape_manual(values  = SHAPES,  breaks = KEYS, name = NULL, drop = FALSE) +
  # One merged key per class. fill carries the size/alpha override; colour and shape
  # only have to agree on the layout for ggplot to fold all three into that one key.
  guides(fill   = guide_legend(ncol = 1, override.aes = list(size = PT_HL * 2, alpha = 1)),
         colour = guide_legend(ncol = 1),
         shape  = guide_legend(ncol = 1)) +
  scale_discrete_identity(aesthetics = "fontface") +   # bold = inside the gate
  coord_fixed(ratio = 1, xlim = c(-LIM_A, LIM_A), ylim = c(-LIM_A, LIM_A)) +
  labs(
    title    = paste("Warming to 39 °C changes the same genes",
                      "with and without cGAS — apart from 23", sep = "\n"),
    # Three blocks, hard-separated: the question this panel opens the sequence with
    # and the answer it gives, then the axes, then the notation the counts block
    # uses. Δ is bound to the quantity the axis line just defined rather than
    # re-stating it, so the counts block can be pure numbers.
    subtitle = paste(
      wrap_at(sprintf(paste0(
        "Question — does 39 °C produce a coherent transcriptional response, and does that ",
        "response require cGAS? Answer — coherent, and largely shared: %s of %s genes ",
        "separate the two genotypes, and every one of them lies on the same side of the line."),
        N("n_sig", CO_INT), N("n_genes", "all"))),
      wrap_at(paste0(
        "Each dot is one gene. Both axes are its log2 fold change at 39 versus 37 °C, ",
        "measured in WT (x) and in cGAS-KO (y).")),
      wrap_at(sprintf("Δ = that log2 fold change; sig = adj.P < %.2g.", FDR)),
      sep = "\n"),
    x = sprintf("%s — log2 fold change", contrast_label(CO_WT)),
    y = sprintf("%s — log2 fold change", contrast_label(CO_KO))
  ) +
  project_theme(config = FIG_CFG) +
  ggplot2::theme(legend.position = "bottom")

p_a <- rasterize_axes(p_a, config = FIG_CFG)   # dense dot layer -> raster inside the vector PDF

tbl_a <- gw %>%
  dplyr::select(gene_symbol, ensembl,
                dplyr::all_of(c(paste0("logFC_", CO_WT), paste0("adjP_", CO_WT),
                                paste0("logFC_", CO_KO), paste0("adjP_", CO_KO),
                                paste0("logFC_", CO_INT), paste0("adjP_", CO_INT))),
                cgas_dependent, in_stringent_gate, reverses_without_cgas,
                arm_class, interaction_rank)

save_overview(
  plot      = p_a,
  stage     = STAGE,
  name      = "heat_response_wt_vs_ko",
  table     = tbl_a,
  finding   = sprintf(
    paste0("Warming iTregs to 39 °C changes %s genes in WT and %s in cGAS-KO, and it changes ",
           "them in step: among the %s heat-responsive genes the cGAS-KO response is %.2f× the ",
           "WT response (r = %.2f, Spearman %.2f), so the dots pile onto the identity line. A ",
           "one-sided minority behaves differently. All %s genes that pass the genotype ",
           "comparison of the heat response lie below the line, a weaker heat response without ",
           "cGAS, and %s lie above. %s keep that direction in both genotypes and %s flip sign: ",
           "up with heat in WT, down without cGAS. Claim tier: L3."),
    N("n_sig", CO_WT), N("n_sig", CO_KO), N("n_genes", "all", HR),
    S("ols_slope", KO_ON_WT, HR), S("pearson_r", WT_KO, HR), S("spearman_rho", WT_KO, HR),
    N("n_up", CO_INT), N("n_down", CO_INT),
    N("n_differs_same_direction", CO_INT, ARM), N("n_reverses_without_cgas", CO_INT, ARM)
  ),
  script    = SCRIPT,
  fn        = "geom_point + geom_abline (effect-versus-effect scatter, equal axes, gate as parallels)",
  config_kv = CFG_KV,
  input     = "03_results/03_de/tables/_overview/cgas_dependence_wide.csv + cgas_dependence_stats.csv",
  how_to_read = sprintf(
    paste0("One dot per gene: x = log2 fold change at 39 vs 37 °C in WT, y = the same in ",
           "cGAS-KO, equal scales, so the dashed identity line runs at 45°. Glyphs run pale to ",
           "dark: pale dots have no detectable cGAS-dependence at n=5, vermillion circles fall ",
           "with heat in both genotypes and further without cGAS, black triangles rise with heat ",
           "in WT and fall without it. Highlighted genes pass at adj.P < %.2g. Dotted lines lie ",
           "%g log2 unit either side; the %s highlighted genes beyond the lower one, names in ",
           "bold, are the ones that also clear |logFC| >= %g, %s of them triangles. The gate is ",
           "not the set that travels to human: the frozen 03_results/human_projection/ contract ",
           "exports this arm at BOTH gates, these %s at fdr_logfc and all %s at fdr_only, the ",
           "sensitivity read a 1 df term at n=5 needs, so bold marks the stringent core and not ",
           "the whole export (per-gate ortholog counts: human_projection/manifest.csv). Read ",
           "one-sidedness as a statement about the HIGHLIGHTED arm only — the pale cloud ",
           "straddles the dashed line in both directions. Claim tier: L3. PROVISIONAL sample ",
           "labels; n=5/group."),
    FDR, GATE, N("n_stringent_gate", CO_INT, ARM), GATE,
    N("n_stringent_and_reverses", CO_INT, ARM),
    N("n_stringent_gate", CO_INT, ARM), N("n_sig", CO_INT)),
  config    = FIG_CFG,
  # Square canvas: coord_fixed keeps the identity line at 45°, which is the whole point.
  width = 9.5, height = 9.5
)

# -----------------------------------------------------------------------------
# 6. Captions for the sibling artifacts save_overview does not key on
#    (it captions the .png; the vector companion and the plotted source table
#    each need their own section so every file under 03_results/ is captioned).
#    These sections carry the power and membership prose that used to sit on the
#    canvas as a nine-line figure caption.
# -----------------------------------------------------------------------------
FLOOR_SENTENCE <- sprintf(
  paste0("The genotype comparison of the heat response is a single-degree-of-freedom test at ",
         "n=5 per group, the least-powered contrast in this design, so %s genes is a floor on ",
         "the cGAS-dependent arm rather than its full size. A gene on the identity line has no ",
         "detectable cGAS-dependence at n=5, which is weaker than independence."),
  N("n_sig", CO_INT))

MEMBERSHIP_SENTENCE <- sprintf(
  paste0("Two subsets of the arm both number %s and are not the same genes: %s clear |logFC| >= ",
         "%g, the stringent export gate, %s reverse sign without cGAS, and %s belong to ",
         "both. Of the reversing genes %s are individually significant in the WT heat contrast ",
         "and %s in both genotypes, so the flip is a direction read off fitted effects rather ",
         "than %s independent significant flips."),
  N("n_stringent_gate", CO_INT, ARM), N("n_stringent_gate", CO_INT, ARM), GATE,
  N("n_reverses_without_cgas", CO_INT, ARM), N("n_stringent_and_reverses", CO_INT, ARM),
  N("n_reverses_wt_significant", CO_INT, ARM),
  N("n_reverses_both_significant", CO_INT, ARM),
  N("n_reverses_without_cgas", CO_INT, ARM))

panel_meta <- list(
  list(stem = "heat_response_wt_vs_ko",
       what = "the WT heat response plotted against the cGAS-KO heat response",
       fn   = "geom_point + geom_abline (effect-versus-effect scatter, equal axes, gate as parallels)")
)

for (m in panel_meta) {
  write_caption(
    stage    = STAGE,
    filename = file.path("figures", OVD, paste0(m$stem, ".pdf")),
    finding  = sprintf(paste0("Vector companion to figures/%s/%s.png, %s. Same plot object, ",
                              "with the dense dot layer rasterised so the file stays small and ",
                              "the axes and text stay editable. %s Claim tier: L3."),
                       OVD, m$stem, m$what, FLOOR_SENTENCE),
    script      = SCRIPT, fn = m$fn, config_kv = CFG_KV,
    input       = "03_results/03_de/tables/_overview/cgas_dependence_wide.csv + cgas_dependence_stats.csv",
    how_to_read = sprintf(paste0("%s Labels name the top %d genes by evidence plus every gene ",
                                 "inside the gate, so both subsets are readable by name. ",
                                 "Claim tier: L3. PROVISIONAL sample labels; n=5/group."),
                          MEMBERSHIP_SENTENCE, LBL_TOP),
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
                                 "comparison of the heat response. in_stringent_gate adds |logFC| >= ",
                                 "%g and is the panel's dotted-line set; reverses_without_cgas is the ",
                                 "sign flip drawn as triangles; arm_class folds the two into the three ",
                                 "plotted classes. interaction_rank orders the arm by evidence and is ",
                                 "NA elsewhere. hallmark_ifn flags membership of a curated interferon ",
                                 "set and is annotation, never a selection rule. Claim tier: L3."),
                                 FDR, GATE),
    config = FIG_CFG
  )
}

# -----------------------------------------------------------------------------
# 7. Read-only captions for the compute sibling's two tables
# -----------------------------------------------------------------------------
write_caption(
  stage    = STAGE,
  filename = file.path("tables", OVD, "cgas_dependence_wide.csv"),
  finding  = sprintf(
    paste0("One row per gene with the four heat-relevant contrasts side by side, joined on ",
           "ensembl: the WT and cGAS-KO heat responses, their shared average, and the genotype ",
           "comparison of the heat response. %s genes, all tested in every contrast. This is the ",
           "single source for the scatter panel and the table in which the panel's geometric ",
           "claim is checkable: wt_minus_ko equals logFC_Interaction gene by gene. %s ",
           "Claim tier: L3."),
    N("n_genes", "all"), COMPOSITION_SENTENCE),
  script      = "02_analysis/scripts/02b_cgas_dependence_geometry.R",
  fn          = "inner_join on ensembl (read-only here; produced by the compute sibling)",
  config_kv   = CFG_KV,
  input       = "03_results/objects/02_de_results.rds",
  how_to_read = sprintf(
    paste0("Columns: ensembl, gene_symbol, then logFC_/adjP_/sig_ per contrast (%s, %s, %s, %s). ",
           "sig_ = adj.P.Val < %.2g. heat_responsive = significant for heat in at least one ",
           "genotype. cgas_dependent = passes the genotype comparison of the heat response. ",
           "up_with_heat_in_wt and down_with_heat_in_ko are logFC signs, reverses_without_cgas ",
           "is both at once, in_stringent_gate adds |logFC| >= %g, and gate_class and arm_class ",
           "fold these into the classes the panel draws. Temp_main is a scored contrast the panel ",
           "does not plot. wt_minus_ko = the WT effect minus the ",
           "cGAS-KO effect. interaction_rank and label_rank order the arm and are NA off it. ",
           "Positive logFC = higher at 39 °C, except the genotype comparison, where positive = ",
           "a larger heat response in WT. Claim tier: L3."),
    CO_WT, CO_KO, CO_INT, CO_SHARE, FDR, GATE),
  config = FIG_CFG
)

write_caption(
  stage    = STAGE,
  filename = file.path("tables", OVD, "cgas_dependence_stats.csv"),
  finding  = paste0(
    "The scalars the scatter panel prints, so no number on a figure is computed at draw ",
    "time. Correlations and the regression slope between the two per-genotype heat responses ",
    "are reported twice, over all genes and over the heat-responsive ones only, because ",
    "agreement measured on the whole universe could be carried by unchanged genes sitting at ",
    "the origin. The restricted scope is the one the panel quotes, and it is the higher of the ",
    "two. Also carries per-contrast counts with their up/down split, the anatomy of the ",
    "cGAS-dependent arm, typical effect sizes, axis ranges, and the identity residual. ",
    "Claim tier: L3."),
  script      = "02_analysis/scripts/02b_cgas_dependence_geometry.R",
  fn          = "cor + lm + tallies (read-only here; produced by the compute sibling)",
  config_kv   = CFG_KV,
  input       = "03_results/objects/02_de_results.rds",
  how_to_read = paste0(
    "Long format: metric, scope, subset, value, note, so a lookup is (metric, scope, subset). ",
    "subset is all_genes, heat_responsive or cgas_dependent_arm and must always be read -- ",
    "pearson_r appears under two of them. The arm rows count how many of the significant genes ",
    "fall with heat in cGAS-KO, how many reverse sign, how many clear the stringent |logFC| ",
    "gate, and how far the two nine-gene subsets overlap. axis_lim gives a symmetric half-range ",
    "per contrast and de_logfc the offset the panel draws as dotted lines; y_expansion and the ",
    "Temp_main/Interaction axis_lim rows are unused. ",
    "max_abs_identity_residual is 0 to numerical precision, which is what ",
    "licenses reading distance from the identity line as the genotype difference. Claim tier: L3."),
  config = FIG_CFG
)

# -----------------------------------------------------------------------------
# 8. Structural asserts (LOUD)
# -----------------------------------------------------------------------------
fig_ovw <- overview_path(STAGE, "figures", FIG_CFG)
expected <- c(file.path(fig_ovw, paste0("heat_response_wt_vs_ko.", c("pdf", "png"))),
              file.path(TBL_OVW, "heat_response_wt_vs_ko.csv"))
missing <- expected[!file.exists(expected)]
if (length(missing) > 0L)
  stop("[02b_viz] expected artifact(s) missing:\n  ", paste(missing, collapse = "\n  "))

# The rotated companion panel was dropped (see the header note). Fail loudly if a stale
# artifact is still on disk, so the README and the results tree cannot disagree.
stale <- c(file.path(fig_ovw, paste0("heat_response_shared_vs_cgas_arm.", c("pdf", "png"))),
           file.path(TBL_OVW, "heat_response_shared_vs_cgas_arm.csv"))
present <- stale[file.exists(stale)]
if (length(present) > 0L)
  stop("[02b_viz] stale artifact(s) from the removed rotated panel:\n  ",
       paste(present, collapse = "\n  "))

message("02b_cgas_dependence_geometry_viz.R complete.")
message("  figures/_overview/heat_response_wt_vs_ko.{pdf,png}")
