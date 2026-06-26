# 16_synthesis_viz.R — VIZ
# =============================================================================
# CAPSTONE headline figure (stage 07_synthesis) — the TWO-ARMS panel.
#
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   07_synthesis
#
# ROLE: VIZ ONLY. Reads two_arms_summary.csv / 16_synthesis.rds (written by
#   16_synthesis.R); recomputes NOTHING (no DE, no fgsea, no decoupleR). Figures
#   ONLY via the figure-style contract: project_theme(config=FIG_CFG),
#   save_overview() (dual print+screen + sibling table + caption). No inline
#   ggsave()/theme()/raw hex; colors come from FIG_CFG$colors.
#
# THE HEADLINE FIGURE (_overview/two_arms_panel):
#   Two STACKED TRACKS, one per arm, read top-to-bottom:
#     TOP    = IFN/ISG arm        (cGAS-DEPENDENT: significant + positive Interaction)
#     BOTTOM = HIF/glycolysis arm (NO DETECTABLE cGAS-dependence: Interaction
#                                  near-zero / NS; rises in BOTH heat arms)
#   x = the headline contrasts (WT_heat | KO_heat | Interaction | Temp_main).
#   Within each track, GLYPH ROWS are the methods (GSEA / PROGENy / TF / DE /
#   GATOM): one tile/point per (method-feature x contrast), colored by signed
#   score, sized/marked by significance. The asymmetry is legible at a glance:
#   the IFN track LIGHTS UP in the Interaction column; the HIF track goes FLAT
#   there while staying lit in WT_heat AND KO_heat.
#
# HOUSE CONSTRAINTS honored in EVERY caption (non-negotiable):
#   * NEVER "cGAS-independent" -> "no detectable cGAS-dependence at n=5".
#   * NEVER crown HIF1a/HIF2a as the driver (HIF is the arm LABEL, not a claim
#     about a TF being #1).
#   * Figure claims FLOOR at L3 (DE/enrichment statistics); mechanism (L7) stays
#     in the reply memo prose, never in a figure title.
#   * PROVISIONAL-sample-label caveat stamped on every panel.
#
# Inputs (read-only):
#   03_results/07_synthesis/tables/two_arms_summary.csv
#   03_results/objects/16_synthesis.rds
#
# Outputs (07_synthesis stage):
#   figures/_overview/two_arms_panel.{print.pdf,screen.png}
#   tables/_overview/two_arms_panel.csv          (source table, via save_overview)
#   03_results/07_synthesis/README.md            (caption, via save_overview)
#
# Run from project root:
#   Rscript 02_analysis/scripts/16_synthesis_viz.R
# =============================================================================

source("02_analysis/helpers/figure_style.R")   # loads figure_helpers.R; sets FIG_CFG

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================================
# CONSTANTS (from FIG_CFG / config; never hardcoded)
# ============================================================================

STAGE  <- "07_synthesis"
SCRIPT <- "02_analysis/scripts/16_synthesis_viz.R"

FDR   <- as.numeric(FIG_CFG$thresholds$gsea_fdr %||% 0.05)
CAP   <- as.numeric(FIG_CFG$figures$nes_cap %||% FIG_CFG$figures$z_clamp %||% 3.5)
# Cap categorical (glyph-row) axis per track to top_n (house anti-clutter rule).
TOP_N <- as.integer(FIG_CFG$figures$top_n %||% 20L)

# Diverging triplet from config (down/neutral/up = blue/white/orange).
NEG <- FIG_CFG$colors$diverging$down    %||% "steelblue4"
MID <- FIG_CFG$colors$diverging$neutral %||% "grey97"
POS <- FIG_CFG$colors$diverging$up      %||% "sienna"

# Headline contrasts (config order). Interaction = the cGAS-dependence test.
HEADLINE_CONTRASTS <- c("WT_heat", "KO_heat", "Interaction", "Temp_main")

# Arm ids (must match 16_synthesis.R).
ARM_IFN <- "IFN_ISG"
ARM_HIF <- "HIF_glycolysis"

# Pretty, line-wrapped facet labels for the two tracks. The arm NAME is a label,
# NOT a mechanistic claim (no "HIF1a is the driver"); the cGAS framing uses the
# house wording.
TRACK_LABELS <- c(
  IFN_ISG        = "IFN / ISG arm\n(cGAS-dependent)",
  HIF_glycolysis = "HIF / glycolysis arm\n(no detectable cGAS-dependence at n=5)")

# Pretty, line-wrapped contrast labels for the x-axis come centrally from
# figure_style.R::CONTRAST_LABELS_SHORT (config-driven; design.contrast_labels_short).
# Used below via scale_x_discrete(labels = CONTRAST_LABELS_SHORT).

# Provisional stamp (matches config.R::provisional_caption()).
PROV_STAMP <- paste0(
  "[PROVISIONAL - inferred sample mapping (Hspa1b/Hsph1 thermometer + Cgas), ",
  "pending collaborator sample sheet; n=5/group]")

# ============================================================================
# 1. GUARD — stop only if BOTH the table and the object are missing; otherwise
#    read whichever is present (table preferred — it is the canonical source).
# ============================================================================

two_arms_csv <- file.path("03_results", STAGE, "tables", "two_arms_summary.csv")
synth_rds    <- file.path("03_results", "objects", "16_synthesis.rds")

if (!file.exists(two_arms_csv) && !file.exists(synth_rds))
  stop("[16_synthesis_viz] neither two_arms_summary.csv nor 16_synthesis.rds found - ",
       "run 16_synthesis.R first.")

synth <- if (file.exists(synth_rds)) readRDS(synth_rds) else NULL

if (file.exists(two_arms_csv)) {
  ev <- readr::read_csv(two_arms_csv, show_col_types = FALSE, progress = FALSE)
  message("[16_synthesis_viz] loaded two_arms_summary.csv (", nrow(ev), " rows).")
} else {
  ev <- as.data.frame(synth$two_arms)
  message("[16_synthesis_viz] two_arms_summary.csv absent; using 16_synthesis.rds$two_arms.")
}

stopifnot(all(c("arm", "method", "feature", "contrast", "score",
                "padj", "direction", "significant") %in% colnames(ev)))

if (nrow(ev) == 0L)
  stop("[16_synthesis_viz] evidence table is EMPTY - the compute arms have not been ",
       "run in-container yet (run 05/06/08/09/10 + 16_synthesis.R), nothing to plot.")

# ============================================================================
# 2. PREP — keep the two cGAS-dependence arms (drop the heat-shock context rows
#    from the headline panel; they are reply-memo context, not an arm), clamp
#    the score for color, order tracks + contrasts, build a stable glyph-row key.
# ============================================================================

ev2 <- ev %>%
  dplyr::filter(arm %in% c(ARM_IFN, ARM_HIF),
                contrast %in% HEADLINE_CONTRASTS) %>%
  dplyr::mutate(
    arm        = factor(arm, levels = c(ARM_IFN, ARM_HIF)),
    contrast   = factor(contrast, levels = HEADLINE_CONTRASTS),
    score_clamp = pmin(pmax(score, -CAP), CAP),
    significant = ifelse(is.na(significant), FALSE, as.logical(significant)),
    # A short, unambiguous glyph-row label: method + feature (the per-method
    # corroborating set/TF/gene). Kept compact so the axis stays legible.
    row_label  = paste0(method, "  ", feature))

if (nrow(ev2) == 0L)
  stop("[16_synthesis_viz] no rows in the two cGAS arms x headline contrasts - nothing to plot.")

# Order glyph rows WITHIN each arm by the Interaction score (descending), so the
# strongest cGAS-dependence signal sits at the top of each track and the eye
# reads the asymmetry top-down. Features absent in Interaction sort last.
row_order <- ev2 %>%
  dplyr::filter(contrast == "Interaction") %>%
  dplyr::group_by(arm, row_label) %>%
  dplyr::summarise(ord_score = mean(score, na.rm = TRUE), .groups = "drop")
# Any row_label not seen in Interaction still needs an order key.
all_rows <- ev2 %>% dplyr::distinct(arm, row_label)
row_order <- dplyr::left_join(all_rows, row_order, by = c("arm", "row_label")) %>%
  dplyr::mutate(ord_score = ifelse(is.na(ord_score), -Inf, ord_score)) %>%
  dplyr::arrange(arm, dplyr::desc(ord_score)) %>%
  dplyr::mutate(row_idx = dplyr::row_number())

# Cap each track to TOP_N glyph rows (house anti-clutter rule), ranked by
# |Interaction score| so the most cGAS-dependence-informative rows survive; ties
# fall back to row_idx. Rank WITHIN arm (so neither track is starved).
row_keep <- row_order %>%
  dplyr::left_join(
    ev2 %>% dplyr::filter(contrast == "Interaction") %>%
      dplyr::group_by(arm, row_label) %>%
      dplyr::summarise(abs_int = mean(abs(score), na.rm = TRUE), .groups = "drop"),
    by = c("arm", "row_label")) %>%
  dplyr::mutate(abs_int = ifelse(is.na(abs_int), -Inf, abs_int)) %>%
  dplyr::group_by(arm) %>%
  dplyr::arrange(dplyr::desc(abs_int), row_idx, .by_group = TRUE) %>%
  dplyr::slice_head(n = TOP_N) %>%
  dplyr::ungroup()
n_dropped <- nrow(row_order) - nrow(row_keep)
if (n_dropped > 0L)
  message(sprintf("[16_synthesis_viz] capping to top_n=%d glyph rows/track; %d rows dropped.",
                  TOP_N, n_dropped))

ev2 <- dplyr::semi_join(ev2, row_keep, by = c("arm", "row_label"))

# Apply the row ordering as a factor on row_label. Because the same feature can
# appear under different methods, row_label is unique per (arm, row_label); we
# order the FACTOR LEVELS by row_idx ascending (top of plot = first level, so we
# reverse for ggplot's bottom-up y-axis to put the strongest signal at the top).
level_order <- row_keep %>%
  dplyr::arrange(row_idx) %>%
  dplyr::pull(row_label) %>%
  unique()
ev2 <- ev2 %>%
  dplyr::left_join(dplyr::select(row_keep, arm, row_label, row_idx),
                   by = c("arm", "row_label")) %>%
  dplyr::mutate(row_label = factor(row_label, levels = rev(level_order)))

# ============================================================================
# 3. THE HEADLINE TWO-ARMS PANEL
#    Geometry: a multi-method glyph grid. Rows = method-feature glyph rows
#    (faceted into the two arm tracks, stacked vertically); columns = the
#    headline contrasts. Fill = signed score (orange up / blue down, clamped);
#    a black ring marks padj < FDR. The Interaction column is the read: IFN
#    track lit, HIF track flat.
# ============================================================================

build_two_arms_panel <- function(df) {
  ggplot(df, aes(x = contrast, y = row_label)) +
    geom_tile(aes(fill = score_clamp), color = "grey92", linewidth = 0.3,
              width = 0.95, height = 0.85) +
    # significance ring (padj < FDR): an open black square outline, color-blind /
    # grayscale safe (does not rely on hue).
    geom_tile(data = dplyr::filter(df, significant),
              fill = NA, color = "black", linewidth = 0.9,
              width = 0.95, height = 0.85) +
    # signed-score label inside each tile (the L3 statistic, made explicit).
    geom_text(aes(label = sprintf("%+.1f", score)),
              size = (FIG_CFG$figures$label_size %||% 4) * 0.8,
              color = "grey15") +
    facet_grid(arm ~ ., scales = "free_y", space = "free_y",
               labeller = labeller(arm = TRACK_LABELS)) +
    scale_fill_gradient2(
      low = NEG, mid = MID, high = POS, midpoint = 0,
      limits = c(-CAP, CAP), oob = scales::squish,
      name = "Signed score\n(orange = up in\nnumerator)") +
    scale_x_discrete(labels = CONTRAST_LABELS_SHORT, position = "top") +
    labs(
      title = "Two-arms cGAS-dependence asymmetry across methods",
      subtitle = paste0(
        "IFN/ISG arm (top) is cGAS-dependent: significant, positive Interaction. ",
        "HIF/glycolysis arm (bottom)\nrises in BOTH heat arms but is flat in the ",
        "Interaction - no detectable cGAS-dependence at n=5.\n", PROV_STAMP),
      x = NULL, y = NULL,
      caption = paste0(
        "Tile = signed score; orange = up in numerator, blue = down. ",
        "Black ring = padj < ", FDR, ". Claim tier: L3 (provisional, n=5/group). ",
        "'Flat Interaction' is NOT proven cGAS-independence.")) +
    project_theme(config = FIG_CFG, legend = TRUE) +
    # Structural facet geometry only (NOT styling): the two arm tracks are stacked
    # vertically via facet_grid(arm ~ .), so the y-strip label must be UN-rotated
    # (ggplot defaults strip.text.y to 270 deg, which would make the two-line track
    # names unreadable), and the stacked tracks need breathing room between them.
    # Font/face/grid styling stays owned by project_theme() per the figure contract.
    theme(
      strip.text.y    = element_text(angle = 0),
      panel.spacing.y = unit(0.6, "lines"))
}

# ============================================================================
# 4. SOURCE TABLE for the figure (same-stem neighbor; the data behind it).
# ============================================================================

# feature_kind / pvalue are carried through from the master join; default them
# in if an older two_arms_summary.csv lacks them (forward-compatible read).
if (!"feature_kind" %in% names(ev2)) ev2$feature_kind <- NA_character_
if (!"pvalue" %in% names(ev2))       ev2$pvalue       <- NA_real_

panel_table <- ev2 %>%
  dplyr::transmute(
    arm        = as.character(arm),
    arm_track  = unname(TRACK_LABELS[as.character(arm)]),
    method,
    feature,
    feature_kind,
    contrast   = as.character(contrast),
    score, pvalue, padj, direction, significant) %>%
  dplyr::arrange(arm, method, feature, factor(contrast, levels = HEADLINE_CONTRASTS))

# ============================================================================
# 5. SAVE via save_overview (figure + sibling table + README caption, atomic).
# ============================================================================

message("[16_synthesis_viz] Rendering headline two-arms panel ...")

# Own the namespace: drop any stale two_arms_panel.* before the re-render (house
# stale-file guard; save_figure also purges, but call it explicitly to match the
# contract pattern).
purge_figures(STAGE, "two_arms_panel", overview = TRUE, config = FIG_CFG)

fig <- build_two_arms_panel(ev2)

n_ifn_int_sig <- ev2 %>%
  dplyr::filter(arm == ARM_IFN, contrast == "Interaction", significant) %>% nrow()
n_hif_int_sig <- ev2 %>%
  dplyr::filter(arm == ARM_HIF, contrast == "Interaction", significant) %>% nrow()

save_overview(
  fig, STAGE, "two_arms_panel",
  table   = panel_table,
  finding = paste0(
    "Two-arms cGAS-dependence asymmetry, multi-method: the IFN/ISG arm is ",
    "cGAS-dependent (positive, significant Interaction; ", n_ifn_int_sig,
    " significant method-features in the Interaction column), while the ",
    "HIF/glycolysis arm rises in BOTH WT_heat and KO_heat yet is flat in the ",
    "Interaction (", n_hif_int_sig, " significant; no detectable cGAS-dependence ",
    "at n=5). Convergent across GSEA, PROGENy, decoupleR-TF, and per-gene DE ",
    "(plus GATOM/CoReSh where provisioned). PROVISIONAL; n=5/group; NOT proven ",
    "independence."),
  script    = SCRIPT,
  fn        = "build_two_arms_panel",
  config_kv = paste0("thresholds.gsea_fdr=", FDR, "; figures.nes_cap=", CAP,
                     "; colors.diverging; design.contrasts"),
  input     = paste0("03_results/07_synthesis/tables/two_arms_summary.csv; ",
                     "03_results/objects/16_synthesis.rds"),
  how_to_read = paste0(
    "Two stacked tracks: TOP = IFN/ISG arm (cGAS-dependent), BOTTOM = ",
    "HIF/glycolysis arm (no detectable cGAS-dependence at n=5). Rows within a ",
    "track = method-feature glyph rows (GSEA gene sets, PROGENy pathways, ",
    "decoupleR TFs, DE marker genes; GATOM modules where present), ordered by ",
    "Interaction score. Columns = headline contrasts (WT heat | cGAS-KO heat | ",
    "Interaction | Temp main). Tile fill = signed score (orange = up in the ",
    "numerator condition, blue = down), clamped to +/-", CAP, "; the printed ",
    "number is the score. A BLACK RING means padj < ", FDR, " (significant). ",
    "READ THE ASYMMETRY DOWN THE 'Interaction' COLUMN: the IFN/ISG track lights ",
    "up (positive, ringed) = cGAS-dependent; the HIF/glycolysis track goes flat ",
    "/ unringed there while staying lit in BOTH heat columns = no detectable ",
    "cGAS-dependence at n=5. The arm NAMES are labels, not claims that any one ",
    "TF (e.g. HIF1a/HIF2a) is the driver. Claim tier: L3 (DE/enrichment ",
    "statistics; provisional, n=5/group). 'Flat Interaction' is NOT proven ",
    "cGAS-independence - the 1-df interaction is the lowest-powered comparison."),
  # GEOMETRY OVERRIDE (contract-sanctioned width/height passthrough; NOT raw ggsave):
  # a 40-glyph-row x 4-column faceted heatmap with long method+feature y-labels
  # cannot honour the fixed canvas (10x8 screen / 3.5x3.0 print) - the columns
  # collapse to a right-edge sliver and the multi-line x labels overprint. A tall
  # canvas sized to the row count gives every column real width (resolving the
  # x-label overprint too); the per-variant font FLOOR is still enforced by
  # project_theme regardless of inches, so a larger canvas only ever helps. Vector
  # PDF stays editable at this size. The y-axis carries long method+feature labels
  # (a fixed-width text block), so canvas width is spent mostly on the panel: 18in
  # gives each of the 4 columns enough absolute width that the two-line
  # contrast_labels_short stop overprinting (the §B label-collision fix).
  width  = 18,
  height = 16,
  config = FIG_CFG)

message("[16_synthesis_viz] two_arms_panel saved (print + screen + table + caption).")

# ============================================================================
# 6. STAGE README — ensure it exists and captions the panel.
#    save_overview() already wrote the per-artifact caption into the README via
#    write_caption(); here we make sure the stage header / overview prose exists
#    above it (write_caption seeds a header only if the file was absent).
# ============================================================================

readme_path <- file.path("03_results", STAGE, "README.md")
if (file.exists(readme_path)) {
  hdr <- readLines(readme_path, warn = FALSE)
  has_overview <- any(grepl("^# 07_synthesis", hdr)) &&
                  any(grepl("two-arms", hdr, ignore.case = TRUE))
  if (!has_overview) {
    # Prepend a stage-overview block above the existing captions (idempotent-ish:
    # only prepend if our overview marker is absent).
    overview_block <- c(
      "# 07_synthesis - capstone synthesis (two-arms panel + reply memo)",
      "",
      "Stage 07 assembles the cross-arm evidence for the STING-cGAS standard sweep",
      "into one tidy two-arms table and renders the headline figure. NO new",
      "statistics are computed here - it joins/summarises the masters the earlier",
      "arms produced (GSEA, PROGENy, decoupleR-TF, DE; GATOM/CoReSh where present).",
      "",
      "Headline result (the publication-relevant payoff): a cGAS-dependence",
      "ASYMMETRY. The IFN/ISG arm is cGAS-dependent (positive, significant",
      "Interaction); the HIF/glycolysis arm rises in BOTH heat arms but is flat in",
      "the Interaction - no detectable cGAS-dependence at n=5. This is an",
      "asymmetry, NOT proven cGAS-independence (the 1-df interaction is the",
      "lowest-powered comparison). HIF1a/HIF2a are NOT crowned as drivers; the arm",
      "names are labels. Figure claims floor at L3 (DE/enrichment statistics);",
      "mechanism (L7: pseudohypoxia / Complex-I) lives only in the reply memo.",
      "",
      "PROVISIONAL - inferred sample mapping (Hspa1b/Hsph1 thermometer + Cgas),",
      "pending the collaborator sample sheet; n=5/group.",
      "",
      "Artifacts:",
      "- figures/_overview/two_arms_panel.png + two_arms_panel.pdf - the headline panel",
      "- tables/_overview/two_arms_panel.csv - the source table behind the panel",
      "- tables/two_arms_summary.csv - the full cross-arm evidence table (16_synthesis.R)",
      "",
      "Reply memo: docs/_internal/reports/2026-06-24_claim-evidence-memo.md",
      "")
    writeLines(c(overview_block, "", hdr), readme_path)
    message("[16_synthesis_viz] prepended 07_synthesis stage overview to README.")
  }
} else {
  warning("[16_synthesis_viz] README not found at ", readme_path,
          " - save_overview should have created it; check the figure-style shim.")
}

# ============================================================================
# 7. FINAL SUMMARY
# ============================================================================

n_fig <- length(list.files(file.path("03_results", STAGE, "figures"),
                           pattern = "\\.(pdf|png)$", recursive = TRUE))
message(sprintf("[16_synthesis_viz] COMPLETE: %d figure file(s) under %s/figures/.",
                n_fig, STAGE))
if (n_fig == 0)
  warning("[16_synthesis_viz] No figures produced - check errors above.")
message("[16_synthesis_viz] Run from project root: Rscript 02_analysis/scripts/16_synthesis_viz.R")
