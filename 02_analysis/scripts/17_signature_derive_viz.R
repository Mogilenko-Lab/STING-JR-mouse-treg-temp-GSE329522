# 17_signature_derive_viz.R — VIZ
# =============================================================================
# Projection-signature derivation figures (stage 10_signature).
#
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   10_signature
#
# ROLE: VIZ ONLY. Reads the three _overview tables, plus the checkpoint for labels, that
#   17_signature_derive.R wrote. Figures go through the figure-style contract alone:
#   project_theme(config=FIG_CFG) + save_overview() (dual pdf+png + sibling source table +
#   README caption, atomic), with colors from FIG_CFG$colors.
#
# Figures (each _overview, each carrying its same-stem source table):
#   _overview/signature_sizes           grouped up/down bars per contrast, faceted by gate
#   _overview/updown_overlap            Jaccard heatmap of up- and down-sets across contrasts
#   _overview/ortholog_coverage_preview stacked mapped_1to1 / one_to_many / unmapped bars per contrast
#
# These three figures are the lens for BREAKPOINT 10. They let the human judge which
# contrasts are distinct, which gate to commit to, and how much ortholog mapping costs a
# signature, before stage 18 freezes the contract.
#
# Inputs (read-only):
#   03_results/10_signature/tables/_overview/{signature_sizes,updown_overlap,ortholog_coverage_preview}.csv
#   03_results/objects/17_signature_sets.rds
#
# Run from project root:
#   Rscript 02_analysis/scripts/17_signature_derive_viz.R
# =============================================================================

source("02_analysis/helpers/figure_style.R")   # project_theme, save_overview, purge_figures, FIG_CFG, contrast_label

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

STAGE  <- "10_signature"
SCRIPT <- "02_analysis/scripts/17_signature_derive_viz.R"
OV_DIR <- file.path("03_results", STAGE, "tables", FIG_CFG$figures$overview_dir %||% "_overview")

DE_FDR   <- as.numeric(FIG_CFG$thresholds$de_fdr   %||% 0.05)
DE_LOGFC <- as.numeric(FIG_CFG$thresholds$de_logfc %||% 1.0)
# Heatmap uses the stringent config gate for the plotted panel (the source table keeps both).
HEATMAP_GATE <- "fdr_logfc"

# Diverging triplet from config (down = blue, up = orange).
POS <- FIG_CFG$colors$diverging$up   %||% "#B35806"   # up
NEG <- FIG_CFG$colors$diverging$down %||% "#2166AC"   # down
# Coverage categories — colorblind-safe okabe-ito picks.
COV_COLORS <- c(
  mapped_1to1 = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  one_to_many = FIG_CFG$colors$okabe_ito$orange       %||% "#E69F00",
  unmapped    = "grey65")

# ============================================================================
# 1. GUARD + READ the overview tables (and checkpoint for contrast order).
# ============================================================================

f_sizes <- file.path(OV_DIR, "signature_sizes.csv")
f_ov    <- file.path(OV_DIR, "updown_overlap.csv")
f_cov   <- file.path(OV_DIR, "ortholog_coverage_preview.csv")
rds     <- file.path("03_results", "objects", "17_signature_sets.rds")
for (f in c(f_sizes, f_ov, f_cov))
  if (!file.exists(f)) stop("[17_viz] missing overview table: ", f, " — run 17_signature_derive.R first.")

sizes_df <- readr::read_csv(f_sizes, show_col_types = FALSE, progress = FALSE)
ov_df    <- readr::read_csv(f_ov,    show_col_types = FALSE, progress = FALSE)
cov_df   <- readr::read_csv(f_cov,   show_col_types = FALSE, progress = FALSE)

CONTRASTS <- if (file.exists(rds)) readRDS(rds)$contrasts else unique(sizes_df$contrast)
# Contrast factor order = config/checkpoint order (WT_heat first).
ord <- function(x) factor(x, levels = CONTRASTS)

message("[17_viz] loaded overview tables (sizes=", nrow(sizes_df),
        " overlap=", nrow(ov_df), " coverage=", nrow(cov_df), " rows).")

# ============================================================================
# 2. FIGURE (a): signature_sizes — grouped up/down bars per contrast, faceted by gate.
# ============================================================================

sizes_p <- sizes_df %>%
  dplyr::mutate(contrast = ord(contrast),
                direction = factor(direction, levels = c("up", "down")),
                gate = factor(gate, levels = c("fdr_only", "fdr_logfc")))

fig_sizes <- ggplot(sizes_p, aes(x = contrast, y = n_genes, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(aes(label = n_genes),
            position = position_dodge(width = 0.8), vjust = -0.35,
            size = (FIG_CFG$figures$label_size %||% 4) * 0.8, colour = "grey20") +
  facet_wrap(~ gate, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(up = POS, down = NEG),
                    labels = c(up = "up (39 > 37 / numerator)", down = "down"),
                    name = "DE direction") +
  scale_x_discrete(labels = function(x) contrast_label(x, short = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Projection signature sizes per contrast",
       subtitle = paste0("Binary up/down set sizes (mouse symbols) at both gates. ",
                         "FDR-only (top) vs FDR+|logFC| (bottom, stringent)."),
       x = NULL, y = "genes in set",
       caption = paste0("Gate: FDR = adj.P.Val < ", DE_FDR,
                        "; +|logFC| adds |log2FC| >= ", DE_LOGFC,
                        ". Sets are mouse symbols (pre-ortholog). Claim tier: L3 (DE statistics).")) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 0))

# ============================================================================
# 3. FIGURE (b): updown_overlap — Jaccard heatmap of up- and down-sets (one gate).
# ============================================================================

ov_plot_df <- ov_df %>%
  dplyr::filter(gate == HEATMAP_GATE) %>%
  dplyr::mutate(contrast_a = ord(contrast_a), contrast_b = ord(contrast_b),
                direction = factor(direction, levels = c("up", "down")))

fig_overlap <- ggplot(ov_plot_df, aes(x = contrast_a, y = contrast_b, fill = jaccard)) +
  geom_tile(colour = "grey92", linewidth = 0.3) +
  geom_text(aes(label = ifelse(is.na(jaccard), "", sprintf("%.2f", jaccard))),
            size = (FIG_CFG$figures$label_size %||% 4) * 0.7, colour = "grey15") +
  facet_wrap(~ direction, nrow = 1,
             labeller = labeller(direction = c(up = "UP sets", down = "DOWN sets"))) +
  scale_fill_gradient(low = "grey97", high = NEG, limits = c(0, 1), na.value = "grey90",
                      name = "Jaccard") +
  scale_x_discrete(labels = function(x) contrast_label(x, short = TRUE)) +
  scale_y_discrete(labels = function(x) contrast_label(x, short = TRUE)) +
  coord_equal() +
  labs(title = "Contrast up/down-set overlap (Jaccard)",
       subtitle = paste0("How distinct the candidate signatures are, gate = ", HEATMAP_GATE,
                        ". Off-diagonal near 0 = distinct sets."),
       x = NULL, y = NULL,
       caption = "Jaccard = |A n B| / |A u B| of mouse-symbol sets. Diagonal = 1 (self). Claim tier: L3.") +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank())

# ============================================================================
# 4. FIGURE (c): ortholog_coverage_preview — stacked mapped/many/unmapped per contrast.
# ============================================================================

cov_long <- cov_df %>%
  tidyr::pivot_longer(cols = c(mapped_1to1, one_to_many, unmapped),
                      names_to = "category", values_to = "n") %>%
  dplyr::mutate(contrast = ord(contrast),
                gate = factor(gate, levels = c("fdr_only", "fdr_logfc")),
                category = factor(category, levels = c("mapped_1to1", "one_to_many", "unmapped")))

fig_cov <- ggplot(cov_long, aes(x = contrast, y = n, fill = category)) +
  geom_col(position = "fill", width = 0.72) +
  facet_wrap(~ gate, ncol = 1) +
  scale_fill_manual(values = COV_COLORS,
                    labels = c(mapped_1to1 = "mapped 1:1", one_to_many = "one mouse -> many human",
                               unmapped = "no human ortholog"),
                    name = "ortholog mapping") +
  scale_x_discrete(labels = function(x) contrast_label(x, short = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(title = "Dry-run ortholog coverage of the up+down sets",
       subtitle = paste0("Fraction of each contrast's significant mouse genes that map to human ",
                        "(babelgene, preview only). High 'unmapped' = mapping loss risk."),
       x = NULL, y = "fraction of up+down genes",
       caption = paste0("Preview via ortholog_utils (min_support from decisions.projection). ",
                        "Does NOT feed stage 18 (which re-maps the frozen sets). Claim tier: L3.")) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 0))

# ============================================================================
# 5. SAVE via save_overview (figure + same-stem source table + README caption).
#    Each figure's neighbor table = the overview table it was built from (passed
#    through unchanged; save_overview re-writes it byte-stably).
# ============================================================================

purge_figures(STAGE, "signature_sizes",           overview = TRUE, config = FIG_CFG)
purge_figures(STAGE, "updown_overlap",            overview = TRUE, config = FIG_CFG)
purge_figures(STAGE, "ortholog_coverage_preview", overview = TRUE, config = FIG_CFG)

save_overview(
  fig_sizes, STAGE, "signature_sizes",
  table   = sizes_df,
  finding = paste0("Per-contrast up/down projection-set sizes in mouse symbols, at both gates; ",
                   "the FDR+|logFC| gate is markedly more stringent than FDR-only."),
  script  = SCRIPT, fn = "ggplot(geom_col)",
  config_kv = paste0("thresholds.de_fdr=", DE_FDR, "; thresholds.de_logfc=", DE_LOGFC,
                     "; gsea.rank_metric=t; colors.diverging"),
  input   = file.path(OV_DIR, "signature_sizes.csv"),
  how_to_read = paste0("Grouped bars per contrast; orange = up (higher in the numerator / 39 C), ",
                       "blue = down. TOP facet = FDR-only gate (adj.P.Val < ", DE_FDR,
                       "), BOTTOM = FDR + |log2FC| >= ", DE_LOGFC, " (stringent). Numbers above bars ",
                       "= genes in that set (mouse symbols, pre-ortholog). Use this to pick the gate ",
                       "before the signatures are frozen: FDR+logFC can decimate small contrasts. ",
                       "Claim tier: L3 (DE statistics), n=5/group. ", sample_mapping_caption()),
  # GEOMETRY OVERRIDE (contract-sanctioned passthrough; NOT raw ggsave): 7 contrasts
  # with two-line labels collide on the fixed 8.5in canvas and the caption clips. A
  # wider canvas gives each column real width; project_theme still enforces the font
  # FLOOR, so a larger canvas only helps.
  width = 12, height = 7.5,
  config = FIG_CFG)

save_overview(
  fig_overlap, STAGE, "updown_overlap",
  table   = ov_df,
  finding = "Off-diagonal Jaccard shows how distinct the candidate signatures (WT_heat / Temp_main / Geno_at_39 / Interaction) are.",
  script  = SCRIPT, fn = "ggplot(geom_tile)",
  config_kv = paste0("gate=", HEATMAP_GATE, "; colors.diverging.down"),
  input   = file.path(OV_DIR, "updown_overlap.csv"),
  how_to_read = paste0("Symmetric heatmap; each tile = Jaccard overlap of two contrasts' UP sets (left ",
                       "facet) or DOWN sets (right facet), gate = ", HEATMAP_GATE, ". Darker = more shared ",
                       "genes; diagonal = 1 (self). Off-diagonal near 0 means the signatures are distinct ",
                       "programs (good — comparators add information). The source table carries BOTH gates. ",
                       "Claim tier: L3."),
  width = 12, height = 6.5,
  config = FIG_CFG)

save_overview(
  fig_cov, STAGE, "ortholog_coverage_preview",
  table   = cov_df,
  finding = "Dry-run mouse->human coverage of each contrast's up+down set; a large 'unmapped' fraction flags mapping-loss risk before the freeze.",
  script  = SCRIPT, fn = "ggplot(geom_col, position=fill)",
  config_kv = "decisions.projection.ortholog_ambiguity.min_support; babelgene(offline)",
  input   = file.path(OV_DIR, "ortholog_coverage_preview.csv"),
  how_to_read = paste0("Stacked bars = fraction of a contrast's significant mouse genes (up+down) that map ",
                       "1:1 to a human ortholog (green), one-mouse->many-human (orange), or have no human ",
                       "ortholog (grey, dropped downstream). TOP facet = FDR-only, BOTTOM = FDR+logFC. This ",
                       "is an offline babelgene PREVIEW for judging mapping loss. The frozen human sets ",
                       "are re-mapped from the approved gene lists by 18_projection_export.R, and that ",
                       "run is what 03_results/11_projection/ reports. Claim tier: L3."),
  width = 12, height = 7.5,
  config = FIG_CFG)

# ============================================================================
# 6. FINAL SUMMARY
# ============================================================================

n_fig <- length(list.files(file.path("03_results", STAGE, "figures"),
                           pattern = "\\.(pdf|png)$", recursive = TRUE))
message(sprintf("[17_viz] COMPLETE: %d figure file(s) under %s/figures/.", n_fig, STAGE))
if (n_fig == 0) warning("[17_viz] No figures produced — check errors above.")
message("[17_viz] Next: render 02_analysis/notebooks/17_signature_review/ (Rscript 02_analysis/notebooks/render.R) — BREAKPOINT 10.")
