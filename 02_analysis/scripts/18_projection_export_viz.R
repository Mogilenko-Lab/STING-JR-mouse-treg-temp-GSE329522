# 18_projection_export_viz.R — VIZ
# =============================================================================
# Human projection export figures (stage 11_projection).
#
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   11_projection
#
# ROLE: VIZ ONLY. Reads the stage tables + manifest that 18_projection_export.R
#   wrote; recomputes NOTHING. Figures ONLY via the figure-style contract
#   (project_theme(config=FIG_CFG) + save_overview()). No inline ggsave()/theme()/
#   raw hex; colors from FIG_CFG$colors.
#
# Figures (each _overview, each carrying its same-stem source table):
#   _overview/human_signature_sizes  up/down set size per exported contrast AFTER mapping
#                                     (human-space counterpart of the stage-17 size bars)
#   _overview/mapping_loss           per contrast, mouse-set size -> human-set size waterfall
#                                     (1:1 kept / many-mapped / dropped-no-ortholog): the
#                                     sanity check that mapping did not decimate the signature
#
# Inputs (read-only; produced by 18_projection_export.R):
#   03_results/11_projection/tables/_overview/{human_signature_sizes,mapping_loss}.csv
#   03_results/human_projection/manifest.csv
#
# Run from project root (after 18_projection_export.R):
#   Rscript 02_analysis/scripts/18_projection_export_viz.R
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

STAGE  <- "11_projection"
SCRIPT <- "02_analysis/scripts/18_projection_export_viz.R"
OV_DIR <- file.path("03_results", STAGE, "tables", FIG_CFG$figures$overview_dir %||% "_overview")

POS <- FIG_CFG$colors$diverging$up   %||% "#B35806"   # up
NEG <- FIG_CFG$colors$diverging$down %||% "#2166AC"   # down
LOSS_COLORS <- c(
  n_mapped_1to1 = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  n_many_mapped = FIG_CFG$colors$okabe_ito$orange       %||% "#E69F00",
  n_unmapped    = "grey65")
# Alarm line: flag a contrast/direction that loses > this fraction to no-ortholog drops.
ALARM_FRAC <- 0.5

# ============================================================================
# 1. GUARD + READ.
# ============================================================================

f_sizes <- file.path(OV_DIR, "human_signature_sizes.csv")
f_loss  <- file.path(OV_DIR, "mapping_loss.csv")
for (f in c(f_sizes, f_loss))
  if (!file.exists(f))
    stop("[18_viz] missing table: ", f,
         " — run 18_projection_export.R first (it only runs after BREAKPOINT-10 sign-off).")

sizes_df <- readr::read_csv(f_sizes, show_col_types = FALSE, progress = FALSE)
loss_df  <- readr::read_csv(f_loss,  show_col_types = FALSE, progress = FALSE)

# Exported-contrast order = manifest order (falls back to appearance order).
man_path <- file.path("03_results", "human_projection", "manifest.csv")
CONTRASTS <- if (file.exists(man_path)) {
  unique(readr::read_csv(man_path, show_col_types = FALSE, progress = FALSE)$contrast)
} else {
  unique(sizes_df$contrast)
}
ord <- function(x) factor(x, levels = CONTRASTS)

message("[18_viz] loaded human sizes (", nrow(sizes_df), ") + mapping loss (", nrow(loss_df), ") rows.")

# ============================================================================
# 2. FIGURE (a): human_signature_sizes — up/down human-set size per exported contrast.
# ============================================================================

sizes_p <- sizes_df %>%
  dplyr::mutate(contrast = ord(contrast),
                direction = factor(direction, levels = c("up", "down")))

fig_sizes <- ggplot(sizes_p, aes(x = contrast, y = n_human, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(aes(label = n_human),
            position = position_dodge(width = 0.8), vjust = -0.35,
            size = (FIG_CFG$figures$label_size %||% 4) * 0.85, colour = "grey20") +
  scale_fill_manual(values = c(up = POS, down = NEG),
                    labels = c(up = "up (39 > 37 / numerator)", down = "down"),
                    name = "DE direction") +
  scale_x_discrete(labels = function(x) contrast_label(x, short = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Human projection signature sizes (after ortholog mapping)",
       subtitle = paste0("Frozen human-symbol up/down set sizes per exported contrast, gate = ",
                        unique(sizes_df$gate)[1], "."),
       x = NULL, y = "human genes in set",
       caption = "Human HGNC symbols after mouse->human mapping (babelgene). Claim tier: L3 (DE statistics).") +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 0))

# ============================================================================
# 3. FIGURE (b): mapping_loss — mouse-set -> human-set waterfall per contrast/direction.
# ============================================================================

loss_long <- loss_df %>%
  tidyr::pivot_longer(cols = c(n_mapped_1to1, n_many_mapped, n_unmapped),
                      names_to = "category", values_to = "n") %>%
  dplyr::mutate(contrast = ord(contrast),
                direction = factor(direction, levels = c("up", "down")),
                category = factor(category,
                                  levels = c("n_mapped_1to1", "n_many_mapped", "n_unmapped")))

# human-out annotation + alarm flag per contrast/direction
loss_ann <- loss_df %>%
  dplyr::mutate(contrast = ord(contrast),
                direction = factor(direction, levels = c("up", "down")),
                frac_unmapped = ifelse(n_mouse > 0, n_unmapped / n_mouse, 0),
                alarm = frac_unmapped > ALARM_FRAC)

fig_loss <- ggplot(loss_long, aes(x = contrast, y = n, fill = category)) +
  geom_col(width = 0.72) +
  geom_text(data = loss_ann,
            aes(x = contrast, y = n_mouse, fill = NULL,
                label = sprintf("->%d", n_human)),
            vjust = -0.35, size = (FIG_CFG$figures$label_size %||% 4) * 0.8, colour = "grey15") +
  geom_point(data = dplyr::filter(loss_ann, alarm),
             aes(x = contrast, y = n_mouse, fill = NULL),
             shape = 8, size = 2.4, colour = "red3", show.legend = FALSE) +
  facet_wrap(~ direction, nrow = 1,
             labeller = labeller(direction = c(up = "UP set", down = "DOWN set"))) +
  scale_fill_manual(values = LOSS_COLORS,
                    labels = c(n_mapped_1to1 = "mapped 1:1", n_many_mapped = "one mouse -> many human",
                               n_unmapped = "dropped (no human ortholog)"),
                    name = "mouse gene fate") +
  scale_x_discrete(labels = function(x) contrast_label(x, short = TRUE)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Ortholog mapping loss: mouse set -> human set",
       subtitle = paste0("Each bar = the mouse up/down set decomposed by fate; '->N' = resulting human ",
                        "set size. Red * = > ", round(100 * ALARM_FRAC), "% dropped (alarm)."),
       x = NULL, y = "mouse genes in set",
       caption = paste0("Bar height = mouse-set size; green kept 1:1, orange one->many (each human gets ",
                        "the mouse t), grey dropped. Claim tier: L3. Correlative input, not a causal claim.")) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 0))

# ============================================================================
# 4. SAVE via save_overview (figure + same-stem source table + README caption).
# ============================================================================

purge_figures(STAGE, "human_signature_sizes", overview = TRUE, config = FIG_CFG)
purge_figures(STAGE, "mapping_loss",          overview = TRUE, config = FIG_CFG)

save_overview(
  fig_sizes, STAGE, "human_signature_sizes",
  table   = sizes_df,
  finding = "Frozen human up/down set sizes per exported contrast after mouse->human ortholog mapping.",
  script  = SCRIPT, fn = "ggplot(geom_col)",
  config_kv = "decisions.projection.gate; colors.diverging",
  input   = f_sizes,
  how_to_read = paste0("Grouped bars per exported contrast; orange = up (higher in numerator / 39 C), blue ",
                       "= down; numbers = human genes in the frozen set. This is the human-space counterpart ",
                       "of the stage-17 (mouse) size bars. Claim tier: L3 (DE statistics); provisional ",
                       "sample mapping, n=5/group."),
  width = 11, height = 6.5,
  config = FIG_CFG)

save_overview(
  fig_loss, STAGE, "mapping_loss",
  table   = loss_df,
  finding = "Mapping did not decimate the primary signature: most mouse genes map 1:1; dropped/many-mapped fractions are small (alarm-flagged if > half lost).",
  script  = SCRIPT, fn = "ggplot(geom_col)",
  config_kv = paste0("alarm_frac=", ALARM_FRAC, "; babelgene(offline); decisions.projection.ortholog_ambiguity"),
  input   = f_loss,
  how_to_read = paste0("Per contrast (x) and direction (facet), the mouse up/down set is stacked by fate: ",
                       "green = mapped 1:1, orange = one-mouse->many-human (each human ortholog inherits the ",
                       "mouse t in the ranked list; unioned in the binary set), grey = dropped (no human ",
                       "ortholog). Bar height = mouse-set size; '->N' above = resulting human-set size. A red ",
                       "* marks > ", round(100 * ALARM_FRAC), "% dropped — a decimated signature to demote. ",
                       "Read it as the sanity check that the frozen human sets are not hollow. Claim tier: L3."),
  width = 11, height = 6.5,
  config = FIG_CFG)

# ============================================================================
# 5. FINAL SUMMARY
# ============================================================================

n_fig <- length(list.files(file.path("03_results", STAGE, "figures"),
                           pattern = "\\.(pdf|png)$", recursive = TRUE))
message(sprintf("[18_viz] COMPLETE: %d figure file(s) under %s/figures/.", n_fig, STAGE))
if (n_fig == 0) warning("[18_viz] No figures produced — check errors above.")
