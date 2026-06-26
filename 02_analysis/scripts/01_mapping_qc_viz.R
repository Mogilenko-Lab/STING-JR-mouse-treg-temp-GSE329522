#!/usr/bin/env Rscript
# =============================================================================
# 01_mapping_qc_viz.R  --  PHASE 1: Sample-mapping QC figures (VISUALIZE)
# =============================================================================
# Phase:        1 (SCIENCE GATE -- visualization of the inferred sample mapping)
# Role:         VISUALIZE half of the "normalize-then-visualize" split. Reads the
#               plot-ready tidy tables emitted by 01_mapping_qc.R and renders the
#               figures. Performs NO statistics (no prcomp, no collapse, no
#               thermometer/Cgas derivation) -- only cosmetic reshaping (factor
#               ordering, faceting). Every save routes through the figure-style
#               contract (project_theme / save_figure / write_caption): dual
#               .print.pdf + .screen.png variants, font floors enforced.
#
# Inputs:       03_results/01_qc/tables/fig1a_thermometer_data.csv
#               03_results/01_qc/tables/fig1b_cgas_data.csv
#               03_results/02_eda/tables/fig1c_pca_data.csv
#               03_results/02_eda/tables/fig1c_pca_varexp.csv
#               03_results/01_qc/tables/fig1d_scramble_data.csv
#               03_results/objects/01_eda.rds  (labels only: cgas_symbol)
# Outputs:      03_results/01_qc/figures/fig1a_thermometer.{print.pdf,screen.png}
#               03_results/01_qc/figures/fig1b_cgas.{print.pdf,screen.png}
#               03_results/02_eda/figures/fig1c_pca_2x2.{print.pdf,screen.png}
#               03_results/02_eda/figures/fig1c_pca_scree.{print.pdf,screen.png}
#               03_results/01_qc/figures/fig1d_scramble.{print.pdf,screen.png}
#
# NOTE: these 4 (now 5) stems are FLAT stage-level figures (label-blind mapping
#       QC) -- there are no DE contrasts at Phase 1, so they route through
#       save_figure(..., variant = "both") + explicit write_caption(), NOT
#       save_overview() (which requires contrast= or overview=). The 20-category
#       panels (1a/1b/1d/scree) get explicit width=/height= so the print column
#       does not crush the axis labels; geometry (width/height) is the only
#       per-call knob -- NO raw ggsave(width=), NO inline theme() styling. Short
#       x-tick labels (library number / PC index) keep the 20-category axes
#       legible without an axis rotation (which save_figure would clobber).
#
# Dependencies: figure_style.R (contract); config.R (stage_dir, DIVERGING_COLORS,
#               provisional_caption, DIR_OBJECTS); ggplot2, dplyr, ggrepel, tidyr
# =============================================================================

source("02_analysis/helpers/figure_style.R")   # FIG_CFG, project_theme(),
                                                # save_figure(), write_caption(),
                                                # purge_figures()
source("02_analysis/config/config.R")           # stage_dir(), DIVERGING_COLORS,
                                                # provisional_caption(), DIR_OBJECTS

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(tidyr)
})

# -----------------------------------------------------------------------------
# 0. Style-contract constants + provenance strings (caption metadata)
# -----------------------------------------------------------------------------
SCRIPT <- "02_analysis/scripts/01_mapping_qc_viz.R"
CFG_KV <- "figures.width_column=3.5; figures.base_size_column=9; figures.base_size=16"

# -----------------------------------------------------------------------------
# 1. Read the plot-ready tidy tables (no recomputation)
# -----------------------------------------------------------------------------
qc_tables  <- stage_dir("01_qc", "tables")
eda_tables <- stage_dir("02_eda", "tables")

fig1a_data   <- read.csv(file.path(qc_tables,  "fig1a_thermometer_data.csv"), check.names = FALSE, stringsAsFactors = FALSE)
fig1b_data   <- read.csv(file.path(qc_tables,  "fig1b_cgas_data.csv"),        check.names = FALSE, stringsAsFactors = FALSE)
fig1c_data   <- read.csv(file.path(eda_tables, "fig1c_pca_data.csv"),         check.names = FALSE, stringsAsFactors = FALSE)
fig1c_varexp <- read.csv(file.path(eda_tables, "fig1c_pca_varexp.csv"),       check.names = FALSE, stringsAsFactors = FALSE)
fig1d_data   <- read.csv(file.path(qc_tables,  "fig1d_scramble_data.csv"),    check.names = FALSE, stringsAsFactors = FALSE)

# Labels only (no stats): Cgas symbol used in fig1b title.
eda <- readRDS(file.path(DIR_OBJECTS, "01_eda.rds"))
cgas_symbol <- eda$cgas_symbol

# Canonical library order (temperature-major) = order rows appear in fig1b_data.
lib_levels <- fig1b_data$lib

# Shared plotting scaffolding (categorical data-as-scale palettes -- kept) -----
temp_cols <- c("37" = DIVERGING_COLORS$negative, "39" = DIVERGING_COLORS$positive)
geno_cols <- c("WT" = "#1B7837", "cGASKO" = "#762A83")

# Read variance-explained + n_top (labels only) from the varexp table.
percentVar <- setNames(fig1c_varexp$pct_var, fig1c_varexp$PC)
n_top <- fig1c_varexp$n_top[1]

# NOTE on axis-label legibility: save_figure() re-applies project_theme() to the
# plot object per variant (figure_helpers.R), which CLOBBERS any post-theme
# axis.text override (rotation/size). So the 20-category panels are kept legible
# by (a) a roomy explicit width= per call and (b) SHORT x-tick labels (the bare
# library number 021..040 / the PC index 1..20), not by an axis rotation that the
# contract would silently drop. The deposited library ids strip to 2-3 chars and
# sit horizontally without collision at the chosen widths.

# Purge any stale single-.pdf stems from the old off-contract run so this run owns
# its figure namespace (save_figure also purges per-stem; belt-and-suspenders).
purge_figures("01_qc", "fig1", config = FIG_CFG)
purge_figures("02_eda", "fig1c", config = FIG_CFG)

# =============================================================================
# A) FIG 1a -- THERMOMETER (label-blind temperature validation)
# =============================================================================
message("\n==================  A) FIG 1a THERMOMETER  ==================")

thermo_plot <- fig1a_data
thermo_plot$temp <- factor(as.character(thermo_plot$inferred_temp), levels = c("37", "39"))
thermo_plot$lib  <- factor(thermo_plot$lib, levels = lib_levels)
thermo_plot$gene <- factor(thermo_plot$marker, levels = unique(fig1a_data$marker))

fig1a <- ggplot(thermo_plot, aes(x = lib, y = value, fill = temp)) +
  geom_col() +
  facet_wrap(~ gene, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = temp_cols, name = "Inferred temp (C)") +
  labs(
    title = "Fig 1a -- Heat-shock thermometer (label-blind)",
    subtitle = "Per-library log2(CPM+0.5); libraries 031-040 (inferred 39C) should be systematically higher",
    x = "Library (12630-RS-)", y = "log2(CPM + 0.5)"
  ) +
  project_theme(config = FIG_CFG)

# 20-library x-axis x 4 stacked facets: give the print variant a tall, roomy
# canvas so the short library labels stay legible at the 9pt print floor.
save_figure(fig1a, "01_qc", "fig1a_thermometer",
            variant = "both", config = FIG_CFG, width = 11, height = 10)
message("[FIG] saved fig1a_thermometer (print+screen)")

write_caption(
  stage      = "01_qc",
  filename   = "figures/fig1a_thermometer.png",
  finding    = paste0(
    "Heat-shock thermometer (label-blind temperature validation): per-library log2(CPM+0.5) ",
    "for the 4 heat-shock marker genes (Hspa1b, Hsph1, Hspa1a, Dnajb1), one facet per gene, ",
    "all 20 libraries (12630-RS-021..040) on the x-axis, colored by INFERRED temperature. ",
    "Libraries 031-040 (inferred 39C, red) run systematically higher than 021-030 (inferred 37C, blue) ",
    "in every marker, confirming the deposited CPM column order is temperature-major. ",
    "Label-blind: the temperature is data-derived (marker expression), not read off the deposited labels. ",
    "Claim tier: PROVISIONAL (pending collaborator sample sheet)."),
  script     = SCRIPT,
  fn         = "geom_col + facet_wrap (project_theme + save_figure)",
  config_kv  = CFG_KV,
  input      = "03_results/01_qc/tables/fig1a_thermometer_data.csv",
  how_to_read = paste0(
    "x = library (12630-RS-021..040, temperature-major order); y = log2(CPM+0.5). ",
    "Fill color = INFERRED temperature: blue = 37C (negative pole of the diverging map), ",
    "red = 39C (positive pole). One facet per heat-shock marker (free y). ",
    "Read: the hot half (031-040) bars should tower over the cool half (021-030) in each facet. ",
    "The labels are inferred from marker expression (label-blind), NOT the deposited sample labels. ",
    "Claim tier: PROVISIONAL, pending the collaborator sample sheet."),
  config     = FIG_CFG
)

# =============================================================================
# B) FIG 1b -- CGAS GENOTYPE CHECK
# =============================================================================
message("\n==================  B) FIG 1b CGAS  ==================")

cgas_df <- data.frame(
  library_id = fig1b_data$library_id,
  lib        = factor(fig1b_data$lib, levels = lib_levels),
  log2cpm    = fig1b_data$cgas_value,
  genotype   = fig1b_data$inferred_genotype,
  temp       = as.character(fig1b_data$inferred_temp),
  stringsAsFactors = FALSE
)

fig1b <- ggplot(cgas_df, aes(x = lib, y = log2cpm, fill = genotype)) +
  geom_col() +
  facet_grid(. ~ temp, scales = "free_x", space = "free_x",
             labeller = labeller(temp = function(x) paste0(x, "C (inferred)"))) +
  scale_fill_manual(values = geno_cols, name = "Inferred genotype") +
  labs(
    title = sprintf("Fig 1b -- %s genotype check (label-blind within temp half)", cgas_symbol),
    subtitle = "WT should exceed cGAS-KO within EACH temperature half",
    x = "Library (12630-RS-)", y = "log2(CPM + 0.5)"
  ) +
  project_theme(config = FIG_CFG)

# 20 libraries split across 2 temperature facets: widen the print canvas so the
# short library labels do not collide at the 9pt floor.
save_figure(fig1b, "01_qc", "fig1b_cgas",
            variant = "both", config = FIG_CFG, width = 10, height = 5.5)
message("[FIG] saved fig1b_cgas (print+screen)")

write_caption(
  stage      = "01_qc",
  filename   = "figures/fig1b_cgas.png",
  finding    = sprintf(paste0(
    "%s genotype check (label-blind within each temperature half): per-library log2(CPM+0.5) of %s, ",
    "all 20 libraries split into the inferred-37C and inferred-39C facets, colored by INFERRED genotype. ",
    "Within EACH temperature half, the WT libraries (green) exceed the cGAS-KO libraries (purple), ",
    "the expected signature of a Cgas knockout -- confirming the genotype axis is data-recoverable. ",
    "Label-blind: genotype is inferred from %s expression, not the deposited labels. ",
    "Claim tier: PROVISIONAL (pending collaborator sample sheet)."),
    cgas_symbol, cgas_symbol, cgas_symbol),
  script     = SCRIPT,
  fn         = "geom_col + facet_grid (project_theme + save_figure)",
  config_kv  = CFG_KV,
  input      = "03_results/01_qc/tables/fig1b_cgas_data.csv",
  how_to_read = sprintf(paste0(
    "x = library (12630-RS-, faceted into inferred-37C | inferred-39C halves); y = %s log2(CPM+0.5). ",
    "Fill color = INFERRED genotype: green = WT, purple = cGASKO. ",
    "Read WITHIN each facet (each temperature half): WT bars should sit above cGAS-KO bars. ",
    "The comparison is intentionally within-temperature so the heat effect does not confound the genotype call. ",
    "Genotype is inferred from %s expression (label-blind). Claim tier: PROVISIONAL."),
    cgas_symbol, cgas_symbol),
  config     = FIG_CFG
)

# =============================================================================
# C) FIG 1c -- PCA (label-blind, overlay inferred labels)
# =============================================================================
message("\n==================  C) FIG 1c PCA 2x2  ==================")

pca_df <- fig1c_data
pca_df$temp     <- as.character(pca_df$temp)
pca_df$genotype <- as.character(pca_df$genotype)

fig1c <- ggplot(pca_df, aes(x = PC1, y = PC2, color = temp, shape = genotype)) +
  geom_point(size = 4, stroke = 1.1) +
  ggrepel::geom_text_repel(aes(label = lib), color = "grey30",
                           max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = temp_cols, name = "Inferred temp (C)") +
  scale_shape_manual(values = c("WT" = 16, "cGASKO" = 17), name = "Inferred genotype") +
  labs(
    title = sprintf("Fig 1c -- PCA on top %d variable genes (label-blind, labels overlaid)", n_top),
    subtitle = "Expect PC1/PC2 to form the 2x2: temperature (primary) and genotype (secondary)",
    x = sprintf("PC1: %.1f%% var", percentVar[["PC1"]]),
    y = sprintf("PC2: %.1f%% var", percentVar[["PC2"]])
  ) +
  project_theme(config = FIG_CFG)

# Repel labels need a roomy square canvas so they fit without a zero-dimension
# viewport at the print column.
save_figure(fig1c, "02_eda", "fig1c_pca_2x2",
            variant = "both", config = FIG_CFG, width = 9, height = 7.5)
message("[FIG] saved fig1c_pca_2x2 (print+screen)")

write_caption(
  stage      = "02_eda",
  filename   = "figures/fig1c_pca_2x2.png",
  finding    = sprintf(paste0(
    "PCA on the top %d variable genes (label-blind): each library plotted on PC1 (%.1f%% var) ",
    "vs PC2 (%.1f%% var), colored by INFERRED temperature and shaped by INFERRED genotype, with ",
    "library ids overlaid. The libraries split into the expected 2x2: temperature dominates PC1 ",
    "(37C vs 39C separate left/right) and genotype tracks the secondary axis -- confirming both ",
    "experimental factors are recoverable from expression alone. ",
    "Claim tier: PROVISIONAL (pending collaborator sample sheet)."),
    n_top, percentVar[["PC1"]], percentVar[["PC2"]]),
  script     = SCRIPT,
  fn         = "geom_point + geom_text_repel (project_theme + save_figure)",
  config_kv  = CFG_KV,
  input      = "03_results/02_eda/tables/fig1c_pca_data.csv",
  how_to_read = paste0(
    "x = PC1 (% var in axis title); y = PC2. Point COLOR = inferred temperature ",
    "(blue = 37C, red = 39C); point SHAPE = inferred genotype (circle = WT, triangle = cGASKO). ",
    "Grey text = library id. Read: temperature should separate along PC1 (the big axis), genotype ",
    "along PC2 -- the canonical 2x2 of a temperature x genotype design. ",
    "Labels are inferred (label-blind), NOT deposited. Companion scree panel: fig1c_pca_scree. ",
    "Claim tier: PROVISIONAL."),
  config     = FIG_CFG
)

# --- FIG 1c (scree) -- variance-explained per PC (viz-only; data already emitted)
message("\n==================  C') FIG 1c SCREE  ==================")

scree_df <- fig1c_varexp[order(-fig1c_varexp$pct_var), ]
# SHORT x-tick labels (bare PC index 1..20): the contract re-applies project_theme
# inside save_figure and would clobber any axis-text rotation, so a wide "PC1..PC20"
# tick set would collide. Numeric indices stay horizontal + legible (mirrors the
# 2-char library ticks on 1a/1b/1d).
scree_df$pc_idx  <- factor(seq_len(nrow(scree_df)),
                           levels = seq_len(nrow(scree_df)))
scree_df$is_lead <- scree_df$PC %in% c("PC1", "PC2")

fig1c_scree <- ggplot(scree_df, aes(x = pc_idx, y = pct_var, fill = is_lead)) +
  geom_col() +
  scale_fill_manual(values = c(`TRUE` = DIVERGING_COLORS$positive,
                               `FALSE` = "grey70"),
                    name = NULL, labels = c(`TRUE` = "PC1/PC2 (2x2 axes)", `FALSE` = "PC3+")) +
  labs(
    title = sprintf("Fig 1c (scree) -- variance explained per PC (top %d variable genes)", n_top),
    subtitle = sprintf("PC1 dominates (%.1f%%); PC2 (%.1f%%) carries the secondary genotype axis",
                       percentVar[["PC1"]], percentVar[["PC2"]]),
    x = "Principal component (rank order)", y = "Variance explained (%)"
  ) +
  project_theme(config = FIG_CFG)

save_figure(fig1c_scree, "02_eda", "fig1c_pca_scree",
            variant = "both", config = FIG_CFG, width = 8, height = 5)
message("[FIG] saved fig1c_pca_scree (print+screen)")

write_caption(
  stage      = "02_eda",
  filename   = "figures/fig1c_pca_scree.png",
  finding    = sprintf(paste0(
    "PCA scree / variance-explained panel (top %d variable genes): %% variance per principal ",
    "component, PC1/PC2 highlighted. PC1 alone captures %.1f%% of the variance (the temperature ",
    "axis) and PC2 a further %.1f%% (the secondary genotype axis); PC3+ contribute little, ",
    "confirming the 2x2 structure seen in fig1c_pca_2x2 is concentrated in the leading two axes. ",
    "Claim tier: PROVISIONAL (pending collaborator sample sheet)."),
    n_top, percentVar[["PC1"]], percentVar[["PC2"]]),
  script     = SCRIPT,
  fn         = "geom_col (project_theme + save_figure)",
  config_kv  = CFG_KV,
  input      = "03_results/02_eda/tables/fig1c_pca_varexp.csv",
  how_to_read = paste0(
    "x = principal component (descending variance); y = % variance explained. ",
    "Orange bars = PC1/PC2 (the two axes that carry the temperature x genotype 2x2); grey = PC3+. ",
    "Read: a steep drop after PC2 means the experimental structure lives in the leading axes ",
    "(the fig1c_pca_2x2 scatter is therefore a faithful 2D summary). ",
    "Companion scatter: fig1c_pca_2x2. Claim tier: PROVISIONAL."),
  config     = FIG_CFG
)

# =============================================================================
# D) FIG 1d -- SCRAMBLE EXHIBIT (competence exhibit)
# =============================================================================
message("\n==================  D) FIG 1d SCRAMBLE  ==================")

scram <- data.frame(
  lib              = factor(fig1d_data$lib, levels = lib_levels),
  inferred_group   = fig1d_data$inferred_condition,
  naive_join_group = fig1d_data$naive_positional_condition,
  discordant       = as.logical(fig1d_data$discordant),
  stringsAsFactors = FALSE
)
n_disc <- sum(scram$discordant)

# Tile heatmap: rows = the two label sources, cols = libraries; fill = group.
scram_long <- scram %>%
  dplyr::select(lib, inferred_group, naive_join_group, discordant) %>%
  tidyr::pivot_longer(c(inferred_group, naive_join_group),
                      names_to = "source", values_to = "assigned_group")
scram_long$source <- factor(scram_long$source,
  levels = c("inferred_group", "naive_join_group"),
  labels = c("Marker-derived (CPM column order, temp-major)",
             "Naive positional GSM->column join (genotype-major)"))

group_levels <- c("WT_37", "cGASKO_37", "WT_39", "cGASKO_39")
group_fill <- c("WT_37" = "#92C5DE", "cGASKO_37" = "#2166AC",
                "WT_39" = "#F4A582", "cGASKO_39" = "#B2182B")
scram_long$assigned_group <- factor(scram_long$assigned_group, levels = group_levels)

# Per-library discordance marker for the x-axis strip.
disc_strip <- data.frame(lib = scram$lib, discordant = scram$discordant)

fig1d <- ggplot(scram_long, aes(x = lib, y = source, fill = assigned_group)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_tile(data = subset(scram_long, discordant),
            aes(x = lib, y = source), fill = NA, color = "black", linewidth = 0.9) +
  geom_point(data = disc_strip, inherit.aes = FALSE,
             aes(x = lib, y = 0.4, shape = discordant),
             size = 2, color = "black") +
  scale_shape_manual(values = c(`TRUE` = 4, `FALSE` = NA),
                     name = "Discordant", labels = c(`TRUE` = "mislabeled", `FALSE` = "")) +
  scale_fill_manual(values = group_fill, name = "Assigned condition", drop = FALSE) +
  labs(
    title = "Fig 1d -- Scramble exhibit: positional GSM join mislabels libraries 026-035",
    subtitle = sprintf("%d/20 libraries discordant (black outline). CPM columns are temperature-major; GSM accessions are genotype-major.", n_disc),
    x = "Library (12630-RS-)", y = NULL
  ) +
  project_theme(config = FIG_CFG)

# Widest panel: 2-row tile x 20-library x-axis. Give the print variant generous
# width + a low height so the 20 rotated labels and the 2 tile rows stay legible.
save_figure(fig1d, "01_qc", "fig1d_scramble",
            variant = "both", config = FIG_CFG, width = 13, height = 4.8)
message("[FIG] saved fig1d_scramble (print+screen)")

write_caption(
  stage      = "01_qc",
  filename   = "figures/fig1d_scramble.png",
  finding    = sprintf(paste0(
    "Scramble exhibit (mapping-competence demonstration): a 2-row tile per library comparing the ",
    "marker-derived condition assignment (top row) against a naive positional GSM->column join ",
    "(bottom row), all 20 libraries on the x-axis. %d/20 libraries (026-035) are DISCORDANT ",
    "(black outline + x marker): the deposited CPM columns are temperature-major while the GEO GSM ",
    "accessions are genotype-major, so a positional join silently swaps KO-37 <-> WT-39 in the ",
    "middle block. This is why the mapping MUST be marker-derived, not accession-positional. ",
    "Claim tier: PROVISIONAL (pending collaborator sample sheet)."),
    n_disc),
  script     = SCRIPT,
  fn         = "geom_tile + geom_point (project_theme + save_figure)",
  config_kv  = CFG_KV,
  input      = "03_results/01_qc/tables/fig1d_scramble_data.csv",
  how_to_read = paste0(
    "x = library (12630-RS-021..040); y = the two label SOURCES (top = marker-derived, ",
    "bottom = naive positional GSM->column join). Tile FILL = assigned condition ",
    "(WT_37/cGASKO_37/WT_39/cGASKO_39, a blue->red diverging-by-condition map). ",
    "A black tile outline + an x glyph on the strip = the two sources DISAGREE for that library ",
    "(a mislabeled column under the naive join). Read: the discordant block (026-035) is exactly ",
    "where the temperature-major vs genotype-major orderings cross. Claim tier: PROVISIONAL."),
  config     = FIG_CFG
)

# -----------------------------------------------------------------------------
# E) Final structural assert (dual variants present; LOUD, no tryCatch)
# -----------------------------------------------------------------------------
expected <- c(
  file.path("03_results", "01_qc", "figures", "fig1a_thermometer.png"),
  file.path("03_results", "01_qc", "figures", "fig1a_thermometer.pdf"),
  file.path("03_results", "01_qc", "figures", "fig1b_cgas.png"),
  file.path("03_results", "01_qc", "figures", "fig1b_cgas.pdf"),
  file.path("03_results", "01_qc", "figures", "fig1d_scramble.png"),
  file.path("03_results", "01_qc", "figures", "fig1d_scramble.pdf"),
  file.path("03_results", "02_eda", "figures", "fig1c_pca_2x2.png"),
  file.path("03_results", "02_eda", "figures", "fig1c_pca_2x2.pdf"),
  file.path("03_results", "02_eda", "figures", "fig1c_pca_scree.png"),
  file.path("03_results", "02_eda", "figures", "fig1c_pca_scree.pdf")
)
missing <- expected[!file.exists(expected)]
if (length(missing) > 0) {
  stop(sprintf("[01_viz] %d expected figure variant(s) missing:\n  %s",
               length(missing), paste(missing, collapse = "\n  ")))
}

message(sprintf(
  "\n[DONE] 01_mapping_qc_viz.R complete: 5 figure stems, dual print+screen variants.\n  01_qc:  fig1a_thermometer, fig1b_cgas, fig1d_scramble\n  02_eda: fig1c_pca_2x2, fig1c_pca_scree (scree panel added; viz-only from existing varexp table)"
))
