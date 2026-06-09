#!/usr/bin/env Rscript
# =============================================================================
# 01_mapping_qc_viz.R  --  PHASE 1: Sample-mapping QC figures (VISUALIZE)
# =============================================================================
# Phase:        1 (SCIENCE GATE -- visualization of the inferred sample mapping)
# Inputs:       03_results/01_qc/tables/fig1a_thermometer_data.csv
#               03_results/01_qc/tables/fig1b_cgas_data.csv
#               03_results/02_eda/tables/fig1c_pca_data.csv
#               03_results/02_eda/tables/fig1c_pca_varexp.csv
#               03_results/01_qc/tables/fig1d_scramble_data.csv
#               03_results/objects/01_eda.rds  (labels only: cgas_symbol)
# Outputs:      03_results/01_qc/figures/fig1a_thermometer.pdf
#               03_results/01_qc/figures/fig1b_cgas.pdf
#               03_results/02_eda/figures/fig1c_pca_2x2.pdf
#               03_results/01_qc/figures/fig1d_scramble.pdf
# Dependencies: config.R; ggplot2, dplyr, ggrepel, custom_minimal_theme.R
#
# NORMALIZE-THEN-VISUALIZE: this script performs NO statistics (no prcomp, no
# collapse, no thermometer/Cgas derivation). It reads the plot-ready tidy tables
# emitted by 01_mapping_qc.R and renders the 4 figures, doing only cosmetic
# reshaping (factor ordering, faceting). Run it AFTER 01_mapping_qc.R.
# =============================================================================

source("02_analysis/config/config.R")
load_packages(extra = c("tidyr"))

source("01_modules/RNAseq-toolkit/scripts/custom_minimal_theme.R")

# --- read the plot-ready tidy tables (no recomputation) ---
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

# Shared plotting scaffolding ------------------------------------------------
temp_cols <- c("37" = DIVERGING_COLORS$negative, "39" = DIVERGING_COLORS$positive)
geno_cols <- c("WT" = "#1B7837", "cGASKO" = "#762A83")
cap <- provisional_caption()

# Read variance-explained + n_top (labels only) from the varexp table.
percentVar <- setNames(fig1c_varexp$pct_var, fig1c_varexp$PC)
n_top <- fig1c_varexp$n_top[1]

# =============================================================================
# B) FIG 1a -- THERMOMETER (label-blind temperature validation)
# =============================================================================
message("\n==================  B) FIG 1a THERMOMETER  ==================")

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
    x = "Library (12630-RS-)", y = "log2(CPM + 0.5)",
    caption = cap
  ) +
  custom_minimal_theme_with_grid() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        plot.caption = element_text(face = "italic", size = 7))

ggsave(file.path(stage_dir("01_qc", "figures"), "fig1a_thermometer.pdf"),
       fig1a, width = 9, height = 8)
message("[FIG] saved fig1a_thermometer.pdf")

# =============================================================================
# C) FIG 1b -- CGAS GENOTYPE CHECK
# =============================================================================
message("\n==================  C) FIG 1b CGAS  ==================")

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
    x = "Library (12630-RS-)", y = "log2(CPM + 0.5)",
    caption = cap
  ) +
  custom_minimal_theme_with_grid() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        plot.caption = element_text(face = "italic", size = 7))

ggsave(file.path(stage_dir("01_qc", "figures"), "fig1b_cgas.pdf"),
       fig1b, width = 8, height = 5)
message("[FIG] saved fig1b_cgas.pdf")

# =============================================================================
# D) FIG 1c -- PCA (label-blind, overlay inferred labels)
# =============================================================================
message("\n==================  D) FIG 1c PCA 2x2  ==================")

pca_df <- fig1c_data
pca_df$temp     <- as.character(pca_df$temp)
pca_df$genotype <- as.character(pca_df$genotype)

fig1c <- ggplot(pca_df, aes(x = PC1, y = PC2, color = temp, shape = genotype)) +
  geom_point(size = 4, stroke = 1.1) +
  ggrepel::geom_text_repel(aes(label = lib), size = 2.6, color = "grey30",
                           max.overlaps = 20, show.legend = FALSE) +
  scale_color_manual(values = temp_cols, name = "Inferred temp (C)") +
  scale_shape_manual(values = c("WT" = 16, "cGASKO" = 17), name = "Inferred genotype") +
  labs(
    title = sprintf("Fig 1c -- PCA on top %d variable genes (label-blind, labels overlaid)", n_top),
    subtitle = "Expect PC1/PC2 to form the 2x2: temperature (primary) and genotype (secondary)",
    x = sprintf("PC1: %.1f%% var", percentVar[["PC1"]]),
    y = sprintf("PC2: %.1f%% var", percentVar[["PC2"]]),
    caption = cap
  ) +
  custom_minimal_theme_with_grid() +
  theme(plot.caption = element_text(face = "italic", size = 7))

ggsave(file.path(stage_dir("02_eda", "figures"), "fig1c_pca_2x2.pdf"),
       fig1c, width = 8, height = 6.5)
message("[FIG] saved fig1c_pca_2x2.pdf")

# =============================================================================
# E) FIG 1d -- SCRAMBLE EXHIBIT (competence exhibit)
# =============================================================================
message("\n==================  E) FIG 1d SCRAMBLE  ==================")

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
    x = "Library (12630-RS-)", y = NULL,
    caption = paste0(cap,
      "\nThe deposited CPM column order is temperature-major (proven by the thermometer: 031-040 hot); the GEO GSM accessions are genotype-major.",
      "\nA positional GSM->column join therefore mislabels the middle block (026-035): KO-37 columns get stamped WT-39 and vice versa.")
  ) +
  custom_minimal_theme_with_grid() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
        plot.caption = element_text(face = "italic", size = 6.5),
        legend.position = "right")

ggsave(file.path(stage_dir("01_qc", "figures"), "fig1d_scramble.pdf"),
       fig1d, width = 11, height = 5)
message("[FIG] saved fig1d_scramble.pdf")

message("\n[DONE] 01_mapping_qc_viz.R complete (4 figures rendered from tidy tables).")
