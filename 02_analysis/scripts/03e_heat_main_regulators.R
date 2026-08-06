#!/usr/bin/env Rscript
# =============================================================================
# 03e_heat_main_regulators.R - PHASE 3 (stage 04_tf) COMPUTE: the HSF1 gap.
#   ULM-native heat-MAIN TF activity for the heat-shock / HIF / IFN axes.
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# Role:  COMPUTE half of the "normalize-then-visualize" split for fig3n. Runs run_ulm (the
#        one approved new statistic for this figure) on the heat-MAIN t-stat matrix
#        {Temp_main, WT_heat, KO_heat} against the cached CollecTRI net, using the same
#        estimator and call as 03_decoupler_tf.R. fig3n therefore sits on the same
#        decoupleR-ULM axis as fig3a/fig3c. Selects the heat-shock / HIF / IFN axis TFs and
#        emits one tidy plot-ready CSV.
#
# WHY ULM-NATIVE HERE (BEAT 5 / fig3n): the design spec recommended reusing the GSEA-NES
#   from master_tf_activities.csv to avoid a new statistic. Per Anton's instruction the
#   figure is built ULM-native so the deck stays on one estimator, and run_ulm is the
#   approved new statistic.
#
# Inputs:
#   - 03_results/objects/02_de_results.rds        (per-gene limma t-stat per contrast)
#   - 03_results/objects/net_collectri_mouse.rds  (source/target/mor; cached, since
#                                                   get_collectri() fails here)
#
# Output:
#   - 03_results/04_tf/tables/fig3n_heat_main_regulators_data.csv
#       cols: tf, contrast, score, padj, direction, axis
#
# Method: t-statistic input, per gsea.rank_metric; run_ulm(.mor="mor", minsize=5); BH within
#         contrast — identical to 03_decoupler_tf.R's primary call.
#
# Dependencies: decoupleR, dplyr, tidyr, readr
# =============================================================================

source("02_analysis/config/config.R")
load_packages(extra = c("decoupleR", "tidyr", "readr"))

set.seed(GSEA_SEED %||% 123)

TBL_DIR <- stage_dir("04_tf", "tables")

# Heat-MAIN contrasts only (the axis the literature says HSF1 dominates at 39C).
HEAT_CONTRASTS <- c("Temp_main", "WT_heat", "KO_heat")

# ----------------------------------------------------------------------------
# AXIS LOOKUP -- single source of truth (define ONCE here).
#   heatshock = Hsf1 (the gap we are closing)
#   HIF       = Hif1a / Epas1
#   IFN       = Irf3 / Stat1 / Stat2 / Nfkb1 / Rela
# Keys must match HEAT_AXIS_COLORS levels in config.R: heatshock / HIF / IFN / other
# ----------------------------------------------------------------------------
TF_AXIS <- c(
  Hsf1  = "heatshock",
  Hif1a = "HIF",   Epas1 = "HIF",
  Irf3  = "IFN",   Stat1 = "IFN", Stat2 = "IFN", Nfkb1 = "IFN", Rela = "IFN"
)
SELECT_TFS <- names(TF_AXIS)

cat("=================================================================\n")
cat("PHASE 3 (fig3n): ULM-native heat-MAIN regulators (the HSF1 gap)\n")
cat("=================================================================\n\n")

# =============================================================================
# 1. LOAD INPUTS + BUILD heat-MAIN t-STAT MATRIX (genes x 3 contrasts; NA -> 0)
# =============================================================================

cat("[1] Loading DE results + CollecTRI net ...\n")
de <- readRDS(file.path(DIR_OBJECTS, "02_de_results.rds"))
stopifnot(all(HEAT_CONTRASTS %in% names(de)))
de <- de[HEAT_CONTRASTS]                  # enforce column order

genes <- de[[1]]$gene_symbol
stopifnot(!any(duplicated(genes)))
mat <- sapply(de, function(d) d$t)        # t-statistic input (NEVER logFC)
rownames(mat) <- genes
mat[is.na(mat)] <- 0
cat(sprintf("  heat-MAIN t-stat matrix: %d genes x %d contrasts (%s)\n",
            nrow(mat), ncol(mat), paste(colnames(mat), collapse = ", ")))

net_ct <- readRDS(file.path(DIR_OBJECTS, "net_collectri_mouse.rds"))
stopifnot(all(SELECT_TFS %in% net_ct$source))
cat(sprintf("  CollecTRI net: %d edges, %d TFs (pre-minsize); all %d target TFs present\n\n",
            nrow(net_ct), length(unique(net_ct$source)), length(SELECT_TFS)))

# =============================================================================
# 2. run_ulm on heat-MAIN -- THE SAME CALL AS 03_decoupler_tf.R, BH within contrast
# =============================================================================

cat("[2] CollecTRI ULM on heat-MAIN contrasts + BH within contrast ...\n")
tf_ulm <- run_ulm(mat = mat, net = net_ct,
                  .source = "source", .target = "target", .mor = "mor",
                  minsize = 5)
tf_ulm <- tf_ulm %>%
  filter(condition %in% HEAT_CONTRASTS) %>%
  group_by(condition) %>%
  mutate(padj = p.adjust(p_value, "BH")) %>%
  ungroup()
cat(sprintf("  ULM TFs surviving minsize=5: %d; score range [%.2f, %.2f]\n\n",
            length(unique(tf_ulm$source)), min(tf_ulm$score), max(tf_ulm$score)))

# =============================================================================
# 3. SELECT the axis TFs + assemble the tidy table
# =============================================================================

cat("[3] Selecting heat-shock / HIF / IFN axis TFs ...\n")
fig3n <- tf_ulm %>%
  filter(source %in% SELECT_TFS) %>%
  transmute(
    tf        = source,
    contrast  = condition,
    score,
    padj,
    direction = ifelse(score > 0, "Up", "Down"),
    axis      = unname(TF_AXIS[source])
  ) %>%
  mutate(
    contrast = factor(contrast, levels = HEAT_CONTRASTS),
    axis     = factor(axis, levels = c("heatshock", "HIF", "IFN", "other"))
  ) %>%
  arrange(contrast, axis, desc(score))

stopifnot(nrow(fig3n) == length(SELECT_TFS) * length(HEAT_CONTRASTS))

# ---- HONESTY REPORT: does ULM agree with the NES co-elevation story? --------
cat("\n  --- HSF1 GAP: ULM scores + padj on heat-MAIN (Hsf1 / Hif1a / Epas1) ---\n")
report <- fig3n %>%
  filter(tf %in% c("Hsf1", "Hif1a", "Epas1")) %>%
  select(tf, contrast, score, padj, direction) %>%
  arrange(contrast, match(tf, c("Hsf1", "Hif1a", "Epas1")))
print(as.data.frame(report), row.names = FALSE, digits = 4)

hsf1_main <- fig3n %>% filter(tf == "Hsf1", contrast == "Temp_main")
hsf1_wt   <- fig3n %>% filter(tf == "Hsf1", contrast == "WT_heat")
co_elev_main <- hsf1_main$score > 0 && hsf1_main$padj < 0.05
co_elev_wt   <- hsf1_wt$score   > 0 && hsf1_wt$padj   < 0.05
cat(sprintf(
  "\n  VERDICT (ULM, heat-MAIN): Hsf1 positive & significant?  Temp_main: %s (score %.2f, padj %.3g) | WT_heat: %s (score %.2f, padj %.3g)\n",
  co_elev_main, hsf1_main$score, hsf1_main$padj,
  co_elev_wt,   hsf1_wt$score,   hsf1_wt$padj))
cat(sprintf("  => ULM %s the NES co-elevation story (Hsf1 co-elevated with HIF on heat-MAIN).\n\n",
            ifelse(co_elev_main, "AGREES WITH", "DISAGREES WITH")))

# =============================================================================
# 4. WRITE the tidy plot-ready table
# =============================================================================

cat("[4] Writing tidy table -> ", TBL_DIR, " ...\n", sep = "")
out <- file.path(TBL_DIR, "fig3n_heat_main_regulators_data.csv")
write.csv(fig3n, out, row.names = FALSE)
cat(sprintf("  [SAVE] %-42s %4d rows x %d cols\n",
            "fig3n_heat_main_regulators_data.csv", nrow(fig3n), ncol(fig3n)))

cat("\n[DONE] fig3n COMPUTE complete. Run 03e_heat_main_regulators_viz.R for the figure.\n")
