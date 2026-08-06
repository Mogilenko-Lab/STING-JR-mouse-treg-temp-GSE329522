#!/usr/bin/env Rscript
# =============================================================================
# 03g_nonidentifiability.R - PHASE 3 (stage 04_tf) COMPUTE: the non-identifiability
#   triptych. Why Hif1a's high heat-MAIN "activity" is unidentifiable.
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# Role:  COMPUTE half of the "normalize-then-visualize" split for fig3p/3q/3r. Runs run_ulm
#        (the one approved statistic, exactly as 03e_heat_main_regulators.R runs it — same
#        .mor="mor", minsize=5 call) on the heat-MAIN t-stat matrix, plus pure
#        set-membership COUNTING over CollecTRI regulons. Emits three tidy plot-ready CSVs.
#
# The honest message (leave every co-regulator and HIF2a uncrowned):
#   On heat-MAIN, Hif1a's "activity" is non-identifiable. Its target set is dominated by
#   genes that the network's most PROMISCUOUS TFs also regulate, so (a) Hif1a ranks #9 among
#   co-elevated stress/IEG/NF-kB TFs and (b) those same promiscuous TFs are the ones sharing
#   its targets. The claim is NON-IDENTIFIABILITY.
#
# Inputs (all cached; get_collectri() is broken here, so no network calls):
#   - 03_results/objects/02_de_results.rds        (per-contrast limma t-stats)
#   - 03_results/objects/net_collectri_mouse.rds  (source/target/mor)
#   - 03_results/04_tf/tables/fig3e_mlm_collinearity_data.csv  (353 Hif1a targets:
#                                                   target, t_wt, n_other_TFs)
#   - 03_results/04_tf/tables/fig3g_target_decomposition_data.csv  (signed contrib)
#
# Outputs (03_results/04_tf/tables/):
#   - fig3p_heatmain_ranking_data.csv   (rank, tf, score, p_value, family,
#                                        is_hif, is_hsf)  -- FULL ranking
#   - fig3q_coregulators_data.csv       (tf, shared_targets, pct_of_hif1a_set,
#                                        family)  -- sorted desc
#   - fig3r_membership_data.csv         (gene, tf, in_regulon, gene_heat_t,
#                                        tf_heatmain_score, gene_contrib, tf_family)
#
# Method: t-statistic input, per gsea.rank_metric; run_ulm(.mor="mor", minsize=5) on the
#   Temp_main column — the same estimator as 03e/03_decoupler_tf.R, so this sits on the same
#   ULM axis as the rest of the deck. Co-regulation is set-membership counting over
#   CollecTRI target lists.
#
# Dependencies: decoupleR, dplyr, tidyr, readr
# =============================================================================

source("02_analysis/config/config.R")
load_packages(extra = c("decoupleR", "tidyr", "readr"))

set.seed(GSEA_SEED %||% 123)

TBL_DIR <- stage_dir("04_tf", "tables")

HEAT_MAIN <- "Temp_main"   # the heat-MAIN contrast (name contains "Temp")

# ----------------------------------------------------------------------------
# TF FAMILY / ROLE LOOKUP -- curated, single source of truth (documented here).
#   This drives BOTH the family coloring across fig3p/3q/3r AND the fig3q "family"
#   column. The point of the coloring is to show that Hif1a's co-elevated peers
#   and target-sharers are the network's generic stress / immediate-early /
#   NF-kB / housekeeping regulators -- a broad stress cohort.
# Buckets (anything not listed -> "other"):
#   housekeeping            : Sp1, Sp3
#   stress/senescence       : Trp53
#   NF-kB                   : Nfkb, Nfkb1, Rela, Rel
#   AP-1 / immediate-early  : Jun, Jund, Fos, Fosb, Fosl2, Ap1, Egr1
#   proliferation/metabolism: Myc
#   nuclear receptor        : Esr1, Ar, Pparg, Ppara, Nr3c1, Nr1h4
#   cytokine/STAT           : Stat1, Stat3, Stat5a
#   growth-factor/other-stress: Creb1, Cebpb, Smad3, Foxo1, Foxo3, Ets1, Runx2,
#                               Ncoa1, Ep300
#   HIF axis                : Hif1a, Epas1
#   heat-shock              : Hsf1, Hsf2
# ----------------------------------------------------------------------------
FAMILY_LOOKUP <- list(
  "housekeeping"               = c("Sp1", "Sp3"),
  "stress/senescence"          = c("Trp53"),
  "NF-kB"                      = c("Nfkb", "Nfkb1", "Rela", "Rel"),
  "AP-1/immediate-early"       = c("Jun", "Jund", "Fos", "Fosb", "Fosl2", "Ap1", "Egr1"),
  "proliferation/metabolism"   = c("Myc"),
  "nuclear receptor"           = c("Esr1", "Ar", "Pparg", "Ppara", "Nr3c1", "Nr1h4"),
  "cytokine/STAT"              = c("Stat1", "Stat3", "Stat5a"),
  "growth-factor/other-stress" = c("Creb1", "Cebpb", "Smad3", "Foxo1", "Foxo3",
                                   "Ets1", "Runx2", "Ncoa1", "Ep300"),
  "HIF axis"                   = c("Hif1a", "Epas1"),
  "heat-shock"                 = c("Hsf1", "Hsf2")
)
fam_of <- function(tf) {
  for (f in names(FAMILY_LOOKUP)) if (tf %in% FAMILY_LOOKUP[[f]]) return(f)
  "other"
}
family_vec <- function(tfs) vapply(tfs, fam_of, character(1), USE.NAMES = FALSE)

cat("=================================================================\n")
cat("PHASE 3 (fig3p/3q/3r): non-identifiability of Hif1a on heat-MAIN\n")
cat("=================================================================\n\n")

# =============================================================================
# 1. LOAD INPUTS
# =============================================================================
cat("[1] Loading inputs ...\n")
de <- readRDS(file.path(DIR_OBJECTS, "02_de_results.rds"))
stopifnot(HEAT_MAIN %in% names(de))
net_ct <- readRDS(file.path(DIR_OBJECTS, "net_collectri_mouse.rds"))

fig3e <- read.csv(file.path(TBL_DIR, "fig3e_mlm_collinearity_data.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
fig3g <- read.csv(file.path(TBL_DIR, "fig3g_target_decomposition_data.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
cat(sprintf("  DE contrasts: %d | CollecTRI edges: %d (%d TFs) | fig3e targets: %d | fig3g rows: %d\n\n",
            length(de), nrow(net_ct), length(unique(net_ct$source)),
            nrow(fig3e), nrow(fig3g)))

# Heat-MAIN gene t-vector (named) -- single source for gene_heat_t below.
de_main <- de[[HEAT_MAIN]]
heat_t <- setNames(de_main$t, de_main$gene_symbol)
heat_t[is.na(heat_t)] <- 0

# =============================================================================
# 2. STEP 1 -- run_ulm on heat-MAIN vs FULL CollecTRI (THE approved statistic)
#    EXACTLY as 03e: .mor="mor", minsize=5, t-stat input, NA->0.
# =============================================================================
cat("[2] run_ulm on heat-MAIN (Temp_main) vs full CollecTRI (.mor='mor', minsize=5) ...\n")
mat <- matrix(heat_t, ncol = 1, dimnames = list(names(heat_t), HEAT_MAIN))

tf_ulm <- run_ulm(mat = mat, net = net_ct,
                  .source = "source", .target = "target", .mor = "mor",
                  minsize = 5)
tf_ulm <- tf_ulm %>% filter(condition == HEAT_MAIN)

ranking <- tf_ulm %>%
  arrange(desc(score)) %>%
  mutate(
    rank   = row_number(),
    tf     = source,
    family = family_vec(source),
    is_hif = source == "Hif1a",
    is_hsf = source == "Hsf1"
  ) %>%
  select(rank, tf, score, p_value, family, is_hif, is_hsf)

cat(sprintf("  ULM TFs surviving minsize=5: %d\n", nrow(ranking)))
cat("  Top 10 by score:\n")
print(as.data.frame(head(ranking, 10)), row.names = FALSE, digits = 4)

# --- ASSERTIONS (verified numbers; STOP if they don't reproduce) -------------
hif_row <- ranking[ranking$tf == "Hif1a", ]
hsf_row <- ranking[ranking$tf == "Hsf1",  ]
cat(sprintf("\n  Hif1a: rank #%d, score %.3f, p=%.3g | Hsf1: rank #%d, score %.3f\n",
            hif_row$rank, hif_row$score, hif_row$p_value, hsf_row$rank, hsf_row$score))
stopifnot(hif_row$rank == 9)
stopifnot(abs(hif_row$score - 5.14) < 0.05)
stopifnot(hsf_row$rank == 50)
stopifnot(abs(hsf_row$score - 3.20) < 0.05)
cat("  [ASSERT OK] Hif1a #9 ~5.14 ; Hsf1 #50 ~3.20\n")

# =============================================================================
# 3. STEP 2 -- co-regulator overlap (pure set-membership counting)
#    The 353 Hif1a targets come from fig3e; for every other TF count how many of
#    those 353 it also regulates in CollecTRI. NO statistic -- just counting.
# =============================================================================
cat("\n[3] Co-regulator overlap over the 353 Hif1a targets ...\n")
hif_targets <- fig3e$target
n_set <- length(hif_targets)
stopifnot(n_set == 353)

# target -> set of TFs (from CollecTRI), then per-other-TF shared count.
net_sub <- net_ct[net_ct$target %in% hif_targets, c("source", "target")]
net_sub <- unique(net_sub)

coreg <- net_sub %>%
  filter(source != "Hif1a") %>%
  group_by(source) %>%
  summarise(shared_targets = n_distinct(target), .groups = "drop") %>%
  transmute(
    tf               = source,
    shared_targets,
    pct_of_hif1a_set = shared_targets / n_set * 100,
    family           = family_vec(source)
  ) %>%
  arrange(desc(shared_targets), tf)

# Re-derive + assert the 92.1% / 325 / mean-22.2 (independently of fig3e cols).
n_other_by_target <- net_sub %>%
  filter(source != "Hif1a") %>%
  group_by(target) %>%
  summarise(n_other = n_distinct(source), .groups = "drop")
# targets with NO other TF contribute 0; align to full 353.
n_other_full <- setNames(rep(0L, n_set), hif_targets)
n_other_full[n_other_by_target$target] <- n_other_by_target$n_other
n_coreg   <- sum(n_other_full >= 1)
pct_coreg <- n_coreg / n_set * 100
mean_other <- mean(n_other_full)

cat(sprintf("  co-regulated (>=1 other TF): %d/%d = %.1f%% | mean other-TFs/target = %.1f\n",
            n_coreg, n_set, pct_coreg, mean_other))
cat("  Top co-regulators by shared count:\n")
print(as.data.frame(head(coreg, 12)), row.names = FALSE, digits = 4)

stopifnot(n_coreg == 325)
stopifnot(abs(pct_coreg - 92.1) < 0.2)
stopifnot(abs(mean_other - 22.2) < 0.2)
cat("  [ASSERT OK] 325/353 = 92.1% co-regulated ; mean other-TFs/target ~22.2\n")

# Sanity vs fig3e's own n_other_TFs column (should agree).
if ("n_other_TFs" %in% names(fig3e)) {
  agree <- all.equal(unname(n_other_full[fig3e$target]), as.integer(fig3e$n_other_TFs))
  cat(sprintf("  cross-check vs fig3e$n_other_TFs: %s\n",
              ifelse(isTRUE(agree), "MATCH", paste("DIFF:", agree))))
}

# =============================================================================
# 4. STEP 3 -- unifier membership matrix for fig3r (long tidy table)
#    rows = top ~20 Hif1a targets by POSITIVE signed contrib (fig3g);
#    cols = Hif1a + top ~12 co-regulators by shared count (step 2).
#    Each (gene, tf): in_regulon = is gene in that TF's CollecTRI regulon?
# =============================================================================
cat("\n[4] Building fig3r membership matrix (top-20 heat-driving genes x Hif1a+12 co-regs) ...\n")
N_GENES <- 20
N_COREG <- 12

top_genes <- fig3g %>%
  arrange(desc(contrib)) %>%
  slice_head(n = N_GENES)
gene_set <- top_genes$target

top_coreg <- head(coreg$tf, N_COREG)
tf_cols   <- c("Hif1a", top_coreg)   # Hif1a first; viz will order cols by score

# gene contrib lookup (Hif1a signed contrib, from fig3g)
contrib_of <- setNames(fig3g$contrib, fig3g$target)
# TF heat-MAIN ULM score lookup (from step 1 ranking; TFs not scored -> NA)
score_of   <- setNames(ranking$score, ranking$tf)

membership <- expand.grid(gene = gene_set, tf = tf_cols,
                          stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(
    in_regulon        = as.integer(any(net_ct$source == tf & net_ct$target == gene)),
    gene_heat_t       = unname(heat_t[gene]),
    tf_heatmain_score = unname(score_of[tf]),
    gene_contrib      = unname(contrib_of[gene]),
    tf_family         = fam_of(tf)
  ) %>%
  ungroup()

cat(sprintf("  membership rows: %d (%d genes x %d TFs); Hif1a-own coverage: %d/%d genes\n",
            nrow(membership), length(gene_set), length(tf_cols),
            sum(membership$in_regulon[membership$tf == "Hif1a"]), length(gene_set)))

# =============================================================================
# 5. WRITE the three tidy tables
# =============================================================================
cat("\n[5] Writing tidy tables -> ", TBL_DIR, " ...\n", sep = "")
w <- function(df, fname) {
  write.csv(df, file.path(TBL_DIR, fname), row.names = FALSE)
  cat(sprintf("  [SAVE] %-38s %4d rows x %d cols\n", fname, nrow(df), ncol(df)))
}
w(ranking,    "fig3p_heatmain_ranking_data.csv")
w(coreg,      "fig3q_coregulators_data.csv")
w(membership, "fig3r_membership_data.csv")

cat("\n[DONE] 03g COMPUTE complete. Run 03g_nonidentifiability_viz.R for fig3p/3q/3r.\n")
