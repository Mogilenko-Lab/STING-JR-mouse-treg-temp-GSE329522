#!/usr/bin/env Rscript
# =============================================================================
# 03_decoupler_tf.R - PHASE 3 (stage 04_tf) COMPUTE: decoupleR TF activity +
#                     HIF1a cross-method rank-cascade forensics
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# Role:    COMPUTE half of the "normalize-then-visualize" split. Runs all statistics and
#          writes the forensics checkpoint, the master accumulator, and one plot-ready tidy
#          CSV per downstream figure into stage_dir("04_tf","tables"). Figures live in
#          03_decoupler_tf_viz.R, which reads only the tidy tables emitted here.
#
# Inputs:
#   - 03_results/objects/02_de_results.rds          (7 contrasts; t-stat per gene)
#   - 03_results/objects/net_collectri_mouse.rds    (primary TF net; source/target/mor)
#   - 03_results/objects/net_dorothea_mouse_ABC.rds (DoRothEA A/B/C; the net Phylo used)
#   - 00_data/references/gene_sets/lombardi2022_hif_consensus_mouse.rds
#       ($Lombardi2022_HIF = mouse HIF consensus; length read at RUNTIME)
#
# Outputs:
#   - 03_results/objects/03_tf_collectri.rds        (CollecTRI ULM, all 7 contrasts, BH)
#   - 03_results/objects/03_tf_forensics.rds        (forensics bundle: all prior blocks
#                                                    + regulon_swap, mlm_collinearity,
#                                                    consensus_zmix, top_tf_by_contrast,
#                                                    refreshed lombardi_vs_phylo)
#   - 03_results/master/master_tf_activities.csv    (13-col; CollecTRI all + DoRothEA WT_heat)
#   - 03_results/04_tf/tables/fig3*_*_data.csv       (one tidy CSV per downstream figure)
#
# Method: t-statistic input, per gsea.rank_metric; run_ulm(.mor="mor", minsize=5); BH within
#         contrast; consensus via decouple(...consensus_score=TRUE).
#
# NETWORK NOTE: decoupleR::get_collectri()/get_progeny() fail here (OmnipathR fallback
#   errors), so this stage loads the CACHED nets built by 00c_prepare_networks.R.
#
# Dependencies: decoupleR, dplyr, tidyr, readr, tibble, purrr
# =============================================================================

source("02_analysis/config/config.R")
load_packages(extra = c("decoupleR", "tidyr", "readr", "tibble", "purrr"))

set.seed(GSEA_SEED %||% 123)

TBL_DIR <- stage_dir("04_tf", "tables")   # all tidy CSVs land here

CONTRASTS <- c("WT_heat", "KO_heat", "Interaction",
               "Geno_at_39", "Geno_at_37", "Temp_main", "Geno_main")

# IFN/NFkB axis vs HIF axis (KEY_TFS from config = the union we report)
HIF_AXIS  <- c("Hif1a", "Epas1")
IFN_AXIS  <- c("Irf1", "Irf3", "Irf7", "Stat1", "Stat2", "Nfkb1", "Rela")

# How many TFs to surface per direction in the contrast-illustration figure.
TOPN_CONTRAST <- 12

cat("=================================================================\n")
cat("PHASE 3: decoupleR TF activity + HIF1a rank-cascade forensics\n")
cat("=================================================================\n\n")

# =============================================================================
# 1. LOAD INPUTS + BUILD t-STAT MATRIX (genes x 7 contrasts; NA -> 0)
# =============================================================================

cat("[1] Loading DE results + networks ...\n")
de <- readRDS(file.path(DIR_OBJECTS, "02_de_results.rds"))
stopifnot(all(CONTRASTS %in% names(de)))
de <- de[CONTRASTS]                       # enforce column order

genes <- de[[1]]$gene_symbol
stopifnot(!any(duplicated(genes)))        # collapse decided in Phase 1; verify here
mat <- sapply(de, function(d) d$t)
rownames(mat) <- genes
mat[is.na(mat)] <- 0
universe <- rownames(mat)
cat(sprintf("  t-stat matrix: %d genes x %d contrasts\n", nrow(mat), ncol(mat)))

net_ct  <- readRDS(file.path(DIR_OBJECTS, "net_collectri_mouse.rds"))
net_dor <- readRDS(file.path(DIR_OBJECTS, "net_dorothea_mouse_ABC.rds"))
cat(sprintf("  CollecTRI net: %d edges, %d TFs (pre-minsize)\n",
            nrow(net_ct), length(unique(net_ct$source))))
cat(sprintf("  DoRothEA ABC net: %d edges, %d TFs (pre-minsize)\n",
            nrow(net_dor), length(unique(net_dor$source))))

# NEW Lombardi mouse set (Agent 1 re-derivation). Read length at RUNTIME.
lombardi <- readRDS("00_data/references/gene_sets/lombardi2022_hif_consensus_mouse.rds")
lombardi_genes <- lombardi[["Lombardi2022_HIF"]]
n_lombardi_full <- length(lombardi_genes)
lombardi_prov <- attr(lombardi_genes, "provenance")
cat(sprintf("  Lombardi-2022 HIF consensus (NEW set): %d mouse genes\n", n_lombardi_full))
if (!is.null(lombardi_prov)) cat(sprintf("    provenance: %s\n", paste(utils::head(lombardi_prov, 3), collapse=" | ")))
cat("\n")

# helper: BH within contrast
bh_within <- function(df) {
  df %>% group_by(condition) %>%
    mutate(padj = p.adjust(p_value, "BH")) %>% ungroup()
}

# =============================================================================
# 2. PRIMARY: CollecTRI ULM on all 7 contrasts, BH within contrast
# =============================================================================

cat("[2] CollecTRI ULM (all 7 contrasts) + BH within contrast ...\n")
tf_ct <- load_or_compute(
  "03_tf_collectri.rds",
  desc = "CollecTRI ULM TF activities (7 contrasts, BH-within-contrast)",
  compute_fn = function() {
    res <- run_ulm(mat = mat, net = net_ct,
                   .source = "source", .target = "target", .mor = "mor",
                   minsize = 5)
    res <- res %>% filter(condition %in% CONTRASTS)
    bh_within(res)
  }
)
n_tf_ct <- length(unique(tf_ct$source))
cat(sprintf("  CollecTRI TFs surviving minsize=5: %d  (acceptance: 100-1500)\n", n_tf_ct))
cat(sprintf("  Activity score range: [%.2f, %.2f]; mean ~ %.3f (acceptance: centered ~0)\n",
            min(tf_ct$score), max(tf_ct$score), mean(tf_ct$score)))
cat(sprintf("  BH applied? padj != p_value for some rows: %s\n\n",
            any(abs(tf_ct$padj - tf_ct$p_value) > 1e-12)))

# ---- THE KEY TEST: Interaction contrast (cGAS-dependence) -------------------
key_test <- tf_ct %>%
  filter(condition == "Interaction", source %in% c(HIF_AXIS, IFN_AXIS)) %>%
  mutate(axis = ifelse(source %in% HIF_AXIS, "HIF", "IFN/NFkB")) %>%
  arrange(axis, p_value) %>%
  select(axis, TF = source, score, p_value, padj)
cat("  --- KEY TEST: TF activity on Interaction (HIF axis vs IFN/NFkB axis) ---\n")
print(as.data.frame(key_test), row.names = FALSE)

# gene-level Interaction cross-reference for ISG markers (regulon-activity dilution check)
isg_gene_int <- de[["Interaction"]] %>%
  filter(gene_symbol %in% ISG_MARKERS) %>%
  arrange(adj.P.Val) %>%
  select(gene = gene_symbol, logFC, t, P.Value, adj.P.Val)
cat("\n  --- gene-level Interaction stats for ISG_MARKERS (corroborates IFN cGAS-dependence) ---\n")
print(as.data.frame(isg_gene_int), row.names = FALSE)

# ---- HIF1a heat response in BOTH genotypes ----------------------------------
hif_heat <- tf_ct %>%
  filter(source == "Hif1a", condition %in% c("WT_heat", "KO_heat", "Temp_main")) %>%
  select(condition, score, p_value, padj)
cat("\n  --- HIF1a TF activity across heat contrasts (driven by heat in BOTH genotypes) ---\n")
print(as.data.frame(hif_heat), row.names = FALSE)
cat("\n")

# =============================================================================
# 3. HIF1a RANK-CASCADE FORENSICS  (#1 -> #12 -> #142 -> #8)
# =============================================================================

cat("[3] HIF1a rank-cascade forensics ...\n")

wt_vec <- mat[, "WT_heat", drop = FALSE]   # single-contrast input for forensics

## 3.1 Reproduce Phylo: DoRothEA ABC ULM on WT_heat -- Hif1a rank ------------
dor_wt <- run_ulm(mat = wt_vec, net = net_dor,
                  .source = "source", .target = "target", .mor = "mor",
                  minsize = 5) %>%
  arrange(desc(score)) %>%
  mutate(rank = row_number())
n_dor <- nrow(dor_wt)
hif_dor <- dor_wt %>% filter(source == "Hif1a")
cat(sprintf("  [3.1] DoRothEA WT_heat: Hif1a rank #%d of %d, score=%.2f, p=%.2g\n",
            hif_dor$rank, n_dor, hif_dor$score, hif_dor$p_value))

## 3.2 Network swap: CollecTRI ULM on WT_heat -- Hif1a rank -----------------
ct_wt <- tf_ct %>% filter(condition == "WT_heat") %>%
  arrange(desc(score)) %>% mutate(rank = row_number())
n_ct <- nrow(ct_wt)
hif_ct <- ct_wt %>% filter(source == "Hif1a")
cat(sprintf("  [3.2] NETWORK SWAP -- CollecTRI WT_heat: Hif1a rank #%d of %d, score=%.2f, p=%.2g\n",
            hif_ct$rank, n_ct, hif_ct$score, hif_ct$p_value))

## 3.3 De-confounding: MLM + consensus on WT_heat (CollecTRI) ---------------
# MLM is multivariate -> penalizes TFs whose targets are shared with other regulons.
mlm_wt <- run_mlm(mat = wt_vec, net = net_ct,
                  .source = "source", .target = "target", .mor = "mor",
                  minsize = 5) %>%
  arrange(desc(score)) %>% mutate(rank = row_number())
hif_mlm <- mlm_wt %>% filter(source == "Hif1a")

# decouple() with ulm + mlm + wsum (decoupleR also emits norm_wsum + corr_wsum),
# then consensus. KEEP the per-statistic rows -- block (C) consumes them.
dec <- decouple(mat = wt_vec, network = net_ct,
                .source = "source", .target = "target",
                statistics = c("ulm", "mlm", "wsum"),
                args = list(
                  ulm  = list(.mor = "mor", minsize = 5),
                  mlm  = list(.mor = "mor", minsize = 5),
                  wsum = list(.mor = "mor", minsize = 5)
                ),
                consensus_score = TRUE, minsize = 5, show_toy_call = FALSE)
cons_wt <- dec %>% filter(statistic == "consensus") %>%
  arrange(desc(score)) %>% mutate(rank = row_number())
hif_cons <- cons_wt %>% filter(source == "Hif1a")

cat(sprintf("  [3.3] DE-CONFOUNDING on WT_heat (CollecTRI):  ULM #%d (%.2f) -> MLM #%d (%.2f) -> consensus #%d (%.2f)\n",
            hif_ct$rank, hif_ct$score, hif_mlm$rank, hif_mlm$score, hif_cons$rank, hif_cons$score))

## 3.4 Target decomposition: HIF-specific vs shared-glycolytic vs other ------
hif_targets <- net_ct %>% filter(source == "Hif1a")
present <- hif_targets %>% filter(target %in% universe)
present$t_wt  <- mat[present$target, "WT_heat"]
present$contrib <- sign(present$mor) * present$t_wt    # signed t-stat * mor sign

# RECONCILED CLASS VOCABULARY (audit item 4 -- fig3g <-> fig3l must not contradict).
# fig3l's finer module decomposition shows the old coarse "HIF-specific" lump was
# misleading: its UP members (Vegfa/Slc2a1/Egln3) are shared HIF1/HIF2 angio-glucose +
# feedback targets (least HIF-diagnostic), while only the REPRESSED Pdk1/Bnip3/Bnip3l/Car9
# are the diagnostic hypoxic-HIF core. So fig3g's coarse 3-class scheme is now:
#   "hypoxic HIF core"  = Pdk1/Bnip3/Bnip3l/Car9 (REPRESSED; the HIF-diagnostic core)
#   "shared/glycolytic" = the UP HIF targets (Vegfa/Slc2a1/Egln3) + the generic glycolytic set
#   "other"             = everything else (the 340-member lump; fig3l carves heat-shock out of it)
# This is a RELABEL of class names + membership boundary only -- contrib values are unchanged
# (sign(mor) x t_wt per gene). "hypoxic HIF core" sum (-8.46) now equals fig3l's
# hif1a_hypoxic_core, so "the HIF core" denotes ONE set across both figures.
HYPOXIC_HIF_CORE <- c("Pdk1", "Bnip3", "Bnip3l", "Car9")  # diagnostic; repressed
SHARED_GLYCO     <- c("Vegfa", "Slc2a1", "Egln3",          # up HIF targets (shared/feedback)
                      "Ldha", "Pgk1", "Aldoa", "Eno1", "Pkm", "Hk2", "Tpi1")

present$class <- ifelse(present$target %in% HYPOXIC_HIF_CORE, "hypoxic HIF core",
                 ifelse(present$target %in% SHARED_GLYCO, "shared/glycolytic", "other"))

decomp <- present %>%
  group_by(class) %>%
  summarise(n_targets = n(),
            sum_contrib = sum(contrib),
            mean_contrib = mean(contrib), .groups = "drop") %>%
  arrange(desc(sum_contrib))
total_contrib <- sum(present$contrib)
decomp$pct_of_total <- 100 * decomp$sum_contrib / total_contrib

cat(sprintf("  [3.4] TARGET DECOMP (Hif1a regulon, %d present): hypoxic HIF core %.1f%% / shared+glyco %.1f%% / other %.1f%%\n",
            nrow(present),
            decomp$pct_of_total[decomp$class=="hypoxic HIF core"],
            decomp$pct_of_total[decomp$class=="shared/glycolytic"],
            decomp$pct_of_total[decomp$class=="other"]))

# =============================================================================
# 3M. MECHANISM BLOCKS -- WHY Hif1a goes #1 -> #12 -> #142 -> #8
# =============================================================================

cat("\n[3M] Mechanism blocks (the rank cascade) ...\n")

## (A) regulon_swap : #1 (DoRothEA-ULM) -> #12 (CollecTRI-ULM) --------------
# DoRothEA's smaller, curated Hif1a regulon is enriched for high-t heat targets;
# CollecTRI's larger regulon dilutes with low-t targets -> ULM score drops.
dor_tgt <- net_dor %>% filter(source == "Hif1a", target %in% universe) %>% pull(target) %>% unique()
ct_tgt  <- net_ct  %>% filter(source == "Hif1a", target %in% universe) %>% pull(target) %>% unique()
both_tgt <- union(dor_tgt, ct_tgt)

regulon_swap_pergene <- tibble(
  target       = both_tgt,
  in_dorothea  = both_tgt %in% dor_tgt,
  in_collectri = both_tgt %in% ct_tgt,
  t_wt         = mat[both_tgt, "WT_heat"]
) %>%
  mutate(membership = case_when(
    in_dorothea & in_collectri ~ "both",
    in_dorothea                ~ "dorothea_only",
    TRUE                       ~ "collectri_only"
  )) %>%
  arrange(desc(abs(t_wt)))

# per-network t_wt distribution summary (the dilution story, plot-ready)
net_t_summary <- function(net_name, tgts) {
  tv <- mat[tgts, "WT_heat"]
  tibble(network = net_name, n = length(tv),
         mean_abs_t = mean(abs(tv)), median_t = median(tv),
         pct_abs_t_gt2 = 100 * mean(abs(tv) > 2))
}
regulon_swap_summary <- bind_rows(
  net_t_summary("DoRothEA", dor_tgt),
  net_t_summary("CollecTRI", ct_tgt)
) %>%
  mutate(n_intersection = length(intersect(dor_tgt, ct_tgt)),
         n_dorothea_only = length(setdiff(dor_tgt, ct_tgt)),
         n_collectri_only = length(setdiff(ct_tgt, dor_tgt)))

regulon_swap <- list(per_gene = regulon_swap_pergene, summary = regulon_swap_summary)
cat(sprintf("  (A) regulon_swap: DoRothEA |%d| vs CollecTRI |%d|; intersect %d; DoRothEA mean|t|=%.2f vs CollecTRI %.2f\n",
            length(dor_tgt), length(ct_tgt), length(intersect(dor_tgt, ct_tgt)),
            regulon_swap_summary$mean_abs_t[regulon_swap_summary$network=="DoRothEA"],
            regulon_swap_summary$mean_abs_t[regulon_swap_summary$network=="CollecTRI"]))

## (B) mlm_collinearity : #12 (ULM) -> #142 (MLM) ---------------------------
# For each Hif1a CollecTRI target present, count OTHER CollecTRI TFs regulating it.
# Shared targets are de-confounded away under the multivariate MLM fit.
ct_present <- net_ct %>% filter(target %in% universe)
target_tf_count <- ct_present %>%
  distinct(source, target) %>%
  count(target, name = "n_tfs_total")          # incl. Hif1a itself
mlm_pergene <- tibble(target = ct_tgt) %>%
  left_join(target_tf_count, by = "target") %>%
  mutate(n_other_TFs = n_tfs_total - 1L,        # exclude Hif1a
         t_wt        = mat[target, "WT_heat"]) %>%
  select(target, n_other_TFs, t_wt) %>%
  arrange(desc(n_other_TFs))

mlm_collinearity_summary <- tibble(
  n_targets             = nrow(mlm_pergene),
  frac_co_regulated     = mean(mlm_pergene$n_other_TFs >= 1),
  mean_other_TFs        = mean(mlm_pergene$n_other_TFs),
  median_other_TFs      = median(mlm_pergene$n_other_TFs)
)
# ULM vs MLM Hif1a WT_heat score collapse (5.110 -> 1.134), 2-row table
score_collapse <- tibble(
  method = factor(c("ULM", "MLM"), levels = c("ULM", "MLM")),
  score  = c(hif_ct$score, hif_mlm$score),
  rank   = c(hif_ct$rank, hif_mlm$rank),
  n_tfs  = c(n_ct, nrow(mlm_wt))
)
mlm_collinearity <- list(per_gene = mlm_pergene,
                         summary = mlm_collinearity_summary,
                         score_collapse = score_collapse)
cat(sprintf("  (B) mlm_collinearity: %.0f%% of Hif1a targets co-regulated by >=1 other TF (mean %.1f others); ULM %.2f -> MLM %.2f\n",
            100 * mlm_collinearity_summary$frac_co_regulated,
            mlm_collinearity_summary$mean_other_TFs, hif_ct$score, hif_mlm$score))

## (C) consensus_zmix : #142 (MLM) -> #8 (consensus) ------------------------
# Reproduce decoupleR's consensus = per-statistic sign-symmetric ("folded")
# standardisation across all TFs, then unweighted mean across statistics.
# run_consensus() folds each statistic about 0 (mirror + scale by symmetric sd)
# then averages -- the two univariate families (ulm, wsum/norm_wsum/corr_wsum)
# sit high, mlm sits low, the unweighted mean re-inflates Hif1a.
nonc <- dec %>% filter(statistic != "consensus")
folded <- nonc %>%
  group_by(statistic, condition) %>% group_split() %>%
  map(function(d) {
    pos <- d %>% filter(score > 0) %>% rbind(., mutate(., score = -score)) %>%
      mutate(score = score / sd(score)) %>% filter(score > 0)
    neg <- d %>% filter(score <= 0) %>% rbind(., mutate(., score = -score)) %>%
      mutate(score = score / sd(score)) %>% filter(score <= 0)
    rbind(pos, neg)
  }) %>% bind_rows() %>%
  rename(z_score = score)

# Hif1a folded-z per statistic + reproduced consensus (assert exact)
hif_zrows <- folded %>% filter(source == "Hif1a") %>%
  select(statistic, z_score) %>% arrange(statistic)
hif_cons_repro <- mean(hif_zrows$z_score)
stopifnot(abs(hif_cons_repro - hif_cons$score) < 1e-6)   # folded-z mean reproduces consensus

# tidy table: ulm/mlm/wsum z for Hif1a + the consensus row (canonical 3 + consensus).
# Keep wsum variants too, flagged, so the viz can show univariate double-counting.
consensus_zmix <- bind_rows(
  hif_zrows %>% mutate(source = "Hif1a"),
  tibble(statistic = "consensus", z_score = hif_cons$score, source = "Hif1a")
) %>%
  mutate(
    family = case_when(
      statistic == "mlm"                                   ~ "multivariate",
      statistic == "consensus"                             ~ "consensus",
      TRUE                                                 ~ "univariate"
    ),
    canonical = statistic %in% c("ulm", "mlm", "wsum", "consensus")
  ) %>%
  select(source, statistic, z_score, family, canonical) %>%
  arrange(match(statistic, c("ulm", "mlm", "wsum", "norm_wsum", "corr_wsum", "consensus")))
cat(sprintf("  (C) consensus_zmix: Hif1a folded-z ulm=%.2f mlm=%.2f wsum=%.2f -> mean over %d stats = consensus %.3f (repro diff %.1e)\n",
            hif_zrows$z_score[hif_zrows$statistic=="ulm"],
            hif_zrows$z_score[hif_zrows$statistic=="mlm"],
            hif_zrows$z_score[hif_zrows$statistic=="wsum"],
            nrow(hif_zrows), hif_cons$score, abs(hif_cons_repro - hif_cons$score)))

## (D) top_tf_by_contrast : the contrast-illustration figure -----------------
axis_of <- function(tf) ifelse(tf %in% HIF_AXIS, "HIF",
                        ifelse(tf %in% IFN_AXIS, "IFN", "other"))
top_tf_by_contrast <- bind_rows(lapply(
  c("WT_heat", "KO_heat", "Temp_main", "Interaction"),
  function(cc) {
    d <- tf_ct %>% filter(condition == cc) %>% arrange(desc(score))
    up   <- d %>% head(TOPN_CONTRAST) %>% mutate(direction = "Up")
    down <- d %>% arrange(score) %>% head(TOPN_CONTRAST) %>% mutate(direction = "Down")
    bind_rows(up, down) %>%
      group_by(direction) %>%
      mutate(rank = ifelse(direction == "Up",
                           rank(-score, ties.method = "first"),
                           rank(score,  ties.method = "first"))) %>%
      ungroup() %>%
      transmute(contrast = cc, source, score, padj,
                direction, rank,
                axis = axis_of(source),
                key_tf = source %in% KEY_TFS)
  }
))
cat(sprintf("  (D) top_tf_by_contrast: %d rows (top/bottom %d x 4 contrasts), %d key-TF flags\n",
            nrow(top_tf_by_contrast), TOPN_CONTRAST, sum(top_tf_by_contrast$key_tf)))

# =============================================================================
# 4. LOMBARDI vs PHYLO-16 ORTHOGONAL HIF CHECK  (block E -- NEW gene set)
# =============================================================================

cat("\n[4] Lombardi vs Phylo-16 orthogonal HIF signature check (NEW Lombardi set) ...\n")

# Phylo's 16-gene HIF set mapped to mouse (Phylo-16 kept AS-IS).
phylo16 <- c("Hif1a", "Egln3", "Ldha", "Slc2a1", "Pkm", "Pgk1", "Eno1", "Bnip3",
             "Bnip3l", "Ddit4", "Vegfa", "Pdk1", "Pfkfb3", "Slc16a3", "Cd274", "Cxcr4")
phylo16_present   <- intersect(phylo16, universe)
lombardi_present  <- intersect(lombardi_genes, universe)
n_lombardi <- length(lombardi_present)   # dynamic count for figure labels
n_phylo    <- length(phylo16_present)
cat(sprintf("  Phylo-16 present in data: %d/%d\n", n_phylo, length(phylo16)))
cat(sprintf("  Lombardi-%d present in data: %d/%d\n", n_lombardi_full, n_lombardi, n_lombardi_full))

mk_net <- function(name, genes) tibble(source = name, target = genes, mor = 1)
sig_net <- bind_rows(
  mk_net("Lombardi2022_HIF", lombardi_present),
  mk_net("Phylo16_HIF",       phylo16_present)
)

sig_res <- run_ulm(mat = mat, net = sig_net,
                   .source = "source", .target = "target", .mor = "mor",
                   minsize = 5) %>%
  filter(condition %in% CONTRASTS) %>%
  bh_within()

sig_tab <- sig_res %>%
  filter(condition %in% c("WT_heat", "KO_heat", "Interaction")) %>%
  select(signature = source, condition, score, p_value, padj) %>%
  arrange(signature, factor(condition, levels = c("WT_heat", "KO_heat", "Interaction"))) %>%
  mutate(n_genes = ifelse(signature == "Lombardi2022_HIF", n_lombardi, n_phylo),
         sig_label = ifelse(signature == "Lombardi2022_HIF",
                            paste0("Lombardi-", n_lombardi),
                            paste0("Phylo-", n_phylo)),
         sig = padj < 0.05)

# Egln3 gene-level (HIF-autofeedback diagnostic)
egln3_gene <- bind_rows(lapply(c("WT_heat", "KO_heat", "Interaction"), function(cc) {
  de[[cc]] %>% filter(gene_symbol == "Egln3") %>%
    transmute(contrast = cc, logFC, t, P.Value, adj.P.Val)
}))

# ---- OLD vs NEW side-by-side (lombardi_vs_phylo IS expected to change) -------
old_forensics <- if (file.exists(file.path(DIR_OBJECTS, "03_tf_forensics.rds"))) {
  readRDS(file.path(DIR_OBJECTS, "03_tf_forensics.rds"))
} else NULL
cat("\n  --- lombardi_vs_phylo OLD vs NEW (Lombardi gene set changed; Phylo-16 unchanged) ---\n")
if (!is.null(old_forensics) && !is.null(old_forensics$lombardi_vs_phylo)) {
  old_lp <- old_forensics$lombardi_vs_phylo %>%
    transmute(signature, condition, score_OLD = score, padj_OLD = padj)
  cmp_lp <- sig_tab %>%
    transmute(signature, condition, score_NEW = score, padj_NEW = padj) %>%
    left_join(old_lp, by = c("signature", "condition")) %>%
    select(signature, condition, score_OLD, score_NEW, padj_OLD, padj_NEW)
  print(as.data.frame(cmp_lp), row.names = FALSE, digits = 4)
} else {
  cat("  (no prior forensics found; printing NEW only)\n")
  print(as.data.frame(sig_tab), row.names = FALSE, digits = 4)
}
# Does the conclusion hold? both signatures up with heat in BOTH genotypes; Interaction NS.
lp <- sig_tab
both_up_heat <- all(lp$score[lp$condition %in% c("WT_heat", "KO_heat")] > 0) &&
  all(lp$padj[lp$condition %in% c("WT_heat", "KO_heat")] < 0.05)
inter_ns <- all(lp$padj[lp$condition == "Interaction"] >= 0.05)
cat(sprintf("  CONCLUSION HOLDS? both sigs UP+sig with heat in both genotypes: %s | Interaction NS: %s\n",
            both_up_heat, inter_ns))

# =============================================================================
# 5. OUTPUTS: master_tf_activities.csv + forensics bundle
# =============================================================================

cat("\n[5] Writing master table + forensics bundle ...\n")

# regulon gene string (core_enrichment) per TF, top-50 by |mor|
regulon_str <- function(network, tf, universe) {
  network %>% filter(source == tf, target %in% universe) %>%
    arrange(desc(abs(mor))) %>% head(50) %>% pull(target) %>%
    paste(collapse = "/")
}
ct_set_size <- net_ct %>% filter(target %in% universe) %>%
  count(source, name = "set_size")

# COLUMN-NAME CAVEAT (do NOT misread): the master table reuses the legacy
# GSEA-style schema (column names nes/pvalue/core_enrichment) as a shared
# accumulator format, but the value written to `nes` here is the decoupleR ULM
# `score` -- the decoupleR ULM score. fig3n (03e_heat_main_regulators.R)
# runs the SAME run_ulm(.mor="mor", minsize=5) call, so its scores are byte-identical
# to this `nes` column and ARE axis-comparable to fig3a/fig3c. The `nes` name is a
# misnomer kept for master-schema compatibility; treat it as the ULM score everywhere.
master_ct <- tf_ct %>%
  left_join(ct_set_size, by = "source") %>%
  rowwise() %>%
  mutate(core_enrichment = regulon_str(net_ct, source, universe)) %>%
  ungroup() %>%
  transmute(
    pathway_id   = paste0("TF_", source),
    pathway_name = source,
    database     = "CollecTRI",
    contrast     = condition,
    nes          = score,   # legacy column name; this is the decoupleR ULM score (see caveat above)
    pvalue       = p_value,
    padj         = padj,
    set_size     = set_size,
    core_enrichment = core_enrichment,
    direction    = ifelse(score > 0, "Up", "Down")
  )

dor_set_size <- net_dor %>% filter(target %in% universe) %>%
  count(source, name = "set_size")
master_dor <- dor_wt %>%
  mutate(padj = p.adjust(p_value, "BH")) %>%
  left_join(dor_set_size, by = "source") %>%
  rowwise() %>%
  mutate(core_enrichment = regulon_str(net_dor, source, universe)) %>%
  ungroup() %>%
  transmute(
    pathway_id   = paste0("TF_", source),
    pathway_name = source,
    database     = "DoRothEA_ABC",
    contrast     = "WT_heat",
    nes          = score,   # legacy column name; this is the decoupleR ULM score (see caveat above)
    pvalue       = p_value,
    padj         = padj,
    set_size     = set_size,
    core_enrichment = core_enrichment,
    direction    = ifelse(score > 0, "Up", "Down")
  )

master_tf <- bind_rows(master_ct, master_dor)
dir.create(DIR_MASTER, recursive = TRUE, showWarnings = FALSE)
write_csv(master_tf, file.path(DIR_MASTER, "master_tf_activities.csv"))
cat(sprintf("  master_tf_activities.csv: %d rows (CollecTRI %d + DoRothEA_ABC %d)\n",
            nrow(master_tf), nrow(master_ct), nrow(master_dor)))

# forensics bundle (ALL prior elements preserved + new mechanism blocks)
forensics <- list(
  dorothea_wt_heat = dor_wt,
  collectri_wt_heat = ct_wt,
  mlm_wt_heat = mlm_wt,
  consensus_wt_heat = cons_wt,
  hif1a_rank_summary = tibble(
    method = c("DoRothEA_ULM", "CollecTRI_ULM", "CollecTRI_MLM", "CollecTRI_consensus"),
    rank   = c(hif_dor$rank, hif_ct$rank, hif_mlm$rank, hif_cons$rank),
    n_tfs  = c(n_dor, n_ct, nrow(mlm_wt), nrow(cons_wt)),
    score  = c(hif_dor$score, hif_ct$score, hif_mlm$score, hif_cons$score)
  ),
  target_decomposition = list(summary = decomp, per_gene = present),
  # NEW mechanism blocks
  regulon_swap     = regulon_swap,
  mlm_collinearity = mlm_collinearity,
  consensus_zmix   = list(per_statistic = consensus_zmix,
                          hif_folded_z = hif_zrows,
                          consensus_reproduced = hif_cons_repro,
                          consensus_stored = hif_cons$score),
  top_tf_by_contrast = top_tf_by_contrast,
  # refreshed signature comparison (NEW Lombardi set)
  lombardi_vs_phylo = sig_tab,
  lombardi_genes_present = lombardi_present,
  phylo16_genes_present = phylo16_present,
  n_lombardi = n_lombardi,
  n_phylo = n_phylo,
  egln3_gene_level = egln3_gene,
  key_test_interaction = key_test,
  isg_gene_interaction = isg_gene_int,
  hif1a_heat = hif_heat
)
saveRDS(forensics, file.path(DIR_OBJECTS, "03_tf_forensics.rds"))
cat("  03_tf_forensics.rds written (all prior blocks + regulon_swap, mlm_collinearity, consensus_zmix, top_tf_by_contrast).\n")

# =============================================================================
# 6. ASSERTION: refactor must NOT move carried-over numbers
# =============================================================================

cat("\n[6] Assertion -- carried-over numbers unchanged vs reference ...\n")
rs <- forensics$hif1a_rank_summary
ds <- forensics$target_decomposition$summary
kt <- forensics$key_test_interaction
pct <- function(cl) ds$pct_of_total[ds$class == cl]
padj_of <- function(tf) kt$padj[kt$TF == tf]

stopifnot(
  # rank cascade ranks 1, 12, 142, 8
  rs$rank[rs$method == "DoRothEA_ULM"]        == 1,
  rs$rank[rs$method == "CollecTRI_ULM"]       == 12,
  rs$rank[rs$method == "CollecTRI_MLM"]       == 142,
  rs$rank[rs$method == "CollecTRI_consensus"] == 8,
  # rank cascade scores 6.087, 5.110, 1.134, 2.890 (tol 1e-2)
  abs(rs$score[rs$method == "DoRothEA_ULM"]        - 6.087) < 1e-2,
  abs(rs$score[rs$method == "CollecTRI_ULM"]       - 5.110) < 1e-2,
  abs(rs$score[rs$method == "CollecTRI_MLM"]       - 1.134) < 1e-2,
  abs(rs$score[rs$method == "CollecTRI_consensus"] - 2.890) < 1e-2,
  # target decomposition (reconciled vocab; item 4). Contrib VALUES unchanged:
  #   "other" % is identical to before (91.8%); the curated genes are merely re-bucketed.
  #   "hypoxic HIF core" sum_contrib == fig3l's hif1a_hypoxic_core (-8.46) -> "the HIF core"
  #   now denotes ONE set across fig3g and fig3l.
  abs(pct("other")             - 91.8) < 0.2,
  abs(ds$sum_contrib[ds$class == "hypoxic HIF core"] - (-8.46)) < 0.1,
  abs(ds$sum_contrib[ds$class == "shared/glycolytic"] - 39.80) < 0.1,
  abs(pct("other")             - 91.8) < 0.2,
  # key_test_interaction directionality of significance
  padj_of("Hif1a") >= 0.05,
  padj_of("Epas1") >= 0.05,
  padj_of("Irf3")  < 0.05,
  padj_of("Stat1") < 0.05,
  padj_of("Stat2") < 0.05
)
cat("  ASSERTION PASSED -- all carried-over numbers preserved.\n")

# =============================================================================
# 7. TIDY CSVs -- one self-contained, plot-ready table per downstream figure
# =============================================================================

cat("\n[7] Writing plot-ready tidy tables -> ", TBL_DIR, " ...\n", sep = "")
written <- character(0)
emit <- function(df, fname) {
  p <- file.path(TBL_DIR, fname)
  write.csv(df, p, row.names = FALSE)
  written <<- c(written, fname)
  cat(sprintf("  [SAVE] %-42s %4d rows x %d cols\n", fname, nrow(df), ncol(df)))
}

# fig3a -- TFs on the Interaction contrast (HIF + IFN axis + top up/down)
inter_all <- tf_ct %>% filter(condition == "Interaction")
inter_top <- bind_rows(
  inter_all %>% arrange(desc(score)) %>% head(TOPN_CONTRAST),
  inter_all %>% arrange(score)       %>% head(TOPN_CONTRAST),
  inter_all %>% filter(source %in% c(HIF_AXIS, IFN_AXIS))
) %>% distinct(source, .keep_all = TRUE)
fig3a <- inter_top %>%
  transmute(source, score, padj, p_value,
            axis = ifelse(source %in% HIF_AXIS, "HIF",
                   ifelse(source %in% IFN_AXIS, "IFN", "other")),
            sig = padj < 0.05,
            direction = ifelse(score > 0, "Up", "Down")) %>%
  arrange(desc(score))
emit(fig3a, "fig3a_tf_interaction_axes_data.csv")

# fig3b -- top TF by contrast (block D)
fig3b <- top_tf_by_contrast %>%
  mutate(contrast = factor(contrast,
                           levels = c("WT_heat", "KO_heat", "Temp_main", "Interaction")),
         sig = padj < 0.05) %>%
  arrange(contrast, direction, rank)
emit(fig3b, "fig3b_top_tf_by_contrast_data.csv")

# fig3c -- HIF1a rank cascade (method factor order = the cascade order)
fig3c <- forensics$hif1a_rank_summary %>%
  mutate(method = factor(method,
                         levels = c("DoRothEA_ULM", "CollecTRI_ULM",
                                    "CollecTRI_MLM", "CollecTRI_consensus")),
         method_order = as.integer(method),
         pct_rank = 100 * rank / n_tfs) %>%
  arrange(method_order)
emit(fig3c, "fig3c_hif1a_rank_cascade_data.csv")

# fig3d -- regulon swap (per-gene + summary)
fig3d <- regulon_swap$per_gene %>%
  mutate(abs_t = abs(t_wt), high_t = abs(t_wt) > 2)
emit(fig3d, "fig3d_regulon_swap_data.csv")
emit(regulon_swap$summary, "fig3d_regulon_swap_summary.csv")

# fig3e -- mlm collinearity (per-gene + score collapse)
fig3e <- mlm_collinearity$per_gene %>%
  mutate(co_regulated = n_other_TFs >= 1)
emit(fig3e, "fig3e_mlm_collinearity_data.csv")
emit(mlm_collinearity$summary, "fig3e_mlm_collinearity_summary.csv")
emit(mlm_collinearity$score_collapse, "fig3e_score_collapse.csv")

# fig3f -- consensus z-mix (per statistic for Hif1a, family-flagged)
emit(consensus_zmix, "fig3f_consensus_zmix_data.csv")

# fig3g -- target decomposition (per-gene contrib + per-class summary)
fig3g_pergene <- present %>%
  select(source, target, mor, t_wt, contrib, class) %>%
  arrange(class, desc(contrib))
emit(fig3g_pergene, "fig3g_target_decomposition_data.csv")
fig3g_summary <- decomp %>%
  mutate(class = factor(class,
                        levels = c("hypoxic HIF core", "shared/glycolytic", "other")))
emit(fig3g_summary, "fig3g_target_decomposition_summary.csv")

# fig3h -- lombardi vs phylo (block E; carries n_lombardi/n_phylo + dynamic label)
fig3h <- sig_tab %>%
  mutate(condition = factor(condition,
                            levels = c("WT_heat", "KO_heat", "Interaction")),
         n_lombardi = n_lombardi, n_phylo = n_phylo)
emit(fig3h, "fig3h_lombardi_vs_phylo_data.csv")

# fig3i -- Interaction primer: 4 group means + 3 contrast values for Ifit1 & Vegfa
# Assembles from ALREADY-COMPUTED values only (no new statistics).
# Source A: four-group means from 03_results/03_de/tables/fig2_marker_means.csv
# Source B: WT_heat / KO_heat / Interaction logFC + adj.P.Val from 02_de_results.rds
# Gene choice: Ifit1 (IFN arm, cGAS-dependent; slopes DIVERGE, Interaction adjP 0.003)
#              Vegfa (HIF arm, no detectable cGAS-dependence; slopes PARALLEL, Interaction adjP 1.00)
primer_genes <- c("Ifit1", "Vegfa")
marker_means_path <- file.path(PROJECT_ROOT, "03_results/03_de/tables/fig2_marker_means.csv")
stopifnot(file.exists(marker_means_path))
primer_means <- read.csv(marker_means_path, check.names = FALSE, stringsAsFactors = FALSE)

# Keep only rows for the two primer genes and select the four-group-mean columns we need.
primer_means <- primer_means[primer_means$gene %in% primer_genes,
                             c("gene", "arm", "group", "genotype", "temp", "mean_log2cpm")]

# Pull contrast values per gene from the DE checkpoint (WT_heat, KO_heat, Interaction).
get_contrast_val <- function(gene_sym, contrast, col) {
  row <- de[[contrast]][de[[contrast]]$gene_symbol == gene_sym, ]
  if (nrow(row) == 0) return(NA_real_)
  row[[col]][1]
}

primer_contrasts <- do.call(rbind, lapply(primer_genes, function(g) {
  data.frame(
    gene       = g,
    wt_heat    = get_contrast_val(g, "WT_heat",    "logFC"),
    ko_heat    = get_contrast_val(g, "KO_heat",    "logFC"),
    interaction = get_contrast_val(g, "Interaction", "logFC"),
    inter_adjP = get_contrast_val(g, "Interaction", "adj.P.Val"),
    stringsAsFactors = FALSE
  )
}))

# Assign arm label matching the goal spec + the project framing constraint
# (never "cGAS-independent" -- the n=5 Interaction null is "no detectable cGAS-dependence").
arm_map <- c("Ifit1" = "cGAS-dependent", "Vegfa" = "no-detectable-cGAS-dependence")

fig3i <- merge(primer_means, primer_contrasts, by = "gene") %>%
  mutate(
    arm    = arm_map[gene],
    temp   = sub("C$", "", temp),   # "37C" -> "37" (viz uses ordered factor)
    gene   = factor(gene, levels = primer_genes)
  ) %>%
  arrange(gene, genotype, temp) %>%
  select(gene, arm, genotype, temp, mean_log2cpm, wt_heat, ko_heat, interaction, inter_adjP)

# Assertion: both genes must be present with expected diverge/parallel pattern.
stopifnot(all(c("Ifit1", "Vegfa") %in% fig3i$gene))
ifit1_row <- fig3i[fig3i$gene == "Ifit1", ][1, ]
vegfa_row  <- fig3i[fig3i$gene == "Vegfa",  ][1, ]
stopifnot(abs(ifit1_row$interaction) > 0.5)   # Ifit1 must diverge (|Interaction| > 0.5)
stopifnot(abs(vegfa_row$interaction)  < 0.1)   # Vegfa must be parallel (|Interaction| < 0.1)

emit(fig3i, "fig3i_interaction_primer_data.csv")

cat(sprintf("\n  Tidy-table inventory (%d files):\n    %s\n",
            length(written), paste(written, collapse = "\n    ")))

cat("\n[DONE] Phase 3 COMPUTE complete. Run 03_decoupler_tf_viz.R for figures.\n")
