# 16_synthesis.R — COMPUTE
# =============================================================================
# CAPSTONE synthesis (stage 07_synthesis) — assemble the unified, cross-arm
# TWO-ARMS evidence table for the STING-cGAS standard sweep.
#
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   07_synthesis
#
# ROLE: COMPUTE ONLY. This script runs NO new statistics — no DE re-fit, no
#   fgsea, no decoupleR, no permutation. It only JOINS / FILTERS / SUMMARISES the
#   masters the earlier sweep arms already produced, into one tidy cross-arm
#   table. There is NO ggplot/ggsave here; the headline panel lives in the
#   sibling 16_synthesis_viz.R, which reads only what this script writes.
#
# THE QUESTION THIS TABLE ANSWERS (the publication-relevant payoff):
#   Across every method in the sweep (GSEA Hallmark/Reactome, the custom DBs,
#   PROGENy, decoupleR-TF, GATOM/CoReSh where present), is the heat response
#   cGAS-DEPENDENT? The genotype x temperature INTERACTION contrast is the
#   cGAS-dependence test (1 df, the lowest-powered comparison). The expectation,
#   to be read out gene-set-wide and multi-method:
#     * IFN/ISG arm (cGAS-DEPENDENT hypothesis): Interaction significant + positive
#       (present in WT_heat, lost/reversed in KO_heat).
#     * HIF/glycolysis arm (NO-DETECTABLE-cGAS-DEPENDENCE hypothesis): rises in BOTH
#       WT_heat AND KO_heat, Interaction near-zero / non-significant.
#   We frame this as an ASYMMETRY, never as "cGAS-independent" (n=5; underpowered).
#
# Inputs (read-only MASTERS; each guarded — warn+skip if absent):
#   03_results/master/master_gsea_table.csv          (GSEA: MSigDB collections +
#                                                     custom DBs + CoReSh_derived +
#                                                     PROGENy + GATOM accumulate here)
#   03_results/master/master_progeny_activities.csv  (PROGENy MLM activities)
#   03_results/master/master_tf_activities.csv        (decoupleR TF activities; CollecTRI)
#   03_results/master/master_de_genes.csv             (limma-trend per-gene DE — marker rollup)
#   03_results/master/master_gatom_modules.csv        (GATOM metabolic-module pseudo-sets)
#
# Outputs:
#   03_results/07_synthesis/tables/two_arms_summary.csv  (the tidy cross-arm table)
#   03_results/objects/16_synthesis.rds                  (the assembled object)
#
# IDEMPOTENT + BYTE-STABLE: a pure read->join->round->write of fixed masters.
#   Re-running yields the same bytes (round_numeric_cols, 9 sig figs).
#
# Dependencies: dplyr/tibble/readr (base project stack via de_gsea_helpers.R).
#   No heavy enrichment packages are attached (this is a join-only stage).
# =============================================================================

source("02_analysis/config/config.R")            # PROJECT_ROOT, YAML_CONFIG, DIR_*, %||%, RANK_METRIC
source("02_analysis/helpers/de_gsea_helpers.R")    # round_numeric_cols, append_master_table (compute idiom)

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
})
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================================
# CONSTANTS (from config; never hardcoded — AGENTS.md rule 1)
# ============================================================================

STAGE  <- "07_synthesis"
SCRIPT <- "02_analysis/scripts/16_synthesis.R"
FDR    <- as.numeric(YAML_CONFIG$thresholds$gsea_fdr %||% 0.05)

# Headline contrasts (config order). WT_heat | KO_heat | Interaction | Temp_main.
# Interaction is the cGAS-dependence test (difference-of-differences, 1 df).
HEADLINE_CONTRASTS <- c("WT_heat", "KO_heat", "Interaction", "Temp_main")

# The two arms. Each arm = a hypothesis label + the corroborating gene-set /
# pathway / TF / gene names to pull from each master. The arm definitions are
# the single source of truth that maps method-specific rows onto the two-arms axis.
#
# IFN/ISG arm  -> the cGAS-DEPENDENT hypothesis (positive significant Interaction).
# HIF/glycolysis arm -> the NO-DETECTABLE-cGAS-DEPENDENCE hypothesis (flat Interaction;
#   rises in both heat arms). NEVER "cGAS-independent".
ARM_IFN <- "IFN_ISG"          # cGAS-dependent
ARM_HIF <- "HIF_glycolysis"   # no detectable cGAS-dependence at n=5

ARM_LABELS <- c(
  IFN_ISG        = "IFN/ISG (cGAS-dependent hypothesis)",
  HIF_glycolysis = "HIF/glycolysis (no-detectable-cGAS-dependence hypothesis)"
)

# ---- Arm membership patterns, per method --------------------------------------
# These are case-insensitive regex on the method's NAME column (gene_set / pathway /
# TF / gene). They are deliberately specific so an arm only collects rows that
# genuinely belong to it; everything else is left out of the two-arms table.

# GSEA gene-set name patterns (matched against pathway_name OR pathway_id).
GSEA_IFN_PAT <- paste0(
  "INTERFERON|_IFN|TYPE_I_IFN|ISG|JAK_STAT|JAK-STAT|",
  "CGAS|STING|DNA_SENSING|RIG_I|ANTIVIRAL|RESPONSE_TO_(TYPE|INTERFERON)|TNFA_SIGNALING_VIA_NFKB|NFKB")
GSEA_HIF_PAT <- paste0(
  "HYPOXIA|HIF|GLYCOLY|GLUCOSE_(METABOL|TRANSPORT)|",
  "PYRUVATE|LACTATE|LOMBARDI|PSEUDOHYPOXIA|VEGF|ANGIOGENESIS")

# PROGENy pathway names (the footprint method; no HIF regulon).
PROGENY_IFN <- c("JAK-STAT", "NFkB", "TNFa")
PROGENY_HIF <- c("Hypoxia")

# decoupleR TF axes (the CollecTRI ULM forensics axis families).
TF_IFN <- c("Irf1", "Irf3", "Irf7", "Stat1", "Stat2", "Nfkb1", "Rela")
TF_HIF <- c("Hif1a", "Epas1")
# Hsf1 is the true thermal sensor (heat-shock arm); it is the contamination that
# the bespoke 16-gene HIF score rode. Tracked separately as context, NOT crowned
# and NOT folded into the HIF arm (it is heat-shock, not hypoxia).
TF_HEATSHOCK <- c("Hsf1")

# Marker genes for the per-gene DE rollup (from config.R; the single source of truth).
DE_IFN_MARKERS <- ISG_MARKERS                                  # cGAS-dependent ISG core
DE_HIF_MARKERS <- HIF_GLYCO_MARKERS                            # HIF-specific + shared-glycolytic

# GATOM is the metabolic-network corroboration; its modules read as the
# HIF/glycolysis (metabolic) arm by construction (Complex-I / glycolytic shift).
GATOM_ARM <- ARM_HIF

# Master-table paths (bare filenames resolved under DIR_MASTER).
MASTER_PATHS <- list(
  gsea    = file.path(DIR_MASTER, "master_gsea_table.csv"),
  progeny = file.path(DIR_MASTER, "master_progeny_activities.csv"),
  tf      = file.path(DIR_MASTER, "master_tf_activities.csv"),
  de      = file.path(DIR_MASTER, "master_de_genes.csv"),
  gatom   = file.path(DIR_MASTER, "master_gatom_modules.csv")
)

# ============================================================================
# DIRECTORIES
# ============================================================================

tbl_dir <- stage_dir(STAGE, "tables")          # 03_results/07_synthesis/tables/; creates dir
dir.create(DIR_OBJECTS, recursive = TRUE, showWarnings = FALSE)

message("=================================================================")
message("16_synthesis: assemble two-arms cross-arm evidence table (07_synthesis)")
message("  COMPUTE ONLY — joins/summarises existing masters; no new statistics.")
message("=================================================================")

# ============================================================================
# 0. GUARDED MASTER READER — warn+skip if a master is absent (e.g. CoReSh /
#    GSEA / PROGENy / GATOM not yet provisioned in-container). Returns NULL on
#    a miss so each arm builder can degrade gracefully. Never fabricates rows.
# ============================================================================

read_master <- function(key) {
  path <- MASTER_PATHS[[key]]
  if (is.null(path) || !file.exists(path)) {
    warning(sprintf("[16_synthesis] master '%s' ABSENT (%s) — SKIPPING this arm source.",
                    key, path %||% "<unset>"), call. = FALSE)
    return(NULL)
  }
  df <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) {
      warning(sprintf("[16_synthesis] failed to read master '%s' (%s): %s — SKIPPING.",
                      key, path, conditionMessage(e)), call. = FALSE)
      NULL
    })
  if (!is.null(df))
    message(sprintf("  [read] %-8s %s  (%d rows)", key, basename(path), nrow(df)))
  df
}

masters <- lapply(setNames(names(MASTER_PATHS), names(MASTER_PATHS)), read_master)
absent  <- names(Filter(is.null, masters))
present <- names(Filter(Negate(is.null), masters))
if (length(absent))
  message("  [guard] absent masters (skipped): ", paste(absent, collapse = ", "))
message("  [guard] present masters: ", paste(present, collapse = ", "))

# ============================================================================
# 1. THE UNIFIED ROW SHAPE
#    Every arm builder emits rows in this single tidy schema so methods stack
#    into one table. The load-bearing columns the plan names:
#      method, gene_set/pathway/TF (-> `feature`), contrast, nes/score (-> `score`),
#      padj, direction, significant. Plus arm + arm_label + hypothesis context.
# ============================================================================

EVIDENCE_COLS <- c(
  "arm", "arm_label", "arm_hypothesis",
  "method", "feature", "feature_kind",
  "contrast", "score", "pvalue", "padj", "direction", "significant")

# Empty, well-typed evidence frame (so binds stay clean when an arm is skipped).
empty_evidence <- function() {
  data.frame(
    arm = character(), arm_label = character(), arm_hypothesis = character(),
    method = character(), feature = character(), feature_kind = character(),
    contrast = character(), score = numeric(), pvalue = numeric(),
    padj = numeric(), direction = character(), significant = logical(),
    stringsAsFactors = FALSE)
}

# Normalise a method's rows into the evidence schema. `arm` is assigned by the
# caller (which already filtered to this arm's features). `significant` is
# recomputed uniformly here from padj < FDR so the flag means the same thing
# across methods (single source of truth; not trusting per-arm `significant`).
to_evidence <- function(df, arm, method, feature_kind,
                        feature_col, score_col, padj_col, pval_col = NULL,
                        direction_col = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(empty_evidence())
  feat  <- as.character(df[[feature_col]])
  score <- suppressWarnings(as.numeric(df[[score_col]]))
  padj  <- suppressWarnings(as.numeric(df[[padj_col]]))
  pval  <- if (!is.null(pval_col) && pval_col %in% names(df))
             suppressWarnings(as.numeric(df[[pval_col]])) else NA_real_
  dirn  <- if (!is.null(direction_col) && direction_col %in% names(df))
             as.character(df[[direction_col]])
           else ifelse(is.na(score), NA_character_, ifelse(score > 0, "Up", "Down"))
  data.frame(
    arm            = arm,
    arm_label      = unname(ARM_LABELS[arm]),
    arm_hypothesis = if (identical(arm, ARM_IFN)) "cGAS-dependent"
                     else "no detectable cGAS-dependence at n=5",
    method         = method,
    feature        = feat,
    feature_kind   = feature_kind,
    contrast       = as.character(df[["contrast"]]),
    score          = score,
    pvalue         = pval,
    padj           = padj,
    direction      = dirn,
    significant    = !is.na(padj) & padj < FDR,
    stringsAsFactors = FALSE)
}

# ============================================================================
# 2. ARM BUILDERS — one per method. Each filters its master to the headline
#    contrasts and to this arm's features, then maps onto the evidence schema.
# ============================================================================

# ---- 2a. GSEA (MSigDB Hallmark/Reactome/... + custom DBs + CoReSh_derived) -----
# The GSEA master accumulates ALL gene-set arms under one schema (PROGENy and
# GATOM also append here, but we read those from their dedicated masters below,
# so exclude database %in% {PROGENy, GATOM} here to avoid double-counting).
build_gsea_arm <- function() {
  df <- masters$gsea
  if (is.null(df)) return(empty_evidence())
  if (!all(c("pathway_name", "database", "nes", "padj", "contrast") %in% names(df))) {
    warning("[16_synthesis] master_gsea_table missing expected columns — skipping GSEA arm.",
            call. = FALSE)
    return(empty_evidence())
  }
  df <- df %>%
    dplyr::filter(contrast %in% HEADLINE_CONTRASTS,
                  !(.data$database %in% c("PROGENy", "GATOM")))
  # Match arm patterns on BOTH the pathway_name and pathway_id (ids carry the
  # collection prefix, names the cleaned form — catch either spelling).
  name_blob <- paste(df$pathway_name %||% "", df$pathway_id %||% "")
  is_ifn <- grepl(GSEA_IFN_PAT, name_blob, ignore.case = TRUE)
  is_hif <- grepl(GSEA_HIF_PAT, name_blob, ignore.case = TRUE)
  # A set could match neither (left out) or, rarely, both — if both, the more
  # specific IFN signalling label wins only when it is an explicit IFN/ISG term;
  # otherwise prefer HIF (hypoxia/glycolysis names are the more specific here).
  ifn_df <- df[is_ifn & !is_hif, , drop = FALSE]
  hif_df <- df[is_hif & !is_ifn, , drop = FALSE]
  # Ambiguous (both) -> resolve by which pattern the NAME (not id) hits first;
  # default to HIF for hypoxia-glycolysis dominance, keep IFN only for clear IFN.
  both_df <- df[is_ifn & is_hif, , drop = FALSE]
  if (nrow(both_df)) {
    clearly_ifn <- grepl("INTERFERON|ISG|JAK_STAT|JAK-STAT|CGAS|STING|ANTIVIRAL",
                         paste(both_df$pathway_name, both_df$pathway_id), ignore.case = TRUE)
    ifn_df <- dplyr::bind_rows(ifn_df, both_df[clearly_ifn, , drop = FALSE])
    hif_df <- dplyr::bind_rows(hif_df, both_df[!clearly_ifn, , drop = FALSE])
  }
  dplyr::bind_rows(
    to_evidence(ifn_df, ARM_IFN, method = paste0("GSEA:", ifn_df$database %||% character()),
                feature_kind = "gene_set", feature_col = "pathway_name",
                score_col = "nes", padj_col = "padj", pval_col = "pvalue",
                direction_col = "direction"),
    to_evidence(hif_df, ARM_HIF, method = paste0("GSEA:", hif_df$database %||% character()),
                feature_kind = "gene_set", feature_col = "pathway_name",
                score_col = "nes", padj_col = "padj", pval_col = "pvalue",
                direction_col = "direction"))
}

# ---- 2b. PROGENy (footprint method; orthogonal — uses NO HIF regulon) ----------
build_progeny_arm <- function() {
  df <- masters$progeny
  if (is.null(df)) return(empty_evidence())
  if (!all(c("pathway_name", "nes", "padj", "contrast") %in% names(df))) {
    warning("[16_synthesis] master_progeny_activities missing expected columns — skipping PROGENy arm.",
            call. = FALSE)
    return(empty_evidence())
  }
  df <- dplyr::filter(df, contrast %in% HEADLINE_CONTRASTS)
  ifn_df <- dplyr::filter(df, pathway_name %in% PROGENY_IFN)
  hif_df <- dplyr::filter(df, pathway_name %in% PROGENY_HIF)
  dplyr::bind_rows(
    to_evidence(ifn_df, ARM_IFN, method = "PROGENy", feature_kind = "pathway",
                feature_col = "pathway_name", score_col = "nes",
                padj_col = "padj", pval_col = "pvalue", direction_col = "direction"),
    to_evidence(hif_df, ARM_HIF, method = "PROGENy", feature_kind = "pathway",
                feature_col = "pathway_name", score_col = "nes",
                padj_col = "padj", pval_col = "pvalue", direction_col = "direction"))
}

# ---- 2c. decoupleR TF activity (CollecTRI ULM) ---------------------------------
build_tf_arm <- function() {
  df <- masters$tf
  if (is.null(df)) return(empty_evidence())
  if (!all(c("pathway_name", "database", "nes", "padj", "contrast") %in% names(df))) {
    warning("[16_synthesis] master_tf_activities missing expected columns — skipping TF arm.",
            call. = FALSE)
    return(empty_evidence())
  }
  # CollecTRI is the de-confounded primary network (run_consensus/MLM in the
  # forensics); use it as the TF read for the two-arms table.
  df <- df %>%
    dplyr::filter(contrast %in% HEADLINE_CONTRASTS) %>%
    dplyr::filter(.data$database == "CollecTRI" |
                  !"CollecTRI" %in% unique(.data$database))   # fall back if CollecTRI absent
  ifn_df <- dplyr::filter(df, pathway_name %in% TF_IFN)
  hif_df <- dplyr::filter(df, pathway_name %in% TF_HIF)
  dplyr::bind_rows(
    to_evidence(ifn_df, ARM_IFN, method = "TF:CollecTRI", feature_kind = "TF",
                feature_col = "pathway_name", score_col = "nes",
                padj_col = "padj", pval_col = "pvalue", direction_col = "direction"),
    to_evidence(hif_df, ARM_HIF, method = "TF:CollecTRI", feature_kind = "TF",
                feature_col = "pathway_name", score_col = "nes",
                padj_col = "padj", pval_col = "pvalue", direction_col = "direction"))
}

# ---- 2d. GATOM metabolic modules (network corroboration; HIF/metabolic arm) -----
build_gatom_arm <- function() {
  df <- masters$gatom
  if (is.null(df)) return(empty_evidence())
  if (!all(c("pathway_name", "nes", "padj", "contrast") %in% names(df))) {
    warning("[16_synthesis] master_gatom_modules missing expected columns — skipping GATOM arm.",
            call. = FALSE)
    return(empty_evidence())
  }
  df <- dplyr::filter(df, contrast %in% HEADLINE_CONTRASTS)
  # GATOM pseudo-NES = mean module-edge log2FC; padj = geo-mean edge pval. Maps
  # onto the HIF/glycolysis (metabolic) arm by construction.
  to_evidence(df, GATOM_ARM, method = "GATOM", feature_kind = "metabolic_module",
              feature_col = "pathway_name", score_col = "nes",
              padj_col = "padj", pval_col = "pvalue", direction_col = "direction")
}

# ---- 2e. per-gene DE marker rollup (limma-trend t / FDR) ------------------------
# Not a gene-set method, but the concrete per-gene anchor of each arm: the ISG
# core (cGAS-dependent) vs the HIF/glycolysis markers. We emit the marker GENES
# directly so the table carries the per-gene evidence behind the gene-set claims.
build_de_arm <- function() {
  df <- masters$de
  if (is.null(df)) return(empty_evidence())
  need <- c("gene_symbol", "t", "adj.P.Val", "P.Value", "contrast", "direction")
  if (!all(need %in% names(df))) {
    warning("[16_synthesis] master_de_genes missing expected columns — skipping DE arm.",
            call. = FALSE)
    return(empty_evidence())
  }
  df <- dplyr::filter(df, contrast %in% HEADLINE_CONTRASTS)
  ifn_df <- dplyr::filter(df, gene_symbol %in% DE_IFN_MARKERS)
  hif_df <- dplyr::filter(df, gene_symbol %in% DE_HIF_MARKERS)
  # DE score = the limma t-statistic (the rank metric); padj = adj.P.Val.
  dplyr::bind_rows(
    to_evidence(ifn_df, ARM_IFN, method = "DE:limma-trend", feature_kind = "gene",
                feature_col = "gene_symbol", score_col = "t",
                padj_col = "adj.P.Val", pval_col = "P.Value", direction_col = "direction"),
    to_evidence(hif_df, ARM_HIF, method = "DE:limma-trend", feature_kind = "gene",
                feature_col = "gene_symbol", score_col = "t",
                padj_col = "adj.P.Val", pval_col = "P.Value", direction_col = "direction"))
}

# Heat-shock context row(s) — Hsf1 from the TF master, kept as a separate
# feature_kind so the reply can say "the thermal sensor is heat-shock, not
# hypoxia" without folding HSF1 into either arm's significance tally.
build_heatshock_context <- function() {
  df <- masters$tf
  if (is.null(df)) return(empty_evidence())
  if (!all(c("pathway_name", "nes", "padj", "contrast") %in% names(df)))
    return(empty_evidence())
  hs <- df %>%
    dplyr::filter(contrast %in% HEADLINE_CONTRASTS, pathway_name %in% TF_HEATSHOCK)
  if (!"database" %in% names(hs)) return(empty_evidence())
  hs <- dplyr::filter(hs, .data$database == "CollecTRI" |
                          !"CollecTRI" %in% unique(.data$database))
  ev <- to_evidence(hs, ARM_HIF, method = "TF:CollecTRI", feature_kind = "heatshock_context",
                    feature_col = "pathway_name", score_col = "nes",
                    padj_col = "padj", pval_col = "pvalue", direction_col = "direction")
  # Override the arm framing: HSF1 is CONTEXT, neither arm's hypothesis.
  if (nrow(ev)) {
    ev$arm            <- "heatshock_context"
    ev$arm_label      <- "Heat-shock context (HSF1; thermal sensor, not hypoxia)"
    ev$arm_hypothesis <- "context (not a cGAS-dependence arm)"
  }
  ev
}

# ============================================================================
# 3. ASSEMBLE THE LONG EVIDENCE TABLE
# ============================================================================

message("[1] Building per-method arm evidence ...")
evidence <- dplyr::bind_rows(
  build_gsea_arm(),
  build_progeny_arm(),
  build_tf_arm(),
  build_gatom_arm(),
  build_de_arm(),
  build_heatshock_context()
)

# Stable factor order for contrasts + arms (so the table + the figure read in
# the headline order regardless of master row order).
evidence <- evidence %>%
  dplyr::filter(!is.na(contrast)) %>%
  dplyr::mutate(
    contrast = factor(contrast, levels = HEADLINE_CONTRASTS),
    arm      = factor(arm, levels = c(ARM_IFN, ARM_HIF, "heatshock_context"))
  ) %>%
  dplyr::arrange(arm, method, feature, contrast) %>%
  dplyr::mutate(contrast = as.character(contrast),
                arm      = as.character(arm))

n_methods <- length(unique(evidence$method))
n_feat    <- length(unique(evidence$feature))
message(sprintf("  evidence rows: %d  (%d methods, %d features, %d contrasts)",
                nrow(evidence), n_methods, n_feat,
                length(unique(evidence$contrast))))

if (nrow(evidence) == 0L)
  warning("[16_synthesis] NO evidence rows assembled — all arm masters absent/empty. ",
          "Run the compute arms (05/06/08/09/10) in-container, then re-run this script.")

# ============================================================================
# 4. THE INTERACTION READOUT — the cGAS-dependence test, per arm.
#    This is the single most publication-relevant summary: does the Interaction
#    contrast separate the two arms? (IFN expected significant+positive; HIF
#    expected NS / near-zero.) We summarise per (arm, method) on the Interaction
#    contrast specifically, plus the WT_heat/KO_heat heat-program read.
# ============================================================================

message("[2] Summarising the Interaction (cGAS-dependence) readout per arm ...")

summarise_contrast <- function(ev, co) {
  ev %>%
    dplyr::filter(contrast == co, arm %in% c(ARM_IFN, ARM_HIF)) %>%
    dplyr::group_by(arm, arm_label, arm_hypothesis) %>%
    dplyr::summarise(
      contrast      = co,
      n_features    = dplyr::n(),
      n_significant = sum(significant, na.rm = TRUE),
      frac_sig      = ifelse(dplyr::n() > 0, mean(significant, na.rm = TRUE), NA_real_),
      mean_score    = mean(score, na.rm = TRUE),
      median_score  = stats::median(score, na.rm = TRUE),
      n_up          = sum(direction == "Up", na.rm = TRUE),
      n_down        = sum(direction == "Down", na.rm = TRUE),
      .groups = "drop")
}

arm_contrast_summary <- dplyr::bind_rows(lapply(HEADLINE_CONTRASTS, function(co)
  summarise_contrast(evidence, co)))
if (nrow(arm_contrast_summary))
  arm_contrast_summary <- arm_contrast_summary %>%
    dplyr::arrange(factor(arm, levels = c(ARM_IFN, ARM_HIF)),
                   factor(contrast, levels = HEADLINE_CONTRASTS))

# Console read of the asymmetry on the Interaction contrast.
int_summary <- dplyr::filter(arm_contrast_summary, contrast == "Interaction")
if (nrow(int_summary)) {
  message("  --- Interaction contrast (the cGAS-dependence test) ---")
  for (i in seq_len(nrow(int_summary))) {
    r <- int_summary[i, ]
    message(sprintf("    %-16s  n=%2d  n_sig=%2d (%.0f%%)  mean_score=%+.2f  [%s]",
                    r$arm, r$n_features, r$n_significant, 100 * (r$frac_sig %||% 0),
                    r$mean_score %||% NA_real_, r$arm_hypothesis))
  }
}

# ============================================================================
# 5. ASSEMBLE THE FINAL OBJECT + WRITE
# ============================================================================

message("[3] Writing two_arms_summary.csv + 16_synthesis.rds ...")

# The flat, plot-ready cross-arm table (the deliverable the viz reads).
two_arms <- evidence %>%
  dplyr::select(dplyr::all_of(EVIDENCE_COLS))

# Byte-stable write (9 sig figs on doubles; characters untouched).
out_csv <- file.path(tbl_dir, "two_arms_summary.csv")
readr::write_csv(round_numeric_cols(two_arms), out_csv)
message(sprintf("  two_arms_summary.csv: %d rows x %d cols", nrow(two_arms), ncol(two_arms)))

# The assembled object: the long table + the per-arm Interaction summary +
# provenance (which masters were present/absent, the arm definitions, the FDR).
synthesis <- list(
  two_arms              = as.data.frame(round_numeric_cols(two_arms)),
  arm_contrast_summary  = as.data.frame(round_numeric_cols(arm_contrast_summary)),
  arm_labels            = ARM_LABELS,
  headline_contrasts    = HEADLINE_CONTRASTS,
  interaction_contrast  = "Interaction",
  fdr                   = FDR,
  arm_definitions = list(
    gsea_ifn_pattern = GSEA_IFN_PAT, gsea_hif_pattern = GSEA_HIF_PAT,
    progeny_ifn = PROGENY_IFN, progeny_hif = PROGENY_HIF,
    tf_ifn = TF_IFN, tf_hif = TF_HIF, tf_heatshock = TF_HEATSHOCK,
    de_ifn_markers = DE_IFN_MARKERS, de_hif_markers = DE_HIF_MARKERS,
    gatom_arm = GATOM_ARM),
  masters_present       = present,
  masters_absent        = absent,
  provisional_caveat    = provisional_caption(),
  built_at              = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")
)
out_rds <- file.path(DIR_OBJECTS, "16_synthesis.rds")
saveRDS(synthesis, out_rds)
message("  16_synthesis.rds written (assembled object).")

# ============================================================================
# 6. ACCEPTANCE CHECKS (structural; no statistics)
# ============================================================================

message("[4] Acceptance checks ...")
stopifnot(
  "two_arms_summary.csv missing" = file.exists(out_csv),
  "16_synthesis.rds missing"     = file.exists(out_rds),
  "evidence schema drift"        = all(EVIDENCE_COLS %in% colnames(two_arms))
)
# If any gene-set arm IS present, the table must be non-empty (catch a silent join bug).
if (length(intersect(present, c("gsea", "progeny", "tf", "de", "gatom"))) > 0 &&
    nrow(two_arms) == 0L)
  warning("[16_synthesis] masters present but evidence table is empty — check arm filters.")

message("=================================================================")
message("16_synthesis COMPLETE")
message("  Table:  ", out_csv)
message("  Object: ", out_rds)
message("  Rows:   ", nrow(two_arms), "  (present masters: ",
        paste(present, collapse = ", "),
        if (length(absent)) paste0("; absent: ", paste(absent, collapse = ", ")) else "",
        ")")
message("  Run 16_synthesis_viz.R for the headline two-arms panel.")
message("=================================================================")
