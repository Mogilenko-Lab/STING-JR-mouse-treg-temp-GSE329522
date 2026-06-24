# 09_activity_progeny.R — COMPUTE
# =============================================================================
# PROGENy pathway activity (14 pathways) via decoupleR MLM on the per-contrast
# limma-trend t-statistic matrix (7 contrasts).
#
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   05_progeny
#
# Role: COMPUTE ONLY. No ggplot/ggsave; figures live in 09_activity_progeny_viz.R.
#       Mirrors the PROGENy-MLM half of 14839/09_activity_progeny_tf.R and
#       matches the load_or_compute + master-table + tidy-CSV idiom of STING's
#       own 03_decoupler_tf.R so PROGENy and TF outputs are consistent.
#
# Key scientific payoff:
#   PROGENy-Hypoxia expected UP in WT_heat AND KO_heat (flat Interaction)
#     => the pseudohypoxia/glycolysis program is cGAS-independent.
#   PROGENy-JAK-STAT expected UP in WT_heat, diminished in KO_heat, positive
#     Interaction => the IFN program is cGAS-dependent.
#   This corroborates the two-arms result by a footprint method that uses NO
#   HIF regulon — the strongest cross-method convergence point.
#
# Inputs:
#   - 03_results/objects/02_de_results.rds       (7 limma topTables; col `t` used)
#   - 03_results/objects/net_progeny_mouse.rds   (PROGENy net; source/target/weight;
#                                                 14 pathways, ~500 genes each;
#                                                 built by 00c_prepare_networks.R via
#                                                 progeny::getModel("Mouse", top=500)
#                                                 then pivot-to-edge-list with filter(weight!=0))
#
# Outputs:
#   - 03_results/objects/09_progeny_activity.rds            (raw decoupleR result list)
#   - 03_results/05_progeny/tables/progeny_activity.csv     (tidy: pathway x contrast)
#   - 03_results/master/master_progeny_activities.csv       (master accumulator; schema-pinned)
#
# Network shape assumed (verified by 00c_prepare_networks.R + stopifnot below):
#   net_progeny: data.frame with columns source (pathway name), target (gene symbol,
#   mouse), weight (continuous float, the PROGENy model coefficient). 14 distinct
#   source values. Genes are already mouse symbols (Title-case MGI). The `weight`
#   column is passed to run_mlm(.mor="weight"), NOT "mor" — PROGENy uses a
#   continuous weight; CollecTRI uses signed ±1 mor. Do NOT swap.
#
# Dependencies: decoupleR (lazy; loaded only when actually needed)
# =============================================================================

source("02_analysis/config/config.R")
source("02_analysis/helpers/de_gsea_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(readr)
  library(purrr)
})

# ============================================================================
# CONSTANTS (from config; never hardcoded — AGENTS.md rule 1)
# ============================================================================

STAGE   <- "05_progeny"
MINSIZE <- as.integer(YAML_CONFIG$thresholds$gsea_min_size %||% 5L)
PADJM   <- YAML_CONFIG$thresholds$padj_method %||% "BH"
RANK    <- RANK_METRIC                                  # "t" from config.R

# PROGENy n=14 pathways; no MT correction warranted at n=14 — padj == pvalue.
# This matches the 14839 template convention: with 14 pathways the expected
# number of false discoveries at α=0.05 is 0.7 — BH still applied for
# transparency but the note is kept for callers reading the CSV.
FDR     <- as.numeric(YAML_CONFIG$thresholds$gsea_fdr %||% 0.05)

# Key PROGENy pathways for the cGAS-dependence asymmetry sanity checks.
KEY_PROGENY <- KEY_PROGENY %||% c("Hypoxia", "JAK-STAT", "NFkB", "TNFa")

# All 7 contrasts declared in config; PROGENy runs on every one.
CONTRASTS <- vapply(YAML_CONFIG$design$contrasts, function(x) x$name, character(1))

# ============================================================================
# DIRECTORIES
# ============================================================================

tbl_dir <- stage_dir(STAGE, "tables")    # 03_results/05_progeny/tables/; also creates dir
dir.create(DIR_MASTER, recursive = TRUE, showWarnings = FALSE)

message("=================================================================")
message("09_activity_progeny: PROGENy MLM pathway activity (stage 05_progeny)")
message("=================================================================")

# ============================================================================
# 1. LOAD DE RESULTS + BUILD GENE x CONTRAST t-STAT MATRIX
#    Mirrors 03_decoupler_tf.R §1 exactly: mat is genes x 7 contrasts; NA -> 0;
#    rownames are gene symbols (the Symbol-rownames contract from load_de_results).
# ============================================================================

message("[1] Loading DE results -> t-stat matrix ...")
de <- load_de_results()                           # from de_gsea_helpers.R; stops loudly on Ensembl rownames
stopifnot("Some declared contrasts missing from 02_de_results.rds — re-run 02_de_limma_trend.R" =
            all(CONTRASTS %in% names(de)))
de <- de[CONTRASTS]                               # enforce config order

genes <- rownames(de[[1]])
stopifnot("Duplicate rownames in topTable — collapse decided in Phase 1 must have run" =
            !anyDuplicated(genes))

# Build gene x contrast t-statistic matrix (NA -> 0, matching 03_decoupler_tf.R).
# decoupleR::run_mlm requires a numeric matrix: genes in rows, samples/contrasts in cols.
mat <- vapply(de, function(d) {
  v <- d[[RANK]]
  v[is.na(v)] <- 0
  v
}, numeric(length(genes)))
rownames(mat) <- genes
stopifnot(is.matrix(mat), nrow(mat) > 0, ncol(mat) == length(CONTRASTS))
message(sprintf("  t-stat matrix: %d genes x %d contrasts", nrow(mat), ncol(mat)))

universe <- rownames(mat)    # gene universe for footprint annotation below

# ============================================================================
# 2. LOAD PROGENY NETWORK
#    net_progeny_mouse.rds: source (pathway) / target (mouse symbol) / weight
#    (continuous float, PROGENy model coefficient).
#    Built by 00c_prepare_networks.R via:
#      prog_mat <- progeny::getModel("Mouse", top=500)   # gene x pathway matrix
#      pivot_longer(-target, ...) %>% filter(weight!=0)
#    => columns: source, target, weight. MUST use .mor="weight" in run_mlm.
# ============================================================================

message("[2] Loading PROGENy network ...")
net_progeny_path <- file.path(DIR_OBJECTS, "net_progeny_mouse.rds")
if (!file.exists(net_progeny_path))
  stop("net_progeny_mouse.rds not found: ", net_progeny_path,
       " — run 00c_prepare_networks.R first.")
net_progeny <- readRDS(net_progeny_path)

# Structural asserts — these must hold before we touch run_mlm
stopifnot("source" %in% colnames(net_progeny),
          "target" %in% colnames(net_progeny),
          "weight" %in% colnames(net_progeny))
n_pathways <- length(unique(net_progeny$source))
stopifnot("PROGENy net must have exactly 14 pathways" = n_pathways == 14L)

n_overlap <- sum(unique(net_progeny$target) %in% universe)
if (n_overlap < 100)
  warning("Only ", n_overlap, " PROGENy target genes overlap the DE universe — ",
          "check mouse-symbol casing (expected >2000).")
message(sprintf("  PROGENy net: %d pathways, %d edges, %d targets, %d overlap DE universe",
                n_pathways, nrow(net_progeny), length(unique(net_progeny$target)), n_overlap))

# ============================================================================
# 3. RUN_MLM (PROGENy) — cached; one run across all contrasts
#    .mor="weight" (continuous PROGENy coefficients — NOT "mor").
#    mat is genes x contrasts; decoupleR treats each column as a condition.
#    The result is a long data.frame with columns:
#      statistic, source, condition, score, p_value
#    (condition == contrast name; source == pathway name).
# ============================================================================

message("[3] Running decoupleR::run_mlm (PROGENy, all ", ncol(mat), " contrasts) ...")
if (!requireNamespace("decoupleR", quietly = TRUE))
  stop("Package 'decoupleR' is required. Install with BiocManager::install('decoupleR').")

progeny_raw <- load_or_compute(
  "09_progeny_activity.rds",
  compute_fn = function() {
    res <- decoupleR::run_mlm(
      mat     = mat,
      net     = net_progeny,
      .source = "source",
      .target = "target",
      .mor    = "weight",    # PROGENy uses a continuous weight — never "mor"
      minsize = MINSIZE
    )
    # Sanity: result must cover all declared pathways and contrasts
    message(sprintf("  run_mlm returned: %d rows, %d unique pathways, %d unique contrasts",
                    nrow(res),
                    length(unique(res$source)),
                    length(unique(res$condition))))
    res
  }
)

# Save the raw decoupleR result (the 09_progeny_activity.rds checkpoint).
# load_or_compute already saves it; also write explicitly so the output path
# matches the documented contract even on a cache-hit re-run.
saveRDS(progeny_raw,
        file.path(DIR_OBJECTS, "09_progeny_activity.rds"))
message("  09_progeny_activity.rds written (decoupleR MLM result, all contrasts).")

# ============================================================================
# 4. BUILD TIDY PATHWAY x CONTRAST TABLE
#    Apply BH padj across pathways WITHIN each contrast (matching 03_decoupler_tf.R
#    BH-within-contrast convention). With n=14 pathways per contrast, padj and
#    pvalue will be close but BH is still applied for schema consistency.
#    direction: "Up" if score > 0 (pathway more active in numerator condition).
# ============================================================================

message("[4] Building tidy progeny_activity table (BH within contrast) ...")

# progeny_raw columns: statistic / source / condition / score / p_value
progeny_tidy <- progeny_raw %>%
  dplyr::filter(!is.na(score), condition %in% CONTRASTS) %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(padj = stats::p.adjust(p_value, method = PADJM)) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(
    pathway_id      = paste0("PROGENY_", source),
    pathway_name    = source,
    database        = "PROGENy",
    contrast        = condition,
    nes             = score,      # "nes" = lowercase, schema-pinned master convention
    pvalue          = p_value,
    padj            = padj,
    direction       = ifelse(score > 0, "Up", "Down")
  ) %>%
  dplyr::arrange(contrast, pvalue)

# Annotate set_size (footprint genes present in our DE universe, per pathway).
# This is the footprint membership intersection — contrast-invariant (same gene set
# for every contrast, as MLM has no leading edge).
footprint_sizes <- net_progeny %>%
  dplyr::filter(target %in% universe) %>%
  dplyr::count(source, name = "set_size") %>%
  dplyr::rename(pathway_name = source)

progeny_tidy <- dplyr::left_join(progeny_tidy, footprint_sizes,
                                 by = "pathway_name")

# Annotate core_enrichment — the full footprint membership (slash-joined) as a
# contrast-invariant gene-set string. Used by downstream viz and schema validation.
footprint_genes <- net_progeny %>%
  dplyr::filter(target %in% universe) %>%
  dplyr::group_by(source) %>%
  dplyr::summarise(core_enrichment = paste(sort(unique(target)), collapse = "/"),
                   .groups = "drop") %>%
  dplyr::rename(pathway_name = source)

progeny_tidy <- dplyr::left_join(progeny_tidy, footprint_genes,
                                 by = "pathway_name")

# Final column order (consistent with master_gsea_table schema from config):
#   pathway_id, pathway_name, database, nes, pvalue, padj,
#   set_size, core_enrichment, contrast, direction
progeny_tidy <- progeny_tidy %>%
  dplyr::select(pathway_id, pathway_name, database,
                nes, pvalue, padj, set_size, core_enrichment,
                contrast, direction)

message(sprintf("  progeny_tidy: %d rows (%d pathways x %d contrasts)",
                nrow(progeny_tidy),
                length(unique(progeny_tidy$pathway_name)),
                length(unique(progeny_tidy$contrast))))

# ============================================================================
# 5. WRITE STAGE TABLE  — 03_results/05_progeny/tables/progeny_activity.csv
#    Tidy format: pathway x contrast, one row per (pathway, contrast).
#    This is the plot-ready input for 09_activity_progeny_viz.R.
# ============================================================================

message("[5] Writing 03_results/05_progeny/tables/progeny_activity.csv ...")
readr::write_csv(round_numeric_cols(progeny_tidy),
                 file.path(tbl_dir, "progeny_activity.csv"))
message(sprintf("  Saved: %d rows x %d cols",
                nrow(progeny_tidy), ncol(progeny_tidy)))

# ============================================================================
# 6. MASTER TABLE — master_progeny_activities.csv
#    Matches the master_gsea_table schema (YAML_CONFIG$schemas$master_gsea_table).
#    append_master_table() (from de_gsea_helpers.R) is idempotent: it drops prior
#    "PROGENy" rows and re-writes, so re-running this script is safe.
# ============================================================================

message("[6] Updating master/master_progeny_activities.csv (idempotent) ...")

# Validate schema requirements before writing
req_cols <- YAML_CONFIG$schemas$master_gsea_table$required_columns
if (!is.null(req_cols)) {
  missing <- setdiff(req_cols, colnames(progeny_tidy))
  if (length(missing) > 0)
    stop("progeny_tidy missing required schema cols: ", paste(missing, collapse = ", "),
         " — fix the column-assembly above.")
}
# Enforce lowercase nes (schema-pinned; explorer contract)
if ("NES" %in% colnames(progeny_tidy))
  stop("progeny_tidy has uppercase NES — rename to nes before writing.")

# Write the standalone PROGENy master (sole writer; overwrite each run)
readr::write_csv(round_numeric_cols(progeny_tidy),
                 file.path(DIR_MASTER, "master_progeny_activities.csv"))
message(sprintf("  master_progeny_activities.csv: %d rows", nrow(progeny_tidy)))

# Also idempotent-append into the canonical master_gsea_table accumulator under
# database="PROGENy" so the synthesis can join GSEA + PROGENy in one table (the
# per-database key keeps PROGENy rows disjoint from the MSigDB/custom/CoReSh rows).
append_master_table(progeny_tidy,
                    master_path = "master_gsea_table.csv",
                    key_col     = "database")
message("  master_gsea_table.csv updated (PROGENy rows replaced).")

# ============================================================================
# 7. SANITY CHECKS  — cGAS-dependence asymmetry positive controls
#    (load-bearing biology: must not be silent on reversal)
# ============================================================================

message("[7] Sanity checks (PROGENy biology positive controls) ...")

# 7a. PROGENy-Hypoxia in WT_heat and KO_heat (expected UP) and Interaction (expected flat/NS).
#     This is the key claim: pseudohypoxia/glycolysis program is cGAS-INDEPENDENT at n=5.
check_pathway <- function(pw, contrast) {
  r <- progeny_tidy %>%
    dplyr::filter(pathway_name == pw, contrast == !!contrast)
  if (nrow(r) == 0) return(NULL)
  r[1, ]
}

hyp_wt   <- check_pathway("Hypoxia", "WT_heat")
hyp_ko   <- check_pathway("Hypoxia", "KO_heat")
hyp_int  <- check_pathway("Hypoxia", "Interaction")
jak_wt   <- check_pathway("JAK-STAT", "WT_heat")
jak_ko   <- check_pathway("JAK-STAT", "KO_heat")
jak_int  <- check_pathway("JAK-STAT", "Interaction")

if (!is.null(hyp_wt))
  message(sprintf("  PROGENy Hypoxia  | WT_heat:     nes=%+.3f  pvalue=%.3e  padj=%.3e  direction=%s",
                  hyp_wt$nes, hyp_wt$pvalue, hyp_wt$padj, hyp_wt$direction))
if (!is.null(hyp_ko))
  message(sprintf("  PROGENy Hypoxia  | KO_heat:     nes=%+.3f  pvalue=%.3e  padj=%.3e  direction=%s",
                  hyp_ko$nes, hyp_ko$pvalue, hyp_ko$padj, hyp_ko$direction))
if (!is.null(hyp_int))
  message(sprintf("  PROGENy Hypoxia  | Interaction: nes=%+.3f  pvalue=%.3e  padj=%.3e  direction=%s",
                  hyp_int$nes, hyp_int$pvalue, hyp_int$padj, hyp_int$direction))

if (!is.null(jak_wt))
  message(sprintf("  PROGENy JAK-STAT | WT_heat:     nes=%+.3f  pvalue=%.3e  padj=%.3e  direction=%s",
                  jak_wt$nes, jak_wt$pvalue, jak_wt$padj, jak_wt$direction))
if (!is.null(jak_ko))
  message(sprintf("  PROGENy JAK-STAT | KO_heat:     nes=%+.3f  pvalue=%.3e  padj=%.3e  direction=%s",
                  jak_ko$nes, jak_ko$pvalue, jak_ko$padj, jak_ko$direction))
if (!is.null(jak_int))
  message(sprintf("  PROGENy JAK-STAT | Interaction: nes=%+.3f  pvalue=%.3e  padj=%.3e  direction=%s",
                  jak_int$nes, jak_int$pvalue, jak_int$padj, jak_int$direction))

# Warn (not stop) on direction reversals — let the scientist decide, but be loud.
if (!is.null(hyp_wt) && hyp_wt$nes < 0)
  warning("SANITY: PROGENy Hypoxia nes < 0 in WT_heat — expected UP (heat-induced glycolytic/HIF overlap).")
if (!is.null(jak_wt) && jak_wt$nes < 0)
  warning("SANITY: PROGENy JAK-STAT nes < 0 in WT_heat — expected UP (STING -> IFN positive control).")
if (!is.null(jak_int) && jak_int$nes < 0)
  warning("SANITY: PROGENy JAK-STAT nes < 0 in Interaction — expected positive (IFN arm is cGAS-dependent).")

# 7b. Full summary by contrast and significance
message("")
message("  --- PROGENy activity summary (n significant at padj<", FDR, ") ---")
summary_df <- progeny_tidy %>%
  dplyr::group_by(contrast) %>%
  dplyr::summarise(
    n_tested   = dplyr::n(),
    n_sig      = sum(padj < FDR, na.rm = TRUE),
    n_up_sig   = sum(direction == "Up"   & padj < FDR, na.rm = TRUE),
    n_down_sig = sum(direction == "Down" & padj < FDR, na.rm = TRUE),
    top_pw     = pathway_name[which.max(abs(nes))],
    top_nes    = nes[which.max(abs(nes))],
    .groups    = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_sig))

for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  message(sprintf("  %-14s  n_tested=%d  n_sig=%d  (up=%d,dn=%d)  top=%-10s nes=%+.2f",
                  r$contrast, r$n_tested, r$n_sig,
                  r$n_up_sig, r$n_down_sig, r$top_pw, r$top_nes))
}

# 7c. KEY_PROGENY watchlist across headline contrasts
message("")
message("  --- Key pathways (KEY_PROGENY) across headline contrasts ---")
key_rows <- progeny_tidy %>%
  dplyr::filter(pathway_name %in% KEY_PROGENY,
                contrast %in% c("WT_heat", "KO_heat", "Interaction", "Temp_main")) %>%
  dplyr::arrange(pathway_name,
                 factor(contrast, levels = c("WT_heat", "KO_heat", "Interaction", "Temp_main")))
if (nrow(key_rows) > 0) {
  for (i in seq_len(nrow(key_rows))) {
    r <- key_rows[i, ]
    message(sprintf("  %-10s | %-12s  nes=%+.3f  padj=%.3e  %s",
                    r$pathway_name, r$contrast, r$nes, r$padj,
                    ifelse(r$padj < FDR, paste0("SIG(", r$direction, ")"), "NS")))
  }
} else {
  warning("SANITY: KEY_PROGENY pathways (", paste(KEY_PROGENY, collapse = ", "),
          ") not found in progeny_tidy — check pathway name spelling vs net_progeny$source.")
}

# ============================================================================
# 8. ACCEPTANCE CHECKS
# ============================================================================

message("[8] Acceptance checks ...")

# Output files exist
stopifnot(
  "09_progeny_activity.rds missing" =
    file.exists(file.path(DIR_OBJECTS, "09_progeny_activity.rds")),
  "progeny_activity.csv missing" =
    file.exists(file.path(tbl_dir, "progeny_activity.csv")),
  "master_progeny_activities.csv missing" =
    file.exists(file.path(DIR_MASTER, "master_progeny_activities.csv"))
)

# Row count: 14 pathways x 7 contrasts (minus any dropped by minsize, but PROGENy
# footprints are pre-curated and all 14 should survive minsize=5)
n_expected <- n_pathways * length(CONTRASTS)
if (nrow(progeny_tidy) != n_expected)
  warning("Expected ", n_expected, " rows (", n_pathways, " pathways x ",
          length(CONTRASTS), " contrasts); got ", nrow(progeny_tidy),
          " — some pathways may have been dropped by minsize=", MINSIZE, ".")

# Schema: no uppercase NES
prog_cols <- colnames(readr::read_csv(
  file.path(DIR_MASTER, "master_progeny_activities.csv"),
  n_max = 1, show_col_types = FALSE))
if ("NES" %in% prog_cols)
  stop("master_progeny_activities.csv has uppercase NES — must be lowercase nes.")

message("  All acceptance checks passed.")
message("=================================================================")
message("09_activity_progeny COMPLETE")
message("  Objects: ", file.path(DIR_OBJECTS, "09_progeny_activity.rds"))
message("  Table:   ", file.path(tbl_dir, "progeny_activity.csv"))
message("  Master:  ", file.path(DIR_MASTER, "master_progeny_activities.csv"))
message("  Rows:    ", nrow(progeny_tidy), " (", n_pathways, " pathways x ",
        length(CONTRASTS), " contrasts)")
message("  Run 09_activity_progeny_viz.R for figures.")
message("=================================================================")
