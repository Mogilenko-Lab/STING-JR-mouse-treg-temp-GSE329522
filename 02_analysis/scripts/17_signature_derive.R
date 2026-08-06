# 17_signature_derive.R — COMPUTE
# =============================================================================
# Projection-signature derivation (stage 10_signature) — reshape the finished per-contrast
# DE master into projectable mouse-symbol gene sets + ranked lists.
#
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   10_signature   (script 17; stage id differs from script number, matching 12–16 → 06/07)
#
# ROLE: COMPUTE ONLY. Reshapes the published master_de_genes.csv into the projection
#   primitives Phase 0 freezes: per-contrast binary up/down sets (at BOTH the fdr_only and
#   fdr_logfc gates) and full signed-t ranked lists, all in MOUSE symbols. It also runs a
#   DISPLAY-ONLY dry-run ortholog-coverage preview (one babelgene call) so the human can
#   judge mapping loss at BREAKPOINT 10 before stage 18 freezes the contract. Figures live
#   in the sibling 17_signature_derive_viz.R, which reads only what this script writes.
#
# WHY THIS STAGE: the human compartments score
#   against a single frozen ortholog-mapped signature artifact, in place of a 7-contrast
#   mouse-symbol master. Stage 17 derives the DEFINITIONS in mouse symbols, the human
#   approves them at BREAKPOINT 10, and stage 18 freezes the human contract. The split is
#   deliberate: derive → approve → freeze.
#
# SIGN CONVENTION (preserved end to end): up/down split by DE direction (sign of logFC);
#   the ranked list ranks on the SIGNED limma-trend t-statistic (RANK_METRIC = "t", per
#   gsea.rank_metric). Correlative framing throughout — this is the INPUT the human phases
#   test for presence, and fever/HIF/STING causality is a separate claim.
#
# Input (read-only, the SOLE input):
#   03_results/master/master_de_genes.csv   (7 contrasts x 19,679 genes; cols gene_symbol,
#     ensembl, logFC, t, P.Value, adj.P.Val, contrast, significant, direction). The binary
#     sets are RE-DERIVED from the raw columns against the config thresholds. The master's
#     own `significant` flag is FDR-only, so the logFC gate is applied explicitly here and
#     stays auditable.
#
# Params (from config):
#   thresholds.de_fdr (0.05), thresholds.de_logfc (1.0); RANK_METRIC = "t".
#   decisions.projection.role_primary (which contrast is the headline; role column).
#   decisions.projection.ortholog_ambiguity.min_support (babelgene vote floor; preview only).
#
# Outputs:
#   03_results/objects/17_signature_sets.rds                 (the checkpoint; see shape below)
#   03_results/10_signature/tables/_overview/signature_sizes.csv           (contrast x role x direction x gate -> n_genes)
#   03_results/10_signature/tables/_overview/updown_overlap.csv            (pairwise Jaccard of up- and down-sets)
#   03_results/10_signature/tables/_overview/ortholog_coverage_preview.csv (contrast x gate x {mapped_1to1, one_to_many, unmapped})
#
# IDEMPOTENT + BYTE-STABLE: read -> reshape -> round -> write over a fixed master, plus an
#   offline babelgene lookup. Re-running yields the same bytes (round_numeric_cols).
#
# Run from project root:
#   Rscript 02_analysis/scripts/17_signature_derive.R
# =============================================================================

source("02_analysis/config/config.R")            # PROJECT_ROOT, YAML_CONFIG, DIR_*, %||%, RANK_METRIC, DE_FDR, DE_LOGFC
source("02_analysis/helpers/de_gsea_helpers.R")   # round_numeric_cols (compute idiom; NO plotting)
source("02_analysis/helpers/ortholog_utils.R")    # build_ortholog_map, ortholog_coverage (dry-run preview)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================================
# CONSTANTS (from config; never hardcoded — AGENTS.md rule 1)
# ============================================================================

STAGE  <- "10_signature"
SCRIPT <- "02_analysis/scripts/17_signature_derive.R"

DE_FDR_   <- as.numeric(YAML_CONFIG$thresholds$de_fdr   %||% 0.05)
DE_LOGFC_ <- as.numeric(YAML_CONFIG$thresholds$de_logfc %||% 1.0)
RANK_MET  <- RANK_METRIC %||% (YAML_CONFIG$gsea$rank_metric %||% "t")   # "t"

# Contrast display order = config order (WT_heat first — the primary axis).
CONTRASTS <- vapply(YAML_CONFIG$design$contrasts, function(c) c$name, character(1))

# Role assignment, config-driven. role=primary is the
# single headline axis; every other contrast is a comparator. This is the manifest
# `role` column stage 18 inherits — derived here for the display tables.
ROLE_PRIMARY <- unlist(YAML_CONFIG$decisions$projection$role_primary %||% list("WT_heat"))
role_of <- function(co) ifelse(co %in% ROLE_PRIMARY, "primary", "comparator")

# The two gates the notebook shows side by side. fdr_logfc is stringent (de_logfc=1.0
# on log2 CPM) and may decimate small contrasts; fdr_only keeps the |logFC| gate off.
GATES <- c("fdr_only", "fdr_logfc")

# babelgene vote floor for the DRY-RUN preview only (does NOT feed stage 18, which
# re-derives the map on the frozen sets at the decision-time min_support).
MIN_SUPPORT <- as.integer(YAML_CONFIG$decisions$projection$ortholog_ambiguity$min_support %||% 3L)

tbl_dir <- stage_dir(STAGE, "tables")             # 03_results/10_signature/tables/ (registered stage)
ov_dir  <- file.path(tbl_dir, YAML_CONFIG$figures$overview_dir %||% "_overview")
dir.create(ov_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(DIR_OBJECTS, recursive = TRUE, showWarnings = FALSE)

message("=================================================================")
message("17_signature_derive: reshape master_de_genes.csv -> projection sets (10_signature)")
message("  COMPUTE ONLY — no new statistics. Gates: ", paste(GATES, collapse = " / "),
        sprintf("  (de_fdr=%.3g, de_logfc=%.3g, rank='%s')", DE_FDR_, DE_LOGFC_, RANK_MET))
message("=================================================================")

# ============================================================================
# 1. READ the sole input (read-only) and validate the schema.
# ============================================================================

de_path <- file.path(DIR_MASTER, "master_de_genes.csv")
if (!file.exists(de_path))
  stop("[17] master_de_genes.csv not found at ", de_path,
       " — run the DE arm (02_de_limma_trend.R -> 16_synthesis.R) first.")

de <- readr::read_csv(de_path, show_col_types = FALSE, progress = FALSE)
need <- c("gene_symbol", "logFC", "t", "P.Value", "adj.P.Val", "contrast")
miss <- setdiff(need, colnames(de))
if (length(miss)) stop("[17] master_de_genes.csv missing columns: ", paste(miss, collapse = ", "))

de$gene_symbol <- as.character(de$gene_symbol)
de <- de[!is.na(de$gene_symbol) & de$gene_symbol != "", , drop = FALSE]

present_contrasts <- intersect(CONTRASTS, unique(de$contrast))
missing_contrasts <- setdiff(CONTRASTS, unique(de$contrast))
if (length(missing_contrasts))
  warning("[17] contrasts in config but ABSENT from master: ",
          paste(missing_contrasts, collapse = ", "), call. = FALSE)
message(sprintf("  [read] master_de_genes.csv: %d rows, %d contrasts (%s).",
                nrow(de), length(present_contrasts), paste(present_contrasts, collapse = ", ")))

# ============================================================================
# 2. RESHAPE — per-contrast binary up/down sets (both gates) + signed-t ranked list.
#    NO statistics: binary membership is a threshold on the existing adj.P.Val/logFC;
#    the ranked list is the existing signed t. `significant` from the master is NOT
#    used (it is FDR-only) — the gate is re-derived explicitly here.
# ============================================================================

# Binary set for one contrast subframe at one gate. Direction is the sign of logFC;
# fdr_only keeps the |logFC| gate off, fdr_logfc adds |logFC| >= de_logfc.
binary_sets <- function(d, gate) {
  sig <- !is.na(d$adj.P.Val) & d$adj.P.Val < DE_FDR_
  if (identical(gate, "fdr_logfc")) {
    up   <- sig & !is.na(d$logFC) & d$logFC >=  DE_LOGFC_
    down <- sig & !is.na(d$logFC) & d$logFC <= -DE_LOGFC_
  } else {                                   # fdr_only: split by direction, no |logFC| gate
    up   <- sig & !is.na(d$logFC) & d$logFC > 0
    down <- sig & !is.na(d$logFC) & d$logFC < 0
  }
  list(up   = sort(unique(d$gene_symbol[up])),
       down = sort(unique(d$gene_symbol[down])))
}

# Full signed-t ranked list for one contrast subframe: (gene_symbol, t) descending,
# NA/inf t dropped, duplicate symbols deduped keeping the most extreme |t|.
ranked_list <- function(d) {
  v <- suppressWarnings(as.numeric(d[[RANK_MET]]))
  keep <- !is.na(v) & is.finite(v)
  r <- data.frame(gene_symbol = d$gene_symbol[keep], t = v[keep], stringsAsFactors = FALSE)
  if (anyDuplicated(r$gene_symbol)) {
    r <- r[order(-abs(r$t)), , drop = FALSE]
    r <- r[!duplicated(r$gene_symbol), , drop = FALSE]
  }
  r[order(-r$t), , drop = FALSE]
}

sets <- list()
for (co in present_contrasts) {
  d <- de[de$contrast == co, , drop = FALSE]
  gate_sets <- lapply(setNames(GATES, GATES), function(g) binary_sets(d, g))
  sets[[co]] <- list(
    up     = lapply(gate_sets, `[[`, "up"),     # up[[gate]]
    down   = lapply(gate_sets, `[[`, "down"),   # down[[gate]]
    ranked = ranked_list(d)                     # data.frame(gene_symbol, t)
  )
  message(sprintf("  [%s] up/down  fdr_only=%d/%d  fdr_logfc=%d/%d  ranked=%d",
                  co,
                  length(sets[[co]]$up$fdr_only),  length(sets[[co]]$down$fdr_only),
                  length(sets[[co]]$up$fdr_logfc), length(sets[[co]]$down$fdr_logfc),
                  nrow(sets[[co]]$ranked)))
}

# ============================================================================
# 3. OVERVIEW TABLE (a): signature_sizes — contrast x role x direction x gate -> n_genes.
# ============================================================================

size_rows <- list()
for (co in present_contrasts) for (g in GATES) {
  size_rows[[length(size_rows) + 1L]] <- data.frame(
    contrast = co, role = role_of(co), direction = "up",   gate = g,
    n_genes = length(sets[[co]]$up[[g]]),   stringsAsFactors = FALSE)
  size_rows[[length(size_rows) + 1L]] <- data.frame(
    contrast = co, role = role_of(co), direction = "down", gate = g,
    n_genes = length(sets[[co]]$down[[g]]), stringsAsFactors = FALSE)
}
signature_sizes <- dplyr::bind_rows(size_rows)
signature_sizes$contrast <- factor(signature_sizes$contrast, levels = present_contrasts)
signature_sizes <- signature_sizes %>%
  dplyr::arrange(contrast, gate, direction) %>%
  dplyr::mutate(contrast = as.character(contrast))

# ============================================================================
# 4. OVERVIEW TABLE (b): updown_overlap — pairwise Jaccard of up-sets and of
#    down-sets across contrasts (per gate). Shows how DISTINCT WT_heat / Temp_main /
#    Geno_at_39 / Interaction are. Full symmetric matrix (self-pairs = 1) so the viz
#    can render a clean heatmap.
# ============================================================================

jaccard <- function(a, b) {
  if (length(a) == 0L && length(b) == 0L) return(NA_real_)
  length(intersect(a, b)) / length(union(a, b))
}
ov_rows <- list()
for (g in GATES) for (dir in c("up", "down")) {
  getset <- function(co) sets[[co]][[dir]][[g]]
  for (a in present_contrasts) for (b in present_contrasts) {
    sa <- getset(a); sb <- getset(b)
    ov_rows[[length(ov_rows) + 1L]] <- data.frame(
      contrast_a = a, contrast_b = b, direction = dir, gate = g,
      n_a = length(sa), n_b = length(sb),
      n_intersect = length(intersect(sa, sb)),
      jaccard = jaccard(sa, sb), stringsAsFactors = FALSE)
  }
}
updown_overlap <- dplyr::bind_rows(ov_rows)

# ============================================================================
# 5. OVERVIEW TABLE (c): ortholog_coverage_preview — DRY RUN (display only, does NOT
#    feed stage 18). One babelgene map over the whole DE universe, then per
#    (contrast x gate) coverage on the union of the up+down significant genes: how
#    many map 1:1, one-mouse->many-human, unmapped. Lets the human judge mapping
#    LOSS before the contract is frozen.
# ============================================================================

message("[dry-run] building ortholog map over the DE universe (babelgene min_support=",
        MIN_SUPPORT, ") — preview only, does NOT feed stage 18 ...")
universe <- sort(unique(de$gene_symbol))
omap_universe <- build_ortholog_map(universe, min_support = MIN_SUPPORT)
message(sprintf("  [dry-run] universe=%d mouse symbols; map edges=%d; distinct mapped mouse genes=%d.",
                length(universe), nrow(omap_universe),
                length(unique(omap_universe$mouse_symbol))))

cov_rows <- list()
for (co in present_contrasts) for (g in GATES) {
  updown <- union(sets[[co]]$up[[g]], sets[[co]]$down[[g]])
  cov <- ortholog_coverage(updown, omap = omap_universe, min_support = MIN_SUPPORT)
  cov_rows[[length(cov_rows) + 1L]] <- data.frame(
    contrast = co, role = role_of(co), gate = g,
    n_input = cov$n_input, mapped_1to1 = cov$mapped_1to1,
    one_to_many = cov$one_to_many, unmapped = cov$unmapped,
    stringsAsFactors = FALSE)
}
ortholog_coverage_preview <- dplyr::bind_rows(cov_rows) %>%
  dplyr::arrange(factor(contrast, levels = present_contrasts), gate)

# ============================================================================
# 6. WRITE the three overview tables (byte-stable) + the checkpoint.
# ============================================================================

message("[write] overview tables + 17_signature_sets.rds ...")
readr::write_csv(round_numeric_cols(signature_sizes),
                 file.path(ov_dir, "signature_sizes.csv"))
readr::write_csv(round_numeric_cols(updown_overlap),
                 file.path(ov_dir, "updown_overlap.csv"))
readr::write_csv(round_numeric_cols(ortholog_coverage_preview),
                 file.path(ov_dir, "ortholog_coverage_preview.csv"))

signature_sets <- list(
  contrasts       = present_contrasts,
  roles           = setNames(vapply(present_contrasts, role_of, character(1)), present_contrasts),
  gates           = GATES,
  thresholds      = list(de_fdr = DE_FDR_, de_logfc = DE_LOGFC_, rank_metric = RANK_MET),
  min_support     = MIN_SUPPORT,
  sets            = sets,                       # per-contrast up/down (both gates) + ranked
  signature_sizes = as.data.frame(signature_sizes),
  updown_overlap  = as.data.frame(updown_overlap),
  coverage_preview = as.data.frame(ortholog_coverage_preview),
  source_master   = "03_results/master/master_de_genes.csv",
  built_by        = SCRIPT,
  built_at        = format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC")
)
out_rds <- file.path(DIR_OBJECTS, "17_signature_sets.rds")
saveRDS(signature_sets, out_rds)

# ============================================================================
# 7. ACCEPTANCE CHECKS (structural; no statistics) + console read for the breakpoint.
# ============================================================================

stopifnot(
  "17_signature_sets.rds missing"    = file.exists(out_rds),
  "signature_sizes.csv missing"      = file.exists(file.path(ov_dir, "signature_sizes.csv")),
  "updown_overlap.csv missing"       = file.exists(file.path(ov_dir, "updown_overlap.csv")),
  "coverage_preview.csv missing"     = file.exists(file.path(ov_dir, "ortholog_coverage_preview.csv"))
)

# Console read: the candidate PRIMARY (WT_heat) set sizes + face-validity heat markers.
prim <- ROLE_PRIMARY[1]
if (prim %in% present_contrasts) {
  thermo <- c("Hspa1b", "Hsph1", "Hspa1a", "Dnajb1")
  up_logfc <- sets[[prim]]$up$fdr_logfc
  message(sprintf("  [primary=%s] fdr_logfc up=%d down=%d | fdr_only up=%d down=%d",
                  prim, length(up_logfc), length(sets[[prim]]$down$fdr_logfc),
                  length(sets[[prim]]$up$fdr_only), length(sets[[prim]]$down$fdr_only)))
  hit <- intersect(thermo, up_logfc)
  message(sprintf("  [face-validity] heat markers in %s up (fdr_logfc): %s",
                  prim, if (length(hit)) paste(hit, collapse = ", ") else "(none — inspect at breakpoint)"))
}

message("=================================================================")
message("17_signature_derive COMPLETE")
message("  Checkpoint: ", out_rds)
message("  Tables:     ", ov_dir, "/{signature_sizes,updown_overlap,ortholog_coverage_preview}.csv")
message("  Next: 17_signature_derive_viz.R, then render notebooks/17_signature_review/ (BREAKPOINT 10).")
message("  Do NOT run 18_projection_export.R until decisions.projection.status = APPROVED.")
message("=================================================================")
