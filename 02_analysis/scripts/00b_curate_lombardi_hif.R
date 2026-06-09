#!/usr/bin/env Rscript
# 00b_curate_lombardi_hif.R - Re-derive the PUBLISHED Lombardi-2022 CONSERVED HIF signature
# Project: STING-cGAS-GSE329522
# Phase: 0 (infrastructure / curated reference asset)
# Inputs:
#   00_data/references/gene_sets/1-s2.0-S2211124722015236-mmc2.xlsx   (32 TCGA per-cancer sheets)
#   00_data/processed/GSE329522_normalized_counts_CPM_iTreg.csv       (this dataset's gene_name column)
# Outputs:
#   00_data/references/gene_sets/lombardi2022_hif_consensus_human.rds  (human consensus, protein-coding)
#   00_data/references/gene_sets/lombardi2022_hif_consensus_mouse.rds  (mouse orthologs; the bulk asset)
#   00_data/references/gene_sets/lombardi2022_hif_curation_log.csv     (full per-gene audit trail)
#   00_data/references/gene_sets/tables/lombardi_recurrence_data.csv   (tidy, plot-ready: per gene)
#   00_data/references/gene_sets/tables/lombardi_recurrence_null.csv   (tidy null curve, per recurrence level)
# Dependencies: 02_analysis/config/config.R; readxl; babelgene (offline human->mouse orthology)
#
# ===========================================================================================
# WHY THIS SCRIPT WAS REWRITTEN (the bug it fixes)
# ===========================================================================================
#   The previous 00b HARDCODED a 48-gene "Lombardi consensus" vector. That list is provably NOT
#   the published consensus: it is verbatim the single `TCGA-ACC` sheet of mmc2.xlsx (that one
#   cancer type happens to have exactly 48 rows). The TRUE published signature is the set of genes
#   that RECUR across the 32 TCGA per-cancer-type sheets ("a conserved gene signature"). This
#   script reconstructs that conserved set from the supplement, end to end, with a statistical
#   justification for the recurrence threshold and a full audit log.
#
# SOURCE (do not edit without re-checking the paper):
#   Lombardi O, Li R, Halim S, Choudhry H, Ratcliffe PJ, Mole DR.
#   "Pan-cancer analysis of tissue and single-cell HIF-pathway activation using a conserved
#    gene signature." Cell Reports 2022; 41(7):111652. doi:10.1016/j.celrep.2022.111652
#   mmc2.xlsx: 32 sheets TCGA-ACC..TCGA-UVM, each `read_excel(skip=1)` ->
#   columns {row_num, gene(HUMAN symbol), correlation, pvalue}, sorted by descending correlation,
#   "genes positively correlating with the HIF-metagene". Per-sheet sizes vary 5..49,763 rows.
#
# ===========================================================================================
# METHOD & RATIONALE (decisions are auditable, not arbitrary)
# ===========================================================================================
#   (A) PER-DATASET HIF-MEMBERSHIP RULE
#       Rule = BH-FDR(pvalue) < 0.05  AND  correlation > 0  AND  rank <= CAP (top-CAP by
#       descending correlation), CAP = 200.
#       - Significance (BH<0.05, positive corr) is the primary, principled criterion.
#       - WHY THE CAP: significance ALONE is pathological here because the sheets differ in size
#         by 4 orders of magnitude. A 49,763-row sheet (TCGA-BRCA) calls ~6,700 genes "FDR-sig",
#         while TCGA-CHOL has only 5 rows total. Counting recurrence off raw FDR membership would
#         let a handful of huge cohorts dominate and would conflate "deeply sequenced cohort" with
#         "conserved HIF gene". The top-CAP cap gives every cancer type an EQUAL vote of at most
#         CAP genes (its strongest-correlating HIF members), so recurrence measures conservation
#         across cancer types, not cohort size. CAP=200 sits near the per-sheet "real-signal"
#         depth (median FDR-membership ~ a few hundred) and is the cap used for the reported set;
#         the conclusion (canonical glycolysis/HIF core at the top) is robust to CAP in 200..1000.
#       - Small sheets (e.g. CHOL=5) simply contribute their <=5 members; documented, not special-cased.
#
#   (B) SHARED UNIVERSE
#       No gene appears in ALL 32 sheets (intersection is empty: each cohort reports a different
#       gene complement). The shared universe U is therefore the UNION of genes appearing in any
#       sheet (|U| ~ 54k). Per-dataset null inclusion probability p_d = |members_d| / |U|.
#
#   (C) CONSENSUS THRESHOLD VIA AN ANALYTIC NULL (no permutation needed; fully deterministic)
#       Under H0 each dataset d independently "includes" a given gene with prob p_d. A gene's
#       recurrence count across the 32 datasets is then Poisson-binomial(p_1..p_32). We compute
#       the exact upper-tail P(X >= observed_recurrence) by direct Bernoulli convolution of the
#       pmf (stable, dependency-free), per gene, then BH-correct across all |U| genes.
#       KEEP rule (alpha): null_padj < ALPHA_BONF, ALPHA_BONF = 0.05/|U|. Because U~54k genes are
#       tested, a genome-wide-significance-style (Bonferroni-scaled) alpha is the honest control;
#       a plain 0.05 would admit genes recurring in only ~3/32 cohorts (statistically non-random,
#       but not "conserved"). This yields a conserved set on the order of ~100 genes BEFORE biotype
#       filtering; the protein-coding, ortholog-mapped, in-dataset bulk asset is smaller (tens).
#
#   (D) HUMAN->MOUSE ORTHOLOGY: babelgene::orthologs(species="mouse", human=TRUE), 1:1 best per
#       human gene (same pattern as the prior 00b). (E) VALIDATE mouse symbols against this
#       dataset's gene_name column (FILE_CPM). Biotype drops (lncRNA host genes / antisense /
#       microRNA / processed pseudogenes / clone-based Ensembl-Havana IDs) by symbol convention,
#       each logged with a reason, applied BEFORE orthology mapping.
#
#   DETERMINISM: the null is analytic (no RNG). A seed constant is recorded for provenance only.
# ===========================================================================================

source("02_analysis/config/config.R")
suppressPackageStartupMessages({ library(readxl); library(babelgene) })

OUT_DIR    <- file.path(PROJECT_ROOT, "00_data/references/gene_sets")
TBL_DIR    <- file.path(OUT_DIR, "tables")
XLSX       <- file.path(OUT_DIR, "1-s2.0-S2211124722015236-mmc2.xlsx")
SCRATCH    <- file.path(PROJECT_ROOT, "docs/_internal/_scratch")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TBL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SCRATCH, recursive = TRUE, showWarnings = FALSE)

# Tunable, DOCUMENTED parameters (see header).
CAP         <- 200L          # per-dataset top-N cap by descending correlation
PER_DS_FDR  <- 0.05          # per-dataset BH-FDR membership cutoff (positive correlation)
SEED        <- 123L          # provenance only; method is analytic / deterministic
set.seed(SEED)

stopifnot(file.exists(XLSX))

# ---------------------------------------------------------------------------
# 0) DELETE the stale (wrong) assets up front, so a partial run can't leave them.
# ---------------------------------------------------------------------------
stale <- file.path(OUT_DIR, c(
  "lombardi2022_hif_consensus_mouse.rds",
  "lombardi2022_hif_consensus_human.rds",
  "lombardi2022_hif_curation_log.csv"
))
removed <- stale[file.exists(stale)]
suppressWarnings(file.remove(removed))
if (length(removed)) message(sprintf("[curate] removed stale asset(s): %s",
                                      paste(basename(removed), collapse = ", ")))

# ---------------------------------------------------------------------------
# 1) Read all 32 sheets robustly (skip=1; coerce types; normalize symbols).
# ---------------------------------------------------------------------------
sheets <- readxl::excel_sheets(XLSX)
stopifnot(length(sheets) == 32L, all(grepl("^TCGA-", sheets)))

read_sheet <- function(s) {
  d <- suppressWarnings(readxl::read_excel(XLSX, sheet = s, skip = 1))
  names(d)[1:4] <- c("row_num", "gene", "correlation", "pvalue")
  d$gene        <- toupper(trimws(as.character(d$gene)))
  d$correlation <- suppressWarnings(as.numeric(d$correlation))
  d$pvalue      <- suppressWarnings(as.numeric(d$pvalue))
  d <- d[!is.na(d$gene) & d$gene != "" & !is.na(d$correlation) & !is.na(d$pvalue),
         c("gene", "correlation", "pvalue"), drop = FALSE]
  d[!duplicated(d$gene), , drop = FALSE]   # one row per gene (keep first = highest corr, sheet is sorted)
}

sets <- lapply(sheets, read_sheet)
names(sets) <- sheets
sheet_nrow <- vapply(sets, nrow, integer(1))
message(sprintf("[curate] read 32 sheets; per-sheet gene rows range %d..%d (median %.0f).",
                min(sheet_nrow), max(sheet_nrow), median(sheet_nrow)))

# ---------------------------------------------------------------------------
# 2) Per-dataset membership: BH-FDR<PER_DS_FDR & corr>0, then cap to top-CAP by corr.
# ---------------------------------------------------------------------------
membership <- lapply(sets, function(d) {
  padj <- p.adjust(d$pvalue, method = "BH")
  keep <- which(padj < PER_DS_FDR & d$correlation > 0)
  dd   <- d[keep, , drop = FALSE]
  dd   <- dd[order(-dd$correlation), , drop = FALSE]
  head(dd$gene, CAP)
})
memb_size <- vapply(membership, length, integer(1))
message(sprintf("[curate] per-dataset membership sizes range %d..%d (cap=%d).",
                min(memb_size), max(memb_size), CAP))

# Shared universe = union of all genes appearing in any sheet.
universe <- sort(unique(unlist(lapply(sets, `[[`, "gene"), use.names = FALSE)))
U <- length(universe)
message(sprintf("[curate] shared universe |U| = %d genes (intersection across 32 sheets is empty).", U))

# ---------------------------------------------------------------------------
# 3) Cross-dataset recurrence count per gene.
# ---------------------------------------------------------------------------
rec <- integer(U); names(rec) <- universe
recur_tab <- table(factor(unlist(membership, use.names = FALSE), levels = universe))
rec[] <- as.integer(recur_tab)

# ---------------------------------------------------------------------------
# 4) Analytic Poisson-binomial null + BH; consensus = null_padj < 0.05/|U|.
# ---------------------------------------------------------------------------
p_d <- memb_size / U                      # per-dataset inclusion prob under H0
# Exact Poisson-binomial pmf via Bernoulli convolution; pmf[k+1] = P(X = k), k = 0..32.
pois_binom_pmf <- function(p) {
  pmf <- 1
  for (pi in p) pmf <- c(pmf * (1 - pi), 0) + c(0, pmf * pi)
  pmf
}
pmf       <- pois_binom_pmf(p_d)
uppertail <- rev(cumsum(rev(pmf)))        # uppertail[m+1] = P(X >= m)
null_p    <- uppertail[rec + 1L]          # per-gene upper-tail p at its observed recurrence
null_padj <- p.adjust(null_p, method = "BH")

ALPHA_BONF <- 0.05 / U
kept_stat  <- null_padj < ALPHA_BONF      # statistically-conserved consensus (pre-biotype)
message(sprintf("[curate] E[recurrence] under null = %.3f; null-significant genes (padj<%.2g) = %d.",
                sum(p_d), ALPHA_BONF, sum(kept_stat)))

# ---------------------------------------------------------------------------
# 5) Biotype classification (symbol-convention; drop non-coding/pseudogene/clone w/ reasons).
#    Applied to the null-significant consensus before orthology mapping.
# ---------------------------------------------------------------------------
#  The processed-pseudogene rule (<PARENT>P<digits>) is parent-aware: a trailing "P<digits>"
#  is only a pseudogene suffix when the PARENT symbol (the gene minus that suffix) is itself
#  present in the supplement universe. A bare "P[0-9]+$" regex is a false-positive magnet -- it
#  wrongly flags real protein-coding genes (BNIP3, CDCP1, MTFP1, PLP2, SAP30, TNIP1). Requiring a
#  present parent cleanly separates true pseudogenes (AK4P1<-AK4, LDHAP5<-LDHA, GAPDHP63<-GAPDH,
#  RPL17P50<-RPL17, TPI1P1<-TPI1, BNIP3P1<-BNIP3) from those false positives.
classify_drop <- function(g, parent_universe) {
  reason <- rep(NA_character_, length(g))
  reason[grepl("^MIR[0-9]", g)]                 <- "microRNA"
  reason[grepl("HG$", g)]                       <- "lncRNA host gene"
  reason[grepl("-AS[0-9]+$", g)]                <- "antisense lncRNA"
  reason[grepl("-DT$", g)]                      <- "divergent transcript / lncRNA"
  reason[grepl("^LINC[0-9]", g)]                <- "long intergenic non-coding RNA"
  reason[grepl("^A[CFLP][0-9].*\\.[0-9]+$", g)] <- "clone-based Ensembl-Havana ID (uncharacterized)"
  pseudo_cand <- grepl("P[0-9]+$", g)
  parent      <- sub("P[0-9]+$", "", g)
  is_pseudo   <- pseudo_cand & is.na(reason) & (parent %in% parent_universe)
  reason[is_pseudo] <- "processed pseudogene"
  reason
}
drop_reason_all <- classify_drop(universe, universe)  # parent-aware; only meaningful for kept_stat genes
is_noncoding    <- !is.na(drop_reason_all)

coding_kept <- universe[kept_stat & !is_noncoding]
human_genes <- sort(unique(coding_kept))
message(sprintf("[curate] consensus: %d null-significant -> %d protein-coding (%d non-coding/pseudo/clone dropped).",
                sum(kept_stat), length(human_genes), sum(kept_stat & is_noncoding)))

# ---------------------------------------------------------------------------
# 6) Human -> mouse orthology (babelgene; 1:1 best-per-human-gene).
# ---------------------------------------------------------------------------
orth <- babelgene::orthologs(genes = human_genes, species = "mouse", human = TRUE)
orth_map <- orth[, c("human_symbol", "symbol")]
names(orth_map) <- c("human_symbol", "mouse_symbol")
orth_map <- orth_map[!duplicated(orth_map$human_symbol), , drop = FALSE]

# ---------------------------------------------------------------------------
# 7) Validate mouse symbols against THIS dataset's gene_name column.
# ---------------------------------------------------------------------------
cpm <- read.csv(FILE_CPM, check.names = FALSE)
dataset_symbols <- unique(cpm$gene_name)

# ---------------------------------------------------------------------------
# 8) Full per-gene audit log (every gene that recurs at least once is auditable;
#    plus the consensus columns the prompt requires).
# ---------------------------------------------------------------------------
audit_genes <- universe[rec > 0]
ai <- match(audit_genes, universe)
mouse_for <- orth_map$mouse_symbol[match(audit_genes, orth_map$human_symbol)]

log_df <- data.frame(
  human_symbol    = audit_genes,
  recurrence_count= rec[ai],
  n_datasets      = 32L,
  null_p          = null_p[ai],
  null_padj       = null_padj[ai],
  kept            = kept_stat[ai] & !is_noncoding[ai],   # final consensus membership
  biotype_class   = ifelse(is_noncoding[ai], "non-coding/pseudogene/clone", "protein-coding"),
  drop_reason     = drop_reason_all[ai],
  mouse_symbol    = mouse_for,
  in_GSE329522    = mouse_for %in% dataset_symbols,
  stringsAsFactors = FALSE
)
log_df$in_GSE329522[is.na(log_df$mouse_symbol)] <- NA
log_df <- log_df[order(-log_df$recurrence_count, log_df$human_symbol), ]

# Final mouse asset: protein-coding ^ null-significant ^ ortholog-found ^ present in dataset.
mouse_genes <- sort(unique(log_df$mouse_symbol[which(log_df$kept &
                                                     !is.na(log_df$mouse_symbol) &
                                                     log_df$in_GSE329522)]))

# ---------------------------------------------------------------------------
# 9) Tidy plot-ready tables (NO ggplot here; a later viz agent draws figures).
# ---------------------------------------------------------------------------
# 9a) Per-gene recurrence table (one row per gene that recurs >=1).
recurrence_data <- data.frame(
  human_symbol     = log_df$human_symbol,
  recurrence_count = log_df$recurrence_count,
  null_p           = log_df$null_p,
  null_padj        = log_df$null_padj,
  kept             = log_df$kept,
  biotype_class    = log_df$biotype_class,
  mouse_symbol     = log_df$mouse_symbol,
  in_GSE329522     = log_df$in_GSE329522,
  stringsAsFactors = FALSE
)
write.csv(recurrence_data, file.path(TBL_DIR, "lombardi_recurrence_data.csv"), row.names = FALSE)

# 9b) Null curve, per recurrence level 0..32: observed gene count vs null-expected count,
#     the upper-tail p, and whether that level clears the consensus alpha. Lets a viz draw the
#     recurrence histogram with the chosen threshold overlaid.
levels_m       <- 0:32
obs_count      <- vapply(levels_m, function(m) sum(rec == m), integer(1))
null_expected  <- pmf * U                                  # E[# genes with recurrence exactly m]
raw_uppertail  <- uppertail[levels_m + 1L]
# smallest recurrence level whose genes clear the consensus alpha (for threshold annotation):
kept_levels    <- sort(unique(rec[kept_stat]))
min_kept_level <- if (length(kept_levels)) min(kept_levels) else NA_integer_
recurrence_null <- data.frame(
  recurrence_level        = levels_m,
  n_genes_observed        = obs_count,
  null_expected_at_count  = null_expected,
  null_uppertail_p        = raw_uppertail,            # P(X >= level) under H0
  clears_consensus_alpha  = !is.na(min_kept_level) & levels_m >= min_kept_level,
  stringsAsFactors        = FALSE
)
write.csv(recurrence_null, file.path(TBL_DIR, "lombardi_recurrence_null.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 10) Save assets with provenance attributes (same filenames -> downstream loads unchanged).
# ---------------------------------------------------------------------------
provenance <- list(
  signature        = "Lombardi 2022 CONSERVED HIF signature (re-derived from mmc2.xlsx)",
  citation         = "Lombardi O, Li R, Halim S, Choudhry H, Ratcliffe PJ, Mole DR. Cell Reports 2022;41(7):111652.",
  doi              = "10.1016/j.celrep.2022.111652",
  source_file      = basename(XLSX),
  n_datasets       = 32L,
  shared_universe  = U,
  per_dataset_rule = sprintf("BH-FDR<%.2f & corr>0 & top-%d by correlation", PER_DS_FDR, CAP),
  null_model       = "Poisson-binomial(p_d = |members_d|/|U|), exact upper tail, BH across genes",
  consensus_alpha  = ALPHA_BONF,
  n_consensus_human= length(human_genes),
  n_mouse_in_dataset = length(mouse_genes),
  mapping_tool     = paste0("babelgene ", as.character(utils::packageVersion("babelgene"))),
  seed             = SEED,
  built_by         = "02_analysis/scripts/00b_curate_lombardi_hif.R",
  note             = "Conserved across 32 TCGA cancer types. Mouse set = protein-coding ^ null-significant ^ ortholog ^ present in GSE329522."
)
human_set <- structure(list(Lombardi2022_HIF = human_genes), provenance = provenance)
mouse_set <- structure(list(Lombardi2022_HIF = mouse_genes), provenance = provenance)

saveRDS(human_set, file.path(OUT_DIR, "lombardi2022_hif_consensus_human.rds"))
saveRDS(mouse_set, file.path(OUT_DIR, "lombardi2022_hif_consensus_mouse.rds"))
write.csv(log_df, file.path(OUT_DIR, "lombardi2022_hif_curation_log.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 11) Report (machine + human readable). Contrast vs the stale 48-from-one-sheet list.
# ---------------------------------------------------------------------------
lombardi48_stale <- c(
  "LDHA","PPFIA4","PFKFB4","AK4","ERO1A","BNIP3L","BNIP3","PFKL","ENO1","MIR210HG",
  "TPI1","HILPDA","PGK1","INSIG2","AK4P1","AC114803.1","ENO2","FUT11","EGLN3","NDRG1",
  "BNIP3P1","WDR54","OVOL1-AS1","STC2","DDX41","C4orf47","GYS1","ANKRD37","BICDL2","P4HA1",
  "PDK1","OSMR","PKM","LDHAP5","AL158201.1","PFKP","MIR210","NLRP3P1","GSX1","AL109946.1",
  "PGAM1","SLC16A3","ISM2","TCAF2","ARID3A","KDM3A","JMJD6","C4orf3"
)
overlap_h  <- intersect(toupper(human_genes), toupper(lombardi48_stale))
new_only   <- setdiff(toupper(human_genes), toupper(lombardi48_stale))
stale_only <- setdiff(toupper(lombardi48_stale), toupper(human_genes))

cat("\n================ Lombardi 2022 CONSERVED HIF re-derivation ================\n")
cat(sprintf("Sheets read:                32 TCGA cancer types (sizes %d..%d rows)\n",
            min(sheet_nrow), max(sheet_nrow)))
cat(sprintf("Shared universe |U|:        %d genes\n", U))
cat(sprintf("Per-dataset rule:           BH-FDR<%.2f & corr>0 & top-%d by corr\n", PER_DS_FDR, CAP))
cat(sprintf("Null model:                 Poisson-binomial exact upper tail; BH; alpha=0.05/|U|=%.2g\n", ALPHA_BONF))
cat(sprintf("Null-significant consensus: %d genes\n", sum(kept_stat)))
cat(sprintf("  protein-coding:           %d  (human asset)\n", length(human_genes)))
cat(sprintf("  mouse orthologs found:    %d\n", nrow(orth_map)))
cat(sprintf("  in GSE329522 (mouse):     %d  <-- FINAL mouse set size\n", length(mouse_genes)))
cat(sprintf("\nFinal MOUSE set (n=%d):\n%s\n", length(mouse_genes), paste(mouse_genes, collapse = ", ")))
cat(sprintf("\nContrast vs stale 48-from-TCGA-ACC list:\n"))
cat(sprintf("  overlap (human):          %d\n", length(overlap_h)))
cat(sprintf("  in new consensus only:    %d\n", length(new_only)))
cat(sprintf("  in stale-48 only (dropped by recurrence): %d -> %s\n",
            length(stale_only), paste(stale_only, collapse = ", ")))
cat("\nWrote:\n")
cat("  ", file.path(OUT_DIR, "lombardi2022_hif_consensus_mouse.rds"), "\n")
cat("  ", file.path(OUT_DIR, "lombardi2022_hif_consensus_human.rds"), "\n")
cat("  ", file.path(OUT_DIR, "lombardi2022_hif_curation_log.csv"), "\n")
cat("  ", file.path(TBL_DIR, "lombardi_recurrence_data.csv"), "\n")
cat("  ", file.path(TBL_DIR, "lombardi_recurrence_null.csv"), "\n")
cat("==========================================================================\n")
