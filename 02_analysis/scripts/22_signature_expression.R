# 22_signature_expression.R — COMPUTE
# =============================================================================
# Signature-gene expression across the 2x2 design (stage 14_signature_expression) —
# read the exported UP arms back out gene by gene, in the four design cells.
#
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   14_signature_expression   (script 22; stage id differs from script number,
#          matching 17/18 -> 10/11 and 20/21 -> 13)
#
# ROLE: COMPUTE ONLY. Every p-value it reports was computed in 02_de_limma_trend.R, and
#   every expression value is a group MEAN of the frozen 01_eda log2 CPM matrix. The
#   figure lives in the sibling 23_signature_expression_viz.R, which reads only what this
#   script writes.
#
# WHY THIS STAGE. Stages 10/11 report a signature as a COUNT: WT_heat_up is "213 mouse
#   symbols -> 199 human". That is the right unit for a projection contract. The bench
#   question is per gene: how strongly is this gene expressed at all, in which of the four
#   design cells is it high, and does removing cGAS change that picture. This stage
#   produces that substrate — a gene x design-cell table with group means, dispersion and
#   an across-group z, plus the DE statistics and ortholog fate each gene carries — so a
#   dot plot and its interactive twin can be drawn from tables alone.
#
# SCOPE: UP ARMS ONLY. Four signatures, all UP, all in MOUSE symbol space:
#   WT_heat_up             (fdr_logfc)  the full thermal response, cGAS-competent
#   KO_heat_up             (fdr_logfc)  the same response with cGAS removed
#   Interaction_up         (fdr_logfc)  the cGAS-dependent slice, 1 df, thin
#   Interaction_up_fdrOnly (fdr_only)   the SAME contrast at the looser gate
#
#   The DOWN arms sit outside this stage — out of the tables, the figure and the README
#   prose — as a scope decision recorded in analysis_config.yaml (signature_expression:).
#   The bench question is about the induced arm, and a two-headed dot plot invites reading
#   up and down as one program. Interaction_up_fdrOnly is flagged in EVERY row and facet
#   label as a gate-sensitivity read on the same 1-df contrast.
#
# Naming discipline (umbrella AGENTS.md): every signature is named by HOW IT WAS DERIVED
#   (contrast + direction + gate). "WT_heat_up" is checkable; a name taken from a
#   hoped-for mechanism smuggles in the conclusion.
#
# Inputs (read-only; each existence-checked, and a missing one is reported):
#   03_results/objects/01_eda.rds              $logcpm_mat (19,679 x 20) + $metadata (20 x 7)
#   03_results/objects/17_signature_sets.rds   $sets$<contrast>$up$<gate> (mouse symbols)
#   03_results/objects/02_de_results.rds       per-contrast limma topTables (logFC, t, adj.P.Val)
#   03_results/master/master_de_genes.csv      the published DE master — cross-checked against
#                                              02_de_results.rds on the signature genes; a
#                                              disagreement is a hard stop, reported as such.
#   03_results/human_projection/ortholog_map.tsv  the APPLIED frozen mouse->human map
#   03_results/human_projection/manifest.csv      the frozen projection ledger — every set size
#                                              and ortholog fate is cross-checked against it.
#
# Params (from config):
#   signature_expression.signatures[]  (name, source_contrast, direction, gate)
#   signature_expression.top_n_genes   (plotted rows per facet; the tables stay complete)
#   design.samples_per_group, design.groups, thresholds.de_fdr/de_logfc, gsea.rank_metric
#
# Outputs (tables only):
#   03_results/14_signature_expression/tables/signature_gene_group_expression.csv
#   03_results/14_signature_expression/tables/signature_gene_stats.csv
#   03_results/14_signature_expression/tables/signature_expression_summary.csv
#   03_results/14_signature_expression/tables/_overview/signature_dotplot.csv
#
# IDEMPOTENT + BYTE-STABLE: read -> group-mean -> join -> round -> write over frozen inputs.
#   Re-running yields the same bytes (round_numeric_cols at 9 significant figures).
#
# Run from project root:
#   Rscript 02_analysis/scripts/22_signature_expression.R
# =============================================================================

source("02_analysis/config/config.R")           # PROJECT_ROOT, YAML_CONFIG, DIR_*, %||%, RANK_METRIC
source("02_analysis/helpers/de_gsea_helpers.R")  # round_numeric_cols (compute idiom; NO plotting)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================================
# CONSTANTS (from config; never hardcoded — AGENTS.md rule 1)
# ============================================================================

STAGE  <- "14_signature_expression"
SCRIPT <- "02_analysis/scripts/22_signature_expression.R"

SE_CFG <- YAML_CONFIG$signature_expression
if (is.null(SE_CFG))
  stop("[22] analysis_config.yaml::signature_expression is missing — register the stage's ",
       "signature list / top_n_genes there before running (config, not hardcoding).")

# The signature roster. Each entry: name, source_contrast, direction, gate.
SIGS <- lapply(SE_CFG$signatures, function(s) list(
  name            = as.character(s$name),
  source_contrast = as.character(s$source_contrast),
  direction       = as.character(s$direction),
  gate            = as.character(s$gate)))
if (!length(SIGS)) stop("[22] signature_expression.signatures is empty.")

# UP-ARM-ONLY GUARD. The scope decision is enforced in code: if a future
# edit adds a down arm to the config roster this script refuses to run rather than
# silently widening a deliberately narrow deliverable.
bad_dir <- vapply(SIGS, function(s) !identical(s$direction, "up"), logical(1))
if (any(bad_dir))
  stop("[22] stage 14 is UP ARMS ONLY by scope decision, but signature_expression.signatures ",
       "contains direction != 'up': ",
       paste(vapply(SIGS[bad_dir], `[[`, character(1), "name"), collapse = ", "),
       ". Remove it, or open a new stage for the down arms — do not widen this one silently.")

TOP_N       <- as.integer(SE_CFG$top_n_genes %||% 25L)
EXPR_SLOT   <- as.character(SE_CFG$expression_matrix %||% "logcpm_mat")
N_PER_GROUP <- as.integer(YAML_CONFIG$design$samples_per_group %||% 5L)
CFG_GROUPS  <- unlist(YAML_CONFIG$design$groups %||% list())
DE_FDR_     <- as.numeric(YAML_CONFIG$thresholds$de_fdr   %||% 0.05)
DE_LOGFC_   <- as.numeric(YAML_CONFIG$thresholds$de_logfc %||% 1.0)
RANK_MET    <- RANK_METRIC %||% (YAML_CONFIG$gsea$rank_metric %||% "t")   # "t"

p_root <- function(rel) file.path(PROJECT_ROOT, rel)
F_EDA      <- p_root(SE_CFG$expression_object   %||% "03_results/objects/01_eda.rds")
F_SETS     <- p_root(SE_CFG$source_object       %||% "03_results/objects/17_signature_sets.rds")
F_DE_RDS   <- p_root("03_results/objects/02_de_results.rds")
F_DE_CSV   <- file.path(DIR_MASTER, "master_de_genes.csv")
F_OMAP     <- p_root(SE_CFG$ortholog_map        %||% "03_results/human_projection/ortholog_map.tsv")
F_MANIFEST <- p_root(SE_CFG$projection_manifest %||% "03_results/human_projection/manifest.csv")

tbl_dir <- stage_dir(STAGE, "tables")           # 03_results/14_signature_expression/tables/
ov_dir  <- file.path(tbl_dir, YAML_CONFIG$figures$overview_dir %||% "_overview")
dir.create(ov_dir, recursive = TRUE, showWarnings = FALSE)

message("=================================================================")
message("22_signature_expression: exported UP arms -> per-gene x design-cell expression")
message("  COMPUTE ONLY — no new statistics. UP ARMS ONLY (down arms out of scope by design).")
message("  Signatures: ", paste(vapply(SIGS, function(s)
        sprintf("%s [%s/%s]", s$name, s$source_contrast, s$gate), character(1)), collapse = "  "))
message(sprintf("  top_n_genes=%d (figure only; tables carry every member)", TOP_N))
message("=================================================================")

# ============================================================================
# 1. GUARD every input. A missing input is REPORTED with its path and stops the run.
# ============================================================================

inputs <- c(eda = F_EDA, signature_sets = F_SETS, de_results = F_DE_RDS,
            master_de = F_DE_CSV, ortholog_map = F_OMAP, projection_manifest = F_MANIFEST)
missing_inputs <- inputs[!file.exists(inputs)]
if (length(missing_inputs))
  stop("[22] required input(s) absent — REPORT this rather than improvising a substitute:\n",
       paste(sprintf("    %s: %s", names(missing_inputs), missing_inputs), collapse = "\n"))
message("  [read] all 6 inputs present.")

# ============================================================================
# 2. EXPRESSION + DESIGN. The four design cells are DERIVED from metadata (genotype x
#    temperature), read from the metadata. Group display order is
#    genotype-major, temperature ascending, with the genotype reference level first
#    (the order metadata itself presents), so a dot plot reads WT 37 -> WT 39 ->
#    KO 37 -> KO 39 and the heat step within each genotype is adjacent.
# ============================================================================

eda <- readRDS(F_EDA)
if (is.null(eda[[EXPR_SLOT]]))
  stop("[22] 01_eda.rds has no '", EXPR_SLOT, "' slot (found: ",
       paste(names(eda), collapse = ", "), ").")
expr <- as.matrix(eda[[EXPR_SLOT]])
meta <- as.data.frame(eda$metadata)
if (is.null(meta) || !all(c("genotype", "temp", "group") %in% colnames(meta)))
  stop("[22] 01_eda.rds$metadata must carry genotype/temp/group columns (found: ",
       paste(colnames(meta), collapse = ", "), ").")

# sample alignment: metadata rows must be the expression columns, in the same order.
if (!identical(rownames(meta), colnames(expr)))
  stop("[22] metadata rownames do not match ", EXPR_SLOT, " colnames — sample alignment is ",
       "load-bearing for every group mean. Inspect 01_eda.rds before proceeding.")
if (anyNA(expr))
  stop("[22] ", EXPR_SLOT, " contains NA — group means would silently change meaning. ",
       "01_eda.rds records na_policy='", eda$na_policy %||% "?", "'; inspect it.")

# as.character() FIRST: metadata$temp is a factor (levels "37"/"39" set in
# 00_setup_metadata.R), and as.numeric() on a factor returns level INDICES (1, 2) —
# which would silently build the cell keys "WT_1"/"WT_2" and break every join.
geno_levels <- unique(as.character(meta$genotype))              # reference level first
temp_chr    <- as.character(meta$temp)
temp_levels <- as.character(sort(unique(suppressWarnings(as.numeric(temp_chr)))))
GROUP_LEVELS <- as.vector(t(outer(geno_levels, temp_levels,
                                  function(g, t) paste0(g, "_", t))))
meta$group <- as.character(meta$group)
if (!setequal(GROUP_LEVELS, unique(meta$group)))
  stop("[22] derived design cells {", paste(GROUP_LEVELS, collapse = ", "),
       "} do not match metadata$group {", paste(unique(meta$group), collapse = ", "),
       "} — the genotype_temp key convention changed; fix the derivation, do not hardcode.")
if (length(CFG_GROUPS) && !setequal(CFG_GROUPS, GROUP_LEVELS))
  stop("[22] design.groups in config {", paste(CFG_GROUPS, collapse = ", "),
       "} disagrees with the groups derived from metadata {",
       paste(GROUP_LEVELS, collapse = ", "), "}.")

group_n <- table(factor(meta$group, levels = GROUP_LEVELS))
if (any(group_n != N_PER_GROUP))
  stop("[22] design is not balanced at design.samples_per_group=", N_PER_GROUP, ": ",
       paste(sprintf("%s=%d", names(group_n), as.integer(group_n)), collapse = ", "))
# One (genotype, temperature) pair per design cell — the labels the tables carry.
group_key <- data.frame(
  group       = GROUP_LEVELS,
  genotype    = sub("_[^_]+$", "", GROUP_LEVELS),
  temperature = as.numeric(sub("^.*_", "", GROUP_LEVELS)),
  n_samples   = as.integer(group_n[GROUP_LEVELS]),
  stringsAsFactors = FALSE)

message(sprintf("  [design] %d samples, %d cells x n=%d: %s",
                ncol(expr), length(GROUP_LEVELS), N_PER_GROUP,
                paste(GROUP_LEVELS, collapse = ", ")))

# Group-wise mean / sd over ALL genes once (cheap at 19,679 x 20), then subset per
# signature. sd is the SAMPLE sd (n-1) across the 5 libraries in the cell.
grp_idx <- lapply(setNames(GROUP_LEVELS, GROUP_LEVELS), function(g) which(meta$group == g))
mean_mat <- vapply(grp_idx, function(ix) rowMeans(expr[, ix, drop = FALSE]),
                   numeric(nrow(expr)))
sd_mat   <- vapply(grp_idx, function(ix) apply(expr[, ix, drop = FALSE], 1L, stats::sd),
                   numeric(nrow(expr)))
dimnames(mean_mat) <- list(rownames(expr), GROUP_LEVELS)
dimnames(sd_mat)   <- list(rownames(expr), GROUP_LEVELS)

# z ACROSS THE FOUR GROUP MEANS (not across samples): each gene's four cell means are
# centred and scaled by their own sample sd, so colour is comparable between a
# high-expressed and a low-expressed gene. A k-cell z is arithmetically bounded at
# +/-(k-1)/sqrt(k) — with k=4 that is +/-1.5, the value the figure uses as its scale
# limit. A gene flat across all four cells (sd ~ 0) is reported as z = 0, not NA, so a
# flat gene reads as "no cell stands out", distinct from missing data.
z_across <- function(m) {
  mu <- rowMeans(m)
  s  <- apply(m, 1L, stats::sd)
  s[!is.finite(s) | s < 1e-12] <- NA_real_
  z <- (m - mu) / s
  z[is.na(z)] <- 0
  z
}
z_mat <- z_across(mean_mat)
Z_BOUND <- (length(GROUP_LEVELS) - 1) / sqrt(length(GROUP_LEVELS))

# ============================================================================
# 3. DE STATISTICS. 02_de_results.rds is the primary limma checkpoint; the published
#    master_de_genes.csv is derived from it. Both are read and CROSS-CHECKED on the
#    signature genes: a disagreement means one of the two is stale, which is a finding
#    to report — never something to reconcile by picking whichever number is convenient.
# ============================================================================

de_list <- readRDS(F_DE_RDS)
de_csv  <- readr::read_csv(F_DE_CSV, show_col_types = FALSE, progress = FALSE)
need_csv <- c("gene_symbol", "logFC", "t", "adj.P.Val", "contrast")
miss_csv <- setdiff(need_csv, colnames(de_csv))
if (length(miss_csv))
  stop("[22] master_de_genes.csv missing columns: ", paste(miss_csv, collapse = ", "))

de_of <- function(co) {
  d <- de_list[[co]]
  if (is.null(d)) stop("[22] 02_de_results.rds has no contrast '", co, "' (have: ",
                       paste(names(de_list), collapse = ", "), ").")
  need <- c("gene_symbol", "logFC", RANK_MET, "adj.P.Val")
  miss <- setdiff(need, colnames(d))
  if (length(miss)) stop("[22] 02_de_results.rds$", co, " missing columns: ",
                         paste(miss, collapse = ", "))
  data.frame(gene_symbol = as.character(d$gene_symbol),
             logFC = as.numeric(d$logFC),
             t = as.numeric(d[[RANK_MET]]),
             adj_P_Val = as.numeric(d$adj.P.Val),
             stringsAsFactors = FALSE)
}

# ============================================================================
# 4. SIGNATURE MEMBERSHIP. Read straight from the frozen mouse-symbol sets — this
#    script does NOT re-derive a gate; the gate arithmetic belongs to stage 10.
# ============================================================================

sig_sets <- readRDS(F_SETS)
manifest <- readr::read_csv(F_MANIFEST, show_col_types = FALSE, progress = FALSE)

genes_of <- function(s) {
  node <- sig_sets$sets[[s$source_contrast]]
  if (is.null(node))
    stop("[22] 17_signature_sets.rds has no contrast '", s$source_contrast, "'.")
  g <- node[[s$direction]][[s$gate]]
  if (is.null(g))
    stop("[22] 17_signature_sets.rds$sets$", s$source_contrast, "$", s$direction,
         " has no gate '", s$gate, "' (gates: ",
         paste(sig_sets$gates, collapse = ", "), ").")
  sort(unique(as.character(g)))
}

# Frozen-ledger cross-check on set size: manifest.csv n_mouse for the same
# (contrast, direction, gate) MUST equal what we read from the checkpoint.
manifest_row <- function(s) {
  r <- manifest[manifest$contrast == s$source_contrast &
                manifest$direction == s$direction &
                manifest$gate == s$gate, , drop = FALSE]
  if (nrow(r) != 1L)
    stop("[22] projection manifest has ", nrow(r), " rows for (", s$source_contrast, ", ",
         s$direction, ", ", s$gate, ") — expected exactly 1.")
  r[1L, , drop = FALSE]
}

# ============================================================================
# 5. ORTHOLOG FATE. The APPLIED frozen map is read from disk (no new babelgene call),
#    so this stage cannot drift from the published contract. Per MOUSE gene:
#      mapped_to_human : has >= 1 human ortholog edge
#      human_symbol    : the ortholog(s); a one2many gene keeps ALL of them, ';'-joined
#      mapping_type    : the gene's edge class, as classified over the whole applied map
#                        (one2many > many2one > one2one), so many2one reflects a human
#                        symbol shared with another mouse gene ANYWHERE in the map.
# ============================================================================

omap <- readr::read_tsv(F_OMAP, show_col_types = FALSE, progress = FALSE)
need_om <- c("mouse_symbol", "human_symbol", "mapping_type")
miss_om <- setdiff(need_om, colnames(omap))
if (length(miss_om))
  stop("[22] ortholog_map.tsv missing columns: ", paste(miss_om, collapse = ", "))

ortholog_fate <- omap %>%
  dplyr::group_by(gene_symbol = .data$mouse_symbol) %>%
  dplyr::summarise(
    human_symbol = paste(sort(unique(.data$human_symbol)), collapse = ";"),
    mapping_type = if (any(.data$mapping_type == "one2many")) "one2many"
                   else if (any(.data$mapping_type == "many2one")) "many2one"
                   else "one2one",
    n_human_symbols = dplyr::n_distinct(.data$human_symbol),
    .groups = "drop")

# ============================================================================
# 6. ASSEMBLE the three per-gene tables plus the summary, signature by signature.
# ============================================================================

expr_rows  <- list()
stat_rows  <- list()
summ_rows  <- list()

for (s in SIGS) {
  genes <- genes_of(s)
  mrow  <- manifest_row(s)
  if (!identical(as.integer(mrow$n_mouse), length(genes)))
    stop("[22] MISMATCH vs the frozen projection ledger for ", s$name, ": ",
         "17_signature_sets.rds gives ", length(genes), " mouse genes, manifest.csv n_mouse = ",
         mrow$n_mouse, ". Report this — do not adjust either side to agree.")

  absent <- setdiff(genes, rownames(expr))
  if (length(absent))
    stop("[22] ", length(absent), " gene(s) of ", s$name, " are absent from ", EXPR_SLOT,
         " (e.g. ", paste(utils::head(absent, 5), collapse = ", "),
         ") — the signature and the expression matrix are on different gene universes.")

  # --- per-gene DE statistics (from the primary checkpoint) + master cross-check ------
  de_all  <- de_of(s$source_contrast)
  if (anyDuplicated(de_all$gene_symbol))
    stop("[22] 02_de_results.rds$", s$source_contrast, " has duplicated gene_symbol — the ",
         "rank would be ambiguous. Inspect the DE arm.")
  de_sig  <- de_all[match(genes, de_all$gene_symbol), , drop = FALSE]
  if (anyNA(de_sig$t))
    stop("[22] ", sum(is.na(de_sig$t)), " gene(s) of ", s$name,
         " have no ", RANK_MET, " in contrast ", s$source_contrast, ".")

  chk <- de_csv[de_csv$contrast == s$source_contrast, , drop = FALSE]
  chk <- chk[match(genes, chk$gene_symbol), , drop = FALSE]
  dmax <- max(abs(chk$t - de_sig$t), abs(chk$logFC - de_sig$logFC),
              abs(chk$adj.P.Val - de_sig$adj_P_Val), na.rm = TRUE)
  if (!is.finite(dmax) || dmax > 1e-6)
    stop("[22] 02_de_results.rds and master_de_genes.csv DISAGREE on ", s$name,
         " (max abs difference in logFC/t/adj.P.Val = ", format(dmax),
         "). One of the two is stale — report this, do not pick the convenient one.")

  # Gate re-assertion, a read over the frozen sets: every member must satisfy the gate its
  # own name advertises. This catches a stale checkpoint.
  if (!all(de_sig$adj_P_Val < DE_FDR_) || !all(de_sig$logFC > 0))
    stop("[22] a member of ", s$name, " violates its advertised gate (adj.P.Val < ",
         DE_FDR_, " and logFC > 0). The checkpoint is stale relative to the DE master.")
  if (identical(s$gate, "fdr_logfc") && !all(de_sig$logFC >= DE_LOGFC_))
    stop("[22] a member of ", s$name, " has |logFC| < thresholds.de_logfc=", DE_LOGFC_,
         " despite the fdr_logfc gate. The checkpoint is stale.")

  # rank 1 = most extreme positive t of the source contrast; ties broken by symbol so
  # the rank (and therefore the plotted top-N) is deterministic across machines.
  ord  <- order(-de_sig$t, de_sig$gene_symbol)
  stat <- data.frame(
    signature       = s$name,
    gene_symbol     = de_sig$gene_symbol[ord],
    source_contrast = s$source_contrast,
    gate            = s$gate,
    logFC           = de_sig$logFC[ord],
    t               = de_sig$t[ord],
    adj_P_Val       = de_sig$adj_P_Val[ord],
    rank_within_signature = seq_along(ord),
    stringsAsFactors = FALSE) %>%
    dplyr::left_join(ortholog_fate, by = "gene_symbol") %>%
    dplyr::mutate(
      mapped_to_human = !is.na(.data$human_symbol),
      human_symbol    = ifelse(is.na(.data$human_symbol), NA_character_, .data$human_symbol),
      mapping_type    = ifelse(is.na(.data$mapping_type), NA_character_, .data$mapping_type)) %>%
    dplyr::select("signature", "gene_symbol", "source_contrast", "gate",
                  "logFC", "t", "adj_P_Val", "rank_within_signature",
                  "human_symbol", "mapping_type", "mapped_to_human")
  stat_rows[[s$name]] <- stat

  # --- per-gene x design-cell expression (long) --------------------------------------
  gm <- mean_mat[genes, , drop = FALSE]
  gs <- sd_mat[genes, , drop = FALSE]
  gz <- z_mat[genes, , drop = FALSE]
  long <- data.frame(
    signature   = s$name,
    gene_symbol = rep(genes, times = length(GROUP_LEVELS)),
    group       = rep(GROUP_LEVELS, each = length(genes)),
    mean_logcpm = as.vector(gm),
    sd_logcpm   = as.vector(gs),
    z_across_groups = as.vector(gz),
    stringsAsFactors = FALSE) %>%
    dplyr::left_join(group_key, by = "group") %>%
    dplyr::mutate(se_logcpm = .data$sd_logcpm / sqrt(.data$n_samples)) %>%
    dplyr::arrange(match(.data$gene_symbol, stat$gene_symbol),
                   match(.data$group, GROUP_LEVELS)) %>%
    dplyr::select("signature", "gene_symbol", "group", "temperature", "genotype",
                  "n_samples", "mean_logcpm", "sd_logcpm", "se_logcpm", "z_across_groups")
  expr_rows[[s$name]] <- long

  # --- one summary row per signature -------------------------------------------------
  # TWO senses of "mapped" are reported side by side because they differ and both are
  # load-bearing downstream:
  #   n_mapped_to_human       = MOUSE genes with >= 1 human ortholog. Equals
  #                             n_genes_mouse - n_dropped_no_ortholog, and equals
  #                             n_one2one + n_one2many + n_many2one. It is the sum of the
  #                             per-gene `mapped_to_human` flag in signature_gene_stats.csv.
  #   n_distinct_human_symbols = distinct HUMAN symbols the set maps onto. This is the
  #                             frozen ledger's `n_human`, and it is SMALLER whenever two
  #                             mouse genes share one human ortholog (many2one collapse).
  # For WT_heat_up these are 201 and 199: the two-gene gap is the collapse.
  n_mouse   <- length(genes)
  fate      <- stat[, c("gene_symbol", "mapped_to_human", "mapping_type", "human_symbol")]
  n_mapped  <- sum(fate$mapped_to_human)
  human_all <- unique(unlist(strsplit(stat$human_symbol[stat$mapped_to_human], ";", fixed = TRUE)))
  summ <- data.frame(
    signature   = s$name,
    gate        = s$gate,
    source_contrast = s$source_contrast,
    n_genes_mouse = n_mouse,
    n_mapped_to_human = n_mapped,
    n_distinct_human_symbols = length(human_all),
    # Every gene of the arm with no edge in the map, which is what this stage can recompute.
    # The frozen ledger SPLITS that total in two, and the split is not recomputable here:
    # babelgene keys on current MGI symbols while this matrix carries GENCODE vM25's older
    # vintage, so a symbol babelgene could not key at all used to be counted as having no
    # ortholog. Those two are carried through from the manifest, and
    # the cross-check below is against the TOTAL, which is the quantity both sides own.
    n_dropped_unmapped_total = n_mouse - n_mapped,
    n_dropped_no_ortholog = as.integer(mrow$n_dropped_no_ortholog),
    n_dropped_stale_query_symbol = as.integer(mrow$n_dropped_stale_query_symbol),
    n_query_symbol_normalised = as.integer(mrow$n_query_symbol_normalised),
    n_one2one   = sum(fate$mapping_type == "one2one",  na.rm = TRUE),
    n_one2many  = sum(fate$mapping_type == "one2many", na.rm = TRUE),
    n_many2one  = sum(fate$mapping_type == "many2one", na.rm = TRUE),
    median_abs_t = stats::median(abs(stat$t)),
    n_genes_plotted = min(TOP_N, n_mouse),
    stringsAsFactors = FALSE)

  # Frozen-ledger cross-check on ortholog fate. manifest.csv columns:
  #   n_human                  -> distinct human symbols   (our n_distinct_human_symbols)
  #   n_dropped_unmapped_total -> mouse genes with no edge (our n_dropped_unmapped_total)
  #   n_many_mapped            -> mouse genes mapping to MANY human (our n_one2many)
  # Plus the ledger's own closure, so a manifest whose split does not add up to its total is
  # caught here, upstream of this stage's summary.
  if (!identical(as.integer(mrow$n_dropped_no_ortholog) +
                   as.integer(mrow$n_dropped_stale_query_symbol),
                 as.integer(mrow$n_dropped_unmapped_total)))
    stop("[22] the frozen ledger's drop split does not close for ", s$name,
         ": no_ortholog + stale_query_symbol != unmapped_total. Report this.")
  led <- c(n_human = as.integer(mrow$n_human),
           n_dropped = as.integer(mrow$n_dropped_unmapped_total),
           n_many = as.integer(mrow$n_many_mapped))
  ours <- c(n_human = as.integer(summ$n_distinct_human_symbols),
            n_dropped = as.integer(summ$n_dropped_unmapped_total),
            n_many = as.integer(summ$n_one2many))
  if (!identical(led, ours))
    stop("[22] MISMATCH vs the frozen projection ledger for ", s$name, ":\n",
         "    manifest.csv : ", paste(sprintf("%s=%d", names(led), led), collapse = ", "), "\n",
         "    recomputed   : ", paste(sprintf("%s=%d", names(ours), ours), collapse = ", "), "\n",
         "  Report this mismatch — do NOT adjust the code to force agreement.")
  summ_rows[[s$name]] <- summ

  message(sprintf("  [%s] %s/%s: %d mouse genes | %d mapped -> %d human symbols (%d dropped: %d no ortholog + %d symbol not keyable; %d arrived only after normalisation) | median|t|=%.2f",
                  s$name, s$source_contrast, s$gate, n_mouse, n_mapped,
                  summ$n_distinct_human_symbols, summ$n_dropped_unmapped_total,
                  summ$n_dropped_no_ortholog, summ$n_dropped_stale_query_symbol,
                  summ$n_query_symbol_normalised, summ$median_abs_t))
}

SIG_LEVELS <- vapply(SIGS, `[[`, character(1), "name")
bind_ord <- function(l) dplyr::bind_rows(l) %>%
  dplyr::arrange(match(.data$signature, SIG_LEVELS))

signature_gene_group_expression <- bind_ord(expr_rows)
signature_gene_stats            <- bind_ord(stat_rows)
signature_expression_summary    <- bind_ord(summ_rows)

# ============================================================================
# 7. THE FIGURE'S SOURCE TABLE — exactly the rows 23_signature_expression_viz.R draws:
#    the top `top_n_genes` genes per signature by descending source-contrast t, crossed
#    with the four design cells, self-contained (no join needed at viz time). Where a
#    signature has fewer than top_n_genes members it contributes all of them, and the
#    per-signature member total travels in the table so the caption can say so.
# ============================================================================

plot_keys <- signature_gene_stats %>%
  dplyr::filter(.data$rank_within_signature <= TOP_N) %>%
  dplyr::select("signature", "gene_symbol", "source_contrast", "gate", "logFC", "t",
                "adj_P_Val", "rank_within_signature", "human_symbol", "mapped_to_human")

signature_dotplot <- plot_keys %>%
  dplyr::inner_join(signature_gene_group_expression, by = c("signature", "gene_symbol")) %>%
  dplyr::left_join(signature_expression_summary %>%
                     dplyr::select("signature", n_genes_in_signature = "n_genes_mouse"),
                   by = "signature") %>%
  dplyr::mutate(z_bound = Z_BOUND) %>%
  dplyr::arrange(match(.data$signature, SIG_LEVELS), .data$rank_within_signature,
                 match(.data$group, GROUP_LEVELS)) %>%
  dplyr::select("signature", "gate", "source_contrast", "gene_symbol",
                "rank_within_signature", "t", "logFC", "adj_P_Val",
                "group", "genotype", "temperature", "n_samples",
                "mean_logcpm", "sd_logcpm", "se_logcpm", "z_across_groups",
                "human_symbol", "mapped_to_human", "n_genes_in_signature", "z_bound")

# ============================================================================
# 8. WRITE (byte-stable: 9 significant figures, no row names).
# ============================================================================

message("[write] four tables under ", file.path("03_results", STAGE, "tables"), " ...")
readr::write_csv(round_numeric_cols(signature_gene_group_expression),
                 file.path(tbl_dir, "signature_gene_group_expression.csv"))
readr::write_csv(round_numeric_cols(signature_gene_stats),
                 file.path(tbl_dir, "signature_gene_stats.csv"))
readr::write_csv(round_numeric_cols(signature_expression_summary),
                 file.path(tbl_dir, "signature_expression_summary.csv"))
readr::write_csv(round_numeric_cols(signature_dotplot),
                 file.path(ov_dir, "signature_dotplot.csv"))

# ============================================================================
# 9. ACCEPTANCE CHECKS (structural + arithmetic; no statistics) and the console read.
# ============================================================================

stopifnot(
  "signature_gene_group_expression.csv missing" =
    file.exists(file.path(tbl_dir, "signature_gene_group_expression.csv")),
  "signature_gene_stats.csv missing" =
    file.exists(file.path(tbl_dir, "signature_gene_stats.csv")),
  "signature_expression_summary.csv missing" =
    file.exists(file.path(tbl_dir, "signature_expression_summary.csv")),
  "_overview/signature_dotplot.csv missing" =
    file.exists(file.path(ov_dir, "signature_dotplot.csv")),
  "long table row count != sum(n_genes) x n_groups" =
    nrow(signature_gene_group_expression) ==
      sum(signature_expression_summary$n_genes_mouse) * length(GROUP_LEVELS),
  "dotplot row count != sum(n_plotted) x n_groups" =
    nrow(signature_dotplot) ==
      sum(signature_expression_summary$n_genes_plotted) * length(GROUP_LEVELS),
  "mapped + dropped != n_genes_mouse" =
    all(signature_expression_summary$n_mapped_to_human +
        signature_expression_summary$n_dropped_unmapped_total ==
        signature_expression_summary$n_genes_mouse),
  "the drop split does not close against its total" =
    all(signature_expression_summary$n_dropped_no_ortholog +
        signature_expression_summary$n_dropped_stale_query_symbol ==
        signature_expression_summary$n_dropped_unmapped_total),
  "one2one + one2many + many2one != n_mapped_to_human" =
    all(signature_expression_summary$n_one2one + signature_expression_summary$n_one2many +
        signature_expression_summary$n_many2one ==
        signature_expression_summary$n_mapped_to_human),
  "z across the design cells exceeds its arithmetic bound" =
    max(abs(signature_gene_group_expression$z_across_groups)) <= Z_BOUND + 1e-9,
  "no DOWN-arm row may exist in this stage" =
    !any(grepl("_down|_dn", signature_gene_stats$signature, ignore.case = TRUE)))

print(signature_expression_summary)

message("=================================================================")
message("22_signature_expression COMPLETE")
message("  Tables: ", file.path("03_results", STAGE, "tables"),
        "/{signature_gene_group_expression,signature_gene_stats,signature_expression_summary}.csv")
message("          ", file.path("03_results", STAGE, "tables", "_overview"), "/signature_dotplot.csv")
message(sprintf("  Every set size and ortholog fate reproduced the frozen projection ledger (%s).",
                "03_results/human_projection/manifest.csv"))
message("  Next: Rscript 02_analysis/scripts/23_signature_expression_viz.R")
message("=================================================================")
