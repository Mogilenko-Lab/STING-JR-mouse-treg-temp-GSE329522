# 23_signature_expression_viz.R — VIZ
# =============================================================================
# Signature-gene expression dot plot (stage 14_signature_expression).
#
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   14_signature_expression
#
# ROLE: VIZ ONLY. Reads the _overview source table that 22_signature_expression.R wrote.
#   Every row plotted is a row of that CSV, which keeps the CSV and the PNG in agreement.
#   Figures go through the figure-style contract alone: project_theme(config = FIG_CFG) +
#   save_overview() (dual pdf+png + sibling source table + README caption, atomic), with
#   colors from FIG_CFG$colors.
#
# Figure (one, _overview, carrying its same-stem source table):
#   _overview/signature_dotplot   gene x design-cell dot plot, faceted by signature.
#                                 dot FILL = across-cell z of the group means (diverging)
#                                 dot SIZE = group mean log2 CPM (expression level)
#
# This is the wet-lab readback of the exported UP arms: which member genes are abundant, in
# which of the four design cells each one peaks, and whether that peak moves when cGAS is
# removed. It is the per-gene view of sets that stages 10/11 report as counts. UP ARMS ONLY
# — the down arms sit outside this stage's scope by design (see
# analysis_config.yaml::signature_expression and the compute sibling's header).
#
# Inputs (read-only):
#   03_results/14_signature_expression/tables/_overview/signature_dotplot.csv
#
# Run from project root:
#   Rscript 02_analysis/scripts/23_signature_expression_viz.R
# =============================================================================

source("02_analysis/helpers/figure_style.R")   # project_theme, save_overview, purge_figures, FIG_CFG
source("02_analysis/config/config.R")          # sample_mapping_stamp()/_caption() — the ONE sample-provenance source

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ============================================================================
# CONSTANTS (from FIG_CFG / config; never hardcoded)
# ============================================================================

STAGE  <- "14_signature_expression"
SCRIPT <- "02_analysis/scripts/23_signature_expression_viz.R"
OV_DIR <- file.path("03_results", STAGE, "tables", FIG_CFG$figures$overview_dir %||% "_overview")

SE_CFG <- FIG_CFG$signature_expression %||% list()
TOP_N  <- as.integer(SE_CFG$top_n_genes %||% 15L)
NCOL   <- as.integer(SE_CFG$facet_ncol  %||% 2L)
# Signature facet order = config roster order (WT_heat_up first — the primary axis).
SIG_LEVELS <- vapply(SE_CFG$signatures %||% list(), function(s) as.character(s$name), character(1))
DOT_RANGE  <- as.numeric(unlist(SE_CFG$dot_size_range %||% c(1.2, 7.0)))

# Diverging triplet from config (low = blue, mid = near-white, high = orange).
POS <- FIG_CFG$colors$diverging$up      %||% "#B35806"
MID <- FIG_CFG$colors$diverging$neutral %||% "#F7F7F7"
NEG <- FIG_CFG$colors$diverging$down    %||% "#2166AC"

DE_FDR   <- as.numeric(FIG_CFG$thresholds$de_fdr   %||% 0.05)
DE_LOGFC <- as.numeric(FIG_CFG$thresholds$de_logfc %||% 1.0)

# GEOMETRY. The figure is laid out 2x2, over 4-across, for a legibility reason worth
# stating: what gets normalised when someone views the PNG is the canvas WIDTH, so a wide
# canvas is a SMALL figure. Four panels side by side need ~17in; displayed at the notebook's
# ~1150px that is 68 px/in, which puts the gene labels near 9px — unreadable. Two panel
# columns fit the project DEFAULT width (8.5in, so no width override at all), giving 135
# px/in and ~19px labels. Only the height is overridden, and only because two stacked panel
# rows of top-N gene rows need the room; project_theme still enforces the font floor.
FIG_W <- as.numeric(FIG_CFG$figures$width %||% 8.5)   # config default — NOT an override
FIG_H <- 10.5                                          # the one geometry override
# Text wrap: conservative characters-per-inch per text tier (the 19_hsr_decomposition_viz
# idiom). Measured-enough, and deliberately pessimistic — an under-filled line costs
# nothing, a clipped one costs the sentence.
CHARS_PER_INCH <- c(title = 8.5, subtitle = 11, caption = 13)
fit_text <- function(..., width_in = FIG_W, tier = "caption") {
  w <- max(20L, as.integer(floor(width_in * CHARS_PER_INCH[[tier]])))
  paste(strwrap(paste0(...), width = w), collapse = "\n")
}

# ============================================================================
# 1. GUARD + READ the source table (the ONLY input).
# ============================================================================

f_dot <- file.path(OV_DIR, "signature_dotplot.csv")
if (!file.exists(f_dot))
  stop("[23_viz] missing source table: ", f_dot,
       " — run 22_signature_expression.R first.")
dot <- readr::read_csv(f_dot, show_col_types = FALSE, progress = FALSE)

need <- c("signature", "gate", "gene_symbol", "rank_within_signature", "t",
          "group", "genotype", "temperature", "mean_logcpm", "z_across_groups",
          "n_genes_in_signature", "z_bound")
miss <- setdiff(need, colnames(dot))
if (length(miss)) stop("[23_viz] signature_dotplot.csv missing columns: ",
                       paste(miss, collapse = ", "))
if (length(SIG_LEVELS) == 0L) SIG_LEVELS <- unique(dot$signature)

message("[23_viz] loaded signature_dotplot.csv (", nrow(dot), " rows, ",
        dplyr::n_distinct(dot$gene_symbol), " genes, ",
        dplyr::n_distinct(dot$signature), " signatures).")

# ============================================================================
# 2. PLOT-READY RESHAPE (labels and factor order only — no numbers are touched).
#    Facets carry different gene sets, so the y axis is free per facet; the row key is
#    therefore (signature, gene) and the label is stripped back to the gene at draw time.
#    Levels run from the WORST rank to rank 1 so rank 1 lands at the TOP of each panel.
# ============================================================================

# Design-cell tick labels are DERIVED from the table's genotype/temperature columns,
# not from a hand-written lookup, and ordered genotype-major / temperature-ascending so
# the heat step within a genotype is the adjacent pair.
cells <- dot %>%
  dplyr::distinct(.data$group, .data$genotype, .data$temperature) %>%
  dplyr::arrange(match(.data$genotype, unique(dot$genotype)), .data$temperature) %>%
  dplyr::mutate(cell_label = paste0(.data$genotype, "\n", .data$temperature, " °C"))

Z_BOUND <- max(abs(dot$z_bound), na.rm = TRUE)

# Facet strip: what the signature IS (contrast/gate) and how much of it is on screen.
# THREE lines: the longest signature name (Interaction_up_fdrOnly) plus its gate
# on one line overruns the panel and the strip clips mid-word. The third line is load-bearing
# too — it is where a capped arm admits it is capped ("top 15 of 213" vs "all 9 of 9"), so a
# truncated panel can never be mistaken for a complete set.
facet_labs <- dot %>%
  dplyr::group_by(.data$signature) %>%
  dplyr::summarise(gate = .data$gate[1], total = .data$n_genes_in_signature[1],
                   shown = dplyr::n_distinct(.data$gene_symbol), .groups = "drop") %>%
  dplyr::mutate(lab = sprintf("%s\n[gate: %s]\n%s of %d genes",
                              .data$signature, .data$gate,
                              ifelse(.data$shown < .data$total,
                                     paste0("top ", .data$shown), paste0("all ", .data$shown)),
                              .data$total))
LAB <- setNames(facet_labs$lab, facet_labs$signature)

# Row key = signature + gene, joined by a delimiter that cannot occur in a gene symbol.
# The facets carry DIFFERENT gene sets, so the y factor has to be unique per facet row;
# scale_y_discrete() strips the signature prefix back off for the visible tick label.
ROW_SEP <- "@@"

dot_p <- dot %>%
  dplyr::mutate(
    signature = factor(.data$signature, levels = SIG_LEVELS),
    group     = factor(.data$group, levels = cells$group,
                       labels = cells$cell_label),
    row_key   = factor(paste(.data$signature, .data$gene_symbol, sep = ROW_SEP))) %>%
  dplyr::arrange(.data$signature, dplyr::desc(.data$rank_within_signature))
dot_p$row_key <- factor(dot_p$row_key, levels = unique(as.character(dot_p$row_key)))

# ============================================================================
# 3. FIGURE: gene x design-cell dot plot, faceted by signature.
#    Dots are shape 21 (fill + thin stroke) on purpose: the diverging MID colour is
#    near-white, so an unfilled point at z ~ 0 would disappear against the panel. The
#    stroke keeps "expressed but flat across the four cells" visible as a pale dot
#    distinct from missing data.
# ============================================================================

fig_dot <- ggplot(dot_p, aes(x = group, y = row_key)) +
  geom_point(aes(fill = z_across_groups, size = mean_logcpm),
             shape = 21, colour = "grey35", stroke = 0.3) +
  facet_wrap(~ signature, ncol = NCOL, scales = "free_y",
             labeller = labeller(signature = LAB)) +
  scale_y_discrete(labels = function(x) sub(paste0("^.*", ROW_SEP), "", x)) +
  scale_fill_gradient2(low = NEG, mid = MID, high = POS, midpoint = 0,
                       limits = c(-Z_BOUND, Z_BOUND),
                       breaks = c(-1, 0, 1),
                       name = "z of the four cell means") +
  scale_size_continuous(range = DOT_RANGE, name = "mean log2 CPM") +
  # Legend BELOW: in a 2-column layout a right-hand legend would eat ~20% of the
  # 8.5in width and hand it back as narrower panels. Underneath, the panels keep the full width.
  guides(fill = guide_colourbar(order = 1, title.position = "top",
                                barwidth = grid::unit(7, "lines"),
                                barheight = grid::unit(0.6, "lines")),
         size = guide_legend(order = 2, title.position = "top", nrow = 1,
                             override.aes = list(fill = "grey75"))) +
  labs(title = "Signature genes across the 2x2 design (up arms)",
       subtitle = fit_text(
         "Each exported UP arm read back gene by gene in the four design cells, ordered by the ",
         "arm's own source-contrast t. Down arms are outside this stage's scope.",
         tier = "subtitle"),
       x = NULL, y = NULL,
       caption = fit_text(
         "Group mean of log2 CPM over n=5 libraries per design cell (01_eda). Fill = the gene's ",
         "four cell means z-scored across those cells (bounded at ±", signif(Z_BOUND, 3),
         " for four cells); size = mean log2 CPM. Gate: adj.P.Val < ", DE_FDR,
         " (fdr_only), or additionally |log2FC| >= ", DE_LOGFC, " (fdr_logfc). Panels show each ",
         "arm's top ", TOP_N, " genes by t, or all of them where an arm has fewer — the strip says ",
         "which, and the stage tables carry every member. Mouse symbols. Claim tier: L3 ",
         "(n=5/group). ", sample_mapping_stamp(), ".")) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 0),
        axis.text.y = element_text(size = (FIG_CFG$figures$axis_text_size %||% 11) * 0.9),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(linewidth = 0.2, colour = "grey92"),
        # The caption is laid out against the PANEL by default; "plot" recovers the full
        # canvas width the wrap above was computed against.
        plot.caption.position = "plot",
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.justification = "left")

# ============================================================================
# 4. SAVE via save_overview (figure + same-stem source table + README caption).
#    The neighbour table is the table this figure was built from, passed through
#    unchanged; save_overview re-writes it byte-stably.
# ============================================================================

purge_figures(STAGE, "signature_dotplot", overview = TRUE, config = FIG_CFG)

save_overview(
  fig_dot, STAGE, "signature_dotplot",
  table   = dot,
  # The finding is stated over EVERY member, past the plotted top-N, and is checkable from the
  # stage's own full table (signature_gene_group_expression.csv): 213/213 and 239/239 for the
  # heat arms, 9/9 and 23/23 for the two interaction gates. An earlier draft added "from a much
  # lower abundance base", which held on the top-25 plotted subset and does NOT hold over all
  # members (median peak log2 CPM 3.6/3.0 for the heat arms vs 1.8 for Interaction_up but 3.7
  # for Interaction_up_fdrOnly), so that clause is gone.
  finding = paste0("Across every member, not only the plotted ones, WT_heat_up and KO_heat_up ",
                   "genes have a 39 °C cell as their high cell and a 37 °C cell as their low cell ",
                   "in BOTH genotypes, whereas cGAS-KO 39 °C is the LOW cell for all 9 ",
                   "Interaction_up and all 23 Interaction_up_fdrOnly genes."),
  script  = SCRIPT, fn = "ggplot(geom_point, shape=21)",
  config_kv = paste0("signature_expression.top_n_genes=", TOP_N,
                     "; signature_expression.dot_size_range=[",
                     paste(DOT_RANGE, collapse = ", "), "]",
                     "; thresholds.de_fdr=", DE_FDR, "; thresholds.de_logfc=", DE_LOGFC,
                     "; colors.diverging"),
  input   = file.path(OV_DIR, "signature_dotplot.csv"),
  how_to_read = paste0(
    "One row per gene, one column per design cell (genotype x temperature, n=5). Dot FILL = ",
    "the gene's four cell means z-scored across those cells: orange = its high cell, blue = ",
    "its low cell, pale = flat; the scale runs to ±", signif(Z_BOUND, 3),
    ". Dot SIZE = mean log2 CPM. Colour and size are independent channels. ",
    "Rows run by the arm's own source-contrast t, rank 1 at top; each strip ",
    "gives the on-screen count. The interaction panels' blue column is cGAS-KO 39 °C — a ",
    "deficit in one cell, the arithmetic content of an interaction term. Interaction_up and ",
    "Interaction_up_fdrOnly are one 1-df contrast at two gates, not two signatures. Up arms ",
    "only. Claim tier: L3 (n=5/group). ", sample_mapping_caption()),
  # GEOMETRY OVERRIDE (contract-sanctioned passthrough; NOT raw ggsave) — see FIG_W/FIG_H
  # above, which the text wrap is also computed against.
  width = FIG_W, height = FIG_H,
  config = FIG_CFG)

# ============================================================================
# 5. CAPTIONS for the compute-owned tables (READ-ONLY here; produced by
#    22_signature_expression.R). save_overview captions the figure and rewrites its
#    same-stem neighbour, but the other three tables need their own entries or the
#    stage README would document only the figure. write_caption() is idempotent and
#    keyed on the filename, so the whole README regenerates from this script alone.
# ============================================================================

write_caption(
  stage    = STAGE,
  filename = "tables/_overview/signature_dotplot.csv",
  finding  = paste0("The ", nrow(dot), " rows the dot plot draws: the top ", TOP_N,
                    " members of each up arm by descending source-contrast t (or all of them ",
                    "where an arm has fewer), crossed with the four design cells."),
  script   = "02_analysis/scripts/22_signature_expression.R",
  fn       = "top-N rank cut + inner_join (rewritten byte-stably by save_overview)",
  config_kv = paste0("signature_expression.top_n_genes=", TOP_N),
  input    = "03_results/14_signature_expression/tables/signature_gene_stats.csv (+ signature_gene_group_expression.csv)",
  how_to_read = paste0(
    "Same-stem neighbour of figures/_overview/signature_dotplot.png and deliberately ",
    "self-contained: the figure joins nothing at draw time, so a row here is a dot there. ",
    "z_across_groups is the colour channel, mean_logcpm the size channel, ",
    "rank_within_signature the row order (1 = top). n_genes_in_signature is the arm's FULL ",
    "member count, so ", TOP_N, "-of-213 reads as a cap rather than a set size; z_bound carries ",
    "the four-cell z limit (", signif(Z_BOUND, 3), ") the colour scale uses. Downstream ",
    "notebooks discover this table by its same-stem pairing with the PNG — keep the pair ",
    "intact. Claim tier: L3 (n=5/group). ", sample_mapping_caption()),
  config   = FIG_CFG)

write_caption(
  stage    = STAGE,
  filename = "tables/signature_gene_group_expression.csv",
  finding  = paste0("Every member of all four up arms in each of the four design cells: ",
                    "group-mean log2 CPM with its dispersion and across-cell z."),
  script   = "02_analysis/scripts/22_signature_expression.R",
  fn       = "group-mean/sd reshape + z_across()",
  config_kv = paste0("signature_expression.expression_matrix=logcpm_mat; ",
                     "design.samples_per_group=", FIG_CFG$design$samples_per_group %||% 5),
  input    = "03_results/objects/01_eda.rds, 03_results/objects/17_signature_sets.rds",
  how_to_read = paste0(
    "Long format, one row per (signature, gene_symbol, group). group is the genotype x ",
    "temperature design cell, also split into genotype and temperature; n_samples is 5 ",
    "everywhere. mean_logcpm / sd_logcpm are the mean and SAMPLE sd over those 5 libraries, ",
    "se_logcpm is sd/sqrt(n_samples). z_across_groups z-scores a gene's FOUR CELL MEANS against ",
    "each other, not its samples, so 'which cell is this gene's high cell' reads the same for an ",
    "abundant and a scarce gene; it is bounded at ±", signif(Z_BOUND, 3), ", and a gene flat ",
    "across all four cells is reported as 0, not NA. One gene can appear under several ",
    "signatures, so signature is part of the key. Up arms only. Claim tier: L3."),
  config   = FIG_CFG)

write_caption(
  stage    = STAGE,
  filename = "tables/signature_gene_stats.csv",
  finding  = paste0("One row per (signature, gene): the DE statistics the gene was selected on, ",
                    "its rank inside its own arm, and what happens to it on the way to human."),
  script   = "02_analysis/scripts/22_signature_expression.R",
  fn       = "de_of() + ortholog_fate left-join",
  config_kv = paste0("gsea.rank_metric=t; thresholds.de_fdr=", DE_FDR,
                     "; thresholds.de_logfc=", DE_LOGFC),
  input    = paste0("03_results/objects/02_de_results.rds, 03_results/master/master_de_genes.csv, ",
                    "03_results/human_projection/ortholog_map.tsv"),
  how_to_read = paste0(
    "logFC, t and adj_P_Val are read verbatim from source_contrast's limma table; nothing is ",
    "re-fitted, and the two DE sources (02_de_results.rds, master_de_genes.csv) are cross-checked ",
    "against each other before any row is written. rank_within_signature = 1 is the arm's most ",
    "extreme positive t, ties broken on symbol so the top-N cut is deterministic. ",
    "mapped_to_human is FALSE when the gene has no ortholog in the frozen applied map, and is ",
    "therefore dropped from the human projection; human_symbol lists all orthologs, ';'-separated. ",
    "mapping_type is the gene's edge class over the WHOLE applied map, so many2one means its ",
    "ortholog is shared with another mouse gene somewhere in the map, not inside this signature. ",
    "Claim tier: L3."),
  config   = FIG_CFG)

write_caption(
  stage    = STAGE,
  filename = "tables/signature_expression_summary.csv",
  finding  = paste0("One row per up arm: its size in mouse symbols, its ortholog fate, and the ",
                    "typical effect size of its members — every value reproduces the frozen ",
                    "projection ledger."),
  script   = "02_analysis/scripts/22_signature_expression.R",
  fn       = "per-signature summary + manifest_row() cross-check",
  config_kv = paste0("signature_expression.signatures; signature_expression.top_n_genes=", TOP_N),
  input    = paste0("03_results/objects/17_signature_sets.rds, ",
                    "03_results/human_projection/{ortholog_map.tsv,manifest.csv}"),
  how_to_read = paste0(
    "TWO senses of 'mapped' sit side by side because they differ and both matter downstream. ",
    "n_mapped_to_human counts MOUSE genes with at least one human ortholog, so n_genes_mouse - ",
    "n_dropped_unmapped_total = n_mapped_to_human = n_one2one + n_one2many + n_many2one; it is ",
    "the sum of the per-gene mapped_to_human flag in signature_gene_stats.csv. ",
    "n_distinct_human_symbols counts the HUMAN symbols the arm lands on — what ",
    "03_results/human_projection/manifest.csv calls n_human — and is smaller whenever two mouse ",
    "genes collapse onto one ortholog.\n\n",
    "The drop total then splits, and the split is the point. n_dropped_no_ortholog is the ",
    "orthology source knowing the symbol and having no human counterpart. ",
    "n_dropped_stale_query_symbol is it not being able to key the symbol at all, because this ",
    "matrix carries an older MGI vintage than the source expects — a vocabulary result that used ",
    "to be reported as the first kind. n_query_symbol_normalised counts genes that arrived ONLY ",
    "because their symbol was lifted to its current form before the source was asked.\n\n",
    "median_abs_t is 10.8/10.3 for the heat arms against 6.7/5.8 for the interaction arm, the ",
    "power gap between a 10-vs-10 contrast and a 1-df interaction term. n_genes_plotted records ",
    "the figure's cap. Tier: L3."),
  config   = FIG_CFG)

# ============================================================================
# 6. FINAL SUMMARY
# ============================================================================

n_fig <- length(list.files(file.path("03_results", STAGE, "figures"),
                           pattern = "\\.(pdf|png)$", recursive = TRUE))
message(sprintf("[23_viz] COMPLETE: %d figure file(s) under %s/figures/.", n_fig, STAGE))
if (n_fig == 0) warning("[23_viz] No figures produced — check errors above.")
