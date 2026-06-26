#!/usr/bin/env Rscript
# =============================================================================
# 03g_nonidentifiability_viz.R - PHASE 3 (stage 04_tf) VIZ: the non-identifiability
#   triptych (fig3p / fig3q / fig3r).
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# Role:  VISUALIZE half of the "normalize-then-visualize" split. Reads ONLY the
#        tidy tables emitted by 03g_nonidentifiability.R and renders three
#        single-claim figures. Performs NO statistics (no run_ulm/run_mlm,
#        no p.adjust/prcomp/cor, no counting); scores/counts/membership are read
#        as-is. Runs STANDALONE after the compute script.
#
# THE HONEST MESSAGE (binding framing constraints):
#   Hif1a's high heat-MAIN "activity" is NOT identifiable. The triptych shows
#   WHY without crowning anyone:
#     fig3p - Hif1a is #9 in a crowd of co-elevated stress/IEG/NF-kB TFs (the
#             gaps are tiny, p~1e-7); Hsf1 is far down at #50. NO clean winner.
#     fig3q - 92% of Hif1a's 353 targets are shared (mean 22 other TFs each);
#             the sharers are the network's most PROMISCUOUS regulators, none
#             hypoxia-specific.
#     fig3r - the same heat-driven genes populate many TFs' regulons -- so none
#             is identifiable from this contrast alone.
#   NEVER "Hif1a the #1/master TF" as truth; NEVER crown Jun/AP-1 or any
#   co-regulator as the driver; NEVER crown HIF2a/Epas1. Data-driven captions.
#
# Inputs (03_results/04_tf/tables/):
#   - fig3p_heatmain_ranking_data.csv   (rank, tf, score, p_value, family, is_hif, is_hsf)
#   - fig3q_coregulators_data.csv       (tf, shared_targets, pct_of_hif1a_set, family)
#   - fig3r_membership_data.csv         (gene, tf, in_regulon, gene_heat_t,
#                                        tf_heatmain_score, gene_contrib, tf_family)
#
# Outputs (figure-style contract, dual variants):
#   - 03_results/04_tf/figures/_overview/fig3{p,q,r}_*.{print.pdf,screen.png}
#   - 03_results/04_tf/tables/_overview/fig3{p,q,r}_*.csv
#   - 03_results/04_tf/README.md captions (via save_overview)
#
# Dependencies: config.R (palettes); figure_style.R (theme+save); ggplot2, dplyr, patchwork
# =============================================================================

source("02_analysis/config/config.R")
source("02_analysis/helpers/figure_style.R")       # contract: project_theme + save_overview; FIG_CFG
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })

STAGE  <- "04_tf"
SCRIPT <- "02_analysis/scripts/03g_nonidentifiability_viz.R"
TBL_DIR <- stage_dir("04_tf", "tables")

PROV <- provisional_caption()
rd  <- function(f) read.csv(file.path(TBL_DIR, f), check.names = FALSE, stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# FAMILY PALETTE -- single source of truth shared across fig3p/3q/3r so a family
# reads the SAME color in all three. HIF axis keeps the deck's signature orange
# (AXIS_COLORS["HIF"]); heat-shock keeps the module-stress purple
# (HEAT_AXIS_COLORS["heatshock"]); "other" is the deck grey. The promiscuous
# generic families get distinct, non-orange hues so the reader sees Hif1a is
# surrounded by NON-hypoxia families.
# -----------------------------------------------------------------------------
FAMILY_COLORS <- c(
  "HIF axis"                   = unname(AXIS_COLORS["HIF"]),        # signature orange
  "heat-shock"                 = unname(HEAT_AXIS_COLORS["heatshock"]), # stress purple
  "NF-kB"                      = "#2166AC",  # deck IFN/NFkB blue
  "AP-1/immediate-early"       = "#D6604D",  # warm red
  "housekeeping"               = "#4D4D4D",  # dark grey
  "stress/senescence"          = "#1B7837",  # green
  "proliferation/metabolism"   = "#C51B7D",  # magenta -- kept clearly OFF the heat-shock purple
  "nuclear receptor"           = "#80CDC1",  # teal
  "cytokine/STAT"              = "#3690C0",  # mid blue
  "growth-factor/other-stress" = "#BF812D",  # tan
  "other"                      = "grey70"
)
FAMILY_LEVELS <- names(FAMILY_COLORS)
fam_factor <- function(x) factor(x, levels = FAMILY_LEVELS)

# =============================================================================
# FIG 3p -- heat-MAIN ranking: Hif1a is #9 in a crowd; Hsf1 is far down at #50.
# CLAIM: no clean winner. The top ~9 scores are separated by tiny gaps
# (all p~1e-7), Hif1a sits at #9 among co-elevated stress/IEG/NF-kB TFs, and the
# heat-shock TF Hsf1 -- the one you might expect to top a HEAT contrast -- is at
# #50. Read horizontally as a ranked dotplot colored by family.
# =============================================================================
cat("[fig3p] heat-MAIN ranking ...\n")
p <- rd("fig3p_heatmain_ranking_data.csv")

N_TOP <- 22
top <- head(p[order(p$rank), ], N_TOP)
hsf <- p[p$is_hsf == TRUE | p$tf == "Hsf1", ][1, ]

# Build a plotting frame: top-22 PLUS a single Hsf1 row appended below a visual
# gap. We give Hsf1 a y-position one slot below the top-22 block so the break is
# obvious; a dashed rule marks the discontinuity.
top$ypos  <- N_TOP - top$rank + 1            # rank1 at top
hsf$ypos  <- 0                               # below the block (gap = rows 0..1)
plot_p <- rbind(
  data.frame(tf = top$tf, score = top$score, family = top$family,
             ypos = top$ypos, is_hif = top$is_hif, rank = top$rank,
             stringsAsFactors = FALSE),
  data.frame(tf = hsf$tf, score = hsf$score, family = hsf$family,
             ypos = hsf$ypos, is_hif = FALSE, rank = hsf$rank,
             stringsAsFactors = FALSE)
)
plot_p$family <- fam_factor(plot_p$family)
# Label: rank tag for Hif1a (#9) and Hsf1 (#50); plain TF elsewhere.
plot_p$lab <- ifelse(plot_p$tf == "Hif1a", "Hif1a  #9",
               ifelse(plot_p$tf == "Hsf1", "Hsf1  ...#50", plot_p$tf))

hif_score <- plot_p$score[plot_p$tf == "Hif1a"]
xmax <- max(plot_p$score) * 1.18

fig3p <- ggplot(plot_p, aes(x = score, y = ypos, color = family)) +
  # discontinuity rule between the top-22 block (ypos>=1) and the appended Hsf1 (ypos=0)
  geom_hline(yintercept = 0.5, linetype = "22", color = "grey60", linewidth = 0.4) +
  geom_segment(aes(xend = 0, yend = ypos), linewidth = 0.6) +
  # Hif1a emphasis: a larger outlined ring under the point (darker + thicker for print contrast)
  geom_point(data = subset(plot_p, tf == "Hif1a"),
             color = "grey5", size = 6.0, shape = 21, fill = NA, stroke = 1.4) +
  geom_point(size = 3.2) +
  geom_point(data = subset(plot_p, tf == "Hif1a"), size = 3.2) +
  geom_text(aes(label = lab), hjust = 0, nudge_x = max(plot_p$score) * 0.02,
            size = 4.0, color = "grey15") +
  scale_color_manual(values = FAMILY_COLORS, drop = TRUE, name = "TF family/role") +
  scale_x_continuous(limits = c(0, xmax), expand = expansion(mult = c(0, 0.02))) +
  scale_y_continuous(breaks = NULL) +
  labs(
    title = "Fig 3p. No clean winner on heat-MAIN: Hif1a is #9 in a crowd of co-elevated TFs",
    subtitle = sprintf(paste0(
      "decoupleR-ULM (CollecTRI, .mor='mor', minsize=5) on the Temp_main t-stat vector; %d TFs scored. The top scores are\n",
      "separated by tiny gaps (Jun 6.06 down to Nfkb1 5.12), all at p~1e-7 -- Hif1a (%.2f, p=%.0e) sits at #9 among generic\n",
      "stress / immediate-early / NF-kB regulators, not atop a hypoxia-specific peak. Hsf1 -- the canonical heat-shock TF --\n",
      "is far down at #50 (3.20, appended below the dashed break). Colored by curated TF family; this ranking crowns no TF."),
      nrow(p), hif_score, p$p_value[p$tf == "Hif1a"]),
    x = "heat-MAIN (Temp_main) ULM activity score", y = NULL
  ) +
  project_theme(config = FIG_CFG)

save_overview(
  fig3p, STAGE, "fig3p_heatmain_ranking",
  table = plot_p[, c("rank", "tf", "score", "family", "is_hif")],
  finding = paste0(
    "On heat-MAIN there is no clean winner: Hif1a is #9 in a crowd of co-elevated ",
    "stress / immediate-early / NF-kB TFs separated by tiny gaps (all p~1e-7), and ",
    "the canonical heat-shock TF Hsf1 is far down at #50. This ranking crowns NO ",
    "TF -- not Hif1a, not Jun/AP-1, not Epas1/HIF2a. ", PROV),
  script = SCRIPT, fn = "save_overview",
  config_kv = "figures.top_n=22; figures.base_size=16; figures.base_size_column=9",
  input = "03_results/04_tf/tables/fig3p_heatmain_ranking_data.csv",
  how_to_read = paste0(
    "Ranked dotplot: each row is a TF, x = heat-MAIN ULM activity score, colored by ",
    "curated TF family. Top-22 shown; Hsf1 appended below the dashed discontinuity ",
    "(it is #50). The dark ring marks Hif1a (#9). Claim tier: descriptive -- the ",
    "tiny gaps mean no single TF is identifiable as the driver."),
  config = FIG_CFG)

# =============================================================================
# FIG 3q -- "Hif1a's targets belong to everyone."
# CLAIM: 92% of Hif1a's 353 targets are co-regulated (mean 22 other TFs each);
# the top sharers are the network's most promiscuous regulators (Sp1, Trp53,
# NF-kB, AP-1, Myc, ...), none hypoxia-specific. Horizontal bars by % of the set.
# =============================================================================
cat("[fig3q] co-regulators ...\n")
q <- rd("fig3q_coregulators_data.csv")

N_Q <- 15
qd <- head(q[order(-q$pct_of_hif1a_set), ], N_Q)
qd$family <- fam_factor(qd$family)
qd$tf <- factor(qd$tf, levels = rev(qd$tf))   # highest at top
qd$lab <- sprintf("%.0f%%  (%d)", qd$pct_of_hif1a_set, qd$shared_targets)

fig3q <- ggplot(qd, aes(x = pct_of_hif1a_set, y = tf, fill = family)) +
  geom_col(width = 0.72) +
  geom_text(aes(label = lab), hjust = -0.08, size = 4.0, color = "grey15") +
  scale_fill_manual(values = FAMILY_COLORS, drop = TRUE, name = "TF family/role") +
  scale_x_continuous(limits = c(0, max(qd$pct_of_hif1a_set) * 1.18),
                     expand = expansion(mult = c(0, 0.02)),
                     labels = function(x) paste0(x, "%")) +
  labs(
    title = "Fig 3q. Hif1a's targets belong to everyone",
    subtitle = paste0(
      "Of Hif1a's 353 analyzed CollecTRI targets, 325 (92.1%) are also regulated by >=1 other TF (mean 22.2 other TFs per\n",
      "target). The top sharers are the network's most promiscuous regulators -- housekeeping (Sp1), stress/senescence\n",
      "(Trp53), NF-kB, AP-1/immediate-early, proliferation (Myc) -- none hypoxia-specific. Bars = % of the 353-target set the\n",
      "TF also regulates (shared count in parentheses); this is why Hif1a's heat-MAIN signal cannot be attributed to Hif1a alone."),
    x = "share of Hif1a's 353 targets also regulated by this TF",
    y = NULL
  ) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank())

save_overview(
  fig3q, STAGE, "fig3q_coregulators",
  table = qd[, c("tf", "shared_targets", "pct_of_hif1a_set", "family")],
  finding = paste0(
    "92% of Hif1a's 353 CollecTRI targets are co-regulated (mean 22 other TFs ",
    "each); the top sharers are the network's most promiscuous regulators ",
    "(Sp1, Trp53, NF-kB, AP-1, Myc), none hypoxia-specific -- so Hif1a's heat-MAIN ",
    "signal cannot be attributed to Hif1a alone, and no co-regulator is crowned ",
    "the driver. ", PROV),
  script = SCRIPT, fn = "save_overview",
  config_kv = "figures.top_n=15; figures.base_size=16; figures.base_size_column=9",
  input = "03_results/04_tf/tables/fig3q_coregulators_data.csv",
  how_to_read = paste0(
    "Horizontal bars (top-15 sharers): x = % of Hif1a's 353 targets that this TF ",
    "also regulates (shared count in parentheses), colored by TF family. None of ",
    "the top sharers is hypoxia-specific. Claim tier: descriptive -- targets ",
    "'belong to everyone', so the signal is not Hif1a-specific."),
  config = FIG_CFG)

# =============================================================================
# FIG 3r -- THE unifier: shared ownership grid.
# CLAIM: the same heat-driven genes populate many TFs' regulons -- so no single
# TF (Hif1a included) is identifiable from this contrast. gene(row) x TF(col)
# grid of regulon membership, with a LEFT marginal bar (gene heat-MAIN t) and a
# TOP marginal bar (TF heat-MAIN ULM score) so the reader sees BOTH margins read
# high. Rows ordered by Hif1a signed contrib; columns by TF heat-MAIN score.
#
# PATCHWORK NOTE: each sub-panel carries project_theme(config = FIG_CFG) so the
# composite is already on-contract BEFORE save_figure's per-variant re-theme.
# The marginal panels then blank their redundant shared axis (structural layout,
# not a styling-floor override) so the three panels align. save_figure adds
# project_theme(variant=) to the assembled patchwork; we verify both .print and
# .screen render filled (not blank) after the round trip.
# =============================================================================
cat("[fig3r] shared-ownership grid ...\n")
m <- rd("fig3r_membership_data.csv")

# Column order = TF heat-MAIN score desc (Hif1a flagged in label).
tf_ord <- m %>%
  dplyr::distinct(tf, tf_heatmain_score) %>%
  dplyr::arrange(dplyr::desc(tf_heatmain_score))
tf_levels <- tf_ord$tf
# Row order = gene contrib desc (top contributor at TOP of the grid).
gene_ord <- m %>%
  dplyr::distinct(gene, gene_contrib, gene_heat_t) %>%
  dplyr::arrange(dplyr::desc(gene_contrib))
gene_levels <- rev(gene_ord$gene)   # rev so highest contrib sits at the top

m$tf   <- factor(m$tf, levels = tf_levels)
m$gene <- factor(m$gene, levels = gene_levels)

# Column-label styling: Hif1a bold + orange so the eye finds it; others grey.
col_face  <- ifelse(tf_levels == "Hif1a", "bold", "plain")
col_color <- ifelse(tf_levels == "Hif1a", unname(AXIS_COLORS["HIF"]), "grey25")

MEMBER_FILL <- "#01665E"   # single clean accent (deck teal) for membership

# --- center grid: tile heatmap of in_regulon -----------------------------
grid <- ggplot(m, aes(x = tf, y = gene)) +
  geom_tile(aes(fill = factor(in_regulon)), color = "white", linewidth = 0.6) +
  scale_fill_manual(values = c("0" = "grey92", "1" = MEMBER_FILL), guide = "none") +
  scale_x_discrete(position = "bottom") +
  labs(x = NULL, y = NULL) +
  project_theme(config = FIG_CFG, legend = FALSE) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, face = col_face,
                                   color = col_color),
        plot.margin = margin(2, 4, 2, 2))

# --- LEFT marginal: gene heat-MAIN t (one bar per row, aligned to grid y) --
gdf <- gene_ord
gdf$gene <- factor(gdf$gene, levels = gene_levels)
left <- ggplot(gdf, aes(x = gene_heat_t, y = gene)) +
  geom_col(fill = "grey45", width = 0.7) +
  scale_x_reverse(expand = expansion(mult = c(0.05, 0))) +
  labs(x = "gene t\n(heat-MAIN)", y = NULL) +
  project_theme(config = FIG_CFG, legend = FALSE) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_blank(),       # structural: shares grid's y
        plot.margin = margin(2, 2, 2, 2))

# --- TOP marginal: TF heat-MAIN ULM score (one bar per column, aligned x) --
tdf <- tf_ord
tdf$tf <- factor(tdf$tf, levels = tf_levels)
top_fill <- ifelse(tdf$tf == "Hif1a", unname(AXIS_COLORS["HIF"]), "grey55")
topbar <- ggplot(tdf, aes(x = tf, y = tf_heatmain_score)) +
  geom_col(fill = top_fill, width = 0.7) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(x = NULL, y = "TF score\n(heat-MAIN)") +
  project_theme(config = FIG_CFG, legend = FALSE) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.text.x = element_blank(),       # structural: shares grid's x
        plot.margin = margin(2, 4, 2, 2))

# --- assemble: blank | topbar  /  left | grid  ---------------------------
blank <- patchwork::plot_spacer()
comp <- (blank + topbar + plot_layout(widths = c(1, 4))) /
        (left  + grid   + plot_layout(widths = c(1, 4))) +
        plot_layout(heights = c(1, 5))

fig3r <- comp + plot_annotation(
  title = "Fig 3r. The same heat-driven genes populate many TFs' regulons -- why none is identifiable",
  subtitle = paste0(
    "Top-20 heat-MAIN-driving genes of Hif1a's regulon (rows, ordered by signed contribution) x Hif1a + its 12 largest target-\n",
    "sharers (columns, ordered by heat-MAIN ULM score; Hif1a in orange). Filled = the gene is in that TF's CollecTRI regulon.\n",
    "LEFT bar = each gene's heat-MAIN t (all high -> the rows are genuinely heat-driven); TOP bar = each TF's heat-MAIN score\n",
    "(all high -> every column reads as 'active'). Because the same heat-driven genes sit in many regulons, the contrast cannot\n",
    "single out Hif1a -- or any one TF -- as the driver."))

# fig3r table neighbor: the membership grid the figure draws.
fig3r_tbl <- m[, c("gene", "tf", "in_regulon", "gene_heat_t",
                   "tf_heatmain_score", "gene_contrib")]

save_overview(
  fig3r, STAGE, "fig3r_shared_ownership",
  table = fig3r_tbl,
  finding = paste0(
    "The same heat-driven genes populate many TFs' CollecTRI regulons: every gene ",
    "(LEFT bar) and every TF (TOP bar) reads high on heat-MAIN, so the contrast ",
    "cannot single out Hif1a -- or any one TF -- as the driver. None is identifiable. ",
    PROV),
  script = SCRIPT, fn = "save_overview",
  config_kv = "figures.base_size=16; figures.base_size_column=9; patchwork=TRUE",
  input = "03_results/04_tf/tables/fig3r_membership_data.csv",
  how_to_read = paste0(
    "Center grid: rows = top-20 heat-MAIN-driving genes of Hif1a's regulon (signed ",
    "contribution order), columns = Hif1a + its 12 largest target-sharers (heat-MAIN ",
    "ULM score order; Hif1a orange + bold). Teal tile = the gene is in that TF's ",
    "regulon. LEFT bar = gene heat-MAIN t; TOP bar = TF heat-MAIN score (both read ",
    "high). Claim tier: descriptive non-identifiability -- shared ownership, no ",
    "single driver."),
  width = 10, height = 7.6, config = FIG_CFG)

cat("[DONE] 03g_nonidentifiability_viz.R -- fig3p/3q/3r rendered from tidy tables (no statistics).\n")
