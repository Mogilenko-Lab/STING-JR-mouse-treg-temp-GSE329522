# 18_projection_export_viz.R — VIZ
# =============================================================================
# Human projection export figures (stage 11_projection).
#
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
# Stage:   11_projection
#
# ROLE: VIZ ONLY. Reads the stage tables + manifest that 18_projection_export.R
#   wrote; recomputes NOTHING. Figures ONLY via the figure-style contract
#   (project_theme(config=FIG_CFG) + save_overview()). No inline ggsave()/theme()/
#   raw hex; colors from FIG_CFG$colors.
#
# GATE-AWARE BY CONSTRUCTION. The contract exports each contrast at one or more
#   THRESHOLD GATES (the primary gate for every contrast, plus a looser secondary
#   gate for the thin comparators named in decisions.projection.secondary_gate_contrasts).
#   The plotted unit is therefore the (contrast, gate) pair, ordered as manifest.csv
#   orders it, so a contrast exported at two gates gets two bars rather than one of
#   its two rows being silently dropped. No gate NAME is written into this script:
#   the gate levels come from the tables/manifest, and the primary/secondary ROLE
#   wording is read from decisions.projection.{gate,secondary_gate}.
#
# Figures (each _overview, each carrying its same-stem source table):
#   _overview/human_signature_sizes  up/down set size per exported (contrast, gate) AFTER
#                                     mapping (human-space counterpart of the stage-17 size bars)
#   _overview/mapping_loss           per (contrast, gate) and direction, mouse-set size ->
#                                     human-set size waterfall (1:1 kept / many-mapped /
#                                     dropped-no-ortholog): the sanity check that mapping
#                                     did not decimate the signature
#   _overview/projection_overlap_ledger
#                                     overlap ledger for the frozen projected human sets:
#                                     WT/KO heat-set sharing plus both Interaction gates
#   _overview/conversion_ledger      per frozen UP arm at the primary gate, the mouse gene
#                                     count decomposed into carried / lost-to-paralog-collapse /
#                                     dropped-no-ortholog. mapping_loss shows the human count
#                                     only as a '->N' annotation, so the collapse loss has no
#                                     segment there; this panel gives it one.
#
# Inputs (read-only; produced by 18_projection_export.R and 17_signature_derive.R):
#   03_results/11_projection/tables/_overview/{human_signature_sizes,mapping_loss,projection_overlap_ledger}.csv
#   03_results/human_projection/manifest.csv
#   03_results/human_projection/ortholog_map.tsv          (conversion_ledger only)
#   03_results/human_projection/signatures/<contrast>/    (conversion_ledger only)
#   03_results/objects/17_signature_sets.rds              (conversion_ledger only)
#   00_data/references/gene_sets/temp_hsr_lens/temp_hsr_mouse_lens.rds  (optional annotation)
#
# Run from project root (after 18_projection_export.R):
#   Rscript 02_analysis/scripts/18_projection_export_viz.R
# =============================================================================

source("02_analysis/helpers/figure_style.R")   # project_theme, save_overview, purge_figures, FIG_CFG, contrast_label

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
})
options(stringsAsFactors = FALSE)

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

STAGE  <- "11_projection"
SCRIPT <- "02_analysis/scripts/18_projection_export_viz.R"
OV_DIR <- file.path("03_results", STAGE, "tables", FIG_CFG$figures$overview_dir %||% "_overview")

POS <- FIG_CFG$colors$diverging$up   %||% "#B35806"   # up
NEG <- FIG_CFG$colors$diverging$down %||% "#2166AC"   # down
# The grey bucket used to be ONE segment called "dropped (no human ortholog)", and for 146
# genes of the projection background that label was wrong: babelgene keys on CURRENT MGI
# symbols and this matrix carries GENCODE vM25's vintage, so a symbol babelgene could not key
# at all was recorded as having no ortholog. Those are now two segments, because one is a
# statement about orthology and the other about vocabulary, and drawing them in one colour is
# the defect this figure exists to expose.
LOSS_COLORS <- c(
  n_mapped_1to1                = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  n_many_mapped                = FIG_CFG$colors$okabe_ito$orange       %||% "#E69F00",
  n_dropped_no_ortholog        = "grey55",
  n_dropped_stale_query_symbol = FIG_CFG$colors$okabe_ito$reddish_purple %||% "#CC79A7")
# Alarm line: flag a contrast/direction that loses > this fraction to no-ortholog drops.
ALARM_FRAC <- 0.5

# Gate ROLE wording comes from the decision block, never from a literal gate name here.
.dcn <- FIG_CFG$decisions$projection %||% list()
GATE_PRIMARY   <- as.character(.dcn$gate           %||% NA_character_)[1]
GATE_SECONDARY <- as.character(.dcn$secondary_gate %||% NA_character_)[1]

gate_role <- function(g) {
  ifelse(!is.na(GATE_PRIMARY)   & g == GATE_PRIMARY,   "primary gate",
  ifelse(!is.na(GATE_SECONDARY) & g == GATE_SECONDARY, "gate-sensitivity", "gate"))
}

# ============================================================================
# 1. GUARD + READ.
# ============================================================================

f_sizes <- file.path(OV_DIR, "human_signature_sizes.csv")
f_loss  <- file.path(OV_DIR, "mapping_loss.csv")
f_ledger <- file.path(OV_DIR, "projection_overlap_ledger.csv")
for (f in c(f_sizes, f_loss, f_ledger))
  if (!file.exists(f))
    stop("[18_viz] missing table: ", f,
         " — run 18_projection_export.R first (it only runs after the signature freeze is signed off).")

sizes_df <- readr::read_csv(f_sizes, show_col_types = FALSE, progress = FALSE)
loss_df  <- readr::read_csv(f_loss,  show_col_types = FALSE, progress = FALSE)
ledger_df <- readr::read_csv(f_ledger, show_col_types = FALSE, progress = FALSE)

# ============================================================================
# 2. PLOTTING UNIT = (contrast, gate), ordered as the manifest orders it.
#    The manifest is the contract index; the stage tables must cover every one of
#    its up/down rows. If they do not, the compute sibling and this viz have drifted
#    apart — refuse to render a figure that quietly omits an exported gate.
# ============================================================================

unit_of <- function(contrast, gate) paste(contrast, gate, sep = "@")

man_path <- file.path("03_results", "human_projection", "manifest.csv")
man <- if (file.exists(man_path))
  readr::read_csv(man_path, show_col_types = FALSE, progress = FALSE) else NULL

if (!is.null(man)) {
  expected <- man %>%
    dplyr::filter(direction %in% c("up", "down")) %>%
    dplyr::transmute(key = paste(contrast, gate, direction, sep = "@")) %>%
    dplyr::pull(key)
  for (nm in c("human_signature_sizes", "mapping_loss")) {
    df  <- if (nm == "human_signature_sizes") sizes_df else loss_df
    got <- paste(df$contrast, df$gate, df$direction, sep = "@")
    gap <- setdiff(expected, got)
    if (length(gap))
      stop("[18_viz] ", nm, ".csv is missing manifest row(s): ", paste(gap, collapse = ", "),
           " — re-run 18_projection_export.R so the stage tables cover every exported gate.")
  }
}

UNITS <- if (!is.null(man)) {
  unique(unit_of(man$contrast, man$gate))
} else {
  unique(unit_of(sizes_df$contrast, sizes_df$gate))
}
GATES <- if (!is.null(man)) unique(man$gate) else unique(sizes_df$gate)

# A contrast exported at more than one gate needs the gate spelled out on its tick;
# a single-gate contrast keeps the plain display label.
.unit_contrast <- sub("@[^@]*$", "", UNITS)
MULTI_GATE <- names(which(table(.unit_contrast) > 1L))

unit_label <- function(u) {
  co  <- sub("@[^@]*$", "", u)
  ga  <- sub("^.*@", "", u)
  lab <- contrast_label(co, short = TRUE)
  ifelse(co %in% MULTI_GATE, paste0(lab, "\n[", ga, "]"), lab)
}

as_unit <- function(df) df %>%
  dplyr::mutate(unit = factor(unit_of(contrast, gate), levels = UNITS),
                direction = factor(direction, levels = c("up", "down"))) %>%
  droplevels()

GATE_NOTE <- paste(sprintf("%s (%s)", GATES, gate_role(GATES)), collapse = " + ")

message("[18_viz] loaded human sizes (", nrow(sizes_df), ") + mapping loss (", nrow(loss_df),
        ") rows over ", length(UNITS), " exported (contrast, gate) unit(s): ", GATE_NOTE, ".")

# ============================================================================
# 3. FIGURE (a): human_signature_sizes — up/down human-set size per exported unit.
# ============================================================================

sizes_p <- as_unit(sizes_df)

fig_sizes <- ggplot(sizes_p, aes(x = unit, y = n_human, fill = direction)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_text(aes(label = n_human),
            position = position_dodge(width = 0.8), vjust = -0.35,
            size = (FIG_CFG$figures$label_size %||% 4) * 0.85, colour = "grey20") +
  scale_fill_manual(values = c(up = POS, down = NEG),
                    labels = c(up = "up (39 > 37 / numerator)", down = "down"),
                    name = "DE direction") +
  scale_x_discrete(labels = unit_label) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Human projection signature sizes (after ortholog mapping)",
       subtitle = paste0("Frozen human-symbol up/down set sizes per exported contrast; gates = ",
                        GATE_NOTE, "."),
       x = NULL, y = "human genes in set",
       caption = "Human HGNC symbols after mouse->human mapping (babelgene). Claim tier: L3 (DE statistics).") +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 0))

# ============================================================================
# 4. FIGURE (b): mapping_loss — mouse-set -> human-set waterfall per unit/direction.
# ============================================================================

loss_long <- as_unit(loss_df) %>%
  tidyr::pivot_longer(cols = c(n_mapped_1to1, n_many_mapped, n_dropped_no_ortholog,
                               n_dropped_stale_query_symbol),
                      names_to = "category", values_to = "n") %>%
  dplyr::mutate(category = factor(category, levels = names(LOSS_COLORS)))

# human-out annotation + alarm flag per unit/direction. The alarm still fires on the TOTAL
# drop, because a set losing half its genes is worth flagging whichever bucket took them.
loss_ann <- as_unit(loss_df) %>%
  dplyr::mutate(frac_unmapped = ifelse(n_mouse > 0, n_unmapped / n_mouse, 0),
                alarm = frac_unmapped > ALARM_FRAC)

fig_loss <- ggplot(loss_long, aes(x = unit, y = n, fill = category)) +
  geom_col(width = 0.72) +
  geom_text(data = loss_ann,
            aes(x = unit, y = n_mouse, fill = NULL,
                label = sprintf("->%d", n_human)),
            vjust = -0.35, size = (FIG_CFG$figures$label_size %||% 4) * 0.8, colour = "grey15") +
  geom_point(data = dplyr::filter(loss_ann, alarm),
             aes(x = unit, y = n_mouse, fill = NULL),
             shape = 8, size = 2.4, colour = "red3", show.legend = FALSE) +
  facet_wrap(~ direction, nrow = 1,
             labeller = labeller(direction = c(up = "UP set", down = "DOWN set"))) +
  scale_fill_manual(values = LOSS_COLORS,
                    labels = c(n_mapped_1to1 = "mapped 1:1",
                               n_many_mapped = "one mouse -> many human",
                               n_dropped_no_ortholog = "dropped: keyed, no human ortholog",
                               n_dropped_stale_query_symbol = "dropped: symbol not keyable"),
                    name = "mouse gene fate") +
  scale_x_discrete(labels = unit_label) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Ortholog mapping loss: mouse set -> human set",
       subtitle = paste0("Each bar = the mouse up/down set decomposed by fate; '->N' = resulting human ",
                        "set size.\nGates = ", GATE_NOTE, ". Red * = > ",
                        round(100 * ALARM_FRAC), "% dropped (alarm)."),
       x = NULL, y = "mouse genes in set",
       caption = paste0("Bar height = mouse-set size; green kept 1:1, orange one->many (each human gets ",
                        "the mouse t). The two drop segments are different statements: grey means the ",
                        "orthology source knew the symbol and had no human counterpart, purple means it ",
                        "could not key the symbol at all because this matrix carries an older MGI ",
                        "vintage — a vocabulary result, not a biological one. Claim tier: L3. ",
                        "Correlative input, not a causal claim.")) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 0))

# ============================================================================
# 5. SAVE via save_overview (figure + same-stem source table + README caption).
# ============================================================================

purge_figures(STAGE, "human_signature_sizes", overview = TRUE, config = FIG_CFG)
purge_figures(STAGE, "mapping_loss",          overview = TRUE, config = FIG_CFG)

save_overview(
  fig_sizes, STAGE, "human_signature_sizes",
  table   = sizes_df,
  finding = "Frozen human up/down set sizes per exported contrast and threshold gate, after mouse->human ortholog mapping.",
  script  = SCRIPT, fn = "ggplot(geom_col)",
  config_kv = "decisions.projection.gate; decisions.projection.secondary_gate; colors.diverging",
  input   = f_sizes,
  how_to_read = paste0("Grouped bars per exported contrast; orange = up (higher in numerator / 39 °C), blue ",
                       "= down; numbers = human genes in the frozen set. A contrast carried at two threshold ",
                       "gates gets one bar pair per gate, with the gate named in square brackets under its ",
                       "tick — the looser gate is the sensitivity read, not a second result. This is the ",
                       "human-space counterpart of the mouse-side size bars. Claim tier: L3 (DE statistics), ",
                       "n=5/group. ", sample_mapping_caption()),
  width = 11, height = 6.5,
  config = FIG_CFG)

save_overview(
  fig_loss, STAGE, "mapping_loss",
  table   = loss_df,
  finding = "Mapping did not decimate the primary signature: most mouse genes map 1:1; dropped/many-mapped fractions are small (alarm-flagged if > half lost).",
  script  = SCRIPT, fn = "ggplot(geom_col)",
  config_kv = paste0("alarm_frac=", ALARM_FRAC, "; babelgene(offline); decisions.projection.ortholog_ambiguity; ",
                     "decisions.projection.secondary_gate"),
  input   = f_loss,
  how_to_read = paste0("Per exported contrast (x) and direction (facet), the mouse up/down set is stacked by ",
                       "fate. Green = mapped 1:1. Orange = one mouse to many human, each ortholog inheriting ",
                       "the mouse t. Grey = the orthology source keyed the symbol and had no human ",
                       "counterpart. Purple = it could not key the symbol at all.\n\n",
                       "Grey and purple were one segment until now, which reported a vocabulary loss as a ",
                       "biological one. A contrast carried at two gates gets one bar per gate, the gate named ",
                       "in brackets under its tick. Bar height = mouse-set size; '->N' = human-set size. A ",
                       "red * marks > ", round(100 * ALARM_FRAC), "% dropped, a decimated signature to ",
                       "demote. Claim tier: L3."),
  # Two direction facets x every exported (contrast, gate) unit: the canvas widens with
  # the unit count so the full contrast+gate ticks stay legible and untruncated.
  width = max(11, 3.0 + 1.6 * 2 * length(UNITS)), height = 6.5,
  config = FIG_CFG)

# ============================================================================
# 5c. FIGURE (c): projection_overlap_ledger — Panel 1D.
# ============================================================================

ledger_groups <- unique(ledger_df$display_group)
ledger_p <- ledger_df %>%
  dplyr::mutate(
    display_group = factor(display_group, levels = ledger_groups),
    direction = factor(direction, levels = c("up", "down"),
                       labels = c("up arm", "down arm")),
    component = factor(component,
                       levels = c("WT_heat only", "WT_heat ∩ KO_heat", "KO_heat only",
                                  "Interaction")),
    label = ifelse(glyph == "structural_empty", "0 structural", as.character(n_human))
  )
ledger_bars <- ledger_p %>% dplyr::filter(glyph == "bar")
ledger_empty <- ledger_p %>% dplyr::filter(glyph == "structural_empty")

ledger_cols <- c(
  "WT_heat only" = POS,
  "WT_heat ∩ KO_heat" = FIG_CFG$colors$okabe_ito$bluish_green,
  "KO_heat only" = NEG,
  "Interaction" = FIG_CFG$colors$okabe_ito$reddish_purple
)

up_shared <- ledger_df %>%
  dplyr::filter(direction == "up", component == "WT_heat ∩ KO_heat") %>%
  dplyr::slice(1)
down_shared <- ledger_df %>%
  dplyr::filter(direction == "down", component == "WT_heat ∩ KO_heat") %>%
  dplyr::slice(1)
cgas_note <- ledger_df %>%
  dplyr::filter(direction == "up", component == "WT_heat only") %>%
  dplyr::slice(1)
cgas_sentence <- if (!is.na(cgas_note$wt_only_up_cgas_dependent_n[1])) {
  sprintf("The WT-only up-arm slice has %d projected human genes; %d map from mouse genes flagged cgas_dependent in cgas_dependence_wide.csv.",
          cgas_note$wt_only_up_n[1], cgas_note$wt_only_up_cgas_dependent_n[1])
} else {
  "No cgas_dependent flag was available in the projection tables, so no WT-only cGAS-dependence claim is made."
}

# The reader arrives here having just been shown that the 39 °C response is
# largely shared between genotypes in mouse. The question this panel closes is
# whether thresholding and ortholog projection recover a separation the mouse
# geometry did not have, so the title states the answer and the subtitle names
# the question it answers.
LEDGER_W <- 12
ledger_wrap <- function(..., chars_per_inch) {
  paste(strwrap(paste0(...), width = as.integer(LEDGER_W * chars_per_inch)), collapse = "\n")
}

fig_ledger <- ggplot(ledger_bars, aes(x = display_group, y = n_human, fill = component)) +
  geom_col(width = 0.68) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5),
            size = (FIG_CFG$figures$label_size %||% 4) * 0.95, fontface = "bold",
            colour = "grey10") +
  geom_point(data = ledger_empty, aes(x = display_group, y = 0, shape = glyph),
             fill = "white", colour = "grey20", size = 3.4, stroke = 1.1,
             inherit.aes = FALSE) +
  geom_text(data = ledger_empty, aes(x = display_group, y = 13, label = label),
            size = (FIG_CFG$figures$label_size %||% 4) * 0.85, colour = "grey20",
            inherit.aes = FALSE) +
  facet_wrap(~ direction, nrow = 1) +
  scale_fill_manual(values = ledger_cols, name = "Projected set slice",
                    labels = c("WT_heat only", "shared WT_heat/KO_heat",
                               "KO_heat only", "Interaction")) +
  scale_shape_manual(values = c(structural_empty = 5), name = "Zero glyph",
                     labels = c(structural_empty = "structural empty")) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  labs(
    title = ledger_wrap(
      sprintf("The projected WT and cGAS-KO up arms share %d of %d human genes",
              up_shared$heat_shared_n[1], up_shared$heat_union_n[1]),
      chars_per_inch = 8.5),
    subtitle = ledger_wrap(
      sprintf(paste("Question — does thresholding and ortholog projection separate the heat response",
                    "by genotype? Answer — no. The two up arms overlap at Jaccard %.3f and the two",
                    "down arms at %.3f, so their enrichment scores from one ranked list are not",
                    "independent. The Interaction gates are disjoint from both."),
              up_shared$heat_jaccard[1], down_shared$heat_jaccard[1]),
      chars_per_inch = 11),
    x = NULL,
    y = "human genes",
    caption = ledger_wrap(
      "Counts are read from the frozen one-symbol-per-line signature files. Open diamond marks a ",
      "set that is structurally empty rather than measured and small. Claim tier: a direct count ",
      "over frozen sets.",
      chars_per_inch = 13)
  ) +
  project_theme(config = FIG_CFG) +
  theme(plot.caption.position = "plot")

purge_figures(STAGE, "projection_overlap_ledger", overview = TRUE, config = FIG_CFG)
save_overview(
  fig_ledger, STAGE, "projection_overlap_ledger",
  table = ledger_df,
  finding = sprintf(
    "The projected WT_heat and KO_heat up arms share %d of %d human genes, while both Interaction gates remain disjoint from the heat-set union.",
    up_shared$heat_shared_n[1], up_shared$heat_union_n[1]
  ),
  script = SCRIPT, fn = "ggplot(geom_col + structural-empty glyph)",
  config_kv = "decisions.projection.gate; decisions.projection.secondary_gate; colors.diverging; colors.okabe_ito",
  input = f_ledger,
  how_to_read = sprintf(
    paste0("Facets separate up and down arms. The left stacked bar partitions the WT_heat/KO_heat union into WT_heat-only, ",
           "shared, and KO_heat-only slices; the two right bars show Interaction at the primary and secondary gates. ",
           "Counts are printed on the bars; open diamond means structurally empty. ",
           "The up-arm heat sets share %d of %d genes (Jaccard %.3f), so their enrichment scores from the same ranked list are not independent; ",
           "seeing both move together is close to guaranteed and is not corroboration. The down-arm heat sets share %d of %d genes (Jaccard %.3f). %s ",
           "The separate HSR/TCR lens-membership question is shown by `hsr_lens_membership`."),
    up_shared$heat_shared_n[1], up_shared$heat_union_n[1], up_shared$heat_jaccard[1],
    down_shared$heat_shared_n[1], down_shared$heat_union_n[1], down_shared$heat_jaccard[1],
    cgas_sentence
  ),
  width = LEDGER_W, height = 6.8,
  config = FIG_CFG)

# ============================================================================
# 5d. FIGURE (d): conversion_ledger: where each frozen mouse UP arm's genes went.
#
#   mapping_loss (section 4) decomposes the MOUSE set by fate and reports the human
#   set size only as a "->N" annotation, so the genes lost when several mouse paralogs
#   collapse onto one human symbol are folded into that annotation. Every human-side
#   denominator turns on that number, so here it gets a segment of its own: bar length
#   is the mouse arm at the primary gate, and the three segments are carried, lost to
#   collapse, and dropped for want of an accepted ortholog.
#
#   COUNTS ARE RECOMPUTED FROM THE FROZEN CONTRACT.
#   Three independent routes to the carried count must agree before anything is
#   plotted: (1) the length of the published one-symbol-per-line signature file,
#   (2) the manifest's n_human, and (3) re-mapping the frozen mouse arm through the
#   applied ortholog map. A disagreement halts the script, because a wrong ortholog
#   ledger would misstate every human denominator downstream of it.
#
#   Reading the frozen contract and 17_signature_sets.rds from a viz script follows
#   the same pattern as 19_hsr_decomposition_viz.R and 23_signature_expression_viz.R:
#   both are published, byte-stable compute outputs, and nothing new is estimated here.
# ============================================================================

PROJ_DIR <- file.path("03_results", "human_projection")
f_omap   <- file.path(PROJ_DIR, "ortholog_map.tsv")
f_sig    <- file.path("03_results", "objects", "17_signature_sets.rds")

if (is.null(man))       stop("[18_viz] conversion_ledger needs ", man_path, ".")
if (!file.exists(f_omap)) stop("[18_viz] conversion_ledger needs ", f_omap, ".")
if (!file.exists(f_sig))  stop("[18_viz] conversion_ledger needs ", f_sig, ".")

omap     <- readr::read_tsv(f_omap, show_col_types = FALSE, progress = FALSE)
sig_sets <- readRDS(f_sig)

BABELGENE_VER <- as.character(unique(omap$babelgene_version))[1]
MIN_SUPPORT   <- (.dcn$ortholog_ambiguity %||% list())$min_support %||% NA

# Map topology: how many partners each side of an edge has. This is what defines the
# three arrival routes, and it is a property of the MAP, not of any one arm.
.named_count <- function(df, key, val) {
  d <- df %>% dplyr::group_by(.data[[key]]) %>%
    dplyr::summarise(n = dplyr::n_distinct(.data[[val]]), .groups = "drop")
  stats::setNames(d$n, d[[key]])
}
LED_N_HUMAN_PER_MOUSE <- .named_count(omap, "mouse_symbol", "human_symbol")
LED_N_MOUSE_PER_HUMAN <- .named_count(omap, "human_symbol", "mouse_symbol")
LED_MOUSE_OF_HUMAN    <- {
  d <- omap %>% dplyr::distinct(human_symbol, .keep_all = TRUE)
  stats::setNames(d$mouse_symbol, d$human_symbol)
}

# How did this human gene arrive: alone, by collapse, or by split?
arrival_route <- function(h) {
  n_mouse_side <- unname(LED_N_MOUSE_PER_HUMAN[h])
  first_mouse  <- unname(LED_MOUSE_OF_HUMAN[h])
  n_human_side <- unname(LED_N_HUMAN_PER_MOUSE[first_mouse])
  out <- rep("one-to-one", length(h))
  out[!is.na(n_human_side) & n_human_side > 1] <- "one mouse split across several human"
  out[!is.na(n_mouse_side) & n_mouse_side > 1] <- "several mouse collapsed onto it"
  out[is.na(n_mouse_side)] <- "not in map"
  out
}

read_frozen_set <- function(rel) {
  p <- file.path(PROJ_DIR, rel)
  if (!file.exists(p)) stop("[18_viz] conversion_ledger: missing frozen signature file ", p, ".")
  x <- trimws(readLines(p, warn = FALSE))
  unique(x[nzchar(x)])
}

# The arms drawn are the frozen primary contrasts, UP direction, at the PRIMARY gate.
# The manifest also carries the down arms and the secondary gate; both are out of scope
# for this panel, and the contrast list and gate come from the decision block rather
# than from literals here.
LED_ARMS <- as.character(.dcn$contrasts_primary %||% character(0))
if (!length(LED_ARMS) || is.na(GATE_PRIMARY))
  stop("[18_viz] conversion_ledger needs decisions.projection.{contrasts_primary,gate}.")

led_man <- man %>%
  dplyr::filter(contrast %in% LED_ARMS, direction == "up", gate == GATE_PRIMARY) %>%
  dplyr::mutate(.ord = match(contrast, LED_ARMS)) %>%
  dplyr::arrange(.ord) %>%
  dplyr::select(-.ord)
if (nrow(led_man) != length(LED_ARMS))
  stop("[18_viz] conversion_ledger: manifest has ", nrow(led_man), " up rows at gate ",
       GATE_PRIMARY, " for contrasts ", paste(LED_ARMS, collapse = "/"), "; expected ",
       length(LED_ARMS), ".")

led_rows <- lapply(seq_len(nrow(led_man)), function(i) {
  r         <- led_man[i, ]
  human_set <- read_frozen_set(r$file)
  routes    <- arrival_route(human_set)
  mouse_arm <- unique(as.character(sig_sets$sets[[r$contrast]]$up[[r$gate]]))
  mouse_arm <- mouse_arm[!is.na(mouse_arm) & nzchar(mouse_arm)]
  arm_edges <- omap %>% dplyr::filter(mouse_symbol %in% mouse_arm)
  # Independent re-derivation of the carried count: push the frozen mouse arm back
  # through the applied map and count distinct human symbols.
  remapped_n <- dplyr::n_distinct(arm_edges$human_symbol)
  # The mouse paralog groups that actually collided INSIDE this arm, which is where
  # the collapse loss comes from.
  cg <- arm_edges %>%
    dplyr::group_by(human_symbol) %>%
    dplyr::summarise(n_mice = dplyr::n_distinct(mouse_symbol),
                     mice   = paste(sort(unique(mouse_symbol)), collapse = "+"),
                     .groups = "drop") %>%
    dplyr::filter(n_mice > 1) %>%
    dplyr::arrange(human_symbol)
  data.frame(
    arm                        = paste0(r$contrast, "_up"),
    contrast                   = r$contrast,
    direction                  = "up",
    gate                       = r$gate,
    n_mouse_at_gate            = as.integer(r$n_mouse),
    n_human_carried            = length(human_set),
    # The collapse loss is what is left after the human genes carried and EVERY drop, so it
    # nets against the unmapped total rather than against one of the two drop classes.
    n_lost_to_collapse         = as.integer(r$n_mouse) -
                                 as.integer(r$n_dropped_unmapped_total) - length(human_set),
    n_dropped_no_ortholog      = as.integer(r$n_dropped_no_ortholog),
    n_dropped_stale_query_symbol = as.integer(r$n_dropped_stale_query_symbol),
    n_dropped_unmapped_total   = as.integer(r$n_dropped_unmapped_total),
    n_query_symbol_normalised  = as.integer(r$n_query_symbol_normalised),
    pct_of_mouse_arm_carried   = round(100 * length(human_set) / as.integer(r$n_mouse), 1),
    n_arrived_one_to_one       = sum(routes == "one-to-one"),
    n_arrived_by_collapse      = sum(routes == "several mouse collapsed onto it"),
    n_arrived_by_split         = sum(routes == "one mouse split across several human"),
    n_arrived_not_in_map       = sum(routes == "not in map"),
    collapsed_within_arm       = if (nrow(cg))
                                   paste(sprintf("%s <- %s", cg$human_symbol, cg$mice),
                                         collapse = "; ") else "",
    n_mouse_arm_object         = length(mouse_arm),
    n_human_remapped_from_arm  = remapped_n,
    manifest_n_human           = as.integer(r$n_human),
    stringsAsFactors           = FALSE)
})
conv_ledger <- dplyr::bind_rows(led_rows)

# ---- THE CROSS-CHECK. Halt on any disagreement rather than plotting either number. ----
led_bad <- conv_ledger %>%
  dplyr::filter(n_human_carried != manifest_n_human |
                n_human_carried != n_human_remapped_from_arm |
                n_mouse_at_gate != n_mouse_arm_object |
                n_lost_to_collapse < 0)
if (nrow(led_bad))
  stop("[18_viz] conversion_ledger cross-check FAILED. The frozen signature files, the ",
       "manifest, and re-mapping the frozen mouse arm through ortholog_map.tsv do not agree, ",
       "so no ledger is plotted:\n",
       paste(utils::capture.output(print(as.data.frame(
         led_bad[, c("arm", "n_mouse_at_gate", "n_mouse_arm_object", "n_human_carried",
                     "manifest_n_human", "n_human_remapped_from_arm", "n_lost_to_collapse")]))),
         collapse = "\n"))

message("[18_viz] conversion_ledger cross-check passed for ", nrow(conv_ledger),
        " up arm(s) at gate ", GATE_PRIMARY, "; map holds ", nrow(omap), " edges over ",
        dplyr::n_distinct(omap$mouse_symbol), " mouse and ",
        dplyr::n_distinct(omap$human_symbol), " human symbols.")

# Curated-lens annotation for the finding line: which genes of the curated HSR core the
# mouse arms carry, and what they become in human. Membership comes from the curated set.
# A symbol-prefix rule would be wrong here: WT_heat_up carries HSPG2, a proteoglycan.
HSR_LENS_PATH <- "00_data/references/gene_sets/temp_hsr_lens/temp_hsr_mouse_lens.rds"
led_hs_sentence <- ""
if (file.exists(HSR_LENS_PATH)) {
  hsr_core <- unique(as.character(readRDS(HSR_LENS_PATH)$HSR_core))
  hs_rows <- lapply(LED_ARMS, function(co) {
    mouse_arm <- unique(as.character(sig_sets$sets[[co]]$up[[GATE_PRIMARY]]))
    hit_m <- sort(intersect(mouse_arm, hsr_core))
    hit_h <- sort(unique(omap$human_symbol[omap$mouse_symbol %in% hit_m]))
    data.frame(contrast = co, n_m = length(hit_m), n_h = length(hit_h),
               m = paste(hit_m, collapse = ", "), h = paste(hit_h, collapse = ", "),
               stringsAsFactors = FALSE)
  })
  hs_df <- dplyr::bind_rows(hs_rows) %>% dplyr::filter(n_m > n_h)
  if (nrow(hs_df)) {
    hs_collapse <- omap %>%
      dplyr::filter(human_symbol %in% strsplit(hs_df$h[1], ", ")[[1]]) %>%
      dplyr::group_by(human_symbol) %>%
      dplyr::summarise(mice = paste(sort(unique(mouse_symbol)), collapse = " and "),
                       n_mice = dplyr::n_distinct(mouse_symbol), .groups = "drop") %>%
      dplyr::filter(n_mice > 1) %>%
      dplyr::slice(1)
    led_hs_sentence <- sprintf(
      paste0("One collapse lands on the heat-shock genes: mouse %s both map to human %s, so the ",
             "%d curated HSR core genes in the %s mouse %s (%s) arrive as %d in human (%s)."),
      hs_collapse$mice[1], hs_collapse$human_symbol[1],
      hs_df$n_m[1], paste(paste0(hs_df$contrast, "_up"), collapse = " and "),
      if (nrow(hs_df) > 1) "arms" else "arm", hs_df$m[1], hs_df$n_h[1], hs_df$h[1])
  }
}

# ---- geometry: one horizontal stacked bar per arm, WT at the top ----
# Four segments, not three. The old grey segment carried two different statements at once:
# a gene the orthology source knew and had no human counterpart for, and a gene whose symbol
# it could not key at all because this matrix carries GENCODE vM25's older MGI vintage. The
# second is a vocabulary result, and 146 genes of the projection background sat in it
# labelled as orthology losses. They get their own colour so the label cannot lie again.
LED_SEGS <- c(n_human_carried              = "human genes carried",
              n_lost_to_collapse           = "lost when paralogs collapsed",
              n_dropped_no_ortholog        = "dropped: keyed, no human ortholog",
              n_dropped_stale_query_symbol = "dropped: symbol not keyable")
LED_FILL <- c("human genes carried"          = FIG_CFG$colors$okabe_ito$blue   %||% "#0072B2",
              "lost when paralogs collapsed" = FIG_CFG$colors$okabe_ito$orange %||% "#E69F00",
              "dropped: keyed, no human ortholog" = "grey55",
              "dropped: symbol not keyable" =
                FIG_CFG$colors$okabe_ito$reddish_purple %||% "#CC79A7")
# Text sitting ON a segment: white reads on the dark carried fill, near-black on the light
# ones. The leader lines keep the segment's own colour so a parked number stays keyed to its
# sliver, except the grey fill which is darkened to stay visible as a line.
LED_INK  <- c("human genes carried"          = "white",
              "lost when paralogs collapsed" = "grey10",
              "dropped: keyed, no human ortholog" = "grey10",
              "dropped: symbol not keyable" = "grey10")
LED_LEAD <- c("human genes carried"           = unname(LED_FILL[1]),
              "lost when paralogs collapsed"  = unname(LED_FILL[2]),
              "dropped: keyed, no human ortholog" = "grey45",
              "dropped: symbol not keyable"   = unname(LED_FILL[4]))

conv_ledger$y    <- rev(seq_len(nrow(conv_ledger)))
conv_ledger$tick <- sprintf("%s\n%d mouse genes",
                            contrast_label(conv_ledger$contrast, short = TRUE),
                            conv_ledger$n_mouse_at_gate)

LED_BARH <- 0.44
led_seg <- conv_ledger %>%
  dplyr::select(arm, y, n_mouse_at_gate, dplyr::all_of(names(LED_SEGS))) %>%
  tidyr::pivot_longer(cols = dplyr::all_of(names(LED_SEGS)),
                      names_to = "key", values_to = "n") %>%
  dplyr::mutate(segment = factor(unname(LED_SEGS[key]), levels = unname(LED_SEGS))) %>%
  dplyr::group_by(y) %>%
  dplyr::arrange(segment, .by_group = TRUE) %>%
  dplyr::mutate(xmax = cumsum(n), xmin = xmax - n, xmid = (xmin + xmax) / 2) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(ymin = y - LED_BARH / 2, ymax = y + LED_BARH / 2)

# A 2-gene segment on a 213-gene bar is a sliver, so its count cannot sit inside it.
# Anything under LED_FIT_FRAC of the widest bar gets its number parked to the right of
# that bar with a leader line back to the sliver, in a fixed order so two slivers on one
# bar never land on top of each other. The threshold is set above the widest loss segment
# (19 of 239) so both loss counts are parked on every bar and the placement rule reads the
# same way on each row.
LED_AXIS_MAX <- max(conv_ledger$n_mouse_at_gate)
LED_FIT_FRAC <- 0.085
LED_W        <- 11
led_lab <- led_seg %>% dplyr::filter(n > 0) %>%
  dplyr::mutate(fits = n / LED_AXIS_MAX >= LED_FIT_FRAC)
led_in  <- led_lab %>% dplyr::filter(fits) %>%
  dplyr::mutate(ink = unname(LED_INK[as.character(segment)]))
led_out <- led_lab %>% dplyr::filter(!fits) %>%
  dplyr::group_by(y) %>%
  dplyr::arrange(segment, .by_group = TRUE) %>%
  dplyr::mutate(lx = n_mouse_at_gate + LED_AXIS_MAX * (0.055 + 0.095 * (dplyr::row_number() - 1))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(lead = unname(LED_LEAD[as.character(segment)]))
LED_XHI <- max(LED_AXIS_MAX * 1.06, if (nrow(led_out)) max(led_out$lx) * 1.10 else 0)

led_wrap <- function(txt, chars_per_inch)
  paste(strwrap(txt, width = as.integer(LED_W * chars_per_inch)), collapse = "\n")

fig_conv <- ggplot() +
  geom_rect(data = dplyr::filter(led_seg, n > 0),
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = segment)) +
  geom_segment(data = led_out,
               aes(x = xmid, xend = lx, y = ymax, yend = y + 0.28, colour = lead),
               linewidth = 0.45, show.legend = FALSE) +
  geom_text(data = led_in,
            aes(x = xmid, y = y, label = n, colour = ink),
            size = (FIG_CFG$figures$label_size %||% 4) * 0.95, fontface = "bold",
            show.legend = FALSE) +
  geom_text(data = led_out,
            aes(x = lx, y = y + 0.32, label = n),
            size = (FIG_CFG$figures$label_size %||% 4) * 0.9, fontface = "bold",
            colour = "grey15", vjust = 0, hjust = 0.5) +
  scale_fill_manual(values = LED_FILL, name = "fate of the mouse gene", drop = FALSE) +
  scale_colour_identity(guide = "none") +
  scale_y_continuous(breaks = conv_ledger$y, labels = conv_ledger$tick,
                     expand = expansion(add = c(0.42, 0.72))) +
  scale_x_continuous(limits = c(0, LED_XHI), expand = expansion(mult = c(0.004, 0))) +
  labs(title = "Fate of each mouse up arm under the ortholog map",
       subtitle = led_wrap(sprintf(
         "Bar length = mouse genes in the arm at the %s gate; segments = where those genes went.",
         GATE_PRIMARY), chars_per_inch = 11),
       x = "genes", y = NULL,
       caption = led_wrap(paste0("Counts recomputed from the frozen signature files and the ",
                                 "applied ortholog map, then cross-checked against the projection ",
                                 "manifest. Claim tier: a direct count over frozen sets."),
                          chars_per_inch = 13)) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank(),
        plot.caption.position = "plot")

LED_FINDING <- paste0(
  "Each frozen mouse up arm loses genes three ways: no human counterpart, a symbol the ",
  "orthology source could not key, and paralogs collapsing onto one human symbol. ",
  paste(sprintf("%s carries %d of %d (%d no ortholog, %d not keyable, %d collapsed, %d recovered by normalisation)",
                conv_ledger$arm, conv_ledger$n_human_carried, conv_ledger$n_mouse_at_gate,
                conv_ledger$n_dropped_no_ortholog,
                conv_ledger$n_dropped_stale_query_symbol, conv_ledger$n_lost_to_collapse,
                conv_ledger$n_query_symbol_normalised),
        collapse = "; "),
  ".\n\nOne-to-one arrivals: ",
  paste(sprintf("%d of %d", conv_ledger$n_arrived_one_to_one, conv_ledger$n_human_carried),
        collapse = ", "),
  " in the same order. ", led_hs_sentence)

purge_figures(STAGE, "conversion_ledger", overview = TRUE, config = FIG_CFG)
save_overview(
  fig_conv, STAGE, "conversion_ledger",
  table   = conv_ledger %>% dplyr::select(-y, -tick),
  finding = LED_FINDING,
  script  = SCRIPT, fn = "ggplot(geom_rect + leader-line sliver labels)",
  config_kv = paste0("decisions.projection.gate=", GATE_PRIMARY,
                     "; decisions.projection.contrasts_primary; ",
                     "decisions.projection.ortholog_ambiguity.min_support=", MIN_SUPPORT,
                     "; babelgene=", BABELGENE_VER, "; colors.okabe_ito"),
  input   = paste0(PROJ_DIR, "/{manifest.csv,ortholog_map.tsv,signatures/}; ", f_sig),
  how_to_read = paste0(
    "One horizontal bar per frozen up arm. Bar length is the mouse genes that passed the ",
    GATE_PRIMARY, " gate; the segments are their fate. Blue is carried into human, orange ",
    "is lost when several mouse paralogs collapsed onto one human symbol, grey is dropped ",
    "for want of an accepted ortholog, and purple is dropped because the orthology source ",
    "could not key the symbol at all. Those last two used to share one grey segment, which ",
    "read every vocabulary loss as a biological one; purple is the honest half of that bucket. ",
    "A count sits inside its segment where it fits, and is otherwise parked to the right of ",
    "the bar on a leader line back to the sliver.\n\n",
    "The map is babelgene ", BABELGENE_VER, ", queried mouse to human at min_support = ",
    MIN_SUPPORT, ", so an edge is accepted when at least three source databases agree; it ",
    "holds ", format(nrow(omap), big.mark = ","), " edges over ",
    format(dplyr::n_distinct(omap$mouse_symbol), big.mark = ","), " mouse and ",
    format(dplyr::n_distinct(omap$human_symbol), big.mark = ","), " human symbols. A binary ",
    "set takes the union both ways: several mouse onto one human shrinks it, one mouse ",
    "across several human grows it. Every human denominator downstream is the carried set ",
    "rather than the mouse set.\n\n",
    "Per-arm arrival routes (one-to-one, several mouse collapsed onto it, one mouse split ",
    "across several human) and the paralog groups that collided inside each arm are in the ",
    "same-stem CSV. Claim tier: a direct count over frozen sets."),
  width = LED_W, height = 5.0,
  config = FIG_CFG)

# ============================================================================
# 5e. CAPTION the one stage table that has no figure.
# ============================================================================
# query_normalisation_ledger.csv is a decision record rather than a plot: one row per matrix
# symbol babelgene could not key, with what org.Mm.eg.db proposed and whether the guards took
# it. It gets a caption here because 18_projection_export.R is compute-only and this stage's
# README is owned by its viz sibling. No figure is drawn from it deliberately — the numbers
# that matter already appear as the purple segment of conversion_ledger and mapping_loss, and
# a per-symbol table is read, not looked at.

f_qnl <- file.path(OV_DIR, "query_normalisation_ledger.csv")
if (file.exists(f_qnl)) {
  qnl <- readr::read_csv(f_qnl, show_col_types = FALSE, progress = FALSE)
  n_acc <- sum(qnl$resolution == "accepted")
  n_rec <- sum(qnl$mapped_after_normalisation)
  write_caption(
    STAGE, "tables/_overview/query_normalisation_ledger.csv",
    finding = sprintf(paste0("Of %d matrix symbols the orthology source could not key, %d ",
                             "resolved one-to-one to a current MGI symbol and %d of those then ",
                             "mapped to a human ortholog — genes that had been counted as ",
                             "having no human ortholog when the truth was that their name had ",
                             "changed. %d pairs were withheld by a guard and %d by human ",
                             "decision."),
                     nrow(qnl), n_acc, n_rec,
                     sum(!qnl$resolution %in% c("accepted", "flagged_for_review")),
                     sum(qnl$resolution == "flagged_for_review")),
    script = SCRIPT, fn = "helpers/ortholog_utils.R::normalise_mouse_query (written by 18_projection_export.R)",
    config_kv = "symbol_alias.ortholog_query_flagged_for_review; org.Mm.eg.db",
    input = "03_results/objects/17_signature_sets.rds (the modelled symbol universe)",
    how_to_read = paste0(
      "One row per candidate. `matrix_symbol` is the name the data carries, `current_symbol` ",
      "what org.Mm.eg.db calls the same Entrez id today, `resolution` what happened.\n\n",
      "`accepted` means the symbol resolved to exactly one Entrez id, the current symbol is ",
      "that id's own official symbol and no other gene's, no second candidate claimed it, and ",
      "it is not already a separate row here. `flagged_for_review` is withheld by decision, not ",
      "by a guard: Ndufb1-ps because handing a pseudogene-named row the real gene's ortholog ",
      "judges the vM25 annotation, H2aw because histone renaming is many-to-many.\n\n",
      "`rejected_current_symbol_already_in_vocabulary` is the guard that matters most quietly. ",
      "Normalising there would give one gene two rows and let the many-mouse-to-one-human rule ",
      "silently replace an existing gene's ranked metric.\n\n",
      "`mapped_after_normalisation` splits the two outcomes an accepted pair can have. TRUE is ",
      "a recovery. FALSE means the symbol was re-asked correctly and still has no human ",
      "ortholog at this min_support — a real orthology result, not a vocabulary loss. Claim ",
      "tier: a direct count over an annotation database, no statistics."),
    config = FIG_CFG)
  message(sprintf("[18_viz] captioned query_normalisation_ledger.csv (%d candidates, %d recovered).",
                  nrow(qnl), n_rec))
}

# ============================================================================
# 6. FINAL SUMMARY
# ============================================================================

n_fig <- length(list.files(file.path("03_results", STAGE, "figures"),
                           pattern = "\\.(pdf|png)$", recursive = TRUE))
message(sprintf("[18_viz] COMPLETE: %d figure file(s) under %s/figures/.", n_fig, STAGE))
if (n_fig == 0) warning("[18_viz] No figures produced — check errors above.")
