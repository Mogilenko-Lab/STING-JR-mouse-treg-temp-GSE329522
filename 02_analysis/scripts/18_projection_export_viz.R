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
#
# Inputs (read-only; produced by 18_projection_export.R):
#   03_results/11_projection/tables/_overview/{human_signature_sizes,mapping_loss,projection_overlap_ledger}.csv
#   03_results/human_projection/manifest.csv
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
LOSS_COLORS <- c(
  n_mapped_1to1 = FIG_CFG$colors$okabe_ito$bluish_green %||% "#009E73",
  n_many_mapped = FIG_CFG$colors$okabe_ito$orange       %||% "#E69F00",
  n_unmapped    = "grey65")
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
  tidyr::pivot_longer(cols = c(n_mapped_1to1, n_many_mapped, n_unmapped),
                      names_to = "category", values_to = "n") %>%
  dplyr::mutate(category = factor(category,
                                  levels = c("n_mapped_1to1", "n_many_mapped", "n_unmapped")))

# human-out annotation + alarm flag per unit/direction
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
                    labels = c(n_mapped_1to1 = "mapped 1:1", n_many_mapped = "one mouse -> many human",
                               n_unmapped = "dropped (no human ortholog)"),
                    name = "mouse gene fate") +
  scale_x_discrete(labels = unit_label) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Ortholog mapping loss: mouse set -> human set",
       subtitle = paste0("Each bar = the mouse up/down set decomposed by fate; '->N' = resulting human ",
                        "set size.\nGates = ", GATE_NOTE, ". Red * = > ",
                        round(100 * ALARM_FRAC), "% dropped (alarm)."),
       x = NULL, y = "mouse genes in set",
       caption = paste0("Bar height = mouse-set size; green kept 1:1, orange one->many (each human gets ",
                        "the mouse t), grey dropped. Claim tier: L3. Correlative input, not a causal claim.")) +
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
                       "human-space counterpart of the mouse-side size bars. Claim tier: L3 (DE statistics); ",
                       "provisional sample mapping, n=5/group."),
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
                       "fate: green = mapped 1:1, orange = one-mouse->many-human (each human ortholog inherits ",
                       "the mouse t in the ranked list; unioned in the binary set), grey = dropped (no human ",
                       "ortholog). A contrast carried at two threshold gates gets one bar per gate, with the ",
                       "gate named in square brackets under its tick. Bar height = mouse-set size; '->N' above ",
                       "= resulting human-set size. A red * marks > ", round(100 * ALARM_FRAC),
                       "% dropped — a decimated signature to demote. Read it as the sanity check that the ",
                       "frozen human sets are not hollow. Claim tier: L3."),
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
# 6. FINAL SUMMARY
# ============================================================================

n_fig <- length(list.files(file.path("03_results", STAGE, "figures"),
                           pattern = "\\.(pdf|png)$", recursive = TRUE))
message(sprintf("[18_viz] COMPLETE: %d figure file(s) under %s/figures/.", n_fig, STAGE))
if (n_fig == 0) warning("[18_viz] No figures produced — check errors above.")
