#!/usr/bin/env Rscript
# 27_arm_composition_viz.R -- VIZ ONLY (no computing)
# =============================================================================
# Draws the six overview figures for 16_arm_composition from the frozen object
# written by 26_arm_composition.R. Every number on every face is read from that
# object; nothing here recomputes a share, a p, or a null.
#
#   arm_composition           the composition itself, per arm, both accountings
#   arm_composition_variants  `unpinned` beside `hypoxia_anchored`, pins marked
#   arm_hypoxia_sources       the three sets carrying hypoxia hits, and their union
#   arm_remainder             the unclaimed genes, split by reason
#   arm_composition_null      observed against the depth-matched null
#   interaction_up_genes      the 7-gene arm, drawn gene by gene
#
# Two accountings are always drawn together and neither is reduced to the other.
# `fractional` gives a gene in k selected sets 1/k to each. `winner_take_all`
# gives each gene to its single best selected set. Both sum to 1.0 over the arm,
# so where the two disagree the disagreement is the quantity on the page.
#
# Guards, and the silent failures they exist for. Each one fired here for real, and none of
# them produced an error on its own. Read this before loosening any of them.
#
#   VARIANTS + the per-frame check     A factor level renamed upstream. A filter on the absent
#                                      level returns ZERO ROWS, not an error, so the figure
#                                      renders as a convincing empty panel and a reader
#                                      filtering the published table gets nothing back. Frames
#                                      this script draws from hard-stop; the rest warn loudly,
#                                      because blocking a correct figure over a table this
#                                      script never touches helps nobody.
#   dup_free()                         A duplicated key. Every caption number is pulled by a
#                                      logical mask expecting one row; a duplicate returns
#                                      length 2, sprintf recycles it in silence, and it
#                                      surfaces far downstream as a coercion error inside the
#                                      caption writer.
#   nz()                               A per-arm subset that comes back empty must stop before
#                                      it reaches ggplot.
#   INPUT carries the object md5       A partial re-run once left five figures from a new
#                                      object, one from the old, and a README from the new,
#                                      with no error and a directory that looked complete.
#                                      Only timestamps disagreed. The stamp makes the check
#                                      one `md5sum` against the caption.
#   `%in%` on the null level           `df[NA, ]` returns an all-NA ROW rather than dropping
#                                      it, which then ships in a sibling table while the geoms
#                                      quietly drop it from the panel.
#   row_order()                        Exact ties in the sort key, broken by input order,
#                                      make vertical position imply a rank the data does not
#                                      support. There are two such ties in WT_heat_up.
#   Captions derive, never name        Two hand-written captions went stale within an hour.
#                                      The leader is computed, and a tied group is reported
#                                      as a group.
#
# A longer write-up of the panel design sits beside the other reasoning notes, in a directory
# this repository does not track, so this header is the tracked copy of anything load-bearing.
#
# Run from project root:
#   Rscript 02_analysis/scripts/27_arm_composition_viz.R
# =============================================================================

source("02_analysis/config/config.R")
source("02_analysis/helpers/figure_style.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})
options(stringsAsFactors = FALSE)

STAGE  <- "16_arm_composition"
SCRIPT <- "02_analysis/scripts/27_arm_composition_viz.R"
OBJ    <- "03_results/objects/26_arm_composition.rds"
# Stamp the source object's md5 into every caption. Twice now a figure has sat beside a table
# regenerated from a newer object, with matching-looking content and no error anywhere; the
# only tell was a file timestamp. A hash on the page makes the check one command instead:
# compare it against `md5sum` of the object.
INPUT  <- sprintf("%s (md5 %s)", OBJ, unname(tools::md5sum(OBJ)))
invisible(set_paper_style(config = FIG_CFG))

obj        <- readRDS(OBJ)
SHARES     <- obj$shares
ASSIGN     <- obj$assignment
REMAINDER  <- obj$remainder
HYP        <- obj$hypoxia
HYP_OVL    <- obj$hypoxia_overlap
NULLS      <- obj$null_summary
ARMS       <- obj$arms
RESIDUAL   <- obj$residual
PROV       <- setNames(obj$provenance$value, obj$provenance$key)

OI  <- FIG_CFG$colors$okabe_ito
LAB <- (FIG_CFG$figures$label_size %||% 4)
PT  <- (FIG_CFG$figures$point_size %||% 2.4)

COLLECTIONS <- c("Hallmark", "GO_BP", "GO_MF", "KEGG", "Reactome")
ROLLUP_KEY  <- "__rollup__"
OTHER_ROW   <- "all other selected sets"
ARM_ORDER   <- names(ARMS)

# The two variant names are the compute stage's vocabulary, held in one place so a rename
# upstream is a one-line change here. The guard below matters more than the constants: a
# filter on a renamed level returns ZERO ROWS rather than an error, and a figure drawn from
# an empty frame renders as an empty panel with a clean exit code. Stop instead.
VARIANT_MAIN <- "unpinned"
VARIANT_PIN  <- "hypoxia_anchored"
VARIANTS     <- c(VARIANT_MAIN, VARIANT_PIN)
# Check EVERY frame in the object that carries a `variant` column, not only the ones drawn
# from. A partial rename leaves the published tables disagreeing with each other, and a reader
# who filters one of them on the other's label gets zero rows and no error. Frames this script
# consumes are a hard stop. Frames it does not consume are reported loudly and let the run
# finish, because the figures are still correct and blocking them helps nobody.
CONSUMED <- c("shares", "assignment", "remainder")
for (nm in names(obj)) {
  fr <- obj[[nm]]
  if (!is.data.frame(fr) || !"variant" %in% names(fr)) next
  seen <- sort(unique(fr$variant))
  if (setequal(seen, VARIANTS)) next
  msg <- sprintf("variant vocabulary in `%s` is {%s}, this stage's is {%s}",
                 nm, paste(seen, collapse = ", "), paste(sort(VARIANTS), collapse = ", "))
  if (nm %in% CONSUMED) stop(msg, call. = FALSE)
  warning(msg, call. = FALSE, immediate. = TRUE)
}
# Same failure mode one level down: a per-arm subset that comes back empty must stop here
# rather than reach ggplot.
nz <- function(df, what) {
  if (!nrow(df)) stop("no rows for ", what, call. = FALSE)
  df
}

# And the mirror of it. Every number quoted in a caption below is pulled out of a frame by a
# logical mask that is expected to match exactly one row. A duplicated key returns a length-2
# vector instead, sprintf recycles it into a length-2 string without complaining, and the
# failure only surfaces further downstream as a coercion error inside the caption writer. The
# grains are therefore asserted once, here, where the diagnostic names the frame.
dup_free <- function(df, keys, what) {
  d <- anyDuplicated(df[, keys, drop = FALSE])
  if (d) stop(sprintf("%s is not unique on (%s): first repeat at row %d",
                      what, paste(keys, collapse = ", "), d), call. = FALSE)
  invisible(TRUE)
}
dup_free(SHARES, c("arm", "variant", "accounting", "ID"), "shares")
dup_free(REMAINDER, c("arm", "variant"), "remainder")
dup_free(NULLS, c("arm", "collection", "null"), "null_summary")
dup_free(HYP, c("arm", "collection", "term"), "hypoxia")
dup_free(HYP_OVL, c("arm", "term_a", "term_b"), "hypoxia_overlap")
# Interaction_up carries 7 genes. A five-collection decomposition of 7 genes is not a
# composition, so it gets a gene-level panel of its own and no share bar anywhere.
BAR_ARMS    <- setdiff(ARM_ORDER, "Interaction_up")
TOP_N       <- 10L

CONFIG_KV <- sprintf(paste0("arm_composition.collections = %s, p_matched_cutoff = %s, ",
                            "max_null_recurrence = %s, prune_hit_jaccard = %s, ",
                            "top_n_per_collection = %s, n_null = %s, seed = %s, ",
                            "winner_take_all_tiebreak = %s"),
                     PROV[["collections"]], PROV[["p_matched_cutoff"]],
                     PROV[["max_null_recurrence"]], PROV[["prune_hit_jaccard"]],
                     PROV[["top_n_per_collection"]], PROV[["n_null"]], PROV[["seed"]],
                     PROV[["winner_take_all_tiebreak"]])

## ---------------------------------------------------------------------------
## Shared vocabulary
## ---------------------------------------------------------------------------

# Category colour is the collection the set came from. The two book-keeping categories
# that are never a set (the residue, and the roll-up of the sets below the drawn cap)
# take greys so they read as book-keeping. KEGG takes the purple rather than the second
# orange: Hallmark and KEGG bands sit adjacent in the stacked panel and two oranges
# there are one glyph, not two.
CAT_COL <- c(Hallmark = OI$vermillion, GO_BP = OI$blue, GO_MF = OI$sky_blue,
             KEGG = OI$reddish_purple, Reactome = OI$bluish_green)
CAT_COL[[RESIDUAL]] <- "grey72"
CAT_COL[["further selected sets"]] <- "grey38"
CAT_LEVELS <- names(CAT_COL)

ACC_SHAPE <- c(fractional = 21L, winner_take_all = 24L)
ACC_ORDER <- c("fractional", "winner_take_all")

arm_n <- vapply(ARMS, function(a) a$n_in_background, integer(1))
arm_strip  <- function(a) sprintf("%s   %d genes", a, unname(arm_n[as.character(a)]))
arm_strip2 <- function(a) sprintf("%s\n%d genes", a, unname(arm_n[as.character(a)]))

# Wrap a long MSigDB identifier onto several lines WITHOUT dropping any of it. The axis
# label has to stay the checkable name, so this breaks on the separators the identifier
# already carries rather than cutting it at a character count.
#
# The underscore is KEPT at the end of the broken line. Dropping it renders
# GOMF_DNA_BINDING_TRANSCRIPTION_ACTIVATOR_ACTIVITY as two lines whose join is a string
# that does not exist, so a reader copying the label off the page gets a name that matches
# nothing. Only a trailing space is stripped.
wrap_id <- function(x, width = 46L) {
  vapply(x, function(s) {
    if (is.na(s) || !nzchar(s)) return("")
    toks <- regmatches(s, gregexpr("[^_ ]+[_ ]?", s))[[1]]
    out <- character(0); cur <- ""
    for (tk in toks) {
      cand <- paste0(cur, tk)
      if (nchar(trimws(cand)) > width && nzchar(cur)) {
        out <- c(out, sub(" $", "", cur)); cur <- tk
      } else cur <- cand
    }
    paste(c(out, sub(" $", "", cur)), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}

pct  <- function(x, acc = 0.1) scales::percent(x, accuracy = acc)
num3 <- function(x) formatC(x, format = "f", digits = 3)

wide_shares <- function(df) {
  df %>%
    dplyr::select("arm", "variant", "ID", "collection", "accounting", "share", "n_genes") %>%
    tidyr::pivot_wider(names_from = "accounting",
                       values_from = c("share", "n_genes"), values_fill = 0)
}

n_cats <- function(arm_id, variant_id, acc)
  sum(SHARES$arm == arm_id & SHARES$variant == variant_id & SHARES$accounting == acc)

## ---------------------------------------------------------------------------
## Figure 1 -- the composition itself
##
## Panel A keeps every gene on the page at collection resolution: five collections plus
## the residue, six bands, exactly 1.0 per bar. Panels B to D go to set resolution, where
## there are up to 42 categories per arm. A 42-band stack is unreadable, so the drawn
## rows are the top ten by the larger of the two accountings and everything below that cap
## is carried as one roll-up row that states how many sets it holds. Nothing is dropped:
## the roll-up plus the residue plus the drawn rows still sum to 1.0, and the full
## per-set list is the sibling table.
##
## The set-level marks are paired points rather than bars because 12 of WT_heat_up's 40
## categories take a share of exactly zero under one accounting. A bar of height zero is
## an axis label with nothing beside it. A point at zero is a mark a reader can see, and
## the segment joining the pair is the disagreement between the two accountings.
## ---------------------------------------------------------------------------

shares_u <- SHARES %>% dplyr::filter(.data$variant == VARIANT_MAIN)

collA <- shares_u %>%
  dplyr::filter(.data$arm %in% BAR_ARMS) %>%
  dplyr::group_by(.data$arm, .data$accounting, .data$collection) %>%
  dplyr::summarise(share = sum(.data$share), .groups = "drop") %>%
  dplyr::mutate(
    collection = factor(.data$collection, levels = CAT_LEVELS),
    accounting = factor(.data$accounting, levels = rev(ACC_ORDER)),
    arm = factor(arm_strip(.data$arm), levels = vapply(BAR_ARMS, arm_strip, character(1))))

# The roll-up category exists only in panels B to D, and panel A carries the one key
# collected for the whole figure. Without a row of its own that key comes out as a label
# with no swatch beside it, so a zero-height row is added to give it one.
collA <- dplyr::bind_rows(
  collA,
  dplyr::mutate(collA[1, ],
                collection = factor("further selected sets", levels = CAT_LEVELS),
                share = 0))

pA <- ggplot(collA, aes(x = .data$share, y = .data$accounting, fill = .data$collection)) +
  geom_col(width = 0.62, colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(.data$share >= 0.055, pct(.data$share, 1), "")),
            position = position_stack(vjust = 0.5), size = LAB * 0.85,
            fontface = "bold", colour = "grey12", show.legend = FALSE) +
  facet_wrap(~ .data$arm, ncol = 1, strip.position = "top") +
  scale_fill_manual(values = CAT_COL, name = NULL, drop = FALSE, limits = CAT_LEVELS) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), expand = expansion(mult = c(0, 0.01))) +
  labs(title = "A. Share of each arm by collection",
       subtitle = sprintf(paste0("Interaction_up, %d genes, carries no share bar anywhere ",
                                 "and is drawn gene by gene in interaction_up_genes."),
                          unname(arm_n[["Interaction_up"]])),
       x = "share of the arm", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE))

# Rows sort on the larger of the two shares. Exact ties happen and they are not a rounding
# artefact: two categories can hold the same number of genes under winner-take-all, in which
# case `peak` is identical to the last bit. Left to `order(-peak)` alone the tie would be
# broken by whatever order pivot_wider happened to emit, and a reader takes vertical position
# for rank. The tie-break is therefore stated and reproducible: fractional share, then the
# identifier. It decides position only, never a share.
row_order <- function(d)
  order(-d$peak, -d$share_fractional, d$ID)

set_rows <- function(arm_id, variant_id = VARIANT_MAIN, top_n = TOP_N) {
  w <- SHARES %>%
    dplyr::filter(.data$arm == arm_id, .data$variant == variant_id) %>%
    nz(sprintf("arm %s, variant %s", arm_id, variant_id)) %>%
    wide_shares() %>%
    dplyr::mutate(peak = pmax(.data$share_fractional, .data$share_winner_take_all))
  res <- w[w$ID == RESIDUAL, , drop = FALSE]
  sets <- w[w$ID != RESIDUAL, , drop = FALSE]
  sets <- sets[row_order(sets), , drop = FALSE]
  drawn <- utils::head(sets, top_n)
  held <- if (nrow(sets) > top_n) sets[(top_n + 1L):nrow(sets), , drop = FALSE] else sets[0, ]
  roll <- if (nrow(held)) data.frame(
    arm = arm_id, variant = variant_id, ID = ROLLUP_KEY,
    collection = "further selected sets",
    share_fractional = sum(held$share_fractional),
    share_winner_take_all = sum(held$share_winner_take_all),
    n_genes_fractional = NA_integer_, n_genes_winner_take_all = NA_integer_,
    peak = max(held$share_fractional, held$share_winner_take_all),
    n_held = nrow(held)) else NULL
  out <- dplyr::bind_rows(
    dplyr::mutate(drawn, n_held = NA_integer_),
    roll,
    dplyr::mutate(res, n_held = NA_integer_))
  out$label <- ifelse(out$ID == ROLLUP_KEY,
                      sprintf("%d further selected sets", out$n_held), out$ID)
  out[row_order(out), , drop = FALSE]
}

SET_ROWS <- dplyr::bind_rows(lapply(BAR_ARMS, set_rows))

# How many categories of the primary arm take a share of exactly zero under one of the two
# accountings. This is the number that rules out bars for the set-level panels.
wt_wide <- wide_shares(SHARES[SHARES$arm == "WT_heat_up" & SHARES$variant == VARIANT_MAIN, ])
wt_wide <- wt_wide[wt_wide$ID != RESIDUAL, ]
N_CAT_WT  <- nrow(wt_wide)
N_ZERO_WT <- sum(wt_wide$share_fractional == 0 | wt_wide$share_winner_take_all == 0)

X_HI  <- 0.80
X_TXT <- 0.605

set_panel <- function(arm_id, tag) {
  d <- SET_ROWS[SET_ROWS$arm == arm_id, , drop = FALSE]
  d$label <- factor(as.character(d$label), levels = rev(as.character(d$label)))
  long <- d %>%
    dplyr::select("label", "collection", "share_fractional", "share_winner_take_all") %>%
    tidyr::pivot_longer(cols = c("share_fractional", "share_winner_take_all"),
                        names_to = "accounting", values_to = "share") %>%
    dplyr::mutate(accounting = factor(sub("^share_", "", .data$accounting),
                                      levels = ACC_ORDER),
                  collection = factor(.data$collection, levels = CAT_LEVELS))
  d$pair <- sprintf("%s  /  %s", num3(d$share_fractional), num3(d$share_winner_take_all))

  ggplot(long, aes(y = .data$label)) +
    geom_segment(data = d,
                 aes(x = .data$share_fractional, xend = .data$share_winner_take_all,
                     y = .data$label, yend = .data$label),
                 colour = "grey60", linewidth = 1.1) +
    geom_point(aes(x = .data$share, fill = .data$collection, shape = .data$accounting),
               size = PT * 1.7, colour = "grey15", stroke = 0.8) +
    geom_text(data = d, aes(x = X_TXT, y = .data$label, label = .data$pair),
              hjust = 0, size = LAB * 0.85, colour = "grey25", family = "mono") +
    scale_fill_manual(values = CAT_COL, drop = FALSE, limits = CAT_LEVELS,
                      guide = "none") +
    scale_shape_manual(values = ACC_SHAPE, name = NULL, drop = FALSE,
                       limits = ACC_ORDER) +
    scale_y_discrete(labels = function(x) wrap_id(x, 46L)) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, X_HI), breaks = seq(0, 0.5, 0.1),
                       expand = expansion(mult = c(0.01, 0))) +
    labs(title = sprintf("%s. %s, %d genes, %d categories fractional and %d winner-take-all",
                         tag, arm_id, unname(arm_n[[arm_id]]),
                         n_cats(arm_id, VARIANT_MAIN, "fractional"),
                         n_cats(arm_id, VARIANT_MAIN, "winner_take_all")),
         x = "share of the arm", y = NULL) +
    project_theme(config = FIG_CFG) +
    theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
          axis.text.y = element_text(lineheight = 0.92)) +
    guides(shape = guide_legend(
      override.aes = list(fill = "grey55", colour = "grey15", size = PT * 1.7)))
}

pB <- set_panel(BAR_ARMS[1], "B")
pC <- set_panel(BAR_ARMS[2], "C")
pD <- set_panel(BAR_ARMS[3], "D")

fig1 <- (pA / pB / pC / pD) +
  patchwork::plot_layout(heights = c(0.80, 1, 1, 1), guides = "collect") +
  patchwork::plot_annotation(
    title = "Pathway shares of the heat-contrast arms, fractional beside winner-take-all",
    theme = project_theme(config = FIG_CFG)) &
  theme(legend.position = "bottom", legend.box = "vertical")

fig1_table <- SET_ROWS %>%
  dplyr::transmute(arm = .data$arm, variant = .data$variant,
                   drawn_row = .data$label, category = .data$ID,
                   collection = .data$collection,
                   share_fractional = .data$share_fractional,
                   share_winner_take_all = .data$share_winner_take_all,
                   n_genes_fractional = .data$n_genes_fractional,
                   n_genes_winner_take_all = .data$n_genes_winner_take_all,
                   n_sets_rolled_up = .data$n_held)

# Which category leads is exactly what the winner-take-all tie-break decides, so it is read
# off the data rather than named here. A hard-coded leader outlives the re-run that changes it.
wt_named <- SET_ROWS[SET_ROWS$arm == "WT_heat_up" &
                       !SET_ROWS$ID %in% c(RESIDUAL, ROLLUP_KEY), , drop = FALSE]
lead_frac <- wt_named[which.max(wt_named$share_fractional), , drop = FALSE]
top_wta   <- max(wt_named$share_winner_take_all)
lead_wta  <- wt_named[abs(wt_named$share_winner_take_all - top_wta) < 1e-12, , drop = FALSE]
lead_wta  <- lead_wta[order(lead_wta$ID), , drop = FALSE]
# Equal winner-take-all share means equal gene count, because that accounting is a gene count
# over the arm size. Assert it rather than assume it before quoting one number for the group.
stopifnot(length(unique(lead_wta$n_genes_winner_take_all)) == 1L)

wta_clause <- if (nrow(lead_wta) > 1L) {
  sprintf(paste0("%d categories tie exactly at the head of `winner_take_all`, %s, each on %d ",
                 "genes and each taking %s of the arm, and they come from collections of very ",
                 "different size, so that accounting names no single leading pathway here"),
          nrow(lead_wta), paste(lead_wta$ID, collapse = " and "),
          lead_wta$n_genes_winner_take_all[1], num3(top_wta))
} else {
  sprintf("`winner_take_all` is led by %s on %d genes, %s of the arm",
          lead_wta$ID[1], lead_wta$n_genes_winner_take_all[1], num3(top_wta))
}

save_overview(
  fig1, STAGE, "arm_composition",
  table = fig1_table,
  finding = sprintf(paste0(
    "Under the `%s` variant %d of %d genes of WT_heat_up sit in no selected set. The largest ",
    "named category under `fractional` is %s at %s of the arm, and the same set takes %s under ",
    "`winner_take_all`. %s. Which pathway leads the composition therefore depends on how a gene ",
    "held by several sets is counted, which is why both readings are drawn and neither is ",
    "called primary. KO_heat_up leaves %d of %d unclaimed and Interaction_fdrOnly_up %d of %d."),
    VARIANT_MAIN,
    REMAINDER$n_unclaimed[REMAINDER$arm == "WT_heat_up" & REMAINDER$variant == VARIANT_MAIN],
    unname(arm_n[["WT_heat_up"]]),
    lead_frac$ID, num3(lead_frac$share_fractional),
    num3(lead_frac$share_winner_take_all), wta_clause,
    REMAINDER$n_unclaimed[REMAINDER$arm == "KO_heat_up" & REMAINDER$variant == VARIANT_MAIN],
    unname(arm_n[["KO_heat_up"]]),
    REMAINDER$n_unclaimed[REMAINDER$arm == "Interaction_fdrOnly_up" &
                            REMAINDER$variant == VARIANT_MAIN],
    unname(arm_n[["Interaction_fdrOnly_up"]])),
  script = SCRIPT, fn = "set_panel() / top-level (pA-pD)",
  config_kv = CONFIG_KV, input = INPUT,
  how_to_read = sprintf(paste0(
    "Panel A is every gene at collection resolution: six bands, two bars per arm, each bar ",
    "summing to 1.0. Panels B to D go to set resolution. Each row is one category with two ",
    "marks, a circle for `fractional` where a gene in k selected sets gives 1/k to each, and a ",
    "triangle for `winner_take_all` where the gene goes whole to its single best set. The grey ",
    "segment between them is the size of the disagreement and both readings are reported. A mark ",
    "on the zero line is a real zero under that accounting, which is why the set-level panels use ",
    "paired points and not bars: %d of WT_heat_up's %d categories take a share of exactly zero ",
    "under one of the two. The pair of numbers on the right repeats each row as fractional / ",
    "winner-take-all, and fill is the collection the set came from. Rows are ordered by the ",
    "larger of the two shares, so the residue and the roll-up sort in among the sets. Only the ",
    "ten largest sets per arm are drawn by name, and the roll-up row states how many further sets ",
    "it holds and carries their summed share, so the drawn rows plus the roll-up plus the residue ",
    "still total 1.0. Every rolled-up set is named in composition_shares.csv. Claim tier: ",
    "descriptive composition of a derived gene set, hypothesis-generating."), N_ZERO_WT, N_CAT_WT),
  config = FIG_CFG, wide = TRUE, height = 22.5)

## ---------------------------------------------------------------------------
## Figure 2 -- `unpinned` beside `hypoxia_anchored`
##
## Anchoring pins the hypoxia sets whether or not they enriched, so a set that never
## reached significance can take a share. Those rows carry NA significance and are
## labelled as such on the panel, so a mark there cannot be read as a finding.
##
## Panel A is drawn on a 0 to 7% axis because that is the range the pinned shares live
## in; panel B is drawn on a 0 to 100% axis because the residue and the rest of the
## composition live there. Splitting the two is what makes both readable at once.
## ---------------------------------------------------------------------------

# The pinned sets the engine could actually test. The fifth pin (the HIF1A signalling GO
# BP set) has 7 background genes, below the 10-gene floor, so it was never testable and
# is reported by figure 3 panel C as an absence of opportunity rather than as a zero here.
PIN_TESTED <- HYP %>%
  dplyr::filter(!is.na(.data$term), grepl("^tested", .data$status)) %>%
  dplyr::distinct(.data$term, .data$collection)

grab_share <- function(w, id, col) { v <- w[[col]][w$ID == id]; if (length(v)) v[1] else 0 }

pin_rows <- dplyr::bind_rows(lapply(BAR_ARMS, function(a) {
  u <- wide_shares(SHARES[SHARES$arm == a & SHARES$variant == VARIANT_MAIN, ])
  h <- wide_shares(SHARES[SHARES$arm == a & SHARES$variant == "hypoxia_anchored", ])
  dplyr::bind_rows(lapply(seq_len(nrow(PIN_TESTED)), function(i) {
    id <- PIN_TESTED$term[i]
    st <- HYP$status[HYP$arm == a & !is.na(HYP$term) & HYP$term == id][1]
    data.frame(arm = a, row = id, kind = "pinned set", status = st,
               fractional_unpinned = grab_share(u, id, "share_fractional"),
               fractional_anchored = grab_share(h, id, "share_fractional"),
               winner_take_all_unpinned = grab_share(u, id, "share_winner_take_all"),
               winner_take_all_anchored = grab_share(h, id, "share_winner_take_all"),
               n_other = NA_integer_)
  }))
}))

bulk_rows <- dplyr::bind_rows(lapply(BAR_ARMS, function(a) {
  u <- wide_shares(SHARES[SHARES$arm == a & SHARES$variant == VARIANT_MAIN, ])
  h <- wide_shares(SHARES[SHARES$arm == a & SHARES$variant == "hypoxia_anchored", ])
  oth <- setdiff(union(u$ID, h$ID), c(PIN_TESTED$term, RESIDUAL))
  data.frame(
    arm = a, row = c(RESIDUAL, OTHER_ROW), kind = "book-keeping",
    status = c("the residue", "every category anchoring does not pin"),
    fractional_unpinned = c(grab_share(u, RESIDUAL, "share_fractional"),
                            sum(u$share_fractional[u$ID %in% oth])),
    fractional_anchored = c(grab_share(h, RESIDUAL, "share_fractional"),
                            sum(h$share_fractional[h$ID %in% oth])),
    winner_take_all_unpinned = c(grab_share(u, RESIDUAL, "share_winner_take_all"),
                                 sum(u$share_winner_take_all[u$ID %in% oth])),
    winner_take_all_anchored = c(grab_share(h, RESIDUAL, "share_winner_take_all"),
                                 sum(h$share_winner_take_all[h$ID %in% oth])),
    n_other = c(NA_integer_, length(oth)))
}))

VAR_ROWS <- dplyr::bind_rows(pin_rows, bulk_rows)
VAR_ROWS$arm_f <- factor(arm_strip2(VAR_ROWS$arm),
                         levels = vapply(BAR_ARMS, arm_strip2, character(1)))
VAR_ROWS$never_enriched <- VAR_ROWS$status == "tested, not enriched"

# The `_unpinned` / `_anchored` strings below are suffixes of column names this script builds
# for the pivot, not the compute stage's vocabulary. The variant guard at the top is what
# protects against a rename upstream.
to_long <- function(df) {
  df %>%
    tidyr::pivot_longer(cols = c("fractional_unpinned", "fractional_anchored",
                                 "winner_take_all_unpinned", "winner_take_all_anchored"),
                        names_to = "key", values_to = "share") %>%
    dplyr::mutate(accounting = ifelse(grepl("^fractional", .data$key),
                                      "fractional", "winner_take_all"),
                  variant = factor(ifelse(grepl("unpinned$", .data$key),
                                          "unpinned", "hypoxia_anchored"),
                                   levels = c("unpinned", "hypoxia_anchored")))
}

to_seg <- function(df) {
  df %>%
    dplyr::select("arm_f", "row", "fractional_unpinned", "fractional_anchored",
                  "winner_take_all_unpinned", "winner_take_all_anchored") %>%
    tidyr::pivot_longer(cols = -c("arm_f", "row"), names_to = "key", values_to = "share") %>%
    dplyr::mutate(accounting = ifelse(grepl("^fractional", .data$key),
                                      "fractional", "winner_take_all"),
                  variant = ifelse(grepl("unpinned$", .data$key),
                                   "unpinned", "hypoxia_anchored")) %>%
    dplyr::select(-"key") %>%
    tidyr::pivot_wider(names_from = "variant", values_from = "share")
}

VAR_COL <- c(unpinned = "grey35", hypoxia_anchored = OI$reddish_purple)

variant_panel <- function(df, x_hi, x_txt, breaks, tag, ttl, tag_fun, wrap_w = 30L) {
  long <- to_long(df); seg <- to_seg(df)
  # The label has to follow the accounting, or the winner-take-all column would be
  # annotated with the fractional pair sitting beside a different point.
  lab <- dplyr::left_join(
    seg, df[, c("arm_f", "row", "never_enriched", "status", "n_other")],
    by = c("arm_f", "row"))
  lab$tag <- tag_fun(lab)
  ggplot(long, aes(y = .data$row)) +
    geom_segment(data = seg, aes(x = .data$unpinned, xend = .data$hypoxia_anchored,
                                 y = .data$row, yend = .data$row),
                 colour = "grey70", linewidth = 1.2) +
    geom_point(aes(x = .data$share, colour = .data$variant, shape = .data$variant),
               size = PT * 1.6, stroke = 1.2, fill = "white") +
    geom_text(data = lab, aes(x = x_txt, y = .data$row, label = .data$tag,
                              colour = ifelse(.data$never_enriched,
                                              "hypoxia_anchored", "unpinned")),
              hjust = 0, vjust = 0.5, size = LAB * 0.80, family = "mono",
              lineheight = 1.05, show.legend = FALSE) +
    facet_grid(.data$arm_f ~ .data$accounting) +
    scale_colour_manual(values = VAR_COL, name = NULL) +
    scale_shape_manual(values = c(unpinned = 21L, hypoxia_anchored = 19L), name = NULL) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                       limits = c(0, x_hi), breaks = breaks,
                       expand = expansion(mult = c(0.015, 0))) +
    scale_y_discrete(labels = function(x) wrap_id(x, wrap_w)) +
    labs(title = sprintf("%s. %s", tag, ttl), x = "share of the arm", y = NULL) +
    project_theme(config = FIG_CFG) +
    theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
          strip.text.y = element_text(angle = 0),
          axis.text.y = element_text(lineheight = 0.92))
}

pins <- VAR_ROWS[VAR_ROWS$kind == "pinned set", ]
pin_ord <- pins %>%
  dplyr::group_by(.data$row) %>%
  dplyr::summarise(peak = max(.data$fractional_anchored), .groups = "drop") %>%
  dplyr::arrange(.data$peak)
pins$row <- factor(pins$row, levels = pin_ord$row)
pV1 <- variant_panel(
  pins, x_hi = 0.088, x_txt = 0.042, breaks = c(0, 0.01, 0.02, 0.03), tag = "A",
  ttl = "The pinned hypoxia sets, on a 0 to 3% axis",
  tag_fun = function(d) sprintf("%s -> %s\n%s", num3(d$unpinned),
                                num3(d$hypoxia_anchored),
                                ifelse(d$never_enriched, "not enriched, NA", "enriched")))

bulk <- VAR_ROWS[VAR_ROWS$kind == "book-keeping", ]
bulk$row <- factor(bulk$row, levels = c(OTHER_ROW, RESIDUAL))
pV2 <- variant_panel(
  bulk, x_hi = 1.72, x_txt = 1.04, breaks = seq(0, 1, 0.25), tag = "B", wrap_w = 42L,
  ttl = "The residue and everything anchoring leaves alone, on a 0 to 100% axis",
  tag_fun = function(d) ifelse(
    is.na(d$n_other),
    sprintf("%s -> %s\nthe residue", num3(d$unpinned), num3(d$hypoxia_anchored)),
    sprintf("%s -> %s\n%d sets", num3(d$unpinned), num3(d$hypoxia_anchored),
            d$n_other)))

fig2 <- (pV1 / pV2) +
  patchwork::plot_layout(heights = c(1.55, 1), guides = "collect") +
  patchwork::plot_annotation(
    title = "The hypoxia-anchored variant drawn beside the unconstrained one",
    theme = project_theme(config = FIG_CFG)) &
  theme(legend.position = "bottom")

fig2_table <- dplyr::bind_rows(lapply(BAR_ARMS, function(a) {
  u <- wide_shares(SHARES[SHARES$arm == a & SHARES$variant == VARIANT_MAIN, ])
  h <- wide_shares(SHARES[SHARES$arm == a & SHARES$variant == "hypoxia_anchored", ])
  m <- dplyr::full_join(
    u[, c("ID", "collection", "share_fractional", "share_winner_take_all")],
    h[, c("ID", "collection", "share_fractional", "share_winner_take_all")],
    by = c("ID", "collection"), suffix = c("_unpinned", "_anchored"))
  for (cl in names(m)) if (is.numeric(m[[cl]])) m[[cl]][is.na(m[[cl]])] <- 0
  m$arm <- a
  m$pinned <- m$ID %in% HYP$term
  st <- HYP[HYP$arm == a & !is.na(HYP$term), c("term", "status", "n_arm_hits")]
  m <- dplyr::left_join(m, st, by = c("ID" = "term"))
  m$delta_fractional <- m$share_fractional_anchored - m$share_fractional_unpinned
  m$delta_winner_take_all <- m$share_winner_take_all_anchored -
    m$share_winner_take_all_unpinned
  m[order(-abs(m$delta_fractional)), ]
}))

wt_res <- VAR_ROWS[VAR_ROWS$arm == "WT_heat_up" & VAR_ROWS$row == RESIDUAL, ]
wt_hyp <- VAR_ROWS[VAR_ROWS$arm == "WT_heat_up" & VAR_ROWS$row == "HALLMARK_HYPOXIA", ]

save_overview(
  fig2, STAGE, "arm_composition_variants",
  table = fig2_table,
  finding = sprintf(paste0(
    "Pinning the hypoxia sets takes WT_heat_up from %d to %d categories fractional and %d to %d ",
    "winner-take-all. Both additions carry NA significance, because ",
    "GOBP_RESPONSE_TO_OXYGEN_LEVELS and GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS were tested and ",
    "did not reach the bar. The residue falls by %s, from %s to %s, which is the whole of what ",
    "anchoring recovers. HALLMARK_HYPOXIA gives fractional weight to the pins it shares genes ",
    "with, %s down to %s, while its winner-take-all share holds at %s. ",
    "REACTOME_CELLULAR_RESPONSE_TO_HYPOXIA was tested in every arm and holds none of their genes, ",
    "so it sits at zero under both variants."),
    n_cats("WT_heat_up", VARIANT_MAIN, "fractional"),
    n_cats("WT_heat_up", "hypoxia_anchored", "fractional"),
    n_cats("WT_heat_up", VARIANT_MAIN, "winner_take_all"),
    n_cats("WT_heat_up", "hypoxia_anchored", "winner_take_all"),
    num3(wt_res$fractional_unpinned - wt_res$fractional_anchored),
    num3(wt_res$fractional_unpinned), num3(wt_res$fractional_anchored),
    num3(wt_hyp$fractional_unpinned), num3(wt_hyp$fractional_anchored),
    num3(wt_hyp$winner_take_all_anchored)),
  script = SCRIPT, fn = "variant_panel() / top-level (pV1/pV2)",
  config_kv = CONFIG_KV, input = INPUT,
  how_to_read = paste0(
    "Every row is drawn twice: a hollow grey circle for the `unpinned` variant and a filled ",
    "purple circle for `hypoxia_anchored`, joined by a grey segment. The monospaced text on the ",
    "right of each row gives that column's pair as `unpinned -> anchored` on the first line and ",
    "the set's status in that arm on the second, printed purple wherever that status reads `not ",
    "enriched, NA`. A purple row is a set that entered because the configuration named it, was ",
    "measured, and did not reach the bar: its share is a gene count and supports no enrichment ",
    "claim. Panel A runs to 3% and panel B to 100%, because the pinned ",
    "shares and the residue differ by more than an order of magnitude and one axis would flatten ",
    "one of them. Panel B's `all other selected sets` row aggregates every category anchoring ",
    "leaves alone; the per-category deltas are in the sibling table. Claim tier: descriptive, and ",
    "the purple rows support no enrichment claim at all."),
  config = FIG_CFG, wide = TRUE, height = 15)

## ---------------------------------------------------------------------------
## Figure 3 -- how much of an arm carries a hypoxia annotation
##
## The sets carrying hypoxia hits overlap, so their counts cannot be added. Panel A draws
## the membership gene by gene, which makes the union countable and the sum uncountable.
## Panel B states both numbers. Panel C is the opportunity audit: two of the five
## collections carry no hypoxia set at all and a third pinned set falls below the tested
## size window, and none of those cells is a zero share.
## ---------------------------------------------------------------------------

HYP_ARMS <- HYP %>%
  dplyr::filter(!is.na(.data$term), .data$n_arm_hits > 0) %>%
  dplyr::pull("arm") %>% unique()

memb <- dplyr::bind_rows(lapply(HYP_ARMS, function(a) {
  h <- HYP[HYP$arm == a & !is.na(HYP$term) & HYP$n_arm_hits > 0, ]
  genes <- sort(unique(unlist(strsplit(h$hits, "/", fixed = TRUE))))
  grid <- expand.grid(gene = genes, term = h$term, stringsAsFactors = FALSE)
  grid$arm <- a
  grid$member <- mapply(function(g, tm)
    g %in% strsplit(h$hits[h$term == tm], "/", fixed = TRUE)[[1]],
    grid$gene, grid$term)
  grid$collection <- h$collection[match(grid$term, h$term)]
  grid
}))
k_per_gene <- memb %>% dplyr::filter(.data$member) %>%
  dplyr::count(.data$arm, .data$gene, name = "k")
memb <- dplyr::left_join(memb, k_per_gene, by = c("arm", "gene"))
# Ordering has to be per arm, and one shared discrete scale cannot hold two orders, so
# the axis key carries the arm and the label strips it back off.
memb$gene_key <- paste(memb$arm, memb$gene)
gene_lv <- k_per_gene %>%
  dplyr::arrange(match(.data$arm, HYP_ARMS), -.data$k, .data$gene) %>%
  dplyr::mutate(key = paste(.data$arm, .data$gene))
memb$gene_key <- factor(memb$gene_key, levels = gene_lv$key)

term_lab <- HYP %>%
  dplyr::filter(!is.na(.data$term), .data$n_arm_hits > 0) %>%
  dplyr::mutate(lab = sprintf("%s\n%d of the arm, %s", .data$term, .data$n_arm_hits,
                              .data$status)) %>%
  dplyr::select("arm", "term", "lab")
memb <- dplyr::left_join(memb, term_lab, by = c("arm", "term"))
memb$arm_f <- factor(arm_strip(memb$arm), levels = vapply(HYP_ARMS, arm_strip, character(1)))

pH1 <- ggplot(memb, aes(x = .data$gene_key, y = .data$lab)) +
  geom_tile(fill = "grey96", colour = "white", linewidth = 0.7) +
  geom_tile(data = memb[memb$member, ], aes(fill = .data$collection),
            colour = "white", linewidth = 0.7) +
  facet_wrap(~ .data$arm_f, ncol = 1, scales = "free") +
  scale_fill_manual(values = CAT_COL, name = NULL, limits = c("Hallmark", "GO_BP")) +
  scale_x_discrete(labels = function(x) sub("^\\S+ ", "", x)) +
  labs(title = "A. Which arm gene each hypoxia set holds", x = NULL, y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        axis.text.y = element_text(lineheight = 1.3),
        panel.grid = element_blank(), legend.position = "bottom")

union_tbl <- dplyr::bind_rows(lapply(HYP_ARMS, function(a) {
  h <- HYP[HYP$arm == a & !is.na(HYP$term) & HYP$n_arm_hits > 0, ]
  u <- length(unique(unlist(strsplit(h$hits, "/", fixed = TRUE))))
  data.frame(arm = a,
             quantity = c("the three set counts added up", "genes in at least one set"),
             n = c(sum(h$n_arm_hits), u),
             n_arm = unname(arm_n[[a]]))
}))
union_tbl$frac <- union_tbl$n / union_tbl$n_arm
union_tbl$quantity <- factor(union_tbl$quantity,
                             levels = c("the three set counts added up",
                                        "genes in at least one set"))
union_tbl$arm_f <- factor(arm_strip(union_tbl$arm),
                          levels = vapply(HYP_ARMS, arm_strip, character(1)))

pH2 <- ggplot(union_tbl, aes(x = .data$n, y = .data$quantity, fill = .data$quantity)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%d genes, %s of the arm", .data$n, pct(.data$frac))),
            hjust = -0.06, size = LAB * 0.95, fontface = "bold", colour = "grey15") +
  facet_wrap(~ .data$arm_f, ncol = 1) +
  scale_fill_manual(values = c("the three set counts added up" = "grey72",
                               "genes in at least one set" = OI$blue), guide = "none") +
  scale_x_continuous(limits = c(0, max(union_tbl$n) * 1.75),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "B. Distinct genes, against the sum of the three set counts",
       x = "genes", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank())

status_tbl <- HYP %>%
  dplyr::mutate(
    row = ifelse(is.na(.data$term), sprintf("%s (no hypoxia set in this collection)",
                                            .data$collection), .data$term),
    cell = dplyr::case_when(
      is.na(.data$term) ~ "no set in the collection",
      grepl("^below|^above", .data$status) ~ "outside the tested size window",
      .data$status == "tested, enriched" ~ "tested, enriched",
      TRUE ~ "tested, not enriched"),
    txt = dplyr::case_when(
      is.na(.data$term) ~ "no set to test",
      grepl("^below", .data$status) ~ sprintf("%d genes, below the floor",
                                              .data$size_in_background),
      grepl("^above", .data$status) ~ sprintf("%d genes, above the ceiling",
                                              .data$size_in_background),
      TRUE ~ sprintf("%d arm genes", .data$n_arm_hits)))
status_tbl$row <- factor(status_tbl$row,
                         levels = rev(unique(status_tbl$row[order(
                           match(status_tbl$collection, COLLECTIONS),
                           is.na(status_tbl$term), status_tbl$term)])))
status_tbl$arm <- factor(status_tbl$arm, levels = ARM_ORDER)

pH3 <- ggplot(status_tbl, aes(x = .data$arm, y = .data$row, fill = .data$cell)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = .data$txt), size = LAB * 0.78, colour = "grey12",
            fontface = "bold") +
  scale_fill_manual(values = c("tested, enriched" = OI$orange,
                               "tested, not enriched" = OI$sky_blue,
                               "outside the tested size window" = "grey78",
                               "no set in the collection" = "grey90"), name = NULL) +
  scale_y_discrete(labels = function(x) wrap_id(x, 40L)) +
  labs(title = "C. What each collection had available to find", x = NULL, y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        axis.text.y = element_text(lineheight = 0.95),
        axis.text.x = element_text(angle = 18, hjust = 1))

fig3 <- (pH1 / pH2 / pH3) +
  patchwork::plot_layout(heights = c(1.15, 0.72, 1.3)) +
  patchwork::plot_annotation(
    title = "The hypoxia sets holding arm genes, their overlap, and their union",
    theme = project_theme(config = FIG_CFG))

# Two grains live in this figure: one row for each pinned set (panels A and C) and one row per
# pair of sets (the overlap the finding quotes). Joining them into one flat frame duplicates
# HALLMARK_HYPOXIA once per pair and makes sum(n_arm_hits) read 51 where the per-set sum is
# 33. So `row_grain` is explicit and no column carries a meaning that depends on it: the
# per-set columns are NA on pair rows and the pair columns are NA on set rows.
# union_tbl only covers the arms that hold hypoxia genes. For the other two the union is a
# MEASURED zero: their pinned sets were tested and hold none of their genes. NA would read
# as "not computed", which is the state reserved for the collections carrying no set at all.
union_per_arm <- union_tbl %>%
  dplyr::filter(.data$quantity == "genes in at least one set") %>%
  dplyr::transmute(arm = .data$arm, n_union = .data$n, union_share_of_arm = .data$frac)
union_per_arm <- dplyr::bind_rows(
  union_per_arm,
  data.frame(arm = setdiff(ARM_ORDER, union_per_arm$arm),
             n_union = 0L, union_share_of_arm = 0))

fig3_set_rows <- HYP %>%
  dplyr::transmute(row_grain = "set", arm = .data$arm, collection = .data$collection,
                   term = .data$term, status = .data$status,
                   size_in_background = .data$size_in_background,
                   n_arm_hits = .data$n_arm_hits, hits = .data$hits,
                   p.adjust = .data$p.adjust, p_matched = .data$p_matched,
                   term_a = NA_character_, term_b = NA_character_,
                   n_a = NA_integer_, n_b = NA_integer_,
                   n_shared = NA_integer_, jaccard = NA_real_)

fig3_pair_rows <- HYP_OVL %>%
  dplyr::transmute(row_grain = "set_pair", arm = .data$arm,
                   collection = NA_character_, term = NA_character_,
                   status = NA_character_, size_in_background = NA_integer_,
                   n_arm_hits = NA_integer_, hits = NA_character_,
                   p.adjust = NA_real_, p_matched = NA_real_,
                   term_a = .data$term_a, term_b = .data$term_b,
                   n_a = .data$n_a, n_b = .data$n_b,
                   n_shared = .data$n_shared, jaccard = .data$jaccard)

fig3_table <- dplyr::bind_rows(fig3_set_rows, fig3_pair_rows) %>%
  dplyr::left_join(union_per_arm, by = "arm") %>%
  dplyr::arrange(match(.data$arm, ARM_ORDER), .data$row_grain)

# The per-set column must sum to the grey bar of panel B, and the union must not move.
stopifnot(
  sum(fig3_table$n_arm_hits[fig3_table$row_grain == "set" &
                              fig3_table$arm == "WT_heat_up"], na.rm = TRUE) ==
    union_tbl$n[union_tbl$arm == "WT_heat_up" &
                  union_tbl$quantity == "the three set counts added up"],
  length(unique(fig3_table$n_union[fig3_table$arm == "WT_heat_up"])) == 1L)

save_overview(
  fig3, STAGE, "arm_hypoxia_sources",
  table = fig3_table,
  finding = sprintf(paste0(
    "Three sets hold hypoxia-annotated genes of WT_heat_up: HALLMARK_HYPOXIA with 18, ",
    "GOBP_RESPONSE_TO_OXYGEN_LEVELS with 11 and GOBP_CELLULAR_RESPONSE_TO_OXYGEN_LEVELS with 4. ",
    # The arm size is read from the data rather than written here: it moved from 199 to 202
    # when the ortholog step stopped counting stale mouse symbols as having no human ortholog.
    "Those counts add to 33 and the union is %d genes, %s of the %d. Hallmark and GO BP share 3 ",
    "genes, a Jaccard of %s, so the two collections read largely different genes under the same ",
    "word. Only HALLMARK_HYPOXIA reached significance. GO MF and KEGG carry no hypoxia set at ",
    "all, and a fourth pinned GO BP set has 7 background genes against the 10-gene floor, so ",
    "three cells of panel C record an absence of opportunity."),
    union_tbl$n[union_tbl$arm == "WT_heat_up" &
                  union_tbl$quantity == "genes in at least one set"],
    pct(union_tbl$frac[union_tbl$arm == "WT_heat_up" &
                         union_tbl$quantity == "genes in at least one set"]),
    union_tbl$n_arm[union_tbl$arm == "WT_heat_up" &
                      union_tbl$quantity == "genes in at least one set"],
    num3(HYP_OVL$jaccard[HYP_OVL$arm == "WT_heat_up" &
                           HYP_OVL$term_a == "HALLMARK_HYPOXIA" &
                           HYP_OVL$term_b == "GOBP_RESPONSE_TO_OXYGEN_LEVELS"])),
  script = SCRIPT, fn = "top-level (pH1/pH2/pH3)",
  config_kv = CONFIG_KV, input = INPUT,
  how_to_read = paste0(
    "Panel A is a membership grid: one column per arm gene, one row per set, a filled cell where ",
    "the set holds the gene. Columns are ordered by how many sets hold them, so the shared genes ",
    "stand on the left and the union is the width of the panel. Panel B puts the two candidate ",
    "numbers side by side, the sum of the three set counts in grey and the count of distinct ",
    "genes in blue. The blue one answers how much of the arm carries a hypoxia annotation. Panel ",
    "C is the opportunity audit, one cell per set per arm: orange reached significance, blue was ",
    "tested and did not, and the two greys are a set below the tested size window and a ",
    "collection with no hypoxia set to test. A grey cell is an absence of opportunity and is ",
    "drawn as a labelled cell rather than as a bar of height zero. Panels A and B cover the two ",
    "arms holding hypoxia genes, and the two interaction arms hold none, which panel C records. ",
    "Claim tier: descriptive gene-content accounting. It measures annotation overlap and says ",
    "nothing about causal structure."),
  config = FIG_CFG, wide = TRUE, height = 18)

## ---------------------------------------------------------------------------
## Figure 4 -- the residue, split by reason
##
## The three classes in panel A partition the unclaimed genes exactly. The count in panel
## B is a SUBSET of panel A's `annotated in a testable set` band, so it is drawn against
## the total rather than stacked into it.
## ---------------------------------------------------------------------------

rem_u <- REMAINDER[REMAINDER$variant == VARIANT_MAIN, ]
rem_u$arm_f <- factor(arm_strip(rem_u$arm),
                      levels = rev(vapply(ARM_ORDER, arm_strip, character(1))))
REASON_LV <- c("annotated in a testable set",
               "annotated only at an untestable set size (0 for every arm here)",
               "no annotation in any of the five collections")
rem_long <- rem_u %>%
  dplyr::select("arm_f", "n_unclaimed", "n_arm_in_background",
                a = "n_unclaimed_annotated_testable",
                b = "n_unclaimed_annotated_only_untestable_size",
                c = "n_unclaimed_no_annotation_at_all") %>%
  tidyr::pivot_longer(cols = c("a", "b", "c"), names_to = "k", values_to = "n") %>%
  dplyr::mutate(reason = factor(REASON_LV[match(.data$k, c("a", "b", "c"))],
                                levels = rev(REASON_LV)))

pR1 <- ggplot(rem_long, aes(x = .data$n, y = .data$arm_f, fill = .data$reason)) +
  geom_col(width = 0.6, colour = "white", linewidth = 0.5) +
  geom_text(data = rem_u,
            aes(x = .data$n_unclaimed, y = .data$arm_f,
                label = sprintf("  %d of %d, %s", .data$n_unclaimed,
                                .data$n_arm_in_background,
                                pct(.data$n_unclaimed / .data$n_arm_in_background))),
            inherit.aes = FALSE, hjust = 0, size = LAB * 0.95, fontface = "bold",
            colour = "grey15") +
  scale_fill_manual(values = setNames(c(OI$sky_blue, OI$yellow, "grey45"), REASON_LV),
                    name = NULL, breaks = REASON_LV) +
  scale_x_continuous(limits = c(0, max(rem_u$n_unclaimed) * 1.45),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "A. Genes no selected set claims, split by why", x = "genes", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom") +
  guides(fill = guide_legend(ncol = 1))

pR2 <- ggplot(rem_u, aes(y = .data$arm_f)) +
  geom_col(aes(x = .data$n_unclaimed), width = 0.6, fill = "grey88") +
  geom_col(aes(x = .data$n_unclaimed_in_an_enriched_set), width = 0.34,
           fill = OI$vermillion) +
  geom_text(aes(x = .data$n_unclaimed_in_an_enriched_set,
                label = sprintf("  %d of %d", .data$n_unclaimed_in_an_enriched_set,
                                .data$n_unclaimed)),
            hjust = 0, size = LAB * 0.95, fontface = "bold", colour = "grey15") +
  scale_x_continuous(limits = c(0, max(rem_u$n_unclaimed) * 1.35),
                     expand = expansion(mult = c(0, 0))) +
  labs(title = "B. Of those, how many sit in a set that enriched and lost selection",
       x = "genes", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank())

fig4 <- (pR1 / pR2) +
  patchwork::plot_layout(heights = c(1, 0.82)) +
  patchwork::plot_annotation(
    title = "The residue of each arm, split by why a gene lands in it",
    theme = project_theme(config = FIG_CFG))

save_overview(
  fig4, STAGE, "arm_remainder",
  table = rem_u %>% dplyr::select(-"arm_f"),
  finding = sprintf(paste0(
    "WT_heat_up leaves %d of %d genes in no selected set, %s of the arm. Of those, %d carry an ",
    "annotation in a set that was testable, %d carry none in any of the five collections, and %d ",
    "of the %d sit in a set that did reach significance and lost its place to the redundancy ",
    "prune or the ten-per-collection cap. KO_heat_up leaves %d of %d and Interaction_fdrOnly_up ",
    "%d of %d. Interaction_up leaves 0 of 7, which is what a 7-gene list with one dominant ",
    "annotation looks like."),
    rem_u$n_unclaimed[rem_u$arm == "WT_heat_up"], unname(arm_n[["WT_heat_up"]]),
    pct(rem_u$n_unclaimed[rem_u$arm == "WT_heat_up"] / arm_n[["WT_heat_up"]]),
    rem_u$n_unclaimed_annotated_testable[rem_u$arm == "WT_heat_up"],
    rem_u$n_unclaimed_no_annotation_at_all[rem_u$arm == "WT_heat_up"],
    rem_u$n_unclaimed_in_an_enriched_set[rem_u$arm == "WT_heat_up"],
    rem_u$n_unclaimed[rem_u$arm == "WT_heat_up"],
    rem_u$n_unclaimed[rem_u$arm == "KO_heat_up"], unname(arm_n[["KO_heat_up"]]),
    rem_u$n_unclaimed[rem_u$arm == "Interaction_fdrOnly_up"],
    unname(arm_n[["Interaction_fdrOnly_up"]])),
  script = SCRIPT, fn = "top-level (pR1/pR2)",
  config_kv = CONFIG_KV, input = INPUT,
  how_to_read = paste0(
    "Panel A stacks three mutually exclusive classes, so the bar length is the whole residue and ",
    "the text beside it repeats that as a count and a share of the arm. The yellow class is in ",
    "the key because it completes the partition and it is zero for all four arms, which is why no ",
    "yellow band is drawn. Panel B asks a different question of the same total: grey is the ",
    "residue again and red counts the genes inside it that do sit in a set that reached ",
    "significance and then lost its place to the redundancy prune or the ten-per-collection cap. ",
    "The red count is a subset of panel A's blue band, so the two panels are read against each ",
    "other and never added. An arm with a residue of zero carries its zero as a label rather than ",
    "an empty row. The ten WT_heat_up genes with no annotation anywhere are named in ",
    "composition_remainder.csv. Claim tier: descriptive coverage accounting."),
  config = FIG_CFG, wide = TRUE, height = 11.5)

## ---------------------------------------------------------------------------
## Figure 5 -- the null frame
##
## A gene set drawn at random but matched on annotation depth already reaches significance
## in GO BP a great many times. A term count printed without that frame reads as a result,
## so both the term count and the gene coverage are drawn against the null they have to
## beat.
## ---------------------------------------------------------------------------

NULL_KEY <- sprintf("depth-matched null, median to 95th percentile of %s draws",
                    format(as.integer(PROV[["n_null"]]), big.mark = ","))
OBS_KEY  <- "observed for the arm"

# `%in%` rather than `==`: a base-R logical subscript propagates NA into an all-NA row,
# which silently lengthens every lookup made against this frame.
# `%in%` rather than `==` on purpose: `null` carries a third level for a cell the engine could
# not run, and if that column ever held NA instead, `df[NA, ]` returns an all-NA ROW rather than
# dropping it. That row would ship in the sibling table while the geoms quietly drop it.
nd <- NULLS[NULLS$null %in% "depth_matched", , drop = FALSE]
NOT_TESTED <- NULLS[!is.na(NULLS$not_tested_reason),
                    c("arm", "collection", "not_tested_reason"), drop = FALSE]
null_long <- dplyr::bind_rows(
  nd %>% dplyr::transmute(arm = .data$arm, collection = .data$collection,
                          measure = "sets reaching significance",
                          med = .data$median_n_signif, q95 = .data$q95_n_signif,
                          obs = .data$n_observed_signif),
  nd %>% dplyr::transmute(arm = .data$arm, collection = .data$collection,
                          measure = "arm genes covered by a significant set",
                          med = .data$median_covered, q95 = .data$q95_covered,
                          obs = .data$covered_obs))
# An arm x collection the engine could not test must not appear as an empty axis row.
full <- expand.grid(arm = ARM_ORDER, collection = COLLECTIONS,
                    measure = unique(null_long$measure), stringsAsFactors = FALSE)
null_long <- dplyr::left_join(full, null_long, by = c("arm", "collection", "measure"))
null_long <- dplyr::left_join(null_long, NOT_TESTED, by = c("arm", "collection"))
null_long$untested <- is.na(null_long$obs)
# The reason a cell was not tested is read from the compute stage and never inferred here. A
# blank cell carrying no recorded reason is a gap to report, not one to label from the source.
stopifnot(!any(null_long$untested & is.na(null_long$not_tested_reason)))
null_long$arm_f <- factor(arm_strip2(null_long$arm),
                          levels = vapply(ARM_ORDER, arm_strip2, character(1)))
null_long$collection <- factor(null_long$collection, levels = rev(COLLECTIONS))
null_long$measure <- factor(null_long$measure,
                            levels = c("sets reaching significance",
                                       "arm genes covered by a significant set"))

fig5 <- ggplot(null_long, aes(y = .data$collection)) +
  geom_segment(aes(x = .data$med, xend = .data$q95, y = .data$collection,
                   yend = .data$collection, colour = NULL_KEY),
               linewidth = 2.6, lineend = "round", na.rm = TRUE) +
  geom_point(aes(x = .data$med), shape = 124, size = PT * 2.2, colour = "grey25",
             na.rm = TRUE) +
  geom_point(aes(x = .data$q95), shape = 124, size = PT * 2.2, colour = "grey25",
             na.rm = TRUE) +
  geom_point(aes(x = .data$obs, fill = .data$collection, shape = OBS_KEY),
             size = PT * 1.9, colour = "grey10", stroke = 0.9, na.rm = TRUE) +
  geom_text(aes(x = .data$obs, label = ifelse(is.na(.data$obs), "", .data$obs)),
            vjust = -1.2, size = LAB * 0.82, fontface = "bold", colour = "grey15",
            na.rm = TRUE) +
  geom_point(data = null_long[null_long$untested, ], aes(x = 0), shape = 4,
             size = PT * 1.5, colour = "grey35", stroke = 1.1) +
  geom_text(data = null_long[null_long$untested, ],
            aes(x = 0, label = paste0("    no test run: ", .data$not_tested_reason)),
            hjust = 0, size = LAB * 0.8, colour = "grey35") +
  facet_grid(.data$arm_f ~ .data$measure, scales = "free_x", switch = "y") +
  scale_fill_manual(values = CAT_COL, guide = "none", drop = FALSE) +
  scale_colour_manual(values = setNames("grey62", NULL_KEY), name = NULL) +
  scale_shape_manual(values = setNames(21L, OBS_KEY), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.16))) +
  labs(title = "Sets found and genes covered, against a depth-matched random gene set",
       x = "count", y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom",
        strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, lineheight = 1.05)) +
  guides(shape = guide_legend(override.aes = list(fill = "grey45", size = PT * 1.9)),
         colour = guide_legend(override.aes = list(linewidth = 2.6)))

save_overview(
  fig5, STAGE, "arm_composition_null",
  table = NULLS %>%
    dplyr::filter(.data$null %in% c("depth_matched", "not_tested")) %>%
    dplyr::select("arm", "collection", "null", "not_tested_reason", "n_replicates",
                  "n_query_annotated", "n_terms_testable",
                  "median_n_signif", "q95_n_signif", "n_observed_signif",
                  "median_covered", "q95_covered", "covered_obs",
                  "p_covered", "max_genes_borrowed"),
  finding = sprintf(paste0(
    "For WT_heat_up a gene set drawn at random and matched on annotation depth reaches ",
    "significance in a median of %s GO BP sets and covers a median of %s of the drawn genes, ",
    "against %s sets and %s genes observed. In Reactome, GO MF and KEGG the same random draw ",
    "reaches a median of zero significant sets and covers a median of zero genes, against %s, %s ",
    "and %s genes observed. Annotation depth does most of the work in GO BP, and the three ",
    "smaller collections separate an observed result from a depth artefact far more sharply. ",
    "Every coverage p is %s, the floor of a %s-replicate permutation."),
    nd$median_n_signif[nd$arm == "WT_heat_up" & nd$collection == "GO_BP"],
    nd$median_covered[nd$arm == "WT_heat_up" & nd$collection == "GO_BP"],
    nd$n_observed_signif[nd$arm == "WT_heat_up" & nd$collection == "GO_BP"],
    nd$covered_obs[nd$arm == "WT_heat_up" & nd$collection == "GO_BP"],
    nd$covered_obs[nd$arm == "WT_heat_up" & nd$collection == "Reactome"],
    nd$covered_obs[nd$arm == "WT_heat_up" & nd$collection == "GO_MF"],
    nd$covered_obs[nd$arm == "WT_heat_up" & nd$collection == "KEGG"],
    format(signif(min(nd$p_covered), 2), scientific = FALSE),
    format(as.integer(PROV[["n_null"]]), big.mark = ",")),
  script = SCRIPT, fn = "top-level (fig5)",
  config_kv = CONFIG_KV, input = INPUT,
  how_to_read = paste0(
    "Each row is one collection. The thick grey bar spans the median to the 95th percentile of ",
    "the depth-matched null, with a tick at each end, and the filled circle with the number above ",
    "it is what the arm achieved. A row where the grey bar collapses to a point sits at a null ",
    "median and 95th percentile of zero, which is the strongest frame a collection can offer. ",
    "Left column counts sets reaching significance, right column counts genes of the arm covered ",
    "by at least one such set. A cross at zero with a label means the engine had under two of ",
    "the arm's genes annotated in that collection and ran no test there, which is the case for ",
    "Interaction_up in KEGG. The uniform null is in composition_null_summary.csv and is left off ",
    "the panel, because the gate that produced these arms admits well-studied genes and only the ",
    "depth-matched draw carries the comparison. Claim tier: this is the calibration frame for ",
    "every count elsewhere in this directory, and by itself supports no biological claim."),
  config = FIG_CFG, wide = TRUE, height = 16)

## ---------------------------------------------------------------------------
## Figure 6 -- Interaction_up, gene by gene
##
## 7 genes. No composition bar is drawn for this arm anywhere in the stage.
## ---------------------------------------------------------------------------

ia <- ASSIGN %>%
  dplyr::filter(.data$arm == "Interaction_up", .data$variant == VARIANT_MAIN,
                .data$ID != RESIDUAL)
ia_grid <- expand.grid(gene = sort(unique(ia$gene)), ID = sort(unique(ia$ID)),
                       stringsAsFactors = FALSE)
ia_grid <- dplyr::left_join(ia_grid, ia[, c("gene", "ID", "is_winner")],
                            by = c("gene", "ID"))
ia_grid$member <- !is.na(ia_grid$is_winner)
ia_grid$collection <- ia$collection[match(ia_grid$ID, ia$ID)]
k_lab <- ia %>% dplyr::distinct(.data$gene, .data$n_sets_for_gene) %>%
  dplyr::mutate(lab = sprintf("%s\n%d sets", .data$gene, .data$n_sets_for_gene))
ia_grid$gene_lab <- k_lab$lab[match(ia_grid$gene, k_lab$gene)]
set_ord <- ia %>% dplyr::count(.data$ID, name = "n") %>% dplyr::arrange(.data$n, .data$ID)
ia_grid$ID <- factor(ia_grid$ID, levels = set_ord$ID)
ia_grid$gene_lab <- factor(ia_grid$gene_lab,
                           levels = k_lab$lab[order(-k_lab$n_sets_for_gene, k_lab$gene)])

fig6 <- ggplot(ia_grid, aes(x = .data$gene_lab, y = .data$ID)) +
  geom_tile(fill = "grey96", colour = "white", linewidth = 0.8) +
  geom_tile(data = ia_grid[ia_grid$member, ], aes(fill = .data$collection),
            colour = "white", linewidth = 0.8) +
  geom_point(data = ia_grid[which(ia_grid$is_winner), ],
             aes(shape = "winner_take_all takes the gene here"),
             size = PT * 1.7, fill = "white", colour = "grey10", stroke = 1.2) +
  scale_fill_manual(values = CAT_COL, name = NULL,
                    limits = intersect(CAT_LEVELS, unique(ia$collection))) +
  scale_shape_manual(values = c("winner_take_all takes the gene here" = 21L), name = NULL) +
  scale_y_discrete(labels = function(x) wrap_id(x, 44L)) +
  labs(title = "Interaction_up, all 7 genes and every selected set holding one",
       x = NULL, y = NULL) +
  project_theme(config = FIG_CFG) +
  theme(panel.grid = element_blank(), legend.position = "bottom",
        legend.box = "vertical",
        axis.text.y = element_text(lineheight = 0.92))

ia_table <- ia %>%
  dplyr::select("arm", "gene", "ID", "collection", "n_sets_for_gene",
                "weight_fractional", "winner_id", "is_winner", "p_matched", "p.adjust") %>%
  dplyr::arrange(.data$gene, .data$ID)

save_overview(
  fig6, STAGE, "interaction_up_genes",
  table = ia_table,
  finding = local({
    # Computed, not named. The winner-take-all split for this arm is decided by the same
    # tie-break as every other arm, so a hand-written split is one re-run from being wrong.
    win <- ia[which(ia$is_winner), ]
    tab <- sort(table(win$ID), decreasing = TRUE)
    parts <- sprintf("%d to %s", as.integer(tab), names(tab))
    sprintf(paste0(
      "Interaction_up holds %d genes, %s, and each of them lands in at least two selected sets, ",
      "so its residue is 0 of %d. Under `winner_take_all` they go %s, which is why that ",
      "accounting has %d categories for this arm against %d fractional. Seven genes over five ",
      "collections is a gene list, and this stage draws it as one."),
      length(unique(ia$gene)), paste(sort(unique(ia$gene)), collapse = ", "),
      length(unique(ia$gene)),
      paste(parts, collapse = " and "),
      n_cats("Interaction_up", VARIANT_MAIN, "winner_take_all"),
      n_cats("Interaction_up", VARIANT_MAIN, "fractional"))
  }),
  script = SCRIPT, fn = "top-level (fig6)",
  config_kv = CONFIG_KV, input = INPUT,
  how_to_read = paste0(
    "One column per gene, one row per selected set, a filled cell where the set holds the gene, ",
    "and fill giving the collection. The white ring marks the set that gene goes to under ",
    "`winner_take_all`. The column header repeats how many selected sets hold the gene, which is ",
    "the k that `fractional` divides its weight by. Rows are ordered by how many of the seven ",
    "genes they hold. Claim tier: a gene list drawn in full, hypothesis-generating, and no share ",
    "statistic is computed from 7 genes anywhere in this directory."),
  config = FIG_CFG, wide = FALSE, width = 11.5, height = 10)

message("27_arm_composition_viz: VIZ DONE.")
