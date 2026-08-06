#!/usr/bin/env Rscript
# 27_arm_composition_viz.R -- VIZ ONLY (no computing)
# =============================================================================
# Draws the two overview figures for 16_arm_composition from the frozen object written by
# 26_arm_composition.R. Every number on every face is read from that object.
#
#   arm_composition           the composition itself, per arm, both accountings
#   arm_composition_variants  `unpinned` beside `hypoxia_anchored`, pins marked
#
# Two accountings are always drawn together and both are kept. `fractional` gives a gene in
# k selected sets 1/k to each. `winner_take_all` gives each gene to its single best selected
# set. Both sum to 1.0 over the arm, so where the two disagree the disagreement is the
# quantity on the page.
#
# Guards, and the silent failures they exist for. Each one fired here for real, and each was
# silent on its own. Read this before loosening any of them.
#
#   VARIANTS + the per-frame check     A factor level renamed upstream. A filter on the absent
#                                      level returns ZERO ROWS silently, so the figure renders
#                                      as a convincing empty panel and a reader filtering the
#                                      published table gets nothing back. Frames this script
#                                      draws from hard-stop; the rest warn loudly, keeping a
#                                      correct figure unblocked over a table this script never
#                                      touches.
#   dup_free()                         A duplicated key. Every caption number is pulled by a
#                                      logical mask expecting one row; a duplicate returns
#                                      length 2, sprintf recycles it silently, and it surfaces
#                                      far downstream as a coercion error inside the caption
#                                      writer.
#   nz()                               An empty per-arm subset stops here, before ggplot.
#   INPUT carries the object md5       A partial re-run once left figures from a new object
#                                      beside one from the old, and a README from the new,
#                                      silently, in a directory that looked complete. Only
#                                      timestamps disagreed. The stamp makes the check one
#                                      `md5sum` against the caption.
#   `%in%` on the null level           `df[NA, ]` returns an all-NA ROW, which then ships in a
#                                      sibling table while the geoms drop it from the panel.
#   row_order()                        Exact ties in the sort key, broken by input order, make
#                                      vertical position imply a rank the data leaves open.
#                                      WT_heat_up carries two such ties.
#   Captions derive from the table    Two hand-written captions went stale within an hour.
#                                      The leader is computed, and a tied group is reported as
#                                      a group.
#
# A longer write-up of the panel design sits with the other reasoning notes, in a directory
# this repository leaves untracked, so this header is the tracked copy of anything
# load-bearing.
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
REMAINDER  <- obj$remainder
HYP        <- obj$hypoxia
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
# filter on a renamed level returns ZERO ROWS silently, and a figure drawn from
# an empty frame renders as an empty panel with a clean exit code. Stop instead.
VARIANT_MAIN <- "unpinned"
VARIANT_PIN  <- "hypoxia_anchored"
VARIANTS     <- c(VARIANT_MAIN, VARIANT_PIN)
# Check EVERY frame in the object that carries a `variant` column, past the ones drawn
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
# before it reaches ggplot.
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
dup_free(HYP, c("arm", "collection", "term"), "hypoxia")
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
# that stand for book-keeping (the residue, and the roll-up of the sets below the drawn cap)
# take greys. KEGG takes the purple, over the second
# orange: Hallmark and KEGG bands sit adjacent in the stacked panel and two oranges
# there are one glyph.
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
# already carries, over a character count.
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
## The set-level marks are paired points, over bars, because 12 of WT_heat_up's 40
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
       subtitle = sprintf(paste0("Interaction_up, %d genes, carries no share bar anywhere: ",
                                 "a five-collection decomposition of %d genes is not a composition."),
                          unname(arm_n[["Interaction_up"]]),
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
# identifier. It decides position alone.
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
# off the data. A hard-coded leader outlives the re-run that changes it.
wt_named <- SET_ROWS[SET_ROWS$arm == "WT_heat_up" &
                       !SET_ROWS$ID %in% c(RESIDUAL, ROLLUP_KEY), , drop = FALSE]
lead_frac <- wt_named[which.max(wt_named$share_fractional), , drop = FALSE]
top_wta   <- max(wt_named$share_winner_take_all)
lead_wta  <- wt_named[abs(wt_named$share_winner_take_all - top_wta) < 1e-12, , drop = FALSE]
lead_wta  <- lead_wta[order(lead_wta$ID), , drop = FALSE]
# Equal winner-take-all share means equal gene count, because that accounting is a gene count
# over the arm size. Assert it before quoting one number for the group.
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
# is reported by figure 3 panel C as an absence of opportunity, distinct from a zero here.
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
# for the pivot, downstream of the compute stage's vocabulary. The variant guard at the top is what
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

message("27_arm_composition_viz: VIZ DONE.")
