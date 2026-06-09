#!/usr/bin/env Rscript
# =============================================================================
# 03c_hif_program_attribution.R - PHASE 3 (stage 04_tf) COMPUTE: fig3l
#   "HIF program" attribution by curated biological MODULE
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# Role:    COMPUTE half of the "normalize-then-visualize" split for the
#          attribution CENTREPIECE (fig3l). Contains NO ggplot()/ggsave() and
#          NO new statistics whatsoever -- it JOINS the already-decomposed
#          Hif1a regulon (fig3g_target_decomposition_data.csv: source/target/
#          mor/t_wt/contrib/class) to a curated, literature-grounded
#          module/isoform lookup. t_wt and contrib are COPIED verbatim from the
#          fig3g table -- never recomputed. The whole point of the figure is to
#          re-BUCKET the SAME signed contributions, not to recompute them.
#
# Inputs (no recomputation):
#   - 03_results/04_tf/tables/fig3g_target_decomposition_data.csv
#       (353 Hif1a regulon members with signed contribution = sign(mor)*t_wt)
#
# Outputs (03_results/04_tf/tables/):
#   - fig3l_hif_attribution_data.csv   (per-member: module/isoform/direction)
#   - fig3l_module_summary.csv         (per-module aggregates; no statistics)
#
# Object name (HARD GUARDRAIL, design-spec sec 0): the thing this figure
#   describes is "a heat-induced glycolytic/stress program partially
#   overlapping HIF targets" -- NEVER "the HIF program" / "HIF1a the TF".
#
# Method: this is a JOIN + simple aggregation only. No lmFit/eBayes/run_ulm/
#   run_mlm/p.adjust/prcomp anywhere. Figures live in the _viz.R partner.
#
# -----------------------------------------------------------------------------
# THE 5 MODULE BUCKETS (single source of truth; defined ONCE below).
# Member genes are verbatim from the design spec sec 2 (fig3l block); the
# literature provenance for the curation:
#   - Heat-shock / stress axis at fever-range 39 C in T cells:
#       Gabai & Calderwood 2007 (PMID 18056375).
#   - HIF1a-selective hypoxic-survival / mitophagy / acidosis core:
#       Hu 2003 (HIF1 vs HIF2 target specificity, PMID 14645546);
#       Keith, Johnson & Simon 2011 (HIF1a vs HIF2a, Nat Rev Cancer);
#       Bellot 2009 (BNIP3/BNIP3L mitophagy);
#       Kaluz, Kaluzova & Stanbridge 2009 (CA9/Car9 regulation, PMID 19344680).
#   - Egln3 / PHD3 negative-feedback (HIF-induced brake):
#       del Peso 2003; Pescador 2005 (PHD3 is a HIF target -> feedback).
#   - Shared HIF1/HIF2 angiogenic / glucose-transport targets (least
#       isoform-diagnostic): Hu 2003; Keith 2011.
# =============================================================================

source("02_analysis/config/config.R")
load_packages(extra = c("readr"))   # dplyr + ggplot2 pulled in; only dplyr used

TBL_DIR <- stage_dir("04_tf", "tables")

# -----------------------------------------------------------------------------
# MODULE LOOKUP -- single source of truth. Each entry: members + isoform tag.
# Facets are ordered (in viz) heatshock_stress -> shared_angio_glucose ->
# autoreg_feedback -> hif1a_hypoxic_core so the punchline (hypoxic core left of
# zero) lands at the bottom. Members NOT listed here fall to other_unclassified.
# -----------------------------------------------------------------------------
MODULE_DEFS <- list(
  heatshock_stress = list(
    # the contamination currently hidden inside fig3g's "other" class
    genes  = c("Hspa1a", "Timp1", "Sdc1", "Cdkn1a", "Serpine1", "Eno2", "Spp1"),
    isoform = "non-HIF (stress)"
  ),
  shared_angio_glucose = list(
    # shared HIF1/HIF2 angiogenic + glucose transport -> least HIF-diagnostic
    genes  = c("Vegfa", "Slc2a1"),
    isoform = "shared HIF1/HIF2"
  ),
  autoreg_feedback = list(
    # Egln3/PHD3 -- HIF-induced negative-feedback brake
    genes  = c("Egln3"),
    isoform = "HIF1a-preferential-feedback"
  ),
  hif1a_hypoxic_core = list(
    # HIF1a-selective survival / mitophagy / acidosis -- the diagnostic core
    genes  = c("Pdk1", "Bnip3", "Bnip3l", "Car9"),
    isoform = "HIF1a-selective"
  )
)

# Flatten the lookup to a per-gene tibble (module + isoform_attribution).
module_lookup <- do.call(rbind, lapply(names(MODULE_DEFS), function(m) {
  data.frame(
    target              = MODULE_DEFS[[m]]$genes,
    module              = m,
    isoform_attribution = MODULE_DEFS[[m]]$isoform,
    stringsAsFactors    = FALSE
  )
}))

# -----------------------------------------------------------------------------
# JOIN fig3g decomposition -> module lookup. COPY t_wt + contrib (no recompute).
# -----------------------------------------------------------------------------
fig3g <- read.csv(file.path(TBL_DIR, "fig3g_target_decomposition_data.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("source", "target", "mor", "t_wt", "contrib", "class") %in% names(fig3g)))

# Sanity: every curated member must exist exactly once in the fig3g table.
missing <- setdiff(module_lookup$target, fig3g$target)
if (length(missing)) {
  stop(sprintf("Curated module genes absent from fig3g table: %s",
               paste(missing, collapse = ", ")))
}

attr_df <- merge(fig3g[, c("source", "target", "t_wt", "contrib")],
                 module_lookup, by = "target", all.x = TRUE)

# Members not in any curated bucket -> the faint context bucket.
attr_df$module[is.na(attr_df$module)]                           <- "other_unclassified"
attr_df$isoform_attribution[is.na(attr_df$isoform_attribution)] <- "unclassified"

# direction = sign(contrib): Up (>=0) vs Down (<0). contrib copied, not recomputed.
attr_df$direction <- ifelse(attr_df$contrib >= 0, "Up", "Down")

# Final column order per the design-spec schema.
attr_df <- attr_df[, c("source", "target", "t_wt", "contrib",
                       "module", "isoform_attribution", "direction")]
# Stable ordering: curated facets in punchline order, then by contrib desc.
module_levels <- c("heatshock_stress", "shared_angio_glucose",
                   "autoreg_feedback", "hif1a_hypoxic_core", "other_unclassified")
attr_df$module <- factor(attr_df$module, levels = module_levels)
attr_df <- attr_df[order(attr_df$module, -attr_df$contrib), ]
attr_df$module <- as.character(attr_df$module)

out_data <- file.path(TBL_DIR, "fig3l_hif_attribution_data.csv")
write.csv(attr_df, out_data, row.names = FALSE)
cat(sprintf("  wrote %s  (%d regulon members)\n", basename(out_data), nrow(attr_df)))

# -----------------------------------------------------------------------------
# MODULE SUMMARY -- aggregations only (no statistics). Lets the README quote
# "heat-shock bucket > net HIF core" without any recomputation.
# -----------------------------------------------------------------------------
total_abs <- sum(abs(attr_df$contrib))
summ <- do.call(rbind, lapply(module_levels, function(m) {
  sub <- attr_df[attr_df$module == m, , drop = FALSE]
  data.frame(
    module      = m,
    n_targets   = nrow(sub),
    sum_contrib = sum(sub$contrib),
    # SIGNED pct (matches the signed convention in fig3g_target_decomposition_summary):
    # numerator is the signed module sum so the sign tracks sum_contrib
    # (hif1a_hypoxic_core sum_contrib = -8.46 -> pct_of_total = -0.56%).
    pct_of_total = 100 * sum(sub$contrib) / total_abs,
    n_up        = sum(sub$direction == "Up"),
    n_down      = sum(sub$direction == "Down"),
    stringsAsFactors = FALSE
  )
}))

out_summ <- file.path(TBL_DIR, "fig3l_module_summary.csv")
write.csv(summ, out_summ, row.names = FALSE)
cat(sprintf("  wrote %s\n", basename(out_summ)))

# -----------------------------------------------------------------------------
# Console verification (the numbers the handoff must report).
# -----------------------------------------------------------------------------
net_core <- summ$sum_contrib[summ$module == "hif1a_hypoxic_core"]
stress   <- summ$sum_contrib[summ$module == "heatshock_stress"]
cat("\n[VERIFY] per-module sum_contrib / n_up / n_down:\n")
print(summ[, c("module", "n_targets", "sum_contrib", "n_up", "n_down")], row.names = FALSE)
cat(sprintf("\n[VERIFY] heatshock_stress sum_contrib = %+.2f  (should be strongly POSITIVE)\n", stress))
cat(sprintf("[VERIFY] hif1a_hypoxic_core sum_contrib = %+.2f  (should be NEGATIVE = repressed)\n", net_core))
cat(sprintf("[VERIFY] stress (%+.2f) %s net HIF core (%+.2f)\n",
            stress, ifelse(stress > net_core, "EXCEEDS", "does NOT exceed"), net_core))

cat("[DONE] 03c_hif_program_attribution.R -- fig3l tidy tables emitted (join only, no statistics).\n")
