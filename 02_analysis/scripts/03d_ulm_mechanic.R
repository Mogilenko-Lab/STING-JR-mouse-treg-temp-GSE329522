#!/usr/bin/env Rscript
# =============================================================================
# 03d_ulm_mechanic.R - PHASE 3 (stage 04_tf) COMPUTE: fig3m
#   "How decoupleR-ULM nominates a TF" -- the SCORING MECHANIC, pedagogical.
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# Role:    COMPUTE half of the "normalize-then-visualize" split for the
#          attribution-arc MECHANIC figure (fig3m), matched partner of
#          fig3g-expanded (landscape) and fig3l (biology). Contains NO ggplot()/
#          ggsave() and NO new statistics whatsoever. It only
#            (i)  COPIES the already-decomposed Hif1a regulon
#                 (fig3g_target_decomposition_data.csv: source/target/mor/t_wt/
#                 contrib/class) -- the promiscuous worked example; and
#            (ii) DESCRIPTIVELY JOINS one SPECIFIC comparator TF's CollecTRI
#                 regulon to the EXISTING WT_heat limma t-vector (a lookup, NOT
#                 a model fit -- no run_ulm). aligned_contrib = sign(mor)*t_wt is
#                 an arithmetic relabel of existing t-statistics.
#
# COMPARATOR-REGULON OPTION (design-spec sec 2 fig3m note): OPTION (b) is used.
#   Rationale: the figure's shared x-axis is the PER-TARGET aligned contribution,
#   and Panel B must show a beeswarm of the comparator regulon's individual
#   members on that same axis. Option (a) (reuse fig3j/master aggregate values)
#   only yields ONE aggregate score per TF -- it has NO per-target distribution,
#   so it cannot populate the beeswarm and is too thin to make the
#   specific-vs-promiscuous contrast legible. Option (b) -- a DESCRIPTIVE join of
#   the existing CollecTRI Stat2 regulon to the existing WT_heat t-vector in
#   02_de_results.rds -- is therefore used. This is a lookup only; NO run_ulm,
#   NO run_mlm, NO p.adjust, NO model refit. FLAG FOR ANTON: option (a) was too
#   thin (aggregate-only) to render the mechanic, so the spec's descriptive
#   fallback (b) was taken.
#
# Comparator TF = Stat2 (SPECIFIC): a small, canonical type-I-IFN-restricted
#   regulon (n~32 members present in WT_heat). On the WT_heat (fever-heat) main
#   contrast its members do NOT coherently move, so the regulon does NOT pile up
#   -- the visual foil to Hif1a's broad regulon that accumulates generic high-|t|
#   stress genes regardless of HIF biology. This figure is about the SCORING
#   MECHANIC ONLY -- it makes NO biology claim about Stat2, IFN, HIF1a, or HIF2a.
#
# Inputs (NO recomputation):
#   - 03_results/04_tf/tables/fig3g_target_decomposition_data.csv  (Hif1a; copy)
#   - 03_results/objects/net_collectri_mouse.rds                   (Stat2 regulon)
#   - 03_results/objects/02_de_results.rds  (WT_heat limma t-vector; lookup only)
#
# Outputs (03_results/04_tf/tables/):
#   - fig3m_ulm_mechanic_data.csv      (per-member: tf/regulon_class/target/mor/
#                                        t_wt/aligned_contrib/is_stress_contaminant)
#   - fig3m_ulm_mechanic_summary.csv   (per-tf aggregates; no statistics)
#
# Object name (HARD GUARDRAIL, design-spec sec 0): the thing the Hif1a panel
#   describes is "a heat-induced glycolytic/stress program partially overlapping
#   HIF targets" -- NEVER "the HIF program" / "HIF1a the TF". Titles/subtitles
#   keep object-name discipline; the comparator NEVER implies HIF2a anything.
#
# Method: COPY + descriptive JOIN + simple aggregation only. No lmFit/eBayes/
#   run_ulm/run_mlm/p.adjust/prcomp anywhere. Figures live in the _viz.R partner.
# =============================================================================

source("02_analysis/config/config.R")
load_packages(extra = c("readr"))   # only base + dplyr used; no statistics

TBL_DIR <- stage_dir("04_tf", "tables")

# -----------------------------------------------------------------------------
# SINGLE SOURCE OF TRUTH (script header constants)
# -----------------------------------------------------------------------------
TF_PROMISCUOUS <- "Hif1a"   # worked example: broad regulon -> score = pile-up
TF_SPECIFIC    <- "Stat2"   # comparator: small IFN-restricted regulon (foil)

# Stress-contaminant genes that "drag the Hif1a score up" (design-spec fig3m).
STRESS_CONTAMINANTS <- c("Hspa1a", "Timp1", "Sdc1", "Cdkn1a",
                         "Serpine1", "Eno2", "Spp1")

# -----------------------------------------------------------------------------
# PANEL A -- Hif1a (PROMISCUOUS). COPY the existing fig3g decomposition verbatim.
#   aligned_contrib := the existing `contrib` column (= sign(mor)*t_wt already);
#   recomputing it from scratch is forbidden, so we COPY it.
# -----------------------------------------------------------------------------
fig3g <- read.csv(file.path(TBL_DIR, "fig3g_target_decomposition_data.csv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(all(c("source", "target", "mor", "t_wt", "contrib", "class") %in% names(fig3g)))
stopifnot(all(fig3g$source == TF_PROMISCUOUS))

hif <- data.frame(
  tf                   = TF_PROMISCUOUS,
  regulon_class        = "promiscuous",
  target               = fig3g$target,
  mor                  = as.integer(fig3g$mor),
  t_wt                 = fig3g$t_wt,
  aligned_contrib      = fig3g$contrib,                 # COPIED (= sign(mor)*t_wt)
  is_stress_contaminant = fig3g$target %in% STRESS_CONTAMINANTS,
  stringsAsFactors     = FALSE
)

# Belt-and-braces: confirm the copied contrib equals sign(mor)*t_wt (no drift).
chk <- max(abs(hif$aligned_contrib - sign(hif$mor) * hif$t_wt))
stopifnot(chk < 1e-8)

# -----------------------------------------------------------------------------
# PANEL B -- Stat2 (SPECIFIC). DESCRIPTIVE JOIN: CollecTRI Stat2 regulon targets
#   looked up against the EXISTING WT_heat limma t-vector. NO run_ulm.
#   aligned_contrib := sign(mor)*t_wt (arithmetic relabel of existing t-stats).
# -----------------------------------------------------------------------------
net <- readRDS(file.path(DIR_OBJECTS, "net_collectri_mouse.rds"))
de  <- readRDS(file.path(DIR_OBJECTS, "02_de_results.rds"))
wt  <- de[["WT_heat"]]
stopifnot(all(c("gene_symbol", "t") %in% names(wt)))

tvec <- setNames(wt$t, wt$gene_symbol)             # existing WT_heat t-statistics
reg  <- net[net$source == TF_SPECIFIC, c("source", "target", "mor")]
reg$t_wt <- tvec[reg$target]                       # lookup only
reg  <- reg[!is.na(reg$t_wt), , drop = FALSE]      # members measured in WT_heat

spec <- data.frame(
  tf                   = TF_SPECIFIC,
  regulon_class        = "specific",
  target               = reg$target,
  mor                  = as.integer(reg$mor),
  t_wt                 = as.numeric(reg$t_wt),
  aligned_contrib      = sign(reg$mor) * as.numeric(reg$t_wt),  # relabel, not a stat
  is_stress_contaminant = FALSE,   # contaminant flag scoped to the Hif1a regulon
  stringsAsFactors     = FALSE
)

# -----------------------------------------------------------------------------
# COMBINE + emit per-member data (design-spec schema).
# -----------------------------------------------------------------------------
mech <- rbind(hif, spec)
mech <- mech[, c("tf", "regulon_class", "target", "mor",
                 "t_wt", "aligned_contrib", "is_stress_contaminant")]

out_data <- file.path(TBL_DIR, "fig3m_ulm_mechanic_data.csv")
write.csv(mech, out_data, row.names = FALSE)
cat(sprintf("  wrote %s  (%d members: %s=%d promiscuous, %s=%d specific)\n",
            basename(out_data), nrow(mech),
            TF_PROMISCUOUS, nrow(hif), TF_SPECIFIC, nrow(spec)))

# -----------------------------------------------------------------------------
# SUMMARY -- per-tf aggregations only (no statistics). The mean aligned
# contribution is the regulon's weighted center -- where the viz draws the
# "score is the pile-up" diamond. pct_contrib_from_stress is computed for the
# Hif1a regulon only (NA for the comparator).
# -----------------------------------------------------------------------------
summ <- do.call(rbind, lapply(split(mech, mech$tf), function(d) {
  total_abs <- sum(abs(d$aligned_contrib))
  pct_stress <- if (any(d$is_stress_contaminant)) {
    100 * sum(abs(d$aligned_contrib[d$is_stress_contaminant])) / total_abs
  } else NA_real_
  mean_ac <- mean(d$aligned_contrib)
  data.frame(
    tf                     = d$tf[1],
    regulon_class          = d$regulon_class[1],
    n_targets              = nrow(d),
    mean_aligned_contrib   = mean_ac,
    pct_contrib_from_stress = pct_stress,
    ulm_score_sign         = ifelse(mean_ac >= 0, "positive", "negative"),
    stringsAsFactors       = FALSE
  )
}))
# Order: promiscuous worked example first, then specific comparator.
summ <- summ[order(summ$regulon_class != "promiscuous"), ]

out_summ <- file.path(TBL_DIR, "fig3m_ulm_mechanic_summary.csv")
write.csv(summ, out_summ, row.names = FALSE)
cat(sprintf("  wrote %s\n", basename(out_summ)))

# -----------------------------------------------------------------------------
# Console verification (numbers the handoff reports).
# -----------------------------------------------------------------------------
cat("\n[VERIFY] per-tf summary:\n")
print(summ, row.names = FALSE)
hif_pct <- summ$pct_contrib_from_stress[summ$tf == TF_PROMISCUOUS]
cat(sprintf("\n[VERIFY] %s regulon: %d members, mean aligned contrib %+.3f (%s ULM sign)\n",
            TF_PROMISCUOUS,
            summ$n_targets[summ$tf == TF_PROMISCUOUS],
            summ$mean_aligned_contrib[summ$tf == TF_PROMISCUOUS],
            summ$ulm_score_sign[summ$tf == TF_PROMISCUOUS]))
cat(sprintf("[VERIFY] %.1f%% of |aligned contribution| comes from 7 generic stress contaminants\n", hif_pct))
cat(sprintf("[VERIFY] %s regulon: %d members, mean aligned contrib %+.3f (%s ULM sign)\n",
            TF_SPECIFIC,
            summ$n_targets[summ$tf == TF_SPECIFIC],
            summ$mean_aligned_contrib[summ$tf == TF_SPECIFIC],
            summ$ulm_score_sign[summ$tf == TF_SPECIFIC]))

cat("[DONE] 03d_ulm_mechanic.R -- fig3m tidy tables emitted (copy + descriptive join; no statistics).\n")
