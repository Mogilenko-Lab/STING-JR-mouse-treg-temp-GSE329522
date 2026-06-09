#!/usr/bin/env Rscript
# =============================================================================
# 03f_experiment_tiers.R - PHASE 3 (stage 04_tf) COMPUTE: fig3o
#   The tiered "what to run next" experiment table (CONFIG / LOOKUP ONLY).
# Project: GSE329522 STING/cGAS Hyperthermia iTreg (2x2 genotype x temperature)
#
# Role:    COMPUTE half of the "normalize-then-visualize" split for fig3o, BUT
#          this is a hand-authored experimental-design CONFIG TABLE -- it is the
#          single source of truth for the tier ladder. It contains NO statistics
#          whatsoever (no lmFit/eBayes/run_ulm/run_mlm/p.adjust/prcomp, no DE, no
#          data read). It emits ONE tidy CSV; the figure lives in the _viz.R
#          partner, which renders the conceptual ladder from this table.
#
# WHY THIS FIGURE: the deck corrects the attribution ("HIF1a #1" is a regulon-
#   aggregation artifact; the rank is method-fragile) and then must hand the wet
#   lab a PLAN, not just caveats. The tiered design says: exclude the cheap
#   confound first, then ask "is it HIF at all?", and only then run the expensive
#   isoform genetics. The guardrail it persists: DO NOT lead with floxed-Hif1a vs
#   floxed-Epas1 genetics.
#
# Output (03_results/04_tf/tables/):
#   - fig3o_experiment_tiers_data.csv
#       cols: tier, tier_label, experiment, what_it_excludes_or_tests, readout,
#             cost_tier, gate, row_in_tier
#
# Framing constraints (bind this table): necessity worded as "some activity the
#   drug perturbs is required at 39C", NEVER "HIF1a is required". belzutifan is an
#   alternative the drug cannot exclude (HIF2a-SELECTIVE tool), never a crowning of
#   HIF2a. Tier 2 genetics is GATED, not led-with.
# =============================================================================

source("02_analysis/config/config.R")
load_packages(extra = c("readr"))   # dplyr pulled in; NO statistics used

TBL_DIR <- stage_dir("04_tf", "tables")

cat("=================================================================\n")
cat("PHASE 3 (fig3o): tiered experiment design table (config only, no stats)\n")
cat("=================================================================\n\n")

# -----------------------------------------------------------------------------
# THE TIER LADDER -- single source of truth. Hand-authored from the
# constructive/adversarial synthesis. One row per experiment arm.
#   tier 0 = exclude the heat x drug-cytotoxicity confound (cheap)
#   tier 1 = "is it HIF at all?"                          (medium)
#   tier 2 = the isoform mind-changer                     (expensive; gated)
# -----------------------------------------------------------------------------
tiers <- tibble::tribble(
  ~tier, ~tier_label,                                  ~experiment,                                  ~what_it_excludes_or_tests,                                                            ~readout,                                              ~cost_tier,

  # ---- TIER 0: exclude the heat x drug-cytotoxicity confound -----------------
  0L, "Tier 0 - exclude the confound",                "Inactive structural analog / vehicle @ 39C", "Excludes: the effect is heat x drug-cytotoxicity, not a specific perturbed activity", "iTreg viability / suppressive function at 39C",       "cheap",
  0L, "Tier 0 - exclude the confound",                "Dose-response of the drug @ 39C",            "Tests: is the 39C effect on-target and dose-dependent, or a cytotoxic cliff?",        "dose vs iTreg readout curve at 39C",                 "cheap",
  0L, "Tier 0 - exclude the confound",                "Tconv comparator @ 39C (HIF-agnostic)",      "Excludes: the drug also kills/impairs Tconv at 39C (generic heat x drug toxicity)",  "Tconv viability/proliferation at 39C vs iTreg",      "cheap",

  # ---- TIER 1: is it HIF at all? ---------------------------------------------
  1L, "Tier 1 - is it HIF at all?",                    "Chaperone / proteostasis inhibitor @ 39C",   "Tests: is the required activity HSP90/HSP70 chaperone dependence (the heat-shock module that dominates fig3l)?", "iTreg readout with HSP90/HSP70 block (e.g. 17-AAG / VER-155008) at 39C", "medium",
  1L, "Tier 1 - is it HIF at all?",                    "Glycolysis inhibitor @ 39C",                 "Tests: whether the shared glycolytic members (one part of the program, NOT the dominant module) carry the dependency", "iTreg readout with 2-DG / glycolysis block at 39C",  "medium",
  1L, "Tier 1 - is it HIF at all?",                    "HSF1 readout / knockdown @ 39C",             "Tests: is the heat-shock axis (HSF1; fig3n co-elevated) the driver rather than HIF?", "iTreg readout vs HSF1 activity / knockdown at 39C",  "medium",
  1L, "Tier 1 - is it HIF at all?",                    "Belzutifan (HIF2a-selective) cross-check",   "Tests: an alternative the drug cannot exclude -- is a HIF2a-selective block sufficient?", "iTreg readout under HIF2a-selective inhibition",     "medium",

  # ---- TIER 2: the isoform mind-changer (gated) ------------------------------
  2L, "Tier 2 - isoform genetics (gated)",             "Floxed-Hif1a iTreg @ 39C",                   "Tests: does deleting HIF1a (specifically) abolish the 39C effect?",                    "iTreg readout in Hif1a-cKO vs control at 39C",       "expensive",
  2L, "Tier 2 - isoform genetics (gated)",             "Floxed-Epas1 (HIF2a) iTreg @ 39C",           "Tests: does deleting HIF2a instead abolish it -- which isoform, if either?",           "iTreg readout in Epas1-cKO vs control at 39C",       "expensive"
)

# -----------------------------------------------------------------------------
# GATE column -- what must hold UPSTREAM before each arm is licensed. This
# encodes the discipline "don't lead with the cross": Tier 2 is gated on Tiers 0-1.
# -----------------------------------------------------------------------------
gate_by_tier <- c(
  "0" = "Start here: licensed immediately (the cheap exclusion the deck demands first).",
  "1" = "Run only AFTER Tier 0 excludes the heat x drug-cytotoxicity confound.",
  "2" = "Run ONLY after Tiers 0-1 -- do NOT lead with floxed-Hif1a vs floxed-Epas1 genetics."
)
tiers$gate <- gate_by_tier[as.character(tiers$tier)]

# row index within tier (for ladder layout in the viz; not a statistic).
tiers <- tiers %>%
  group_by(tier) %>%
  mutate(row_in_tier = row_number()) %>%
  ungroup() %>%
  arrange(tier, row_in_tier)

# -----------------------------------------------------------------------------
# Light schema checks (NOT statistics): three tiers present, costs ascend.
# -----------------------------------------------------------------------------
stopifnot(all(c(0L, 1L, 2L) %in% tiers$tier))
cost_rank <- c("cheap" = 1L, "medium" = 2L, "expensive" = 3L)
stopifnot(all(tiers$cost_tier %in% names(cost_rank)))
# cost is non-decreasing with tier (the ladder ascends in cost).
agg_cost <- tapply(cost_rank[tiers$cost_tier], tiers$tier, max)
stopifnot(all(diff(agg_cost) >= 0))

out <- file.path(TBL_DIR, "fig3o_experiment_tiers_data.csv")
readr::write_csv(tiers, out)
cat(sprintf("  [SAVE] %-42s %4d rows x %d cols\n",
            basename(out), nrow(tiers), ncol(tiers)))

cat("\n[VERIFY] tier ladder (tier / cost / experiment):\n")
print(as.data.frame(tiers[, c("tier", "cost_tier", "experiment")]), row.names = FALSE)

cat("\n[DONE] fig3o COMPUTE complete (config table only). Run 03f_experiment_tiers_viz.R.\n")
