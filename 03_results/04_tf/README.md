# 04_tf: TF activity (decoupleR) and the Hif1a rank cascade

Transcription-factor activity inferred by decoupleR (CollecTRI ULM, primary; DoRothEA-ABC ULM, comparator; CollecTRI MLM and consensus, forensics) from limma t-statistics across seven contrasts of the GSE329522 2×2 cGAS×temperature iTreg factorial design (murine iTreg, WT vs cGAS-KO × 37 °C vs 39 °C fever-heat, n=5/group). The central question this phase answers is whether the HIF/glycolytic program shows cGAS-dependent activation: the data do not support that conclusion at n=5. The IFN members Irf3/Stat2/Stat1 reach significance on the Interaction contrast (padj < 0.05), while the HIF axis (Hif1a/Epas1) is flat and non-significant on the same contrast — asymmetry, not independence.

A secondary forensic thread explains why Hif1a's rank varies from #1 to #142 across inference methods and networks, and concludes that its apparent activation is substantially a heat-generic stress/glycolytic phenomenon whose rank instability is fully accounted for by regulon size, target collinearity, and consensus aggregation arithmetic. A three-figure attribution arc (fig3m → fig3g → fig3l) forms the analytical spine: mechanic (HOW ULM nominates a TF) → landscape (WHERE Hif1a's score comes from: the full ranked 353-member contribution distribution) → biology (WHAT KIND of genes and in which direction: heat-shock/stress up, HIF-selective hypoxic core down). This arc establishes that the program revealed here is **a heat-induced glycolytic/stress program partially overlapping HIF targets**, not a canonical hypoxic-HIF output.

**Sample mapping**: confirmed against the owner's sample sheet on 2026-07-22, 20 of 20 libraries concordant with the label-blind marker call.

---

## Generating Scripts

| Script | Role | Outputs |
|---|---|---|
| `02_analysis/scripts/03_decoupler_tf.R` | COMPUTE — all statistics, no ggplot | `03_results/objects/03_tf_collectri.rds`, `03_results/objects/03_tf_forensics.rds`, `03_results/master/master_tf_activities.csv`, all `tables/fig3a–fig3k_*.csv` |
| `02_analysis/scripts/03_decoupler_tf_viz.R` | VIZ — reads tidy tables, no recomputation | `figures/fig3a–fig3i_*.pdf`; redesigned `figures/fig3g_target_decomposition.pdf` (expanded full-landscape view) |
| `02_analysis/scripts/03b_decoupler_method_comparison.R` | COMPUTE — six-statistic method comparison, no ggplot | `tables/fig3j_allmethods_topTF_data.csv`, `tables/fig3k_method_rank_divergence_data.csv`, `tables/fig3k_method_rank_spearman.csv` |
| `02_analysis/scripts/03b_decoupler_method_comparison_viz.R` | VIZ — reads method-comparison tidy tables, no recomputation | `figures/fig3j_topTF_allmethods_WT_heat.pdf`, `figures/fig3j_topTF_allmethods_Interaction.pdf`, `figures/fig3k_method_rank_divergence.pdf` |
| `02_analysis/scripts/03c_hif_program_attribution.R` | COMPUTE — join of fig3g table to curated module/isoform lookup; no new statistics | `tables/fig3l_hif_attribution_data.csv`, `tables/fig3l_module_summary.csv` |
| `02_analysis/scripts/03c_hif_program_attribution_viz.R` | VIZ — reads fig3l tables, no recomputation | `figures/fig3l_hif_attribution.pdf` |
| `02_analysis/scripts/03d_ulm_mechanic.R` | COMPUTE — join/copy of fig3g contributions + comparator regulon; no run_ulm | `tables/fig3m_ulm_mechanic_data.csv`, `tables/fig3m_ulm_mechanic_summary.csv` |
| `02_analysis/scripts/03d_ulm_mechanic_viz.R` | VIZ — reads fig3m tables, no recomputation | `figures/fig3m_ulm_mechanic.pdf` |
| `02_analysis/scripts/03e_heat_main_regulators.R` | COMPUTE — `run_ulm` on a fresh heat-MAIN t-stat matrix (same call as `03_decoupler_tf.R`); no ggplot | `tables/fig3n_heat_main_regulators_data.csv` |
| `02_analysis/scripts/03e_heat_main_regulators_viz.R` | VIZ — reads fig3n table, no recomputation | `figures/fig3n_heat_main_regulators.pdf` |
| `02_analysis/scripts/03g_nonidentifiability.R` | COMPUTE — `run_ulm` on heat-MAIN t-stat vector + CollecTRI set-membership counting; no ggplot | `tables/fig3p_heatmain_ranking_data.csv`, `tables/fig3q_coregulators_data.csv`, `tables/fig3r_membership_data.csv` |
| `02_analysis/scripts/03g_nonidentifiability_viz.R` | VIZ — reads the three fig3p/3q/3r tables, no recomputation; shared `FAMILY_COLORS` palette | `figures/fig3p_heatmain_ranking.{pdf,png}`, `figures/fig3q_coregulators.{pdf,png}`, `figures/fig3r_shared_ownership.{pdf,png}` |

---

## Narrative deck sequence

```
AFFIRM ──────────────▶ fig3i  (real phenomenon; shared HIF/glycolytic target Vegfa up in both genotypes; IFN Ifit1 cGAS-dependent)

ATTRIBUTION ARC ─────▶ fig3m  (HOW ULM nominates a TF — regulon-weighted aggregation; promiscuous vs specific)
                       fig3g* (WHERE Hif1a's score comes from — full ranked 353-member landscape; stress genes dominate)
                       fig3l  (WHAT KIND — module-bucketed sign-split; heat-shock/stress UP, HIF-selective hypoxic core DOWN)

LIBRARY STACK ───────▶ fig3d  (WHY #1: tight DoRothEA regulon → #1; fuller CollecTRI regulon dilutes → #12)
                       fig3e  (WHY not robust: 92% targets shared, mean 22 TFs/target; MLM re-attributes, #12 → #142)

NON-IDENTIFIABLE ────▶ fig3c  (rank cascade recap #1 → #12 → #142 → #8)
                       fig3k  (MLM lone de-confounder; reshuffles Rela, Nfkb1, Hif1a alike)

TWO ARMS ────────────▶ fig3a  (IFN interaction-significant / HIF flat — asymmetry, not independence, at n=5)

HSF1 GAP ────────────▶ fig3n  (heat-shock Hsf1 co-elevated on heat-MAIN; same CollecTRI-ULM estimator as fig3a/fig3c)
```

`*` fig3g = expanded viz (full ranked landscape), same underlying table `fig3g_target_decomposition_data.csv`.

### Non-identifiability cluster (fig3p / fig3q / fig3r) — deck placement PENDING

A three-figure triptych making one combined claim: on the heat-main contrast the activation signal **cannot be attributed to Hif1a — or to any single TF — as the driver**. fig3p shows Hif1a is only #9 in a dense crowd of generic stress/immediate-early/NF-κB TFs (no clean winner); fig3q shows 92% of Hif1a's CollecTRI targets are shared by the network's most promiscuous, non-hypoxia-specific regulators; fig3r overlays the two — the same heat-driven genes populate many TFs' regulons. Deck placement is **PENDING** (the conductor will slot the cluster — likely after fig3e, whose target-sharing it generalises — but no final deck-spine position is asserted here). All three share a single `FAMILY_COLORS` palette defined in `03g_nonidentifiability_viz.R` (heat-shock = purple `#762A83`; proliferation/metabolism = magenta `#C51B7D`; HIF axis = orange; NF-κB = blue; AP-1/immediate-early = red; etc.).

```
NON-IDENTIFIABILITY ─▶ fig3p  (heat-MAIN ranking: no clean winner — Hif1a #9 in a generic-stress crowd; canonical Hsf1 far down at #50)
   (placement PENDING)   fig3q  (Hif1a's 353 targets are 92% shared by the most promiscuous non-hypoxia regulators — Sp1/Trp53/NF-κB/AP-1/Myc)
                         fig3r  (shared ownership: top-20 Hif1a targets × its 12 largest sharers — the contrast cannot single out any one TF)
```

### Final wet-lab deck spine (serendipity-resolution ordering)

`fig3i → fig3m → fig3g → fig3l → fig3d → fig3e → fig3c → fig3k → fig3a → fig3n`

(affirm what's real → HOW ULM nominates → WHERE the score comes from → WHAT KIND of genes + direction → **WHY #1 happened (regulon dilution, #1→#12)** → **WHY it isn't robust (target collinearity, #12→#142)** → rank cascade recap → MLM reshuffle → cGAS asymmetry → HSF1 co-elevation.)

**fig3d and fig3e sit in the deck spine.** They carry the account of how the collaborator's DoRothEA-ULM produced HIF1α #1 (#1→#12 regulon dilution) and why CollecTRI-MLM does not hold it (#12→#142 collinearity redistribution). They sit between fig3l (what the genes are) and fig3c (the full cascade) so the mechanism precedes the cascade recap, which lets fig3c read as a recap.


**Archive-only (retained in phase directory, NOT part of the deck):** fig3b, fig3f, the two fig3j variants, and **lombardi_recurrence** — canonical PDFs and tables are retained; they are simply not placed in the deck. See DROP notes under each legend below.

---

## Figure Legends

---

### fig3a_tf_interaction_axes.pdf

**The IFN axis members Irf3 (+5.93, padj = 1.04×10⁻⁶), Stat2 (+5.85, padj = 1.10×10⁻⁶), and Stat1 (+5.30, padj = 1.88×10⁻⁵) are significantly interaction-positive — their activation is cGAS-dependent — at BH padj < 0.05 on the cGAS×heat Interaction contrast; NFkB members Nfkb1 (+2.16, padj = 0.598) and Rela (+2.03, padj = 0.631) lean positive but are non-significant at n=5; and the HIF axis (Hif1a −0.19, padj = 0.97; Epas1 −0.86, padj = 0.91) is flat with no detectable cGAS-dependence. The two arms are asymmetric: the IFN arm is cGAS-dependent, and the HIF arm carries no detectable cGAS-dependence at n=5.**

**Method.** CollecTRI ULM (`run_ulm(.mor="mor", minsize=5)`) applied to limma t-statistics from the 2×2 Interaction contrast (captures cGAS-dependent temperature response); BH correction within contrast; TFs classified as HIF axis {Hif1a, Epas1}, IFN axis {Irf1, Irf3, Irf7, Stat1, Stat2, Nfkb1, Rela}, or other; lollipop plot on ULM score axis, significant TFs (padj < 0.05) starred.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf_viz.R` |
| Function | `run_ulm(.source="source", .target="target", .mor="mor", minsize=5)` + BH within Interaction contrast |
| Input | `03_results/04_tf/tables/fig3a_tf_interaction_axes_data.csv` |

---

### fig3b_top_tf_by_contrast.pdf

> **ARCHIVE ONLY — NOT in narrative deck.** The "heat drives broad stress/glycolytic program" point is made more sharply and on-message by the expanded fig3g + fig3l, which name the stress genes explicitly; the "Interaction isolates IFN" point belongs to fig3a. The 96-row 4-facet top/bottom-12 grid is dense and off the single-message spine. Retained as a forensic backup.

**The top-12 activated and top-12 suppressed TFs per contrast (CollecTRI ULM, BH) reveal that heat drives a broad stress-response and glycolytic program in both genotypes (WT_heat and KO_heat) while the Interaction contrast selectively isolates the IFN/Irf9/Stat2 axis; the HIF axis shows no detectable signal on the Interaction contrast.**

**Method.** For each of four contrasts (WT_heat, KO_heat, Temp_main, Interaction), the top and bottom 12 TFs by ULM score from the CollecTRI run are extracted and faceted on a shared score axis; HIF-axis and IFN-axis members annotated; all facets use the same x-axis to enable direct cross-contrast comparison.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf_viz.R` |
| Function | `run_ulm(.mor="mor", minsize=5)` — scores read from pre-computed table; top/bottom 12 per contrast (`TOPN_CONTRAST = 12`) |
| Input | `03_results/04_tf/tables/fig3b_top_tf_by_contrast_data.csv` |

---

### fig3c_hif1a_rank_cascade.pdf

**Hif1a's activation rank on the WT_heat contrast is highly method/network-dependent — #1 of 263 (DoRothEA-ULM, score = 6.09), #12 of 658 (CollecTRI-ULM, score = 5.11), #142 of 658 (CollecTRI-MLM, score = 1.13), #8 of 658 (CollecTRI-consensus, score = 2.89) — so the apparent rank tracks the inference choice, where a stable biological signal would hold across all four; its flat, non-significant Interaction score (ULM = −0.19, padj = 0.97) confirms no detectable cGAS-dependence regardless of which method is used.**

**Method.** Four inferences run on the WT_heat t-statistic vector: DoRothEA-ABC ULM (`run_ulm`, net_dorothea_mouse_ABC), CollecTRI ULM (`run_ulm`, net_collectri_mouse), CollecTRI MLM (`run_mlm`), and CollecTRI consensus (`decouple(statistics=c("ulm","mlm","wsum"), consensus_score=TRUE)`); Hif1a rank extracted from each ranked list; percentage rank computed as rank/n_tfs×100.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf_viz.R` |
| Function | `run_ulm` / `run_mlm` / `decouple(..., consensus_score=TRUE)` on WT_heat t-stat vector |
| Input | `03_results/04_tf/tables/fig3c_hif1a_rank_cascade_data.csv` |

---

### fig3d_regulon_swap.pdf  *(REINSTATED — library-stack spine, "why #1 happened")*

**The drop from rank #1 (DoRothEA-ULM) to #12 (CollecTRI-ULM) is explained by regulon expansion: CollecTRI's Hif1a regulon (353 targets, mean|t| = 4.26) is 2.7× larger than DoRothEA's (131 targets, mean|t| = 4.70), with 252 CollecTRI-only targets that carry high |t| but mixed sign (median signed t ≈ 0.3), adding magnitude with little net directional signal, which dilutes the mor-weighted ULM score.**

SERENDIPITY ROLE: This is the regulon-definition effect that mechanically explains the collaborator's "#1" call. DoRothEA's smaller, tighter Hif1α regulon concentrates the high-t heat targets, producing the #1 rank under ULM; CollecTRI's fuller regulon dilutes that concentration to #12. The "#1" is therefore a property of which regulon library was used, not a biological fact — explained without invoking any biology.

**Method.** DoRothEA and CollecTRI Hif1a target lists are intersected against the expression universe; per-gene t-statistics (WT_heat) retrieved; targets classified as `both` (n = 101), `dorothea_only` (n = 30), or `collectri_only` (n = 252); per-network mean|t|, median t, and %|t|>2 computed.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf_viz.R` |
| Function | Set arithmetic on filtered `net_dorothea_mouse_ABC` and `net_collectri_mouse`; `mean(abs(t_wt))` per network |
| Input | `03_results/04_tf/tables/fig3d_regulon_swap_data.csv`, `03_results/04_tf/tables/fig3d_regulon_swap_summary.csv` |

---

### fig3e_mlm_collinearity.pdf  *(REINSTATED — library-stack spine, "why it isn't robust")*

**The drop from rank #12 (CollecTRI-ULM, score = 5.11) to #142 (CollecTRI-MLM, score = 1.13) is explained by target collinearity: 92% (frac_co_regulated = 0.9207) of Hif1a's 353 CollecTRI targets are co-regulated by at least one other TF in the same network (mean = 22.23 co-regulators per target, median = 11), so the multivariate MLM re-assigns shared signal away from Hif1a.**

SERENDIPITY ROLE: This is the collinearity-redistribution half of the "#1" explanation. Because 92% of Hif1α's CollecTRI targets are co-regulated by a mean of 22 other TFs, the multivariate MLM distributes that shared signal across the competing regulons rather than crediting it all to Hif1α, and the score collapses 5.11/#12 (ULM) → 1.13/#142 (MLM). Together with fig3d (regulon-definition effect), this accounts for the collaborator's "#1" call mechanically — without invoking any biology. The broad, defensible statement is that Hif1α's signal is shared across many stress-responsive TFs, not uniquely its own.

**Method.** For each Hif1a CollecTRI target present in the expression universe, all other CollecTRI TFs also targeting that gene are counted; co-regulated flag set when n_other_TFs ≥ 1; fraction and mean computed across the 353-target regulon; two-bar score-collapse panel shows ULM score 5.11 → MLM score 1.13 for the same network and data.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf_viz.R` |
| Function | `count(target, name="n_tfs_total")` across CollecTRI edges; `run_mlm(.mor="mor", minsize=5)` |
| Input | `03_results/04_tf/tables/fig3e_mlm_collinearity_data.csv`, `03_results/04_tf/tables/fig3e_mlm_collinearity_summary.csv`, `03_results/04_tf/tables/fig3e_score_collapse.csv` |

---

### fig3f_consensus_zmix.pdf

> **ARCHIVE ONLY — NOT in narrative deck.** The #142→#8 arithmetic (4:1 univariate z-vote drowning MLM) is a methods-internal detail; including it in the deck invites the "consensus recovers HIF1a" mis-read. Archive only; cite in prose if a methods reviewer presses.

**The partial recovery from rank #142 (MLM) to #8 (consensus, score = 2.89) occurs because the decoupleR consensus is an unweighted mean of folded-z scores across all statistics: the four univariate-family statistics (ulm z = 2.44, wsum z = 4.13, norm_wsum z = 2.72, corr_wsum z = 4.31) numerically outvote the single multivariate statistic (mlm z = 0.85), re-inflating Hif1a's rank while inheriting univariate double-counting.**

**Method.** `decouple(statistics=c("ulm","mlm","wsum"), consensus_score=TRUE)` run on the WT_heat vector; non-consensus rows folded (mirror + scale by symmetric sd) to produce z-scores per statistic; consensus = mean of folded-z across all emitted statistics; Hif1a's per-statistic z-scores tabulated; algebraic mean verified to reproduce the stored consensus score (diff < 1×10⁻⁶).

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf_viz.R` |
| Function | `decouple(..., consensus_score=TRUE)`; folded-z via sign-symmetric sd standardisation per statistic then `mean(z_score)` |
| Input | `03_results/04_tf/tables/fig3f_consensus_zmix_data.csv` |

---

### fig3g_target_decomposition.pdf  *(REDESIGNED — expanded full-landscape view)*

**Hif1a's positive CollecTRI-ULM score on the WT_heat contrast is built almost entirely from its 340 "other" regulon members — generic stress/ECM/proliferation genes that happen to sit in the regulon. The six highest-contributing stress/remodelling genes (Timp1 +24.92, Sdc1 +22.89, Cdkn1a +12.33, Serpine1 +11.21, Eno2 +11.00, Hspa1a +10.81) dominate the entire ranked landscape. Across the full 353-member regulon (coarse 3-class view): 91.8% of the total signed contribution comes from the 340 "other" targets (sum = +351.43); the 9 `shared/glycolytic` members (Hk2, Vegfa, Egln3, Slc2a1, Aldoa, Pgk1, Ldha, Eno1, Pkm) sum to +39.81; and the 4-member `hypoxic HIF core` (Pdk1 −4.16, Bnip3l −1.81, Bnip3 −1.43, Car9 −1.06) is uniformly repressed, summing to −8.46. The score is a regulon-aggregation dominated by non-HIF-diagnostic genes. Vegfa/Slc2a1/Egln3 are `shared/glycolytic`, NOT HIF-specific.**

Reconciling the two views of the same regulon (fig3g and fig3l). fig3g is the coarse three-class view of the regulon — `hypoxic HIF core` (Pdk1/Bnip3/Bnip3l/Car9, repressed, sum = −8.46) / `shared/glycolytic` (the UP angio-glucose targets Vegfa/Slc2a1/Egln3 plus the glycolytic set) / `other`. fig3l is the FINER module decomposition of the SAME regulon: it carves 7 heat-shock genes out of fig3g's 340-member 91.8% "other" lump into a named `heatshock_stress` bucket, leaving a 339-member 90.7% residual "other." "The HIF core" denotes ONE set across both figures — Pdk1/Bnip3/Bnip3l/Car9, sum −8.46. The two figures cannot contradict; the % "other" differs (91.8% coarse vs 90.7% finer) by construction, not by error. Per-gene contribution values are identical between the figures; only class names/membership boundaries differ.

**Method.** The redesigned viz reads the existing `fig3g_target_decomposition_data.csv` (no recomputation of statistics) and renders the full ranked 353-member contribution distribution on a single signed-contribution axis (x = `contrib = sign(mor) × t_wt`). Each gene is plotted at its rank-ordered contribution position, coloured by `class` (`hypoxic HIF core` / `shared/glycolytic` / `other`) with a purple `#762A83` heat-shock/stress highlight overlaid on the six stress contaminants (Timp1, Sdc1, Cdkn1a, Serpine1, Eno2, Hspa1a) to visually connect this landscape to the module-level view in fig3l. The top stress drivers and the curated-class members are text-labelled; the coarse 91.8% / +39.81 / −8.46 provenance split is annotated from the existing `fig3g_target_decomposition_summary.csv`. Zero-line drawn. A viz-local label-subset file (`fig3g_target_landscape_labels.csv`) identifies the labelled genes by rank and class; it is a `head/tail` slice of the sorted table, not a statistic.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf_viz.R` (fig3g block, expanded viz) |
| Function | `contrib = sign(mor) * t_wt` computed in `03_decoupler_tf.R`; coarse 3-class lookup (`hypoxic HIF core` / `shared/glycolytic` / `other`); viz reads pre-computed table |
| Input | `03_results/04_tf/tables/fig3g_target_decomposition_data.csv`, `03_results/04_tf/tables/fig3g_target_decomposition_summary.csv`, `03_results/04_tf/tables/fig3g_target_landscape_labels.csv` |

---

### fig3i_interaction_primer.pdf

**Fever-heat raises the shared HIF/glycolytic target Vegfa in both genotypes — parallel genotype slopes (WT heat logFC = +0.970, KO heat logFC = +0.954, Interaction logFC = +0.016, Interaction adjP ~1.000) signal no detectable cGAS-dependence at n=5 — while the IFN target Ifit1 rises only in WT (WT heat logFC = +0.497, KO heat logFC = −0.511, Interaction logFC = +1.007, Interaction adjP = 0.003; slopes diverge). The phenomenon of heat-responsive gene activation is real and shared across both genotypes; the two arms differ in cGAS-dependence.**

**Method.** Four group-mean log2CPM values per gene (WT/cGASKO × 37 °C/39 °C) drawn from `fig2_marker_means.csv`, plotted as WT (solid) vs cGAS-KO (dashed) across temperature, faceted by gene with free y-axis scales. Contrast logFCs and BH-adjusted Interaction p-value read directly from `02_de_results.rds` and placed in the facet strip. A positive Interaction value means heat raises the gene more in WT than in cGAS-KO (cGAS-dependent); a near-zero Interaction means heat raises the gene similarly in both genotypes (no detectable cGAS-dependence). The Vegfa facet strip now reads "Vegfa [shared HIF/glycolytic target]" (was "[HIF arm]"), because Vegfa is the LEAST HIF-diagnostic target (fig3l: `shared_angio_glucose`); Ifit1's "[IFN arm]" tag is unchanged. Vegfa is the exemplar because it is the highest-t upregulated member of the shared/glycolytic class on the WT_heat contrast (t = 9.39 in fig3g) and its parallel slopes illustrate the absence of a detectable cGAS effect. No new statistics computed.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf_viz.R` |
| Function | Line-plot of pre-computed group means + DE contrast annotations (no recomputation) |
| Input | `03_results/04_tf/tables/fig3i_interaction_primer_data.csv` |

---

### fig3j_topTF_allmethods_WT_heat.pdf

> **RECOMMENDED ARCHIVE / SUPPLEMENT — NOT in primary narrative deck.** fig3c and fig3k already carry method-fragility cleanly on single axes; fig3j's value is the complete per-method audit, which belongs in the supplement. Flagged for Anton's decision.

**The top-12 activated and top-12 suppressed TFs on the WT_heat contrast rendered under all six decoupleR statistics show that the HIF axis is consistently high-ranked under the five univariate-family statistics (Hif1a rank #3–#12: norm_wsum #3, consensus #8, ulm #12, wsum #7, corr_wsum #7) but collapses to rank #142 under MLM, and that the IFN members Irf3/Stat2 are bottom-ranked (ranks ~522–636) across all six methods on the WT_heat contrast, consistent with their heat-agnostic, cGAS-dependent profile.**

**Method.** `decouple(statistics=c("ulm","mlm","wsum"), consensus_score=TRUE)` run on the WT_heat limma t-statistic vector using CollecTRI (minsize=5); top-12 activated and top-12 suppressed TFs per statistic extracted and faceted with free x-axes (absolute scores are not comparable across statistics — only rank and identity are meaningful across facets).

| | |
|---|---|
| Script | `02_analysis/scripts/03b_decoupler_method_comparison_viz.R` |
| Function | `decouple(statistics=c("ulm","mlm","wsum"), consensus_score=TRUE)` on WT_heat t-vector; top/bottom 12 per statistic; facet_wrap by statistic, free x-axis scales |
| Input | `03_results/04_tf/tables/fig3j_allmethods_topTF_data.csv` |

---

### fig3j_topTF_allmethods_Interaction.pdf

> **RECOMMENDED ARCHIVE / SUPPLEMENT — NOT in primary narrative deck.** See note above for fig3j_WT_heat.

**The same six-method top-TF view applied to the Interaction contrast shows the HIF axis uniformly low-ranked across all six statistics (Hif1a: ulm #423, norm_wsum #411, mlm #248, consensus #307, corr_wsum #167, wsum #123), establishing that the IFN-significant / HIF-flat asymmetry is robust to inference method choice.**

**Method.** Identical pipeline as fig3j_topTF_allmethods_WT_heat.pdf, applied to the Interaction t-statistic vector; top-12 activated and top-12 suppressed TFs per statistic extracted and faceted with free x-axes.

| | |
|---|---|
| Script | `02_analysis/scripts/03b_decoupler_method_comparison_viz.R` |
| Function | `decouple(statistics=c("ulm","mlm","wsum"), consensus_score=TRUE)` on Interaction t-vector; top/bottom 12 per statistic; facet_wrap by statistic, free x-axis scales |
| Input | `03_results/04_tf/tables/fig3j_allmethods_topTF_data.csv` |

---

### fig3k_method_rank_divergence.pdf

**MLM is the lone outlier among decoupleR statistics: Spearman rank-correlation vs ULM on the WT_heat contrast is 0.619 for MLM versus 0.968–0.996 for every other statistic (norm_wsum 0.996, corr_wsum 0.975, consensus 0.969, wsum 0.968; n = 658 TFs). Rank instability is an "MLM vs the entire univariate family" multivariate de-confounding effect, not general noise across methods. MLM reshuffles both axes: Rela drops from rank #3 (ULM) to #350 (MLM), Nfkb1 from #6 to #86, and Hif1a from #12 to #142. Hif1a's rank collapse is therefore one instance of a general MLM de-confounding reshuffle of high-collinearity TFs, not a HIF-specific quirk.**

**Method.** Per-statistic TF ranks tabulated for the WT_heat contrast across all 658 TFs surviving minsize=5; Spearman rank-correlation of each statistic vs ULM computed over the full 658-TF rank vector; rank heatmap over a focused TF set with per-method rho annotated on the method axis; Hif1a highlighted; equivalent drops for Rela and Nfkb1 annotated to frame MLM de-confounding as a broad re-ranking of high-collinearity TFs.

| | |
|---|---|
| Script | `02_analysis/scripts/03b_decoupler_method_comparison_viz.R` |
| Function | `cor(rank_ulm, rank_x, method="spearman")` over 658-TF rank vectors; heatmap of per-TF ranks per statistic |
| Input | `03_results/04_tf/tables/fig3k_method_rank_divergence_data.csv`, `03_results/04_tf/tables/fig3k_method_rank_spearman.csv` |

---

### fig3l_hif_attribution.pdf  *(NEW — attribution arc, panel 3 of 3: WHAT KIND)*

**Of the 353 members of the Hif1a CollecTRI regulon, the four curated modules carry 9.3% of the regulon's contribution magnitude and the 339 unassigned members carry 90.7%, so the curated partition accounts for under a tenth of the score and is a reading aid rather than an explanation. Within that tenth: the 7 `heatshock_stress` members (Timp1, Sdc1, Spp1, Cdkn1a, Serpine1, Eno2, Hspa1a) all rise and sum to +108.54, while the four genes diagnostic of canonical HIF1a output — `hif1a_hypoxic_core`: Pdk1, Bnip3, Bnip3l, Car9 — are all four repressed, sum −8.46. The shared angiogenic/glycolytic pair (Vegfa/Slc2a1, +15.93) and the autoregulatory feedback member (Egln3, +6.81) both rise, and neither is diagnostic for HIF1a specifically. The directional result is the one that does not depend on how the buckets were drawn: the genes whose induction is the signature of HIF1a activity move the other way in the contrast that is supposed to induce them. What the remaining 90.7% represents is not established here, and the module sums are not a partition of the score.**

Reconciling the two views of the same regulon (fig3l and fig3g). fig3l is the finer module decomposition of the same 353-member Hif1a regulon shown coarsely in fig3g. "The HIF core" denotes ONE set across both figures — Pdk1/Bnip3/Bnip3l/Car9, sum = −8.46 (repressed). The 7 `heatshock_stress` genes are carved OUT of fig3g's 340-member 91.8% "other" lump, leaving the 339-member 90.7% residual `other_unclassified` here; the % "other" differs (91.8% vs 90.7%) by construction (coarse vs finer bucketing of the identical regulon), not by error. Vegfa/Slc2a1/Egln3 are `shared_angio_glucose`/`autoreg_feedback` (shared/glycolytic), NOT HIF-specific. Per-gene contribution values are unchanged from fig3g; only class names/membership boundaries changed.

**Method.** `03c_hif_program_attribution.R` reads `fig3g_target_decomposition_data.csv` and joins it to a curated module/isoform lookup (module assignment per a curated lookup; definitions and refs in the script header). No statistical recomputation — `t_wt` and `contrib` are copied from the existing fig3g table. Emits `fig3l_hif_attribution_data.csv` (per-gene with `module`, `isoform_attribution`, `direction` columns) and `fig3l_module_summary.csv` (per-module aggregates). The viz renders a horizontal sign-aware lollipop on the signed-contribution (x = `contrib`) axis, faceted by module; the four curated buckets are drawn and `other_unclassified` is left off the panel (its 339 members stay in `fig3l_hif_attribution_data.csv`). Facet order top→bottom: `heatshock_stress` → `shared_angio_glucose` → `autoreg_feedback` → `hif1a_hypoxic_core`, so the module whose members sit left of zero reads last.

Module palette (single source of truth in `config.R::MODULE_COLORS`): `heatshock_stress` = `#762A83` (purple); `shared_angio_glucose` = `#E08214` (orange); `autoreg_feedback` = `#FDB863` (pale orange); `hif1a_hypoxic_core` = `#01665E` (dark teal — repalette; the repressed core is teal, not a third orange, so it separates cleanly from the orange up-set and the purple heat-shock bucket); `other_unclassified` = `grey80`. Bucket **names** come from the companion `config.R::MODULE_LABELS` (added 2026-08-05, same keys as `MODULE_COLORS`), so the human counterpart of this panel can reuse one bucket definition — colour *and* name — through the frozen ortholog map instead of retyping the labels. The viz script appends only each bucket's `n_targets` and `sum_contrib`, read from `fig3l_module_summary.csv`, into that bucket's facet strip.

Drawn text (rebuilt 2026-08-05, canvas 11 × 7.5 in): the title names the panel, the subtitle is two descriptive lines giving what a point is and what x means, and each facet strip names its bucket by curated membership plus that bucket's n and sum. The reading of the pattern lives in the caption below, not on the canvas. A `fits_canvas()` guard in the viz script measures every drawn title/subtitle line against the canvas width and errors out rather than letting `ggsave` emit a silently clipped line — the previous render shipped with its title cut off mid-word and three subtitle lines cut off at the right edge.

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03c_hif_program_attribution.R` |
| Script (viz) | `02_analysis/scripts/03c_hif_program_attribution_viz.R` |
| Function | Join of `fig3g_target_decomposition_data.csv` to curated module lookup; `group_by(module) %>% summarise(...)` aggregations only |
| Input | `03_results/04_tf/tables/fig3g_target_decomposition_data.csv` + curated module/isoform lookup defined in script header |
| Output tables | `03_results/04_tf/tables/fig3l_hif_attribution_data.csv`, `03_results/04_tf/tables/fig3l_module_summary.csv` |

---

### fig3m_ulm_mechanic.pdf  *(NEW — attribution arc, panel 1 of 3: HOW)*

**decoupleR-ULM scores a TF by the regulon-weighted alignment of its targets' t-statistics to their modes of regulation (`t_gene ~ mor_TF`; signed contribution = `sign(mor) × t_wt`). Because the score is a regulon-weighted aggregation, a promiscuous regulon picks up any high-|t| members regardless of biological specificity. Hif1a (353 members, weighted contribution center = +1.08) scores positive in the WT_heat contrast primarily because seven generic heat-shock/stress genes that happen to reside in the CollecTRI regulon (Timp1 +24.92, Sdc1 +22.89, Cdkn1a +12.33, Serpine1 +11.21, Eno2 +11.00, Hspa1a +10.81, Spp1 +15.36) pile the distribution to the right. By contrast, a specific IFN regulon (Stat2, 32 members, weighted center = −0.62) does not pile up on this contrast: its targets are distributed across both sides of zero, reflecting the absence of a coherent heat-MAIN signal for the IFN arm.**

**Method.** `03d_ulm_mechanic.R` reads `fig3g_target_decomposition_data.csv` and copies the Hif1a `contrib` column directly (no recomputation). For the Stat2 comparator regulon, it joins the existing CollecTRI network RDS targets to the WT_heat t-stat vector in `02_de_results.rds` (a descriptive lookup only, no `run_ulm`), computes `aligned_contrib = sign(mor) × t_wt` as an arithmetic relabel of existing t-stats, and flags the seven stress contaminants via `is_stress_contaminant`. Emits `fig3m_ulm_mechanic_data.csv` (353 Hif1a + 32 Stat2 rows) and `fig3m_ulm_mechanic_summary.csv` (per-tf aggregates: n_targets, mean_aligned_contrib, pct_contrib_from_stress, ulm_score_sign). The viz renders two panels sharing the same signed-contribution x-axis (licensing direct comparison): Panel A = Hif1a regulon as a strip/beeswarm with the aggregate weighted center marked as a labelled diamond and the seven stress contaminants highlighted in `#762A83` purple; Panel B = Stat2 regulon on the same axis showing the dispersed, near-zero distribution. DECLUTTER (this round): the score label, the dashed center rule, and the `t_gene ~ mor_TF` formula were moved OFF the swarm — the formula now sits in a footer band — and the Stat2 panel gained an in-panel cue "no pile-up of aligned targets → low score." No statistics changed; the score values, the weighted centers, and the stress-contaminant set are unchanged.

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03d_ulm_mechanic.R` |
| Script (viz) | `02_analysis/scripts/03d_ulm_mechanic_viz.R` |
| Function | Copy of `contrib` from fig3g table (Hif1a); descriptive join of CollecTRI Stat2 targets to WT_heat t-vector + `sign(mor)*t_wt` (no run_ulm); `group_by(tf) %>% summarise(mean_aligned_contrib=mean(...), ...)` |
| Input | `03_results/04_tf/tables/fig3g_target_decomposition_data.csv`, `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/02_de_results.rds` |
| Output tables | `03_results/04_tf/tables/fig3m_ulm_mechanic_data.csv`, `03_results/04_tf/tables/fig3m_ulm_mechanic_summary.csv` |

---

### fig3n_heat_main_regulators.pdf  *(NEW — HSF1 gap)*

**On the heat-main contrast (Temp_main, averaging across both genotypes), Hsf1 is significantly co-elevated (decoupleR-ULM score = +3.20, padj = 0.015) alongside Hif1a (score = +5.14, padj = 2.0×10⁻⁵) and Epas1 (score = +4.17, padj = 8.9×10⁻⁴). Hsf1 was not foregrounded earlier in this analysis. The same direction holds in the genotype-stratified contrasts: Hsf1 is significantly elevated in both WT_heat (score = +2.83, padj = 0.034) and KO_heat (score = +3.48, padj = 0.008), so its co-elevation is present in both genotypes. This figure closes the analytical gap but does NOT claim Hsf1 causes or outranks HIF axis activity; the message is co-elevation and prior neglect, naming HSF1 as an alternative the design can now interrogate. The elevated Hsf1 activity does NOT mean Hif1α's inflated ULM score is an "HSF1 capture": the Hif1α∩Hsf1 regulon overlap is thin (≈24–32 shared targets, ~7–12% of Hif1α's signed contribution), and of the named heat-shock drivers only Hspa1a/Cdkn1a/Serpine1 are shared with Hsf1 while the biggest drivers (Timp1, Sdc1, Spp1, Eno2) are Hif1α-only. The defensible statement is the BROAD one already carried by fig3e — Hif1α's signal is shared across many stress-responsive TFs, not HSF1 specifically — not a localized HSF1 hijack. The Hsf1 cGAS×heat interaction term (raw p = 0.022) is non-significant after BH correction (padj = 0.515) and must not be cited as evidence of cGAS-dependence.**

**Method.** `03e_heat_main_regulators.R` builds a heat-main t-statistic matrix (genes × {Temp_main, WT_heat, KO_heat}) directly from `02_de_results.rds` and runs `run_ulm(mat, net_collectri, .source="source", .target="target", .mor="mor", minsize=5)` — the IDENTICAL decoupleR-ULM call as the primary `03_decoupler_tf.R` — then applies BH within contrast and selects the heat-shock / HIF / IFN axis TFs (`axis`: `heatshock` = Hsf1; `HIF` = Hif1a/Epas1; `IFN` = Irf3/Stat1/Stat2/Nfkb1/Rela). It does NOT filter `master_tf_activities.csv`. fig3n scores are therefore decoupleR-ULM and ARE axis-comparable to fig3a/fig3c/fig3k; they come out byte-identical to the master table's `nes` column only because that column is itself the ULM `score` written under a legacy schema name (`nes = score` in `03_decoupler_tf.R`). The viz renders a lollipop on the decoupleR-ULM score axis coloured by axis (heatshock = `#762A83` purple, matching the stress highlight in fig3g/fig3l for cross-figure consistency; HIF = `#B35806`; IFN = `#2166AC`), grouped by contrast, with padj < 0.05 stars matching fig3a's convention.

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03e_heat_main_regulators.R` |
| Script (viz) | `02_analysis/scripts/03e_heat_main_regulators_viz.R` |
| Function | `run_ulm(.source="source", .target="target", .mor="mor", minsize=5)` on a fresh heat-MAIN t-stat matrix + BH within contrast — same call as `03_decoupler_tf.R` |
| Input | `03_results/objects/02_de_results.rds`, `03_results/objects/net_collectri_mouse.rds` |
| Output table | `03_results/04_tf/tables/fig3n_heat_main_regulators_data.csv` |

---

### fig3p_heatmain_ranking.pdf  *(NEW — non-identifiability cluster, 1 of 3; deck placement PENDING)*

**On the heat-main contrast (Temp_main) there is no clean winner. Hif1a is only #9 of 658 scored TFs (ULM score = 5.14, p = 2.8×10⁻⁷), embedded in a dense crowd of generic stress / immediate-early / NF-κB regulators separated by tiny gaps — Jun #1 (6.06), Egr1 #2 (5.54), Fos #3 (5.46), Sp1 #4 (5.29), Stat5a #5 (5.27), Rela #6 (5.25), Ncoa1 #7 (5.25), Hoxd3 #8 (5.20), then Hif1a #9 (5.14), Nfkb1 #10 (5.12), Jund #11 (5.10), Nfkb #12 (4.81). The canonical heat-shock TF Hsf1 sits far down at #50 (3.20). Hif1a's nine-place position among interchangeable promiscuous regulators, and the gap of only ~0.9 score units spanning ranks #1–#12, mean this contrast does not nominate Hif1a (or anything else) as a unique driver.**

**Method.** `03g_nonidentifiability.R` runs `run_ulm(net_collectri, .source="source", .target="target", .mor="mor", minsize=5)` on the Temp_main limma t-statistic vector from `02_de_results.rds` (the same CollecTRI-ULM call as `03_decoupler_tf.R`), scoring 658 TFs and ranking by descending score; raw p-values stored. The viz renders a horizontal lollipop of the top ~22 TFs by score with Hsf1 appended below a visual break, each TF coloured by curated family via the shared `FAMILY_COLORS` palette (HIF axis = orange; AP-1/immediate-early = red; NF-κB = blue; heat-shock = purple `#762A83`; housekeeping/STAT/etc.). Hif1a's #9 placement is read directly off the ranked axis; no per-TF crowning — the figure's message is the crowd, not any member.

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03g_nonidentifiability.R` |
| Script (viz) | `02_analysis/scripts/03g_nonidentifiability_viz.R` |
| Function | `run_ulm(.source="source", .target="target", .mor="mor", minsize=5)` on the Temp_main t-stat vector; `arrange(desc(score)) %>% mutate(rank = row_number())` |
| Input | `03_results/04_tf/tables/fig3p_heatmain_ranking_data.csv` (from `03_results/objects/02_de_results.rds` + `net_collectri_mouse.rds`) |

---

### fig3q_coregulators.pdf  *(NEW — non-identifiability cluster, 2 of 3; deck placement PENDING)*

**92% of Hif1a's 353 CollecTRI targets (≈325/353; mean 22 other TFs per target — fig3e) are also targets of the network's most promiscuous, non-hypoxia-specific regulators, so the heat-MAIN signal cannot be attributed to Hif1a alone. The largest target-sharers are generic stress/proliferation/inflammation factors, not hypoxia specialists: Sp1 shares 171/353 (48.4%), Trp53 143 (40.5%), the NF-κB family (Nfkb 127 = 36.0%, Rela 102 = 28.9%, Nfkb1 81 = 22.9%), the AP-1/immediate-early family (Jun 118 = 33.4%, Ap1 100 = 28.3%, Fos 66 = 18.7%), Myc 110 (31.2%), nuclear receptors (Esr1 96 = 27.2%, Ar 80 = 22.7%), and STATs (Stat3 90 = 25.5%). For reference the other HIF-axis member Epas1 shares only 46 (13.0%) and heat-shock Hsf1 only 24 (6.8%). This shared-target collinearity is the regulon basis for the ULM 5.11/#12 → MLM 1.13/#142 collapse documented in fig3e — the same genes credited to Hif1a are competed for by dozens of equally plausible regulators.**

**Method.** `03g_nonidentifiability.R` performs set-membership counting of every CollecTRI regulon against the fixed 353-member Hif1a target set: for each other TF, `shared_targets` = |its regulon ∩ Hif1a's 353 targets| and `pct_of_hif1a_set` = shared_targets / 353 × 100. No statistical inference — pure set arithmetic over CollecTRI edges. The viz renders the top sharers as horizontal bars (length = shared/353), coloured by curated family via `FAMILY_COLORS`. The figure foregrounds that the heaviest sharers are non-hypoxia-specific (Sp1/Trp53/NF-κB/AP-1/Myc), establishing shared ownership of the target set rather than crowning any of them as the driver.

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03g_nonidentifiability.R` |
| Script (viz) | `02_analysis/scripts/03g_nonidentifiability_viz.R` |
| Function | Set-membership counting: `count(target)` of each CollecTRI regulon ∩ Hif1a's 353-target set; `pct_of_hif1a_set = shared_targets / 353 * 100` |
| Input | `03_results/04_tf/tables/fig3q_coregulators_data.csv` (from `03_results/objects/net_collectri_mouse.rds`) |

---

### fig3r_shared_ownership.pdf  *(NEW — non-identifiability cluster, 3 of 3; deck placement PENDING)*

**The same heat-driven genes populate many TFs' regulons, so the contrast cannot single out Hif1a — or any one TF — as the driver. The top-20 Hif1a targets by signed heat-MAIN contribution (Timp1 contrib = 24.92, Sdc1 22.89, Itga5 18.46, Glb1 17.39, Lgals3 16.80, Bcl2l1 16.08, Spp1 15.36, … down to Rora 10.95) are each also members of several of Hif1a's 12 largest target-sharers' regulons (Sp1, Trp53, Nfkb, Jun, Myc, Rela, Ap1, Esr1, Stat3, Egr1, Nfkb1, Ar — the heaviest sharers from fig3q). Because these high-|t| genes are jointly owned, the heat-MAIN ULM score they generate is claimed simultaneously by Hif1a and by a dozen generic regulators that all score in the same #1–#12 band (fig3p). This is the unifier of fig3e (targets shared) and fig3p (Hif1a only #9 among those same sharers): non-identifiability, not attribution.**

**Method.** `03g_nonidentifiability.R` selects the top-20 Hif1a targets by signed contribution (`gene_contrib`, copied from fig3g — `sign(mor) × t_wt` on the heat-MAIN t-stat) and crosses them against {Hif1a + its 12 largest target-sharers from fig3q}; `in_regulon` is the CollecTRI set-membership indicator (1/0) for each gene×TF cell. Per-gene `gene_heat_t` (heat-MAIN limma t-statistic) and per-TF `tf_heatmain_score` (the fig3p heat-MAIN ULM score) are carried for the marginals. The viz renders a tile heatmap (tile = CollecTRI membership), with a LEFT marginal showing each gene's heat-MAIN t-statistic and a TOP marginal showing each TF's heat-MAIN ULM score; Hif1a is flagged orange (HIF-axis `FAMILY_COLORS`) but not privileged in the ranking — it is one column among thirteen. No statistics recomputed; all values are joins/copies of fig3g + fig3p.

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03g_nonidentifiability.R` |
| Script (viz) | `02_analysis/scripts/03g_nonidentifiability_viz.R` |
| Function | Cross of top-20 Hif1a targets (by `gene_contrib` from fig3g) × {Hif1a + 12 largest sharers}; `in_regulon` = CollecTRI membership; marginals = `gene_heat_t`, `tf_heatmain_score` (no run_ulm) |
| Input | `03_results/04_tf/tables/fig3r_membership_data.csv` (from `net_collectri_mouse.rds` + `02_de_results.rds` + fig3g/fig3p tables) |

---

## Table Legends

---

### fig3a_tf_interaction_axes_data.csv

**All TFs represented on fig3a: ULM scores and BH-adjusted p-values for the Interaction contrast, with axis classification and significance flag. Key values: Irf3 score = +5.93 padj = 1.04×10⁻⁶ (sig); Stat2 +5.85 padj = 1.10×10⁻⁶ (sig); Stat1 +5.30 padj = 1.88×10⁻⁵ (sig); Nfkb1 +2.16 padj = 0.598 (NS); Rela +2.03 padj = 0.631 (NS); Hif1a −0.19 padj = 0.97 (NS); Epas1 −0.86 padj = 0.91 (NS).**

29 rows × 7 columns. `source` = TF gene symbol; `score` = CollecTRI ULM activity score on the Interaction contrast; `padj` = BH-adjusted p-value within Interaction contrast; `p_value` = raw p-value; `axis` = HIF | IFN | other; `sig` = logical, padj < 0.05; `direction` = Up | Down.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | `run_ulm(.mor="mor", minsize=5)` + `p.adjust("BH")` within Interaction contrast |
| Input | `03_results/objects/03_tf_collectri.rds` |

---

### fig3b_top_tf_by_contrast_data.csv

**Top-12 activated and top-12 suppressed TFs (by ULM score) for each of four contrasts (WT_heat, KO_heat, Temp_main, Interaction), with axis and key-TF annotation.**

96 rows × 9 columns. `contrast` = one of WT_heat | KO_heat | Temp_main | Interaction; `source` = TF symbol; `score` = ULM score; `padj` = BH-adjusted p-value within contrast; `direction` = Up | Down; `rank` = rank within direction within contrast (1 = most extreme); `axis` = HIF | IFN | other; `key_tf` = logical; `sig` = logical, padj < 0.05.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | `run_ulm(.mor="mor", minsize=5)` + BH; top/bottom `TOPN_CONTRAST = 12` per contrast |
| Input | `03_results/objects/03_tf_collectri.rds` |

---

### fig3c_hif1a_rank_cascade_data.csv

**Hif1a's rank, score, number of TFs assessed, and percentile rank across four inference method/network combinations on the WT_heat contrast.**

4 rows × 5 columns. `method` = DoRothEA_ULM | CollecTRI_ULM | CollecTRI_MLM | CollecTRI_consensus; `rank` = integer rank by descending score (1 = most activated); `n_tfs` = total TFs surviving minsize=5; `score` = Hif1a activity score; `pct_rank` = 100 × rank / n_tfs.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | `run_ulm` / `run_mlm` / `decouple(..., consensus_score=TRUE)` on WT_heat; `row_number()` after `arrange(desc(score))` |
| Input | `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/net_dorothea_mouse_ABC.rds`, `03_results/objects/02_de_results.rds` |

---

### fig3d_regulon_swap_data.csv

**Per-target t-statistics and network membership for all genes in the union of the DoRothEA and CollecTRI Hif1a regulons present in the expression universe.**

384 rows × 7 columns. `target` = gene symbol; `in_dorothea` = logical; `in_collectri` = logical; `t_wt` = WT_heat limma t-statistic; `membership` = both | dorothea_only | collectri_only; `abs_t` = |t_wt|; `high_t` = logical, |t_wt| ≥ 2.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | Set intersection/difference of filtered network target lists; `mat[target, "WT_heat"]` |
| Input | `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/net_dorothea_mouse_ABC.rds`, `03_results/objects/02_de_results.rds` |

---

### fig3d_regulon_swap_summary.csv

**Per-network summary statistics for the Hif1a regulon t-distribution, quantifying the dilution of signal when switching from DoRothEA (mean|t| = 4.70) to CollecTRI (mean|t| = 4.26).**

2 rows × 8 columns. `network` = DoRothEA | CollecTRI; `n` = number of targets present in universe; `mean_abs_t`; `median_t`; `pct_abs_t_gt2`; `n_intersection` = 101; `n_dorothea_only` = 30; `n_collectri_only` = 252.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | `net_t_summary()` helper: `mean(abs(tv))`, `median(tv)`, `mean(abs(tv)>2)` |
| Input | Derived from `fig3d_regulon_swap_data.csv` |

---

### fig3e_mlm_collinearity_data.csv

**Per-target count of co-regulating TFs and WT_heat t-statistic for all 353 Hif1a CollecTRI targets, demonstrating that 92% are shared with at least one other TF regulon (median co-regulator count = 11, mean ~22).**

353 rows × 4 columns. `target` = gene symbol; `n_other_TFs` = number of other CollecTRI TFs also targeting this gene (Hif1a excluded); `t_wt` = WT_heat limma t-statistic; `co_regulated` = logical, n_other_TFs ≥ 1.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | `distinct(source, target) %>% count(target, name="n_tfs_total")`; `n_other_TFs = n_tfs_total − 1` |
| Input | `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/02_de_results.rds` |

---

### fig3e_mlm_collinearity_summary.csv

**Regulon-level aggregates: 92.1% of Hif1a's 353 CollecTRI targets co-regulated by ≥1 other TF; mean 22.2 (median 11) other regulators per target.**

1 row × 4 columns. `n_targets` = 353; `frac_co_regulated`; `mean_other_TFs`; `median_other_TFs`.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | `summarise(frac_co_regulated = mean(n_other_TFs ≥ 1), ...)` |
| Input | `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/02_de_results.rds` |

---

### fig3e_score_collapse.csv

**Two-row table recording Hif1a's ULM score (5.11, rank #12) and MLM score (1.13, rank #142) on the WT_heat CollecTRI run, quantifying the 4.5-fold score collapse attributable to multivariate de-confounding.**

2 rows × 4 columns. `method` = ULM | MLM; `score`; `rank`; `n_tfs`.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | `run_ulm(.mor="mor", minsize=5)` and `run_mlm(.mor="mor", minsize=5)` on WT_heat vector |
| Input | `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/02_de_results.rds` |

---

### fig3f_consensus_zmix_data.csv

**Per-statistic folded-z scores and the consensus score for Hif1a (WT_heat, CollecTRI): ulm z = 2.44, wsum z = 4.13, norm_wsum z = 2.72, corr_wsum z = 4.31, mlm z = 0.85, consensus z = 2.89.**

6 rows × 5 columns. `source` = "Hif1a"; `statistic` = ulm | mlm | wsum | norm_wsum | corr_wsum | consensus; `z_score`; `family` = univariate | multivariate | consensus; `canonical` = logical.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | `decouple(statistics=c("ulm","mlm","wsum"), consensus_score=TRUE)`; folded-z via sign-symmetric sd standardisation; `consensus = mean(z_score)` verified vs stored value (diff < 1×10⁻⁶) |
| Input | `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/02_de_results.rds` |

---

### fig3g_target_decomposition_data.csv

**Per-target signed contributions (`contrib = sign(mor) × t_wt`) for all 353 Hif1a CollecTRI targets in the WT_heat expression universe, classified by biological specificity. This is the primary input table for all three attribution arc figures (fig3m, fig3g, fig3l) — no arc figure recomputes any statistics from it.**

353 rows × 6 columns. `source` = "Hif1a"; `target` = gene symbol; `mor` = mode of regulation (±1, from CollecTRI); `t_wt` = WT_heat limma t-statistic; `contrib` = sign(mor) × t_wt; `class` = coarse 3-class label `hypoxic HIF core` (4: Pdk1/Bnip3/Bnip3l/Car9) | `shared/glycolytic` (9: Hk2/Vegfa/Egln3/Slc2a1/Aldoa/Pgk1/Ldha/Eno1/Pkm) | `other` (340). fig3l decomposes the `other` lump finer (carving out 7 heat-shock genes).

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | `contrib = sign(mor) * t_wt`; coarse-class lookup against the `hypoxic HIF core` and `shared/glycolytic` member sets |
| Input | `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/02_de_results.rds` |

---

### fig3g_target_decomposition_summary.csv

**Three-row coarse-class aggregation of signed contributions: `other` 91.8% (n = 340, sum = +351.43), `shared/glycolytic` (n = 9, sum = +39.81, 10.4% of total), `hypoxic HIF core` (n = 4, sum = −8.46, −2.21% of total — repressed). The `shared/glycolytic` and `hypoxic HIF core` percentages sum against the positive total, hence the negative core %.**

3 rows × 5 columns. `class` = other | shared/glycolytic | hypoxic HIF core; `n_targets`; `sum_contrib`; `mean_contrib`; `pct_of_total`.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | `group_by(class) %>% summarise(n_targets=n(), sum_contrib=sum(contrib), ...)` |
| Input | Derived from `fig3g_target_decomposition_data.csv` |

---

### fig3g_target_landscape_labels.csv

**Viz-local label subset: the genes text-labelled in the redesigned fig3g ranked contribution landscape (the curated-class members — `hypoxic HIF core` and `shared/glycolytic` — plus the highest stress-contaminant entries). A `head/tail` slice of `fig3g_target_decomposition_data.csv` sorted by `contrib`, augmented with `rank` and `legend_class` columns for plot annotation. Not a statistic — a convenience slice for the viz script.**

8 columns. `source`, `target`, `mor`, `t_wt`, `contrib`, `class` (from upstream table); `rank` = integer position in the 353-gene contribution ranking (1 = most negative); `legend_class` = re-label for the legend (curated-class / heat-shock-stress).

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf_viz.R` (computed inline as a viz-local slice) |
| Function | `arrange(contrib) %>% mutate(rank = row_number())` on the upstream table; `head()`/`tail()` to select label candidates |
| Input | `03_results/04_tf/tables/fig3g_target_decomposition_data.csv` |

---

### fig3h_lombardi_vs_phylo_data.csv

**ULM scores, p-values, BH-adjusted p-values, and significance flags for two HIF signatures (Lombardi-100 and Phylo-16) across three contrasts (WT_heat, KO_heat, Interaction). Both signatures heat-activated in both genotypes (padj ≤ 0.015); Interaction NS for both (Lombardi padj = 0.34, Phylo padj = 0.38) — no detectable cGAS-dependence at n=5 for either signature.**

6 rows × 9 columns. `signature`; `condition`; `score`; `p_value`; `padj`; `n_genes`; `sig_label`; `sig` = logical, padj < 0.05; `n_lombardi` = 100; `n_phylo` = 16.

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | Unit-mor signature networks scored by `run_ulm(.mor="mor", minsize=5)` across 7-contrast matrix; BH within contrast |
| Input | `00_data/references/gene_sets/lombardi2022_hif_consensus_mouse.rds`, Phylo-16 hard-coded vector, `03_results/objects/02_de_results.rds` |

---

### fig3i_interaction_primer_data.csv

**The eight group-mean log2CPM values (Ifit1 and Vegfa × genotype × temperature) with the three pre-computed contrast logFCs and Interaction BH-adjusted p-value visualised in fig3i.**

8 rows × 9 columns. `gene`; `arm` = cGAS-dependent | no-detectable-cGAS-dependence; `genotype` = WT | cGASKO; `temp` = 37 | 39 (°C); `mean_log2cpm`; `wt_heat` = WT heat logFC (Ifit1: +0.497; Vegfa: +0.970); `ko_heat` = cGASKO heat logFC (Ifit1: −0.511; Vegfa: +0.954); `interaction` = Interaction logFC (Ifit1: +1.007; Vegfa: +0.016); `inter_adjP` = BH-adjusted Interaction p-value (Ifit1: 0.003; Vegfa: ~1.000).

| | |
|---|---|
| Script | `02_analysis/scripts/03_decoupler_tf.R` |
| Function | Group means joined from `fig2_marker_means.csv`; contrast values extracted from `02_de_results.rds` |
| Input | `03_results/objects/02_de_results.rds`, `03_results/03_de/tables/fig2_marker_means.csv` |

---

### fig3j_allmethods_topTF_data.csv

**Top-12 activated and top-12 suppressed TFs per contrast × per decoupleR statistic (4 contrasts × 6 statistics = 576 rows); tidy source for the fig3j faceted lollipop figures. Scores are NOT comparable across statistics — only rank and identity are meaningful across facets.**

576 rows × 8 columns. `contrast`; `statistic` = ulm | mlm | wsum | norm_wsum | corr_wsum | consensus; `source`; `score`; `rank`; `direction`; `axis`; `key_tf`.

| | |
|---|---|
| Script | `02_analysis/scripts/03b_decoupler_method_comparison.R` |
| Function | `decouple(statistics=c("ulm","mlm","wsum"), consensus_score=TRUE)` on each contrast; top/bottom 12 per statistic per contrast |
| Input | `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/02_de_results.rds` |

---

### fig3k_method_rank_divergence_data.csv

**Per-TF ranks across all six decoupleR statistics for WT_heat and Interaction contrasts, restricted to a focused TF set for the rank-divergence heatmap. Key values (WT_heat): Rela ulm #3 → mlm #350; Nfkb1 ulm #6 → mlm #86; Hif1a ulm #12 → mlm #142; Irf3/Stat2 bottom-ranked in all methods on WT_heat.**

348 rows × 5 columns. `contrast`; `source`; `statistic`; `rank` (among all 658 TFs surviving minsize=5, 1 = most activated); `score`.

| | |
|---|---|
| Script | `02_analysis/scripts/03b_decoupler_method_comparison.R` |
| Function | `decouple(statistics=c("ulm","mlm","wsum"), consensus_score=TRUE)`; `row_number()` after `arrange(desc(score))` per statistic per contrast; filtered to focused TF set |
| Input | `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/02_de_results.rds` |

---

### fig3k_method_rank_spearman.csv

**Spearman rank-correlation of each decoupleR statistic vs ULM over all 658 TFs on the WT_heat contrast: MLM rho = 0.619 (lone outlier); all others rho = 0.968–0.996 (norm_wsum 0.996; corr_wsum 0.975; consensus 0.969; wsum 0.968).**

6 rows × 3 columns. `statistic`; `spearman_vs_ulm`; `n_tfs` = 658 for all rows.

| | |
|---|---|
| Script | `02_analysis/scripts/03b_decoupler_method_comparison.R` |
| Function | `cor(rank_ulm, rank_x, method="spearman")` over 658-TF rank vectors per statistic on WT_heat |
| Input | Derived from full-rank tables in `03_results/objects/03_tf_collectri.rds` and method-comparison decouple run |

---

### fig3l_hif_attribution_data.csv

**Per-target signed contributions for the 15 curated Hif1a regulon members (7 heatshock_stress + 2 shared_angio_glucose + 1 autoreg_feedback + 4 hif1a_hypoxic_core) plus a pointer to the 339 other_unclassified members. Values are copied directly from `fig3g_target_decomposition_data.csv` — no recomputation.**

15 curated rows (see note on other_unclassified) × 7 columns. `source` = "Hif1a"; `target`; `t_wt` = WT_heat t-statistic (copied from fig3g); `contrib` = sign(mor)×t_wt (copied from fig3g); `module` = heatshock_stress | shared_angio_glucose | autoreg_feedback | hif1a_hypoxic_core | other_unclassified; `isoform_attribution` = non-HIF (stress) | shared HIF1/HIF2 | HIF1a-preferential-feedback | HIF1a-selective | unclassified; `direction` = Up | Down.

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03c_hif_program_attribution.R` |
| Function | Left join of `fig3g_target_decomposition_data.csv` to curated module lookup; no statistical functions |
| Input | `03_results/04_tf/tables/fig3g_target_decomposition_data.csv` |

---

### fig3l_module_summary.csv

**Per-module aggregation of signed contributions. Two percentage columns, and they answer different questions — quote the one you mean. `pct_of_total_magnitude` is each module's share of the regulon's total contribution magnitude; it sums to 100 and it is the column to use when the claim is how much of the score a module accounts for: other_unclassified 90.7%, heatshock_stress 7.2%, shared_angio_glucose 1.1%, hif1a_hypoxic_core 0.6%, autoreg_feedback 0.5%. `pct_of_total` keeps the signed convention `fig3g_target_decomposition_summary` uses — a net module sum over the same gross denominator — so its sign tracks `sum_contrib` and it reads −0.56% for the repressed core. It does not sum to 100, and it understates any module whose members cancel: other_unclassified is 189 up against 150 down, so it reads 17.3% there against a true magnitude share of 90.7%. The two columns are identical for every module of uniform sign, which is why they diverge on exactly the two modules where the distinction matters. Sums: heatshock_stress +108.54 (7 genes, all up); shared_angio_glucose +15.93 (2, all up); autoreg_feedback +6.81 (1, up); hif1a_hypoxic_core −8.46 (4, all down); other_unclassified +259.97 (339, 189 up / 150 down).**

5 rows × 6 columns. `module`; `n_targets`; `sum_contrib`; `pct_of_total`; `n_up`; `n_down`.

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03c_hif_program_attribution.R` |
| Function | `group_by(module) %>% summarise(n_targets=n(), sum_contrib=sum(contrib), pct_of_total=..., n_up=sum(direction=="Up"), n_down=sum(direction=="Down"))` |
| Input | `03_results/04_tf/tables/fig3l_hif_attribution_data.csv` |

---

### fig3m_ulm_mechanic_data.csv

**Per-target aligned contributions for the Hif1a (n = 353) and Stat2 (n = 32) regulons, with a flag identifying the seven stress contaminants within the Hif1a regulon. The primary data source for the algorithmic mechanic figure (fig3m). No statistics computed — Hif1a `contrib` values are copied from fig3g; Stat2 `aligned_contrib = sign(mor) × t_wt` is an arithmetic relabel of existing t-statistics (descriptive lookup only).**

385 rows × 7 columns. `tf`; `regulon_class` = promiscuous | specific; `target`; `mor` = ±1; `t_wt` = WT_heat limma t-statistic; `aligned_contrib` = sign(mor) × t_wt; `is_stress_contaminant` = logical (TRUE for {Hspa1a, Timp1, Sdc1, Cdkn1a, Serpine1, Eno2, Spp1} within the Hif1a regulon).

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03d_ulm_mechanic.R` |
| Function | Copy of `contrib` from fig3g table (Hif1a); join of CollecTRI Stat2 targets to WT_heat t-vector + arithmetic `sign(mor)*t_wt` (no run_ulm) |
| Input | `03_results/04_tf/tables/fig3g_target_decomposition_data.csv`, `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/02_de_results.rds` |

---

### fig3m_ulm_mechanic_summary.csv

**Per-TF regulon-level summary: Hif1a (n_targets = 353, mean_aligned_contrib = +1.08, pct_contrib_from_stress = 7.2%, ulm_score_sign = positive); Stat2 (n_targets = 32, mean_aligned_contrib = −0.62, pct_contrib_from_stress = NA, ulm_score_sign = negative).**

2 rows × 6 columns. `tf`; `regulon_class`; `n_targets`; `mean_aligned_contrib`; `pct_contrib_from_stress`; `ulm_score_sign`.

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03d_ulm_mechanic.R` |
| Function | `group_by(tf) %>% summarise(mean_aligned_contrib = mean(aligned_contrib), ...)` |
| Input | `03_results/04_tf/tables/fig3m_ulm_mechanic_data.csv` |

---

### fig3n_heat_main_regulators_data.csv

**decoupleR-ULM scores and BH-adjusted p-values for the heat-shock / HIF / IFN axis TFs {Hsf1, Hif1a, Epas1, Irf3, Stat1, Stat2, Nfkb1, Rela} on contrasts {Temp_main, WT_heat, KO_heat}, used to build fig3n. ESTIMATOR NOTE: scores are decoupleR-ULM — the SAME estimator as fig3a/fig3c/fig3k (`run_ulm(.mor="mor", minsize=5)`), computed on a fresh heat-MAIN t-stat matrix — and ARE axis-comparable to those figures. (They equal the `master_tf_activities.csv` `nes` column only because that column is itself the ULM `score` under a legacy schema name; `nes` is a misnomer, not a GSEA score.) Key values: Hsf1 Temp_main +3.20 padj = 0.015 (sig); WT_heat +2.83 padj = 0.034 (sig); KO_heat +3.48 padj = 0.008 (sig); Hif1a Temp_main +5.14 padj = 2.0×10⁻⁵; Epas1 Temp_main +4.17 padj = 8.9×10⁻⁴.**

24 rows × 6 columns. `tf`; `contrast` = Temp_main | WT_heat | KO_heat; `score` = decoupleR-ULM activity score; `padj` = BH-adjusted p-value within contrast; `direction` = Up | Down; `axis` = heatshock | HIF | IFN.

| | |
|---|---|
| Script (compute) | `02_analysis/scripts/03e_heat_main_regulators.R` |
| Function | `run_ulm(.source="source", .target="target", .mor="mor", minsize=5)` on a fresh heat-MAIN t-stat matrix + BH within contrast + axis lookup — same call as `03_decoupler_tf.R` |
| Input | `03_results/objects/02_de_results.rds`, `03_results/objects/net_collectri_mouse.rds` |

---

## Reproduction

```bash
# From project root — in pipeline order:
Rscript 02_analysis/scripts/03_decoupler_tf.R
Rscript 02_analysis/scripts/03_decoupler_tf_viz.R
Rscript 02_analysis/scripts/03b_decoupler_method_comparison.R
Rscript 02_analysis/scripts/03b_decoupler_method_comparison_viz.R
Rscript 02_analysis/scripts/03c_hif_program_attribution.R
Rscript 02_analysis/scripts/03c_hif_program_attribution_viz.R
Rscript 02_analysis/scripts/03d_ulm_mechanic.R
Rscript 02_analysis/scripts/03d_ulm_mechanic_viz.R
Rscript 02_analysis/scripts/03e_heat_main_regulators.R
Rscript 02_analysis/scripts/03e_heat_main_regulators_viz.R
```

All COMPUTE scripts (`03_`, `03b_`, `03c_`, `03d_`, `03e_`) contain no `ggplot`/`ggsave`. All VIZ scripts (`_viz.R`) contain no `lmFit`/`eBayes`/`run_ulm`/`run_mlm`/`p.adjust`/`prcomp`. Networks pre-built by `02_analysis/scripts/00c_prepare_networks.R` — do NOT call `decoupleR::get_collectri()` or `get_progeny()` in this environment (OmniPath unreachable; see `01_modules/SciAgent-toolkit/skills/bulk-rnaseq-activity-inference/references/known-issues.md`).

## Dependencies

- R with packages: decoupleR, limma, dplyr, ggplot2, readr
- Pre-built network objects: `03_results/objects/net_collectri_mouse.rds`, `03_results/objects/net_dorothea_mouse_ABC.rds`
- Upstream objects: `03_results/objects/02_de_results.rds`, `03_results/master/master_tf_activities.csv`
- Reference gene sets: `00_data/references/gene_sets/lombardi2022_hif_consensus_mouse.rds`
- Config: `02_analysis/config/config.R` (palette, theme, `provisional_caption()`)

## Notes

- **Sample mapping**: every figure reads `sample_mapping_caption()` from `design.sample_mapping` in the config, which records the mapping as owner-confirmed on 2026-07-22, 20 of 20 libraries concordant with the label-blind marker call.
- **cGAS framing**: the HIF arm shows "no detectable cGAS-dependence at n=5" — this is an underpowered n=5 null, not evidence of independence. Never rephrase as "cGAS-independent," "parallel," or "cascade."
- **Object name**: the program revealed here is "a heat-induced glycolytic/stress program partially overlapping HIF targets." Never "the HIF program," "HIF1a the TF," or "canonical HIF1a program."
- **HIF2a/Epas1**: elevated on all heat contrasts (Temp_main decoupleR-ULM score = +4.17, padj = 8.9×10⁻⁴) but named only as an alternative the HIF inhibitor cannot exclude — not crowned as a driver.
- **fig3n estimator**: fig3n scores are decoupleR-ULM — the SAME `run_ulm(.mor="mor", minsize=5)` estimator as fig3a/fig3c/fig3k, computed on a fresh heat-MAIN t-stat matrix — and ARE axis-comparable to those figures. They equal the `master_tf_activities.csv` `nes` column only because that column is the ULM `score` written under a legacy schema name (`nes = score`); the `nes` name is a misnomer kept for master-schema compatibility, NOT a GSEA score. The Hsf1 cGAS×heat interaction raw p (0.022) is non-significant after BH correction (padj = 0.515) and must not be cited as evidence of cGAS-dependence.
- **fig3d + fig3e in the deck spine**: they carry how the collaborator got #1 and why it does not hold, regulon dilution (#1→#12) and collinearity redistribution (#12→#142), and they sit between fig3l and fig3c.
- **Archive figures (not in the deck; canonical PDFs/tables retained)**: fig3b, fig3f, fig3j×2, and **lombardi_recurrence** are forensic/supplement backups (fig3b's broad-program point is made more sharply by fig3g/fig3l; fig3f risks the "consensus recovers HIF1a" mis-read; fig3j belongs in the supplement as the full per-method audit). Canonical artifacts kept, just not on the deck spine.

## tables/_overview/fig3a_tf_interaction_axes.csv

Plot-ready TF activity data for fig3a written by `save_overview()`: CollecTRI ULM scores and BH-adjusted p-values for the Interaction contrast panel (29 TFs), factor-ordered by score for lollipop rendering.

**How to read:** Each row = one TF. `source` = TF gene symbol; `score` = ULM activity on the Interaction contrast (Interaction = WT_heat − KO_heat; positive = more active when cGAS is present); `padj` = BH-adjusted p-value within the Interaction contrast; `axis` = HIF / IFN / other; `sig` = logical padj < 0.05; `direction` = Up / Down. Rows ordered by ascending score for bottom-to-top lollipop layout. Key values: Irf3 +5.93 padj = 1.04×10⁻⁶ (sig); Hif1a −0.19 padj = 0.97 (NS). Claim tier L3 (n=5, underpowered null).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `colors.diverging; design.axis_colors` | `03_results/04_tf/tables/fig3a_tf_interaction_axes_data.csv` |

## tables/_overview/fig3b_top_tf_by_contrast.csv

Plot-ready data for fig3b: top-N activated and suppressed TFs (CollecTRI ULM) per contrast, filtered to at most `figures.top_n` = 20 per direction per contrast by |score|, with a `row_key` integer for within-facet ordering.

**How to read:** Each row = one TF × contrast combination. `contrast` = WT_heat / KO_heat / Temp_main / Interaction; `source` = TF symbol; `score` = ULM; `padj` = BH within contrast; `direction` = Up / Down; `rank` = rank within direction within contrast; `axis` = HIF / IFN / other; `key_tf` = logical (HIF/IFN watchlist); `sig` = logical padj < 0.05; `row_key` = plot ordering integer. Claim tier L3 (n=5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `figures.top_n=20; design.axis_colors` | `03_results/04_tf/tables/fig3b_top_tf_by_contrast_data.csv` |

## tables/_overview/fig3c_hif1a_rank_cascade.csv

Plot-ready data for fig3c: Hif1a rank, score, n_tfs, and pct_rank across four inference method/network combinations on the WT_heat contrast, documenting the #1 → #12 → #142 → #8 cascade.

**How to read:** 4 rows × 5 columns. `method` = DoRothEA_ULM / CollecTRI_ULM / CollecTRI_MLM / CollecTRI_consensus; `rank` = integer rank among all scored TFs by descending score (1 = most activated); `n_tfs` = total TFs surviving minsize=5; `score` = Hif1a activity estimate; `pct_rank` = 100 × rank / n_tfs. Read left to right as the cascade: the rank is not a stable biological signal, it is a method/network artefact. Claim tier L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `colors.diverging` | `03_results/04_tf/tables/fig3c_hif1a_rank_cascade_data.csv` |

## tables/_overview/fig3d_regulon_swap.csv

Plot-ready data for fig3d: per-target WT_heat limma t-statistics for all genes in the union of the DoRothEA and CollecTRI Hif1a regulons, classified by network membership.

**How to read:** Each row = one gene (384 total in union). `target` = gene symbol; `in_dorothea` / `in_collectri` = logical membership; `t_wt` = WT_heat limma t-statistic; `membership` = both (n = 101) / dorothea_only (n = 30) / collectri_only (n = 252); `abs_t` = |t_wt|; `high_t` = logical |t_wt| ≥ 2. The CollecTRI-only targets straddle zero (mixed sign), diluting the mor-weighted ULM mean and explaining the #1 → #12 rank drop. Claim tier L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `design.module_colors` | `03_results/04_tf/tables/fig3d_regulon_swap_data.csv` |

## tables/_overview/fig3e_mlm_collinearity.csv

Plot-ready data for fig3e: per-target co-regulator count for all 353 Hif1a CollecTRI targets on the WT_heat contrast, demonstrating that 92% are shared with ≥1 other CollecTRI TF (mean ~22 co-regulators per target).

**How to read:** Each row = one gene (353 total). `target` = gene symbol; `n_other_TFs` = number of other CollecTRI TFs also targeting this gene (Hif1a excluded); `t_wt` = WT_heat limma t-statistic; `co_regulated` = logical, n_other_TFs ≥ 1. The high co-regulation count explains why the multivariate MLM re-distributes Hif1a's shared signal, collapsing its rank from #12 (ULM) to #142 (MLM). Claim tier L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `colors.diverging` | `03_results/04_tf/tables/fig3e_mlm_collinearity_data.csv` |

## tables/_overview/fig3f_consensus_zmix.csv

Plot-ready data for fig3f: per-statistic folded-z scores for Hif1a on the WT_heat CollecTRI run, showing the five-to-one univariate-vs-multivariate vote that re-inflates Hif1a's consensus rank to #8.

**How to read:** 6 rows × 5 columns. `source` = "Hif1a"; `statistic` = ulm / mlm / wsum / norm_wsum / corr_wsum / consensus; `z_score` = folded z (sign-symmetric sd standardisation per statistic); `family` = univariate / multivariate / consensus; `canonical` = logical. The lone MLM bar (z = 0.85) is outvoted by four univariate-family statistics (z = 2.44–4.31); consensus = unweighted mean = 2.89. Claim tier L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `colors.diverging` | `03_results/04_tf/tables/fig3f_consensus_zmix_data.csv` |

## tables/_overview/fig3g_target_decomposition.csv

Plot-ready data for fig3g: signed contribution (sign(mor) × t_wt) for all 353 Hif1a CollecTRI targets ranked by contribution, classified by biological specificity class, plus a `legend_class` column for the heat-shock/stress highlight overlay.

**How to read:** 353 rows. `source` = "Hif1a"; `target` = gene symbol; `mor` = ±1 (mode of regulation from CollecTRI); `t_wt` = WT_heat limma t-statistic; `contrib` = sign(mor) × t_wt (positive = pushes ULM score up); `class` = hypoxic HIF core (4 genes: Pdk1/Bnip3/Bnip3l/Car9, repressed, sum = −8.46) / shared/glycolytic (9 genes, sum = +39.81) / other (340 genes, 91.8% of total); `rank` = rank position in the contribution landscape (1 = most negative); `legend_class` = re-label for the heat-shock/stress highlight. Claim tier L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `design.module_colors` | `03_results/04_tf/tables/fig3g_target_decomposition_data.csv` |

## tables/_overview/fig3i_interaction_primer.csv

Plot-ready data for fig3i: eight group-mean log2CPM values (Ifit1 and Vegfa × genotype × temperature) with pre-computed contrast logFCs and the Interaction BH-adjusted p-value.

**How to read:** 8 rows × ~9 columns. `gene` = Ifit1 / Vegfa; `arm` = cGAS-dependent / no-detectable-cGAS-dependence; `genotype` = WT / cGASKO; `temp` = 37 / 39 (°C); `mean_log2cpm`; `wt_heat` = WT heat logFC; `ko_heat` = KO heat logFC; `interaction` = Interaction logFC; `inter_adjP` = BH-adjusted Interaction p-value. Ifit1: Interaction logFC = +1.007, adjP = 0.003 (diverging slopes = cGAS-dependent); Vegfa: Interaction ~0, adjP ~1.000 (parallel slopes = no detectable cGAS-dependence). Claim tier L3 (n=5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `design.module_colors; colors.diverging` | `03_results/04_tf/tables/fig3i_interaction_primer_data.csv` |

## tables/_overview/fig3j_topTF_allmethods_WT_heat.csv

Plot-ready data for fig3j (WT_heat panel): top-12 activated and suppressed TFs per decoupleR statistic on the WT_heat contrast, across all six statistics; scores are NOT comparable across statistics.

**How to read:** Each row = one TF × statistic combination. `contrast` = WT_heat; `statistic` = ulm / mlm / wsum / norm_wsum / corr_wsum / consensus; `source` = TF symbol; `score`; `rank` (within statistic); `direction` = Up / Down; `axis` = HIF / IFN / other; `key_tf` = logical. Hif1a is consistently high-ranked under five univariate-family statistics (#3–#12) but collapses to #142 under MLM. Compare rank and identity across statistics — not raw scores. Claim tier: methodological comparison.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03b_decoupler_method_comparison_viz.R` | `make_fig3j` | `figures.base_size=16; figures.base_size_column=9; statistics` | `03_results/04_tf/tables/fig3j_allmethods_topTF_data.csv` |

## tables/_overview/fig3j_topTF_allmethods_Interaction.csv

Plot-ready data for fig3j (Interaction panel): top-12 activated and suppressed TFs per decoupleR statistic on the Interaction contrast, across all six statistics; HIF axis uniformly low-ranked across all six methods.

**How to read:** Each row = one TF × statistic combination. `contrast` = Interaction; `statistic` = ulm / mlm / wsum / norm_wsum / corr_wsum / consensus; `source`; `score`; `rank`; `direction`; `axis`; `key_tf`. Hif1a is bottom-half across all six statistics on the Interaction contrast (ulm #423, mlm #248, consensus #307), confirming the IFN-significant / HIF-flat asymmetry is method-robust. Scores not comparable across statistics. Claim tier: methodological comparison.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03b_decoupler_method_comparison_viz.R` | `make_fig3j` | `figures.base_size=16; figures.base_size_column=9; statistics` | `03_results/04_tf/tables/fig3j_allmethods_topTF_data.csv` |

## tables/_overview/fig3k_method_rank_divergence.csv

Plot-ready data for fig3k: per-TF ranks across all six decoupleR statistics for a focused TF set on the WT_heat contrast, with Spearman rho of each statistic vs ULM annotated (MLM rho = 0.619; all others 0.968–0.996).

**How to read:** Each row = one TF × statistic combination. `contrast` = WT_heat (primarily); `source` = TF symbol; `statistic`; `rank` = integer rank among 658 TFs (1 = most activated); `score`. Key rank collapses under MLM: Rela ulm #3 → mlm #350; Nfkb1 ulm #6 → mlm #86; Hif1a ulm #12 → mlm #142. Spearman rho values are stored in the sibling `fig3k_method_rank_spearman.csv`. Claim tier: methodological rank divergence.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03b_decoupler_method_comparison_viz.R` | `save_overview` | `figures.base_size=16; figures.base_size_column=9; statistics` | `03_results/04_tf/tables/fig3k_method_rank_divergence_data.csv` |

## tables/_overview/fig3l_hif_attribution.csv

Plot-ready data for fig3l: per-target signed contributions for the 15 curated Hif1a regulon members bucketed into five biological modules, with the 339 other_unclassified genes represented as an aggregate row.

**How to read:** Each curated row = one gene. `source` = "Hif1a"; `target`; `t_wt` = WT_heat t-statistic (copied from fig3g; no recomputation); `contrib` = sign(mor) × t_wt; `module` = heatshock_stress / shared_angio_glucose / autoreg_feedback / hif1a_hypoxic_core / other_unclassified; `isoform_attribution` = non-HIF (stress) / shared HIF1/HIF2 / HIF1a-preferential-feedback / HIF1a-selective / unclassified; `direction` = Up / Down. Key result: heatshock_stress sum = +108.54 (7 genes, all UP); hif1a_hypoxic_core sum = −8.46 (4 genes, ALL DOWN). Claim tier: descriptive target attribution.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03c_hif_program_attribution_viz.R` | `save_overview` | `figures.base_size=16; figures.base_size_column=9; modules` | `03_results/04_tf/tables/fig3l_hif_attribution_data.csv` |

## tables/_overview/fig3m_ulm_mechanic.csv

Plot-ready data for fig3m: aligned contributions for the Hif1a (n = 353) and Stat2 (n = 32) regulons side by side, with a `is_stress_contaminant` flag identifying the seven heat-shock/stress genes driving Hif1a's rightward pile-up.

**How to read:** 385 rows. `tf` = Hif1a / Stat2; `regulon_class` = promiscuous / specific; `target`; `mor` = ±1; `t_wt` = WT_heat limma t-statistic; `aligned_contrib` = sign(mor) × t_wt; `is_stress_contaminant` = logical (TRUE for {Hspa1a, Timp1, Sdc1, Cdkn1a, Serpine1, Eno2, Spp1} within the Hif1a regulon). Values are copies/joins from fig3g and CollecTRI — no run_ulm called. A high ULM score is literally a rightward pile-up; the Stat2 panel shows a dispersed, near-zero distribution. Claim tier: illustrative mechanic.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03d_ulm_mechanic_viz.R` | `save_overview` | `figures.base_size=16; figures.base_size_column=9` | `03_results/04_tf/tables/fig3m_ulm_mechanic_data.csv` |

## tables/_overview/fig3n_heat_main_regulators.csv

Plot-ready data for fig3n: decoupleR-ULM scores and BH padj for the heat-shock / HIF / IFN axis TFs across heat-MAIN contrasts (Temp_main, WT_heat, KO_heat), showing Hsf1 significantly co-elevated alongside Hif1a and Epas1.

**How to read:** 24 rows. `tf` = Hsf1 / Hif1a / Epas1 / Irf3 / Stat1 / Stat2 / Nfkb1 / Rela; `contrast` = Temp_main / WT_heat / KO_heat; `score` = decoupleR-ULM activity (same `run_ulm(.mor="mor", minsize=5)` call as fig3a/fig3c — scores are axis-comparable); `padj` = BH within contrast; `direction`; `axis` = heatshock / HIF / IFN. Key values: Hsf1 Temp_main +3.20 padj = 0.015 (sig); Hif1a Temp_main +5.14 padj = 2.0×10⁻⁵. Co-elevation only — no crowning of any TF as driver. Claim tier L3.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03e_heat_main_regulators_viz.R` | `save_overview` | `figures.base_size=16; figures.base_size_column=9` | `03_results/04_tf/tables/fig3n_heat_main_regulators_data.csv` |

## tables/_overview/fig3p_heatmain_ranking.csv

Plot-ready data for fig3p: CollecTRI-ULM rankings for the top-22 TFs + Hsf1 on the heat-MAIN (Temp_main) contrast, showing Hif1a at #9 in a crowd of stress / immediate-early / NF-kB regulators with no clean winner.

**How to read:** Each row = one TF (~23 rows: top-22 + Hsf1 appended below a visual break). `source` = TF symbol; `score` = CollecTRI ULM on Temp_main; `rank` = integer rank among 658 TFs (1 = most activated); `p_value`; `family` = curated TF family (HIF axis = orange; AP-1/immediate-early = red; NF-kB = blue; heat-shock = purple; etc.; defined by `FAMILY_COLORS` palette). Key: Jun #1 (6.06), Hif1a #9 (5.14), Hsf1 #50 (3.20). The tiny gaps between ranks mean no single TF is identifiable as the driver. Claim tier: descriptive ranking.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03g_nonidentifiability_viz.R` | `save_overview` | `figures.top_n=22; figures.base_size=16; figures.base_size_column=9` | `03_results/04_tf/tables/fig3p_heatmain_ranking_data.csv` |

## tables/_overview/fig3q_coregulators.csv

Plot-ready data for fig3q: set-membership counts showing that the top sharers of Hif1a's 353-member CollecTRI target set are all non-hypoxia-specific regulators, making the heat-MAIN signal non-attributable to Hif1a alone.

**How to read:** Each row = one TF (top-15 sharers). `source` = TF symbol; `shared_targets` = count of Hif1a's 353 targets also in this TF's CollecTRI regulon; `pct_of_hif1a_set` = shared_targets / 353 × 100; `family` = curated TF family. Top sharers: Sp1 48.4% (171), Trp53 40.5% (143), Nfkb 36.0% (127), Jun 33.4% (118), Myc 31.2% (110). For reference: Epas1 13.0% (46); Hsf1 6.8% (24). No statistical inference — pure set arithmetic over CollecTRI edges. Claim tier: descriptive.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03g_nonidentifiability_viz.R` | `save_overview` | `figures.top_n=15; figures.base_size=16; figures.base_size_column=9` | `03_results/04_tf/tables/fig3q_coregulators_data.csv` |

## tables/_overview/fig3r_shared_ownership.csv

Plot-ready data for fig3r: cross-tabulation of the top-20 heat-MAIN-driving Hif1a targets against Hif1a and its 12 largest target-sharing TFs, with CollecTRI membership indicator and gene/TF-level marginal statistics.

**How to read:** Each row = one gene × TF combination (top-20 genes × 13 TFs = 260 rows). `target` = gene symbol; `tf` = Hif1a or one of its 12 largest sharers (from fig3q); `in_regulon` = 0/1 (CollecTRI membership; teal tile in the heatmap when 1); `gene_heat_t` = gene's heat-MAIN limma t-statistic (left marginal bar); `tf_heatmain_score` = TF's heat-MAIN ULM score (top marginal bar); `contrib` = signed contribution of the gene to Hif1a's ULM score (copied from fig3g). Every gene reads high on heat-MAIN AND is claimed by multiple TFs simultaneously — shared ownership, not unique attribution. Claim tier: descriptive non-identifiability.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03g_nonidentifiability_viz.R` | `save_overview` | `figures.base_size=16; figures.base_size_column=9; patchwork=TRUE` | `03_results/04_tf/tables/fig3r_membership_data.csv` |


## figures/_overview/fig3a_tf_interaction_axes.png

On the cGAS x heat Interaction contrast, the HIF axis (Hif1a/Epas1)
is flat/NS -- no detectable cGAS-dependence -- while IFN members
Irf3/Stat2/Stat1 are interaction-positive.

**How to read:** Lollipop x = TF activity (ULM score; Interaction = WT_heat -
KO_heat); color/shape = TF axis (orange triangle = HIF; blue square =
IFN/NFkB; grey = other). Glyph: * = BH padj < 0.05. Asymmetry, not
proven independence; claim tier L3 (n=5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `colors.diverging; design.axis_colors` | `03_results/04_tf/tables/fig3a_tf_interaction_axes_data.csv` |

## figures/_overview/fig3b_top_tf_by_contrast.png

Top activated/suppressed TFs (CollecTRI ULM) per contrast on a shared
ULM-score x-axis: reads out what each heat arm, the temperature main
effect, and the cGAS x heat Interaction capture.

**How to read:** Four facets share the x-axis (TF activity, ULM score). Lollipop color
= TF axis (orange = HIF; blue = IFN/NFkB; grey = other); HIF/IFN
watchlist TFs labelled. Capped to top-20 per direction per contrast
by |score|. Claim tier L3 (n=5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `figures.top_n=20; design.axis_colors` | `03_results/04_tf/tables/fig3b_top_tf_by_contrast_data.csv` |

## figures/_overview/fig3c_hif1a_rank_cascade.png

Hif1a's #1 DoRothEA ranking is method/network-fragile: it collapses
under the CollecTRI swap (#12) and multivariate MLM de-confounding
(#142), then is partially re-inflated by the consensus (#8).

**How to read:** x = inference method/network config (left -> right reading order); y
= Hif1a rank among activated TFs (reversed; rank #1 at top = most
activated). Point color = ULM/consensus score (orange up / blue
down). Labels give rank/total and score. Tier L3 (n=5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `colors.diverging` | `03_results/04_tf/tables/fig3c_hif1a_rank_cascade_data.csv` |

## figures/_overview/fig3d_regulon_swap.png

The DoRothEA -> CollecTRI regulon swap dilutes Hif1a's ULM signal:
the CollecTRI-only targets carry high |t| but mixed sign, so they
cancel and drag the regulon mean down (explains the #1 -> #12 drop).

**How to read:** Box+jitter of per-target WT_heat t-statistic by regulon membership
(DoRothEA-only / shared / CollecTRI-only). The CollecTRI-only box
straddles 0 (mixed sign), diluting the mor-weighted mean. Annotation
gives regulon sizes. Tier L3 (n=5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `design.module_colors` | `03_results/04_tf/tables/fig3d_regulon_swap_data.csv` |

## figures/_overview/fig3e_mlm_collinearity.png

Nearly every Hif1a target is co-regulated by many other CollecTRI
TFs, so the multivariate MLM re-attributes Hif1a's signal away from
it -- collapsing its score (explains the #12 -> #142 drop).

**How to read:** Histogram: for each Hif1a target, how many OTHER CollecTRI TFs also
regulate it. Dashed orange line = mean. Annotation gives the ULM ->
MLM score/rank collapse. High co-regulation is why MLM de-confounds
Hif1a's signal away. Tier L3 (n=5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `colors.diverging` | `03_results/04_tf/tables/fig3e_mlm_collinearity_data.csv` |

## figures/_overview/fig3f_consensus_zmix.png

Only the multivariate MLM gives Hif1a a low folded-z; the four
univariate statistics stay high and the consensus (their mean)
re-inflates Hif1a back to #8 (explains the #142 -> #8 recovery).

**How to read:** Bars = Hif1a folded z per decoupleR statistic; fill = statistic
family (orange univariate / blue multivariate / grey consensus).
Dashed line = consensus mean z. The lone low MLM bar is outvoted by
the univariate family. Tier L3 (n=5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `colors.diverging` | `03_results/04_tf/tables/fig3f_consensus_zmix_data.csv` |

## figures/_overview/fig3g_target_decomposition.png

Hif1a's positive WT_heat ULM score is built almost entirely from its
'other' regulon members -- generic stress/ECM genes tower over the
HIF-specific core, which is repressed. This is a heat-induced
glycolytic/stress program partially overlapping HIF targets, not the
canonical HIF program.

**How to read:** Lollipop landscape: all 353 Hif1a regulon members ranked by signed
contribution (sign(mor) x t_wt; + = pushes the ULM score up). Color =
regulon class (purple heat-shock/stress; orange shared/glycolytic;
teal repressed hypoxic core; grey other). Labelled = top stress
drivers + the 7 HIF-specific members (the source table). Tier L3
(n=5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `design.module_colors` | `03_results/04_tf/tables/fig3g_target_decomposition_data.csv` |

## figures/_overview/fig3i_interaction_primer.png

Reading the cGAS x heat Interaction as the difference of genotype
heat-slopes: Ifit1 (IFN arm) slopes diverge -> large positive
Interaction (cGAS-dependent); Vegfa (shared HIF/glycolytic target)
slopes stay parallel -> Interaction ~0 (no detectable
cGAS-dependence).

**How to read:** Two facets (free y; same metric log2 CPM): x = temperature, lines =
genotype (WT teal solid, cGAS-KO orange dashed). Diverging slopes =
cGAS-dependent; parallel slopes = no detectable cGAS-dependence.
Strip text carries the contrast values (WT heat / KO heat /
Interaction + adjP). Tier L3 (n=5).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03_decoupler_tf_viz.R` | `save_overview` | `design.module_colors; colors.diverging` | `03_results/04_tf/tables/fig3i_interaction_primer_data.csv` |

## figures/_overview/fig3j_topTF_allmethods_WT_heat.png

Method comparison for contrast WT_heat: lollipops show the top
activated / suppressed TFs under all 6 decoupleR statistics. Sample
mapping owner-confirmed (GSE329522 iTreg 2x2 genotype x temperature)

**How to read:** Lollipops = per-TF activity, faceted by statistic (ULM, MLM, wsum,
norm_wsum, corr_wsum, consensus), colored by axis (HIF orange;
IFN/NFkB blue; other grey). Key TFs are labelled. Note that y-axes
differ (scales='free') to compare rank. Claim tier: methodological
comparison (ranking identity).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03b_decoupler_method_comparison_viz.R` | `make_fig3j` | `figures.base_size=16; figures.base_size_column=9; statistics` | `03_results/04_tf/tables/fig3j_allmethods_topTF_data.csv` |

## figures/_overview/fig3j_topTF_allmethods_Interaction.png

Method comparison for contrast Interaction: lollipops show the top
activated / suppressed TFs under all 6 decoupleR statistics. Sample
mapping owner-confirmed (GSE329522 iTreg 2x2 genotype x temperature)

**How to read:** Lollipops = per-TF activity, faceted by statistic (ULM, MLM, wsum,
norm_wsum, corr_wsum, consensus), colored by axis (HIF orange;
IFN/NFkB blue; other grey). Key TFs are labelled. Note that y-axes
differ (scales='free') to compare rank. Claim tier: methodological
comparison (ranking identity).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03b_decoupler_method_comparison_viz.R` | `make_fig3j` | `figures.base_size=16; figures.base_size_column=9; statistics` | `03_results/04_tf/tables/fig3j_allmethods_topTF_data.csv` |

## figures/_overview/fig3k_method_rank_divergence.png

MLM (the only multivariate estimator) reshuffles both axes on
WT_heat: rank heatmap shows MLM is a structural outlier with a
Spearman rho of 0.62 vs ULM, while univariate statistics sit at
0.97-1.00. Top univariate TFs are demoted under MLM across both axes:
Rela 3->350, Nfkb1 6->86, Hif1a 12->142. Hif1a's collapse is one
instance of that general de-confounding reshuffle, not a HIF-specific
quirk. Sample mapping owner-confirmed (GSE329522 iTreg 2x2 genotype x
temperature)

**How to read:** Rank heatmap: rows = key TFs, columns = statistics. Fill = binned
rank quantile (bright yellow = highly activated, dark purple =
demoted/inactive). Cell label = absolute rank. The boxed column marks
MLM (multivariate outlier). Claim tier: methodological rank
divergence.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03b_decoupler_method_comparison_viz.R` | `save_overview` | `figures.base_size=16; figures.base_size_column=9; statistics` | `03_results/04_tf/tables/fig3k_method_rank_divergence_data.csv; fig3k_method_rank_spearman.csv` |

## figures/_overview/fig3l_hif_attribution.png

Regrouping Hif1a's WT_heat target contributions by curated module
splits the positive signal into a heat-shock/stress fraction summing
to +108.54, all seven members up, against a HIF-annotated remainder
of +14.27 (shared angio/glucose +15.93, autoregulatory feedback
+6.81, HIF1a-selective hypoxic core -8.46). Pdk1, Bnip3, Bnip3l and
Car9 are canonical HIF1a-induced targets, yet all four sit left of
zero on the contrast where Hif1a scores positive. The positive score
is thus consistent with a heat-induced stress/glycolytic program
overlapping the HIF1a regulon, not with canonical hypoxic HIF1a
output. Correlative: attribution of an inferred activity score, not
evidence about HIF1a protein. Sample mapping owner-confirmed
(GSE329522 iTreg 2x2 genotype x temperature)

**How to read:** Horizontal lollipop, one row per gene. x = the signed contribution
sign(mor) x t_wt this gene makes to Hif1a's WT_heat ULM score; the
grey rule marks zero, so points right of it are genes the contrast
moves up and points left of it genes it moves down, on a symmetric
scale. Facets are the four curated modules, ordered
heat-shock/stress, shared angio/glucose, autoregulatory feedback,
HIF1a-selective hypoxic core; each strip carries that module's gene
count and summed contribution. Colour encodes module only
(config.R::MODULE_COLORS). Per-gene isoform attributions and the 339
unclassified regulon members sit in fig3l_hif_attribution_data.csv.
Claim tier: descriptive attribution; no statistics are computed here.
The seven heat-shock genes here are carved out of the coarser 91.8%
'other' lump in fig3g, leaving a 339-member 90.7% residual: the two
figures differ by construction, not by error.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03c_hif_program_attribution_viz.R` | `save_overview` | `figures.title_size=16; figures.subtitle_size=11; MODULE_COLORS; MODULE_LABELS` | `03_results/04_tf/tables/fig3l_hif_attribution_data.csv` |

## figures/_overview/fig3m_ulm_mechanic.png

decoupleR-ULM scores a TF by the regulon-weighted pile-up of aligned
target contributions: a promiscuous regulon (Hif1a) scores high off
generic heat-shock/stress contaminants it happens to contain (a
heat-induced glycolytic/stress program partially overlapping HIF
targets), while a small specific regulon (Stat2) does not pile up.
This is the scoring MECHANIC, not a biology claim; the Hif1a panel is
NOT 'the HIF program'/'Hif1a the TF'. Sample mapping owner-confirmed
(GSE329522 iTreg 2x2 genotype x temperature)

**How to read:** Two panels share the x-axis = aligned (signed) target contribution =
sign(mor) x t_gene. Each jittered point is one regulon member; purple
= heat-shock/stress contaminant, orange = HIF-specific member, grey =
other. The grey diamond + dashed rule sit at the regulon's weighted
center (the ULM score); a high score is literally a rightward
pile-up. Pedagogical (claim tier: illustrative mechanic, not a
quantitative result).

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03d_ulm_mechanic_viz.R` | `save_overview` | `figures.base_size=16; figures.base_size_column=9` | `03_results/04_tf/tables/fig3m_ulm_mechanic_{data,summary}.csv` |

## figures/_overview/fig3n_heat_main_regulators.png

On the heat-main contrasts the heat-shock regulator Hsf1 is
significantly co-elevated (Temp_main 3.20, padj 0.015) alongside the
HIF axis (Hif1a +5.14 / Epas1 +4.17) -- co-elevation only, NOT a
claim that Hsf1 causes or outranks HIF, and no single master TF is
crowned. The axis is the decoupleR-ULM score, the same estimator as
fig3a/fig3c and cross-quotable with them. Sample mapping
owner-confirmed (GSE329522 iTreg 2x2 genotype x temperature)

**How to read:** Lollipops = per-TF activity, faceted by heat-MAIN contrast, colored
by axis (heat-shock purple = Hsf1; HIF orange; IFN blue; other grey).
x = decoupleR-ULM score (CollecTRI), the SAME estimator as fig3a/3c
-- cross-quotable, no GSEA-NES caveat. * = BH padj < 0.05. Claim
tier: descriptive co-elevation; the equal visual weight crowns no
master TF.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03e_heat_main_regulators_viz.R` | `save_overview` | `figures.base_size=16; figures.base_size_column=9` | `03_results/04_tf/tables/fig3n_heat_main_regulators_data.csv` |

## figures/_overview/fig3p_heatmain_ranking.png

On heat-MAIN there is no clean winner: Hif1a is #9 in a crowd of
co-elevated stress / immediate-early / NF-kB TFs separated by tiny
gaps (Jun 6.06 down to Nfkb1 5.12, all p~1e-7), it does not sit atop
a hypoxia-specific peak, and the canonical heat-shock TF Hsf1 is far
down at #50. This ranking crowns no TF -- not Hif1a, not Jun/AP-1,
not Epas1/HIF2a. Sample mapping owner-confirmed (GSE329522 iTreg 2x2
genotype x temperature)

**How to read:** Ranked dotplot: each row is a TF, x = heat-MAIN ULM activity score,
colored by curated TF family. Top-22 shown; Hsf1 appended below the
dashed discontinuity (it is #50). The dark ring marks Hif1a (#9).
Claim tier: descriptive -- the tiny gaps mean no single TF is
identifiable as the driver.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03g_nonidentifiability_viz.R` | `save_overview` | `figures.top_n=22; figures.base_size=16; figures.base_size_column=9` | `03_results/04_tf/tables/fig3p_heatmain_ranking_data.csv` |

## figures/_overview/fig3q_coregulators.png

92% of Hif1a's 353 CollecTRI targets are co-regulated (mean 22 other
TFs each); the top sharers are the network's most promiscuous
regulators (Sp1, Trp53, NF-kB, AP-1, Myc), none hypoxia-specific --
so Hif1a's heat-MAIN signal cannot be attributed to Hif1a alone, and
no co-regulator is crowned the driver. Sample mapping owner-confirmed
(GSE329522 iTreg 2x2 genotype x temperature)

**How to read:** Horizontal bars (top-15 sharers): x = % of Hif1a's 353 targets that
this TF also regulates (shared count in parentheses), colored by TF
family. None of the top sharers is hypoxia-specific. Claim tier:
descriptive -- targets 'belong to everyone', so the signal is not
Hif1a-specific.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03g_nonidentifiability_viz.R` | `save_overview` | `figures.top_n=15; figures.base_size=16; figures.base_size_column=9` | `03_results/04_tf/tables/fig3q_coregulators_data.csv` |

## figures/_overview/fig3r_shared_ownership.png

The same heat-driven genes populate many TFs' CollecTRI regulons:
every gene (LEFT bar) and every TF (TOP bar) reads high on heat-MAIN,
so the contrast cannot single out Hif1a -- or any one TF -- as the
driver. None is identifiable. Sample mapping owner-confirmed
(GSE329522 iTreg 2x2 genotype x temperature)

**How to read:** Center grid: rows = top-20 heat-MAIN-driving genes of Hif1a's regulon
(signed contribution order), columns = Hif1a + its 12 largest
target-sharers (heat-MAIN ULM score order; Hif1a orange + bold). Teal
tile = the gene is in that TF's regulon. LEFT bar = gene heat-MAIN t;
TOP bar = TF heat-MAIN score (both read high). Claim tier:
descriptive non-identifiability -- shared ownership, no single
driver.

| Script | Function | Config | Input |
|---|---|---|---|
| `02_analysis/scripts/03g_nonidentifiability_viz.R` | `save_overview` | `figures.base_size=16; figures.base_size_column=9; patchwork=TRUE` | `03_results/04_tf/tables/fig3r_membership_data.csv` |

