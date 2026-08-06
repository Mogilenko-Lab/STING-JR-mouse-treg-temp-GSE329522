# 04_tf — Which transcription factor the heat response nominates

Transcription-factor activity is inferred with decoupleR from the limma moderated-t vectors of
[`../03_de/`](../03_de/). The primary run is a univariate linear model over the CollecTRI mouse
regulons; DoRothEA-ABC, a multivariate model and a consensus statistic run alongside as
comparators.

**An inferred activity is a statistic over target-gene expression.** The model regresses each
gene's contrast statistic on the mode of regulation the network records for that edge and reports
the slope's t-statistic. It describes how the genes a network files under a factor behave in a
contrast. HIF1α protein, its stabilisation and its occupancy sit outside what it measures, and
every artifact here is annotation tier for that reason.

Two threads run through the directory.

**The two arms are asymmetric.** On the interaction contrast the interferon members carry
significance — Irf3 +5.93 (adjusted p 1.0e-06), Stat2 +5.85 (1.1e-06), Stat1 +5.30 (1.9e-05) —
while the HIF axis is flat (Hif1a −0.19 at 0.97, Epas1 −0.86 at 0.91), so the interferon arm is
cGAS-dependent and the HIF arm has no detectable cGAS-dependence at n = 5
(`tables/fig3a_tf_interaction_axes_data.csv`).

**Hif1a's high rank is a property of the method and the network.** On `WT_heat` its rank runs
#1 of 263 under DoRothEA-ULM, #12 of 658 under CollecTRI-ULM, #142 of 658 under CollecTRI-MLM
and #8 of 658 under the consensus (`tables/fig3c_hif1a_rank_cascade_data.csv`). Three
measurements explain that traverse and bound what the score can support: the CollecTRI regulon
is 2.7× larger and dilutes a tighter DoRothEA one; 92% of its 353 targets are shared with a mean
of 22 other factors, which the multivariate fit redistributes; and 91.8% of the signed
contribution comes from 340 members that generic stress and remodelling genes dominate, while
the four genes diagnostic of canonical HIF1α output — Pdk1, Bnip3, Bnip3l, Car9 — are all
repressed, summing to −8.46 (`tables/fig3g_target_decomposition_summary.csv`).

So the score reports a heat-induced stress and glycolytic program that overlaps the HIF1α
regulon. Naming a single driver takes a separate measurement, and `fig3p`, `fig3q` and `fig3r`
are the three panels that show why this contrast supplies none.

**Compute** — `03_decoupler_tf.R` (primary run, all `fig3a`–`fig3i` tables),
`03b_decoupler_method_comparison.R` (`fig3j`, `fig3k`), `03c_hif_program_attribution.R`
(`fig3l`), `03d_ulm_mechanic.R` (`fig3m`), `03e_heat_main_regulators.R` (`fig3n`),
`03g_nonidentifiability.R` (`fig3p`–`fig3r`). Each has a `_viz.R` sibling that reads the tables
and recomputes nothing.

---

## Figures

All seventeen sit in `figures/_overview/` as a `.pdf` and `.png` pair, with the same-stem source
table in `tables/_overview/`. `deck_assets/` holds presentation copies of the same PNGs.

### `fig3a_tf_interaction_axes` — the two arms on the cGAS-dependence test

**The interferon axis carries the interaction and the HIF axis is flat.**
Lollipop over the interaction contrast. x, univariate activity score; y, transcription factor.
Colour and shape give the axis: orange triangles for HIF (Hif1a, Epas1), blue squares for
IFN/NF-κB (Irf1, Irf3, Irf7, Stat1, Stat2, Nfkb1, Rela), grey for everything else. An asterisk
marks BH-adjusted p < 0.05. Nfkb1 (+2.16, p 0.60) and Rela (+2.03, p 0.63) lean positive and
stay above threshold at n = 5.
*Source* `tables/_overview/fig3a_tf_interaction_axes.csv`.

### `fig3c_hif1a_rank_cascade` — the traverse

**Hif1a's rank moves with the inference choice.**
x, the four method-and-network configurations in reading order; y, Hif1a's rank among activated
factors, inverted so rank #1 sits at the top. Point colour gives the score. Labels carry
rank/total and the score behind it. Reading left to right: DoRothEA-ULM #1 of 263 (6.09),
CollecTRI-ULM #12 of 658 (5.11), CollecTRI-MLM #142 of 658 (1.13), CollecTRI-consensus #8 of 658
(2.89). A signal stable across estimator and network would hold its rank; this one traverses 141
places. Its interaction score stays flat in all four.
*Source* `tables/_overview/fig3c_hif1a_rank_cascade.csv`.

### `fig3d_regulon_swap` — why the rank was #1

**Regulon expansion dilutes the score.**
Box and jitter of per-target `WT_heat` moderated t by regulon membership. x, membership class
(DoRothEA-only n = 30, shared n = 101, CollecTRI-only n = 252); y, moderated t. The
CollecTRI-only box straddles zero — high |t| with mixed sign, median signed t ≈ 0.3 — so those
252 targets add magnitude and little net direction, which drags the weighted mean down.
DoRothEA's 353-versus-131 smaller regulon concentrates the high-t heat targets and produces the
#1 rank; CollecTRI's fuller one gives #12. The rank is a property of the library.
*Source* `tables/_overview/fig3d_regulon_swap.csv`, with per-network summaries in
`tables/fig3d_regulon_swap_summary.csv`.

### `fig3e_mlm_collinearity` — why the rank collapses

**Nearly every Hif1a target belongs to many regulons at once.**
Histogram over Hif1a's 353 CollecTRI targets. x, how many other CollecTRI factors also regulate
that target; y, target count. The dashed orange rule marks the mean, 22.23; the median is 11 and
92.07% of targets carry at least one other claimant. A multivariate fit distributes shared
movement across competing regulons, which takes the score from 5.11 at rank #12 to 1.13 at rank
#142 on the same network and the same data.
*Source* `tables/_overview/fig3e_mlm_collinearity.csv`, with the two-bar collapse in
`tables/fig3e_score_collapse.csv` and the summary in `tables/fig3e_mlm_collinearity_summary.csv`.

### `fig3g_target_decomposition` — where the score comes from

**The full 353-member landscape, ranked by what each target contributes.**
Lollipop over every regulon member. x, signed contribution, `sign(mor) × t_wt`, so a point right
of zero pushes the score up; y, rank order. Colour gives the coarse class: purple for the
heat-shock and stress highlight, orange for shared/glycolytic, teal for the repressed hypoxic
core, grey for the rest. The six largest contributors are labelled (Timp1 +24.92, Sdc1 +22.89,
Cdkn1a +12.33, Serpine1 +11.21, Eno2 +11.00, Hspa1a +10.81), as are the curated-class members.

Across the regulon: the 340 "other" targets sum to +351.43, which is 91.8% of the signed total;
the nine shared/glycolytic members (Hk2, Vegfa, Egln3, Slc2a1, Aldoa, Pgk1, Ldha, Eno1, Pkm) sum
to +39.81; the four-member hypoxic core sums to −8.46 with all four repressed. Vegfa, Slc2a1 and
Egln3 are shared targets rather than HIF-diagnostic ones.
*Source* `tables/_overview/fig3g_target_decomposition.csv`, plus
`tables/fig3g_target_decomposition_summary.csv` for the class provenance split and
`tables/fig3g_target_landscape_labels.csv` for the labelled subset.

### `fig3l_hif_attribution` — what kind of genes, and in which direction

**The genes that diagnose HIF1α move against the score.**
Horizontal lollipop, one row per gene, faceted by curated module in the order heat-shock/stress,
shared angio/glucose, autoregulatory feedback, HIF1α-selective hypoxic core. x, signed
contribution on the same axis `fig3g` uses; a grey rule marks zero, so points right of it are
genes the contrast raises. Each facet strip carries its module's gene count and summed
contribution.

The seven heat-shock members (Timp1, Sdc1, Spp1, Cdkn1a, Serpine1, Eno2, Hspa1a) all rise and
sum to +108.54. The shared angiogenic/glycolytic pair (Vegfa, Slc2a1) reaches +15.93 and the
autoregulatory member (Egln3) +6.81. Pdk1, Bnip3, Bnip3l and Car9 are canonical HIF1α-**induced**
targets and all four sit left of zero, summing to −8.46, on the contrast where the factor scores
positive. The four curated modules hold 14 of 353 members and 9.3% of the contribution magnitude,
so the partition is a reading aid over a tenth of the score.

`fig3g` and `fig3l` decompose the same regulon at two grains: `fig3l` carves the seven heat-shock
genes out of `fig3g`'s 340-member lump, leaving a 339-member residual, so the two differ by
construction. Per-gene contributions are identical between them.
*Source* `tables/_overview/fig3l_hif_attribution.csv`, per gene in
`tables/fig3l_hif_attribution_data.csv` (353 rows including the 339 unclassified), module totals
in `tables/fig3l_module_summary.csv`.

### `fig3m_ulm_mechanic` — how the estimator nominates a factor

**A promiscuous regulon scores high when a few large targets pile to one side.**
Two panels on one shared x axis: aligned signed contribution, `sign(mor) × t_gene`. Each jittered
point is one regulon member; purple marks the seven heat-shock contaminants, orange the
HIF-specific members, grey the rest. A grey diamond with a dashed rule sits at the regulon's
weighted centre, which is the score.

Panel A is Hif1a — 353 members, centre +1.08 — and its rightward pile-up. Panel B is Stat2 — 32
members, centre −0.62 — a small specific regulon whose targets spread across both sides. The
scoring mechanic is drawn once here so the landscape panels read.
*Source* `tables/_overview/fig3m_ulm_mechanic.csv`, per gene in
`tables/fig3m_ulm_mechanic_data.csv` (353 Hif1a + 32 Stat2 rows), aggregates in
`tables/fig3m_ulm_mechanic_summary.csv`.

### `fig3n_heat_main_regulators` — the heat-shock factor

**Hsf1 is co-elevated on the heat contrasts alongside the HIF axis.**
Lollipop faceted by heat-main contrast. x, univariate score on the same estimator `fig3a` and
`fig3c` use, so the three panels are cross-quotable. Colour gives the axis: purple for Hsf1,
orange for HIF, blue for IFN, grey for other; an asterisk marks BH-adjusted p < 0.05. On
`Temp_main` Hsf1 reaches +3.20 (p 0.015) beside Hif1a +5.14 (2.0e-05) and Epas1 +4.17 (8.9e-04),
and it holds in both genotypes separately (`WT_heat` +2.83 at p 0.034, `KO_heat` +3.48 at 0.008).

The reading is co-elevation. Hif1a's regulon overlaps Hsf1's thinly — 24 to 32 shared targets,
7–12% of the signed contribution — and of the named heat-shock drivers only Hspa1a, Cdkn1a and
Serpine1 are shared, while Timp1, Sdc1, Spp1 and Eno2 are Hif1a-only. Hsf1's cGAS-by-heat
interaction reaches raw p 0.022 and adjusted p 0.515, so it carries no cGAS-dependence.
*Source* `tables/_overview/fig3n_heat_main_regulators.csv`.

### `fig3p` · `fig3q` · `fig3r` — the non-identifiability triptych

Three panels making one point: on the heat-main contrast this design nominates no single factor.

**`fig3p_heatmain_ranking` — the crowd.** Ranked dotplot of the top 22 factors with Hsf1 appended
below a dashed break. x, heat-main univariate score; y, factor. Colour gives the curated family
(heat-shock purple, AP-1/immediate-early red, NF-κB blue, HIF orange) and a dark ring marks
Hif1a. Hif1a sits #9 of 658 at 5.14 inside a band running Jun 6.06, Egr1 5.54, Fos 5.46, Sp1
5.29, Stat5a 5.27, Rela 5.25, Ncoa1 5.25, Hoxd3 5.20, then Nfkb1 5.12 and Jund 5.10 — about 0.9
score units across ranks #1 to #12, all near p 1e-07. Hsf1 sits #50 at 3.20.
*Source* `tables/_overview/fig3p_heatmain_ranking.csv`.

**`fig3q_coregulators` — shared ownership.** Horizontal bars, top fifteen sharers. x, the
percentage of Hif1a's 353 targets that this factor also regulates, with the shared count in
parentheses; colour gives the family. The largest sharers are generic: Sp1 171 (48.4%), Trp53
143 (40.5%), Jun 118 (33.4%), Nfkb 127 (36.0%), Myc 110 (31.2%), Rela 102 (28.9%). The other HIF
member Epas1 shares 46 (13.0%) and Hsf1 shares 24 (6.8%), so the heaviest sharers are all
outside the hypoxia axis.
*Source* `tables/_overview/fig3q_coregulators.csv`.

**`fig3r_shared_ownership` — the same genes in many regulons.** Tile heatmap with marginals.
Rows, the top twenty Hif1a targets by signed heat-main contribution (Timp1 24.92, Sdc1 22.89,
Itga5 18.46, Glb1 17.39, Lgals3 16.80 … Rora 10.95); columns, Hif1a plus its twelve largest
sharers in heat-main score order, Hif1a flagged orange. A teal tile marks membership. The left
marginal gives each gene's heat-main moderated t and the top marginal each factor's heat-main
score, and both read high, so the same driving genes are claimed by every column.
*Source* `tables/_overview/fig3r_shared_ownership.csv`, membership grid in
`tables/fig3r_membership_data.csv`.

### Supporting panels

| Figure | What it shows |
|---|---|
| `fig3b_top_tf_by_contrast` | Top and bottom twelve factors per contrast on a shared score axis, four facets. Heat raises a broad stress and glycolytic set in both genotypes; the interaction isolates the IFN/Irf9/Stat2 axis. |
| `fig3f_consensus_zmix` | Hif1a's folded z per decoupleR statistic, with fill by statistic family. The four univariate statistics (ulm 2.44, wsum 4.13, norm_wsum 2.72, corr_wsum 4.31) outnumber the single multivariate one (mlm 0.85) in an unweighted mean, which is the arithmetic behind the #142 → #8 recovery. |
| `fig3i_interaction_primer` | The interaction read as the difference of two genotype slopes. x, temperature; lines, genotype; facets, gene with free y. Vegfa slopes stay parallel (WT +0.970, KO +0.954, interaction +0.016 at adjusted p ≈ 1.000); Ifit1 slopes diverge (+0.497, −0.511, +1.007 at adjusted p 0.003). |
| `fig3j_topTF_allmethods_WT_heat` · `fig3j_topTF_allmethods_Interaction` | The full six-statistic audit, faceted by statistic with free x. On `WT_heat` Hif1a holds ranks #3–#12 under the five univariate statistics and #142 under MLM; on the interaction it sits #123–#423 under all six. |
| `fig3k_method_rank_divergence` | Rank heatmap over key factors × six statistics, fill by binned rank quantile, cell label the absolute rank. Spearman against ULM over all 658 factors: MLM 0.619, every other statistic 0.968–0.996. MLM demotes Rela #3 → #350 and Nfkb1 #6 → #86 alongside Hif1a #12 → #142, so the collapse is a general reshuffle of high-collinearity factors. |

---

## Tables

Every figure has a same-stem source table under `tables/_overview/`. The `tables/` root holds
the compute-side outputs those are cut from, plus four files with no panel of their own.

| File | What it holds |
|---|---|
| `fig3a_tf_interaction_axes_data.csv` | Per-factor interaction score, p, adjusted p and axis assignment. |
| `fig3b_top_tf_by_contrast_data.csv` | Top and bottom twelve per contrast, 96 rows. |
| `fig3c_hif1a_rank_cascade_data.csv` | One row per configuration: score, rank, total factors scored, percentage rank. |
| `fig3d_regulon_swap_{data,summary}.csv` | Per-target membership class with its `WT_heat` t; per-network mean\|t\|, median t and %\|t\|>2. |
| `fig3e_mlm_collinearity_{data,summary}.csv`, `fig3e_score_collapse.csv` | Per-target co-regulator count; the 0.9207 fraction and 22.23 mean; the two-bar ULM-to-MLM collapse. |
| `fig3f_consensus_zmix_data.csv` | Hif1a's folded z per statistic; the algebraic mean reproduces the stored consensus to under 1e-06. |
| `fig3g_target_decomposition_{data,summary}.csv`, `fig3g_target_landscape_labels.csv` | All 353 members with `contrib = sign(mor) × t_wt` and coarse class; the per-class sums; the head/tail slice the panel labels. |
| `fig3h_lombardi_vs_phylo_data.csv` | A comparison retained without a panel. The Lombardi 2022 pan-cancer HIF consensus is used elsewhere as a membership lens over symbols and is scored as a gene-set database by no stage in these results. |
| `fig3i_interaction_primer_data.csv` | Four group means per gene plus the contrast logFCs and interaction adjusted p. |
| `fig3j_allmethods_topTF_data.csv` | Top and bottom twelve per statistic for both plotted contrasts. |
| `fig3k_method_rank_divergence_data.csv`, `fig3k_method_rank_spearman.csv` | Per-factor rank per statistic over 658 factors; the five Spearman correlations against ULM. |
| `fig3l_hif_attribution_data.csv`, `fig3l_module_summary.csv` | Per gene with `module`, `isoform_attribution` and `direction`, including the 339 the panel leaves off; per-module counts and sums. |
| `fig3m_ulm_mechanic_{data,summary}.csv` | 353 Hif1a plus 32 Stat2 rows with `aligned_contrib` and the stress-contaminant flag; per-factor aggregates. |
| `fig3n_heat_main_regulators_data.csv` | Axis-selected factors × three heat-main contrasts with score, p and adjusted p. |
| `fig3p_heatmain_ranking_data.csv` | All 658 factors on `Temp_main` with score, rank and raw p. |
| `fig3q_coregulators_data.csv` | Per factor: `shared_targets` against Hif1a's 353, and `pct_of_hif1a_set`. |
| `fig3r_membership_data.csv` | The 20 × 13 membership grid with both marginals carried on it. |

`deck_assets/` holds a PNG copy of each panel for presentation use; the canonical renders are
in `figures/_overview/`.
